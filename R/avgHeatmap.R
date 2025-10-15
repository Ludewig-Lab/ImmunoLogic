#' Average Expression Heatmap for Seurat Objects
#'
#' Generates a heatmap of averaged expression values for selected genes across cell groups or identities in a Seurat object.
#' If no genes are provided, the top variable features will be used.
#'
#' @param seurat A Seurat object.
#' @param selGenes Character vector of genes to plot, or a data.frame with a column named `gene` or `geneID`. Defaults to `NULL`.
#' @param group_by Metadata column in the Seurat object to group cells by. Defaults to `NULL` (uses active identity).
#' @param condition_by Optional metadata column for a second grouping variable (e.g., treatment condition). Default is `NULL`.
#' @param scale_method Scaling method for the heatmap: `"row"`, `"column"`, or `"none"`. Default is `"row"`.
#' @param cluster_rows Logical, whether to cluster rows. Default is `FALSE`.
#' @param cluster_cols Logical, whether to cluster columns. Default is `FALSE`.
#' @param show_rownames Logical, whether to display row names (genes). Default is `TRUE`.
#' @param show_colnames Logical, whether to display column names (groups). Default is `TRUE`.
#' @param cellwidth Numeric, width of each cell in the heatmap. Default is `15`.
#' @param cellheight Numeric, height of each cell in the heatmap. Default is `10`.
#' @param color_palette Color vector for the heatmap. Default is a blue-white-red gradient.
#' @param annotation_colors Named list of colors for annotations. Can contain elements named 'Group' and 'Condition'.
#'   Alternatively, provide a named vector which will be used for the Group annotation.
#' @param condition_colors Named vector of colors for condition annotation. Only used if `condition_by` is specified.
#' @param gaps_row Optional vector specifying rows after which to insert gaps. Default is `NULL`.
#' @param gaps_col Optional vector specifying columns after which to insert gaps. Default is `NULL`.
#' @param n_variable_genes Number of top variable genes to use if `selGenes` is `NULL`. Default is `20`.
#' @param ... Additional arguments passed to [pheatmap::pheatmap()].
#'
#' @return A `pheatmap` object.
#'
#' @author Mechthild Lütge, Roman Stadler
#'
#' @details
#' - Handles gene names as symbols or ENSEMBL IDs.
#' - Automatically removes genes with zero variance across groups.
#' - Orders genes by cluster of maximum expression if `cluster_rows = FALSE`.
#' - Automatically selects color palettes based on group names and number of groups.
#' - When `condition_by` is specified, creates combined groups (e.g., "Cluster0_ConditionA").
#'
#' @importFrom Seurat VariableFeatures Idents GetAssayData
#' @importFrom dplyr mutate left_join select group_by summarise
#' @importFrom pheatmap pheatmap
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' seurat <- YourSeuratObject
#'
#' # Simple heatmap by cell type
#' avgHeatmap(seurat, selGenes = c("GeneA", "GeneB"), group_by = "celltype")
#'
#' # Heatmap with cell type and condition
#' avgHeatmap(seurat,
#'            selGenes = c("GeneA", "GeneB"),
#'            group_by = "celltype",
#'            condition_by = "treatment")
#'
#' # With custom colors for both annotations
#' cluster_cols <- c("Cluster0" = "red", "Cluster1" = "blue")
#' condition_cols <- c("WT" = "forestgreen", "Mutant" = "firebrick")
#' avgHeatmap(seurat,
#'            annotation_colors = cluster_cols,
#'            condition_colors = condition_cols,
#'            condition_by = "treatment")
#' }
#'
#' @export

avgHeatmap <- function(seurat,
                       selGenes = NULL,
                       group_by = NULL,
                       condition_by = NULL,
                       scale_method = "row",
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       cellwidth = 15,
                       cellheight = 10,
                       color_palette = NULL,
                       annotation_colors = NULL,
                       condition_colors = NULL,
                       gaps_row = NULL,
                       gaps_col = NULL,
                       n_variable_genes = 20,
                       # Old argument names for backward compatibility
                       colVecIdent = NULL,
                       colVecCond = NULL,
                       ordVec = NULL,
                       gapVecR = NULL,
                       gapVecC = NULL,
                       cc = NULL,
                       cr = NULL,
                       condCol = NULL,
                       ...) {

  # ===========================================================================
  # BACKWARD COMPATIBILITY: Map old arguments to new ones
  # ===========================================================================

  # Map colVecIdent to annotation_colors
  if (!is.null(colVecIdent)) {
    message("Note: 'colVecIdent' is deprecated. Please use 'annotation_colors' instead.")
    if (is.null(annotation_colors)) {
      annotation_colors <- colVecIdent
    }
  }

  # Map colVecCond to condition_colors
  if (!is.null(colVecCond)) {
    message("Note: 'colVecCond' is deprecated. Please use 'condition_colors' instead.")
    if (is.null(condition_colors)) {
      condition_colors <- colVecCond
    }
  }

  # Map condCol to condition_by (if it's TRUE, try to detect condition column)
  if (!is.null(condCol) && condCol == TRUE && is.null(condition_by)) {
    # Try to find a 'cond' or 'condition' column in metadata
    if ("cond" %in% colnames(seurat@meta.data)) {
      condition_by <- "cond"
      message("Note: 'condCol=TRUE' detected. Using 'cond' column for condition annotation.")
    } else if ("condition" %in% colnames(seurat@meta.data)) {
      condition_by <- "condition"
      message("Note: 'condCol=TRUE' detected. Using 'condition' column for condition annotation.")
    }
  }

  # Map gapVecR to gaps_row
  if (!is.null(gapVecR)) {
    message("Note: 'gapVecR' is deprecated. Please use 'gaps_row' instead.")
    if (is.null(gaps_row)) {
      gaps_row <- gapVecR
    }
  }

  # Map gapVecC to gaps_col
  if (!is.null(gapVecC)) {
    message("Note: 'gapVecC' is deprecated. Please use 'gaps_col' instead.")
    if (is.null(gaps_col)) {
      gaps_col <- gapVecC
    }
  }

  # Map cc to cluster_cols
  if (!is.null(cc)) {
    message("Note: 'cc' is deprecated. Please use 'cluster_cols' instead.")
    cluster_cols <- cc
  }

  # Map cr to cluster_rows
  if (!is.null(cr)) {
    message("Note: 'cr' is deprecated. Please use 'cluster_rows' instead.")
    cluster_rows <- cr
  }

  # ===========================================================================
  # MAIN FUNCTION LOGIC
  # ===========================================================================

  # Handle gene input - can be vector of gene names or data.frame with gene column
  # If no genes provided, use top variable features
  if (is.null(selGenes)) {
    # Get top n highly variable features as default
    if ("SCT" %in% names(seurat@assays)) {
      # If SCT assay exists, use its variable features
      variable_features <- VariableFeatures(seurat, assay = "SCT")
    } else {
      # Otherwise use RNA assay variable features
      variable_features <- VariableFeatures(seurat, assay = "RNA")
    }

    if (length(variable_features) == 0) {
      stop("No variable features found. Please run FindVariableFeatures() first or provide selGenes manually.")
    }

    gene_list <- head(variable_features, n_variable_genes)
    message("No genes provided. Using top ", n_variable_genes, " variable features: ",
            paste(head(gene_list, 5), collapse = ", "), "...")

  } else if (is.data.frame(selGenes)) {
    if ("gene" %in% colnames(selGenes)) {
      gene_list <- selGenes$gene
    } else if ("geneID" %in% colnames(selGenes)) {
      gene_list <- selGenes$geneID
    } else {
      gene_list <- selGenes[, 1]  # Use first column if no standard column names
    }
  } else {
    gene_list <- selGenes
  }

  # Set grouping variable - default to active identity
  if (is.null(group_by)) {
    group_by <- "ident"
    clusterAssigned <- data.frame(
      ident = Idents(seurat),
      cell = names(Idents(seurat))
    )
  } else {
    clusterAssigned <- data.frame(
      ident = seurat@meta.data[[group_by]],
      cell = rownames(seurat@meta.data)
    )
  }

  # Add condition if specified
  if (!is.null(condition_by)) {
    if (!condition_by %in% colnames(seurat@meta.data)) {
      stop("condition_by '", condition_by, "' not found in seurat metadata")
    }
    clusterAssigned$condition <- seurat@meta.data[[condition_by]]
    # Create combined identity: celltype_condition
    clusterAssigned$combined_ident <- paste0(clusterAssigned$ident, "_", clusterAssigned$condition)
  }

  # Get assay data
  seuratDat <- GetAssayData(seurat, assay = "RNA", slot = "data")

  # Find genes in the Seurat object
  # Handle ENSEMBL IDs if present (assumes format: ENSEMBL.SYMBOL)
  genes <- data.frame(gene = rownames(seurat)) %>%
    mutate(geneID = gsub("^.*\\.", "", gene))

  # Match genes more flexibly
  if (all(grepl("^ENSG", gene_list))) {
    # If input genes are ENSEMBL IDs, match directly
    matched_genes <- genes[genes$gene %in% gene_list, ]
  } else {
    # If input genes are symbols, match against cleaned gene symbols
    matched_genes <- genes[genes$geneID %in% gene_list | genes$gene %in% gene_list, ]
  }

  if (nrow(matched_genes) == 0) {
    stop("No matching genes found in the Seurat object. Check gene names.")
  }

  # Create expression matrix averaged by identity (and condition if specified)
  logNormExpres <- as.data.frame(t(as.matrix(
    seuratDat[matched_genes$gene, ]
  )))

  logNormExpres <- logNormExpres %>%
    mutate(cell = rownames(.)) %>%
    left_join(clusterAssigned, by = "cell") %>%
    dplyr::select(-cell)

  # Group by combined identity if condition is specified
  if (!is.null(condition_by)) {
    logNormExpres <- logNormExpres %>%
      group_by(combined_ident) %>%
      summarise_at(vars(-ident, -condition), mean, .groups = 'drop')
  } else {
    logNormExpres <- logNormExpres %>%
      group_by(ident) %>%
      summarise_all(mean, .groups = 'drop')
  }

  # Convert to matrix format for heatmap
  if (!is.null(condition_by)) {
    logNormExpresMa <- logNormExpres %>%
      dplyr::select(-combined_ident) %>%
      as.matrix()
    rownames(logNormExpresMa) <- logNormExpres$combined_ident
  } else {
    logNormExpresMa <- logNormExpres %>%
      dplyr::select(-ident) %>%
      as.matrix()
    rownames(logNormExpresMa) <- logNormExpres$ident
  }

  logNormExpresMa <- t(logNormExpresMa)

  # Clean up gene names (remove ENSEMBL prefix if present)
  rownames(logNormExpresMa) <- gsub("^.*?\\.", "", rownames(logNormExpresMa))

  # Remove genes with zero variance across all groups
  zero_var_genes <- apply(logNormExpresMa, 1, sd) == 0
  if (any(zero_var_genes)) {
    logNormExpresMa <- logNormExpresMa[!zero_var_genes, , drop = FALSE]
    warning(paste("Removed", sum(zero_var_genes), "genes with zero variance across groups"))
  }

  # Apply ordVec if provided (backward compatibility for column ordering)
  if (!is.null(ordVec)) {
    message("Note: 'ordVec' is deprecated but will be applied to reorder columns.")
    # Only reorder columns that exist in the matrix
    ordVec_valid <- ordVec[ordVec %in% colnames(logNormExpresMa)]
    if (length(ordVec_valid) > 0) {
      # Add any columns not in ordVec to the end
      remaining_cols <- setdiff(colnames(logNormExpresMa), ordVec_valid)
      final_order <- c(ordVec_valid, remaining_cols)
      logNormExpresMa <- logNormExpresMa[, final_order, drop = FALSE]
    }
  }

  # Order genes to create "staircase" pattern - genes ordered by which cluster has highest expression
  if (!cluster_rows) {
    # Find which cluster has max expression for each gene
    max_clusters <- apply(logNormExpresMa, 1, function(x) colnames(logNormExpresMa)[which.max(x)])

    # Order genes by their max-expressing cluster, then by max expression level within cluster
    gene_order <- names(sort(factor(max_clusters, levels = colnames(logNormExpresMa))))

    # Within each cluster group, order genes by decreasing max expression
    ordered_genes <- c()
    for (cluster in colnames(logNormExpresMa)) {
      cluster_genes <- names(max_clusters[max_clusters == cluster])
      if (length(cluster_genes) > 0) {
        # Order by decreasing expression in that cluster
        cluster_gene_order <- cluster_genes[order(logNormExpresMa[cluster_genes, cluster], decreasing = TRUE)]
        ordered_genes <- c(ordered_genes, cluster_gene_order)
      }
    }

    # Reorder the matrix
    logNormExpresMa <- logNormExpresMa[ordered_genes, , drop = FALSE]
  }

  # Set default color palette for heatmap
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(50)
  }

  # ===========================================================================
  # CREATE ANNOTATION DATA FRAME
  # ===========================================================================

  if (!is.null(condition_by)) {
    # Parse combined identities back into celltype and condition
    col_names <- colnames(logNormExpresMa)

    # Extract celltype and condition from combined names
    celltype_vec <- gsub("_[^_]*$", "", col_names)  # Remove last underscore part
    condition_vec <- gsub("^.*_", "", col_names)     # Get last underscore part

    annotation_col <- data.frame(
      Celltype = celltype_vec,
      Condition = condition_vec,
      row.names = col_names
    )

    # Get unique values for color assignment
    celltypes_present <- unique(celltype_vec)
    conditions_present <- unique(condition_vec)
    n_celltypes <- length(celltypes_present)
    n_conditions <- length(conditions_present)

  } else {
    # Simple annotation without condition
    annotation_col <- data.frame(
      Group = colnames(logNormExpresMa),
      row.names = colnames(logNormExpresMa)
    )

    groups_present <- unique(annotation_col$Group)
    n_groups <- length(groups_present)
  }

  # ===========================================================================
  # SETUP ANNOTATION COLORS
  # ===========================================================================

  # Define color palettes
  disease_colors <- c("#dfc27d", "#BE3144", "#202547", "#355C7D", "#779d8d")
  fibroblast_colors <- c("#D53E4F", "#f4a582", "#ff7b7b", "#8e0b00", "#FEE08B",
                         "#42090D", "#FF7B00", "#FFF4DF")
  nice_colors <- c("#67001f", "#D53E4F", "#f4a582", "#FEE08B", "#003c30", "#01665e",
                   "#66C2A5", "#3288BD", "#BEAED4", "#c7eae5", "#355C7D", "#202547",
                   "#B45B5C", "#8c510a")
  extended_colors <- c("#fde0dd", "#fa9fb5", "#d95f0e", "#dd1c77", "#D53E4F",
                       "#f4a582", "#FEE08B", "#f03b20", "#ffffcc", "#43a2ca",
                       "#1c9099", "#355C7D", "#3288BD", "#BEAED4", "#756bb1",
                       "#c7eae5")
  large_palette <- c(
    "#fde0dd", "#fa9fb5", "#f768a1", "#dd1c77", "#980043",
    "#f4a582", "#fdae61", "#f46d43", "#d73027", "#a50026",
    "#fee08b", "#ffffbf", "#e6f598", "#99d594", "#66c2a5",
    "#43a2ca", "#1c9099", "#016c59", "#3288bd", "#5e4fa2",
    "#beaed4", "#9e9ac8", "#756bb1", "#542788", "#3f007d",
    "#c7eae5", "#a6bddb", "#74a9cf", "#3690c0", "#045a8d"
  )

  if (!is.null(condition_by)) {
    # Handle colors for dual annotation (celltype + condition)

    # === CELLTYPE COLORS ===
    if (!is.null(annotation_colors)) {
      if (is.vector(annotation_colors) && !is.list(annotation_colors)) {
        if (!is.null(names(annotation_colors))) {
          celltype_colors <- annotation_colors[names(annotation_colors) %in% celltypes_present]
          if (length(celltype_colors) < n_celltypes) {
            missing <- setdiff(celltypes_present, names(celltype_colors))
            default_cols <- rainbow(length(missing))
            names(default_cols) <- missing
            celltype_colors <- c(celltype_colors, default_cols)
          }
        } else {
          if (length(annotation_colors) < n_celltypes) {
            annotation_colors <- c(annotation_colors, rainbow(n_celltypes - length(annotation_colors)))
          }
          celltype_colors <- annotation_colors[1:n_celltypes]
          names(celltype_colors) <- celltypes_present
        }
      } else if (is.list(annotation_colors)) {
        if ("Celltype" %in% names(annotation_colors)) {
          celltype_colors <- annotation_colors$Celltype
        } else if ("Group" %in% names(annotation_colors)) {
          celltype_colors <- annotation_colors$Group
        } else {
          celltype_colors <- annotation_colors[[1]]
        }
        if (!is.null(names(celltype_colors))) {
          celltype_colors <- celltype_colors[names(celltype_colors) %in% celltypes_present]
        }
      }
    } else {
      # Auto-select palette for celltypes
      if (n_celltypes <= 5 && any(grepl("healthy|explant|visit", celltypes_present, ignore.case = TRUE))) {
        celltype_colors <- disease_colors[1:n_celltypes]
      } else if (n_celltypes <= 8 && any(grepl("Fb|Periv|VSMC", celltypes_present, ignore.case = TRUE))) {
        celltype_colors <- fibroblast_colors[1:n_celltypes]
      } else if (n_celltypes <= 14) {
        celltype_colors <- nice_colors[1:n_celltypes]
      } else if (n_celltypes <= 16) {
        celltype_colors <- extended_colors[1:n_celltypes]
      } else if (n_celltypes <= 30) {
        celltype_colors <- large_palette[1:n_celltypes]
      } else {
        celltype_colors <- c(large_palette, rainbow(n_celltypes - length(large_palette)))
      }
      names(celltype_colors) <- celltypes_present
    }

    # === CONDITION COLORS ===
    if (!is.null(condition_colors)) {
      if (!is.null(names(condition_colors))) {
        cond_colors <- condition_colors[names(condition_colors) %in% conditions_present]
        if (length(cond_colors) < n_conditions) {
          missing <- setdiff(conditions_present, names(cond_colors))
          default_cols <- c("forestgreen", "firebrick", "steelblue", "orange")[1:length(missing)]
          names(default_cols) <- missing
          cond_colors <- c(cond_colors, default_cols)
        }
      } else {
        cond_colors <- condition_colors[1:n_conditions]
        names(cond_colors) <- conditions_present
      }
    } else {
      # Default condition colors
      default_cond_colors <- c("forestgreen", "firebrick", "steelblue", "orange", "purple")
      cond_colors <- default_cond_colors[1:n_conditions]
      names(cond_colors) <- conditions_present
    }

    ann_colors_final <- list(
      Celltype = celltype_colors,
      Condition = cond_colors
    )

  } else {
    # Handle colors for single annotation (group only)

    if (is.null(annotation_colors)) {
      # Auto-select palette
      if (n_groups <= 5 && any(grepl("healthy|explant|visit", groups_present, ignore.case = TRUE))) {
        group_colors <- disease_colors[1:n_groups]
      } else if (n_groups <= 8 && any(grepl("Fb|Periv|VSMC", groups_present, ignore.case = TRUE))) {
        group_colors <- fibroblast_colors[1:n_groups]
      } else if (n_groups <= 14) {
        group_colors <- nice_colors[1:n_groups]
      } else if (n_groups <= 16) {
        group_colors <- extended_colors[1:n_groups]
      } else if (n_groups <= 30) {
        group_colors <- large_palette[1:n_groups]
      } else {
        group_colors <- c(large_palette, rainbow(n_groups - length(large_palette)))
      }
      names(group_colors) <- groups_present
      ann_colors_final <- list(Group = group_colors)

    } else {
      # Process user-provided colors
      if (is.vector(annotation_colors) && !is.list(annotation_colors)) {
        if (!is.null(names(annotation_colors))) {
          group_colors <- annotation_colors[names(annotation_colors) %in% groups_present]
          if (length(group_colors) < n_groups) {
            missing <- setdiff(groups_present, names(group_colors))
            default_cols <- rainbow(length(missing))
            names(default_cols) <- missing
            group_colors <- c(group_colors, default_cols)
          }
        } else {
          if (length(annotation_colors) < n_groups) {
            annotation_colors <- c(annotation_colors, rainbow(n_groups - length(annotation_colors)))
          }
          group_colors <- annotation_colors[1:n_groups]
          names(group_colors) <- groups_present
        }
        ann_colors_final <- list(Group = group_colors)

      } else if (is.list(annotation_colors)) {
        if ("Group" %in% names(annotation_colors)) {
          group_colors <- annotation_colors$Group
        } else {
          group_colors <- annotation_colors[[1]]
        }
        if (!is.null(names(group_colors))) {
          group_colors <- group_colors[names(group_colors) %in% groups_present]
        }
        ann_colors_final <- list(Group = group_colors)
      }
    }
  }

  # Generate and return heatmap
  p <- pheatmap(logNormExpresMa,
                scale = scale_method,
                cluster_rows = cluster_rows,
                cluster_cols = cluster_cols,
                color = color_palette,
                annotation_col = annotation_col,
                annotation_colors = ann_colors_final,
                cellwidth = cellwidth,
                cellheight = cellheight,
                show_rownames = show_rownames,
                show_colnames = show_colnames,
                gaps_row = gaps_row,
                gaps_col = gaps_col,
                ...)

  return(p)
}
