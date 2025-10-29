#' Average Expression Heatmap for Seurat Objects with Significance Testing
#'
#' Generates a heatmap of averaged expression values for selected genes across cell groups or identities in a Seurat object.
#' If no genes are provided, the top variable features will be used. Optionally performs statistical testing and displays
#' significance stars on the heatmap.
#'
#' @param seurat A Seurat object.
#' @param selGenes Character vector of genes to plot, or a data.frame with a column named `gene` or `geneID`. Defaults to `NULL`.
#' @param group_by Metadata column in the Seurat object to group cells by. Defaults to `NULL` (uses active identity).
#' @param condition_by Optional metadata column for a second grouping variable (e.g., treatment condition). Default is `NULL`.
#' @param scale_method Scaling method for the heatmap: `"row"`, `"column"`, or `"none"`. Default is `"row"`.
#' @param cluster_rows Logical, whether to cluster rows. Default is `FALSE`.
#' @param cluster_cols Logical, whether to cluster columns. Default is `FALSE`.
#' @param cluster_order Optional character vector to manually set the order of columns (clusters). Overrides `cluster_cols = FALSE`.
#' @param gene_order Optional character vector to manually set the order of rows (genes). Overrides `cluster_rows = FALSE`.
#' @param show_rownames Logical, whether to display row names (genes). Default is `TRUE`.
#' @param show_colnames Logical, whether to display column names (groups). Default is `TRUE`.
#' @param cellwidth Numeric, width of each cell in the heatmap. Default is `15`.
#' @param cellheight Numeric, height of each cell in the heatmap. Default is `10`.
#' @param color_palette Color vector for the heatmap. Default is a blue-white-red gradient.
#' @param annotation_colors Named list of colors for annotations. Can contain elements named 'Group' and 'Condition'.
#'    Alternatively, provide a named vector which will be used for the Group annotation.
#' @param condition_colors Named vector of colors for condition annotation. Only used if `condition_by` is specified.
#' @param gaps_row Optional vector specifying rows after which to insert gaps. Default is `NULL`.
#' @param gaps_col Optional vector specifying columns after which to insert gaps. Default is `NULL`.
#' @param n_variable_genes Number of top variable genes to use if `selGenes` is `NULL`. Default is `20`.
#' @param show_significance Logical, whether to perform statistical testing and show significance stars. Default is `FALSE`.
#' @param significance_direction Direction for showing significance: "higher" (default), "lower", or "both".
#' @param significance_test Statistical test to use: "wilcox" (default) or "t.test".
#' @param pval_cutoffs Named numeric vector of p-value cutoffs. Default is `c("***" = 0.001, "**" = 0.01, "*" = 0.05)`.
#' @param star_size Numeric, font size for significance stars. Default is `8`.
#' @param p_adjust_method Method for p-value adjustment. Default is `"BH"` (Benjamini-Hochberg/FDR).
#'    Use `"none"` for unadjusted p-values. See `?p.adjust` for other options.
#' @param ... Additional arguments passed to [pheatmap::pheatmap()].
#'
#' @return A `pheatmap` object.
#'
#' @author Mechthild Lütge, Roman Stadler
#'
#' @details
#' - Handles gene names as symbols or ENSEMBL IDs.
#' - Automatically removes genes with zero variance across groups.
#' - If `cluster_rows = FALSE` and `gene_order = NULL` (default),
#'    orders genes by cluster of maximum expression to create a "staircase" pattern.
#' - If `gene_order` is provided, genes are kept in the user-provided order.
#' - If `cluster_order` is provided, clusters are kept in the user-provided order.
#'    (When `condition_by` is used, this groups columns based on the cluster prefix).
#' - Automatically selects color palettes based on group names and number of groups.
#' - When `condition_by` is specified, creates combined groups (e.g., "Cluster0_ConditionA").
#' - When `show_significance = TRUE`:
#'   - If `condition_by` is `NULL`: performs test comparing each cluster vs. all other clusters.
#'   - If `condition_by` is set: performs test comparing each group (e.g., "Cluster0_CondA")
#'     vs. all *other groups within that same base cluster* (e.g., "Cluster0_CondB", "Cluster0_CondC").
#'   - P-values are adjusted for multiple testing using the Benjamini-Hochberg method (FDR) by default.
#'
#' @importFrom Seurat VariableFeatures Idents GetAssayData
#' @importFrom dplyr mutate left_join select group_by summarise
#' @importFrom pheatmap pheatmap
#' @importFrom stats sd wilcox.test t.test p.adjust
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' seurat <- YourSeuratObject
#'
#' # Simple heatmap with significance stars (higher than average)
#' avgHeatmap(seurat,
#'            selGenes = c("GeneA", "GeneB", "GeneC"),
#'            group_by = "celltype",
#'            show_significance = TRUE)
#'
#' # Show significance for both higher and lower expression
#' avgHeatmap(seurat,
#'            selGenes = c("GeneA", "GeneB"),
#'            group_by = "celltype",
#'            show_significance = TRUE,
#'            significance_direction = "both")
#'
#' # Custom p-value cutoffs
#' avgHeatmap(seurat,
#'            selGenes = c("GeneA", "GeneB"),
#'            group_by = "celltype",
#'            show_significance = TRUE,
#'            pval_cutoffs = c("***" = 0.0001, "**" = 0.001, "*" = 0.01))
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
                       cluster_order = NULL,
                       gene_order = NULL,
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
                       show_significance = FALSE,
                       significance_direction = "higher",
                       significance_test = "wilcox",
                       pval_cutoffs = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
                       star_size = 8,
                       p_adjust_method = "BH",
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

  if (!is.null(colVecIdent)) {
    message("Note: 'colVecIdent' is deprecated. Please use 'annotation_colors' instead.")
    if (is.null(annotation_colors)) {
      annotation_colors <- colVecIdent
    }
  }

  if (!is.null(colVecCond)) {
    message("Note: 'colVecCond' is deprecated. Please use 'condition_colors' instead.")
    if (is.null(condition_colors)) {
      condition_colors <- colVecCond
    }
  }

  if (!is.null(ordVec)) {
    message("Note: 'ordVec' is deprecated. Please use 'cluster_order' instead.")
    if (is.null(cluster_order)) {
      cluster_order <- ordVec
    }
  }

  if (!is.null(condCol) && condCol == TRUE && is.null(condition_by)) {
    if ("cond" %in% colnames(seurat@meta.data)) {
      condition_by <- "cond"
      message("Note: 'condCol=TRUE' detected. Using 'cond' column for condition annotation.")
    } else if ("condition" %in% colnames(seurat@meta.data)) {
      condition_by <- "condition"
      message("Note: 'condCol=TRUE' detected. Using 'condition' column for condition annotation.")
    }
  }

  if (!is.null(gapVecR)) {
    message("Note: 'gapVecR' is deprecated. Please use 'gaps_row' instead.")
    if (is.null(gaps_row)) {
      gaps_row <- gapVecR
    }
  }

  if (!is.null(gapVecC)) {
    message("Note: 'gapVecC' is deprecated. Please use 'gaps_col' instead.")
    if (is.null(gaps_col)) {
      gaps_col <- gapVecC
    }
  }

  if (!is.null(cc)) {
    message("Note: 'cc' is deprecated. Please use 'cluster_cols' instead.")
    cluster_cols <- cc
  }

  if (!is.null(cr)) {
    message("Note: 'cr' is deprecated. Please use 'cluster_rows' instead.")
    cluster_rows <- cr
  }

  # ===========================================================================
  # MAIN FUNCTION LOGIC
  # ===========================================================================

  # Handle gene input
  if (is.null(selGenes)) {
    if ("SCT" %in% names(seurat@assays)) {
      variable_features <- VariableFeatures(seurat, assay = "SCT")
    } else {
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
      gene_list <- selGenes[, 1]
    }
  } else {
    gene_list <- selGenes
  }

  # Set grouping variable
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
    clusterAssigned$combined_ident <- paste0(clusterAssigned$ident, "_", clusterAssigned$condition)
  }

  # Get assay data
  seuratDat <- GetAssayData(seurat, assay = "RNA", slot = "data")

  # Find genes in the Seurat object
  genes <- data.frame(gene = rownames(seurat)) %>%
    mutate(geneID = gsub("^.*\\.", "", gene))

  # Match genes
  if (all(grepl("^ENSG", gene_list))) {
    matched_genes <- genes[genes$gene %in% gene_list, ]
  } else {
    matched_genes <- genes[genes$geneID %in% gene_list | genes$gene %in% gene_list, ]
  }

  if (nrow(matched_genes) == 0) {
    stop("No matching genes found in the Seurat object. Check gene names.")
  }

  # Create expression matrix averaged by identity
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

  # Clean up gene names
  rownames(logNormExpresMa) <- gsub("^.*?\\.", "", rownames(logNormExpresMa))

  # Remove genes with zero variance
  zero_var_genes <- apply(logNormExpresMa, 1, sd) == 0
  if (any(zero_var_genes)) {
    logNormExpresMa <- logNormExpresMa[!zero_var_genes, , drop = FALSE]
    warning(paste("Removed", sum(zero_var_genes), "genes with zero variance across groups"))
  }

  # ===========================================================================
  # SIGNIFICANCE TESTING (MODIFIED)
  # ===========================================================================

  significance_matrix <- NULL

  if (show_significance) {
    message("Performing significance testing...")

    # Create a matrix to store significance symbols
    significance_matrix <- matrix("", nrow = nrow(logNormExpresMa), ncol = ncol(logNormExpresMa))
    rownames(significance_matrix) <- rownames(logNormExpresMa)
    colnames(significance_matrix) <- colnames(logNormExpresMa)

    # Get the raw expression data for statistical testing
    raw_expression <- as.data.frame(t(as.matrix(seuratDat[matched_genes$gene, ])))
    raw_expression <- raw_expression %>%
      mutate(cell = rownames(.)) %>%
      left_join(clusterAssigned, by = "cell")

    # Clean column names to match
    gene_cols <- setdiff(colnames(raw_expression), c("cell", "ident", "condition", "combined_ident"))
    clean_gene_names <- gsub("^.*?\\.", "", gene_cols)

    # Determine which grouping variable to use
    group_var <- if (!is.null(condition_by)) "combined_ident" else "ident"

    # Store p-values for adjustment
    pval_list <- list()
    test_info <- list()

    # For each cluster and each gene, perform test
    test_counter <- 0
    for (cluster in colnames(logNormExpresMa)) {
      for (i in seq_along(gene_cols)) {
        gene_col <- gene_cols[i]
        gene_name <- clean_gene_names[i]

        # Skip if gene was filtered out
        if (!gene_name %in% rownames(logNormExpresMa)) next

        # Get expression values
        cluster_cells <- raw_expression[[group_var]] == cluster
        in_cluster <- raw_expression[cluster_cells, gene_col]

        # --- MODIFICATION START ---
        # Check if we are in the 'condition_by' mode
        if (!is.null(condition_by)) {
          # New logic: Compare only within the same base cluster (e.g., "Fb1")

          # 1. Get base cluster (e.g., "Fb1" from "Fb1_donorheart")
          base_cluster <- gsub("_[^_]*$", "", cluster)

          # 2. Find all cluster names from the heatmap
          all_heatmap_clusters <- colnames(logNormExpresMa)

          # 3. Find all "siblings" (e.g., "Fb1_donorheart", "Fb1_explant", "Fb1_control")
          sibling_clusters <- grep(paste0("^", base_cluster, "_"), all_heatmap_clusters, value = TRUE)

          # 4. Find the *other* siblings to compare against
          other_sibling_clusters <- setdiff(sibling_clusters, cluster)

          # 5. Get the cells belonging to these other siblings
          out_cells <- raw_expression[[group_var]] %in% other_sibling_clusters
          out_cluster <- raw_expression[out_cells, gene_col]

        } else {
          # Original logic: compare cluster vs. all other clusters
          out_cluster <- raw_expression[!cluster_cells, gene_col]
        }
        # --- MODIFICATION END ---


        # Skip if not enough data
        if (length(in_cluster) < 3 || length(out_cluster) < 3) next

        # Perform statistical test
        tryCatch({
          if (significance_test == "wilcox") {
            test_result <- wilcox.test(in_cluster, out_cluster, alternative = "two.sided")
          } else if (significance_test == "t.test") {
            test_result <- t.test(in_cluster, out_cluster, alternative = "two.sided")
          } else {
            stop("significance_test must be 'wilcox' or 't.test'")
          }

          test_counter <- test_counter + 1
          pval <- test_result$p.value

          # Use median for comparison, as requested for wilcox
          # (mean is also fine for direction, but median is more robust)
          median_in <- median(in_cluster)
          median_out <- median(out_cluster)
          is_higher <- median_in > median_out

          # Store for adjustment
          pval_list[[test_counter]] <- pval
          test_info[[test_counter]] <- list(
            gene = gene_name,
            cluster = cluster,
            is_higher = is_higher
          )

        }, error = function(e) {
          # Skip genes that cause errors
        })
      }
    }

    # Adjust p-values for multiple testing
    if (length(pval_list) > 0) {
      if (p_adjust_method == "none") {
        adjusted_pvals <- unlist(pval_list)
        message("Using unadjusted p-values (no multiple testing correction)")
      } else {
        adjusted_pvals <- p.adjust(unlist(pval_list), method = p_adjust_method)
        message(sprintf("Adjusted %d p-values using %s method", length(pval_list), p_adjust_method))
      }

      # Assign significance stars based on adjusted p-values
      for (j in seq_along(adjusted_pvals)) {
        info <- test_info[[j]]
        pval <- adjusted_pvals[j]

        # Determine if we should show based on direction preference
        should_show <- FALSE
        if (significance_direction == "higher" && info$is_higher) {
          should_show <- TRUE
        } else if (significance_direction == "lower" && !info$is_higher) {
          should_show <- TRUE
        } else if (significance_direction == "both") {
          should_show <- TRUE
        }

        if (should_show) {
          # Sort cutoffs by p-value (most stringent first)
          sorted_cutoffs <- sort(pval_cutoffs)

          for (k in seq_along(sorted_cutoffs)) {
            if (pval < sorted_cutoffs[k]) {
              significance_matrix[info$gene, info$cluster] <- names(sorted_cutoffs)[k]
              break
            }
          }
        }
      }
    }
  }

  # ===========================================================================
  # Column and Row Ordering
  # ===========================================================================

  # Apply cluster_order if provided
  if (!is.null(cluster_order)) {
    if (!is.null(condition_by)) {
      message("Applying manual column order (cluster_order) with condition_by grouping.")
      final_order <- c()
      current_cols <- colnames(logNormExpresMa)

      for (base_cluster in cluster_order) {
        cols_to_add <- sort(grep(paste0("^", base_cluster, "_"), current_cols, value = TRUE))
        if (length(cols_to_add) > 0) {
          final_order <- c(final_order, cols_to_add)
          current_cols <- setdiff(current_cols, cols_to_add)
        }
      }
      final_order <- c(final_order, current_cols)

      if (length(final_order) == length(colnames(logNormExpresMa))) {
        logNormExpresMa <- logNormExpresMa[, final_order, drop = FALSE]
        if (!is.null(significance_matrix)) {
          significance_matrix <- significance_matrix[, final_order, drop = FALSE]
        }
      }
    } else {
      message("Applying manual column order (cluster_order).")
      order_valid <- cluster_order[cluster_order %in% colnames(logNormExpresMa)]
      if (length(order_valid) > 0) {
        remaining_cols <- setdiff(colnames(logNormExpresMa), order_valid)
        final_order <- c(order_valid, remaining_cols)
        logNormExpresMa <- logNormExpresMa[, final_order, drop = FALSE]
        if (!is.null(significance_matrix)) {
          significance_matrix <- significance_matrix[, final_order, drop = FALSE]
        }
      }
    }
  }

  # Order genes
  if (!cluster_rows) {
    if (!is.null(gene_order)) {
      message("Applying manual gene order (gene_order).")
      cleaned_gene_list <- gsub("^.*?\\.", "", gene_order)
      ordered_genes_from_list <- cleaned_gene_list[cleaned_gene_list %in% rownames(logNormExpresMa)]
      remaining_genes <- setdiff(rownames(logNormExpresMa), ordered_genes_from_list)
      final_gene_order <- c(ordered_genes_from_list, remaining_genes)

      if (length(final_gene_order) > 0) {
        logNormExpresMa <- logNormExpresMa[final_gene_order, , drop = FALSE]
        if (!is.null(significance_matrix)) {
          significance_matrix <- significance_matrix[final_gene_order, , drop = FALSE]
        }
      }
    } else {
      message("Ordering genes by cluster of max expression (default).")
      max_clusters <- apply(logNormExpresMa, 1, function(x) colnames(logNormExpresMa)[which.max(x)])
      current_col_order <- colnames(logNormExpresMa)
      ordered_genes <- c()

      for (cluster in current_col_order) {
        cluster_genes <- names(max_clusters[max_clusters == cluster])
        if (length(cluster_genes) > 0) {
          cluster_gene_order <- cluster_genes[order(logNormExpresMa[cluster_genes, cluster], decreasing = TRUE)]
          ordered_genes <- c(ordered_genes, cluster_gene_order)
        }
      }

      if (length(ordered_genes) > 0) {
        logNormExpresMa <- logNormExpresMa[ordered_genes, , drop = FALSE]
        if (!is.null(significance_matrix)) {
          significance_matrix <- significance_matrix[ordered_genes, , drop = FALSE]
        }
      }
    }
  }

  # ===========================================================================
  # Color Setup
  # ===========================================================================

  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(50)
  }

  # Create annotation
  if (!is.null(condition_by)) {
    col_names <- colnames(logNormExpresMa)
    celltype_vec <- gsub("_[^_]*$", "", col_names)
    condition_vec <- gsub("^.*_", "", col_names)

    annotation_col <- data.frame(
      Celltype = celltype_vec,
      Condition = condition_vec,
      row.names = col_names
    )

    celltypes_present <- unique(celltype_vec)
    conditions_present <- unique(condition_vec)
    n_celltypes <- length(celltypes_present)
    n_conditions <- length(conditions_present)
  } else {
    annotation_col <- data.frame(
      Group = colnames(logNormExpresMa),
      row.names = colnames(logNormExpresMa)
    )
    groups_present <- unique(annotation_col$Group)
    n_groups <- length(groups_present)
  }

  # Color palettes
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
    # Handle dual annotation colors
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

    if (!is.null(condition_colors)) {
      if (!is.null(names(condition_colors))) {
        cond_colors <- condition_colors[names(condition_colors) %in% conditions_present]
        if (length(cond_colors) < n_conditions) {
          missing <- setdiff(conditions_present, names(cond_colors))
          default_cols <- c("steelblue", "orange", "purple", "forestgreen", "firebrick",
                            "goldenrod", "turquoise", "violet", "darkolivegreen", "coral",
                            "slateblue", "tomato", "mediumorchid", "darkgoldenrod", "cadetblue")[1:length(missing)]
          names(default_cols) <- missing
          cond_colors <- c(cond_colors, default_cols)
        }
      } else {
        cond_colors <- condition_colors[1:n_conditions]
        names(cond_colors) <- conditions_present
      }
    } else {
      default_cond_colors <- c("steelblue", "orange", "purple", "forestgreen", "firebrick",
                               "goldenrod", "turquoise", "violet", "darkolivegreen", "coral",
                               "slateblue", "tomato", "mediumorchid", "darkgoldenrod", "cadetblue")
      cond_colors <- default_cond_colors[1:n_conditions]
      names(cond_colors) <- conditions_present
    }

    ann_colors_final <- list(
      Celltype = celltype_colors,
      Condition = cond_colors
    )
  } else {
    # Single annotation colors
    if (is.null(annotation_colors)) {
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

  # ===========================================================================
  # Generate Heatmap
  # ===========================================================================

  # Prepare parameters for pheatmap
  pheatmap_params <- list(
    mat = logNormExpresMa,
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
    gaps_col = gaps_col
  )

  # Add significance stars if requested
  if (show_significance && !is.null(significance_matrix)) {
    pheatmap_params$display_numbers <- significance_matrix
    pheatmap_params$number_color <- "black"
    pheatmap_params$fontsize_number <- star_size
  }

  # Add any additional arguments
  pheatmap_params <- c(pheatmap_params, list(...))

  # Generate heatmap
  p <- do.call(pheatmap::pheatmap, pheatmap_params)

  return(p)
}
