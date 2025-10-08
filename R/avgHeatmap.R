#' Average Expression Heatmap for Seurat Objects
#'
#' Generates a heatmap of averaged expression values for selected genes across cell groups or identities in a Seurat object.
#' If no genes are provided, the top variable features will be used.
#'
#' @param seurat A Seurat object.
#' @param selGenes Character vector of genes to plot, or a data.frame with a column named `gene` or `geneID`. Defaults to `NULL`.
#' @param group_by Metadata column in the Seurat object to group cells by. Defaults to `NULL` (uses active identity).
#' @param scale_method Scaling method for the heatmap: `"row"`, `"column"`, or `"none"`. Default is `"row"`.
#' @param cluster_rows Logical, whether to cluster rows. Default is `FALSE`.
#' @param cluster_cols Logical, whether to cluster columns. Default is `FALSE`.
#' @param show_rownames Logical, whether to display row names (genes). Default is `TRUE`.
#' @param show_colnames Logical, whether to display column names (groups). Default is `TRUE`.
#' @param cellwidth Numeric, width of each cell in the heatmap. Default is `15`.
#' @param cellheight Numeric, height of each cell in the heatmap. Default is `10`.
#' @param color_palette Color vector for the heatmap. Default is a blue-white-red gradient.
#' @param annotation_colors Named list or named vector of colors for annotations.
#'   If a named vector is provided, it will be automatically wrapped as `list(Group = your_vector)`.
#'   The vector should be named with group levels, or have length matching the number of groups.
#' @param gaps_row Optional vector specifying rows after which to insert gaps. Default is `NULL`.
#' @param gaps_col Optional vector specifying columns after which to insert gaps. Default is `NULL`.
#' @param n_variable_genes Number of top variable genes to use if `selGenes` is `NULL`. Default is `20`.
#' @param ... Additional arguments passed to [pheatmap::pheatmap()].
#'
#' @return A `pheatmap` object.
#'
#' @details
#' - Handles gene names as symbols or ENSEMBL IDs.
#' - Automatically removes genes with zero variance across groups.
#' - Orders genes by cluster of maximum expression if `cluster_rows = FALSE`.
#' - Automatically selects color palettes based on group names and number of groups.
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
#' # Using default colors
#' avgHeatmap(seurat, selGenes = c("GeneA", "GeneB"), group_by = "celltype")
#'
#' # Using custom colors (named vector)
#' my_colors <- c("Cluster0" = "#FF0000", "Cluster1" = "#00FF00")
#' avgHeatmap(seurat, selGenes = c("GeneA", "GeneB"), annotation_colors = my_colors)
#'
#' # Using custom colors (unnamed vector - will be matched to groups in order)
#' my_colors <- c("#FF0000", "#00FF00", "#0000FF")
#' avgHeatmap(seurat, selGenes = c("GeneA", "GeneB"), annotation_colors = my_colors)
#' }
#'
#' @export

avgHeatmap <- function(seurat,
                       selGenes = NULL,
                       group_by = NULL,
                       scale_method = "row",
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       cellwidth = 15,
                       cellheight = 10,
                       color_palette = NULL,
                       annotation_colors = NULL,
                       gaps_row = NULL,
                       gaps_col = NULL,
                       n_variable_genes = 20,
                       ...) {

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

  # Create expression matrix averaged by identity
  logNormExpres <- as.data.frame(t(as.matrix(
    seuratDat[matched_genes$gene, ]
  )))

  logNormExpres <- logNormExpres %>%
    mutate(cell = rownames(.)) %>%
    left_join(clusterAssigned, by = "cell") %>%
    dplyr::select(-cell) %>%
    group_by(ident) %>%
    summarise_all(mean, .groups = 'drop')

  # Convert to matrix format for heatmap
  logNormExpresMa <- logNormExpres %>%
    dplyr::select(-ident) %>%
    as.matrix()
  rownames(logNormExpresMa) <- logNormExpres$ident
  logNormExpresMa <- t(logNormExpresMa)

  # Clean up gene names (remove ENSEMBL prefix if present)
  rownames(logNormExpresMa) <- gsub("^.*?\\.", "", rownames(logNormExpresMa))

  # Remove genes with zero variance across all groups
  zero_var_genes <- apply(logNormExpresMa, 1, sd) == 0
  if (any(zero_var_genes)) {
    logNormExpresMa <- logNormExpresMa[!zero_var_genes, , drop = FALSE]
    warning(paste("Removed", sum(zero_var_genes), "genes with zero variance across groups"))
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

  # Create annotation for columns (cell types/groups)
  annotation_col <- data.frame(
    Group = colnames(logNormExpresMa)
  )
  rownames(annotation_col) <- colnames(logNormExpresMa)

  # Get the actual groups present in the data (after averaging)
  groups_present <- unique(annotation_col$Group)
  n_groups <- length(groups_present)

  # Set default annotation colors if not provided
  if (is.null(annotation_colors)) {
    # Define nice color palettes
    nice_colors <- c("#67001f", "#D53E4F", "#f4a582", "#FEE08B", "#003c30", "#01665e",
                     "#66C2A5", "#3288BD", "#BEAED4", "#c7eae5", "#355C7D", "#202547",
                     "#B45B5C", "#8c510a")

    disease_colors <- c("#dfc27d", "#BE3144", "#202547", "#355C7D", "#779d8d")

    fibroblast_colors <- c("#D53E4F", "#f4a582", "#ff7b7b", "#8e0b00", "#FEE08B",
                           "#42090D", "#FF7B00", "#FFF4DF")

    extended_colors <- c("#fde0dd", "#fa9fb5", "#d95f0e", "#dd1c77", "#D53E4F",
                         "#f4a582", "#FEE08B", "#f03b20", "#ffffcc", "#43a2ca",
                         "#1c9099", "#355C7D", "#3288BD", "#BEAED4", "#756bb1",
                         "#c7eae5")

    # Choose color palette based on group names and count
    if (n_groups <= 5 && any(grepl("healthy|explant|visit", groups_present, ignore.case = TRUE))) {
      # Disease condition colors
      group_colors <- disease_colors[1:n_groups]
    } else if (n_groups <= 8 && any(grepl("Fb|Periv|VSMC", groups_present, ignore.case = TRUE))) {
      # Fibroblast colors
      group_colors <- fibroblast_colors[1:n_groups]
    } else if (n_groups <= 14) {
      # General nice colors for clusters
      group_colors <- nice_colors[1:n_groups]
    } else if (n_groups <= 16) {
      # Extended colors for 15-16 groups
      group_colors <- extended_colors[1:n_groups]
    } else {
      # For many groups, use extended palette + rainbow
      group_colors <- c(extended_colors, rainbow(n_groups - length(extended_colors)))
    }

    names(group_colors) <- groups_present
    annotation_colors <- list(Group = group_colors)

  } else {
    # User provided annotation_colors - need to process it

    # If user provided a vector, wrap it in a list
    if (is.vector(annotation_colors) && !is.list(annotation_colors)) {

      # Check if the vector has names
      if (!is.null(names(annotation_colors))) {
        # Named vector - filter to only groups present in data
        group_colors <- annotation_colors[names(annotation_colors) %in% groups_present]

        # Check if we have colors for all groups
        if (length(group_colors) < n_groups) {
          missing_groups <- setdiff(groups_present, names(group_colors))
          warning("Some groups are missing from annotation_colors. ",
                  "Missing: ", paste(missing_groups, collapse = ", "))

          # Add default colors for missing groups
          default_colors <- rainbow(length(missing_groups))
          names(default_colors) <- missing_groups
          group_colors <- c(group_colors, default_colors)
        }

      } else {
        # Unnamed vector - assign to groups in order
        if (length(annotation_colors) < n_groups) {
          warning("annotation_colors has ", length(annotation_colors),
                  " colors but there are ", n_groups, " groups. ",
                  "Adding default colors for remaining groups.")
          # Extend with default colors
          extra_colors <- rainbow(n_groups - length(annotation_colors))
          annotation_colors <- c(annotation_colors, extra_colors)
        } else if (length(annotation_colors) > n_groups) {
          message("annotation_colors has more colors than groups. Using first ", n_groups, " colors.")
          annotation_colors <- annotation_colors[1:n_groups]
        }

        # Assign names based on groups present
        group_colors <- annotation_colors[1:n_groups]
        names(group_colors) <- groups_present
      }

      annotation_colors <- list(Group = group_colors)

    } else if (is.list(annotation_colors)) {
      # User provided a list

      # Try to find the right element
      if ("Group" %in% names(annotation_colors)) {
        group_colors <- annotation_colors$Group
      } else if (length(annotation_colors) > 0) {
        # Use first element
        group_colors <- annotation_colors[[1]]
        message("Using first element of annotation_colors list for Group colors")
      } else {
        stop("annotation_colors list is empty")
      }

      # Filter to groups present
      if (!is.null(names(group_colors))) {
        group_colors <- group_colors[names(group_colors) %in% groups_present]

        # Check if we need to add colors for missing groups
        if (length(group_colors) < n_groups) {
          missing_groups <- setdiff(groups_present, names(group_colors))
          default_colors <- rainbow(length(missing_groups))
          names(default_colors) <- missing_groups
          group_colors <- c(group_colors, default_colors)
        }
      } else {
        # Unnamed - assign to groups
        if (length(group_colors) < n_groups) {
          extra_colors <- rainbow(n_groups - length(group_colors))
          group_colors <- c(group_colors, extra_colors)
        }
        names(group_colors) <- groups_present[1:length(group_colors)]
      }

      annotation_colors <- list(Group = group_colors)
    }
  }

  # Generate and return heatmap
  p <- pheatmap(logNormExpresMa,
                scale = scale_method,
                cluster_rows = cluster_rows,
                cluster_cols = cluster_cols,
                color = color_palette,
                annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                cellwidth = cellwidth,
                cellheight = cellheight,
                show_rownames = show_rownames,
                show_colnames = show_colnames,
                gaps_row = gaps_row,
                gaps_col = gaps_col,
                ...)

  return(p)
}
