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
#' @param group_legend_title Optional string to set the title for the main group/cluster legend.
#'   Defaults to "Group" (or "Celltype" if `condition_by` is used).
#' @param condition_legend_title Optional string to set the title for the condition legend.
#'   Only used if `condition_by` is specified. Defaults to "Condition".
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
#'     orders genes by cluster of maximum expression to create a "staircase" pattern.
#' - If `gene_order` is provided, genes are kept in the user-provided order.
#' - If `cluster_order` is provided, clusters are kept in the user-provided order.
#'     (When `condition_by` is used, this groups columns based on the cluster prefix).
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
#' # Example with custom legend titles
#' avgHeatmap(seurat,
#'            selGenes = c("GeneA", "GeneB"),
#'            group_by = "celltype",
#'            condition_by = "treatment",
#'            group_legend_title = "Cell Type",
#'            condition_legend_title = "Treatment")
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
                       group_legend_title = NULL,
                       condition_legend_title = NULL,
                       n_variable_genes = 20,
                       show_significance = FALSE,
                       significance_direction = "higher",
                       significance_test = "wilcox",
                       pval_cutoffs = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
                       star_size = 8,
                       p_adjust_method = "BH",
                       colVecIdent = NULL,
                       colVecCond = NULL,
                       ordVec = NULL,
                       gapVecR = NULL,
                       gapVecC = NULL,
                       cc = NULL,
                       cr = NULL,
                       condCol = NULL,
                       return_ggplot = FALSE,
                       ...) {

  # ===========================================================================
  # BACKWARD COMPATIBILITY
  # ===========================================================================
  if (!is.null(colVecIdent)) {
    message("Note: 'colVecIdent' is deprecated. Please use 'annotation_colors' instead.")
    if (is.null(annotation_colors)) annotation_colors <- colVecIdent
  }
  if (!is.null(colVecCond)) {
    message("Note: 'colVecCond' is deprecated. Please use 'condition_colors' instead.")
    if (is.null(condition_colors)) condition_colors <- colVecCond
  }
  if (!is.null(ordVec)) {
    message("Note: 'ordVec' is deprecated. Please use 'cluster_order' instead.")
    if (is.null(cluster_order)) cluster_order <- ordVec
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
    if (is.null(gaps_row)) gaps_row <- gapVecR
  }
  if (!is.null(gapVecC)) {
    message("Note: 'gapVecC' is deprecated. Please use 'gaps_col' instead.")
    if (is.null(gaps_col)) gaps_col <- gapVecC
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
  # MAIN FUNCTION LOGIC & DATA PREP
  # ===========================================================================
  if (is.null(selGenes)) {
    variable_features <- VariableFeatures(seurat, assay = ifelse("SCT" %in% names(seurat@assays), "SCT", "RNA"))
    if (length(variable_features) == 0) stop("No variable features found. Please run FindVariableFeatures() first.")
    gene_list <- head(variable_features, n_variable_genes)
    message("No genes provided. Using top ", n_variable_genes, " variable features...")
  } else if (is.data.frame(selGenes)) {
    gene_list <- if("gene" %in% colnames(selGenes)) selGenes$gene else if("geneID" %in% colnames(selGenes)) selGenes$geneID else selGenes[, 1]
  } else {
    gene_list <- selGenes
  }

  if (is.null(group_by)) {
    group_by <- "ident"
    clusterAssigned <- data.frame(ident = Idents(seurat), cell = names(Idents(seurat)))
  } else {
    clusterAssigned <- data.frame(ident = seurat@meta.data[[group_by]], cell = rownames(seurat@meta.data))
  }

  if (!is.null(condition_by)) {
    if (!condition_by %in% colnames(seurat@meta.data)) stop("condition_by '", condition_by, "' not found in metadata")
    clusterAssigned$condition <- seurat@meta.data[[condition_by]]
    clusterAssigned$combined_ident <- paste0(clusterAssigned$ident, "_", clusterAssigned$condition)
  }

  seuratDat <- tryCatch(
    GetAssayData(seurat, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(seurat, assay = "RNA", slot = "data")
  )

  genes <- data.frame(gene = rownames(seurat)) %>% mutate(geneID = gsub("^.*\\.", "", gene))
  matched_genes <- if(all(grepl("^ENSG", gene_list))) genes[genes$gene %in% gene_list, ] else genes[genes$geneID %in% gene_list | genes$gene %in% gene_list, ]
  if (nrow(matched_genes) == 0) stop("No matching genes found in the Seurat object.")

  logNormExpres <- as.data.frame(t(as.matrix(seuratDat[matched_genes$gene, ]))) %>%
    mutate(cell = rownames(.)) %>%
    left_join(clusterAssigned, by = "cell") %>%
    dplyr::select(-cell)

  if (!is.null(condition_by)) {
    logNormExpres <- logNormExpres %>% group_by(combined_ident) %>% summarise_at(vars(-ident, -condition), mean, .groups = 'drop')
    logNormExpresMa <- logNormExpres %>% dplyr::select(-combined_ident) %>% as.matrix()
    rownames(logNormExpresMa) <- logNormExpres$combined_ident
  } else {
    logNormExpres <- logNormExpres %>% group_by(ident) %>% summarise_all(mean, .groups = 'drop')
    logNormExpresMa <- logNormExpres %>% dplyr::select(-ident) %>% as.matrix()
    rownames(logNormExpresMa) <- logNormExpres$ident
  }

  logNormExpresMa <- t(logNormExpresMa)
  rownames(logNormExpresMa) <- gsub("^.*?\\.", "", rownames(logNormExpresMa))

  zero_var_genes <- apply(logNormExpresMa, 1, sd) == 0
  if (any(zero_var_genes)) {
    logNormExpresMa <- logNormExpresMa[!zero_var_genes, , drop = FALSE]
    warning(paste("Removed", sum(zero_var_genes), "genes with zero variance"))
  }

  # ===========================================================================
  # SIGNIFICANCE TESTING
  # ===========================================================================
  significance_matrix <- NULL
  if (show_significance) {
    message("Performing significance testing...")
    significance_matrix <- matrix("", nrow = nrow(logNormExpresMa), ncol = ncol(logNormExpresMa))
    rownames(significance_matrix) <- rownames(logNormExpresMa)
    colnames(significance_matrix) <- colnames(logNormExpresMa)

    raw_expression <- as.data.frame(t(as.matrix(seuratDat[matched_genes$gene, ]))) %>%
      mutate(cell = rownames(.)) %>% left_join(clusterAssigned, by = "cell")

    gene_cols <- setdiff(colnames(raw_expression), c("cell", "ident", "condition", "combined_ident"))
    clean_gene_names <- gsub("^.*?\\.", "", gene_cols)
    group_var <- if (!is.null(condition_by)) "combined_ident" else "ident"

    pval_list <- list(); test_info <- list(); test_counter <- 0

    for (cluster in colnames(logNormExpresMa)) {
      for (i in seq_along(gene_cols)) {
        gene_col <- gene_cols[i]; gene_name <- clean_gene_names[i]
        if (!gene_name %in% rownames(logNormExpresMa)) next

        cluster_cells <- raw_expression[[group_var]] == cluster
        in_cluster <- raw_expression[cluster_cells, gene_col]

        if (!is.null(condition_by)) {
          base_cluster <- gsub("_[^_]*$", "", cluster)
          all_heatmap_clusters <- colnames(logNormExpresMa)
          sibling_clusters <- grep(paste0("^", base_cluster, "_"), all_heatmap_clusters, value = TRUE)
          other_sibling_clusters <- setdiff(sibling_clusters, cluster)
          out_cells <- raw_expression[[group_var]] %in% other_sibling_clusters
          out_cluster <- raw_expression[out_cells, gene_col]
        } else {
          out_cluster <- raw_expression[!cluster_cells, gene_col]
        }

        if (length(in_cluster) < 3 || length(out_cluster) < 3) next

        tryCatch({
          if (significance_test == "wilcox") {
            test_result <- wilcox.test(in_cluster, out_cluster, alternative = "two.sided")
          } else {
            test_result <- t.test(in_cluster, out_cluster, alternative = "two.sided")
          }

          test_counter <- test_counter + 1
          pval_list[[test_counter]] <- test_result$p.value
          test_info[[test_counter]] <- list(gene = gene_name, cluster = cluster, is_higher = median(in_cluster) > median(out_cluster))
        }, error = function(e) {})
      }
    }

    if (length(pval_list) > 0) {
      if (p_adjust_method == "none") {
        adjusted_pvals <- unlist(pval_list)
      } else {
        adjusted_pvals <- p.adjust(unlist(pval_list), method = p_adjust_method)
      }

      for (j in seq_along(adjusted_pvals)) {
        info <- test_info[[j]]; pval <- adjusted_pvals[j]
        should_show <- (significance_direction == "both") ||
          (significance_direction == "higher" && info$is_higher) ||
          (significance_direction == "lower" && !info$is_higher)

        if (should_show) {
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

  # Respect factor level ordering when no explicit cluster_order is given
  if (is.null(cluster_order)) {
    if (is.factor(clusterAssigned$ident)) {
      base_levels <- levels(clusterAssigned$ident)
    } else {
      base_levels <- unique(clusterAssigned$ident)
    }

    if (!is.null(condition_by)) {
      if (is.factor(clusterAssigned$condition)) {
        cond_levels <- levels(clusterAssigned$condition)
      } else {
        cond_levels <- unique(clusterAssigned$condition)
      }
      desired_order <- as.vector(t(outer(base_levels, cond_levels, paste, sep = "_")))
      desired_order <- desired_order[desired_order %in% colnames(logNormExpresMa)]
    } else {
      desired_order <- base_levels[base_levels %in% colnames(logNormExpresMa)]
    }

    if (length(desired_order) == ncol(logNormExpresMa)) {
      logNormExpresMa <- logNormExpresMa[, desired_order, drop = FALSE]
      if (!is.null(significance_matrix)) {
        significance_matrix <- significance_matrix[, desired_order, drop = FALSE]
      }
    }
  }

  # ===========================================================================
  # Column and Row Ordering
  # ===========================================================================
  if (!is.null(cluster_order)) {
    if (!is.null(condition_by)) {
      final_order <- c(); current_cols <- colnames(logNormExpresMa)
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
        if (!is.null(significance_matrix)) significance_matrix <- significance_matrix[, final_order, drop = FALSE]
      }
    } else {
      order_valid <- cluster_order[cluster_order %in% colnames(logNormExpresMa)]
      if (length(order_valid) > 0) {
        final_order <- c(order_valid, setdiff(colnames(logNormExpresMa), order_valid))
        logNormExpresMa <- logNormExpresMa[, final_order, drop = FALSE]
        if (!is.null(significance_matrix)) significance_matrix <- significance_matrix[, final_order, drop = FALSE]
      }
    }
  }

  if (!cluster_rows) {
    if (!is.null(gene_order)) {
      cleaned_gene_list <- gsub("^.*?\\.", "", gene_order)
      ordered_genes_from_list <- cleaned_gene_list[cleaned_gene_list %in% rownames(logNormExpresMa)]
      final_gene_order <- c(ordered_genes_from_list, setdiff(rownames(logNormExpresMa), ordered_genes_from_list))
      if (length(final_gene_order) > 0) {
        logNormExpresMa <- logNormExpresMa[final_gene_order, , drop = FALSE]
        if (!is.null(significance_matrix)) significance_matrix <- significance_matrix[final_gene_order, , drop = FALSE]
      }
    } else {
      max_clusters <- apply(logNormExpresMa, 1, function(x) colnames(logNormExpresMa)[which.max(x)])
      ordered_genes <- c()
      for (cluster in colnames(logNormExpresMa)) {
        cluster_genes <- names(max_clusters[max_clusters == cluster])
        if (length(cluster_genes) > 0) {
          ordered_genes <- c(ordered_genes, cluster_genes[order(logNormExpresMa[cluster_genes, cluster], decreasing = TRUE)])
        }
      }
      if (length(ordered_genes) > 0) {
        logNormExpresMa <- logNormExpresMa[ordered_genes, , drop = FALSE]
        if (!is.null(significance_matrix)) significance_matrix <- significance_matrix[ordered_genes, , drop = FALSE]
      }
    }
  }

  # ===========================================================================
  # Color Setup (Refactored to use generate_colors)
  # ===========================================================================
  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(50)
  }

  if (!is.null(condition_by)) {
    col_names <- colnames(logNormExpresMa)
    celltype_vec <- gsub("_[^_]*$", "", col_names)
    condition_vec <- gsub("^.*_", "", col_names)

    final_group_name <- if (!is.null(group_legend_title)) group_legend_title else "Celltype"
    final_cond_name <- if (!is.null(condition_legend_title)) condition_legend_title else "Condition"

    annotation_col <- data.frame(GroupColumn = celltype_vec, ConditionColumn = condition_vec, row.names = col_names)
    colnames(annotation_col) <- c(final_group_name, final_cond_name)

    celltypes_present <- unique(celltype_vec)
    conditions_present <- unique(condition_vec)

    # 1. Resolve Celltype Colors
    if (!is.null(annotation_colors)) {
      celltype_colors <- if (is.list(annotation_colors)) annotation_colors[[1]] else annotation_colors
      celltype_colors <- celltype_colors[names(celltype_colors) %in% celltypes_present]
      if (length(celltype_colors) < length(celltypes_present)) {
        missing <- setdiff(celltypes_present, names(celltype_colors))
        celltype_colors <- c(celltype_colors, generate_colors(missing, "default"))
      }
    } else {
      celltype_colors <- generate_colors(celltypes_present, "default")
    }

    # 2. Resolve Condition Colors
    if (!is.null(condition_colors)) {
      # Handle unnamed vectors: assign names from conditions_present
      if (is.null(names(condition_colors))) {
        if (length(condition_colors) >= length(conditions_present)) {
          names(condition_colors) <- conditions_present[seq_along(conditions_present)]
        }
      }
      cond_colors <- condition_colors[names(condition_colors) %in% conditions_present]
      if (length(cond_colors) < length(conditions_present)) {
        missing <- setdiff(conditions_present, names(cond_colors))
        cond_colors <- c(cond_colors, generate_colors(missing, "condition"))
      }
    }else {
      cond_colors <- generate_colors(conditions_present, "condition")
    }

    ann_colors_final <- list(celltype_colors, cond_colors)
    names(ann_colors_final) <- c(final_group_name, final_cond_name)

  } else {
    final_group_name <- if (!is.null(group_legend_title)) group_legend_title else "Group"
    annotation_col <- data.frame(Group = colnames(logNormExpresMa), row.names = colnames(logNormExpresMa))
    colnames(annotation_col) <- final_group_name
    groups_present <- unique(annotation_col[[final_group_name]])

    # 3. Resolve Group Colors (Single Annotation)
    if (is.null(annotation_colors)) {
      group_colors <- generate_colors(groups_present, "default")
    } else {
      group_colors <- if (is.list(annotation_colors)) annotation_colors[[1]] else annotation_colors
      group_colors <- group_colors[names(group_colors) %in% groups_present]
      if (length(group_colors) < length(groups_present)) {
        missing <- setdiff(groups_present, names(group_colors))
        group_colors <- c(group_colors, generate_colors(missing, "default"))
      }
    }
    ann_colors_final <- list(group_colors)
    names(ann_colors_final) <- final_group_name
  }

  # ===========================================================================
  # Generate Heatmap
  # ===========================================================================
  pheatmap_params <- list(
    mat = logNormExpresMa, scale = scale_method, cluster_rows = cluster_rows,
    cluster_cols = cluster_cols, color = color_palette, annotation_col = annotation_col,
    annotation_colors = ann_colors_final, cellwidth = cellwidth, cellheight = cellheight,
    show_rownames = show_rownames, show_colnames = show_colnames, gaps_row = gaps_row, gaps_col = gaps_col
  )

  if (show_significance && !is.null(significance_matrix)) {
    pheatmap_params$display_numbers <- significance_matrix
    pheatmap_params$number_color <- "black"
    pheatmap_params$fontsize_number <- star_size
  }

  pheatmap_params <- c(pheatmap_params, list(...))
  p <- do.call(pheatmap::pheatmap, pheatmap_params)

  if (return_ggplot) return(ggplotify::as.ggplot(p$gtable)) else return(p)
}
