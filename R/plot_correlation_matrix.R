#' Plot a correlation matrix with significance overlays
#'
#' Creates a detailed correlation plot from a data frame or precomputed
#' correlation/p-value matrices. The plot includes correlation coefficients,
#' significance highlighting, and legends for both correlations and p-values.
#'
#' @param data Optional data frame. If provided, correlations and p-values
#'   are computed automatically.
#' @param cor_matrix Optional correlation matrix if \code{data} is not given.
#' @param p_matrix Optional matrix of p-values.
#' @param correlation_type Type of correlation to compute (\code{"pearson"},
#'   \code{"spearman"}, or \code{"kendall"}).
#' @param method Method passed to \code{corrplot::corrplot()} (e.g. "ellipse").
#' @param sig.level Significance threshold for highlighting.
#' @param order Variable ordering method.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns a data frame with coordinates and p-values.
#' @author Roman Stadler
#' @export
#'
#' @examples
#' \dontrun{
#' data(mtcars)
#' plot_correlation_matrix(mtcars)
#' }

plot_correlation_matrix <- function(data = NULL,
                                    cor_matrix = NULL,
                                    p_matrix = NULL,
                                    correlation_type = "pearson",
                                    method = "ellipse",
                                    sig.level = 0.05,
                                    order = "original",
                                    number_size = 0.9,
                                    number.font = 2,
                                    insig = "blank",
                                    label.col = "black",
                                    impute = F,
                                    verbose = T,
                                    clean_names = F,
                                    legend_text_margin = 0.15,
                                    legend_width = 0.15,
                                    legend_height= 0.6,
                                    legend_x_start=1,
                                    legend_y_start=0.9,
                                    legend_title_offset = 0.3) {


  if (!requireNamespace("Hmisc", quietly = TRUE) ||
      !requireNamespace("corrplot", quietly = TRUE) ||
      !requireNamespace("viridis", quietly = TRUE)) {
    stop("Please install 'Hmisc', 'corrplot' and 'viridis'.")
  }

  # if (!requireNamespace("biostatUZH", quietly = TRUE)) {
  #   stop(
  #     "Package 'biostatUZH' is required for p-value formatting.\n",
  #     "Install it with one of the following:\n",
  #     "  remotes::install_github(\"EBPI-Biostatistics/biostatUZH\")\n",
  #     "  install.packages(\"biostatUZH\", repos = \"http://R-Forge.R-project.org\")"
  #   )
  # }

  # Helper function to clean variable names
  clean_variable_names <- function(names) {
    names <- gsub("_", " ", names)
    split_words <- strsplit(names, " ")
    cleaned <- lapply(split_words, function(words) {
      sapply(words, function(w) {
        if (nchar(w) <= 3) {
          toupper(w)               # Make short words ALL CAPS
        } else {
          tools::toTitleCase(w)    # Title case for longer words
        }
      })
    })
    cleaned <- sapply(cleaned, paste, collapse = " ")
    return(cleaned)
  }

  if (!is.null(data)) {

    # Remove columns that are entirely NA
    all_na_cols <- names(data)[colSums(is.na(data)) == nrow(data)]
    if (length(all_na_cols) > 0) {
      cat("Removing columns with only NA values:", paste(all_na_cols, collapse = ", "), "\n")
    }
    data <- data[, colSums(is.na(data)) < nrow(data), drop = FALSE]
    cat("After removing all-NA columns:", dim(data), "\n")

    # Remove rows that are entirely NA
    all_na_rows <- which(rowSums(is.na(data)) == ncol(data))
    if (length(all_na_rows) > 0) {
      cat("Removing rows with only NA values:", paste(all_na_rows, collapse = ", "), "\n")
    }
    data <- data[rowSums(is.na(data)) < ncol(data), , drop = FALSE]
    cat("After removing all-NA rows:", dim(data), "\n")



    if(impute){
      if (!requireNamespace("mice", quietly = TRUE) ) {
        stop("Please install 'mice' for the imputation.")
      }
      if (verbose) cat("Performing MICE imputation...\n")
      library(mice)
      imp.data <- mice(data = data, m = 50, maxit = 10, seed = 12345, printFlag = F)
      data <- complete(imp.data)
      if (verbose) cat("Imputation completed\n")
    }

    cor_results <- Hmisc::rcorr(as.matrix(data), type=correlation_type)
    cor_matrix <- cor_results$r
    p_matrix <- cor_results$P

    # Clean names if requested
    if (clean_names) {
      cleaned_names <- clean_variable_names(rownames(cor_matrix))
      rownames(cor_matrix) <- cleaned_names
      colnames(cor_matrix) <- cleaned_names
      rownames(p_matrix) <- cleaned_names
      colnames(p_matrix) <- cleaned_names
    }

  } else if (is.null(cor_matrix) || is.null(p_matrix)) {
    stop("If 'data' is not provided, both 'cor_matrix' and 'p_matrix' must be supplied.")
  } else {
    # Clean names for provided matrices if requested
    if (clean_names) {
      cleaned_names <- clean_variable_names(rownames(cor_matrix))
      rownames(cor_matrix) <- cleaned_names
      colnames(cor_matrix) <- cleaned_names
      rownames(p_matrix) <- cleaned_names
      colnames(p_matrix) <- cleaned_names
    }
  }

  cp <- corrplot::corrplot(cor_matrix,
                           p.mat = p_matrix,
                           method = method,
                           type = "full",
                           order = order,
                           insig = "blank",
                           diag = FALSE,
                           tl.col = label.col,
                           number.cex = number_size,
                           number.font = number.font,
                           cl.pos = "n",
                           sig.level = sig.level,
                           col= colorRampPalette(c("blue", "white", "red"))(100))
  cp$corrPos -> p1
  # text(p1$x, p1$y, round(p1$corr, 2), cex = number_size)
  text_colors <- ifelse(!is.na(p1$corr) & p1$corr < -0.5 & p1$p.value < sig.level, "white", "black")

  # Add correlations on bottom
  text(p1$x, p1$y,
       labels =  round(p1$corr, 2),
       cex = number_size,
       col = text_colors,
       font = 1)

  # Extract coordinates and add p-values
  pos <- cp$corrPos
  row_names <- rownames(cor_matrix)
  col_names <- colnames(cor_matrix)


  get_p <- function(xname, yname) {
    i <- match(yname, row_names)
    j <- match(xname, col_names)
    return(p_matrix[i, j])
  }
  pos$pval <- mapply(get_p, pos$xName, pos$yName)
  pos

  # Filter for upper triangle only (where x > y)
  upper_triangle <- pos[pos$x + pos$y -1 > length(unique(pos$y)), ]
  upper_triangle


  library(biostatUZH)
  upper_triangle$pval_text <- sapply(upper_triangle$pval, formatPval)


  # Filter for significant p-values (< 0.05) for colored circles
  significant <- upper_triangle[!is.na(upper_triangle$pval) & upper_triangle$pval < sig.level, ]

  # Create viridis color mapping for significant p-values
  library(viridis)
  if (nrow(significant) > 0) {
    # Map p-values from 0 to 0.05 to viridis colors
    # Lower p-values get more intense colors
    p_range <- c(min(upper_triangle$pval), 0.05)
    # Normalize p-values to 0-1 scale, then invert so lower p-values = higher color intensity
    color_scale <- 1 - (significant$pval - p_range[1]) / (p_range[2] - p_range[1])
    # color_scale <- color_scale**2 # extra color boost possible
    circle_colors <- viridis(100)[pmax(1, pmin(100, round(color_scale * 100)))]
  }

  # Clear the upper triangle by drawing white rectangles
  margin <- 0.49
  for (i in 1:length(unique(pos$y)) ) {
    for (j in 1:length(unique(pos$y)) ) {
      if(i+j-1 > length(unique(pos$y))) {
        # cat(i, j, "\n")
        rect(i - margin, j - margin,
             i + margin, j + margin,
             col = "white", border = NA)
      }
    }
  }

  # Add colored circles for significant p-values
  if (nrow(significant) > 0) {
    circle_radius <- 0.35  # Adjust size as needed
    for (i in 1:nrow(significant)) {
      # Draw filled circle
      theta <- seq(0, 2*pi, length.out = 100)
      x_circle <- significant$x[i] + circle_radius * cos(theta)
      y_circle <- significant$y[i] + circle_radius * sin(theta)
      polygon(x_circle, y_circle, col = circle_colors[i], border = NA)
    }
  }

  # Add p-value text on top with adaptive color
  # Use white text for p-values < 0.03 (dark viridis background), black for others
  text_colors <- ifelse(!is.na(upper_triangle$pval) & upper_triangle$pval > 0.03 & upper_triangle$pval < sig.level, "white", "black")

  # Add p-value text on top
  text(upper_triangle$x, upper_triangle$y,
       labels = upper_triangle$pval_text,
       cex = number_size,
       col = text_colors,
       font = 1)


  # Get plot dimensions for legend positioning
  plot_dims <- par("usr")
  n_vars <- ncol(cor_matrix)

  # Add custom correlation legend (Red-White-Blue)
  legend_height <- n_vars * legend_height
  legend_x_start <- n_vars + legend_x_start
  legend_y_start <- n_vars * legend_y_start

  # Create correlation color gradient
  cor_colors <- c(
    colorRampPalette(c("red", "white"))(50),
    colorRampPalette(c("white", "blue"))(50)
  )

  # Draw correlation legend rectangles
  legend_step <- legend_height / 100
  for (i in 1:100) {
    y_pos <- legend_y_start - (i-1) * legend_step
    rect(legend_x_start, y_pos - legend_step,
         legend_x_start + legend_width, y_pos,
         col = cor_colors[i], border = NA)
  }

  require(stringr)
  # Add correlation legend labels and title
  text(legend_x_start + legend_width/2, legend_y_start + legend_title_offset,
       paste(stringr::str_to_title(correlation_type),"\nCorrelation"), cex = number_size, font = 2, xpd = TRUE)
  text(legend_x_start + legend_width + legend_text_margin, legend_y_start, "1", cex = number_size, xpd = TRUE)
  text(legend_x_start + legend_width + legend_text_margin, legend_y_start - legend_height/2, "0", cex = number_size, xpd = TRUE)
  text(legend_x_start + legend_width + legend_text_margin, legend_y_start - legend_height, "-1", cex = number_size, xpd = TRUE)

  # Add custom p-value legend (Viridis) if there are significant values
  if (nrow(significant) > 0) {
    pval_legend_x_start <- legend_x_start + legend_width + 0.8

    # Create viridis color gradient
    viridis_colors <- viridis(100)

    # Draw p-value legend rectangles
    for (i in 1:100) {
      y_pos <- legend_y_start - (i-1) * legend_step
      rect(pval_legend_x_start, y_pos - legend_step,
           pval_legend_x_start + legend_width, y_pos,
           col = viridis_colors[101-i], border = NA)  # Reverse order so low p-values are at top
    }

    # Add p-value legend labels and title
    text(pval_legend_x_start + legend_width/2, legend_y_start + legend_title_offset,
         "P-value", cex = number_size, font = 2, xpd = TRUE)
    text(pval_legend_x_start + legend_width + legend_text_margin, legend_y_start, "0.0", cex = number_size, xpd = TRUE)
    text(pval_legend_x_start + legend_width + legend_text_margin, legend_y_start - legend_height/2, paste(sig.level/2), cex = number_size, xpd = TRUE)
    text(pval_legend_x_start + legend_width + legend_text_margin, legend_y_start - legend_height, paste(sig.level), cex = number_size, xpd = TRUE)
  }

}

