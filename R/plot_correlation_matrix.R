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
#'  or \code{"spearman"}).
#' @param method Method passed to \code{corrplot::corrplot()} (e.g. "ellipse").
#' @param sig.level Significance threshold for highlighting.
#' @param order Variable ordering method.
#' @param pval_color_low Color for low (significant) p-values. Default is NULL (uses viridis).
#' @param pval_color_high Color for high (less significant) p-values. Default is NULL (uses viridis).
#' @param ... Additional arguments.
#'
#' @return Invisibly returns a data frame with coordinates and p-values.
#' @author Roman Stadler
#' @export
#'
#' @examples
#' \dontrun{
#' data(mtcars)
#' # Default viridis colors
#' plot_correlation_matrix(mtcars)
#'
#' # Custom colors (e.g., white to purple)
#' plot_correlation_matrix(mtcars, pval_color_low = "white", pval_color_high = "purple")
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
                                    legend_title_offset = 0.3,
                                    pval_color_low = NULL,
                                    pval_color_high = NULL) {


  if (!requireNamespace("Hmisc", quietly = TRUE) ||
      !requireNamespace("corrplot", quietly = TRUE)) {
    stop("Please install 'Hmisc' and 'corrplot'.")
  }

  # Check for viridis only if using default colors
  if (is.null(pval_color_low) || is.null(pval_color_high)) {
    if (!requireNamespace("viridis", quietly = TRUE)) {
      stop("Please install 'viridis' or specify custom colors with pval_color_low and pval_color_high.")
    }
  }

  # Helper function to clean variable names
  clean_variable_names <- function(names) {
    new_name <- gsub("k__", "", names) # also clean taxonomic levels
    new_name <- gsub("p__", "", new_name)
    new_name <- gsub("c__", "", new_name)
    new_name <- gsub("o__", "", new_name)
    new_name <- gsub("f__", "", new_name)
    new_name <- gsub("g__", "", new_name)
    new_name <- gsub("s__", "", new_name)

    new_name <- gsub("_", " ", new_name)
    split_words <- strsplit(new_name, " ")
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

    # Warn if there are non-numeric columns (categorical data)
    non_numeric_cols <- names(data)[!sapply(data, is.numeric)]
    if (length(non_numeric_cols) > 0) {
      warning(
        paste0(
          "The following columns are non-numeric and will be excluded from the correlation plot: ",
          paste(non_numeric_cols, collapse = ", "),
          "\nPlease check your data or convert these columns if correlations are desired."
        ),
        call. = FALSE
      )
      data <- data[, sapply(data, is.numeric), drop = FALSE]
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

  # Store original par settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Calculate legend space needed (as fraction of plot width)
  n_vars <- ncol(cor_matrix)
  has_significant <- any(!is.na(p_matrix) & p_matrix < sig.level)
  n_legends <- if(has_significant) 2 else 1

  # Adjust right margin to accommodate legends
  # Each legend needs space for: width + text + gap
  legend_space_per <- legend_width + legend_text_margin + 0.3
  total_legend_space <- n_legends * legend_space_per

  # Set margins: bottom, left, top, right (in margin lines)
  # Minimize left margin, increase right margin for legends
  par(mar = c(1, 0.5, 1, 8), xpd = FALSE)

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

  # Now allow drawing outside plot region for legends
  par(xpd = TRUE)

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

  upper_triangle$pval_text <- sapply(upper_triangle$pval, formatPval)


  # Filter for significant p-values (< 0.05) for colored circles
  significant <- upper_triangle[!is.na(upper_triangle$pval) & upper_triangle$pval < sig.level, ]

  # Create color mapping for significant p-values
  # Use custom colors if provided, otherwise use viridis
  if (nrow(significant) > 0) {
    # Map p-values from 0 to 0.05 to colors
    # Lower p-values get more intense colors
    p_range <- c(min(upper_triangle$pval), 0.05)
    # Normalize p-values to 0-1 scale, then invert so lower p-values = higher color intensity
    color_scale <- 1 - (significant$pval - p_range[1]) / (p_range[2] - p_range[1])
    # color_scale <- color_scale**2 # extra color boost possible

    # Generate color palette based on user input or default to viridis
    if (!is.null(pval_color_low) && !is.null(pval_color_high)) {
      color_palette <- colorRampPalette(c(pval_color_high, pval_color_low))(100)
    } else {
      library(viridis)
      color_palette <- viridis(100)
    }

    circle_colors <- color_palette[pmax(1, pmin(100, round(color_scale * 100)))]
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
       paste(stringr::str_to_title(correlation_type),"\nCorrelation"), cex = number_size, font = 2)
  text(legend_x_start + legend_width + legend_text_margin, legend_y_start, "1", cex = number_size)
  text(legend_x_start + legend_width + legend_text_margin, legend_y_start - legend_height/2, "0", cex = number_size)
  text(legend_x_start + legend_width + legend_text_margin, legend_y_start - legend_height, "-1", cex = number_size)

  # Add custom p-value legend if there are significant values
  if (nrow(significant) > 0) {
    pval_legend_x_start <- legend_x_start + legend_width + 0.8

    # Create color gradient based on user input or default to viridis
    if (!is.null(pval_color_low) && !is.null(pval_color_high)) {
      pval_legend_colors <- rev(colorRampPalette(c(pval_color_high, pval_color_low))(100))
    } else {
      library(viridis)
      pval_legend_colors <- rev(viridis(100))
    }

    # Draw p-value legend rectangles
    for (i in 1:100) {
      y_pos <- legend_y_start - (i-1) * legend_step
      rect(pval_legend_x_start, y_pos - legend_step,
           pval_legend_x_start + legend_width, y_pos,
           col = pval_legend_colors[i], border = NA)
    }

    # Add p-value legend labels and title
    text(pval_legend_x_start + legend_width/2, legend_y_start + legend_title_offset,
         "P-value", cex = number_size, font = 2)
    text(pval_legend_x_start + legend_width + legend_text_margin, legend_y_start, "0.0", cex = number_size)
    text(pval_legend_x_start + legend_width + legend_text_margin, legend_y_start - legend_height/2, paste(sig.level/2), cex = number_size)
    text(pval_legend_x_start + legend_width + legend_text_margin, legend_y_start - legend_height, paste(sig.level), cex = number_size)
  }

  # Crop the plot to remove excess white space
  plot_right_edge <- legend_x_start + legend_width * n_legends + 1
  par(usr = c(plot_dims[1], plot_right_edge, plot_dims[3], plot_dims[4]))

  invisible(pos)
}
