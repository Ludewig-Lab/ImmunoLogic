#' Goodness-of-Fit Diagnostic Plots for MMLT Models
#'
#' Creates diagnostic plots to assess the goodness of fit of a multivariate
#' maximum likelihood transformation (MMLT) model using Rosenblatt transformations.
#' Plots empirical cumulative distribution functions (ECDFs) of transformed marginal
#' and conditional distributions and compares them against the theoretical uniform
#' distribution.
#'
#' @param new_data Data frame containing the variables used for goodness-of-fit testing.
#'   Must contain the same variables as used in fitting the MMLT model.
#' @param mmlt_object Fitted MMLT model object returned by \code{fit_mmlt_flexible} or \code{mmlt}.
#' @param plot_layout Numeric vector of length 2 specifying the plot grid dimensions
#'   (rows, columns). Default is \code{c(3, 3)} for a 3x3 grid.
#' @param exclude_variables Character vector of variable names to exclude from the
#'   goodness-of-fit assessment. Default is empty.
#' @param status Logical. If \code{TRUE}, prints progress messages during computation.
#'   Default is \code{FALSE}.
#' @param M Integer. Number of quasi-random points for multivariate normal probability
#'   calculations. Higher values increase accuracy but computation time. Default is 2000.
#' @param line_color Color for the reference diagonal line (theoretical uniform CDF).
#'   Default is \code{"red"}.
#' @param digits Integer. Number of decimal places for displayed p-values. Default is 3.
#' @param add_names Logical. If \code{TRUE}, adds variable names to plot titles.
#'   Default is \code{FALSE}.
#'
#' @return A list containing:
#'   \item{pvalmin}{Minimum p-value across all tests, adjusted using
#'     the beta distribution.}
#'   \item{pvalprod}{Product of all p-values, adjusted for multiple testing.}
#'   \item{individual_pvals}{Named vector of individual test p-values for each dimension.}
#'   \item{test_types}{Named vector indicating which test was used ("KS" or "CvM") for each dimension.}
#'   \item{transformed_data}{Named list of transformed uniform variates for each dimension.}
#'
#' @details
#' This function implements goodness-of-fit diagnostics based on Rosenblatt transformations.
#' For a J-dimensional model, it creates J plots showing:
#' \itemize{
#'   \item Plot 1: Marginal ECDF of \eqn{F_1(Y_1)}
#'   \item Plot i (2 ≤ i ≤ J): Conditional ECDF of \eqn{F_i(Y_i | Y_1, ..., Y_{i-1})}
#' }
#'
#' Under correct model specification, each transformed variate should follow a
#' uniform distribution on [0,1]. Deviations indicate model misspecification.
#'
#' The function automatically detects ties in the transformed data (which occur with
#' discrete or ordered variables) and uses the Cramér-von Mises test instead of the
#' Kolmogorov-Smirnov test when ties are present. For continuous variables without
#' ties, the standard KS test is used.
#'
#' Two omnibus tests are reported:
#' \itemize{
#'   \item \strong{Minimum p-value test}: Tests whether the minimum p-value is
#'     suspiciously small (adjusted via beta distribution)
#'   \item \strong{Product p-value test}: Combines all p-values using Fisher's method
#' }
#'
#' @note
#' \itemize{
#'   \item Requires the \code{tram} package and its dependencies
#'   \item Computation time increases with M and the number of dimensions
#'   \item For high-dimensional models, consider using a larger plot layout
#' }
#'
#' @author Roman Stadler
#'
#' @seealso \code{\link{fit_mmlt_flexible}}, \code{\link[tram]{mmlt}}
#'
#' @importFrom tram predict.mmlt coef
#' @importFrom stats ecdf ks.test pnorm pbeta
#' @importFrom qrng ghalton
#' @export
#'
#' @examples
#' \dontrun{
#' library(tram)
#' data <- data.frame(
#'   x1 = rnorm(200),
#'   x2 = rnorm(200),
#'   x3 = rnorm(200)
#' )
#'
#' # Fit MMLT model
#' model <- fit_mmlt_flexible(data, automatic = TRUE)
#'
#' # Assess goodness of fit
#' gof_results <- plot_gof(
#'   new_data = data,
#'   mmlt_object = model$object,
#'   plot_layout = c(2, 2),
#'   add_names = TRUE
#' )
#'
#' # Check overall fit
#' print(gof_results$pvalmin)
#' print(gof_results$pvalprod)
#' }

plot_gof <- function(new_data,
                     mmlt_object,
                     plot_layout = c(3, 3),
                     exclude_variables = c(),
                     status = FALSE,
                     M = 2000,
                     line_color = "red",
                     digits = 3,
                     add_names = FALSE) {

  # Load required packages
  if (!requireNamespace("tram", quietly = TRUE)) {
    stop("Package 'tram' is required but not installed.")
  }
  if (!requireNamespace("qrng", quietly = TRUE)) {
    stop("Package 'qrng' is required but not installed.")
  }

  # Input validation
  if (!is.data.frame(new_data)) {
    stop("'new_data' must be a data frame.")
  }
  if (nrow(new_data) < 10) {
    warning("Small sample size may lead to unreliable goodness-of-fit results.")
  }
  if (length(plot_layout) != 2 || any(plot_layout < 1)) {
    stop("'plot_layout' must be a numeric vector of length 2 with positive integers.")
  }
  if (M < 100) {
    warning("'M' is very small; consider using M >= 1000 for reliable results.")
  }

  # Exclude specified variables
  if (length(exclude_variables) > 0) {
    if (status) {
      cat("Excluding variables:", paste(exclude_variables, collapse = ", "), "\n")
    }
    new_data <- new_data[, !names(new_data) %in% exclude_variables, drop = FALSE]
  }

  variable_names <- colnames(new_data)
  J <- ncol(new_data)

  if (J < 2) {
    stop("At least 2 variables are required for goodness-of-fit assessment.")
  }

  if (status) {
    cat("=== Goodness-of-Fit Assessment ===\n")
    cat("Number of variables:", J, "\n")
    cat("Number of observations:", nrow(new_data), "\n")
    cat("Monte Carlo samples (M):", M, "\n\n")
  }

  # Compute marginal transformations
  if (status) cat("Computing marginal transformations...\n")
  for (i in 1:J) {
    lp_col <- paste0("lp", i)
    new_data[[lp_col]] <- predict(mmlt_object,
                                  newdata = new_data,
                                  margin = i,
                                  type = "trafo")
  }

  lcol <- paste0("lp", 1:J)

  # Set up plotting parameters
  par(mar = c(1, 0, 3, 0) + 0.1, oma = c(4, 4, 0.5, 0.5))
  par(mfrow = plot_layout)

  cxa <- 1      # Axis text size
  cxt <- 0.8    # P-value text size
  grd <- seq(0, 1, by = 0.001)
  xlim <- c(0, 1)
  ylim <- c(0, 1)

  # Storage for transformed distributions and p-values
  us <- vector("list", J)
  pvals <- numeric(J)
  test_types <- character(J)  # Track which test was used

  # ============================================================================
  # Plot 1: Marginal distribution of first variable
  # ============================================================================
  if (status) cat("Processing dimension 1 (marginal)...\n")

  u1 <- predict(mmlt_object,
                newdata = new_data,
                margins = 1,
                type = "distribution")
  us[[1]] <- u1
  e1 <- ecdf(u1)

  # Check for ties and use appropriate test
  has_ties_1 <- any(duplicated(u1))
  if (has_ties_1) {
    # Use Cramér-von Mises test for discrete data
    pvals[1] <- cvm_test(u1)
    test_types[1] <- "CvM"
  } else {
    pvals[1] <- suppressWarnings(ks.test(u1, "punif")$p.value)
    test_types[1] <- "KS"
  }

  var_name <- if (add_names) variable_names[1] else NULL
  plot_title <- get_title(1, var_name)
  pvt <- bquote(paste("", p-value == .(round(pvals[1], digits))))

  # Determine axis visibility
  x_index <- ((1 - 1) %/% plot_layout[2]) + 1
  y_index <- (1 - 1) %% plot_layout[2]
  show_x <- (x_index == plot_layout[1]) || (1 + plot_layout[2] > J)
  show_y <- (y_index == 0)

  plot(grd, e1(grd),
       xaxt = ifelse(show_x, "s", "n"),
       yaxt = ifelse(show_y, "s", "n"),
       las = 1,
       type = "s",
       cex.axis = cxa,
       cex.main = cxa,
       xlim = xlim,
       ylim = ylim,
       main = plot_title)

  mtext(pvt, side = 3, outer = FALSE, line = -1.5, adj = 0.03, cex = cxt)
  abline(0, 1, col = line_color)

  # ============================================================================
  # Plots 2 to J-1: Conditional distributions (if J > 2)
  # ============================================================================
  if (J > 2) {
    for (i in 2:(J - 1)) {
      if (status) cat("Processing dimension", i, "(conditional)...\n")

      which_given <- 1:(i - 1)

      # Compute conditional distribution
      cd <- cond_mvnorm(
        invchol = coef(mmlt_object, type = "Lambda"),
        which_given = which_given,
        given = t(as.matrix(new_data[, lcol[which_given], drop = FALSE]))
      )

      # Set up integration bounds
      lower <- upper <- matrix(0, nrow = J - i + 1, ncol = nrow(new_data))
      lower[1, ] <- -Inf
      upper[1, ] <- new_data[[lcol[i]]]

      if ((J - i + 1) > 1) {
        for (k in 2:(J - i + 1)) {
          lower[k, ] <- -Inf
          upper[k, ] <- Inf
        }
      }

      # Compute probabilities using quasi-Monte Carlo
      pmvargs <- list(M = M, w = t(ghalton(M, d = J - i)))
      u <- exp(
        lpmvnorm(
          lower = lower,
          upper = upper,
          mean = cd$mean,
          invchol = cd$invchol,
          w = pmvargs$w,
          logLik = FALSE
        )
      )

      us[[i]] <- u
      e <- ecdf(u)

      # Check for ties and use appropriate test
      has_ties <- any(duplicated(u))
      if (has_ties) {
        pvals[i] <- cvm_test(u)
        test_types[i] <- "CvM"
      } else {
        pvals[i] <- suppressWarnings(ks.test(u, "punif")$p.value)
        test_types[i] <- "KS"
      }

      var_name <- if (add_names) variable_names[i] else NULL
      plot_title <- get_title(i, var_name)
      pvt <- bquote(paste("", p-value == .(round(pvals[i], digits))))

      # Determine axis visibility
      x_index <- ((i - 1) %/% plot_layout[2]) + 1
      y_index <- (i - 1) %% plot_layout[2]
      show_x <- (x_index == plot_layout[1]) || (i + plot_layout[2] > J)
      show_y <- (y_index == 0)

      plot(grd, e(grd),
           xaxt = ifelse(show_x, "s", "n"),
           yaxt = ifelse(show_y, "s", "n"),
           las = 1,
           type = "s",
           cex.axis = cxa,
           cex.main = cxa,
           xlim = xlim,
           ylim = ylim,
           main = plot_title)

      mtext(pvt, side = 3, outer = FALSE, line = -1.5, adj = 0.03, cex = cxt)
      abline(0, 1, col = line_color)
    }
  }

  # ============================================================================
  # Plot J: Final conditional distribution
  # ============================================================================
  if (status) cat("Processing dimension", J, "(final conditional)...\n")

  which_given <- 1:(J - 1)
  cd <- cond_mvnorm(
    invchol = coef(mmlt_object, type = "Lambda"),
    which_given = which_given,
    given = t(as.matrix(new_data[, lcol[which_given], drop = FALSE]))
  )

  # For the last dimension, use normal CDF directly
  un <- pnorm(new_data[[lcol[J]]],
              mean = c(cd$mean),
              sd = sqrt(c(diagonals(invchol2cov(cd$invchol)))))

  us[[J]] <- un
  e <- ecdf(un)

  # Check for ties and use appropriate test
  has_ties_J <- any(duplicated(un))
  if (has_ties_J) {
    pvals[J] <- cvm_test(un)
    test_types[J] <- "CvM"
  } else {
    pvals[J] <- suppressWarnings(ks.test(un, "punif")$p.value)
    test_types[J] <- "KS"
  }

  var_name <- if (add_names) variable_names[J] else NULL
  plot_title <- get_title(J, var_name)
  pvt <- bquote(paste("", p-value == .(round(pvals[J], digits))))

  # Determine axis visibility
  x_index <- ((J - 1) %/% plot_layout[2]) + 1
  y_index <- (J - 1) %% plot_layout[2]
  show_x <- (x_index == plot_layout[1]) || (J + plot_layout[2] > J)
  show_y <- (y_index == 0)

  plot(grd, e(grd),
       xaxt = ifelse(show_x, "s", "n"),
       yaxt = ifelse(show_y, "s", "n"),
       las = 1,
       type = "s",
       cex.axis = cxa,
       cex.main = cxa,
       xlim = xlim,
       ylim = ylim,
       main = plot_title)

  mtext(pvt, side = 3, outer = FALSE, line = -1.5, adj = 0.03, cex = cxt)
  abline(0, 1, col = line_color)

  # ============================================================================
  # Omnibus goodness-of-fit tests
  # ============================================================================

  # Minimum p-value test (adjusted via beta distribution)
  pvalmin <- pbeta(min(pvals), 1, J)

  # Product of p-values test (Fisher's method)
  pprod <- prod(pvals)
  pvalprod <- pprod * sum(sapply(1:J, function(j) {
    ((-1)^(j - 1) / factorial(j - 1)) * log(pprod)^(j - 1)
  }))

  # Add overall plot labels
  mtext("Empirical cumulative distribution function (ECDF)",
        side = 2, outer = TRUE, line = 2.5, cex = 0.9)

  mtext(
    bquote(
      paste("Rosenblatt transformation: ",
            p[min] == .(round(pvalmin, digits)),
            ", ",
            p[prod] == .(round(pvalprod, digits)))
    ),
    side = 1, outer = TRUE, line = 2, cex = 0.9
  )

  if (status) {
    cat("\n=== Goodness-of-Fit Results ===\n")
    cat("Individual test p-values:\n")
    for (i in 1:J) {
      cat(sprintf("  Dimension %d (%s): %.4f\n", i, test_types[i], pvals[i]))
    }
    cat(sprintf("\nMinimum p-value (adjusted): %.4f\n", pvalmin))
    cat(sprintf("Product p-value (adjusted): %.4f\n", pvalprod))
    cat("\nNote: KS = Kolmogorov-Smirnov test (continuous data)\n")
    cat("      CvM = Cramér-von Mises test (discrete/tied data)\n")
    cat("================================\n")
  }

  # Return results
  names(pvals) <- variable_names
  names(test_types) <- variable_names
  names(us) <- variable_names

  invisible(list(
    pvalmin = pvalmin,
    pvalprod = pvalprod,
    individual_pvals = pvals,
    test_types = test_types,
    transformed_data = us
  ))
}

################################################################################
# Helper function for title generation
################################################################################

#' Generate Plot Titles for Goodness-of-Fit Plots
#'
#' Internal helper function to create formatted plot titles for marginal and
#' conditional distribution plots.
#'
#' @param i Integer. Dimension index.
#' @param var_name Character. Optional variable name to include in title.
#'
#' @return An expression object suitable for plot titles.
#' @keywords internal

get_title <- function(i, var_name = NULL) {
  var_label <- if (!is.null(var_name)) {
    paste0(", ", clean_capitalize(var_name))
  } else {
    ""
  }

  titles <- list(
    bquote(paste("ECDF of ", F[1], "(", Y[1], ")", .(var_label))),
    bquote(paste("ECDF of ", F[2], "(", Y[2], " | ", Y[1], ")", .(var_label))),
    bquote(paste("ECDF of ", F[3], "(", Y[3], " | ", Y[2], ", ", Y[1], ")",
                 .(var_label))),
    bquote(paste("ECDF of ", F[4], "(", Y[4], " | ", Y[3], ", ", Y[2], ", ",
                 Y[1], ")", .(var_label))),
    bquote(paste("ECDF of ", F[.(i)], "(", Y[.(i)], " | ", Y[.(i-1)],
                 ", ... , ", Y[1], ")", .(var_label)))
  )

  if (i <= 4) {
    return(titles[[i]])
  } else {
    return(titles[[5]])
  }
}

################################################################################
# Helper function for string cleaning
################################################################################

#' Clean and Capitalize Variable Names
#'
#' Internal helper function to format variable names for display by replacing
#' underscores and periods with spaces and capitalizing each word.
#'
#' @param input_string Character. String to clean and capitalize.
#'
#' @return Character. Formatted string.
#' @keywords internal

clean_capitalize <- function(input_string) {
  clean_string <- gsub("[_.]", " ", input_string)
  capitalized_string <- gsub("\\b(\\w)", "\\U\\1", clean_string, perl = TRUE)
  return(capitalized_string)
}

################################################################################
# Helper function for Cramér-von Mises test
################################################################################

#' Cramér-von Mises Test for Uniform Distribution
#'
#' Internal helper function to test whether data follows a uniform[0,1] distribution
#' using the Cramér-von Mises test. This test is more appropriate than the
#' Kolmogorov-Smirnov test when ties are present in the data (e.g., with discrete
#' or ordered variables).
#'
#' @param x Numeric vector of observations to test.
#'
#' @return Numeric. P-value for the test.
#' @keywords internal
#'
#' @details
#' The Cramér-von Mises statistic is computed as:
#' \deqn{W^2 = \sum_{i=1}^{n} (F_n(x_i) - (2i-1)/(2n))^2 + 1/(12n)}
#' where \eqn{F_n} is the empirical CDF.
#'
#' P-values are approximated using the asymptotic distribution. For small samples
#' (n < 8), results may be unreliable.

cvm_test <- function(x) {
  n <- length(x)

  # Handle edge cases
  if (n < 8) {
    warning("Sample size too small for reliable Cramér-von Mises test (n < 8)")
    return(NA)
  }

  # Sort the data
  x_sorted <- sort(x)

  # Compute CVM statistic
  # W^2 = sum((F_n(x_i) - (2i-1)/(2n))^2) + 1/(12n)
  i <- 1:n
  expected_quantiles <- (2 * i - 1) / (2 * n)

  # CVM statistic
  W2 <- sum((x_sorted - expected_quantiles)^2) + 1 / (12 * n)

  # Adjust for sample size
  W2_adj <- W2 * (1 + 0.5 / n)

  # Approximate p-value using asymptotic distribution
  # Based on Anderson-Darling approximation for CVM statistic
  # Reference: Anderson & Darling (1952), Stephens (1974)

  if (W2_adj < 0.0275) {
    pval <- 1 - exp(-13.953 + 775.5 * W2_adj - 12542.61 * W2_adj^2)
  } else if (W2_adj < 0.051) {
    pval <- 1 - exp(-5.903 + 179.546 * W2_adj - 1515.29 * W2_adj^2)
  } else if (W2_adj < 0.092) {
    pval <- exp(0.886 - 31.62 * W2_adj + 10.897 * W2_adj^2)
  } else if (W2_adj < 0.1591) {
    pval <- exp(1.111 - 34.242 * W2_adj + 12.832 * W2_adj^2)
  } else {
    pval <- exp(-W2_adj * (4.33 + W2_adj * 0.5))
  }

  # Ensure p-value is in [0, 1]
  pval <- max(0, min(1, pval))

  return(pval)
}
