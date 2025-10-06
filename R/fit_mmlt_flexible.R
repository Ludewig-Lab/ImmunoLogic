#' Flexible Multivariate Maximum Likelihood Transformation Model (MMLT)
#'
#' Fits a flexible multivariate maximum likelihood transformation (MMLT) model
#' using the `tram` package. Handles automatic variable classification into
#' continuous or discrete types, missing data imputation, and model post-processing.
#'
#' @param dataframe Data frame containing the variables to be modeled.
#' @param impute Logical. Whether to impute missing data using `mice`. Default is `FALSE`.
#' @param seed Random seed for reproducibility.
#' @param marginal_formula Formula for the marginal models.
#' @param multivariate_formula Formula for the multivariate dependency.
#' @param continuous_vars Character vector of continuous variable names.
#' @param discrete_vars Character vector of discrete variable names.
#' @param exclude_vars Variables to exclude from modeling.
#' @param continuous_formula Type of continuous response handling: `"default"`, `"interval"`, or `"ordered"`.
#' @param order Polynomial order for the Box–Cox basis.
#' @param automatic Logical. If `TRUE`, classify variables automatically.
#' @param auto_threshold Integer. Threshold for automatic discrete classification (default = 5).
#' @author Roman Stadler
#'
#' @return A list containing:
#'   \item{object}{The fitted MMLT model.}
#'   \item{coef_corr}{Estimated correlation matrix.}
#'   \item{coef_lambda}{Model coefficients.}
#'   \item{standard_errors_lambda}{Standard errors of coefficients.}
#'   \item{p_values}{P-values for all coefficients.}
#'   \item{p_matrix}{Matrix of pairwise p-values.}
#'   \item{AIC, BIC, logLik}{Model fit statistics.}
#'   \item{automatic_classification}{Information about auto-classified variables (if used).}
#'
#' @details
#' This function wraps a multi-step modeling pipeline:
#' data cleaning, optional imputation, automatic variable classification,
#' marginal model fitting for each variable, and final multivariate model fitting.
#' It outputs a structured list with results and diagnostics.
#'
#' @importFrom tram BoxCox mmlt mltoptim
#' @importFrom qrng sobol
#' @importFrom stats logLik coef pnorm
#' @export
#'
#' @examples
#' \dontrun{
#' library(tram)
#' data <- data.frame(
#'   x1 = rnorm(100),
#'   x2 = rnorm(100),
#'   x3 = sample(1:5, 100, TRUE)
#' )
#' result <- fit_mmlt_flexible(data, automatic = TRUE)
#' }



fit_mmlt_flexible <- function(dataframe,
                              impute = F,
                              seed=490,
                              marginal_formula="~ 1",
                              multivariate_formula="~ 1",
                              continuous_vars = c(),
                              discrete_vars = c(),
                              exclude_vars = c(), # removes a variable for the responses
                              continuous_formula = "default", # default, interval, ordered
                              order=6,
                              automatic = T, # new parameter
                              auto_threshold = 5 # threshold for automatic classification
){
  require(tram)
  require(qrng)
  set.seed(seed)

  cat("=== MMLT Flexible Model Fitting ===\n")
  cat("Dataset dimensions:", nrow(dataframe), "observations,", ncol(dataframe), "variables\n")
  cat("Seed:", seed, "\n")

  # Remove columns that are entirely NA
  all_na_cols <- names(dataframe)[colSums(is.na(dataframe)) == nrow(dataframe)]
  if (length(all_na_cols) > 0) {
    cat("Removing columns with only NA values:", paste(all_na_cols, collapse = ", "), "\n")
  }
  dataframe <- dataframe[, colSums(is.na(dataframe)) < nrow(dataframe), drop = FALSE]
  cat("After removing all-NA columns:", dim(dataframe), "\n")

  # Remove rows that are entirely NA
  all_na_rows <- which(rowSums(is.na(dataframe)) == ncol(dataframe))
  if (length(all_na_rows) > 0) {
    cat("Removing rows with only NA values:", paste(all_na_rows, collapse = ", "), "\n")
  }
  dataframe <- dataframe[rowSums(is.na(dataframe)) < ncol(dataframe), , drop = FALSE]
  cat("After removing all-NA rows:", dim(dataframe), "\n")


  # ensure data compatibility
  dataframe <- as.data.frame(dataframe)


  # Check if all continuous and discrete variables are in the dataframe
  all_vars <- c(continuous_vars, discrete_vars)
  missing_vars <- setdiff(all_vars, names(dataframe))
  if (length(missing_vars) > 0) {
    warning(paste("Missing variables in dataframe:", paste(missing_vars, collapse = ", "), ", omitting"))
  }

  cat("\n--- Variable Processing ---\n")
  cat("User-specified continuous variables:", ifelse(length(continuous_vars) > 0, paste(continuous_vars, collapse = ", "), "none"), "\n")
  cat("User-specified discrete variables:", ifelse(length(discrete_vars) > 0, paste(discrete_vars, collapse = ", "), "none"), "\n")
  cat("Variables to exclude:", ifelse(length(exclude_vars) > 0, paste(exclude_vars, collapse = ", "), "none"), "\n")

  # make sure to exclude the marginal covariate from the responses, a regression on itself causes errors
  exclusion_variable <- ifelse(trimws(sub("~", "", multivariate_formula)) != "1", trimws(sub("~", "", multivariate_formula)), "")
  nm <- setdiff(names(dataframe), exclusion_variable) # if the argument is 1, the exclusion will be "", which is ignored on this line
  nm <- setdiff(nm, exclude_vars)

  if (exclusion_variable != "") {
    cat("Excluding multivariate covariate from responses:", exclusion_variable, "\n")
  }
  cat("Variables to model:", paste(nm, collapse = ", "), "\n")



  cat("\n--- Missing Values Summary ---\n")
  na_counts <- sapply(dataframe[, nm, drop = FALSE], function(x) sum(is.na(x)))
  na_percent <- round(100 * na_counts / nrow(dataframe), 2)
  na_summary <- data.frame(Variable = names(na_counts),
                           NAs = na_counts,
                           Percent = na_percent)
  print(na_summary)

  if(max(na_summary$Percent)>10 & !impute){
    warning("Data contains large amount of missing entries. Consider setting impute = TRUE")
  }

  if(impute){
    if (!requireNamespace("mice", quietly = TRUE) ) {
      stop("Please install 'mice' for the imputation.")
    }
    if (verbose) cat("Performing MICE imputation...\n")
    library(mice)
    imp.data <- mice(data = dataframe, m = 50, maxit = 10, seed = 12345, printFlag = F)
    dataframe <- complete(imp.data)
    if (verbose) cat("Imputation completed\n")
  }


  # Automatic variable classification if enabled
  if (automatic) {
    cat("\n--- Automatic Variable Classification ---\n")
    cat("Threshold for discrete classification:", auto_threshold, "unique values\n")

    auto_continuous_vars <- c()
    auto_discrete_vars <- c()

    for (var in nm) {
      if (!(var %in% exclude_vars)) {
        # Count unique non-missing values
        unique_vals <- length(unique(dataframe[[var]][!is.na(dataframe[[var]])]))

        if (unique_vals <= auto_threshold) {
          auto_discrete_vars <- c(auto_discrete_vars, var)
          cat("  ", var, ":", unique_vals, "unique values -> DISCRETE\n")
        } else {
          auto_continuous_vars <- c(auto_continuous_vars, var)
          cat("  ", var, ":", unique_vals, "unique values -> CONTINUOUS\n")
        }
      }
    }

    # Combine automatic classification with user-specified variables
    # User-specified variables take precedence
    final_continuous_vars <- unique(c(continuous_vars, auto_continuous_vars))
    final_discrete_vars <- unique(c(discrete_vars, auto_discrete_vars))

    # Remove conflicts (user specification overrides automatic)
    final_continuous_vars <- setdiff(final_continuous_vars, discrete_vars)
    final_discrete_vars <- setdiff(final_discrete_vars, continuous_vars)

    # Print final classification summary
    cat("\nFinal variable classification:\n")
    if (length(final_discrete_vars) > 0) {
      cat("  Discrete variables:", paste(final_discrete_vars, collapse = ", "), "\n")
    }
    if (length(final_continuous_vars) > 0) {
      cat("  Continuous variables:", paste(final_continuous_vars, collapse = ", "), "\n")
    }

    # Report any overrides
    overridden_vars <- intersect(c(continuous_vars, discrete_vars), c(auto_continuous_vars, auto_discrete_vars))
    if (length(overridden_vars) > 0) {
      cat("  User specifications override automatic classification for:", paste(overridden_vars, collapse = ", "), "\n")
    }

  } else {
    # Use original user-specified variables
    final_continuous_vars <- continuous_vars
    final_discrete_vars <- discrete_vars
    cat("\nUsing user-specified variable classifications (automatic = FALSE)\n")
  }

  # Specify Optimizers
  op <- mltoptim()[1:3] # Order: auglag, spg, nloptr

  cat("\n--- Marginal Model Fitting ---\n")
  cat("Marginal formula:", marginal_formula, "\n")
  cat("Continuous formula method:", continuous_formula, "\n")
  cat("Polynomial order:", order, "\n")
  cat("Optimizers:", paste(names(op), collapse = ", "), "\n\n")

  m <- list() # Initialize m as a list of marginal models for the joint call below
  for (i in seq_along(nm)) {
    name <- nm[i]

    if(!(name %in% exclude_vars)){
      cat("Fitting marginal model for:", name)

      if (name %in% final_discrete_vars) {
        cat(" (DISCRETE)\n")

        formula_string <- paste("R(", nm[i], ", as.R.interval= TRUE) ", marginal_formula)
        formula <- as.formula(formula_string)

        m[[i]] <- BoxCox(formula, data = dataframe, order=order)

      } else {
        if(!(name %in% final_continuous_vars)){
          warning(paste(name,"not categorized, including it as continuous"))
        }

        cat(" (CONTINUOUS - ", continuous_formula, ")\n", sep="")

        if(continuous_formula=="interval") {
          # The likelihood is  calculated over intervals
          # trans.fun. is non parametric likelihood

          formula_string <- paste("R(", nm[i], ", as.R.interval= TRUE) ", marginal_formula)
          formula <- as.formula(formula_string)

        } else if(continuous_formula=="ordered") {
          formula_string <- paste("R(", nm[i], ", as.R.ordered= TRUE) ", marginal_formula)
          formula <- as.formula(formula_string)

        } else {
          # log density
          # The likelihood is the continuous density on Bernstein Polynomial

          formula_string <- paste("R(", nm[i], ") ", marginal_formula)
          formula <- as.formula(formula_string)
        }

        m[[i]] <- BoxCox(formula, data = dataframe, order=order) # add the formula from above

      }
    }

  }


  m$data <- dataframe
  m$formula <- as.formula(multivariate_formula)
  m$optim <- op

  cat("\n--- Joint Model Fitting ---\n")
  cat("Multivariate formula:", multivariate_formula, "\n")
  cat("Fitting joint MMLT model with", length(nm), "variables...\n")

  # Always Joint here:
  m$domargins <- TRUE
  m$theta <- coef(m) # Use hot start
  mm <- do.call("mmlt", m)

  cat("Joint model fitting completed successfully!\n")

  log_likelihood <- logLik(mm)


  n_par_total <- length(coef(mm))
  n_par_marginal <- sum(sapply(coef(mm, type = "marginal"), length))
  n_par_multivariate <- n_par_total - n_par_marginal
  n_obs <- nrow(dataframe)

  AIC <- 2*n_par_total - 2*as.numeric(log_likelihood)
  BIC <- n_par_total*log(n_obs) - 2*as.numeric(log_likelihood)

  cat("\n--- Model Summary ---\n")
  cat("Log-likelihood:", round(as.numeric(log_likelihood), 4), "\n")
  cat("Number of parameters:", n_par_total, "(", n_par_marginal, "marginal +", n_par_multivariate, "multivariate)\n")
  cat("AIC:", round(AIC, 4), "\n")
  cat("BIC:", round(BIC, 4), "\n")

  coef_corr <- tryCatch( # If we have stratified correlation matrices, this will fail, only useful for combined correlations
    expr = {
      as.array(coef(mm, type="Corr"))[,,1]
    },error = function(e){  })

  cat("\n--- Post-Processing ---\n")
  cat("Extracting variances and computing p-values...\n")

  coef_lambda <- coef(mm)
  standard_errors_lambda <- sqrt(diag(vcov(mm)))
  p_values <- 2 * pnorm(-abs( coef_lambda/ standard_errors_lambda))


  clean_names <- function(x, pretty = FALSE) {
    x <- as.character(x)

    # remove any "R." tokens
    x <- gsub("R\\.", "", x, perl = TRUE)

    # simplify Bs terms: Bs1(...whatever...) -> Bs1
    x <- gsub("Bs([0-9]+)\\([^)]*\\)", "Bs\\1", x, perl = TRUE)

    # intercept
    x <- gsub("\\(Intercept\\)", "Intercept", x, perl = TRUE)

    # remove "as.R.interval...TRUE" and variants
    x <- gsub("\\.?as\\.R\\.interval\\.*TRUE\\.*", "", x, perl = TRUE)
    x <- gsub("\\.?as\\.interval\\.*TRUE\\.*", "", x, perl = TRUE)
    x <- gsub("\\(as\\.interval\\.TRUE\\)", "", x, perl = TRUE)

    # fix duplicated Species.Species.*
    x <- gsub("^Species\\.Species\\.", "Species.", x, perl = TRUE)

    # collapse multiple dots, remove trailing/leading
    x <- gsub("\\.\\.+", ".", x, perl = TRUE)
    x <- gsub("^\\.|\\.$", "", x, perl = TRUE)

    if (pretty) {
      # pretty formatting: ".BsN" -> " [BsN]" and dots -> spaces
      x <- gsub("\\.Bs([0-9]+)$", " [Bs\\1]", x, perl = TRUE)
      x <- gsub("\\.", " ", x, perl = TRUE)
      x <- gsub("\\s+", " ", x, perl = TRUE)
      x <- trimws(x)
    }

    x
  }


  rownames(coef_corr) <- clean_names(rownames(coef_corr))
  colnames(coef_corr) <- clean_names(colnames(coef_corr))

  # apply to coef_lambda and p_values
  names(coef_lambda) <- clean_names(names(coef_lambda))
  names(p_values)    <- clean_names(names(p_values))


  # create the p value matrix too:
  p_values_marginal <- p_values[n_par_marginal+1:n_par_multivariate]
  var_names <- rownames(coef_corr)
  n_vars <- length(var_names)
  # Initialize a matrix for p-values
  p_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars,
                     dimnames = list(var_names, var_names))
  # Fill diagonal with 0 (or NA, depending on your preference)
  diag(p_matrix) <- 0  # or use NA if you prefer
  # Fill the upper triangle and mirror to lower triangle
  for (i in 1:(n_vars-1)) {
    for (j in (i+1):n_vars) {
      # Create the name pattern that matches your p_values_marginal names
      var1 <- var_names[i]
      var2 <- var_names[j]
      # Look for the p-value name in the format: R.var2..R.var1..(Intercept)
      p_name <- paste0(var2, ".", var1, ".Intercept")
      # print(p_name)

      if (p_name %in% names(p_values_marginal)) {
        p_val <- p_values_marginal[p_name]
        p_matrix[i, j] <- p_val
        p_matrix[j, i] <- p_val  # Mirror to lower triangle
      }
    }
  }

  cat("Creating correlation p-value matrix...\n")
  cat("Model fitting and post-processing completed!\n")
  cat("=====================================\n\n")




  return(
    list(
      object = mm,
      coef_corr =  coef_corr,
      coef_lambda = coef_lambda,
      standard_errors_lambda = standard_errors_lambda,
      p_values = p_values,
      p_matrix = p_matrix,
      logLik = log_likelihood,
      AIC = AIC,
      BIC = BIC,
      n_par_total = n_par_total,
      n_par_marginal = n_par_marginal,
      n_par_multivariate = n_par_multivariate,
      n_obs = n_obs,
      n_var = length(nm),
      # New return values for automatic classification info
      automatic_classification = if(automatic) list(
        continuous_vars = final_continuous_vars,
        discrete_vars = final_discrete_vars,
        threshold = auto_threshold
      ) else NULL
    )
  )
}
