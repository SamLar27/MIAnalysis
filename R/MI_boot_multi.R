#' Bootstrap Analysis for Multiple Imputation Data with Multiple Models
#'
#' @description
#' Performs bootstrap analysis on multiply imputed data to estimate and compare multiple models
#' using the same bootstrap samples. This approach ensures fair comparison by testing all models
#' on identical bootstrap resamples.
#'
#' @param data A data frame containing the imputed datasets.
#' @param outcome_var Character string specifying the name of the outcome variable.
#' @param models_list A named list where each element contains a character vector of predictor
#'   variables for a specific model. Names will be used as model identifiers.
#' @param imp_col Character string specifying the name of the imputation identifier column. Default is ".imp".
#' @param model_type Character string specifying the type of model to fit: "nb" (negative binomial),
#'   "lm" (linear model), or "glm" (Poisson). Default is "nb".
#' @param followup_offset Character string, either "Yes" or "No", indicating whether to include
#'   an offset term for follow-up duration. Default is "No".
#' @param followup_col Character string specifying the name of the follow-up duration column.
#'   Required if followup_offset = "Yes".
#' @param random_intercept Character string, either "Yes" or "No", indicating whether to include
#'   a random intercept term. Default is "No".
#' @param random_intercept_var Character string specifying the grouping variable for random intercepts.
#'   Required if random_intercept = "Yes".
#' @param boot_strata Character string, either "Yes" or "No", indicating whether to perform
#'   stratified bootstrapping. Default is "No".
#' @param strata_var Character string specifying the stratification variable. Required if
#'   boot_strata = "Yes".
#' @param R Integer specifying the number of bootstrap replications. Default is 1000.
#' @param parallel Logical indicating whether to use parallel processing. Default is TRUE.
#' @param n_cores Integer specifying the number of cores to use for parallel processing.
#'   If NULL, uses one less than the available cores. Default is NULL.
#' @param suppress_warnings Logical indicating whether to suppress warnings. Default is TRUE.
#' @param progress_bar Logical indicating whether to display a progress bar. Default is TRUE.
#' @param max_model_time Maximum time in seconds to allow for fitting a single model.
#'   Default is 60.
#' @param compare_models Logical indicating whether to compute formal model comparisons. Default is TRUE.
#'
#' @return A list containing results for each model, plus model comparison statistics if requested.
#'   For each model, results include:
#'   - pooled_bootstrap_matrix: Matrix of pooled bootstrap estimates
#'   - results_table: Parameter estimates with confidence intervals
#'   - exp_results_table: Exponentiated estimates (rate ratios) for count models
#'   - performance_table: Model performance metrics
#'   - model_tracking: Counts of different model fitting outcomes
#'
#'   If compare_models=TRUE, also includes:
#'   - model_comparisons: Pairwise comparisons between all models
#'   - best_model: Summary of which model performed best on different metrics
#'
#' @importFrom dplyr count arrange
#' @importFrom MASS glm.nb
#' @importFrom glmmTMB glmmTMB fixef ranef
#' @importFrom boot boot
#' @importFrom pROC roc
#' @importFrom purrr map2
#' @importFrom pbapply pblapply pboptions
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom stats as.formula BIC AIC coef predict quantile
#' @importFrom utils read.csv
#'
#' @export

MI_boot_multi <- function(
    data,
    outcome_var,
    models_list,
    imp_col = ".imp",
    model_type = "nb",
    followup_offset = "No",
    followup_col = NULL,
    random_intercept = "No",
    random_intercept_var = NULL,
    boot_strata = "No",
    strata_var = NULL,
    R = 1000,
    parallel = TRUE,
    n_cores = NULL,
    suppress_warnings = TRUE,
    progress_bar = TRUE,
    max_model_time = 60,
    compare_models = TRUE
) {
  # Record start time for total runtime
  total_start_time <- Sys.time()

  # Validate models_list
  if (!is.list(models_list) || is.null(names(models_list)) || any(names(models_list) == "")) {
    stop("models_list must be a named list where each element contains predictors for a model")
  }

  # Load required packages
  required_packages <- c("dplyr", "MASS", "glmmTMB", "boot", "pROC", "purrr")
  if (parallel) required_packages <- c(required_packages, "parallel", "foreach", "doParallel", "pbapply")

  # Check if R.utils is available, but don't require it
  has_r_utils <- requireNamespace("R.utils", quietly = TRUE)
  if (!has_r_utils) {
    cat("Note: 'R.utils' package is not installed. Timeout protection will be disabled.\n")
    cat("Consider installing it with: install.packages('R.utils')\n")
  }

  sapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(paste0("Package '", pkg, "' is required."))
  })

  # Create function to expand interaction terms (e.g., "x1*x2" into c("x1", "x2", "x1:x2"))
  expand_interaction_terms <- function(term) {
    if (grepl("\\*", term)) {
      # Split by * to get individual variables
      vars <- unlist(strsplit(term, "\\*"))
      vars <- trimws(vars)

      # Generate all combinations (main effects and interactions)
      result <- c(vars)  # Start with main effects

      if (length(vars) > 1) {
        # Add interaction terms
        for (i in 2:length(vars)) {
          # Get all combinations of i variables
          combs <- combn(vars, i)
          # For each combination, create an interaction term
          for (j in 1:ncol(combs)) {
            result <- c(result, paste(combs[,j], collapse=":"))
          }
        }
      }
      return(result)
    } else {
      # No interaction, return as is
      return(term)
    }
  }

  # Helper function to extract base variable names (non-interactions)
  extract_base_vars <- function(terms) {
    # Split interaction terms by colon
    all_vars <- unique(unlist(strsplit(terms, ":")))
    # Trim whitespace
    all_vars <- trimws(all_vars)
    return(all_vars)
  }

  # Function for timeout handling
  safe_timeout <- function(expr, timeout = 60, onTimeout = "error") {
    if (has_r_utils) {
      R.utils::withTimeout(expr, timeout = timeout, onTimeout = onTimeout)
    } else {
      # If R.utils not available, just evaluate the expression without timeout
      expr
    }
  }

  # Function to compute performance metrics from a model
  compute_performance_metrics <- function(model, data, outcome_var) {
    # If model fitting failed, return NA for all metrics
    if (is.null(model) || inherits(model, "try-error")) {
      perf_names <- c("R2", "AIC", "AICc", "BIC", "BICc", "C_Index", "RMSE", "MAE")
      return(setNames(rep(NA, length(perf_names)), perf_names))
    }

    preds <- try(predict(model, type = "response"), silent = TRUE)
    if (inherits(preds, "try-error")) {
      perf_names <- c("R2", "AIC", "AICc", "BIC", "BICc", "C_Index", "RMSE", "MAE")
      return(setNames(rep(NA, length(perf_names)), perf_names))
    }

    obs <- data[[outcome_var]]
    residuals <- obs - preds
    RMSE <- sqrt(mean(residuals^2, na.rm = TRUE))
    MAE <- mean(abs(residuals), na.rm = TRUE)

    # For R2, use a simpler approach for all model types
    R2 <- tryCatch({
      1 - sum(residuals^2) / sum((obs - mean(obs))^2)
    }, error = function(e) NA)

    # For information criteria, use the model's methods directly
    AIC_val <- tryCatch(AIC(model), error = function(e) NA)
    BIC_val <- tryCatch(BIC(model), error = function(e) NA)

    # For AICc and BICc, we need sample size and parameter count
    n <- nrow(data)
    k <- tryCatch({
      if (inherits(model, "glmmTMB")) {
        length(fixef(model)$cond)
      } else {
        length(coef(model))
      }
    }, error = function(e) NA)

    AICc_val <- if (!is.na(AIC_val) && !is.na(k)) AIC_val + (2 * k^2 + 2 * k) / (n - k - 1) else NA
    BICc_val <- if (!is.na(BIC_val) && !is.na(k)) BIC_val + (k * (k + 1)) / (n - k - 1) else NA

    # For C-index (AUC), use ROC
    cindex <- tryCatch(pROC::roc(obs, preds, quiet = TRUE)$auc, error = function(e) NA)

    return(c(R2 = R2, AIC = AIC_val, AICc = AICc_val, BIC = BIC_val, BICc = BICc_val, C_Index = cindex, RMSE = RMSE, MAE = MAE))
  }

  # Setup for parallel processing
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    cat("🔵 Parallel mode: using", n_cores, "cores.\n")

    parallel::clusterEvalQ(cl, {
      library(glmmTMB)
      library(MASS)
      library(pbapply)
      library(pROC)
      library(dplyr)
      # Load R.utils if available
      if (requireNamespace("R.utils", quietly = TRUE)) library(R.utils)
    })
  }

  # Validate inputs
  if (!imp_col %in% colnames(data)) stop("Missing imputation column.")
  imputations <- unique(data[[imp_col]])
  num_imputations <- length(imputations)

  if (boot_strata == "Yes") {
    if (is.null(strata_var)) stop("Stratification variable must be provided if boot_strata = 'Yes'.")
    if (!strata_var %in% colnames(data)) stop(paste("Stratification variable", strata_var, "not found"))
  }

  if (random_intercept == "Yes") {
    if (is.null(random_intercept_var)) stop("Random intercept variable must be provided if random_intercept = 'Yes'.")
    if (!random_intercept_var %in% colnames(data)) stop(paste("Random intercept variable", random_intercept_var, "not found"))

    # Check if random intercept variable is a factor or convert it
    if (!is.factor(data[[random_intercept_var]])) {
      cat("Converting random intercept variable to factor...\n")
      data[[random_intercept_var]] <- as.factor(data[[random_intercept_var]])
    }

    # Check if there are enough levels
    n_levels <- length(levels(data[[random_intercept_var]]))
    if (n_levels < 5) {
      cat("⚠️ Warning: Random intercept variable has only", n_levels, "levels. Models with random effects typically work better with more levels.\n")
    } else {
      cat("Random intercept variable has", n_levels, "levels.\n")
    }
  }

  # Create a shared progress counter for tracking overall progress if parallel
  if (parallel && progress_bar) {
    progress_file <- tempfile("progress_")
    cat("0", file = progress_file)

    # Function to update progress
    update_progress <- function() {
      old <- as.integer(readLines(progress_file))
      new_val <- old + 1
      cat(as.character(new_val), file = progress_file)
      return(new_val)
    }

    parallel::clusterExport(cl, c("progress_file", "update_progress",
                                  "has_r_utils", "safe_timeout"),
                            envir = environment())
  }

  # Initialize result storage
  all_model_results <- list()

  # Precompute bootstrap indices for each imputation (CRUCIAL for fair model comparisons)
  # This ensures all models are evaluated on EXACTLY the same bootstrap samples
  cat("\n🔹 Precomputing bootstrap indices for fair model comparison...\n")
  bootstrap_indices <- list()

  for (imp_idx in seq_along(imputations)) {
    imp <- imputations[imp_idx]
    imp_data <- data[data[[imp_col]] == imp, ]
    n_obs <- nrow(imp_data)

    # Create bootstrap indices based on stratification option
    if (boot_strata == "Yes") {
      strata_counts <- imp_data %>%
        dplyr::count(!!rlang::sym(strata_var)) %>%
        dplyr::arrange(!!rlang::sym(strata_var))

      bootstrap_indices[[imp_idx]] <- replicate(R, {
        unlist(purrr::map2(
          strata_counts[[strata_var]],
          strata_counts$n,
          function(stratum_name, stratum_n) {
            idx <- which(imp_data[[strata_var]] == stratum_name)
            if (length(idx) == 0) stop(paste0("No subjects found for stratum: ", stratum_name))
            sample(idx, size = stratum_n, replace = TRUE)
          }
        ))
      }, simplify = FALSE)
    } else {
      bootstrap_indices[[imp_idx]] <- replicate(R, sample(seq_len(n_obs), n_obs, replace = TRUE), simplify = FALSE)
    }
  }

  # Process each model using the shared bootstrap samples
  cat("\n🔹 Processing all models using shared bootstrap samples...\n")

  # Main loop: process each model
  for (model_name in names(models_list)) {
    cat("\n🔧 Running bootstrap for model:", model_name, "\n")
    predictor_vars <- models_list[[model_name]]

    # Track model fitting outcomes
    model_counters <- list(
      mixed_success = 0,
      fixed_fallback = 0,
      total_failure = 0
    )

    # Create model tracking file if parallel
    if (parallel) {
      model_tracking_file <- tempfile(paste0("model_tracking_", model_name, "_"))
      cat("", file = model_tracking_file)
    }

    # Function to update model tracking
    update_model_tracking <- function(model_status, imp, boot_idx) {
      if (parallel) {
        # Append to file in atomic operation
        entry <- sprintf("%d,%d,%s\n", imp, boot_idx, model_status)
        cat(entry, file = model_tracking_file, append = TRUE)
      } else {
        # Update in-memory counter
        model_counters[[model_status]] <<- model_counters[[model_status]] + 1
      }
    }

    # For this model, export model tracking file to cluster
    if (parallel) {
      parallel::clusterExport(cl, c("model_tracking_file", "update_model_tracking"),
                              envir = environment())
    }

    # Expand interaction terms in predictor variables
    expanded_terms_list <- lapply(predictor_vars, expand_interaction_terms)
    expanded_predictors <- unique(unlist(expanded_terms_list))

    # Get all base variables
    base_vars <- unique(unlist(lapply(expanded_terms_list, extract_base_vars)))

    # Validate base variable presence in dataset
    missing_vars <- base_vars[!base_vars %in% colnames(data)]
    if (length(missing_vars) > 0) {
      stop(paste("Predictor variables not found in dataset:", paste(missing_vars, collapse=", ")))
    }

    # Notify user of interaction terms detected
    interaction_terms <- expanded_predictors[grepl(":", expanded_predictors)]
    if (length(interaction_terms) > 0) {
      cat("Detected interaction terms:", paste(interaction_terms, collapse = ", "), "\n")
    }

    # Generate formula for this model
    formula_str <- paste(outcome_var, "~", paste(expanded_predictors, collapse = " + "))
    if (followup_offset == "Yes" && !is.null(followup_col)) {
      formula_str <- paste(formula_str, "+ offset(log(", followup_col, "))")
    }
    if (random_intercept == "Yes" && !is.null(random_intercept_var)) {
      formula_str <- paste(formula_str, "+ (1|", random_intercept_var, ")")
    }
    formula_obj <- as.formula(formula_str)
    cat("Formula:", deparse(formula_obj), "\n")

    # Get parameter names from a test model
    get_model_param_names <- function() {
      test_data <- data[data[[imp_col]] == imputations[1], ]

      cat("Running test model to determine parameter names...\n")
      test_model <- tryCatch({
        if (model_type == "nb") {
          if (random_intercept == "Yes") {
            # Use glmmTMB with negative binomial family
            safe_timeout({
              suppressWarnings(glmmTMB(formula_obj, data = test_data, family = nbinom2))
            }, timeout = max_model_time, onTimeout = "error")
          } else {
            suppressWarnings(MASS::glm.nb(formula_obj, data = test_data))
          }
        } else if (model_type == "lm") {
          if (random_intercept == "Yes") {
            safe_timeout({
              suppressWarnings(glmmTMB(formula_obj, data = test_data, family = gaussian))
            }, timeout = max_model_time, onTimeout = "error")
          } else {
            stats::lm(formula_obj, data = test_data)
          }
        } else if (model_type == "glm") {
          if (random_intercept == "Yes") {
            safe_timeout({
              suppressWarnings(glmmTMB(formula_obj, data = test_data, family = poisson))
            }, timeout = max_model_time, onTimeout = "error")
          } else {
            stats::glm(formula_obj, family = poisson, data = test_data)
          }
        }
      }, error = function(e) {
        cat("Error in test model:", conditionMessage(e), "\n")
        cat("Trying alternative approach to determine parameter names...\n")

        # If random intercept model fails, try a simpler fixed-effects model
        if (random_intercept == "Yes") {
          simple_formula <- as.formula(gsub("\\+ \\(1\\|[^)]+\\)", "", formula_str))
          try_simple_model <- try({
            if (model_type == "nb") {
              MASS::glm.nb(simple_formula, data = test_data)
            } else if (model_type == "lm") {
              stats::lm(simple_formula, data = test_data)
            } else if (model_type == "glm") {
              stats::glm(simple_formula, family = poisson, data = test_data)
            }
          }, silent = TRUE)

          if (!inherits(try_simple_model, "try-error")) {
            cat("Using fixed effects model for parameter names.\n")
            return(try_simple_model)
          }
        }

        return(NULL)
      })

      if (is.null(test_model)) {
        cat("Warning: Failed to build test model. Using default parameter names.\n")
        return(c("(Intercept)", expanded_predictors))
      }

      if (random_intercept == "Yes" && inherits(test_model, "glmmTMB")) {
        names(fixef(test_model)$cond)
      } else {
        names(stats::coef(test_model))
      }
    }

    # Get actual parameter names
    actual_param_names <- get_model_param_names()
    cat("Model parameters:", paste(actual_param_names, collapse=", "), "\n")
    n_actual_params <- length(actual_param_names)

    # Define perf_names for consistent output
    perf_names <- c("R2", "AIC", "AICc", "BIC", "BICc", "C_Index", "RMSE", "MAE")

    # Function to run bootstrap for one imputation
    run_one_imputation <- function(imp_idx) {
      imp <- imputations[imp_idx]
      imp_data <- data[data[[imp_col]] == imp, ]
      n_obs <- nrow(imp_data)

      # Track model types for this imputation
      local_model_tracking <- list(
        mixed_success = 0,
        fixed_fallback = 0,
        total_failure = 0
      )

      # Use the precomputed bootstrap indices for this imputation
      boot_results <- lapply(seq_len(R), function(b) {
        # Get bootstrap indices
        boot_indices <- bootstrap_indices[[imp_idx]][[b]]
        d <- imp_data[boot_indices, ]

        model <- tryCatch({
          # Use timeout function to prevent model fitting from taking too long
          safe_timeout({
            if (model_type == "nb") {
              if (random_intercept == "Yes") {
                # Use glmmTMB for negative binomial with random effects
                suppressWarnings(glmmTMB(formula_obj, data = d, family = nbinom2))
              } else {
                suppressWarnings(MASS::glm.nb(formula_obj, data = d))
              }
            } else if (model_type == "lm") {
              if (random_intercept == "Yes") {
                suppressWarnings(glmmTMB(formula_obj, data = d, family = gaussian))
              } else {
                stats::lm(formula_obj, data = d)
              }
            } else if (model_type == "glm") {
              if (random_intercept == "Yes") {
                suppressWarnings(glmmTMB(formula_obj, data = d, family = poisson))
              } else {
                stats::glm(formula_obj, family = poisson, data = d)
              }
            }
          }, timeout = max_model_time, onTimeout = "error")
        }, error = function(e) {
          # If random effect model fails, try a simpler fixed-effects model as fallback
          if (random_intercept == "Yes") {
            simple_formula <- as.formula(gsub("\\+ \\(1\\|[^)]+\\)", "", formula_str))
            try_simple_model <- try({
              if (model_type == "nb") {
                MASS::glm.nb(simple_formula, data = d)
              } else if (model_type == "lm") {
                stats::lm(simple_formula, data = d)
              } else if (model_type == "glm") {
                stats::glm(simple_formula, family = poisson, data = d)
              }
            }, silent = TRUE)

            if (!inherits(try_simple_model, "try-error")) {
              # Track that we used a fixed effects fallback model
              update_model_tracking("fixed_fallback", imp, b)
              if (!parallel) local_model_tracking$fixed_fallback <- local_model_tracking$fixed_fallback + 1

              return(try_simple_model)  # Return the fallback model
            }
          }

          # Track complete failure
          update_model_tracking("total_failure", imp, b)
          if (!parallel) local_model_tracking$total_failure <- local_model_tracking$total_failure + 1

          return(NULL)  # Return NULL if all attempts fail
        })

        # If model successful and we got here, track as mixed effects success
        if (!is.null(model) && random_intercept == "Yes" && inherits(model, "glmmTMB")) {
          update_model_tracking("mixed_success", imp, b)
          if (!parallel) local_model_tracking$mixed_success <- local_model_tracking$mixed_success + 1
        }

        # Extract coefficients or return NAs if model failed
        coefs <- if (!is.null(model)) {
          if (random_intercept == "Yes" && inherits(model, "glmmTMB")) {
            fixef_vals <- fixef(model)$cond
            # Make sure we get all expected coefficients
            result <- setNames(rep(NA, n_actual_params), actual_param_names)
            # Fill in the ones that exist in this model
            result[names(fixef_vals)] <- as.numeric(fixef_vals)
            result
          } else {
            coef_vals <- tryCatch(stats::coef(model), error = function(e) NULL)
            if (is.null(coef_vals)) {
              setNames(rep(NA, n_actual_params), actual_param_names)
            } else {
              # Make sure we get all expected coefficients
              result <- setNames(rep(NA, n_actual_params), actual_param_names)
              # Fill in the ones that exist in this model
              result[names(coef_vals)] <- as.numeric(coef_vals)
              result
            }
          }
        } else {
          setNames(rep(NA, n_actual_params), actual_param_names)
        }

        # Compute performance metrics
        perfs <- compute_performance_metrics(model, d, outcome_var)

        list(coef = coefs, perf = perfs)
      })

      # For parallel mode, update progress counter
      if (parallel && progress_bar) {
        imp_done <- update_progress()
        progress_pct <- 100 * imp_done / (num_imputations * length(models_list))
        cat(sprintf("\rProgress: %.1f%% completed", progress_pct))
      }

      # Extract coefficients ensuring all have the same names
      coef_matrix <- t(sapply(boot_results, function(x) {
        if (is.null(x$coef)) return(rep(NA, n_actual_params))
        as.numeric(x$coef)
      }))

      # Extract performance metrics
      perf_matrix <- t(sapply(boot_results, function(x) {
        if (is.null(x$perf)) return(rep(NA, length(perf_names)))
        as.numeric(x$perf)
      }))

      colnames(coef_matrix) <- actual_param_names
      colnames(perf_matrix) <- perf_names

      # Add the model tracking info to the result
      if (!parallel) {
        attr(coef_matrix, "model_tracking") <- local_model_tracking
      }

      return(list(coef = coef_matrix, perf = perf_matrix))
    }

    # Process all imputations
    if (parallel) {
      parallel::clusterExport(cl, c("data", "outcome_var", "expanded_predictors", "imp_col",
                                    "model_type", "followup_offset", "followup_col",
                                    "random_intercept", "random_intercept_var",
                                    "bootstrap_indices", "formula_str", "formula_obj",
                                    "actual_param_names", "n_actual_params",
                                    "perf_names", "max_model_time"),
                              envir = environment())

      all_boot_results <- parallel::parLapply(cl, seq_along(imputations), run_one_imputation)
    } else {
      all_boot_results <- lapply(seq_along(imputations), run_one_imputation)

      # Aggregate model tracking info from all imputations
      for (i in seq_along(all_boot_results)) {
        if (!is.null(attr(all_boot_results[[i]]$coef, "model_tracking"))) {
          tracking_info <- attr(all_boot_results[[i]]$coef, "model_tracking")
          model_counters$mixed_success <- model_counters$mixed_success + tracking_info$mixed_success
          model_counters$fixed_fallback <- model_counters$fixed_fallback + tracking_info$fixed_fallback
          model_counters$total_failure <- model_counters$total_failure + tracking_info$total_failure
        }
      }
    }

    names(all_boot_results) <- paste0("imp_", imputations)

    # Pool results across bootstrap samples and imputations
    n_bootstrap <- R
    n_metrics <- length(perf_names)
    n_imputations <- length(all_boot_results)

    boot_coef_array <- array(NA, dim = c(n_bootstrap, n_actual_params, n_imputations))
    boot_perf_array <- array(NA, dim = c(n_bootstrap, n_metrics, n_imputations))

    for (i in seq_along(all_boot_results)) {
      boot_coef_array[,,i] <- all_boot_results[[i]]$coef
      boot_perf_array[,,i] <- all_boot_results[[i]]$perf
    }

    pooled_coef_matrix <- apply(boot_coef_array, c(1,2), mean, na.rm = TRUE)
    pooled_perf_matrix <- apply(boot_perf_array, c(1,2), mean, na.rm = TRUE)
    # Create results tables
    results_df <- data.frame(
      Parameter = actual_param_names,
      Estimate = colMeans(pooled_coef_matrix, na.rm = TRUE),
      Lower_CI = apply(pooled_coef_matrix, 2, quantile, probs = 0.025, na.rm = TRUE),
      Upper_CI = apply(pooled_coef_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
    )

    # Create performance metrics summary table
    performance_df <- data.frame(
      Metric = perf_names,
      Estimate = colMeans(pooled_perf_matrix, na.rm = TRUE),
      Lower_CI = apply(pooled_perf_matrix, 2, quantile, probs = 0.025, na.rm = TRUE),
      Upper_CI = apply(pooled_perf_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
    )

    # Create exponentiated estimates table for count models
    if (model_type %in% c("nb", "glm", "poisson")) {
      exp_results_df <- data.frame(
        Parameter = actual_param_names,
        RateRatio = exp(results_df$Estimate),
        Lower_CI = exp(results_df$Lower_CI),
        Upper_CI = exp(results_df$Upper_CI)
      )
    } else {
      exp_results_df <- NULL
    }

    # Compile model tracking results from parallel processing if used
    if (parallel && file.exists(model_tracking_file)) {
      # Read the tracking data
      model_data <- tryCatch({
        read.csv(model_tracking_file, header = FALSE,
                 col.names = c("imp", "boot", "model_type"))
      }, error = function(e) {
        data.frame(imp = numeric(0), boot = numeric(0), model_type = character(0))
      })

      # Compile statistics
      model_counters$mixed_success <- sum(model_data$model_type == "mixed_success")
      model_counters$fixed_fallback <- sum(model_data$model_type == "fixed_fallback")
      model_counters$total_failure <- sum(model_data$model_type == "total_failure")

      # Clean up the tracking file
      file.remove(model_tracking_file)
    }

    # Display model fitting summary
    total_models <- model_counters$mixed_success + model_counters$fixed_fallback + model_counters$total_failure

    cat("\n--- Model Fitting Summary for", model_name, "---\n")
    if (total_models > 0) {
      cat("  - Mixed Effects Models (successful): ", model_counters$mixed_success,
          " (", round(100 * model_counters$mixed_success / total_models, 1), "%)\n", sep = "")
      cat("  - Fixed Effects Models (fallback): ", model_counters$fixed_fallback,
          " (", round(100 * model_counters$fixed_fallback / total_models, 1), "%)\n", sep = "")
      cat("  - Failed Models (NAs): ", model_counters$total_failure,
          " (", round(100 * model_counters$total_failure / total_models, 1), "%)\n", sep = "")
      cat("  - Total Model Fits: ", total_models, "\n", sep = "")
    } else {
      cat("  No model fitting information available.\n")
    }

    # Store results for this model
    all_model_results[[model_name]] <- list(
      model_formula = formula_str,
      predictor_vars = predictor_vars,
      expanded_predictors = expanded_predictors,
      pooled_bootstrap_matrix = pooled_coef_matrix,
      pooled_perf_matrix = pooled_perf_matrix,
      all_bootstrap_results = all_boot_results,
      results_table = results_df,
      exp_results_table = exp_results_df,
      performance_table = performance_df,
      model_type = model_type,
      parameters = actual_param_names,
      bootstrap_samples = R,
      imputations = n_imputations,
      model_tracking = model_counters
    )
  }

  # After all models are processed, perform model comparisons if requested
  if (compare_models && length(models_list) > 1) {
    cat("\n🔍 Performing model comparisons...\n")

    # Initialize container for pairwise comparisons
    model_comparisons <- list()

    # Get pairs of models to compare
    model_names <- names(models_list)
    model_pairs <- combn(model_names, 2, simplify = FALSE)

    for (pair in model_pairs) {
      model1_name <- pair[1]
      model2_name <- pair[2]

      cat("  Comparing:", model1_name, "vs", model2_name, "\n")

      # Extract performance metrics from both models
      model1_perf_matrix <- all_model_results[[model1_name]]$pooled_perf_matrix
      model2_perf_matrix <- all_model_results[[model2_name]]$pooled_perf_matrix

      # Run paired comparison for each metric
      metrics_to_compare <- c("AIC", "BIC", "R2", "C_Index")
      comparison_results <- list()

      for (metric in metrics_to_compare) {
        metric_idx <- which(perf_names == metric)

        if (length(metric_idx) == 0) {
          next
        }

        # Extract metric values for each bootstrap sample
        metric1_values <- model1_perf_matrix[, metric_idx]
        metric2_values <- model2_perf_matrix[, metric_idx]

        # Calculate differences (accounting for direction - for AIC/BIC, lower is better)
        if (metric %in% c("AIC", "BIC", "AICc", "BICc")) {
          # For these metrics, negative differences mean model2 is better
          diff_values <- metric1_values - metric2_values
          p_value <- mean(diff_values > 0, na.rm = TRUE)  # Prob model1 is worse than model2
          favored_model <- ifelse(mean(diff_values, na.rm = TRUE) < 0, model1_name, model2_name)
        } else {
          # For these metrics, positive differences mean model1 is better
          diff_values <- metric1_values - metric2_values
          p_value <- mean(diff_values < 0, na.rm = TRUE)  # Prob model1 is worse than model2
          favored_model <- ifelse(mean(diff_values, na.rm = TRUE) > 0, model1_name, model2_name)
        }

        # Calculate mean difference and confidence interval
        mean_diff <- mean(diff_values, na.rm = TRUE)
        ci_lower <- quantile(diff_values, 0.025, na.rm = TRUE)
        ci_upper <- quantile(diff_values, 0.975, na.rm = TRUE)

        # Store results
        comparison_results[[metric]] <- list(
          metric = metric,
          model1_name = model1_name,
          model2_name = model2_name,
          mean_diff = mean_diff,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          p_value = p_value,
          favored_model = favored_model,
          significant = (p_value < 0.05),
          diff_values = diff_values
        )
      }

      # Create summary table for this pair
      pair_name <- paste(model1_name, "vs", model2_name)
      metrics <- names(comparison_results)

      if (length(metrics) > 0) {
        summary_df <- data.frame(
          Metric = metrics,
          Mean_Difference = sapply(comparison_results, function(x) x$mean_diff),
          CI_Lower = sapply(comparison_results, function(x) x$ci_lower),
          CI_Upper = sapply(comparison_results, function(x) x$ci_upper),
          P_Value = sapply(comparison_results, function(x) x$p_value),
          Favored_Model = sapply(comparison_results, function(x) x$favored_model),
          Significant = sapply(comparison_results, function(x) x$significant)
        )

        # Print summary
        cat("\n  Summary for", pair_name, ":\n")
        print(summary_df)

        # Store results
        model_comparisons[[pair_name]] <- list(
          summary_table = summary_df,
          detailed_results = comparison_results
        )
      }
    }

    # Identify best model for each metric
    best_model_summary <- list()

    for (metric in c("AIC", "BIC", "R2", "C_Index")) {
      all_values <- sapply(model_names, function(model_name) {
        perf_table <- all_model_results[[model_name]]$performance_table
        idx <- which(perf_table$Metric == metric)
        if (length(idx) > 0) perf_table$Estimate[idx] else NA
      })

      if (all(is.na(all_values))) {
        next
      }

      # Determine best model (lower is better for AIC/BIC, higher is better for others)
      if (metric %in% c("AIC", "BIC", "AICc", "BICc")) {
        best_idx <- which.min(all_values)
        best_model <- model_names[best_idx]
      } else {
        best_idx <- which.max(all_values)
        best_model <- model_names[best_idx]
      }

      # Compare best model to each other model
      comparison_p_values <- numeric(length(model_names))
      names(comparison_p_values) <- model_names

      for (i in seq_along(model_names)) {
        if (i != best_idx) {
          # Find the comparison with this model
          pair_name1 <- paste(model_names[best_idx], "vs", model_names[i])
          pair_name2 <- paste(model_names[i], "vs", model_names[best_idx])

          if (pair_name1 %in% names(model_comparisons)) {
            p_value <- model_comparisons[[pair_name1]]$detailed_results[[metric]]$p_value
            comparison_p_values[i] <- p_value
          } else if (pair_name2 %in% names(model_comparisons)) {
            p_value <- model_comparisons[[pair_name2]]$detailed_results[[metric]]$p_value
            comparison_p_values[i] <- p_value
          }
        } else {
          comparison_p_values[i] <- NA  # No comparison with itself
        }
      }

      # Store best model info
      best_model_summary[[metric]] <- list(
        best_model = best_model,
        metric_values = all_values,
        comparison_p_values = comparison_p_values,
        significantly_better = comparison_p_values < 0.05
      )
    }

    # Print best model summary
    cat("\n--- Best Model Summary ---\n")
    for (metric in names(best_model_summary)) {
      cat(metric, ": Best model is", best_model_summary[[metric]]$best_model, "\n")
      cat("  Values for all models: \n")
      for (i in seq_along(model_names)) {
        value <- best_model_summary[[metric]]$metric_values[i]
        p_value <- best_model_summary[[metric]]$comparison_p_values[i]

        if (is.na(p_value)) {
          cat("    ", model_names[i], ": ", format(value, digits = 4), " (best)\n", sep = "")
        } else if (p_value < 0.05) {
          cat("    ", model_names[i], ": ", format(value, digits = 4),
              " (p = ", format(p_value, digits = 3), ", significantly different)\n", sep = "")
        } else {
          cat("    ", model_names[i], ": ", format(value, digits = 4),
              " (p = ", format(p_value, digits = 3), ", not significantly different)\n", sep = "")
        }
      }
      cat("\n")
    }

    # Add comparison results to output
    all_model_results$model_comparisons <- model_comparisons
    all_model_results$best_model_summary <- best_model_summary
  }

  # Calculate total execution time
  total_time <- difftime(Sys.time(), total_start_time, units = "mins")
  cat("\n✅ Analysis completed in", round(as.numeric(total_time), 2), "minutes\n")

  # Clean up parallel cluster if used
  if (parallel) {
    parallel::stopCluster(cl)

    # Clean up progress file if exists
    if (progress_bar && file.exists(progress_file)) {
      file.remove(progress_file)
    }
  }

  # Return all results
  return(all_model_results)
}

