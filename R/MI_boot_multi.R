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
    reference_model = NULL,  # Parameter to specify reference model
    bootstrap_strategy = "separate",  # New parameter: "separate" or "shared"
    R = 1000,
    parallel = TRUE,
    n_cores = NULL,
    suppress_warnings = TRUE,
    progress_bar = TRUE,
    max_model_time = 60,
    compare_models = TRUE,
    calculate_delta_metrics = TRUE  # Parameter to calculate delta metrics per sample
) {
  # Record start time for total runtime
  total_start_time <- Sys.time()

  # Validate bootstrap_strategy
  if (!bootstrap_strategy %in% c("separate", "shared")) {
    stop("bootstrap_strategy must be one of: 'separate' or 'shared'")
  }

  # Validate models_list
  if (!is.list(models_list) || is.null(names(models_list)) || any(names(models_list) == "")) {
    stop("models_list must be a named list where each element contains predictors for a model")
  }

  # Set reference model if not specified (use first model by default)
  if (is.null(reference_model)) {
    reference_model <- names(models_list)[1]
    cat("No reference model specified. Using first model (", reference_model, ") as reference.\n")
  } else if (!reference_model %in% names(models_list)) {
    stop("Reference model '", reference_model, "' not found in models_list.")
  }

  # Load required packages
  required_packages <- c("dplyr", "MASS", "glmmTMB", "boot", "pROC", "purrr")
  if (parallel) required_packages <- c(required_packages, "parallel", "foreach", "doParallel", "pbapply")

  # Check if spline terms are present in any model and load rms if needed
  has_spline_terms <- FALSE
  for (model_name in names(models_list)) {
    if (any(grepl("^rcs\\(", models_list[[model_name]]))) {
      has_spline_terms <- TRUE
      break
    }
  }

  if (has_spline_terms) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("Package 'rms' is needed for restricted cubic splines. Please install it.")
    }
    required_packages <- c(required_packages, "rms")
  }

  # Check if R.utils is available, but don't require it
  has_r_utils <- requireNamespace("R.utils", quietly = TRUE)
  if (!has_r_utils) {
    cat("Note: 'R.utils' package is not installed. Timeout protection will be disabled.\n")
    cat("Consider installing it with: install.packages('R.utils')\n")
  }

  sapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(paste0("Package '", pkg, "' is required."))
  })

  # Helper function to expand interaction terms (e.g., "x1*x2" into c("x1", "x2", "x1:x2"))
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

  # Helper function to extract base variable names from terms including rcs() functions
  extract_base_vars <- function(terms) {
    result <- character(0)

    for (term in terms) {
      if (grepl("^rcs\\(", term)) {
        # Extract the variable name from rcs() function
        var_name <- gsub("^rcs\\(([^,]+),.*$", "\\1", term)
        result <- c(result, trimws(var_name))
      } else if (grepl(":", term)) {
        # Handle interaction terms
        vars <- unlist(strsplit(term, ":"))
        result <- c(result, trimws(vars))
      } else if (grepl("\\*", term)) {
        # Handle interaction terms with *
        vars <- unlist(strsplit(term, "\\*"))
        result <- c(result, trimws(vars))
      } else {
        # Regular variable
        result <- c(result, trimws(term))
      }
    }

    return(unique(result))
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

      # Load rms if spline terms are present
      if (requireNamespace("rms", quietly = TRUE)) library(rms)

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

  # Precompute bootstrap indices
  cat("\n🔹 Precomputing bootstrap indices...\n")

  # For different bootstrapping strategies
  if (bootstrap_strategy == "separate") {
    # Original approach: different bootstrap samples for each imputation
    cat("Using 'separate' bootstrap strategy: different samples for each imputation, consistent across models\n")

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
  } else {
    # New approach: same bootstrap samples across all imputations
    cat("Using 'shared' bootstrap strategy: same bootstrap samples across all imputations\n")

    # Get full data dimensions
    full_data_ids <- unique(data[!duplicated(data[, setdiff(names(data), imp_col)]), ])
    n_original <- nrow(full_data_ids)

    # Create bootstrap indices for original dataset
    if (boot_strata == "Yes") {
      strata_counts <- full_data_ids %>%
        dplyr::count(!!rlang::sym(strata_var)) %>%
        dplyr::arrange(!!rlang::sym(strata_var))

      # Generate R bootstrap samples with stratification
      master_bootstrap_indices <- replicate(R, {
        unlist(purrr::map2(
          strata_counts[[strata_var]],
          strata_counts$n,
          function(stratum_name, stratum_n) {
            idx <- which(full_data_ids[[strata_var]] == stratum_name)
            if (length(idx) == 0) stop(paste0("No subjects found for stratum: ", stratum_name))
            sample(idx, size = stratum_n, replace = TRUE)
          }
        ))
      }, simplify = FALSE)
    } else {
      # Generate R bootstrap samples without stratification
      master_bootstrap_indices <- replicate(R, sample(seq_len(n_original), n_original, replace = TRUE), simplify = FALSE)
    }

    # Map these indices to each imputation
    bootstrap_indices <- list()

    for (imp_idx in seq_along(imputations)) {
      imp <- imputations[imp_idx]
      imp_data <- data[data[[imp_col]] == imp, ]

      # Map the master bootstrap indices to this imputation's row indices
      bootstrap_indices[[imp_idx]] <- lapply(master_bootstrap_indices, function(boot_idx) {
        # Get bootstrapped original data rows
        boot_original_rows <- full_data_ids[boot_idx, ]

        # Find matching rows in this imputation
        # For each bootstrapped row, find its corresponding row in the current imputation
        # This is a simplification - you'll need to implement logic to match rows across imputations
        # This might involve using row IDs or other identifying columns

        # For demonstration, assuming rows can be matched by their position:
        match_indices <- seq_len(nrow(imp_data))
        match_indices[match(boot_idx, seq_len(n_original))]
      })
    }
  }

  # Process each model using the shared bootstrap samples
  cat("\n🔹 Processing all models using", bootstrap_strategy, "bootstrap samples...\n")

  # Define perf_names for consistent output
  perf_names <- c("R2", "AIC", "AICc", "BIC", "BICc", "C_Index", "RMSE", "MAE")

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

    # Extract all base variables (including from rcs terms)
    base_vars <- extract_base_vars(predictor_vars)

    # Validate that all base variables exist in the dataset
    missing_vars <- base_vars[!base_vars %in% colnames(data)]
    if (length(missing_vars) > 0) {
      stop(paste("Predictor variables not found in dataset:", paste(missing_vars, collapse=", ")))
    }

    # Notify user of interaction terms detected
    interaction_terms <- expanded_predictors[grepl(":", expanded_predictors)]
    if (length(interaction_terms) > 0) {
      cat("Detected interaction terms:", paste(interaction_terms, collapse=", "), "\n")
    }

    # Notify user of spline terms detected
    spline_terms <- predictor_vars[grepl("^rcs\\(", predictor_vars)]
    if (length(spline_terms) > 0) {
      cat("Detected restricted cubic spline terms:", paste(spline_terms, collapse=", "), "\n")
    }

    # Generate formula for this model
    # First, determine if we need to modify predictor_vars for rcs terms
    formula_parts <- predictor_vars

    # Create formula string
    formula_str <- paste(outcome_var, "~", paste(formula_parts, collapse = " + "))
    if (followup_offset == "Yes" && !is.null(followup_col)) {
      formula_str <- paste(formula_str, "+ offset(log(", followup_col, "))")
    }
    if (random_intercept == "Yes" && !is.null(random_intercept_var)) {
      formula_str <- paste(formula_str, "+ (1|", random_intercept_var, ")")
    }
    formula_obj <- as.formula(formula_str)

    cat("Formula:", deparse(formula_obj, width.cutoff = 500), "\n")

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

    # Function to run bootstrap for one imputation
    run_one_imputation <- function(imp_idx) {
      imp <- imputations[imp_idx]
      imp_data <- data[data[[imp_col]] == imp, ]
      n_obs <- nrow(imp_data)

      # Ensure rms package is loaded for spline terms
      if (length(spline_terms) > 0) {
        suppressMessages(require(rms))
      }

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
                                    "actual_param_names", "n_actual_params", "bootstrap_strategy",
                                    "perf_names", "max_model_time", "spline_terms"),
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
    # Pool results differently based on bootstrap_strategy
    n_bootstrap <- R
    n_metrics <- length(perf_names)
    n_imputations <- length(all_boot_results)

    if (bootstrap_strategy == "separate") {
      # Original approach: Pool across imputations for each bootstrap
      boot_coef_array <- array(NA, dim = c(n_bootstrap, n_actual_params, n_imputations))
      boot_perf_array <- array(NA, dim = c(n_bootstrap, n_metrics, n_imputations))

      for (i in seq_along(all_boot_results)) {
        boot_coef_array[,,i] <- all_boot_results[[i]]$coef
        boot_perf_array[,,i] <- all_boot_results[[i]]$perf
      }

      # Pool across imputations for each bootstrap sample
      pooled_coef_matrix <- apply(boot_coef_array, c(1,2), mean, na.rm = TRUE)
      pooled_perf_matrix <- apply(boot_perf_array, c(1,2), mean, na.rm = TRUE)
    } else {
      # Shared bootstrap approach: First, organize by bootstrap
      # Each bootstrap sample has results from all imputations
      # We'll use Rubin's rules to pool across imputations for each bootstrap

      # Initialize matrices to store pooled results
      pooled_coef_matrix <- matrix(NA, nrow = n_bootstrap, ncol = n_actual_params)
      colnames(pooled_coef_matrix) <- actual_param_names

      pooled_perf_matrix <- matrix(NA, nrow = n_bootstrap, ncol = n_metrics)
      colnames(pooled_perf_matrix) <- perf_names

      # For each bootstrap sample
      for (b in 1:n_bootstrap) {
        # Collect coefficients across imputations for this bootstrap sample
        coefs_across_imps <- matrix(NA, nrow = n_imputations, ncol = n_actual_params)
        for (i in 1:n_imputations) {
          coefs_across_imps[i,] <- all_boot_results[[i]]$coef[b,]
        }

        # Apply Rubin's rules to pool coefficients
        pooled_coef_matrix[b,] <- colMeans(coefs_across_imps, na.rm = TRUE)

        # Collect performance metrics across imputations for this bootstrap sample
        perf_across_imps <- matrix(NA, nrow = n_imputations, ncol = n_metrics)
        for (i in 1:n_imputations) {
          perf_across_imps[i,] <- all_boot_results[[i]]$perf[b,]
        }

        # Pool performance metrics
        pooled_perf_matrix[b,] <- colMeans(perf_across_imps, na.rm = TRUE)
      }
    }

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
      model_tracking = model_counters,
      bootstrap_strategy = bootstrap_strategy
    )
  }  # End of the model-specific processing loop

  # After all models are processed, calculate delta metrics if requested
  if (calculate_delta_metrics && length(models_list) > 1) {
    cat("\n🔹 Calculating delta metrics for each bootstrap sample...\n")

    # Get metrics names (AIC, BIC, etc.)
    perf_names <- c("R2", "AIC", "AICc", "BIC", "BICc", "C_Index", "RMSE", "MAE")

    # Initialize storage for delta metrics
    delta_metrics <- list()

    # Process each model (except reference)
    for (model_name in setdiff(names(models_list), reference_model)) {
      cat("  Calculating delta metrics for model:", model_name, "vs", reference_model, "\n")

      delta_metrics[[model_name]] <- list()

      # For each imputation
      for (imp_idx in seq_along(imputations)) {
        imp <- imputations[imp_idx]
        imp_name <- paste0("imp_", imp)

        # Get performance metrics for this model and reference model
        model_perf <- all_model_results[[model_name]]$all_bootstrap_results[[imp_name]]$perf
        ref_perf <- all_model_results[[reference_model]]$all_bootstrap_results[[imp_name]]$perf

        # Calculate delta metrics (model - reference)
        delta_perf <- model_perf - ref_perf

        # Store delta metrics for this imputation
        delta_metrics[[model_name]][[imp_name]] <- delta_perf
      }

      # Calculate summary statistics for delta metrics
      delta_summary <- list()

      # For each metric
      for (metric_idx in seq_along(perf_names)) {
        metric_name <- perf_names[metric_idx]

        # Collect all delta values for this metric across all imputations and bootstraps
        all_delta_values <- numeric(0)

        for (imp_name in names(delta_metrics[[model_name]])) {
          all_delta_values <- c(all_delta_values,
                                delta_metrics[[model_name]][[imp_name]][, metric_idx])
        }

        # Calculate summary statistics
        delta_summary[[metric_name]] <- list(
          mean = mean(all_delta_values, na.rm = TRUE),
          median = median(all_delta_values, na.rm = TRUE),
          sd = sd(all_delta_values, na.rm = TRUE),
          ci_lower = quantile(all_delta_values, 0.025, na.rm = TRUE),
          ci_upper = quantile(all_delta_values, 0.975, na.rm = TRUE),
          p_value = if (metric_name %in% c("AIC", "AICc", "BIC", "BICc")) {
            # For these metrics, negative values indicate the model is better than reference
            mean(all_delta_values < 0, na.rm = TRUE)
          } else {
            # For these metrics, positive values indicate the model is better than reference
            mean(all_delta_values > 0, na.rm = TRUE)
          }
        )
      }

      # Store summary in all_model_results for this model
      all_model_results[[model_name]]$delta_metrics <- delta_metrics[[model_name]]
      all_model_results[[model_name]]$delta_summary <- delta_summary
    }

    # Create a summary data frame of delta AIC
    delta_aic_summary <- data.frame(
      Model = setdiff(names(models_list), reference_model),
      Reference = reference_model,
      Mean_Delta_AIC = sapply(setdiff(names(models_list), reference_model),
                              function(model) all_model_results[[model]]$delta_summary$AIC$mean),
      CI_Lower = sapply(setdiff(names(models_list), reference_model),
                        function(model) all_model_results[[model]]$delta_summary$AIC$ci_lower),
      CI_Upper = sapply(setdiff(names(models_list), reference_model),
                        function(model) all_model_results[[model]]$delta_summary$AIC$ci_upper),
      P_Value = sapply(setdiff(names(models_list), reference_model),
                       function(model) all_model_results[[model]]$delta_summary$AIC$p_value),
      Better_Than_Reference = sapply(setdiff(names(models_list), reference_model),
                                     function(model) all_model_results[[model]]$delta_summary$AIC$p_value > 0.95)
    )

    # Print delta AIC summary
    cat("\n📊 Delta AIC Summary (vs. reference model):\n")
    print(delta_aic_summary)

    # Add delta_aic_summary to all_model_results
    all_model_results$delta_aic_summary <- delta_aic_summary
  }

  # After all models are processed, perform model comparisons if requested
  if (compare_models && length(models_list) > 1) {
    cat("\n🔍 Performing model comparisons...\n")

    # Initialize container for pairwise comparisons
    model_comparisons <- list()

    # Reference model performance matrix
    ref_model_perf_matrix <- all_model_results[[reference_model]]$pooled_perf_matrix

    # Calculate delta AIC/BIC values relative to reference model
    for (model_name in names(models_list)) {
      if (model_name != reference_model) {
        cat("  Comparing:", model_name, "vs", reference_model, "(reference)\n")

        # Extract performance metrics from both models
        model_perf_matrix <- all_model_results[[model_name]]$pooled_perf_matrix

        # Run comparison for each criterion
        metrics_to_compare <- c("AIC", "AICc", "BIC", "BICc", "R2", "C_Index")
        comparison_results <- list()

        for (metric in metrics_to_compare) {
          metric_idx <- which(perf_names == metric)

          if (length(metric_idx) == 0) {
            next
          }

          # Extract metric values for each bootstrap sample
          metric_values <- model_perf_matrix[, metric_idx]
          ref_metric_values <- ref_model_perf_matrix[, metric_idx]

          # Calculate differences based on whether lower or higher is better
          if (metric %in% c("AIC", "AICc", "BIC", "BICc")) {
            # For these metrics, negative differences mean current model is better than reference
            diff_values <- metric_values - ref_metric_values
            p_value <- mean(diff_values < 0, na.rm = TRUE)  # Prob current model is better than reference
            favored_model <- ifelse(mean(diff_values, na.rm = TRUE) < 0, model_name, reference_model)
          } else {
            # For these metrics, positive differences mean current model is better than reference
            diff_values <- metric_values - ref_metric_values
            p_value <- mean(diff_values > 0, na.rm = TRUE)  # Prob current model is better than reference
            favored_model <- ifelse(mean(diff_values, na.rm = TRUE) > 0, model_name, reference_model)
          }

          # Calculate mean difference and confidence interval
          mean_diff <- mean(diff_values, na.rm = TRUE)
          ci_lower <- quantile(diff_values, 0.025, na.rm = TRUE)
          ci_upper <- quantile(diff_values, 0.975, na.rm = TRUE)

          # Store results
          comparison_results[[metric]] <- list(
            metric = metric,
            model_name = model_name,
            reference_model = reference_model,
            mean_diff = mean_diff,
            ci_lower = ci_lower,
            ci_upper = ci_upper,
            p_value = p_value,
            favored_model = favored_model,
            significant = (p_value < 0.05 || p_value > 0.95),
            diff_values = diff_values
          )
        }

        # Create summary table for this comparison
        pair_name <- paste(model_name, "vs", reference_model)
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
    }

    # Identify best model for each metric
    best_model_summary <- list()

    for (metric in c("AIC", "AICc", "BIC", "BICc", "R2", "C_Index")) {
      all_values <- sapply(names(models_list), function(model_name) {
        perf_table <- all_model_results[[model_name]]$performance_table
        idx <- which(perf_table$Metric == metric)
        if (length(idx) > 0) perf_table$Estimate[idx] else NA
      })

      if (all(is.na(all_values))) {
        next
      }

      # Determine best model (lower is better for AIC/BIC, higher is better for others)
      if (metric %in% c("AIC", "AICc", "BIC", "BICc")) {
        best_idx <- which.min(all_values)
        best_model <- names(models_list)[best_idx]
      } else {
        best_idx <- which.max(all_values)
        best_model <- names(models_list)[best_idx]
      }

      # Calculate bootstrap p-values for each model vs best model
      comparison_p_values <- numeric(length(names(models_list)))
      names(comparison_p_values) <- names(models_list)

      # Best model's performance matrix
      best_perf_matrix <- all_model_results[[best_model]]$pooled_perf_matrix
      best_metric_values <- best_perf_matrix[, which(perf_names == metric)]

      for (model_name in names(models_list)) {
        if (model_name != best_model) {
          model_perf_matrix <- all_model_results[[model_name]]$pooled_perf_matrix
          model_metric_values <- model_perf_matrix[, which(perf_names == metric)]

          # Calculate differences and p-values
          if (metric %in% c("AIC", "AICc", "BIC", "BICc")) {
            # Lower is better
            diff_values <- model_metric_values - best_metric_values
            p_value <- mean(diff_values < 0, na.rm = TRUE)  # Prob model is better than best
          } else {
            # Higher is better
            diff_values <- model_metric_values - best_metric_values
            p_value <- mean(diff_values > 0, na.rm = TRUE)  # Prob model is better than best
          }

          comparison_p_values[model_name] <- p_value
        } else {
          comparison_p_values[model_name] <- NA  # No comparison with itself
        }
      }

      # Store best model info
      best_model_summary[[metric]] <- list(
        best_model = best_model,
        metric_values = all_values,
        comparison_p_values = comparison_p_values,
        significantly_better = sapply(comparison_p_values, function(p) !is.na(p) && (p < 0.05 || p > 0.95))
      )
    }

    # Print best model summary
    cat("\n--- Best Model Summary ---\n")
    for (metric in names(best_model_summary)) {
      cat(metric, ": Best model is", best_model_summary[[metric]]$best_model, "\n")
      cat("  Values for all models: \n")
      for (i in seq_along(names(models_list))) {
        model_name <- names(models_list)[i]
        value <- best_model_summary[[metric]]$metric_values[i]
        p_value <- best_model_summary[[metric]]$comparison_p_values[model_name]

        if (is.na(p_value)) {
          cat("    ", model_name, ": ", format(value, digits = 4), " (best)\n", sep = "")
        } else if (p_value < 0.05 || p_value > 0.95) {
          cat("    ", model_name, ": ", format(value, digits = 4),
              " (p = ", format(p_value, digits = 3), ", significantly different)\n", sep = "")
        } else {
          cat("    ", model_name, ": ", format(value, digits = 4),
              " (p = ", format(p_value, digits = 3), ", not significantly different)\n", sep = "")
        }
      }
      cat("\n")
    }

    # Add comparison results to output
    all_model_results$model_comparisons <- model_comparisons
    all_model_results$best_model_summary <- best_model_summary
    all_model_results$reference_model <- reference_model
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
