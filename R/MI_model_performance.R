#' Assess Model Performance for Multiple Imputed Data (with Random Effects Support)
#'
#' This function evaluates model performance across multiple imputations with optimizations
#' and supports interaction terms, restricted cubic splines, polynomial terms, and random effects in models.
#'
#' @name MI_model_performance
#' @param data A data frame with imputed datasets.
#' @param outcome_var The dependent variable in the model.
#' @param predictor_vars A vector of predictor variables. Can include interactions using c() and : notation.
#' @param imp_col The imputation column name (default is ".imp").
#' @param followup_offset Whether to include an offset for follow-up duration ("Yes" or "No").
#' @param followup_col Column for follow-up time (optional, required if followup_offset = "Yes").
#' @param trial_factor Whether to include trial as a factor ("Yes" or "No").
#' @param trial_col Column for trial effect (optional).
#' @param imp_n Number of imputations (if NULL, detected automatically).
#' @param model_type Model type: "nb" (negative binomial), "lm", "bin", "poisson", "gamma", "quasipoisson", "quasibinomial", or "cox".
#' @param formula_string Optional custom formula string. If provided, overrides the formula created from predictor_vars.
#' @param spline_terms A named list to specify variables to be modeled with restricted cubic splines. Each element should be a list with 'var' (variable name) and 'knots' (number of knots or explicit knot positions).
#' @param polynomial_terms A named list to specify variables to be modeled with polynomial terms. Each element should be a list with 'var' (variable name) and 'degree' (polynomial degree).
#' @param random_intercept Whether to include a random intercept in the model ("Yes" or "No").
#' @param random_intercept_var The grouping variable for the random intercept (required if random_intercept = "Yes").
#' @param verbose Show progress information (default: FALSE).
#' @param aic_method Method to compute AIC: \"corrected\" (default) or \"pooled\" (Rubin rule pooling of log_link).
#' @param r2_method Method to compute R2: \"pseudo\" (default, log-likelihood based) or \"pooled\" (classic pooled predictions).
#' @return A list containing performance metrics, including degrees of freedom (df).
#' @import dplyr MASS Hmisc mice survival splines
#' @importFrom lme4 lmer glmer fixef VarCorr
#' @importFrom glmmTMB glmmTMB fixef VarCorr
#' @importFrom coxme coxme
#' @export
MI_model_performance <- function(data,
                                 outcome_var,
                                 predictor_vars,
                                 imp_col = ".imp",
                                 followup_offset = "No",
                                 followup_col = NULL,
                                 trial_factor = "No",
                                 trial_col = NULL,
                                 imp_n = NULL,
                                 model_type = "nb",
                                 formula_string = NULL,
                                 spline_terms = NULL,
                                 polynomial_terms = NULL,
                                 random_intercept = "No",   # Added parameter
                                 random_intercept_var = NULL, # Added parameter
                                 custom_model_func = NULL,
                                 verbose = FALSE,
                                 aic_method = "corrected",
                                 r2_method = "pseudo") {

  # Start timing for performance analysis
  start_time <- Sys.time()

  # Import namespaces without attaching
  required_packages <- c("MASS", "Hmisc", "stats", "splines")
  if (model_type == "cox") required_packages <- c(required_packages, "survival")

  # Add package requirements for random effects
  if (random_intercept == "Yes") {
    required_packages <- c(required_packages, "lme4")

    if (model_type == "nb") {
      if (requireNamespace("glmmTMB", quietly = TRUE)) {
        required_packages <- c(required_packages, "glmmTMB")
      } else {
        warning("Package 'glmmTMB' not available. Using lme4 with Poisson family as approximation for negative binomial with random effects.")
      }
    }

    if (model_type == "cox") {
      if (!requireNamespace("coxme", quietly = TRUE)) {
        stop("Package 'coxme' is required for Cox models with random effects.")
      }
      required_packages <- c(required_packages, "coxme")
    }
  }

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Check if rms package is available if spline_terms is provided
  if (!is.null(spline_terms)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      message("Package 'rms' is not available. Using splines package instead.")
      use_rms <- FALSE
    } else {
      suppressPackageStartupMessages(require(rms))
      use_rms <- TRUE
    }
  } else {
    use_rms <- FALSE
  }

  # Fast validation with minimal overhead
  if (!imp_col %in% names(data)) stop("imp_col not found in data.")
  if (!outcome_var %in% names(data) && model_type != "cox") stop("outcome_var not found in data.")

  # Validate random_intercept input
  if (!random_intercept %in% c("Yes", "No")) stop("random_intercept must be either 'Yes' or 'No'.")

  if (random_intercept == "Yes" && is.null(random_intercept_var)) {
    stop("If random_intercept = 'Yes', random_intercept_var must be provided.")
  }

  if (!is.null(random_intercept_var) && !random_intercept_var %in% names(data)) {
    stop("random_intercept_var not found in data.")
  }

  # Validate spline_terms if provided
  if (!is.null(spline_terms)) {
    if (!is.list(spline_terms)) {
      stop("spline_terms must be a list")
    }

    for (i in seq_along(spline_terms)) {
      if (!is.list(spline_terms[[i]]) || is.null(spline_terms[[i]]$var)) {
        stop("Each element in spline_terms must be a list with at least a 'var' element")
      }

      if (!spline_terms[[i]]$var %in% names(data)) {
        stop(paste("Spline variable", spline_terms[[i]]$var, "not found in data"))
      }

      if (is.null(spline_terms[[i]]$knots)) {
        # Default to 3 knots if not specified
        spline_terms[[i]]$knots <- 3
      } else if (!is.numeric(spline_terms[[i]]$knots)) {
        stop("Knots must be numeric")
      }
    }
  }

  # Validate polynomial_terms if provided
  if (!is.null(polynomial_terms)) {
    if (!is.list(polynomial_terms)) {
      stop("polynomial_terms must be a list")
    }

    for (i in seq_along(polynomial_terms)) {
      if (!is.list(polynomial_terms[[i]]) || is.null(polynomial_terms[[i]]$var)) {
        stop("Each element in polynomial_terms must be a list with at least a 'var' element")
      }

      if (!polynomial_terms[[i]]$var %in% names(data)) {
        stop(paste("Polynomial variable", polynomial_terms[[i]]$var, "not found in data"))
      }

      if (is.null(polynomial_terms[[i]]$degree)) {
        # Default to quadratic (degree 2) if not specified
        polynomial_terms[[i]]$degree <- 2
      } else if (!is.numeric(polynomial_terms[[i]]$degree) ||
                 polynomial_terms[[i]]$degree < 2 ||
                 polynomial_terms[[i]]$degree != round(polynomial_terms[[i]]$degree)) {
        stop("Degree must be a positive integer greater than or equal to 2")
      }
    }
  }

  # Check if formula_string is provided, otherwise validate predictor_vars
  if (is.null(formula_string)) {
    # Check for interaction terms in predictor_vars
    has_interactions <- any(grepl(":", predictor_vars))

    if (!has_interactions) {
      # Regular variable check if no interactions
      missing_predictors <- predictor_vars[!predictor_vars %in% names(data)]
      if (length(missing_predictors) > 0) {
        stop(paste("Predictor variables not found in data:", paste(missing_predictors, collapse = ", ")))
      }
    } else {
      # For interactions, we need to parse them more carefully
      # Extract all unique variable names from interaction terms
      all_vars <- character(0)
      for (term in predictor_vars) {
        if (grepl(":", term)) {
          # Split interaction terms
          interaction_vars <- unlist(strsplit(term, ":"))
          # Clean up any whitespace
          interaction_vars <- trimws(interaction_vars)
          all_vars <- c(all_vars, interaction_vars)
        } else {
          all_vars <- c(all_vars, term)
        }
      }
      all_vars <- unique(all_vars)

      # Check if all extracted variables exist
      missing_predictors <- all_vars[!all_vars %in% names(data)]
      if (length(missing_predictors) > 0) {
        stop(paste("Predictor variables not found in data:", paste(missing_predictors, collapse = ", ")))
      }
    }
  }

  # Validate only what's needed based on options
  if (followup_offset == "Yes") {
    if (is.null(followup_col)) stop("If followup_offset = 'Yes', followup_col must be provided.")
    if (!followup_col %in% names(data)) stop("followup_col not found in data.")
  }

  if (trial_factor == "Yes") {
    if (is.null(trial_col)) stop("If trial_factor = 'Yes', trial_col must be provided.")
    if (!trial_col %in% names(data)) stop("trial_col not found in data.")
  }

  # Get unique imputation values efficiently using table()
  imp_table <- table(data[[imp_col]])
  unique_imps <- as.numeric(names(imp_table))

  if (is.null(imp_n)) {
    imp_n <- length(unique_imps)
    if (imp_n < 2) stop("At least 2 imputations required.")
  }

  if (verbose) message(paste("Processing", imp_n, "imputations..."))

  # Precompute formula only once - major performance boost
  trial_term <- if (trial_factor == "Yes") paste("+ as.factor(", trial_col, ")") else ""
  offset_term <- if (followup_offset == "Yes") paste("+ offset(log(", followup_col, "))") else ""

  # Random intercept term for formula construction
  random_term <- if (random_intercept == "Yes") paste("+ (1 | ", random_intercept_var, ")") else ""

  # Process spline terms if provided
  spline_formula_parts <- character(0)
  if (!is.null(spline_terms)) {
    for (i in seq_along(spline_terms)) {
      var_name <- spline_terms[[i]]$var
      knots <- spline_terms[[i]]$knots

      # Remove the original variable from predictor_vars if it's there
      if (!is.null(predictor_vars)) {
        predictor_vars <- predictor_vars[predictor_vars != var_name]
      }

      # Create the spline term
      if (use_rms) {
        # Use rcs() from rms package
        if (length(knots) == 1) {
          # Number of knots specified
          spline_formula_parts <- c(spline_formula_parts,
                                    paste0("rcs(", var_name, ", ", knots, ")"))
        } else {
          # Explicit knot positions specified
          spline_formula_parts <- c(spline_formula_parts,
                                    paste0("rcs(", var_name, ", c(",
                                           paste(knots, collapse = ", "), "))"))
        }
      } else {
        # Use bs() from splines package as a fallback
        if (length(knots) == 1) {
          # For bs(), df = knots + degree - 1 (for cubic splines with degree=3)
          df <- knots + 2  # degree=3, so 3-1=2
          spline_formula_parts <- c(spline_formula_parts,
                                    paste0("bs(", var_name, ", df = ", df, ", degree = 3)"))
        } else {
          # Explicit knot positions
          spline_formula_parts <- c(spline_formula_parts,
                                    paste0("bs(", var_name, ", knots = c(",
                                           paste(knots, collapse = ", "), "), degree = 3)"))
        }
      }
    }
  }

  # Process polynomial terms if provided
  poly_formula_parts <- character(0)
  if (!is.null(polynomial_terms)) {
    for (i in seq_along(polynomial_terms)) {
      var_name <- polynomial_terms[[i]]$var
      degree <- polynomial_terms[[i]]$degree

      # Remove the original variable from predictor_vars if it's there
      if (!is.null(predictor_vars)) {
        predictor_vars <- predictor_vars[predictor_vars != var_name]
      }

      # Create the polynomial term using poly() function
      poly_formula_parts <- c(poly_formula_parts,
                              paste0("poly(", var_name, ", degree = ", degree, ", raw = TRUE)"))
    }
  }

  if (is.null(formula_string)) {
    # Combine regular predictors, spline terms, and polynomial terms
    all_terms <- c(predictor_vars, spline_formula_parts, poly_formula_parts)

    # Construct formula from predictor_vars (handling interactions properly)
    predictors_string <- paste(all_terms, collapse = " + ")

    if (model_type == "cox") {
      formula_string <- paste("Surv(time, event) ~", predictors_string, trial_term, random_term)
    } else {
      formula_string <- paste(outcome_var, "~", predictors_string, offset_term, trial_term, random_term)
    }
  } else {
    # If a custom formula is provided, we need to ensure it's compatible with random effects
    if (random_intercept == "Yes" && !grepl("\\|", formula_string)) {
      # Add random effect if not present in custom formula
      formula_string <- paste(formula_string, random_term)
    }

    # Use the provided formula_string but add trial and offset terms if needed
    # First, check if formula_string already has the outcome_var
    if (!grepl(paste0("^", outcome_var, "\\s*~"), formula_string) && model_type != "cox") {
      formula_string <- paste(outcome_var, "~", formula_string)
    }

    # Add offset and trial terms if they're not already included
    if (followup_offset == "Yes" && !grepl("offset\\(log\\(.*\\)\\)", formula_string)) {
      formula_string <- paste(formula_string, offset_term)
    }

    if (trial_factor == "Yes" && !grepl(paste0("as\\.factor\\(", trial_col, "\\)"), formula_string)) {
      formula_string <- paste(formula_string, trial_term)
    }
  }

  # Store interaction terms for later reference
  interaction_terms <- character(0)
  if (grepl(":", formula_string)) {
    # Extract terms after the "~"
    terms_part <- strsplit(formula_string, "~")[[1]][2]
    # Split by "+"
    terms <- strsplit(terms_part, "\\+")[[1]]
    # Trim whitespace
    terms <- trimws(terms)
    # Find terms containing ":"
    interaction_terms <- terms[grepl(":", terms)]
  }

  # Store spline terms for later reference
  spline_terms_detected <- character(0)
  if (grepl("rcs\\(|bs\\(", formula_string)) {
    # Extract terms after the "~"
    terms_part <- strsplit(formula_string, "~")[[1]][2]
    # Split by "+"
    terms <- strsplit(terms_part, "\\+")[[1]]
    # Trim whitespace
    terms <- trimws(terms)
    # Find terms containing "rcs(" or "bs("
    spline_terms_detected <- terms[grepl("rcs\\(|bs\\(", terms)]
  }

  # Store polynomial terms for later reference
  poly_terms_detected <- character(0)
  if (grepl("poly\\(", formula_string)) {
    # Extract terms after the "~"
    terms_part <- strsplit(formula_string, "~")[[1]][2]
    # Split by "+"
    terms <- strsplit(terms_part, "\\+")[[1]]
    # Trim whitespace
    terms <- trimws(terms)
    # Find terms containing "poly("
    poly_terms_detected <- terms[grepl("poly\\(", terms)]
  }

  model_formula <- stats::as.formula(formula_string)

  # Create null formula preserving random effects
  if (random_intercept == "Yes") {
    null_formula <- stats::as.formula(paste(outcome_var, "~ 1", offset_term, trial_term, random_term))
  } else {
    null_formula <- stats::as.formula(paste(outcome_var, "~ 1", offset_term, trial_term))
  }

  # Pre-allocate result vectors (faster than growing lists)
  logL_full <- numeric(imp_n)
  logL_null <- numeric(imp_n)
  df_values <- numeric(imp_n)
  c_index_values <- numeric(imp_n)
  rss_values <- numeric(imp_n)
  ae_values <- numeric(imp_n)
  r2_values <- numeric(imp_n)

  # Get a reference subset for first imputation
  reference_subset <- data[data[[imp_col]] == unique_imps[1], ]
  n <- nrow(reference_subset)
  subject_n <- NA






  # Pre-build subsetting indices for all imputations
  # This is significantly faster than repeated filtering
  subset_indices <- vector("list", imp_n)
  for (i in 1:imp_n) {
    subset_indices[[i]] <- which(data[[imp_col]] == unique_imps[i])
  }

  # Cache model creation functions outside the loop
  fit_main_model <- function(formula, data_subset) {
    if (random_intercept == "Yes") {
      # Models with random effects
      switch(model_type,
             "lm" = lme4::lmer(formula, data = data_subset),
             "bin" = lme4::glmer(formula, family = stats::binomial(), data = data_subset),
             "poisson" = lme4::glmer(formula, family = stats::poisson(), data = data_subset),
             "nb" = {
               if (requireNamespace("glmmTMB", quietly = TRUE)) {
                 glmmTMB::glmmTMB(formula, family = glmmTMB::nbinom2, data = data_subset)
               } else {
                 # Fallback to Poisson if glmmTMB not available
                 lme4::glmer(formula, family = stats::poisson(), data = data_subset)
               }
             },
             "cox" = coxme::coxme(formula, data = data_subset),
             "gamma" = {
               if (requireNamespace("glmmTMB", quietly = TRUE)) {
                 glmmTMB::glmmTMB(formula, family = stats::Gamma, data = data_subset)
               } else {
                 # Gamma family may not be stable in lme4
                 warning("Gamma family with random effects may not be stable in lme4.")
                 lme4::glmer(formula, family = stats::Gamma, data = data_subset)
               }
             },
             stop("Unsupported model type for random effects.")
      )
    } else {
      # Standard models without random effects (original code)
      switch(model_type,
             "nb" = MASS::glm.nb(formula, data = data_subset),
             "lm" = stats::glm(formula, family = stats::gaussian(), data = data_subset),
             "bin" = stats::glm(formula, family = stats::binomial(), data = data_subset),
             "poisson" = stats::glm(formula, family = stats::poisson(), data = data_subset),
             "gamma" = stats::glm(formula, family = stats::Gamma(), data = data_subset),
             "quasipoisson" = stats::glm(formula, family = stats::quasipoisson(), data = data_subset),
             "quasibinomial" = stats::glm(formula, family = stats::quasibinomial(), data = data_subset),
             "cox" = survival::coxph(formula, data = data_subset),
             stop("Unsupported model type.")
      )
    }
  }

  # Function to extract log-likelihood from different model types
  extract_logLik <- function(model) {
    if (inherits(model, "merMod")) {
      # Use stats::logLik for lmer/glmer models
      return(as.numeric(stats::logLik(model)))
    } else if (inherits(model, "glmmTMB")) {
      # Use logLik for glmmTMB models
      return(as.numeric(stats::logLik(model)))
    } else if (inherits(model, "coxme")) {
      # For coxme models, logLik is accessible but may need special handling
      return(as.numeric(model$loglik[2]))  # The second element is the fitted model log-likelihood
    } else {
      # For standard models, use stats::logLik
      return(as.numeric(stats::logLik(model)))
    }
  }

  # Function to get predictions from different model types
  get_predictions <- function(model, data_subset) {
    if (inherits(model, "merMod")) {
      # For lme4 models, use fitted values
      return(stats::fitted(model))
    } else if (inherits(model, "glmmTMB")) {
      # For glmmTMB models, use fitted values
      return(stats::fitted(model))
    } else if (inherits(model, "coxme")) {
      # For Cox models with random effects, use linear predictor
      return(predict(model, type = "lp"))
    } else {
      # Standard models
      return(stats::fitted(model))
    }
  }

  # Function to get fixed effects degrees of freedom
  get_df <- function(model) {
    if (inherits(model, "merMod")) {
      # For lme4 models, count fixed effects
      return(length(lme4::fixef(model)))
    } else if (inherits(model, "glmmTMB")) {
      # For glmmTMB models, count fixed effects
      return(length(glmmTMB::fixef(model)$cond))
    } else if (inherits(model, "coxme")) {
      # For coxme models, get fixed effects
      return(length(coxme::fixef(model)))
    } else {
      # Standard models
      return(length(stats::coef(model)))
    }
  }

  # Loop through imputations - optimized sequential processing
  for (i in 1:imp_n) {
    if (verbose && (i %% 2 == 0 || i == 1 || i == imp_n)) {
      message(paste("  Processing imputation", i, "of", imp_n))
    }

    # Get data subset using pre-computed indices (very fast)
    data_subset <- data[subset_indices[[i]], ]

    # Ensure required packages are loaded for the modeling
    if (length(spline_terms_detected) > 0) {
      if (any(grepl("rcs\\(", spline_terms_detected)) && use_rms) {
        # For rcs(), we need rms package
        suppressPackageStartupMessages(require(rms))
      }
      if (any(grepl("bs\\(", spline_terms_detected))) {
        # For bs(), we need splines package
        suppressPackageStartupMessages(require(splines))
      }
    }

    # Fit models with error handling for convergence issues
    model <- tryCatch({
      fit_main_model(model_formula, data_subset)
    }, error = function(e) {
      warning(paste("Error fitting full model for imputation", i, ":", conditionMessage(e)))
      return(NULL)
    })

    # Add this block to capture the actual number of subjects used in the model
    if (!is.null(model)) {
      # For the first successful model, extract number of subjects included
      if (is.na(subject_n)) {
        # Different model types store this information differently
        if (inherits(model, "glm") || inherits(model, "negbin")) {
          subject_n <- length(model$fitted.values)
        } else if (inherits(model, "lm")) {
          subject_n <- length(model$fitted.values)
        } else if (inherits(model, "coxph")) {
          subject_n <- model$n
        } else if (inherits(model, "merMod")) {
          subject_n <- nrow(model@frame)
        } else if (inherits(model, "glmmTMB")) {
          subject_n <- length(fitted(model))
        } else if (inherits(model, "coxme")) {
          subject_n <- model$n
        } else {
          # Fallback: try to get number of fitted values
          tryCatch({
            subject_n <- length(fitted(model))
          }, error = function(e) {
            warning("Could not determine number of subjects in model. Using total count.")
            subject_n <- n
          })
        }
      }
    }

    null_model <- tryCatch({
      fit_main_model(null_formula, data_subset)
    }, error = function(e) {
      warning(paste("Error fitting null model for imputation", i, ":", conditionMessage(e)))
      return(NULL)
    })

    # Skip this imputation if either model failed
    if (is.null(model) || is.null(null_model)) {
      logL_full[i] <- NA
      logL_null[i] <- NA
      df_values[i] <- NA
      c_index_values[i] <- NA
      rss_values[i] <- NA
      ae_values[i] <- NA
      next
    }






    # Store metrics (direct assignments, no list manipulation)
    logL_full[i] <- extract_logLik(model)
    logL_null[i] <- extract_logLik(null_model)
    df_values[i] <- get_df(model)

    # Get predictions once and reuse
    pred_values <- get_predictions(model, data_subset)

    # Get observed values
    if (model_type == "cox") {
      obs_values <- with(data_subset, Surv(time, event))
    } else {
      obs_values <- data_subset[[outcome_var]]
    }

    # Calculate R2 for each imputation (only if not Cox)
    if (model_type != "cox") {
      r2_values[i] <- 1 - sum((obs_values - pred_values)^2) / sum((obs_values - mean(obs_values))^2)
    }

    # For outcome_var we need to handle Cox models differently
    if (model_type == "cox") {
      # For Cox models, use survival time and event as observed value
      obs_values <- with(data_subset, Surv(time, event))
    } else {
      obs_values <- data_subset[[outcome_var]]
    }

    # Calculate performance metrics
    if (model_type == "cox") {
      # For Cox models, use a specific concordance approach
      c_index_values[i] <- tryCatch(
        survival::concordance(model)$concordance,
        error = function(e) NA
      )
      # For Cox, we don't calculate RMSE or MAE in the usual way
      rss_values[i] <- NA
      ae_values[i] <- NA
    } else {
      # For other models, use standard metrics
      c_index_values[i] <- tryCatch(
        Hmisc::rcorr.cens(pred_values, obs_values)[1],
        error = function(e) NA
      )
      # Use vectorized operations
      rss_values[i] <- sum((obs_values - pred_values)^2)
      ae_values[i] <- sum(abs(obs_values - pred_values))
    }
  }

  # Pool performance metrics with vectorized operations
  valid_imps <- !is.na(logL_full)
  n_valid <- sum(valid_imps)

  if (n_valid == 0) {
    stop("All models failed to converge. Check your model specification.")
  }
  # Finally, before returning the result at the end, add this check
  # to handle cases where we couldn't determine subject_n
  if (is.na(subject_n)) {
    warning("Could not determine number of subjects in model. Using total count.")
    subject_n <- n
  }
  # Use only valid imputations for calculations
  K <- mean(df_values[valid_imps], na.rm = TRUE)  # Pooled degrees of freedom

  # Between-imputation variance of log-likelihoods
  B_logL <- var(logL_full[valid_imps], na.rm = TRUE)

  # Within-imputation variance (approximated by degrees of freedom)
  W_logL <- K

  # Total variance with multiple imputation correction factor
  T_logL <- W_logL + (1 + 1/n_valid) * B_logL

  # Pooled log-likelihood
  logL_pooled <- mean(logL_full[valid_imps], na.rm = TRUE)

  if (tolower(r2_method) == "pseudo") {
    if (model_type == "cox") {
      R2_pooled <- NA
      R2_method_used <- "NA (Cox model - no standard R2)"
    } else {
      R2_pooled <- 1 - (logL_pooled / mean(logL_null[valid_imps], na.rm = TRUE))
      R2_method_used <- "Pseudo R2 (based on log-likelihoods)"
    }
  } else if (tolower(r2_method) == "pooled") {
    R2_pooled <- mean(r2_values, na.rm = TRUE)
    R2_method_used <- "Pooled R2 (average of per-imputation R2)"
  } else {
    stop("Invalid value for r2_method. Must be 'pseudo' or 'pooled'.")
  }

  if (tolower(aic_method) == "corrected") {
    # AIC using correction for between-imputation variance
    AIC_pooled <- -2 * logL_pooled + 2 * T_logL
    AICc_pooled <- AIC_pooled + (2 * T_logL * (T_logL + 1)) / (n - T_logL - 1)

    BIC_pooled <- -2 * logL_pooled + T_logL * log(n)
    BICc_pooled <- BIC_pooled + (T_logL * (T_logL + 1)) / (n - T_logL - 1)
    AIC_method_used <- "Corrected for between-imputation variance"
  } else if (tolower(aic_method) == "pooled") {
    # AIC based only on pooled logLik and pooled degrees of freedom (K)
    AIC_pooled <- -2 * logL_pooled + 2 * K
    AICc_pooled <- AIC_pooled + (2 * K * (K + 1)) / (n - K - 1)

    BIC_pooled <- -2 * logL_pooled + K * log(n)
    BICc_pooled <- BIC_pooled + (K * (K + 1)) / (n - K - 1)
    AIC_method_used <- "Based on pooled logLik and degrees of freedom"
  } else {
    stop("Invalid value for aic_method. Must be 'corrected' or 'pooled_loglik'.")
  }

  # C-index
  C_index_pooled <- mean(c_index_values, na.rm = TRUE)

  # Pooled RMSE and MAE (only for non-Cox models)
  if (model_type != "cox") {
    # Pooled RMSE
    W <- mean(rss_values[valid_imps], na.rm = TRUE) / n
    B <- stats::var(rss_values[valid_imps] / n, na.rm = TRUE)
    T_val <- W + (1 + 1/n_valid) * B
    RMSE_pooled <- sqrt(T_val)

    # Pooled MAE
    W_mae <- mean(ae_values[valid_imps], na.rm = TRUE) / n
    B_mae <- stats::var(ae_values[valid_imps] / n, na.rm = TRUE)
    T_mae <- W_mae + (1 + 1/n_valid) * B_mae
    MAE_pooled <- T_mae
  } else {
    RMSE_pooled <- NA
    MAE_pooled <- NA
  }

  # Calculate execution time
  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "secs")

  if (verbose) {
    message(paste("MI_model_performance execution completed in",
                  round(as.numeric(execution_time), 2), "seconds"))
  }

  # Return performance metrics with BIC and BICc added, plus interaction/spline/polynomial information
  result <- list(
    Model_Formula = formula_string,
    Model_Type = model_type,
    Number_Subjects = subject_n,
    Number_Imputation = imp_n,
    Number_Valid_Imputation = n_valid,  # Also add this for clarity
    df = K,
    B_logL = B_logL,  # Add between-imputation variance
    T_logL = T_logL,  # Add total variance-adjusted df
    logL = logL_pooled,  # Update to use logL_pooled variable name
    R2 = R2_pooled,
    R2_Method = R2_method_used,
    AIC = AIC_pooled,
    AICc = AICc_pooled,
    BIC = BIC_pooled,
    BICc = BICc_pooled,
    AIC_method = AIC_method_used,
    C_Index = C_index_pooled,
    RMSE = RMSE_pooled,
    MAE = MAE_pooled,
    Has_Interactions = length(interaction_terms) > 0,
    Interaction_Terms = if(length(interaction_terms) > 0) interaction_terms else NULL,
    Has_Splines = length(spline_terms_detected) > 0,
    Spline_Terms = if(length(spline_terms_detected) > 0) spline_terms_detected else NULL,
    Spline_Method = if(use_rms) "rcs" else "bs",
    Has_Polynomials = length(poly_terms_detected) > 0,
    Polynomial_Terms = if(length(poly_terms_detected) > 0) poly_terms_detected else NULL
  )
  return(result)
}
