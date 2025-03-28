#' Assess Model Performance for Multiple Imputed Data (Optimized Sequential Version)
#'
#' This function evaluates model performance across multiple imputations with optimizations
#' and supports interaction terms, restricted cubic splines, and polynomial terms in models.
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
#' @param verbose Show progress information (default: FALSE).
#'
#' @return A list containing performance metrics, including degrees of freedom (df).
#' @import dplyr MASS Hmisc mice survival splines
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
                                 custom_model_func = NULL,
                                 verbose = FALSE) {

  # Start timing for performance analysis
  start_time <- Sys.time()

  # Import namespaces without attaching
  required_packages <- c("MASS", "Hmisc", "stats", "splines")
  if (model_type == "cox") required_packages <- c(required_packages, "survival")

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
    formula_string <- paste(outcome_var, "~", predictors_string, offset_term, trial_term)
  } else {
    # If a custom formula is provided, we need to ensure it's compatible with spline and polynomial terms
    if (!is.null(spline_terms) || !is.null(polynomial_terms)) {
      warning("Using custom formula_string with spline_terms or polynomial_terms. Make sure your formula includes these terms correctly.")
    }

    # Use the provided formula_string but add trial and offset terms if needed
    # First, check if formula_string already has the outcome_var
    if (!grepl(paste0("^", outcome_var, "\\s*~"), formula_string)) {
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
  null_formula <- stats::as.formula(paste(outcome_var, "~ 1", offset_term, trial_term))

  # Pre-allocate result vectors (faster than growing lists)
  logL_full <- numeric(imp_n)
  logL_null <- numeric(imp_n)
  df_values <- numeric(imp_n)
  c_index_values <- numeric(imp_n)
  rss_values <- numeric(imp_n)
  ae_values <- numeric(imp_n)

  # Get a reference subset for first imputation
  reference_subset <- data[data[[imp_col]] == unique_imps[1], ]
  n <- nrow(reference_subset)

  # Pre-build subsetting indices for all imputations
  # This is significantly faster than repeated filtering
  subset_indices <- vector("list", imp_n)
  for (i in 1:imp_n) {
    subset_indices[[i]] <- which(data[[imp_col]] == unique_imps[i])
  }

  # Cache model creation functions outside the loop
  fit_main_model <- function(formula, data_subset) {
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

    # Fit models
    model <- fit_main_model(model_formula, data_subset)
    null_model <- fit_main_model(null_formula, data_subset)

    # Store metrics (direct assignments, no list manipulation)
    logL_full[i] <- as.numeric(stats::logLik(model))
    logL_null[i] <- as.numeric(stats::logLik(null_model))
    df_values[i] <- length(stats::coef(model))

    # Get predictions once and reuse
    pred_values <- stats::fitted(model)
    obs_values <- data_subset[[outcome_var]]

    # Calculate performance metrics
    c_index_values[i] <- tryCatch(
      Hmisc::rcorr.cens(pred_values, obs_values)[1],
      error = function(e) NA
    )

    # Use vectorized operations
    rss_values[i] <- sum((obs_values - pred_values)^2)
    ae_values[i] <- sum(abs(obs_values - pred_values))
  }

  # Pool performance metrics with vectorized operations
  K <- mean(df_values)  # Pooled degrees of freedom

  # McFadden's R²
  R2_pooled <- 1 - (mean(logL_full) / mean(logL_null))

  # AIC and AICc
  AIC_pooled <- -2 * mean(logL_full) + 2 * K
  AICc_pooled <- AIC_pooled + (2 * K * (K + 1)) / (n - K - 1)

  # BIC and BICc calculations
  BIC_pooled <- -2 * mean(logL_full) + K * log(n)
  BICc_pooled <- BIC_pooled + (K * (K + 1)) / (n - K - 1)

  # C-index
  C_index_pooled <- mean(c_index_values, na.rm = TRUE)

  # Pooled RMSE
  W <- mean(rss_values) / n
  B <- stats::var(rss_values / n)
  T_val <- W + (1 + 1/imp_n) * B
  RMSE_pooled <- sqrt(T_val)

  # Pooled MAE
  W_mae <- mean(ae_values) / n
  B_mae <- stats::var(ae_values / n)
  T_mae <- W_mae + (1 + 1/imp_n) * B_mae
  MAE_pooled <- T_mae  # Original had sqrt here, but MAE doesn't need it

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
    Number_Imputation = imp_n,
    df = K,
    logL = mean(logL_full),
    R2 = R2_pooled,
    AIC = AIC_pooled,
    AICc = AICc_pooled,
    BIC = BIC_pooled,
    BICc = BICc_pooled,
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
