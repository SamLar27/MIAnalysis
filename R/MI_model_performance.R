#' Assess Model Performance for Multiple Imputed Data (Predictor/Covariable Style)
#'
#' This function evaluates model performance across multiple imputations with the same
#' structure as MI_estimates: testing individual predictors with fixed covariables.
#' All R2 calculations are compared against an intercept-only model for meaningful comparison.
#' Modified to handle RCS terms directly in predictor_vars and robust null model support.
#'
#' @param data A data frame with imputed datasets.
#' @param outcome_var The dependent variable in the model.
#' @param predictor_vars A vector of predictor variables to be tested individually. Each predictor will be tested separately along with all covariables. Can include "" for null model. Can now include RCS terms like "rcs(FEV1_preBD_PCT_0W,4)".
#' @param covariables A vector of covariable names to be included as fixed effects in all models.
#' @param imp_col The imputation column name (default is ".imp").
#' @param imp_n Number of imputations (if NULL, detected automatically).
#' @param model_type Model type: "nb", "lm", "bin", "poisson", "gamma", "quasipoisson", "quasibinomial", or "cox".
#' @param followup_offset Whether to include an offset for follow-up duration ("Yes" or "No").
#' @param followup_col Column for follow-up time (required if followup_offset = "Yes").
#' @param trial_factor Whether to include trial as a factor ("Yes" or "No").
#' @param trial_col Column for trial effect (required if trial_factor = "Yes").
#' @param time_col The time variable for Cox regression (only required if model_type = "cox").
#' @param event_col The event variable for Cox regression (only required if model_type = "cox").
#' @param formula_string Optional custom formula string. If provided, overrides the formula created from predictor_vars and covariables.
#' @param spline_terms A named list to specify variables to be modeled with restricted cubic splines.
#' @param poly_terms A named list to specify variables to be modeled with polynomial terms.
#' @param random_intercept Whether to include a random intercept in the model ("Yes" or "No").
#' @param random_intercept_var The grouping variable for the random intercept (required if random_intercept = "Yes").
#' @param verbose Show progress information (default: FALSE).
#' @param aic_method Method to compute AIC: "Rubin" (default) or "Simple_pooled".
#' @param r2_method Method to compute R2: "pseudo" (default, log-likelihood based) or "pooled" (classic pooled predictions).
#' @return A list containing performance metrics for each predictor, with combined results if multiple predictors.
#' @import dplyr MASS mice survival splines
#' @importFrom lme4 lmer glmer fixef VarCorr
#' @importFrom glmmTMB glmmTMB fixef VarCorr
#' @importFrom coxme coxme
#' @export
MI_model_performance <- function(data,
                                 outcome_var,
                                 predictor_vars,
                                 covariables = NULL,
                                 imp_col = ".imp",
                                 imp_n = NULL,
                                 model_type = "nb",
                                 followup_offset = "No",
                                 followup_col = NULL,
                                 trial_factor = "No",
                                 trial_col = NULL,
                                 time_col = NULL,
                                 event_col = NULL,
                                 formula_string = NULL,
                                 spline_terms = NULL,
                                 poly_terms = NULL,
                                 random_intercept = "No",
                                 random_intercept_var = NULL,
                                 verbose = FALSE,
                                 aic_method = "Rubin",
                                 r2_method = "pseudo") {

  # Load required packages
  require(dplyr)
  require(MASS)
  require(mice)
  require(rlang)
  require(survival)
  require(splines)

  # Load packages for mixed models if needed
  if (random_intercept == "Yes") {
    require(lme4)
    require(mitools)

    if (model_type == "cox") {
      if (!requireNamespace("coxme", quietly = TRUE)) {
        stop("Package 'coxme' is required for Cox models with random effects.")
      }
      require(coxme)
    }

    if (model_type == "nb") {
      if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        warning("Package 'glmmTMB' not available. Using 'lme4' with Poisson family as approximation.")
      } else {
        require(glmmTMB)
      }
    }
  }

  # Handle empty or NULL predictor_vars - ENHANCED
  if (is.null(predictor_vars) || length(predictor_vars) == 0) {
    predictor_vars <- c("")  # Convert to empty string for null model
    if (verbose) message("No predictor variables provided. Testing null model (covariables only).")
  }

  # Clean predictor_vars - handle various empty cases
  predictor_vars <- as.character(predictor_vars)  # Ensure character vector
  predictor_vars[is.na(predictor_vars)] <- ""     # Convert NAs to empty strings

  # If all predictors are empty after cleaning, keep one empty string
  if (all(predictor_vars == "" | is.na(predictor_vars))) {
    predictor_vars <- c("")
    if (verbose) message("All predictor variables are empty. Testing null model (covariables only).")
  }

  # Check for RCS/spline terms in predictor_vars and load appropriate package
  has_rcs_in_predictors <- any(grepl("rcs\\(|bs\\(", predictor_vars))

  if (has_rcs_in_predictors || !is.null(spline_terms)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      message("Package 'rms' is not available. Using splines package instead.")
      use_rms <- FALSE
    } else {
      require(rms)
      use_rms <- TRUE
    }
  } else {
    use_rms <- FALSE
  }

  # Input validation
  if (!imp_col %in% names(data)) stop("imp_col not found in data.")
  if (!outcome_var %in% names(data) && model_type != "cox") stop("outcome_var not found in data.")

  # Determine number of imputations
  actual_imps <- sort(unique(data[[imp_col]]))
  if (is.null(imp_n)) {
    imp_n <- length(actual_imps)
  } else {
    if (imp_n != length(actual_imps)) {
      warning("imp_n does not match the number of unique imputations found in the data. Using only present imputations.")
      imp_n <- length(actual_imps)
    }
  }

  # Validate inputs
  if (followup_offset == "Yes" && is.null(followup_col)) {
    stop("If followup_offset = 'Yes', followup_col must be provided.")
  }
  if (trial_factor == "Yes" && is.null(trial_col)) {
    stop("If trial_factor = 'Yes', trial_col must be provided.")
  }
  if (random_intercept == "Yes" && is.null(random_intercept_var)) {
    stop("If random_intercept = 'Yes', random_intercept_var must be provided.")
  }
  if (model_type == "cox" && (is.null(time_col) || is.null(event_col))) {
    stop("For Cox regression, time_col and event_col must be provided.")
  }

  # Helper function to extract variable names from RCS/spline terms
  extract_vars_from_rcs <- function(term) {
    if (is.na(term) || term == "") return("")  # Handle empty terms

    if (grepl("^rcs\\(", term)) {
      var_name <- gsub("rcs\\(([^,]+),.*\\)", "\\1", term)
      return(trimws(var_name))
    } else if (grepl("^bs\\(", term)) {
      var_name <- gsub("bs\\(([^,]+),.*\\)", "\\1", term)
      return(trimws(var_name))
    } else if (grepl("^poly\\(", term)) {
      var_name <- gsub("poly\\(([^,]+),.*\\)", "\\1", term)
      return(trimws(var_name))
    } else {
      return(term)
    }
  }

  # Validate predictor_vars and covariables - MODIFIED for RCS support and null model
  if (is.null(formula_string)) {
    # Helper function to extract individual variables from interaction terms
    extract_vars_from_terms <- function(terms) {
      all_vars <- character(0)
      for (term in terms) {
        if (is.na(term) || term == "") next  # Skip empty terms

        if (grepl("\\*", term)) {
          vars_in_interaction <- unlist(strsplit(term, "\\*"))
          vars_in_interaction <- trimws(vars_in_interaction)
          # Extract variables from RCS terms within interactions
          vars_in_interaction <- sapply(vars_in_interaction, extract_vars_from_rcs)
          all_vars <- c(all_vars, vars_in_interaction)
        } else {
          # Extract variable from RCS term or use as is
          var_name <- extract_vars_from_rcs(term)
          all_vars <- c(all_vars, var_name)
        }
      }
      return(unique(all_vars[all_vars != ""]))  # Remove empty strings
    }

    # Extract individual variables from predictor_vars (handling interactions and RCS)
    non_empty_predictors <- predictor_vars[predictor_vars != "" & !is.na(predictor_vars)]
    if (length(non_empty_predictors) > 0) {
      individual_predictor_vars <- extract_vars_from_terms(non_empty_predictors)
      if (length(individual_predictor_vars) > 0 && !all(individual_predictor_vars %in% names(data))) {
        missing_vars <- individual_predictor_vars[!individual_predictor_vars %in% names(data)]
        stop(paste("Predictor variables not found in data:", paste(missing_vars, collapse = ", ")))
      }
    }

    # Check covariables
    if (!is.null(covariables)) {
      individual_covariable_vars <- extract_vars_from_terms(covariables)
      if (length(individual_covariable_vars) > 0 && !all(individual_covariable_vars %in% names(data))) {
        missing_vars <- individual_covariable_vars[!individual_covariable_vars %in% names(data)]
        stop(paste("Covariables not found in data:", paste(missing_vars, collapse = ", ")))
      }
    }
  }

  # Helper function to expand interaction terms - MODIFIED for RCS support
  expand_terms <- function(term) {
    if (is.na(term) || term == "") return(character(0))  # Handle empty terms

    # If it's already an RCS/spline/poly function, don't expand
    if (grepl("^rcs\\(|^bs\\(|^poly\\(", term)) {
      return(term)
    }

    if (grepl("\\*", term)) {
      vars_split <- unlist(strsplit(term, "\\*"))
      vars_split <- trimws(vars_split)
      all_combinations <- lapply(1:length(vars_split), function(k) {
        combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
      })
      unlist(all_combinations)
    } else {
      term
    }
  }

  # Helper function to process spline and polynomial terms
  process_special_terms <- function() {
    spline_formula_parts <- character(0)
    if (!is.null(spline_terms)) {
      for (i in seq_along(spline_terms)) {
        var_name <- spline_terms[[i]]$var
        knots <- spline_terms[[i]]$knots

        if (use_rms) {
          if (length(knots) == 1) {
            spline_formula_parts <- c(spline_formula_parts,
                                      paste0("rcs(", var_name, ", ", knots, ")"))
          } else {
            spline_formula_parts <- c(spline_formula_parts,
                                      paste0("rcs(", var_name, ", c(",
                                             paste(knots, collapse = ", "), "))"))
          }
        } else {
          if (length(knots) == 1) {
            df <- knots + 2
            spline_formula_parts <- c(spline_formula_parts,
                                      paste0("bs(", var_name, ", df = ", df, ", degree = 3)"))
          } else {
            spline_formula_parts <- c(spline_formula_parts,
                                      paste0("bs(", var_name, ", knots = c(",
                                             paste(knots, collapse = ", "), "), degree = 3)"))
          }
        }
      }
    }

    poly_formula_parts <- character(0)
    if (!is.null(poly_terms)) {
      for (i in seq_along(poly_terms)) {
        var_name <- poly_terms[[i]]$var
        degree <- poly_terms[[i]]$degree
        poly_formula_parts <- c(poly_formula_parts,
                                paste0("poly(", var_name, ", degree = ", degree, ", raw = TRUE)"))
      }
    }

    return(list(spline_parts = spline_formula_parts, poly_parts = poly_formula_parts))
  }

  # Get special terms
  special_terms <- process_special_terms()
  spline_formula_parts <- special_terms$spline_parts
  poly_formula_parts <- special_terms$poly_parts

  # Helper function to build formula for a given predictor - ENHANCED for null model
  build_formula_for_predictor <- function(current_predictor) {
    trial_term <- if (trial_factor == "Yes")
      paste("+ as.factor(", trial_col, ")") else ""
    offset_term <- if (followup_offset == "Yes")
      paste("+ offset(log(", followup_col, "))") else ""
    random_term <- if (random_intercept == "Yes")
      paste("+ (1 | ", random_intercept_var, ")") else ""

    if (is.null(formula_string)) {
      # Handle empty predictor (null model) - ENHANCED
      if (is.na(current_predictor) || current_predictor == "" || is.null(current_predictor)) {
        if (!is.null(covariables) && length(covariables) > 0) {
          expanded_covariables <- unlist(lapply(covariables, expand_terms))
          expanded_covariables <- unique(expanded_covariables[expanded_covariables != ""])
        } else {
          expanded_covariables <- character(0)
        }
        all_terms <- c(expanded_covariables, spline_formula_parts, poly_formula_parts)
      } else {
        # Regular predictor case - MODIFIED to handle RCS terms directly
        current_spline_parts <- spline_formula_parts
        current_poly_parts <- poly_formula_parts

        # Check if current_predictor is an RCS/spline/poly term
        is_rcs_term <- grepl("^rcs\\(|^bs\\(|^poly\\(", current_predictor)

        if (is_rcs_term) {
          # If current predictor is already an RCS/spline/poly term, use it directly
          expanded_predictor <- current_predictor

          # Extract the base variable name to avoid conflicts
          base_var <- extract_vars_from_rcs(current_predictor)

          # Remove any spline/poly versions of the base variable from special terms
          if (!is.null(spline_terms)) {
            for (i in seq_along(spline_terms)) {
              if (spline_terms[[i]]$var == base_var) {
                pattern <- if (use_rms) {
                  paste0("rcs\\(", base_var, ",")
                } else {
                  paste0("bs\\(", base_var, ",")
                }
                current_spline_parts <- current_spline_parts[!grepl(pattern, current_spline_parts)]
              }
            }
          }

          if (!is.null(poly_terms)) {
            for (i in seq_along(poly_terms)) {
              if (poly_terms[[i]]$var == base_var) {
                pattern <- paste0("poly\\(", base_var, ",")
                current_poly_parts <- current_poly_parts[!grepl(pattern, current_poly_parts)]
              }
            }
          }

          # Prepare covariables (exclude base variable to avoid duplication)
          if (!is.null(covariables)) {
            current_covariables <- covariables[covariables != base_var]
            expanded_covariables <- unlist(lapply(current_covariables, expand_terms))
            expanded_covariables <- unique(expanded_covariables[expanded_covariables != ""])
          } else {
            expanded_covariables <- character(0)
          }

        } else {
          # Regular variable - existing logic
          # Remove spline/poly versions of current predictor
          if (!is.null(spline_terms)) {
            for (i in seq_along(spline_terms)) {
              if (spline_terms[[i]]$var == current_predictor) {
                pattern <- if (use_rms) {
                  paste0("rcs\\(", current_predictor, ",")
                } else {
                  paste0("bs\\(", current_predictor, ",")
                }
                current_spline_parts <- current_spline_parts[!grepl(pattern, current_spline_parts)]
              }
            }
          }

          if (!is.null(poly_terms)) {
            for (i in seq_along(poly_terms)) {
              if (poly_terms[[i]]$var == current_predictor) {
                pattern <- paste0("poly\\(", current_predictor, ",")
                current_poly_parts <- current_poly_parts[!grepl(pattern, current_poly_parts)]
              }
            }
          }

          # Prepare covariables (exclude current predictor to avoid duplication)
          if (!is.null(covariables)) {
            current_covariables <- covariables[covariables != current_predictor]
            expanded_covariables <- unlist(lapply(current_covariables, expand_terms))
            expanded_covariables <- unique(expanded_covariables[expanded_covariables != ""])
          } else {
            expanded_covariables <- character(0)
          }

          # Expand current predictor for interactions
          expanded_predictor <- expand_terms(current_predictor)
        }

        # Combine all terms
        all_terms <- c(expanded_predictor, expanded_covariables, current_spline_parts, current_poly_parts)
      }

      # Remove empty terms and clean
      all_terms <- all_terms[!is.na(all_terms) & all_terms != ""]
      all_terms <- unique(all_terms)

      # Build formula based on model type - ENHANCED for null model
      if (model_type == "cox") {
        if (random_intercept == "Yes") {
          if (length(all_terms) == 0) {
            formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ 1", trial_term,
                                 "+ (1 | ", random_intercept_var, ")")
          } else {
            formula_str <- paste("Surv(", time_col, ",", event_col, ") ~",
                                 paste(all_terms, collapse = " + "), trial_term,
                                 "+ (1 | ", random_intercept_var, ")")
          }
        } else {
          if (length(all_terms) == 0) {
            formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ 1", trial_term)
          } else {
            formula_str <- paste("Surv(", time_col, ",", event_col, ") ~",
                                 paste(all_terms, collapse = " + "), trial_term)
          }
        }
      } else {
        if (random_intercept == "Yes") {
          if (length(all_terms) == 0) {
            formula_str <- paste(outcome_var, "~ 1",
                                 offset_term, trial_term, random_term)
          } else {
            formula_str <- paste(outcome_var, "~",
                                 paste(all_terms, collapse = " + "),
                                 offset_term, trial_term, random_term)
          }
        } else {
          if (length(all_terms) == 0) {
            formula_str <- paste(outcome_var, "~ 1",
                                 offset_term, trial_term)
          } else {
            formula_str <- paste(outcome_var, "~",
                                 paste(all_terms, collapse = " + "),
                                 offset_term, trial_term)
          }
        }
      }
    } else {
      # Use provided formula_string
      formula_str <- formula_string
      # Add components if needed
      if (model_type == "cox") {
        if (!grepl("^Surv\\(", formula_str)) {
          formula_str <- paste("Surv(", time_col, ",", event_col, ") ~", formula_str)
        }
      } else {
        if (!grepl(paste0("^", outcome_var, "\\s*~"), formula_str)) {
          formula_str <- paste(outcome_var, "~", formula_str)
        }
      }

      if (trial_factor == "Yes" && !grepl(paste0("as\\.factor\\(", trial_col, "\\)"), formula_str)) {
        formula_str <- paste(formula_str, trial_term)
      }

      if (followup_offset == "Yes" && !grepl("offset\\(log\\(.*\\)\\)", formula_str)) {
        formula_str <- paste(formula_str, offset_term)
      }
    }

    return(formula_str)
  }

  # Function to calculate performance for a single predictor
  calculate_performance_for_predictor <- function(current_predictor) {
    formula_string_current <- build_formula_for_predictor(current_predictor)
    model_formula <- as.formula(formula_string_current)

    # Build null formula (covariables only, no predictors)
    null_formula_str <- build_formula_for_predictor("")  # Empty predictor for null model
    null_formula <- as.formula(null_formula_str)

    # Build intercept-only formula (for R2 calculation baseline)
    if (model_type == "cox") {
      if (random_intercept == "Yes") {
        intercept_formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ 1 + (1 | ", random_intercept_var, ")")
      } else {
        intercept_formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ 1")
      }
    } else {
      offset_term_intercept <- if (followup_offset == "Yes") paste("+ offset(log(", followup_col, "))") else ""
      if (random_intercept == "Yes") {
        intercept_formula_str <- paste(outcome_var, "~ 1", offset_term_intercept, "+ (1 | ", random_intercept_var, ")")
      } else {
        intercept_formula_str <- paste(outcome_var, "~ 1", offset_term_intercept)
      }
    }
    intercept_formula <- as.formula(intercept_formula_str)

    if (verbose) {
      # Determine display name for RCS terms and null model
      display_name <- if (is.na(current_predictor) || current_predictor == "") {
        "No_predictor (null model)"
      } else if (grepl("^rcs\\(|^bs\\(|^poly\\(", current_predictor)) {
        current_predictor  # Show the full RCS term
      } else {
        current_predictor
      }

      message(paste("Processing performance for predictor:", display_name))
    }

    # Pre-allocate result vectors
    logL_full <- numeric(imp_n)
    logL_null <- numeric(imp_n)
    logL_intercept <- numeric(imp_n)
    df_values <- numeric(imp_n)
    c_index_values <- numeric(imp_n)
    rss_values <- numeric(imp_n)
    ae_values <- numeric(imp_n)
    r2_values <- numeric(imp_n)

    subject_n <- NA

    # Helper functions for model fitting
    fit_model <- function(formula, data_subset) {
      if (random_intercept == "Yes") {
        switch(model_type,
               "lm" = lme4::lmer(formula, data = data_subset),
               "bin" = lme4::glmer(formula, family = binomial(), data = data_subset),
               "poisson" = lme4::glmer(formula, family = poisson(), data = data_subset),
               "nb" = {
                 if (requireNamespace("glmmTMB", quietly = TRUE)) {
                   glmmTMB::glmmTMB(formula, family = glmmTMB::nbinom2, data = data_subset)
                 } else {
                   lme4::glmer(formula, family = poisson(), data = data_subset)
                 }
               },
               "cox" = coxme::coxme(formula, data = data_subset),
               "gamma" = {
                 if (requireNamespace("glmmTMB", quietly = TRUE)) {
                   glmmTMB::glmmTMB(formula, family = Gamma, data = data_subset)
                 } else {
                   lme4::glmer(formula, family = Gamma, data = data_subset)
                 }
               },
               stop("Unsupported model type for random effects.")
        )
      } else {
        switch(model_type,
               "nb" = MASS::glm.nb(formula, data = data_subset),
               "lm" = glm(formula, family = gaussian(), data = data_subset),
               "bin" = glm(formula, family = binomial(), data = data_subset),
               "poisson" = glm(formula, family = poisson(), data = data_subset),
               "gamma" = glm(formula, family = Gamma(), data = data_subset),
               "quasipoisson" = glm(formula, family = quasipoisson(), data = data_subset),
               "quasibinomial" = glm(formula, family = quasibinomial(), data = data_subset),
               "cox" = survival::coxph(formula, data = data_subset),
               stop("Unsupported model type.")
        )
      }
    }

    # Extract log-likelihood
    extract_logLik <- function(model) {
      if (inherits(model, "merMod") || inherits(model, "glmmTMB")) {
        return(as.numeric(logLik(model)))
      } else if (inherits(model, "coxme")) {
        return(as.numeric(model$loglik[2]))
      } else {
        return(as.numeric(logLik(model)))
      }
    }

    # Get model degrees of freedom
    get_df <- function(model) {
      if (inherits(model, "merMod")) {
        return(length(lme4::fixef(model)))
      } else if (inherits(model, "glmmTMB")) {
        return(length(glmmTMB::fixef(model)$cond))
      } else if (inherits(model, "coxme")) {
        return(length(coxme::fixef(model)))
      } else {
        return(length(coef(model)))
      }
    }

    # Loop through imputations
    for (i in seq_along(actual_imps)) {
      imp_val <- actual_imps[i]
      data_subset <- dplyr::filter(data, !!rlang::sym(imp_col) == imp_val)

      # Fit models with error handling
      model <- tryCatch({
        fit_model(model_formula, data_subset)
      }, error = function(e) {
        if (verbose) warning(paste("Error fitting full model for imputation", i, ":", conditionMessage(e)))
        return(NULL)
      })

      null_model <- tryCatch({
        fit_model(null_formula, data_subset)
      }, error = function(e) {
        if (verbose) warning(paste("Error fitting null model for imputation", i, ":", conditionMessage(e)))
        return(NULL)
      })

      # Fit intercept-only model for R2 baseline
      intercept_model <- tryCatch({
        fit_model(intercept_formula, data_subset)
      }, error = function(e) {
        if (verbose) warning(paste("Error fitting intercept model for imputation", i, ":", conditionMessage(e)))
        return(NULL)
      })

      # Skip if any model failed
      if (is.null(model) || is.null(null_model) || is.null(intercept_model)) {
        logL_full[i] <- NA
        logL_null[i] <- NA
        logL_intercept[i] <- NA
        df_values[i] <- NA
        c_index_values[i] <- NA
        rss_values[i] <- NA
        ae_values[i] <- NA
        r2_values[i] <- NA
        next
      }

      # Get subject count from first successful model
      if (is.na(subject_n)) {
        if (inherits(model, "glm") || inherits(model, "negbin")) {
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
          subject_n <- nrow(data_subset)
        }
      }

      # Store metrics
      logL_full[i] <- extract_logLik(model)
      logL_null[i] <- extract_logLik(null_model)
      logL_intercept[i] <- extract_logLik(intercept_model)
      df_values[i] <- get_df(model)

      # Get subject count from first successful model
      if (is.na(subject_n)) {
        if (inherits(model, "glm") || inherits(model, "negbin")) {
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
          subject_n <- nrow(data_subset)
        }
      }

      # Store metrics
      logL_full[i] <- extract_logLik(model)
      logL_null[i] <- extract_logLik(null_model)
      logL_intercept[i] <- extract_logLik(intercept_model)
      df_values[i] <- get_df(model)

      # Calculate performance metrics
      if (model_type == "cox") {
        c_index_values[i] <- tryCatch(
          survival::concordance(model)$concordance,
          error = function(e) NA
        )
        rss_values[i] <- NA
        ae_values[i] <- NA
        r2_values[i] <- NA
      } else {
        # Get predictions and observed values with better error handling
        tryCatch({
          pred_values <- fitted(model)
          obs_values <- data_subset[[outcome_var]]

          # Debug: Check dimensions
          if (verbose) {
            cat("Debug - Imputation", i, ": pred_length =", length(pred_values),
                ", obs_length =", length(obs_values), "\n")
          }

          # Ensure same length
          min_length <- min(length(pred_values), length(obs_values))
          if (length(pred_values) != length(obs_values)) {
            warning(paste("Length mismatch in imputation", i,
                          ": pred =", length(pred_values),
                          ", obs =", length(obs_values),
                          ". Using first", min_length, "values."))
            pred_values <- pred_values[1:min_length]
            obs_values <- obs_values[1:min_length]
          }

          # Calculate metrics with error handling
          c_index_values[i] <- tryCatch({
            if (requireNamespace("Hmisc", quietly = TRUE)) {
              Hmisc::rcorr.cens(pred_values, obs_values)[1]
            } else {
              cor(pred_values, obs_values, use = "complete.obs")
            }
          }, error = function(e) {
            if (verbose) warning(paste("Error calculating C-index for imputation", i, ":", e$message))
            NA
          })

          rss_values[i] <- tryCatch({
            sum((obs_values - pred_values)^2, na.rm = TRUE)
          }, error = function(e) {
            if (verbose) warning(paste("Error calculating RSS for imputation", i, ":", e$message))
            NA
          })

          ae_values[i] <- tryCatch({
            sum(abs(obs_values - pred_values), na.rm = TRUE)
          }, error = function(e) {
            if (verbose) warning(paste("Error calculating AE for imputation", i, ":", e$message))
            NA
          })

          r2_values[i] <- tryCatch({
            ss_res <- sum((obs_values - pred_values)^2, na.rm = TRUE)
            ss_tot <- sum((obs_values - mean(obs_values, na.rm = TRUE))^2, na.rm = TRUE)
            1 - (ss_res / ss_tot)
          }, error = function(e) {
            if (verbose) warning(paste("Error calculating R2 for imputation", i, ":", e$message))
            NA
          })

        }, error = function(e) {
          if (verbose) warning(paste("Error in performance calculation for imputation", i, ":", e$message))
          c_index_values[i] <- NA
          rss_values[i] <- NA
          ae_values[i] <- NA
          r2_values[i] <- NA
        })
      }
    }

    # Pool performance metrics
    valid_imps <- !is.na(logL_full)
    n_valid <- sum(valid_imps)

    if (verbose) {
      cat("Valid imputations:", n_valid, "out of", imp_n, "\n")
    }

    if (n_valid == 0) {
      warning(paste("All models failed to converge for predictor:", current_predictor))
      # Return a basic result structure instead of NULL
      return(list(
        Model_Formula = formula_string_current,
        Model_Type = model_type,
        Predictor_Tested = ifelse(is.na(current_predictor) || current_predictor == "", "No_predictor", current_predictor),
        Number_Subjects = nrow(data[data[[imp_col]] == actual_imps[1], ]),
        Number_Imputation = imp_n,
        Number_Valid_Imputation = 0,
        df = NA,
        B_logL = NA,
        T_logL = NA,
        logL = NA,
        R2 = NA,
        R2_Method = "Failed to calculate",
        AIC = NA,
        AICc = NA,
        BIC = NA,
        BICc = NA,
        AIC_method = "Failed to calculate",
        C_Index = NA,
        RMSE = NA,
        MAE = NA
      ))
    }

    if (is.na(subject_n)) {
      subject_n <- nrow(data[data[[imp_col]] == actual_imps[1], ])
    }

    # Calculate pooled metrics with better error handling
    K <- mean(df_values[valid_imps], na.rm = TRUE)
    if (is.na(K) || is.infinite(K)) K <- 1  # Fallback

    B_logL <- var(logL_full[valid_imps], na.rm = TRUE)
    if (is.na(B_logL) || is.infinite(B_logL)) B_logL <- 0  # Fallback

    W_logL <- K
    T_logL <- W_logL + (1 + 1/n_valid) * B_logL
    logL_pooled <- mean(logL_full[valid_imps], na.rm = TRUE)

    # R2 calculation - corrected logic
    if (tolower(r2_method) == "pseudo") {
      if (model_type == "cox") {
        R2_pooled <- NA
        R2_method_used <- "NA (Cox model - no standard R2)"
      } else {
        # All models compared against intercept-only model
        intercept_logL_pooled <- mean(logL_intercept[valid_imps], na.rm = TRUE)
        if (is.na(intercept_logL_pooled) || intercept_logL_pooled == 0) {
          R2_pooled <- NA
        } else {
          R2_pooled <- 1 - (logL_pooled / intercept_logL_pooled)
        }
        R2_method_used <- "Pseudo R2 (compared to intercept-only model)"
      }
    } else {
      R2_pooled <- mean(r2_values[valid_imps], na.rm = TRUE)
      R2_method_used <- "Pooled R2 (average of per-imputation R2)"
    }

    # AIC calculation
    if (tolower(aic_method) == "rubin") {
      AIC_pooled <- -2 * logL_pooled + 2 * T_logL
      AICc_pooled <- AIC_pooled + (2 * T_logL * (T_logL + 1)) / (subject_n - T_logL - 1)
      BIC_pooled <- -2 * logL_pooled + T_logL * log(subject_n)
      BICc_pooled <- BIC_pooled + (T_logL * (T_logL + 1)) / (subject_n - T_logL - 1)
      AIC_method_used <- "Corrected for between-imputation variance"
    } else {
      AIC_pooled <- -2 * logL_pooled + 2 * K
      AICc_pooled <- AIC_pooled + (2 * K * (K + 1)) / (subject_n - K - 1)
      BIC_pooled <- -2 * logL_pooled + K * log(subject_n)
      BICc_pooled <- BIC_pooled + (K * (K + 1)) / (subject_n - K - 1)
      AIC_method_used <- "Based on pooled logLik and degrees of freedom"
    }

    # Other pooled metrics
    C_index_pooled <- mean(c_index_values[valid_imps], na.rm = TRUE)

    if (model_type != "cox") {
      # RMSE calculation with proper error handling
      valid_rss <- rss_values[valid_imps & !is.na(rss_values)]
      if (length(valid_rss) > 0) {
        W <- mean(valid_rss, na.rm = TRUE) / subject_n
        B <- var(valid_rss / subject_n, na.rm = TRUE)
        if (is.na(B)) B <- 0
        T_val <- W + (1 + 1/length(valid_rss)) * B
        RMSE_pooled <- sqrt(T_val)
      } else {
        RMSE_pooled <- NA
      }

      # MAE calculation with proper error handling
      valid_ae <- ae_values[valid_imps & !is.na(ae_values)]
      if (length(valid_ae) > 0) {
        W_mae <- mean(valid_ae, na.rm = TRUE) / subject_n
        B_mae <- var(valid_ae / subject_n, na.rm = TRUE)
        if (is.na(B_mae)) B_mae <- 0
        T_mae <- W_mae + (1 + 1/length(valid_ae)) * B_mae
        MAE_pooled <- T_mae
      } else {
        MAE_pooled <- NA
      }
    } else {
      RMSE_pooled <- NA
      MAE_pooled <- NA
    }

    # Determine predictor name for output - MODIFIED for RCS support and null model
    predictor_tested <- if (is.na(current_predictor) || current_predictor == "") {
      "No_predictor"
    } else if (grepl("^rcs\\(|^bs\\(|^poly\\(", current_predictor)) {
      current_predictor  # Keep the full RCS/spline/poly term
    } else {
      current_predictor
    }

    if (verbose) {
      cat("Final calculations:\n")
      cat("  logL_pooled:", logL_pooled, "\n")
      cat("  AIC_pooled:", AIC_pooled, "\n")
      cat("  R2_pooled:", R2_pooled, "\n")
      cat("  C_index_pooled:", C_index_pooled, "\n")
      cat("  predictor_tested:", predictor_tested, "\n")
    }

    # Return performance metrics
    result <- tryCatch({
      list(
        Model_Formula = formula_string_current,
        Model_Type = model_type,
        Predictor_Tested = predictor_tested,
        Number_Subjects = subject_n,
        Number_Imputation = imp_n,
        Number_Valid_Imputation = n_valid,
        df = K,
        B_logL = B_logL,
        T_logL = T_logL,
        logL = logL_pooled,
        R2 = R2_pooled,
        R2_Method = R2_method_used,
        AIC = AIC_pooled,
        AICc = AICc_pooled,
        BIC = BIC_pooled,
        BICc = BICc_pooled,
        AIC_method = AIC_method_used,
        C_Index = C_index_pooled,
        RMSE = RMSE_pooled,
        MAE = MAE_pooled
      )
    }, error = function(e) {
      if (verbose) cat("Error creating result list:", e$message, "\n")
      return(NULL)
    })

    if (verbose) {
      cat("Result created successfully:", !is.null(result), "\n")
      if (!is.null(result)) {
        cat("Result class:", class(result), "\n")
        cat("Result length:", length(result), "\n")
      }
    }

    # Add attributes
    if (!is.null(result)) {
      tryCatch({
        attr(result, "predictor_tested") <- predictor_tested
        attr(result, "formula") <- formula_string_current
        attr(result, "covariables") <- covariables
      }, error = function(e) {
        if (verbose) cat("Error adding attributes:", e$message, "\n")
      })
    }

    return(result)
  }

  # Main execution: calculate performance for each predictor
  results_list <- vector("list", length(predictor_vars))

  # Set proper names for the results list
  result_names <- character(length(predictor_vars))
  for (i in seq_along(predictor_vars)) {
    current_pred <- predictor_vars[i]
    if (is.na(current_pred) || current_pred == "" || is.null(current_pred)) {
      result_names[i] <- "No_predictor"
    } else {
      result_names[i] <- current_pred
    }
  }
  names(results_list) <- result_names

  for (i in seq_along(predictor_vars)) {
    current_predictor <- predictor_vars[i]

    if (verbose) {
      cat("Processing predictor", i, "of", length(predictor_vars), "\n")
    }

    if (is.na(current_predictor) || current_predictor == "" || is.null(current_predictor)) {
      if (verbose) message("Calculating performance for null model")
      result <- calculate_performance_for_predictor(current_predictor)
      results_list[["No_predictor"]] <- result

      if (verbose) {
        cat("Assigned result to No_predictor slot:", !is.null(result), "\n")
        cat("Results list names:", paste(names(results_list), collapse = ", "), "\n")
        cat("Results list length:", length(results_list), "\n")
      }
    } else {
      if (verbose) {
        display_name <- if (grepl("^rcs\\(|^bs\\(|^poly\\(", current_predictor)) {
          current_predictor
        } else {
          current_predictor
        }
        message(paste("Calculating performance for predictor:", display_name))
      }
      result <- calculate_performance_for_predictor(current_predictor)
      results_list[[i]] <- result

      if (verbose) {
        cat("Assigned result to slot", i, ":", !is.null(result), "\n")
      }
    }
  }

  # Create combined results for multiple predictors
  if (length(predictor_vars) > 1) {
    combined_performance <- data.frame()

    for (pred_name in names(results_list)) {
      result <- results_list[[pred_name]]
      if (!is.null(result)) {
        performance_row <- data.frame(
          Predictor = result$Predictor_Tested,
          AIC = result$AIC,
          AICc = result$AICc,
          BIC = result$BIC,
          BICc = result$BICc,
          R2 = result$R2,
          C_Index = result$C_Index,
          RMSE = result$RMSE,
          MAE = result$MAE,
          df = result$df,
          logL = result$logL,
          stringsAsFactors = FALSE
        )
        combined_performance <- rbind(combined_performance, performance_row)
      }
    }

    # Reset row names
    rownames(combined_performance) <- NULL
  }

  # Return results
  if (verbose) {
    cat("Final return section:\n")
    cat("Length of predictor_vars:", length(predictor_vars), "\n")
    cat("Length of results_list:", length(results_list), "\n")
    cat("Names of results_list:", paste(names(results_list), collapse = ", "), "\n")
    cat("Non-null results:", sum(sapply(results_list, function(x) !is.null(x))), "\n")
  }

  if (length(predictor_vars) == 1) {
    # Single predictor case
    result_to_return <- results_list[[1]]

    if (verbose) {
      cat("Returning single result:\n")
      cat("  Result is null:", is.null(result_to_return), "\n")
      if (!is.null(result_to_return)) {
        cat("  Result class:", class(result_to_return), "\n")
        cat("  Result names:", paste(names(result_to_return), collapse = ", "), "\n")
      }
    }

    return(result_to_return)
  } else {
    # Multiple predictors case
    if (verbose) {
      cat("Returning multiple results as list\n")
    }

    class(results_list) <- c("MI_performance_multi", "list")
    attr(results_list, "covariables") <- covariables
    attr(results_list, "predictors") <- predictor_vars
    attr(results_list, "model_type") <- model_type
    attr(results_list, "n_predictors") <- length(predictor_vars)
    attr(results_list, "combined_performance") <- combined_performance
    return(results_list)
  }
}
