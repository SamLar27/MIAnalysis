#' Compute Estimates for Multiple Imputed Data
#'
#' This function fits statistical models to multiply imputed datasets and pools the results using Rubin's Rules.
#' It supports various regression models, including negative binomial, logistic, Poisson, linear, and Cox regression.
#' The function also supports modeling with random intercepts using `glmmTMB`, `lme4`, and `coxme`, and can handle interaction terms, restricted cubic splines (rcs), and polynomial terms.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for GLM models.
#' @param predictor_vars A vector of predictor variables to be tested in the model. Each predictor will be included individually along with all covariables.
#' @param covariables A vector of covariable names to be included as fixed effects in all models. These are always included regardless of which predictor is being tested.
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param imp_n Number of imputations (if NULL, it is detected automatically).
#' @param model_type The model type to be used:
#'   \itemize{
#'     \item "nb" for Negative Binomial regression (using `MASS::glm.nb` or `glmmTMB` if random intercepts are used)
#'     \item "lm" for Linear regression
#'     \item "bin" for Logistic regression (binomial family)
#'     \item "poisson" for Poisson regression
#'     \item "gamma" for Gamma regression
#'     \item "quasipoisson" for Overdispersed Poisson regression
#'     \item "quasibinomial" for Overdispersed logistic regression
#'     \item "cox" for Cox Proportional Hazards regression (using `survival` or `coxme`)
#'   }
#' @param followup_offset Whether to include an offset for follow-up duration ("Yes" or "No").
#' @param followup_col The column representing follow-up duration (required if `followup_offset = "Yes"`).
#' @param trial_factor Whether to include trial as a fixed factor ("Yes" or "No").
#' @param trial_col The column used for trial factor (required if `trial_factor = "Yes"`).
#' @param time_col The time variable for Cox regression (only required if `model_type = "cox"`).
#' @param event_col The event variable for Cox regression (only required if `model_type = "cox"`).
#' @param formula_string Optional custom formula string. If provided, overrides the formula created from predictor_vars and covariables.
#' @param highlight_interactions Whether to flag interaction terms in the output (default: TRUE).
#' @param spline_terms A named list to specify variables to be modeled with restricted cubic splines. Each element should be a list with 'var' (variable name) and 'knots' (number of knots or explicit knot positions).
#' @param poly_terms A named list to specify variables to be modeled with polynomial terms. Each element should be a list with 'var' (variable name) and 'degree' (2 for quadratic, 3 for cubic).
#' @param include_spline_terms Whether to include individual spline components in the output (default: FALSE).
#' @param include_poly_terms Whether to include individual polynomial components in the output (default: TRUE).
#' @param random_intercept Whether to include a random intercept in the model ("Yes" or "No").
#' @param random_intercept_var The grouping variable for the random intercept (required if `random_intercept = "Yes"`).
#'
#' @return A list containing model results for each predictor variable. Each element is a data frame with model estimates, standard errors, confidence intervals, p-values, and exponentiated estimates (if applicable).
#'   Additional attributes for each result include:
#'   \itemize{
#'     \item \code{predictor_tested}: the predictor variable that was tested in this model
#'     \item \code{formula}: the model formula used
#'     \item \code{has_random_effects}: logical flag indicating if random effects were used
#'     \item \code{random_intercept_var}: the grouping variable used for the random intercept
#'     \item \code{variance_components}: the estimated variance of the random intercept (if applicable)
#'     \item \code{ICC}: intraclass correlation coefficient (for linear mixed models only)
#'     \item \code{model_type_used}: the model type used internally
#'   }
#'
#' @import dplyr MASS mice survival splines
#' @importFrom lme4 lmer glmer fixef VarCorr
#' @importFrom glmmTMB glmmTMB fixef VarCorr
#' @importFrom coxme coxme
#' @export
MI_estimates <- function(data,
                         outcome_var,
                         predictor_vars,
                         covariables = NULL,
                         imp_col = ".imp",
                         imp_n = NULL,
                         model_type = "nb",
                         followup_offset = "No",  # Default to "No"
                         followup_col = NULL,
                         trial_factor = "No",  # Default to "No"
                         trial_col = NULL,
                         time_col = NULL,
                         event_col = NULL,
                         formula_string = NULL,
                         highlight_interactions = TRUE,
                         spline_terms = NULL,
                         poly_terms = NULL,
                         include_spline_terms = FALSE,
                         include_poly_terms = TRUE,
                         random_intercept = "No",  # Default to "No"
                         random_intercept_var = NULL) {

  # Load required packages
  require(dplyr)
  require(MASS)         # For Negative Binomial
  require(mice)         # For multiple imputations
  require(rlang)        # For formula manipulation
  require(survival)     # For Cox regression
  require(splines)      # For basic spline functions

  # Load packages for mixed models
  if (random_intercept == "Yes") {
    require(lme4)       # For mixed models
    require(mitools)    # For handling multiple imputations with mixed models

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
        require(glmmTMB)  # For negative binomial mixed models
      }
    }
  }

  # Check if rms package is available if spline_terms is provided
  if (!is.null(spline_terms)) {
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

  # Determine number of imputations if not provided
  actual_imps <- sort(unique(data[[imp_col]]))
  if (is.null(imp_n)) {
    imp_n <- length(actual_imps)
  } else {
    if (imp_n != length(actual_imps)) {
      warning("imp_n does not match the number of unique imputations found in the data. Using only present imputations.")
      imp_n <- length(actual_imps)
    }
  }

  # Create imputation list using only present imputations
  implist <- vector("list", length(actual_imps))
  for (i in seq_along(actual_imps)) {
    imp_val <- actual_imps[i]
    data_i <- dplyr::filter(data, !!rlang::sym(imp_col) == imp_val)
    if (nrow(data_i) == 0) stop(paste("No data for imputation", imp_val))
    implist[[i]] <- data_i
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

  # Validate polynomial terms if provided
  if (!is.null(poly_terms)) {
    if (!is.list(poly_terms)) {
      stop("poly_terms must be a list")
    }

    for (i in seq_along(poly_terms)) {
      if (!is.list(poly_terms[[i]]) || is.null(poly_terms[[i]]$var)) {
        stop("Each element in poly_terms must be a list with at least a 'var' element")
      }

      if (!poly_terms[[i]]$var %in% names(data)) {
        stop(paste("Polynomial variable", poly_terms[[i]]$var, "not found in data"))
      }

      if (is.null(poly_terms[[i]]$degree)) {
        # Default to quadratic if not specified
        poly_terms[[i]]$degree <- 2
      } else if (!is.numeric(poly_terms[[i]]$degree) || !(poly_terms[[i]]$degree %in% c(2, 3))) {
        stop("Degree must be either 2 (quadratic) or 3 (cubic)")
      }
    }
  }

  # Validate predictor_vars and covariables
  if (is.null(formula_string)) {
    # Check that predictor_vars exist in data
    if (!all(predictor_vars %in% names(data))) {
      missing_vars <- predictor_vars[!predictor_vars %in% names(data)]
      stop(paste("Predictor variables not found in data:", paste(missing_vars, collapse = ", ")))
    }

    # Check that covariables exist in data (if provided)
    if (!is.null(covariables)) {
      if (!all(covariables %in% names(data))) {
        missing_vars <- covariables[!covariables %in% names(data)]
        stop(paste("Covariables not found in data:", paste(missing_vars, collapse = ", ")))
      }
    }
  }

  if (!is.null(followup_col) && !followup_col %in% names(data)) stop("followup_col not found in data.")
  if (!is.null(trial_col) && !trial_col %in% names(data)) stop("trial_col not found in data.")
  if (!is.null(time_col) && !time_col %in% names(data)) stop("time_col not found in data.")
  if (!is.null(event_col) && !event_col %in% names(data)) stop("event_col not found in data.")

  # Validate followup_offset input
  if (!followup_offset %in% c("Yes", "No")) stop("followup_offset must be either 'Yes' or 'No'.")

  if (followup_offset == "Yes" && is.null(followup_col)) {
    stop("If followup_offset = 'Yes', followup_col must be provided.")
  }

  if (!is.null(followup_col) && any(data[[followup_col]] <= 0, na.rm = TRUE)) {
    stop("followup_col must be strictly positive for offset.")
  }

  # Validate trial_factor input
  if (!trial_factor %in% c("Yes", "No")) stop("trial_factor must be either 'Yes' or 'No'.")

  if (trial_factor == "Yes" && is.null(trial_col)) {
    stop("If trial_factor = 'Yes', trial_col must be provided.")
  }

  # Validate random_intercept input
  if (!random_intercept %in% c("Yes", "No")) stop("random_intercept must be either 'Yes' or 'No'.")

  if (random_intercept == "Yes" && is.null(random_intercept_var)) {
    stop("If random_intercept = 'Yes', random_intercept_var must be provided.")
  }

  if (!is.null(random_intercept_var) && !random_intercept_var %in% names(data)) {
    stop("random_intercept_var not found in data.")
  }

  # Determine number of imputations if not provided
  if (is.null(imp_n)) {
    imp_n <- length(unique(data[[imp_col]]))
    if (imp_n < 2) stop("At least 2 imputations required.")
  }

  # Helper function to expand interaction terms
  expand_terms <- function(term) {
    if (grepl("^rcs\\(|^bs\\(", term)) {
      # If it is already a spline function, don't expand
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
    # Process spline terms if provided
    spline_formula_parts <- character(0)
    if (!is.null(spline_terms)) {
      for (i in seq_along(spline_terms)) {
        var_name <- spline_terms[[i]]$var
        knots <- spline_terms[[i]]$knots

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
    if (!is.null(poly_terms)) {
      for (i in seq_along(poly_terms)) {
        var_name <- poly_terms[[i]]$var
        degree <- poly_terms[[i]]$degree

        # Create the polynomial term using poly() function
        # The raw=TRUE option creates simple polynomials rather than orthogonal ones
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

  # Helper function to build formula for a given predictor
  build_formula_for_predictor <- function(current_predictor) {
    # Construct formula with optional terms
    trial_term <- if (trial_factor == "Yes")
      paste("+ as.factor(", trial_col, ")") else ""

    offset_term <- if (followup_offset == "Yes")
      paste("+ offset(log(", followup_col, "))") else ""

    # Random intercept term (for formula construction)
    random_term <- if (random_intercept == "Yes")
      paste("+ (1 | ", random_intercept_var, ")") else ""

    if (is.null(formula_string)) {
      # Remove current predictor from spline_terms and poly_terms if present
      current_spline_parts <- spline_formula_parts
      current_poly_parts <- poly_formula_parts

      # Remove spline/poly versions of current predictor
      if (!is.null(spline_terms)) {
        for (i in seq_along(spline_terms)) {
          if (spline_terms[[i]]$var == current_predictor) {
            # Remove this spline term
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
            # Remove this polynomial term
            pattern <- paste0("poly\\(", current_predictor, ",")
            current_poly_parts <- current_poly_parts[!grepl(pattern, current_poly_parts)]
          }
        }
      }

      # Prepare covariables (exclude current predictor to avoid duplication)
      if (!is.null(covariables)) {
        current_covariables <- covariables[covariables != current_predictor]
        # Expand covariables for interactions
        expanded_covariables <- unlist(lapply(current_covariables, expand_terms))
        expanded_covariables <- unique(expanded_covariables)
      } else {
        expanded_covariables <- character(0)
      }

      # Expand current predictor for interactions
      expanded_predictor <- expand_terms(current_predictor)

      # Combine all terms: predictor + covariables + special terms
      all_terms <- c(expanded_predictor, expanded_covariables, current_spline_parts, current_poly_parts)

      if (model_type == "cox") {
        if (is.null(time_col) || is.null(event_col)) {
          stop("For Cox regression, time_col and event_col must be provided.")
        }

        if (random_intercept == "Yes") {
          # For Cox mixed model with coxme
          formula_str <- paste("Surv(", time_col, ",", event_col, ") ~",
                               paste(all_terms, collapse = " + "), trial_term,
                               "+ (1 | ", random_intercept_var, ")")
        } else {
          # Regular Cox model
          formula_str <- paste("Surv(", time_col, ",", event_col, ") ~",
                               paste(all_terms, collapse = " + "), trial_term)
        }
      } else {
        # For GLM/mixed models
        if (random_intercept == "Yes") {
          # For mixed models, we add random effect to fixed effects formula
          formula_str <- paste(outcome_var, "~",
                               paste(all_terms, collapse = " + "),
                               offset_term, trial_term, random_term)
        } else {
          # Regular GLM formula
          formula_str <- paste(outcome_var, "~",
                               paste(all_terms, collapse = " + "),
                               offset_term, trial_term)
        }
      }
    } else {
      # Use the provided formula_string (assuming it includes the current predictor appropriately)
      formula_str <- formula_string

      if (model_type == "cox") {
        if (is.null(time_col) || is.null(event_col)) {
          stop("For Cox regression, time_col and event_col must be provided.")
        }

        # Check if formula already has Surv() term
        if (!grepl("^Surv\\(", formula_str)) {
          formula_str <- paste("Surv(", time_col, ",", event_col, ") ~", formula_str)
        }
      } else {
        # Check if formula already has outcome variable
        if (!grepl(paste0("^", outcome_var, "\\s*~"), formula_str)) {
          formula_str <- paste(outcome_var, "~", formula_str)
        }
      }

      # Add trial term if needed and not already in formula
      if (trial_factor == "Yes" && !grepl(paste0("as\\.factor\\(", trial_col, "\\)"), formula_str)) {
        formula_str <- paste(formula_str, trial_term)
      }

      # Add offset if needed and not already in formula
      if (followup_offset == "Yes" && !grepl("offset\\(log\\(.*\\)\\)", formula_str)) {
        formula_str <- paste(formula_str, offset_term)
      }
    }

    return(formula_str)
  }

  # Function to fit model for a single predictor (existing model fitting logic)
  fit_model_for_predictor <- function(current_predictor) {
    formula_string_current <- build_formula_for_predictor(current_predictor)
    model_formula <- as.formula(formula_string_current)

    # Detect interaction terms
    interaction_terms <- character(0)
    if (grepl(":", formula_string_current)) {
      # Extract terms after the "~"
      terms_part <- strsplit(formula_string_current, "~")[[1]][2]
      # Split by "+"
      terms <- strsplit(terms_part, "\\+")[[1]]
      # Trim whitespace
      terms <- trimws(terms)
      # Find terms containing ":"
      interaction_terms <- terms[grepl(":", terms)]
    }

    # Detect spline terms
    spline_terms_detected <- character(0)
    if (grepl("rcs\\(|bs\\(", formula_string_current)) {
      # Extract terms after the "~"
      terms_part <- strsplit(formula_string_current, "~")[[1]][2]
      # Split by "+"
      terms <- strsplit(terms_part, "\\+")[[1]]
      # Trim whitespace
      terms <- trimws(terms)
      # Find terms containing "rcs(" or "bs("
      spline_terms_detected <- terms[grepl("rcs\\(|bs\\(", terms)]
    }

    # Detect polynomial terms
    poly_terms_detected <- character(0)
    if (grepl("poly\\(", formula_string_current)) {
      # Extract terms after the "~"
      terms_part <- strsplit(formula_string_current, "~")[[1]][2]
      # Split by "+"
      terms <- strsplit(terms_part, "\\+")[[1]]
      # Trim whitespace
      terms <- trimws(terms)
      # Find terms containing "poly("
      poly_terms_detected <- terms[grepl("poly\\(", terms)]
    }

    # Fit models to each imputed dataset (using existing logic)
    if (random_intercept == "Yes") {
      # Approach models with random effects
      models_list <- vector("list", imp_n)

      # Fit the model to each imputation
      for (i in seq_along(actual_imps)) {
        data_i <- implist[[i]]

        # Fit model based on type
        if (model_type == "lm") {
          models_list[[i]] <- lme4::lmer(model_formula, data = data_i)
        } else if (model_type == "bin") {
          models_list[[i]] <- lme4::glmer(model_formula, family = binomial(), data = data_i)
        } else if (model_type == "poisson" || model_type == "quasipoisson") {
          models_list[[i]] <- lme4::glmer(model_formula, family = poisson(), data = data_i)
        } else if (model_type == "nb") {
          if (requireNamespace("glmmTMB", quietly = TRUE)) {
            models_list[[i]] <- glmmTMB::glmmTMB(model_formula, family = nbinom2, data = data_i)
          } else {
            # Fallback to Poisson if glmmTMB not available
            warning("Using Poisson instead of Negative Binomial. Install glmmTMB for proper NB mixed models.")
            models_list[[i]] <- lme4::glmer(model_formula, family = poisson(), data = data_i)
          }
        } else if (model_type == "cox") {
          models_list[[i]] <- coxme::coxme(model_formula, data = data_i)
        } else if (model_type == "gamma") {
          if (requireNamespace("glmmTMB", quietly = TRUE)) {
            models_list[[i]] <- glmmTMB::glmmTMB(model_formula, family = Gamma, data = data_i)
          } else {
            warning("Gamma family with random effects may not be stable in lme4.")
            models_list[[i]] <- lme4::glmer(model_formula, family = Gamma, data = data_i)
          }
        } else {
          stop("Unsupported model type for random effects.")
        }
      }

      # Extract coefficients
      coefs <- lapply(models_list, function(m) {
        if (model_type == "cox") {
          # coxme coefficients
          c_est <- fixef(m)
          c_se <- sqrt(diag(vcov(m)))
        } else if (inherits(m, "glmmTMB")) {
          # glmmTMB coefficients
          c_est <- fixef(m)$cond
          v <- vcov(m)$cond
          c_se <- sqrt(diag(v))
        } else {
          # lme4 coefficients - improved handling
          c_est <- fixef(m)

          # More robust standard error extraction for lmer/glmer objects
          if (inherits(m, "lmerMod")) {
            # For linear mixed models, use summary to get standard errors
            model_summary <- summary(m)
            c_se <- model_summary$coefficients[, "Std. Error"]
          } else {
            # For other lme4 models (glmer), try vcov approach with error handling
            tryCatch({
              vcov_matrix <- vcov(m)
              c_se <- sqrt(diag(vcov_matrix))
            }, error = function(e) {
              # Fallback to summary method
              model_summary <- summary(m)
              c_se <- model_summary$coefficients[, "Std. Error"]
            })
          }
        }

        data.frame(
          term = names(c_est),
          estimate = c_est,
          std.error = c_se,
          stringsAsFactors = FALSE
        )
      })
      # Pool the results using Rubin's rules manually
      # 1. Combine estimates
      all_terms <- unique(unlist(lapply(coefs, function(df) df$term)))

      # 2. Calculate pooled estimates
      pooled_results <- data.frame(term = all_terms,
                                   estimate = NA,
                                   std.error = NA,
                                   stringsAsFactors = FALSE)

      for (term in all_terms) {
        # Get estimates and standard errors for this term from all models
        term_ests <- sapply(coefs, function(df) {
          if (term %in% df$term) df$estimate[df$term == term] else NA
        })
        term_ests <- term_ests[!is.na(term_ests)]

        term_ses <- sapply(coefs, function(df) {
          if (term %in% df$term) df$std.error[df$term == term] else NA
        })
        term_ses <- term_ses[!is.na(term_ses)]

        if (length(term_ests) < 2) next

        # Rubin's rules for pooling
        pooled_est <- mean(term_ests)

        # Within-imputation variance
        w_var <- mean(term_ses^2)

        # Between-imputation variance
        b_var <- sum((term_ests - pooled_est)^2) / (length(term_ests) - 1)

        # Total variance
        total_var <- w_var + (1 + 1/length(term_ests)) * b_var
        pooled_se <- sqrt(total_var)

        # Update results
        pooled_results$estimate[pooled_results$term == term] <- pooled_est
        pooled_results$std.error[pooled_results$term == term] <- pooled_se
      }

      # Calculate confidence intervals and p-values
      pooled_results <- pooled_results %>%
        mutate(
          `2.5 %` = estimate - 1.96 * std.error,
          `97.5 %` = estimate + 1.96 * std.error,
          # More numerically stable p-value calculation
          p.value = 2 * pnorm(-abs(estimate / std.error))
        )

      # Set as the final results
      Results_multivariate_analysis <- pooled_results

      # Store models_list for later use in attributes
      current_models_list <- models_list

    } else {
      # Approach for non-mixed models
      res_comb <- vector("list", length(actual_imps))
      for (i in seq_along(actual_imps)) {
        imp_val <- actual_imps[i]
        data_subset <- dplyr::filter(data, !!rlang::sym(imp_col) == imp_val)
        if (nrow(data_subset) == 0) stop(paste("No data for imputation", imp_val))

        model <- switch(model_type,
                        "nb" = MASS::glm.nb(model_formula, data = data_subset),
                        "lm" = glm(model_formula, family = gaussian(), data = data_subset),
                        "bin" = glm(model_formula, family = binomial(), data = data_subset),
                        "poisson" = glm(model_formula, family = poisson(), data = data_subset),
                        "gamma" = glm(model_formula, family = Gamma(), data = data_subset),
                        "quasipoisson" = glm(model_formula, family = quasipoisson(), data = data_subset),
                        "quasibinomial" = glm(model_formula, family = quasibinomial(), data = data_subset),
                        "cox" = coxph(model_formula, data = data_subset),
                        stop("Unsupported model type."))

        res_comb[[i]] <- model
      }

      # Pool results using mice::pool()
      pooled <- mice::pool(res_comb)
      Results_multivariate_analysis <- summary(pooled, conf.int = TRUE, exp = FALSE)
      current_models_list <- res_comb
    }

    # Determine if exponentiation is needed
    exp_required <- model_type %in% c("nb", "bin", "poisson", "quasipoisson", "quasibinomial", "cox")

    # For mixed models with random effects
    if (random_intercept == "Yes") {
      # We already have the pooled_results from above
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(exp_estimate = exp(estimate),
               exp_CI95_lower = exp(`2.5 %`),
               exp_CI95_upper = exp(`97.5 %`)) %>%
        select(term, estimate, `std.error`, `2.5 %`, `97.5 %`, exp_estimate, exp_CI95_lower, exp_CI95_upper, p.value)
    } else {
      # For standard models - Results_multivariate_analysis already assigned above
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(exp_estimate = exp(estimate),
               exp_CI95_lower = exp(`2.5 %`),
               exp_CI95_upper = exp(`97.5 %`)) %>%
        select(term, estimate, `std.error`, `2.5 %`, `97.5 %`, exp_estimate, exp_CI95_lower, exp_CI95_upper, p.value)
    }

    Results_multivariate_analysis$term <- gsub("poly\\(([^,]+), degree = ([0-9]+), raw = TRUE\\)", "poly(\\1, \\2, raw = TRUE)", Results_multivariate_analysis$term)

    # Filter out trial terms if needed
    if (!is.null(trial_col)) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        filter(!grepl(trial_col, term))
    }

    # Filter out random intercept term if needed
    if (random_intercept == "Yes" && !is.null(random_intercept_var)) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        filter(!grepl(paste0("\\(Intercept\\)|SD\\(", random_intercept_var, "\\)"), term))
    }

    # Filter out individual spline terms if requested
    if (!include_spline_terms && length(spline_terms_detected) > 0) {
      # Create patterns to match individual spline components
      spline_patterns <- sapply(spline_terms_detected, function(x) {
        if (grepl("rcs\\(", x)) {
          # Extract variable name from rcs() term
          var_name <- gsub("rcs\\(([^,]+),.*", "\\1", x)
        } else if (grepl("bs\\(", x)) {
          # Extract variable name from bs() term
          var_name <- gsub("bs\\(([^,]+),.*", "\\1", x)
        } else {
          return("")
        }
        var_name <- trimws(var_name)
        # Pattern to match the variable's spline components
        paste0("^", var_name, "'")
      })

      # Combine patterns
      pattern <- paste(spline_patterns, collapse = "|")

      # Keep the main spline terms but filter out the individual components
      if (pattern != "") {
        Results_multivariate_analysis <- Results_multivariate_analysis %>%
          filter(!grepl(pattern, term))
      }
    }

    # Filter out individual polynomial terms if requested
    if (!include_poly_terms && length(poly_terms_detected) > 0) {
      # Create patterns to match individual polynomial components
      poly_patterns <- sapply(poly_terms_detected, function(x) {
        if (grepl("poly\\(", x)) {
          # Extract variable name from poly() term
          var_name <- gsub("poly\\(([^,]+),.*", "\\1", x)
          var_name <- trimws(var_name)
          # Pattern to match the variable's polynomial components
          return(paste0("^poly\\(", var_name, ".*\\)"))
        } else {
          return("")
        }
      })

      # Combine patterns
      pattern <- paste(poly_patterns, collapse = "|")

      # Remove polynomial terms according to pattern
      if (pattern != "") {
        Results_multivariate_analysis <- Results_multivariate_analysis %>%
          filter(!grepl(pattern, term))
      }
    }

    # Add indicator for interaction terms if requested
    if (highlight_interactions && length(interaction_terms) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_interaction = sapply(term, function(t) any(sapply(interaction_terms, function(i) grepl(i, t, fixed = TRUE)))))
    }

    # Add indicator for spline terms if requested
    if (highlight_interactions && length(spline_terms_detected) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_spline = sapply(term, function(t) {
          is_spline_term <- grepl("rcs\\(|bs\\(", t)
          if (is_spline_term) return(TRUE)

          # Check if it's a component of a spline
          for (spline_term in spline_terms_detected) {
            if (grepl("rcs\\(", spline_term)) {
              var_name <- gsub("rcs\\(([^,]+),.*", "\\1", spline_term)
            } else if (grepl("bs\\(", spline_term)) {
              var_name <- gsub("bs\\(([^,]+),.*", "\\1", spline_term)
            } else {
              next
            }
            var_name <- trimws(var_name)
            if (grepl(paste0("^", var_name, "'"), t)) return(TRUE)
          }
          return(FALSE)
        }))
    }

    # Add indicator for polynomial terms if requested
    if (highlight_interactions && length(poly_terms_detected) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_polynomial = sapply(term, function(t) {
          # Check for direct polynomial term
          is_poly_term <- grepl("poly\\(", t)
          if (is_poly_term) return(TRUE)

          # Check for individual polynomial components (like .1, .2, .3)
          for (poly_term in poly_terms_detected) {
            if (grepl("poly\\(", poly_term)) {
              var_name <- gsub("poly\\(([^,]+),.*", "\\1", poly_term)
              var_name <- trimws(var_name)
              degree_match <- regmatches(poly_term, regexpr("degree\\s*=\\s*[0-9]+", poly_term))
              if (length(degree_match) > 0) {
                degree <- as.numeric(gsub("degree\\s*=\\s*", "", degree_match))
              } else {
                degree <- 2  # Default to degree 2 if not found
              }
              # Check for each possible polynomial component
              for (d in 1:degree) {
                if (grepl(paste0("poly\\(", var_name, ".*\\)", d), t)) return(TRUE)
              }
            }
          }
          return(FALSE)
        }))
    }

    # Add random effects attributes if applicable
    if (random_intercept == "Yes") {
      attr(Results_multivariate_analysis, "has_random_effects") <- TRUE
      attr(Results_multivariate_analysis, "random_intercept_var") <- random_intercept_var

      # Store random effects variance components if available
      if (exists("current_models_list") && length(current_models_list) > 0) {
        # Extract from first model as example
        model1 <- current_models_list[[1]]
        if (inherits(model1, "glmmTMB")) {
          re_var <- tryCatch({
            as.data.frame(glmmTMB::VarCorr(model1)$cond)[1, "vcov"]
          }, error = function(e) NA)

          attr(Results_multivariate_analysis, "variance_components") <- re_var
        }

        if (model_type == "cox") {
          # coxme
          if (!is.null(model1$vcoef)) {
            var_comp <- as.numeric(model1$vcoef)
            names(var_comp) <- "Var(Random Intercept)"
            attr(Results_multivariate_analysis, "variance_components") <- var_comp
          }
        } else {
          # lme4 models
          if (inherits(model1, "merMod")) {
            var_comp <- as.data.frame(lme4::VarCorr(model1))
            var_comp <- setNames(var_comp$vcov, paste0("Var(", var_comp$grp, ")"))
            attr(Results_multivariate_analysis, "variance_components") <- var_comp

            # Calculate ICC for linear models if applicable
            if (model_type == "lm") {
              # ICC = tau^2 / (tau^2 + sigma^2)
              tau2 <- var_comp[1]  # Random intercept variance
              sigma2 <- attr(lme4::VarCorr(model1), "sc")^2  # Residual variance
              ICC <- tau2 / (tau2 + sigma2)
              attr(Results_multivariate_analysis, "ICC") <- ICC
            }
          }
        }
      }
    } else {
      attr(Results_multivariate_analysis, "has_random_effects") <- FALSE
    }

    # Add formula and term information as attributes
    attr(Results_multivariate_analysis, "predictor_tested") <- current_predictor
    attr(Results_multivariate_analysis, "model_type_used") <- model_type
    attr(Results_multivariate_analysis, "formula") <- formula_string_current
    attr(Results_multivariate_analysis, "has_interactions") <- length(interaction_terms) > 0
    attr(Results_multivariate_analysis, "interaction_terms") <- if(length(interaction_terms) > 0) interaction_terms else NULL
    attr(Results_multivariate_analysis, "has_splines") <- length(spline_terms_detected) > 0
    attr(Results_multivariate_analysis, "spline_terms") <- if(length(spline_terms_detected) > 0) spline_terms_detected else NULL
    attr(Results_multivariate_analysis, "spline_method") <- if(use_rms) "rcs" else "bs"
    attr(Results_multivariate_analysis, "has_polynomials") <- length(poly_terms_detected) > 0
    attr(Results_multivariate_analysis, "polynomial_terms") <- if(length(poly_terms_detected) > 0) poly_terms_detected else NULL
    attr(Results_multivariate_analysis, "imputations") <- actual_imps
    attr(Results_multivariate_analysis, "n_imp") <- length(actual_imps)

    return(Results_multivariate_analysis)
  }

  # Main execution: fit models for each predictor variable
  results_list <- vector("list", length(predictor_vars))
  names(results_list) <- predictor_vars

  for (i in seq_along(predictor_vars)) {
    current_predictor <- predictor_vars[i]
    message(paste("Fitting model for predictor:", current_predictor))

    results_list[[i]] <- fit_model_for_predictor(current_predictor)
  }

  # Create combined results dataframe for multiple predictors
  if (length(predictor_vars) > 1) {
    create_combined_results <- function() {
      combined_results <- data.frame()

      for (pred_name in names(results_list)) {
        result <- results_list[[pred_name]]

        # Find rows that correspond to the predictor variable (not covariables, intercept, or other terms)
        predictor_rows <- result[grepl(paste0("^", pred_name), result$term) &
                                   !grepl("Intercept", result$term), ]

        # If no exact match found, try partial matching (for categorical variables)
        if (nrow(predictor_rows) == 0) {
          # For categorical variables, the term might be like "GenderMale", "Treatment_stepStep 4", etc.
          predictor_rows <- result[grepl(pred_name, result$term) &
                                     !grepl("Intercept", result$term), ]

          # Further filter out covariables if any are found
          if (!is.null(covariables) && length(covariables) > 0) {
            covariable_pattern <- paste(covariables, collapse = "|")
            predictor_rows <- predictor_rows[!grepl(covariable_pattern, predictor_rows$term), ]
          }
        }

        # Add predictor information to each row and select only desired columns
        if (nrow(predictor_rows) > 0) {
          # Select only the desired columns
          predictor_rows <- predictor_rows[, c("term", "estimate", "std.error", "2.5 %", "97.5 %",
                                               "exp_estimate", "exp_CI95_lower", "exp_CI95_upper", "p.value")]

          combined_results <- rbind(combined_results, predictor_rows)
        }
      }

      # Reset row names
      rownames(combined_results) <- NULL

      return(combined_results)
    }

    # Create the combined results
    combined_predictors_results <- create_combined_results()
  }

  # If only one predictor, return the result directly (maintaining backward compatibility)
  if (length(predictor_vars) == 1) {
    return(results_list[[1]])
  } else {
    # Return list of results for multiple predictors with combined results
    class(results_list) <- c("MI_estimates_multi", "list")
    attr(results_list, "covariables") <- covariables
    attr(results_list, "predictors") <- predictor_vars
    attr(results_list, "model_type") <- model_type
    attr(results_list, "n_predictors") <- length(predictor_vars)
    attr(results_list, "combined_results") <- combined_predictors_results
    return(results_list)
  }
}
