#' Create Spline Plots from Multiple Imputed Data (with Random Effects Support)
#'
#' This function fits restricted cubic spline models to multiply imputed datasets and creates
#' visualization of the relationships between a continuous predictor and an outcome, with optional
#' stratification by subgroups. The function now supports random effects models.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for the model.
#' @param variable_x The continuous predictor variable to be modeled with splines.
#' @param subgroups Optional variable for stratification (default is NULL).
#' @param covariables Optional vector of covariates to adjust for.
#' @param knot_n Number of knots for the restricted cubic spline (default is 4).
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param model_type The model type: "nb" for Negative Binomial, "poisson" for Poisson, "cox" for Cox, or "logistic" for Logistic.
#' @param followup_offset Whether to include offset for follow-up duration ("Yes" or "No").
#' @param followup_col The column representing follow-up duration (required if followup_offset = "Yes").
#' @param trial_factor Whether to include trial as a factor ("Yes" or "No").
#' @param trial_col The column for trial factor adjustment (required if trial_factor = "Yes").
#' @param time_col The time variable for Cox regression (only required if model_type = "cox").
#' @param event_col The event variable for Cox regression (only required if model_type = "cox").
#' @param random_intercept Whether to include a random intercept in the model ("Yes" or "No").
#' @param random_intercept_var The grouping variable for the random intercept (required if random_intercept = "Yes").
#' @param plot_options List of options for plot customization.
#' @param subgroup_as_factor Whether to convert subgroups to factors (default is TRUE).
#' @param subgroup_labels Optional labels for subgroup levels.
#' @param prediction_range Range for predictions as quantiles (default is c(0.01, 0.99)).
#' @param calculate_derivatives Whether to calculate derivatives of the predictions (default is FALSE).
#' @param derivative_points Optional specific points to calculate derivatives at.
#'
#' @return A list containing predictions, model object, plot, and optionally derivatives and thresholds.
#'
#' @importFrom stats as.formula glm binomial poisson gaussian quasipoisson quasibinomial Gamma predict quantile median plogis pchisq setNames
#' @importFrom MASS glm.nb
#' @importFrom survival Surv coxph
#' @importFrom rms rcs
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon xlab ylab labs scale_x_log10 scale_x_continuous scale_y_continuous coord_cartesian scale_color_manual scale_fill_manual scale_linetype_manual guides guide_legend facet_wrap theme_minimal element_rect element_text margin element_blank element_line unit annotate
#' @importFrom dplyr %>%
#' @importFrom lme4 lmer glmer fixef ranef
#' @importFrom glmmTMB glmmTMB
#' @importFrom coxme coxme
#' @export
MI_spline <- function(data,
                      outcome_var,
                      variable_x,
                      subgroups = NULL,
                      covariables = NULL,
                      knot_n = 4,
                      imp_col = ".imp",
                      model_type = "nb",
                      followup_offset = "No",
                      followup_col = NULL,
                      trial_factor = "No",
                      trial_col = NULL,
                      time_col = NULL,
                      event_col = NULL,
                      random_intercept = "No",  # Added parameter
                      random_intercept_var = NULL,  # Added parameter
                      plot_options = NULL,
                      subgroup_as_factor = TRUE,
                      subgroup_labels = NULL,
                      prediction_range = c(0.01, 0.99),
                      calculate_derivatives = FALSE,
                      derivative_points = NULL)  {      # Optional specific points to calculate derivatives at

  # Load required packages
  requireNamespace("rms", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # Load packages for random effects models if needed
  if (random_intercept == "Yes") {
    requireNamespace("lme4", quietly = TRUE)

    if (model_type == "nb") {
      if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        warning("Package 'glmmTMB' not available. Using lme4 with Poisson family as approximation for negative binomial with random effects.")
      }
    }

    if (model_type == "cox") {
      if (!requireNamespace("coxme", quietly = TRUE)) {
        stop("Package 'coxme' is required for Cox models with random effects.")
      }
    }
  }

  # Input validation
  if (!model_type %in% c("nb", "poisson", "cox", "logistic", "lm")) {
    stop("model_type must be one of 'nb', 'poisson', 'cox', 'logistic', or 'lm'")
  }

  if (followup_offset == "Yes" && is.null(followup_col)) {
    stop("followup_col must be provided when followup_offset = 'Yes'")
  }

  if (trial_factor == "Yes" && is.null(trial_col)) {
    stop("trial_col must be provided when trial_factor = 'Yes'")
  }

  if (model_type == "cox" && (is.null(time_col) || is.null(event_col))) {
    stop("time_col and event_col must be provided for Cox models")
  }

  # Validate random intercept parameters
  if (!random_intercept %in% c("Yes", "No")) {
    stop("random_intercept must be either 'Yes' or 'No'")
  }

  if (random_intercept == "Yes" && is.null(random_intercept_var)) {
    stop("If random_intercept = 'Yes', random_intercept_var must be provided")
  }

  if (random_intercept == "Yes" && !random_intercept_var %in% names(data)) {
    stop("random_intercept_var not found in data")
  }

  # Validate prediction_range parameter
  if (length(prediction_range) != 2) {
    stop("prediction_range must be a vector of length 2")
  }

  if (any(prediction_range < 0) || any(prediction_range > 1)) {
    stop("prediction_range values must be between 0 and 1")
  }

  if (prediction_range[1] >= prediction_range[2]) {
    stop("First value of prediction_range must be less than the second value")
  }

  # Select only the first imputation
  Data_Subset <- subset(data, get(imp_col) == 1)

  # Store original subgroup values before conversion, if available
  original_subgroup_values <- NULL
  if (!is.null(subgroups)) {
    original_subgroup_values <- Data_Subset[[subgroups]]
  }

  # Create a mapping between original values and labels if both exist
  subgroup_mapping <- NULL
  if (!is.null(subgroups) && !is.null(subgroup_labels)) {
    # Identify unique values in original data
    unique_values <- sort(unique(Data_Subset[[subgroups]]))

    # Create mapping
    if (length(unique_values) == length(subgroup_labels)) {
      subgroup_mapping <- setNames(subgroup_labels, unique_values)
    }
  }

  # Convert subgroups to factor if requested and needed
  if (!is.null(subgroups) && subgroup_as_factor && !is.factor(Data_Subset[[subgroups]])) {
    # Convert subgroups column to factor
    Data_Subset[[subgroups]] <- as.factor(Data_Subset[[subgroups]])

    # Apply custom labels if provided
    if (!is.null(subgroup_labels)) {
      levels(Data_Subset[[subgroups]]) <- subgroup_labels
    }
  }

  # Setup formula components
  spline_term <- paste0("rcs(", variable_x, ", ", knot_n, ")")

  # Handle subgroups
  if (!is.null(subgroups)) {
    spline_term <- paste0(spline_term, " * ", subgroups)
  }

  # Expand covariables if interaction '*' is used
  expand_terms <- function(term) {
    if (grepl("^rcs\\(", term) || grepl("^bs\\(", term)) {
      # If the term already starts with rcs( or bs(, do not expand it
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

  expanded_covariables <- NULL
  if (!is.null(covariables)) {
    expanded_covariables <- unlist(lapply(covariables, expand_terms))
    expanded_covariables <- unique(expanded_covariables)

    # Helper function to extract variable names from spline terms
    extract_variables_cov <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^bs\\(", term)) {
        inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
        return(trimws(inside))
      } else {
        return(term)
      }
    }

    # Validate that all base variables exist in data
    all_base_vars_cov <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
    all_base_vars_cov <- trimws(all_base_vars_cov)
    all_base_vars_cov <- sapply(all_base_vars_cov, extract_variables_cov)

    missing_vars_cov <- all_base_vars_cov[!all_base_vars_cov %in% names(Data_Subset)]
    if (length(missing_vars_cov) > 0) {
      stop(paste("Covariates not found in data:", paste(missing_vars_cov, collapse = ", ")))
    }

    covariates_str <- paste0(" + ", paste(expanded_covariables, collapse = " + "))
  } else {
    covariates_str <- ""
  }

  # Add trial factor if requested
  trial_str <- ""
  if (trial_factor == "Yes") {
    trial_str <- paste0(" + as.factor(", trial_col, ")")
  }

  # Add offset for count models if requested
  offset_str <- ""
  if (followup_offset == "Yes") {
    offset_str <- paste0(" + offset(log(", followup_col, "))")
  }

  # Add random intercept if requested
  random_effect_str <- ""
  if (random_intercept == "Yes") {
    random_effect_str <- paste0(" + (1 | ", random_intercept_var, ")")
  }

  # Build the complete formula
  if (model_type == "cox") {
    # Cox models need a different formula structure
    formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                          spline_term, covariates_str, trial_str, random_effect_str)
  } else {
    formula_str <- paste0(outcome_var, " ~ ", spline_term, covariates_str,
                          trial_str, offset_str, random_effect_str)
  }

  formula_obj <- as.formula(formula_str)

  # Fit the model based on the specified type and random effects
  model <- NULL
  if (random_intercept == "Yes") {
    # Models with random effects
    if (model_type == "nb") {
      if (requireNamespace("glmmTMB", quietly = TRUE)) {
        model <- glmmTMB::glmmTMB(formula_obj, family = glmmTMB::nbinom2, data = Data_Subset)
      } else {
        warning("Using Poisson mixed model as approximation for negative binomial.")
        model <- lme4::glmer(formula_obj, family = poisson(link = "log"), data = Data_Subset)
      }
    } else if (model_type == "poisson") {
      model <- lme4::glmer(formula_obj, family = poisson(link = "log"), data = Data_Subset)
    } else if (model_type == "logistic") {
      model <- lme4::glmer(formula_obj, family = binomial(link = "logit"), data = Data_Subset)
    } else if (model_type == "lm") {
      model <- lme4::lmer(formula_obj, data = Data_Subset)
    } else if (model_type == "cox") {
      model <- coxme::coxme(formula_obj, data = Data_Subset)
    }
  } else {
    # Standard models without random effects
    if (model_type == "nb") {
      requireNamespace("MASS", quietly = TRUE)
      model <- MASS::glm.nb(formula_obj, data = Data_Subset)
    } else if (model_type == "poisson") {
      model <- glm(formula_obj, family = poisson(link = "log"), data = Data_Subset)
    } else if (model_type == "logistic") {
      model <- glm(formula_obj, family = binomial(link = "logit"), data = Data_Subset)
    } else if (model_type == "lm") {
      model <- glm(formula_obj, family = gaussian(), data = Data_Subset)
    } else if (model_type == "cox") {
      requireNamespace("survival", quietly = TRUE)
      model <- survival::coxph(formula_obj, data = Data_Subset)
    }
  }

  # Create prediction data
  # 1. Create sequence of x values using the prediction_range parameter
  x_values <- seq(
    from = quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
    to = quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
    length.out = 100
  )

  # 2. Create prediction frame
  if (is.null(subgroups)) {
    # Simple prediction frame without subgroups
    pred_data <- data.frame(x_values)
    colnames(pred_data) <- variable_x
  } else {
    # Prediction frame with subgroups
    subgroup_levels <- if (is.factor(Data_Subset[[subgroups]])) {
      levels(Data_Subset[[subgroups]])
    } else {
      sort(unique(Data_Subset[[subgroups]]))
    }

    pred_data <- expand.grid(
      x_values,
      subgroup_levels
    )
    colnames(pred_data) <- c(variable_x, subgroups)

    # Ensure subgroups column is a factor for prediction/plotting
    if (is.factor(Data_Subset[[subgroups]])) {
      pred_data[[subgroups]] <- factor(pred_data[[subgroups]], levels = levels(Data_Subset[[subgroups]]))
    } else if (subgroup_as_factor) {
      pred_data[[subgroups]] <- factor(pred_data[[subgroups]])
      if (!is.null(subgroup_labels)) {
        levels(pred_data[[subgroups]]) <- subgroup_labels
      }
    }
  }

  # 3. Add median/mode values for covariates
  if (!is.null(expanded_covariables)) {
    # Helper to extract the variable name from rcs(...) or bs(...) if needed
    extract_variables_cov <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^bs\\(", term)) {
        inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
        return(trimws(inside))
      } else {
        return(term)
      }
    }

    covariate_vars_for_prediction <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
    covariate_vars_for_prediction <- trimws(covariate_vars_for_prediction)
    covariate_vars_for_prediction <- sapply(covariate_vars_for_prediction, extract_variables_cov)

    for (cov in covariate_vars_for_prediction) {
      if (cov %in% colnames(Data_Subset)) {
        if (is.factor(Data_Subset[[cov]])) {
          pred_data[[cov]] <- as.factor(names(which.max(table(Data_Subset[[cov]]))))
        } else {
          pred_data[[cov]] <- median(Data_Subset[[cov]], na.rm = TRUE)
        }
      }
    }
  }

  # 4. Add trial factor if needed
  if (trial_factor == "Yes") {
    pred_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
  }

  # 5. Add follow-up time for offset if needed
  if (followup_offset == "Yes") {
    pred_data[[followup_col]] <- 365  # 1 year follow-up for standardized predictions
  }

  # 6. Add random effect variable but set to reference level for population-level predictions
  if (random_intercept == "Yes") {
    # For prediction with mixed models, we use the most common level or an "average" random effect
    # This is for population-level (marginal) predictions
    pred_data[[random_intercept_var]] <- names(which.max(table(Data_Subset[[random_intercept_var]])))
  }

  # Get predictions with standard errors
  preds <- NULL

  # Function to handle predictions from mixed models
  get_mixed_model_predictions <- function(model, newdata, type = "link") {
    if (inherits(model, "merMod")) {
      # For lme4 models, use predict with re.form=NA for population-level predictions
      preds <- predict(model, newdata = newdata, re.form = NA, type = type)
      # For standard errors, we need a more complex approach since predict doesn't provide SEs
      # We'll use a bootstrap approach or approximate method based on the fixed effects variance-covariance matrix

      # Get the fixed-effects design matrix for new data
      mm <- model.matrix(formula(model, fixed.only = TRUE), newdata)

      # Get the variance-covariance matrix of fixed effects
      vcov_fixed <- as.matrix(vcov(model))

      # Calculate standard errors for fixed effects predictions
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))

      return(list(fit = preds, se.fit = se_fixed))
    } else if (inherits(model, "glmmTMB")) {
      # For glmmTMB models
      preds <- predict(model, newdata = newdata, re.form = NA, type = type, se.fit = TRUE)
      return(list(fit = preds$fit, se.fit = preds$se.fit))
    } else if (inherits(model, "coxme")) {
      # For coxme models, predict linear predictor
      preds <- predict(model, newdata = newdata, type = "lp")
      # Standard error not directly available, use fixed effects se as approximation
      mm <- model.matrix(formula(model, fixed.only = TRUE), newdata)
      vcov_fixed <- as.matrix(vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))

      return(list(fit = preds, se.fit = se_fixed))
    }
  }

  # Get predictions based on model type
  if (random_intercept == "Yes") {
    # For mixed models
    if (model_type == "cox") {
      # For Cox mixed models
      preds <- get_mixed_model_predictions(model, pred_data)

      # Transform to hazard ratio
      pred_data$prediction <- exp(preds$fit)
      pred_data$lower_ci <- exp(preds$fit - 1.96 * preds$se.fit)
      pred_data$upper_ci <- exp(preds$fit + 1.96 * preds$se.fit)
    } else {
      # For other mixed models
      preds <- get_mixed_model_predictions(model, pred_data)

      # Transform to response scale
      if (model_type %in% c("nb", "poisson")) {
        # For count models: exponentiate for rate/count
        pred_data$prediction <- exp(preds$fit)
        pred_data$lower_ci <- exp(preds$fit - 1.96 * preds$se.fit)
        pred_data$upper_ci <- exp(preds$fit + 1.96 * preds$se.fit)
      } else if (model_type == "logistic") {
        # For logistic: inverse logit for probability
        pred_data$prediction <- plogis(preds$fit)
        pred_data$lower_ci <- plogis(preds$fit - 1.96 * preds$se.fit)
        pred_data$upper_ci <- plogis(preds$fit + 1.96 * preds$se.fit)
      } else if (model_type == "lm") {
        # For linear models: no transformation
        pred_data$prediction <- preds$fit
        pred_data$lower_ci <- preds$fit - 1.96 * preds$se.fit
        pred_data$upper_ci <- preds$fit + 1.96 * preds$se.fit
      }
    }
  } else {
    # Original code for standard models
    if (model_type == "cox") {
      # For Cox models
      preds <- list()
      preds$fit <- predict(model, newdata = pred_data, type = "lp")
      preds$se.fit <- rep(NA, length(preds$fit))  # Cox models don't provide SE directly

      # Transform to hazard ratio
      pred_data$prediction <- exp(preds$fit)
      pred_data$lower_ci <- NA
      pred_data$upper_ci <- NA
    } else {
      # For GLMs
      preds <- predict(model, newdata = pred_data, type = "link", se.fit = TRUE)

      # Transform to response scale
      if (model_type %in% c("nb", "poisson")) {
        pred_data$prediction <- exp(preds$fit)
        pred_data$lower_ci <- exp(preds$fit - 1.96 * preds$se.fit)
        pred_data$upper_ci <- exp(preds$fit + 1.96 * preds$se.fit)
      } else if (model_type == "logistic") {
        pred_data$prediction <- plogis(preds$fit)
        pred_data$lower_ci <- plogis(preds$fit - 1.96 * preds$se.fit)
        pred_data$upper_ci <- plogis(preds$fit + 1.96 * preds$se.fit)
      } else if (model_type == "lm") {
        pred_data$prediction <- preds$fit
        pred_data$lower_ci <- preds$fit - 1.96 * preds$se.fit
        pred_data$upper_ci <- preds$fit + 1.96 * preds$se.fit
      }
    }
  }

  # ============= CALCULATE DERIVATIVES IF REQUESTED =============
  derivative_data <- NULL

  if (calculate_derivatives) {
    # Load rms package for restricted cubic spline operations
    requireNamespace("rms", quietly = TRUE)

    # Define derivative calculation method
    if (is.null(derivative_points)) {
      # If no specific points requested, use a reasonable number of points across the range
      derivative_points <- seq(
        from = quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
        to = quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
        length.out = 20
      )
    }

    # Create data for derivative prediction
    if (is.null(subgroups)) {
      # Simple derivative frame without subgroups
      deriv_data <- data.frame(derivative_points)
      colnames(deriv_data) <- variable_x
    } else {
      # Derivative frame with subgroups
      subgroup_levels <- if (is.factor(Data_Subset[[subgroups]])) {
        levels(Data_Subset[[subgroups]])
      } else {
        sort(unique(Data_Subset[[subgroups]]))
      }

      deriv_data <- expand.grid(
        derivative_points,
        subgroup_levels
      )
      colnames(deriv_data) <- c(variable_x, subgroups)

      # Ensure subgroups column is a factor
      if (is.factor(Data_Subset[[subgroups]])) {
        deriv_data[[subgroups]] <- factor(deriv_data[[subgroups]], levels = levels(Data_Subset[[subgroups]]))
      } else if (subgroup_as_factor) {
        deriv_data[[subgroups]] <- factor(deriv_data[[subgroups]])
        if (!is.null(subgroup_labels)) {
          levels(deriv_data[[subgroups]]) <- subgroup_labels
        }
      }
    }

    # Add the same covariate values used for prediction
    if (!is.null(covariables)) {
      for (cov in covariables) {
        if (cov %in% colnames(pred_data)) {
          deriv_data[[cov]] <- pred_data[[cov]][1]  # Take the first value from pred_data
        }
      }
    }

    # Add trial factor if needed
    if (trial_factor == "Yes") {
      deriv_data[[trial_col]] <- pred_data[[trial_col]][1]
    }

    # Add follow-up time for offset if needed
    if (followup_offset == "Yes") {
      deriv_data[[followup_col]] <- 365  # Same as in prediction
    }

    # Add random intercept variable if needed
    if (random_intercept == "Yes") {
      deriv_data[[random_intercept_var]] <- pred_data[[random_intercept_var]][1]
    }

    # Use finite differences to estimate derivatives
    # We'll compute this for each subgroup separately

    derivative_results <- list()

    # Helper function to calculate derivatives for a specific subgroup
    calc_derivatives_for_subgroup <- function(subgroup_value = NULL) {
      # Filter data for this subgroup if applicable
      if (!is.null(subgroups) && !is.null(subgroup_value)) {
        # Filter prediction data for this subgroup
        subgroup_pred_data <- pred_data[pred_data[[subgroups]] == subgroup_value, ]
        # Sort by x variable
        subgroup_pred_data <- subgroup_pred_data[order(subgroup_pred_data[[variable_x]]), ]

        # Filter derivative points for this subgroup
        subgroup_deriv_data <- deriv_data[deriv_data[[subgroups]] == subgroup_value, ]
      } else {
        # Use all data if no subgroups
        subgroup_pred_data <- pred_data[order(pred_data[[variable_x]]), ]
        subgroup_deriv_data <- deriv_data
      }

      # Calculate derivatives using finite differences on the prediction curve
      derivatives <- numeric(nrow(subgroup_deriv_data))
      lower_ci_derivatives <- numeric(nrow(subgroup_deriv_data))
      upper_ci_derivatives <- numeric(nrow(subgroup_deriv_data))

      for (i in 1:nrow(subgroup_deriv_data)) {
        point <- subgroup_deriv_data[i, variable_x]

        # Find closest points in prediction data
        idx <- which.min(abs(subgroup_pred_data[[variable_x]] - point))

        # Ensure we have points on both sides for differentiation when possible
        if (idx > 1 && idx < nrow(subgroup_pred_data)) {
          # Use central difference
          h1 <- subgroup_pred_data[[variable_x]][idx] - subgroup_pred_data[[variable_x]][idx-1]
          h2 <- subgroup_pred_data[[variable_x]][idx+1] - subgroup_pred_data[[variable_x]][idx]

          # Central difference formula for uneven spacing
          # f'(x) ≈ [h₁²f(x+h₂) - h₂²f(x-h₁) + (h₂²-h₁²)f(x)] / [h₁h₂(h₁+h₂)]
          denominator <- h1 * h2 * (h1 + h2)

          # For the main prediction
          f_minus <- subgroup_pred_data$prediction[idx-1]
          f_center <- subgroup_pred_data$prediction[idx]
          f_plus <- subgroup_pred_data$prediction[idx+1]

          derivatives[i] <- (h1^2 * f_plus - h2^2 * f_minus + (h2^2 - h1^2) * f_center) / denominator

          # For the lower CI
          f_minus_lower <- subgroup_pred_data$lower_ci[idx-1]
          f_center_lower <- subgroup_pred_data$lower_ci[idx]
          f_plus_lower <- subgroup_pred_data$lower_ci[idx+1]

          lower_ci_derivatives[i] <- (h1^2 * f_plus_lower - h2^2 * f_minus_lower +
                                        (h2^2 - h1^2) * f_center_lower) / denominator

          # For the upper CI
          f_minus_upper <- subgroup_pred_data$upper_ci[idx-1]
          f_center_upper <- subgroup_pred_data$upper_ci[idx]
          f_plus_upper <- subgroup_pred_data$upper_ci[idx+1]

          upper_ci_derivatives[i] <- (h1^2 * f_plus_upper - h2^2 * f_minus_upper +
                                        (h2^2 - h1^2) * f_center_upper) / denominator
        } else if (idx == 1) {
          # Use forward difference at the start
          h <- subgroup_pred_data[[variable_x]][idx+1] - subgroup_pred_data[[variable_x]][idx]

          derivatives[i] <- (subgroup_pred_data$prediction[idx+1] - subgroup_pred_data$prediction[idx]) / h
          lower_ci_derivatives[i] <- (subgroup_pred_data$lower_ci[idx+1] - subgroup_pred_data$lower_ci[idx]) / h
          upper_ci_derivatives[i] <- (subgroup_pred_data$upper_ci[idx+1] - subgroup_pred_data$upper_ci[idx]) / h
        } else if (idx == nrow(subgroup_pred_data)) {
          # Use backward difference at the end
          h <- subgroup_pred_data[[variable_x]][idx] - subgroup_pred_data[[variable_x]][idx-1]

          derivatives[i] <- (subgroup_pred_data$prediction[idx] - subgroup_pred_data$prediction[idx-1]) / h
          lower_ci_derivatives[i] <- (subgroup_pred_data$lower_ci[idx] - subgroup_pred_data$lower_ci[idx-1]) / h
          upper_ci_derivatives[i] <- (subgroup_pred_data$upper_ci[idx] - subgroup_pred_data$upper_ci[idx-1]) / h
        }
      }
      # Ensure lower CI is always lower than upper CI
      for (i in 1:length(derivatives)) {
        if (lower_ci_derivatives[i] > upper_ci_derivatives[i]) {
          # Swap the values
          temp <- lower_ci_derivatives[i]
          lower_ci_derivatives[i] <- upper_ci_derivatives[i]
          upper_ci_derivatives[i] <- temp
        }
      }

      # Create result data frame
      result <- data.frame(
        x_point = subgroup_deriv_data[[variable_x]],
        derivative = derivatives,
        lower_ci_derivative = lower_ci_derivatives,
        upper_ci_derivative = upper_ci_derivatives
      )

      # Add subgroup column if applicable
      if (!is.null(subgroups) && !is.null(subgroup_value)) {
        result[[subgroups]] <- subgroup_value
      }

      return(result)
    }

    # Calculate derivatives for each subgroup or for the whole dataset
    if (!is.null(subgroups)) {
      subgroup_levels <- if (is.factor(Data_Subset[[subgroups]])) {
        levels(Data_Subset[[subgroups]])
      } else {
        sort(unique(Data_Subset[[subgroups]]))
      }

      # Calculate for each subgroup
      all_derivatives <- list()
      for (level in subgroup_levels) {
        all_derivatives[[level]] <- calc_derivatives_for_subgroup(level)
      }

      # Combine all results
      derivative_data <- do.call(rbind, all_derivatives)
    } else {
      # Calculate for the whole dataset
      derivative_data <- calc_derivatives_for_subgroup()
    }

    # Check for significant positive derivatives (where lower CI > 0)
    derivative_data$significant_positive <- derivative_data$lower_ci_derivative > 0

    # Find threshold points where the derivative becomes significantly positive
    if (!is.null(subgroups)) {
      # Create a named vector to store thresholds for each subgroup
      threshold_values <- numeric(length(subgroup_levels))
      names(threshold_values) <- subgroup_levels

      for (level in subgroup_levels) {
        subgroup_deriv <- derivative_data[derivative_data[[subgroups]] == level, ]
        subgroup_deriv <- subgroup_deriv[order(subgroup_deriv$x_point), ]

        # Find first point where derivative is significantly positive
        sig_pos_idx <- which(subgroup_deriv$significant_positive)

        if (length(sig_pos_idx) > 0) {
          threshold_values[level] <- subgroup_deriv$x_point[min(sig_pos_idx)]
        } else {
          threshold_values[level] <- NA
        }
      }

      # Now assign the correct threshold to each row based on its subgroup
      derivative_data$threshold <- NA
      for (level in subgroup_levels) {
        derivative_data$threshold[derivative_data[[subgroups]] == level] <- threshold_values[level]
      }
    }
  }
  # ============= ENHANCED PLOTTING CAPABILITIES =============

  # Set default plot options if not provided
  if (is.null(plot_options)) {
    plot_options <- list()
  }

  # Define default plot labels
  x_lab <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x

  # Process x_lab to handle superscripts if requested
  if (!is.null(plot_options$use_superscript) && plot_options$use_superscript) {
    if (is.character(x_lab)) {
      # Instead of changing to expressions, just replace the "^9" with the unicode superscript 9
      # We'll use simple string replacement for maximum compatibility
      x_lab <- gsub("\\^9", "\u2079", x_lab)  # Unicode superscript 9
      x_lab <- gsub("\\^8", "\u2078", x_lab)  # Unicode superscript 8
      x_lab <- gsub("\\^7", "\u2077", x_lab)  # Unicode superscript 7
      x_lab <- gsub("\\^6", "\u2076", x_lab)  # Unicode superscript 6
      x_lab <- gsub("\\^5", "\u2075", x_lab)  # Unicode superscript 5
      x_lab <- gsub("\\^4", "\u2074", x_lab)  # Unicode superscript 4
      x_lab <- gsub("\\^3", "\u00B3", x_lab)  # Unicode superscript 3
      x_lab <- gsub("\\^2", "\u00B2", x_lab)  # Unicode superscript 2
      x_lab <- gsub("\\^1", "\u00B9", x_lab)  # Unicode superscript 1
      x_lab <- gsub("\\^0", "\u2070", x_lab)  # Unicode superscript 0
      # Handle minus sign in superscript
      x_lab <- gsub("\\^-", "\u207B", x_lab)  # Unicode superscript minus
    }
  }

  y_lab <- if (!is.null(plot_options$y_lab)) plot_options$y_lab else {
    if (model_type %in% c("nb", "poisson")) {
      "Predicted Rate"
    } else if (model_type == "logistic") {
      "Predicted Probability"
    } else {
      "Predicted Hazard Ratio"
    }
  }

  # Process y_lab for superscripts too if needed
  if (!is.null(plot_options$use_superscript) && plot_options$use_superscript) {
    if (is.character(y_lab)) {
      # Apply the same superscript replacements to y-axis
      y_lab <- gsub("\\^9", "\u2079", y_lab)  # Unicode superscript 9
      y_lab <- gsub("\\^8", "\u2078", y_lab)  # Unicode superscript 8
      y_lab <- gsub("\\^7", "\u2077", y_lab)  # Unicode superscript 7
      y_lab <- gsub("\\^6", "\u2076", y_lab)  # Unicode superscript 6
      y_lab <- gsub("\\^5", "\u2075", y_lab)  # Unicode superscript 5
      y_lab <- gsub("\\^4", "\u2074", y_lab)  # Unicode superscript 4
      y_lab <- gsub("\\^3", "\u00B3", y_lab)  # Unicode superscript 3
      y_lab <- gsub("\\^2", "\u00B2", y_lab)  # Unicode superscript 2
      y_lab <- gsub("\\^1", "\u00B9", y_lab)  # Unicode superscript 1
      y_lab <- gsub("\\^0", "\u2070", y_lab)  # Unicode superscript 0
      # Handle minus sign in superscript
      y_lab <- gsub("\\^-", "\u207B", y_lab)  # Unicode superscript minus
    }
  }

  # Initialize the plot
  if (!is.null(subgroups)) {
    # Count observations in each subgroup
    if (!is.null(subgroup_mapping) && !is.null(original_subgroup_values)) {
      # If we have custom labels, count from original values
      subgroup_counts <- table(original_subgroup_values)
      # Map the names to the new labels
      names(subgroup_counts) <- subgroup_mapping[names(subgroup_counts)]
    } else {
      # Otherwise count from the current subgroup values in Data_Subset
      subgroup_counts <- table(Data_Subset[[subgroups]])
    }

    # Set up legend labels
    if (!is.null(plot_options$legend_labels)) {
      legend_labels <- plot_options$legend_labels

      # Include counts in the custom legend labels if requested
      if (!is.null(plot_options$include_counts) && plot_options$include_counts) {
        # Get the levels in the correct order from pred_data
        plot_levels <- levels(factor(pred_data[[subgroups]]))

        # If we have a mapping, we need to associate levels with original counts
        if (!is.null(subgroup_mapping)) {
          # Create counts with new labels as names
          legend_labels_with_counts <- character(length(plot_levels))

          for (i in seq_along(plot_levels)) {
            level <- plot_levels[i]
            label <- plot_options$legend_labels[i]

            # Find count for this level
            count <- subgroup_counts[level]
            if (is.na(count)) count <- 0

            # Add count to label
            legend_labels_with_counts[i] <- paste0(label, " (N=", count, ")")
          }

          legend_labels <- legend_labels_with_counts
        } else {
          # Simple case - just add counts to custom labels
          for (i in seq_along(legend_labels)) {
            count <- subgroup_counts[plot_levels[i]]
            if (is.na(count)) count <- 0
            legend_labels[i] <- paste0(legend_labels[i], " (N=", count, ")")
          }
        }
      }
    } else if (!is.null(plot_options$include_counts) && plot_options$include_counts) {
      # Default labels with counts
      legend_labels <- paste0(names(subgroup_counts), " (N=", subgroup_counts, ")")
    } else {
      # Just use the levels as labels
      legend_labels <- levels(factor(pred_data[[subgroups]]))
    }

    # Create base plot with subgroups
    plot <- ggplot2::ggplot(pred_data, ggplot2::aes_string(x = variable_x, y = "prediction",
                                                           color = subgroups,
                                                           fill = subgroups,
                                                           linetype = subgroups))
  } else {
    # Simple plot without subgroups
    plot <- ggplot2::ggplot(pred_data, ggplot2::aes_string(x = variable_x, y = "prediction"))
  }

  # Add lines and ribbons
  line_size <- if (!is.null(plot_options$line_size)) plot_options$line_size else 1
  ribbon_alpha <- if (!is.null(plot_options$ribbon_alpha)) plot_options$ribbon_alpha else 0.3

  plot <- plot +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                         alpha = ribbon_alpha, color = NA)

  # Add labels
  plot <- plot +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab)

  # Add title, subtitle, and caption if provided
  title_args <- list()
  if (!is.null(plot_options$title)) title_args$title <- plot_options$title
  if (!is.null(plot_options$subtitle)) title_args$subtitle <- plot_options$subtitle
  if (!is.null(plot_options$caption)) title_args$caption <- plot_options$caption

  if (length(title_args) > 0) {
    plot <- plot + do.call(ggplot2::labs, title_args)
  }

  # Apply log scale for x-axis if requested
  if (!is.null(plot_options$use_log_x) && plot_options$use_log_x) {
    if (!is.null(plot_options$x_breaks)) {
      plot <- plot + ggplot2::scale_x_log10(breaks = plot_options$x_breaks)
    } else {
      plot <- plot + ggplot2::scale_x_log10()
    }
  } else if (!is.null(plot_options$x_breaks)) {
    plot <- plot + ggplot2::scale_x_continuous(breaks = plot_options$x_breaks)
  }

  # Apply custom y-axis breaks if provided
  if (!is.null(plot_options$y_breaks)) {
    plot <- plot + ggplot2::scale_y_continuous(breaks = plot_options$y_breaks)
  }

  # Apply axis limits using coord_cartesian
  # Start with empty list for coord_cartesian arguments
  coord_args <- list()

  # Add y-axis limits if provided
  if (!is.null(plot_options$y_limits)) {
    coord_args$ylim <- plot_options$y_limits
  }

  # Add x-axis limits if provided
  if (!is.null(plot_options$x_limits)) {
    coord_args$xlim <- plot_options$x_limits
  }

  # Apply coord_cartesian if any limits are set
  if (length(coord_args) > 0) {
    plot <- plot + do.call(ggplot2::coord_cartesian, coord_args)
  }

  # Apply custom colors if provided
  if (!is.null(plot_options$colors) && !is.null(subgroups)) {
    plot <- plot + ggplot2::scale_color_manual(
      values = plot_options$colors,
      labels = legend_labels,
      name = if (!is.null(plot_options$legend_title)) plot_options$legend_title else subgroups
    )

    # Use fill_colors if provided, otherwise use the same colors as lines
    fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors else plot_options$colors

    plot <- plot + ggplot2::scale_fill_manual(
      values = fill_colors,
      guide = "none"
    )
  }

  # Apply custom line types if provided
  if (!is.null(plot_options$line_types) && !is.null(subgroups)) {
    plot <- plot + ggplot2::scale_linetype_manual(
      values = plot_options$line_types,
      labels = legend_labels,
      name = if (!is.null(plot_options$legend_title)) plot_options$legend_title else subgroups
    )
  }

  # Apply custom guides if provided
  if (!is.null(plot_options$custom_guides)) {
    plot <- plot + plot_options$custom_guides
  } else if (!is.null(subgroups) && !is.null(plot_options$colors)) {
    # Default guides setup for subgroups with custom colors
    fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors else plot_options$colors

    plot <- plot + ggplot2::guides(
      linetype = ggplot2::guide_legend(
        override.aes = list(
          color = plot_options$colors,
          fill = fill_colors,
          shape = NA
        ),
        keywidth = 2
      ),
      color = "none",
      fill = "none"
    )
  }

  # Apply faceting if requested
  if (!is.null(plot_options$facet_var)) {
    facet_scales <- if (!is.null(plot_options$facet_scales)) plot_options$facet_scales else "fixed"
    plot <- plot + ggplot2::facet_wrap(as.formula(paste("~", plot_options$facet_var)),
                                       scales = facet_scales)
  }

  # Apply theme customization
  # Initialize theme settings or use provided ones
  theme_settings <- if (!is.null(plot_options$theme_settings)) plot_options$theme_settings else list()

  # Create basic theme
  plot <- plot + ggplot2::theme_minimal(base_size = 10)

  # Define default professional styling
  default_theme <- list(
    legend.background = ggplot2::element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    legend.title = ggplot2::element_text(face = "bold", size = 11),
    legend.text = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(face = "bold", size = 12),
    axis.text.y = ggplot2::element_text(face = "bold", size = 12),
    axis.title.x = ggplot2::element_text(face = "bold", size = 12),
    axis.title.y = ggplot2::element_text(face = "bold", size = 12),
    plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
    plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", size = 0.5),
    axis.ticks = ggplot2::element_line(size = 0.5, colour = "black"),
    axis.ticks.length = ggplot2::unit(0.2, "cm"),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA)
  )

  # Set legend position if provided
  if (!is.null(plot_options$legend_position)) {
    default_theme$legend.position <- plot_options$legend_position
  }

  # Handle individual theme elements from plot_options
  # This allows direct setting of theme elements without the theme_settings structure
  possible_direct_theme_elements <- c(
    "legend.position", "legend.background", "legend.title", "legend.text",
    "axis.text.x", "axis.text.y", "axis.title.x", "axis.title.y", "axis.title",
    "plot.title", "plot.margin", "panel.grid.major", "panel.grid.minor",
    "axis.line", "axis.ticks", "axis.ticks.length", "panel.background", "plot.background"
  )

  for (element_name in possible_direct_theme_elements) {
    # Look for direct element setting in plot_options
    direct_element <- paste0("theme_", gsub("\\.", "_", element_name))

    if (!is.null(plot_options[[direct_element]])) {
      default_theme[[element_name]] <- plot_options[[direct_element]]
    }
  }

  # Merge user theme settings with defaults, user settings taking precedence
  # This still supports the theme_settings structure for backward compatibility
  for (name in names(theme_settings)) {
    default_theme[[name]] <- theme_settings[[name]]
  }

  # Apply the combined theme settings
  plot <- plot + do.call(ggplot2::theme, default_theme)

  # Apply any additional custom scales
  if (!is.null(plot_options$custom_scales)) {
    for (scale in plot_options$custom_scales) {
      plot <- plot + scale
    }
  }


  # Add vertical lines if specified
  if (!is.null(plot_options$vline)) {
    for (v in plot_options$vline) {
      plot <- plot + ggplot2::geom_vline(
        xintercept = v,
        color = if (!is.null(plot_options$vline_color)) plot_options$vline_color else "black",
        linetype = if (!is.null(plot_options$vline_type)) plot_options$vline_type else "dashed",
        linewidth = if (!is.null(plot_options$vline_size)) plot_options$vline_size else 0.5
      )
    }
  }

  # Add horizontal lines if specified
  if (!is.null(plot_options$hline)) {
    for (h in plot_options$hline) {
      plot <- plot + ggplot2::geom_hline(
        yintercept = h,
        color = if (!is.null(plot_options$hline_color)) plot_options$hline_color else "black",
        linetype = if (!is.null(plot_options$hline_type)) plot_options$hline_type else "dashed",
        linewidth = if (!is.null(plot_options$hline_size)) plot_options$hline_size else 0.5
      )
    }
  }

  # Add annotations if specified
  if (!is.null(plot_options$annotations)) {
    for (annotation in plot_options$annotations) {
      plot <- plot + ggplot2::annotate(
        geom = annotation$geom,
        x = annotation$x,
        y = annotation$y,
        label = annotation$label,
        color = if (!is.null(annotation$color)) annotation$color else "black",
        size = if (!is.null(annotation$size)) annotation$size else 4,
        hjust = if (!is.null(annotation$hjust)) annotation$hjust else 0.5,
        vjust = if (!is.null(annotation$vjust)) annotation$vjust else 0.5,
        angle = if (!is.null(annotation$angle)) annotation$angle else 0,
        fontface = if (!is.null(annotation$fontface)) annotation$fontface else "plain"
      )
    }
  }
  # Add shaded areas if specified
  if (!is.null(plot_options$shaded_areas)) {
    for (area in plot_options$shaded_areas) {
      plot <- plot + ggplot2::annotate(
        geom = "rect",
        xmin = area$xmin,
        xmax = area$xmax,
        ymin = if (!is.null(area$ymin)) area$ymin else -Inf,
        ymax = if (!is.null(area$ymax)) area$ymax else Inf,
        fill = if (!is.null(area$fill)) area$fill else "lightblue",
        alpha = if (!is.null(area$alpha)) area$alpha else 0.2
      )
    }
  }



  # Return results
  return(list(
    predictions = pred_data,
    model = model,
    plot = plot,
    plot_data = pred_data,  # Return the plot data for further customization
    prediction_range_values = c(
      quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
      quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE)
    ),   # Return the actual values used for prediction range
    derivatives = derivative_data,  # Add this line to return the derivatives
    threshold = if (calculate_derivatives && !is.null(subgroups)) {
      # Return thresholds as a named vector
      threshold_values
    } else if (calculate_derivatives) {
      # Return single threshold for non-subgroup case
      threshold
    } else {
      NULL  # Return NULL if derivatives weren't calculated
    }
  ))
}
