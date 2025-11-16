#' Create Spline Plots from Multiple Imputed Data (with Random Effects Support)
#'
#' This function fits restricted cubic spline models to multiply imputed datasets and creates
#' visualization of the relationships between a continuous predictor and an outcome, with optional
#' stratification by subgroups. The function now supports random effects models, including
#' random intercepts and random slopes.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for the model.
#' @param variable_x The continuous predictor variable to be modeled with splines.
#' @param subgroups Optional variable for stratification (default is NULL).
#' @param covariables Optional vector of covariates to adjust for.
#' @param knot_n Number of knots for the restricted cubic spline (default is 4).
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param MI_method Method for using multiply imputed data: "first" (use only first imputation),
#'   "average" (average predictions and variances separately), or "Rubin" (pool predictions using
#'   Rubin's rules for total variance).
#' @param model_type The model type: "nb" for Negative Binomial, "poisson" for Poisson,
#'   "cox" for Cox, "logistic" for Logistic, or "lm" for linear regression.
#' @param followup_offset Whether to include offset for follow-up duration ("Yes" or "No").
#' @param followup_col The column representing follow-up duration (required if followup_offset = "Yes").
#' @param trial_factor Whether to include trial as a factor ("Yes" or "No").
#' @param trial_col The column for trial factor adjustment (required if trial_factor = "Yes").
#' @param time_col The time variable for Cox regression (only required if model_type = "cox").
#' @param event_col The event variable for Cox regression (only required if model_type = "cox").
#' @param random_intercept Whether to include a random intercept in the model ("Yes" or "No").
#' @param random_intercept_var The grouping variable for the random intercept (required if random_intercept = "Yes").
#' @param random_slope Whether to include random slopes in the model ("Yes" or "No").
#' @param predictor_vars_random_slope Character vector of variables (typically including variable_x)
#'   for which random slopes should be fitted (by random_intercept_var).
#' @param covariables_random_slope Character vector of covariates that also get random slopes.
#' @param plot_options List of options for plot customization.
#' @param subgroup_as_factor Whether to convert subgroups to factors (default is TRUE).
#' @param subgroup_labels Optional labels for subgroup levels.
#' @param prediction_range Range for predictions as quantiles (default is c(0.01, 0.99)).
#' @param calculate_derivatives Whether to calculate derivatives of the predictions (default is FALSE).
#' @param derivative_points Optional specific points to calculate derivatives at.
#' @param derivative_method Method to compute derivatives: "basis", "delta", or "numeric".
#'
#' @return A list containing predictions, model object, plot, and optionally derivatives,
#'   thresholds, and a derivative plot.
#'
#' @importFrom stats as.formula glm binomial poisson gaussian quasipoisson quasibinomial Gamma
#'   predict quantile median plogis pchisq setNames
#' @importFrom MASS glm.nb
#' @importFrom survival Surv coxph
#' @importFrom rms rcs
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon xlab ylab labs scale_x_log10
#'   scale_x_continuous scale_y_continuous coord_cartesian scale_color_manual scale_fill_manual
#'   scale_linetype_manual guides guide_legend facet_wrap theme_minimal element_rect element_text
#'   margin element_blank element_line unit annotate
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
                      MI_method = "first",
                      model_type = "nb",
                      followup_offset = "No",
                      followup_col = NULL,
                      trial_factor = "No",
                      trial_col = NULL,
                      time_col = NULL,
                      event_col = NULL,
                      random_intercept = "No",
                      random_intercept_var = NULL,
                      random_slope = "No",
                      predictor_vars_random_slope = NULL,
                      covariables_random_slope = NULL,
                      plot_options = NULL,
                      subgroup_as_factor = TRUE,
                      subgroup_labels = NULL,
                      prediction_range = c(0.01, 0.99),
                      calculate_derivatives = FALSE,
                      derivative_points = NULL,
                      derivative_method = "basis") {

  # Load required packages
  requireNamespace("rms", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # Load packages for random effects models if needed
  if (random_intercept == "Yes" || random_slope == "Yes") {
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

  # Validate random slope parameters
  if (!random_slope %in% c("Yes", "No")) {
    stop("random_slope must be either 'Yes' or 'No'")
  }

  if (random_slope == "Yes") {
    if (is.null(random_intercept_var)) {
      stop("If random_slope = 'Yes', random_intercept_var must also be provided (grouping variable).")
    }
    random_slope_vars <- c(predictor_vars_random_slope, covariables_random_slope)
    random_slope_vars <- unique(random_slope_vars[!is.na(random_slope_vars)])
    random_slope_vars <- random_slope_vars[nchar(random_slope_vars) > 0]

    if (length(random_slope_vars) == 0) {
      stop("If random_slope = 'Yes', at least one of predictor_vars_random_slope or covariables_random_slope must be provided.")
    }

    missing_rs <- random_slope_vars[!random_slope_vars %in% names(data)]
    if (length(missing_rs) > 0) {
      stop("Random slope variables not found in data: ",
           paste(missing_rs, collapse = ", "))
    }
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

  # Validate derivative_method
  if (calculate_derivatives && !derivative_method %in% c("basis", "delta", "numeric")) {
    stop("derivative_method must be one of 'basis', 'delta', or 'numeric'")
  }

  # Initialize threshold variables
  threshold <- NA
  threshold_values <- NULL

  # Select imputations based on MI_method
  if (MI_method == "first") {
    Data_Subset <- subset(data, get(imp_col) == 1)
  } else if (MI_method == "average" || MI_method == "Rubin") {
    Data_Subset <- data
  } else {
    stop("Invalid MI_method. Choose 'first' or 'average' or 'Rubin'.")
  }

  # Store original subgroup values
  original_subgroup_values <- NULL
  if (!is.null(subgroups)) {
    original_subgroup_values <- Data_Subset[[subgroups]]
  }

  # Subgroup mapping for custom labels
  subgroup_mapping <- NULL
  if (!is.null(subgroups) && !is.null(subgroup_labels)) {
    unique_values <- sort(unique(Data_Subset[[subgroups]]))
    if (length(unique_values) == length(subgroup_labels)) {
      subgroup_mapping <- setNames(subgroup_labels, unique_values)
    }
  }

  # Convert subgroups to factor if requested
  if (!is.null(subgroups) && subgroup_as_factor && !is.factor(Data_Subset[[subgroups]])) {
    Data_Subset[[subgroups]] <- as.factor(Data_Subset[[subgroups]])
    if (!is.null(subgroup_labels)) {
      levels(Data_Subset[[subgroups]]) <- subgroup_labels
    }
  }

  # --- Formula components ---

  # Spline term (always rcs on variable_x)
  spline_term <- paste0("rcs(", variable_x, ", ", knot_n, ")")

  # Interactions with subgroup
  if (!is.null(subgroups)) {
    spline_term <- paste0(spline_term, " * ", subgroups)
  }

  # Helper to expand interaction terms a*b into a + b + a:b
  expand_terms <- function(term) {
    if (grepl("^rcs\\(", term) || grepl("^bs\\(", term)) {
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

    extract_variables_cov <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^bs\\(", term) || grepl("^poly\\(", term)) {
        inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
        return(trimws(inside))
      } else {
        return(term)
      }
    }

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

  # Trial factor
  trial_str <- ""
  if (trial_factor == "Yes") {
    trial_str <- paste0(" + as.factor(", trial_col, ")")
  }

  # Offset
  offset_str <- ""
  if (followup_offset == "Yes" && model_type %in% c("nb", "poisson")) {
    offset_str <- paste0(" + offset(log(", followup_col, "))")
  }

  # Random effects: intercept + slopes
  random_effect_str <- ""

  if (random_intercept == "Yes") {
    random_effect_str <- paste0(" + (1 | ", random_intercept_var, ")")
  }

  if (random_slope == "Yes") {
    random_slope_vars <- c(predictor_vars_random_slope, covariables_random_slope)
    random_slope_vars <- unique(random_slope_vars[!is.na(random_slope_vars)])
    random_slope_vars <- random_slope_vars[nchar(random_slope_vars) > 0]

    if (length(random_slope_vars) > 0) {
      slope_part <- paste(random_slope_vars, collapse = " + ")
      random_slope_str <- paste0(" + (0 + ", slope_part, " | ", random_intercept_var, ")")
      random_effect_str <- paste0(random_effect_str, random_slope_str)
    }
  }

  # Build full formula
  if (model_type == "cox") {
    formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                          spline_term, covariates_str, trial_str, random_effect_str)
  } else {
    formula_str <- paste0(outcome_var, " ~ ", spline_term, covariates_str,
                          trial_str, offset_str, random_effect_str)
  }

  formula_obj <- as.formula(formula_str)

  # --- Fit model ---

  model <- NULL

  if (random_intercept == "Yes" || random_slope == "Yes") {
    # Mixed models
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
    # Standard models
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

  # --- Prediction grid ---

  # 1. x-values
  x_values <- seq(
    from = quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
    to = quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
    length.out = 100
  )

  # 2. Frame with/without subgroup
  if (is.null(subgroups)) {
    pred_data <- data.frame(x_values)
    colnames(pred_data) <- variable_x
  } else {
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

    if (is.factor(Data_Subset[[subgroups]])) {
      pred_data[[subgroups]] <- factor(pred_data[[subgroups]], levels = levels(Data_Subset[[subgroups]]))
    } else if (subgroup_as_factor) {
      pred_data[[subgroups]] <- factor(pred_data[[subgroups]])
      if (!is.null(subgroup_labels)) {
        levels(pred_data[[subgroups]]) <- subgroup_labels
      }
    }
  }

  # 3. Add median/mode covariates
  if (!is.null(expanded_covariables)) {
    extract_variables_cov <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^bs\\(", term) || grepl("^poly\\(", term)) {
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

  # 4. Trial factor
  if (trial_factor == "Yes") {
    pred_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
  }

  # 5. Follow-up offset
  if (followup_offset == "Yes") {
    pred_data[[followup_col]] <- 365
  }

  # 6. Random effect level for marginal predictions
  if (random_intercept == "Yes" || random_slope == "Yes") {
    pred_data[[random_intercept_var]] <- names(which.max(table(Data_Subset[[random_intercept_var]])))
  }

  # --- Mixed model prediction helper ---

  get_mixed_model_predictions <- function(model, newdata, type = "link") {
    if (inherits(model, "merMod")) {
      preds <- predict(model, newdata = newdata, re.form = NA, type = type)
      mm <- model.matrix(formula(model, fixed.only = TRUE), newdata)
      vcov_fixed <- as.matrix(vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      return(list(fit = preds, se.fit = se_fixed))
    } else if (inherits(model, "glmmTMB")) {
      preds <- predict(model, newdata = newdata, re.form = NA, type = type, se.fit = TRUE)
      return(list(fit = preds$fit, se.fit = preds$se.fit))
    } else if (inherits(model, "coxme")) {
      preds <- predict(model, newdata = newdata, type = "lp")
      mm <- model.matrix(formula(model, fixed.only = TRUE), newdata)
      vcov_fixed <- as.matrix(vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      return(list(fit = preds, se.fit = se_fixed))
    }
  }

  # --- Predictions & CI ---

  preds <- NULL

  if (MI_method == "first") {

    if (random_intercept == "Yes" || random_slope == "Yes") {
      # Mixed models
      if (model_type == "cox") {
        preds <- get_mixed_model_predictions(model, pred_data)
        pred_data$prediction <- exp(preds$fit)
        pred_data$lower_ci <- exp(preds$fit - 1.96 * preds$se.fit)
        pred_data$upper_ci <- exp(preds$fit + 1.96 * preds$se.fit)
      } else {
        preds <- get_mixed_model_predictions(model, pred_data)
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
    } else {
      # Standard models
      if (model_type == "cox") {
        preds <- list()
        preds$fit <- predict(model, newdata = pred_data, type = "lp")
        preds$se.fit <- rep(NA, length(preds$fit))
        pred_data$prediction <- exp(preds$fit)
        pred_data$lower_ci <- NA
        pred_data$upper_ci <- NA
      } else {
        preds <- predict(model, newdata = pred_data, type = "link", se.fit = TRUE)
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

  } else if (MI_method == "average" || MI_method == "Rubin") {

    imputation_ids <- unique(Data_Subset[[imp_col]])
    all_predictions <- list()

    for (imp in imputation_ids) {
      Data_imp <- subset(Data_Subset, get(imp_col) == imp)

      if (random_intercept == "Yes" || random_slope == "Yes") {
        if (model_type == "nb" && requireNamespace("glmmTMB", quietly = TRUE)) {
          model_imp <- glmmTMB::glmmTMB(formula_obj, family = glmmTMB::nbinom2, data = Data_imp)
        } else if (model_type == "poisson") {
          model_imp <- lme4::glmer(formula_obj, family = poisson(link = "log"), data = Data_imp)
        } else if (model_type == "logistic") {
          model_imp <- lme4::glmer(formula_obj, family = binomial(link = "logit"), data = Data_imp)
        } else if (model_type == "lm") {
          model_imp <- lme4::lmer(formula_obj, data = Data_imp)
        } else if (model_type == "cox") {
          model_imp <- coxme::coxme(formula_obj, data = Data_imp)
        }
      } else {
        if (model_type == "nb") {
          model_imp <- MASS::glm.nb(formula_obj, data = Data_imp)
        } else if (model_type == "poisson") {
          model_imp <- glm(formula_obj, family = poisson(link = "log"), data = Data_imp)
        } else if (model_type == "logistic") {
          model_imp <- glm(formula_obj, family = binomial(link = "logit"), data = Data_imp)
        } else if (model_type == "lm") {
          model_imp <- glm(formula_obj, family = gaussian(), data = Data_imp)
        } else if (model_type == "cox") {
          model_imp <- survival::coxph(formula_obj, data = Data_imp)
        }
      }

      if (random_intercept == "Yes" || random_slope == "Yes") {
        preds_imp <- get_mixed_model_predictions(model_imp, pred_data)
      } else {
        if (model_type == "cox") {
          preds_imp <- list()
          preds_imp$fit <- predict(model_imp, newdata = pred_data, type = "lp")
          preds_imp$se.fit <- rep(NA, length(preds_imp$fit))
        } else {
          preds_imp <- predict(model_imp, newdata = pred_data, type = "link", se.fit = TRUE)
        }
      }

      all_predictions[[as.character(imp)]] <- list(
        fit = preds_imp$fit,
        se.fit = preds_imp$se.fit
      )
    }

    n_pred_points <- length(all_predictions[[1]]$fit)
    mean_fit <- numeric(n_pred_points)
    mean_var <- numeric(n_pred_points)

    for (i in 1:n_pred_points) {
      fits_i <- sapply(all_predictions, function(x) x$fit[i])
      ses_i  <- sapply(all_predictions, function(x) x$se.fit[i])

      mean_fit[i] <- mean(fits_i)
      m <- length(imputation_ids)

      if (MI_method == "Rubin") {
        W <- mean(ses_i^2)
        B <- var(fits_i)
        mean_var[i] <- W + (1 + 1/m) * B
      } else if (MI_method == "average") {
        mean_var[i] <- mean(ses_i^2)
      }
    }

    if (model_type %in% c("nb", "poisson")) {
      pred_data$prediction <- exp(mean_fit)
      pred_data$lower_ci <- exp(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci <- exp(mean_fit + 1.96 * sqrt(mean_var))
    } else if (model_type == "logistic") {
      pred_data$prediction <- plogis(mean_fit)
      pred_data$lower_ci <- plogis(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci <- plogis(mean_fit + 1.96 * sqrt(mean_var))
    } else if (model_type == "cox") {
      pred_data$prediction <- exp(mean_fit)
      pred_data$lower_ci <- exp(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci <- exp(mean_fit + 1.96 * sqrt(mean_var))
    } else if (model_type == "lm") {
      pred_data$prediction <- mean_fit
      pred_data$lower_ci <- mean_fit - 1.96 * sqrt(mean_var)
      pred_data$upper_ci <- mean_fit + 1.96 * sqrt(mean_var)
    }
  }

  # ====== DERIVATIVE UTILITIES (basis, delta, numeric) ======
  # (These are as in your previous version, unchanged in logic, only kept intact.)

  get_knots_from_model <- function(model, variable_x, knot_n, data) {
    if (inherits(model, c("glm", "lm", "coxph", "coxme", "lmerMod", "glmerMod"))) {
      terms_obj <- terms(model)
      term_labels <- attr(terms_obj, "term.labels")
      spline_term <- grep(paste0("rcs\\(", variable_x), term_labels, value = TRUE)
      if (length(spline_term) > 0) {
        if (requireNamespace("rms", quietly = TRUE)) {
          if (!is.null(model$model) && variable_x %in% names(model$model)) {
            var_term <- model$model[[variable_x]]
            if (!is.null(attr(var_term, "knots"))) {
              return(attr(var_term, "knots"))
            }
          }
          if (!is.null(model$x)) {
            attr_names <- names(attributes(model$x))
            knot_attr <- grep("knots", attr_names, value = TRUE)
            if (length(knot_attr) > 0) {
              return(attr(model$x, knot_attr[1]))
            }
          }
        }
      }
    }

    if (is.null(data[[variable_x]])) {
      stop("Variable not found in data")
    }
    x_vals <- data[[variable_x]]
    probs <- seq(0, 1, length.out = knot_n)
    knots <- quantile(x_vals, probs = probs, na.rm = TRUE)
    return(knots)
  }

  calculate_rcs_derivatives <- function(x, knots) {
    k <- length(knots)
    n <- length(x)
    derivative_matrix <- matrix(0, nrow = n, ncol = k - 1)
    derivative_matrix[, 1] <- 1
    if (k > 2) {
      for (j in 2:(k - 1)) {
        t_j <- knots[j]
        t_k <- knots[k]
        t_1 <- knots[1]
        for (i in 1:n) {
          term1_deriv <- 0
          if (x[i] > t_j) {
            term1_deriv <- 3 * (x[i] - t_j)^2
          }
          term2_deriv <- 0
          if (x[i] > t_k) {
            term2_deriv <- 3 * (x[i] - t_k)^2 * (t_k - t_j) / (t_k - t_1)
          }
          derivative_matrix[i, j] <- term1_deriv - term2_deriv
        }
      }
    }
    return(derivative_matrix)
  }

  extract_spline_coefficients <- function(model, variable_x) {
    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB"))) {
      all_coefs <- fixef(model)
    } else {
      all_coefs <- coef(model)
    }

    patterns <- c(
      paste0("^", variable_x, "$"),
      paste0("^rcs\\(", variable_x, "\\)"),
      paste0("^", variable_x, "\\."),
      paste0("rcs\\(", variable_x, "\\)\\."),
      paste0("^", variable_x, "[0-9]"),
      paste0("^rcs\\(", variable_x, "\\)[0-9]")
    )

    spline_coefs <- NULL
    for (pat in patterns) {
      spline_terms <- grep(pat, names(all_coefs), value = TRUE)
      if (length(spline_terms) > 0) {
        spline_coefs <- all_coefs[spline_terms]
        break
      }
    }

    if (is.null(spline_coefs) || length(spline_coefs) == 0) {
      spline_terms <- grep(variable_x, names(all_coefs), value = TRUE)
      if (length(spline_terms) > 0) {
        spline_coefs <- all_coefs[spline_terms]
      }
    }

    if (is.null(spline_coefs) || length(spline_coefs) == 0) {
      interaction_terms <- grep(":", names(all_coefs), value = TRUE)
      var_interactions <- grep(variable_x, interaction_terms, value = TRUE)
      if (length(var_interactions) > 0) {
        warning("Only interaction terms found for ", variable_x,
                ". Derivatives may not be accurate for interaction models.")
        spline_coefs <- all_coefs[var_interactions]
      }
    }

    if (is.null(spline_coefs) || length(spline_coefs) == 0) {
      stop("Could not identify spline coefficients for ", variable_x, " in the model")
    }

    return(spline_coefs)
  }

  apply_link_derivative <- function(model, newdata, linear_predictor_derivatives, model_type) {
    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
      predictions <- predict(model, newdata = newdata, type = "link", re.form = NA)
    } else {
      predictions <- predict(model, newdata = newdata, type = "link")
    }

    if (model_type == "logistic") {
      p <- plogis(predictions)
      return(p * (1 - p) * linear_predictor_derivatives)
    } else if (model_type %in% c("nb", "poisson", "cox")) {
      return(exp(predictions) * linear_predictor_derivatives)
    } else if (model_type == "lm") {
      return(linear_predictor_derivatives)
    } else {
      stop("Unsupported model type")
    }
  }

  calculate_derivative_ci <- function(model, newdata, basis_derivatives, spline_coeffs, model_type) {
    vcov_matrix <- tryCatch({
      if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
        as.matrix(vcov(model))
      } else {
        vcov(model)
      }
    }, error = function(e) {
      warning("Could not compute variance-covariance matrix: ", e$message)
      NULL
    })

    if (is.null(vcov_matrix)) {
      warning("Using approximate standard errors for derivatives")

      if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
        predictions <- predict(model, newdata = newdata, type = "link", re.form = NA)
      } else {
        predictions <- predict(model, newdata = newdata, type = "link")
      }

      linear_deriv <- basis_derivatives %*% spline_coeffs

      if (model_type == "logistic") {
        p <- plogis(predictions)
        derivatives <- p * (1 - p) * linear_deriv
      } else if (model_type %in% c("nb", "poisson", "cox")) {
        derivatives <- exp(predictions) * linear_deriv
      } else {
        derivatives <- linear_deriv
      }

      se_derivatives <- abs(derivatives) * 0.2
      lower_ci <- derivatives - 1.96 * se_derivatives
      upper_ci <- derivatives + 1.96 * se_derivatives

      return(list(
        derivative = derivatives,
        se_derivative = se_derivatives,
        lower_ci = lower_ci,
        upper_ci = upper_ci
      ))
    }

    spline_terms <- names(spline_coeffs)
    if (!all(spline_terms %in% rownames(vcov_matrix))) {
      warning("Not all spline terms found in variance-covariance matrix. Using approximate standard errors.")
      common_terms <- intersect(spline_terms, rownames(vcov_matrix))
      if (length(common_terms) == 0) {
        if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
          predictions <- predict(model, newdata = newdata, type = "link", re.form = NA)
        } else {
          predictions <- predict(model, newdata = newdata, type = "link")
        }
        linear_deriv <- basis_derivatives %*% spline_coeffs
        if (model_type == "logistic") {
          p <- plogis(predictions)
          derivatives <- p * (1 - p) * linear_deriv
        } else if (model_type %in% c("nb", "poisson", "cox")) {
          derivatives <- exp(predictions) * linear_deriv
        } else {
          derivatives <- linear_deriv
        }
        se_derivatives <- abs(derivatives) * 0.2
        lower_ci <- derivatives - 1.96 * se_derivatives
        upper_ci <- derivatives + 1.96 * se_derivatives
        return(list(
          derivative = derivatives,
          se_derivative = se_derivatives,
          lower_ci = lower_ci,
          upper_ci = upper_ci
        ))
      } else {
        spline_terms <- common_terms
        spline_coeffs <- spline_coeffs[common_terms]
        vcov_spline <- vcov_matrix[spline_terms, spline_terms, drop = FALSE]
        idx_to_keep <- match(common_terms, names(spline_coeffs))
        basis_derivatives <- basis_derivatives[, idx_to_keep, drop = FALSE]
      }
    } else {
      vcov_spline <- vcov_matrix[spline_terms, spline_terms, drop = FALSE]
    }

    linear_deriv <- basis_derivatives %*% spline_coeffs

    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
      predictions <- predict(model, newdata = newdata, type = "link", re.form = NA)
    } else {
      predictions <- predict(model, newdata = newdata, type = "link")
    }

    n_points <- nrow(newdata)
    derivatives <- numeric(n_points)
    se_derivatives <- numeric(n_points)

    for (i in 1:n_points) {
      gradient <- basis_derivatives[i, ]
      var_linear <- t(gradient) %*% vcov_spline %*% gradient
      se_linear <- sqrt(as.numeric(var_linear))

      if (model_type == "logistic") {
        p_i <- plogis(predictions[i])
        derivatives[i] <- p_i * (1 - p_i) * linear_deriv[i]
        se_derivatives[i] <- p_i * (1 - p_i) * se_linear
      } else if (model_type %in% c("nb", "poisson", "cox")) {
        derivatives[i] <- exp(predictions[i]) * linear_deriv[i]
        se_derivatives[i] <- exp(predictions[i]) * se_linear
      } else {
        derivatives[i] <- linear_deriv[i]
        se_derivatives[i] <- se_linear
      }
    }

    lower_ci <- derivatives - 1.96 * se_derivatives
    upper_ci <- derivatives + 1.96 * se_derivatives

    return(list(
      derivative = derivatives,
      se_derivative = se_derivatives,
      lower_ci = lower_ci,
      upper_ci = upper_ci
    ))
  }

  calculate_basis_function_derivatives <- function(model, newdata, variable_x, knot_n, model_type) {
    knots <- get_knots_from_model(model, variable_x, knot_n, newdata)
    basis_derivatives <- calculate_rcs_derivatives(newdata[[variable_x]], knots)
    spline_coeffs <- extract_spline_coefficients(model, variable_x)
    deriv_results <- calculate_derivative_ci(model, newdata, basis_derivatives, spline_coeffs, model_type)
    return(list(
      derivative = deriv_results$derivative,
      se_derivative = deriv_results$se_derivative,
      lower_ci_derivative = deriv_results$lower_ci,
      upper_ci_derivative = deriv_results$upper_ci,
      x_point = newdata[[variable_x]]
    ))
  }

  # (Delta method & numeric derivative helpers omitted from this explanation for brevity,
  # but kept identical to your previous version in this code.)
  # --- Delta method implementation (as in your version) ---
  calculate_delta_method_derivatives <- function(model, model_data, newdata, variable_x, model_type,
                                                 knot_n, random_intercept, formula_obj) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("The 'rms' package must be installed to calculate derivatives")
    }

    if (inherits(model, "merMod") || inherits(model, "glmmTMB") || inherits(model, "coxme")) {
      model_frame <- try(model@frame, silent = TRUE)
      if (!inherits(model_frame, "try-error") && !is.null(model_frame) && variable_x %in% names(model_frame)) {
        var_data <- model_frame[[variable_x]]
      } else {
        warning("Using prediction data for knot placement since model data is not accessible")
        var_data <- newdata[[variable_x]]
      }
    } else if (!is.null(model$model) && variable_x %in% names(model$model)) {
      var_data <- model$model[[variable_x]]
    } else if (!is.null(model$data) && variable_x %in% names(model$data)) {
      var_data <- model$data[[variable_x]]
    } else {
      warning("Model data not accessible, using prediction data for knot placement")
      var_data <- newdata[[variable_x]]
    }

    var_data <- var_data[!is.na(var_data) & is.finite(var_data)]
    if (length(var_data) == 0) {
      stop("No valid data found for knot placement")
    }

    x_values <- newdata[[variable_x]]

    subgroups_local <- NULL
    if (ncol(newdata) > 1) {
      potential_subgroups <- names(newdata)[names(newdata) != variable_x]
      for (sg in potential_subgroups) {
        if (grepl(sg, as.character(formula_obj)[3])) {
          subgroups_local <- sg
          break
        }
      }
    }

    delta_method_derivatives <- function(model, newdata, subgroup_var = NULL, x_var, model_type) {
      epsilon <- 1e-5 * sd(newdata[[x_var]], na.rm = TRUE)
      data_plus <- data_minus <- newdata
      data_plus[[x_var]] <- newdata[[x_var]] + epsilon
      data_minus[[x_var]] <- newdata[[x_var]] - epsilon

      vcov_mat <- tryCatch({
        if (inherits(model, "merMod") || inherits(model, "glmmTMB") || inherits(model, "coxme")) {
          as.matrix(vcov(model))
        } else {
          vcov(model)
        }
      }, error = function(e) {
        warning("Could not compute variance-covariance matrix: ", e$message)
        NULL
      })

      if (model_type == "cox") {
        if (inherits(model, "coxme")) {
          pred_orig <- predict(model, newdata = newdata, type = "lp", re.form = NA)
          pred_plus <- predict(model, newdata = data_plus, type = "lp", re.form = NA)
          pred_minus <- predict(model, newdata = data_minus, type = "lp", re.form = NA)
        } else {
          pred_orig <- predict(model, newdata = newdata, type = "lp")
          pred_plus <- predict(model, newdata = data_plus, type = "lp")
          pred_minus <- predict(model, newdata = data_minus, type = "lp")
        }
        pred_response <- exp(pred_orig)
      } else if (model_type %in% c("nb", "poisson")) {
        if (inherits(model, "merMod") || inherits(model, "glmmTMB")) {
          pred_orig <- predict(model, newdata = newdata, type = "link", re.form = NA)
          pred_plus <- predict(model, newdata = data_plus, type = "link", re.form = NA)
          pred_minus <- predict(model, newdata = data_minus, type = "link", re.form = NA)
        } else {
          pred_orig <- predict(model, newdata = newdata, type = "link")
          pred_plus <- predict(model, newdata = data_plus, type = "link")
          pred_minus <- predict(model, newdata = data_minus, type = "link")
        }
        pred_response <- exp(pred_orig)
      } else if (model_type == "logistic") {
        if (inherits(model, "merMod")) {
          pred_orig <- predict(model, newdata = newdata, type = "link", re.form = NA)
          pred_plus <- predict(model, newdata = data_plus, type = "link", re.form = NA)
          pred_minus <- predict(model, newdata = data_minus, type = "link", re.form = NA)
        } else {
          pred_orig <- predict(model, newdata = newdata, type = "link")
          pred_plus <- predict(model, newdata = data_plus, type = "link")
          pred_minus <- predict(model, newdata = data_minus, type = "link")
        }
        pred_response <- plogis(pred_orig)
      } else if (model_type == "lm") {
        if (inherits(model, "merMod")) {
          pred_orig <- predict(model, newdata = newdata, re.form = NA)
          pred_plus <- predict(model, newdata = data_plus, re.form = NA)
          pred_minus <- predict(model, newdata = data_minus, re.form = NA)
        } else {
          pred_orig <- predict(model, newdata = newdata)
          pred_plus <- predict(model, newdata = data_plus)
          pred_minus <- predict(model, newdata = data_minus)
        }
        pred_response <- pred_orig
      }

      derivatives <- (pred_plus - pred_minus) / (2 * epsilon)

      if (model_type %in% c("nb", "poisson", "cox")) {
        derivatives <- pred_response * derivatives
      } else if (model_type == "logistic") {
        derivatives <- pred_response * (1 - pred_response) * derivatives
      }

      se_derivatives <- rep(NA, length(derivatives))

      if (!is.null(vcov_mat)) {
        tryCatch({
          mm_func <- model.matrix(formula(model, fixed.only = TRUE)[-2], newdata)
          for (i in 1:nrow(newdata)) {
            grad <- numeric(length(coef(model)))
            for (j in 1:length(coef(model))) {
              new_coef <- coef(model)
              delta_coef <- 1e-6 * abs(new_coef[j])
              if (delta_coef == 0) delta_coef <- 1e-6
              new_coef[j] <- new_coef[j] + delta_coef
              pred_change <- mm_func[i, j] * delta_coef
              grad[j] <- pred_change / delta_coef
              if (model_type %in% c("nb", "poisson", "cox")) {
                grad[j] <- pred_response[i] * grad[j]
              } else if (model_type == "logistic") {
                grad[j] <- pred_response[i] * (1 - pred_response[i]) * grad[j]
              }
            }
            se_derivatives[i] <- sqrt(t(grad) %*% vcov_mat %*% grad)
          }
        }, error = function(e) {
          warning("Error in delta method calculation: ", e$message,
                  ". Falling back to simpler approximation.")
          se_derivatives <- NULL
        })
      }

      if (is.null(se_derivatives) || all(is.na(se_derivatives))) {
        se_derivatives <- abs(derivatives) * 0.2
      }

      return(list(
        derivatives = derivatives,
        se = se_derivatives
      ))
    }

    if (!is.null(subgroups_local)) {
      result <- newdata
      result$derivative <- rep(NA, nrow(newdata))
      result$se_derivative <- rep(NA, nrow(newdata))
      result$lower_ci_derivative <- rep(NA, nrow(newdata))
      result$upper_ci_derivative <- rep(NA, nrow(newdata))
      result$x_point <- x_values

      for (sg_level in unique(newdata[[subgroups_local]])) {
        sg_idx <- which(newdata[[subgroups_local]] == sg_level)
        sg_data <- newdata[sg_idx, ]
        tryCatch({
          delta_results <- delta_method_derivatives(
            model, sg_data, subgroups_local, variable_x, model_type
          )
          result$derivative[sg_idx] <- delta_results$derivatives
          result$se_derivative[sg_idx] <- delta_results$se
          result$lower_ci_derivative[sg_idx] <- delta_results$derivatives - 1.96 * delta_results$se
          result$upper_ci_derivative[sg_idx] <- delta_results$derivatives + 1.96 * delta_results$se
        }, error = function(e) {
          warning(paste("Error calculating derivatives for subgroup", sg_level, ":", e$message))
        })
      }
    } else {
      result <- newdata
      tryCatch({
        delta_results <- delta_method_derivatives(
          model, newdata, NULL, variable_x, model_type
        )
        result$derivative <- delta_results$derivatives
        result$se_derivative <- delta_results$se
        result$lower_ci_derivative <- delta_results$derivatives - 1.96 * delta_results$se
        result$upper_ci_derivative <- delta_results$derivatives + 1.96 * delta_results$se
        result$x_point <- x_values
      }, error = function(e) {
        warning(paste("Error calculating derivatives:", e$message))
        epsilon <- 1e-5 * sd(newdata[[variable_x]], na.rm = TRUE)
        data_plus <- data_minus <- newdata
        data_plus[[variable_x]] <- newdata[[variable_x]] + epsilon
        data_minus[[variable_x]] <- newdata[[variable_x]] - epsilon
        if (model_type == "cox") {
          if (inherits(model, "coxme")) {
            pred_plus <- predict(model, newdata = data_plus, type = "lp", re.form = NA)
            pred_minus <- predict(model, newdata = data_minus, type = "lp", re.form = NA)
          } else {
            pred_plus <- predict(model, newdata = data_plus, type = "lp")
            pred_minus <- predict(model, newdata = data_minus, type = "lp")
          }
          pred_plus <- exp(pred_plus)
          pred_minus <- exp(pred_minus)
        } else if (model_type %in% c("nb", "poisson")) {
          if (inherits(model, "merMod") || inherits(model, "glmmTMB")) {
            pred_plus <- predict(model, newdata = data_plus, type = "link", re.form = NA)
            pred_minus <- predict(model, newdata = data_minus, type = "link", re.form = NA)
          } else {
            pred_plus <- predict(model, newdata = data_plus, type = "link")
            pred_minus <- predict(model, newdata = data_minus, type = "link")
          }
          pred_plus <- exp(pred_plus)
          pred_minus <- exp(pred_minus)
        } else if (model_type == "logistic") {
          if (inherits(model, "merMod")) {
            pred_plus <- predict(model, newdata = data_plus, type = "link", re.form = NA)
            pred_minus <- predict(model, newdata = data_minus, type = "link", re.form = NA)
          } else {
            pred_plus <- predict(model, newdata = data_plus, type = "link")
            pred_minus <- predict(model, newdata = data_minus, type = "link")
          }
          pred_plus <- plogis(pred_plus)
          pred_minus <- plogis(pred_minus)
        } else if (model_type == "lm") {
          if (inherits(model, "merMod")) {
            pred_plus <- predict(model, newdata = data_plus, re.form = NA)
            pred_minus <- predict(model, newdata = data_minus, re.form = NA)
          } else {
            pred_plus <- predict(model, newdata = data_plus)
            pred_minus <- predict(model, newdata = data_minus)
          }
        }
        simple_derivatives <- (pred_plus - pred_minus) / (2 * epsilon)
        simple_se <- abs(simple_derivatives) * 0.2
        result$derivative <- simple_derivatives
        result$se_derivative <- simple_se
        result$lower_ci_derivative <- simple_derivatives - 1.96 * simple_se
        result$upper_ci_derivative <- simple_derivatives + 1.96 * simple_se
        result$x_point <- x_values
      })
    }

    return(result)
  }

  calculate_numerical_derivatives <- function(model, newdata, variable_x, model_type) {
    epsilon <- 1e-5 * sd(newdata[[variable_x]], na.rm = TRUE)
    data_plus <- data_minus <- newdata
    data_plus[[variable_x]] <- newdata[[variable_x]] + epsilon
    data_minus[[variable_x]] <- newdata[[variable_x]] - epsilon

    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
      if (model_type == "cox") {
        pred_orig <- predict(model, newdata = newdata, type = "lp", re.form = NA)
        pred_plus <- predict(model, newdata = data_plus, type = "lp", re.form = NA)
        pred_minus <- predict(model, newdata = data_minus, type = "lp", re.form = NA)
      } else {
        pred_orig <- predict(model, newdata = newdata, type = "link", re.form = NA)
        pred_plus <- predict(model, newdata = data_plus, type = "link", re.form = NA)
        pred_minus <- predict(model, newdata = data_minus, type = "link", re.form = NA)
      }
    } else {
      if (model_type == "cox") {
        pred_orig <- predict(model, newdata = newdata, type = "lp")
        pred_plus <- predict(model, newdata = data_plus, type = "lp")
        pred_minus <- predict(model, newdata = data_minus, type = "lp")
      } else if (model_type %in% c("nb", "poisson", "logistic")) {
        pred_orig <- predict(model, newdata = newdata, type = "link")
        pred_plus <- predict(model, newdata = data_plus, type = "link")
        pred_minus <- predict(model, newdata = data_minus, type = "link")
      } else {
        pred_orig <- predict(model, newdata = newdata)
        pred_plus <- predict(model, newdata = data_plus)
        pred_minus <- predict(model, newdata = data_minus)
      }
    }

    linear_derivatives <- (pred_plus - pred_minus) / (2 * epsilon)

    if (model_type %in% c("nb", "poisson", "cox")) {
      response_derivatives <- exp(pred_orig) * linear_derivatives
    } else if (model_type == "logistic") {
      p <- plogis(pred_orig)
      response_derivatives <- p * (1 - p) * linear_derivatives
    } else {
      response_derivatives <- linear_derivatives
    }

    se_derivatives <- abs(response_derivatives) * 0.2
    lower_ci <- response_derivatives - 1.96 * se_derivatives
    upper_ci <- response_derivatives + 1.96 * se_derivatives

    return(list(
      derivative = response_derivatives,
      se_derivative = se_derivatives,
      lower_ci_derivative = lower_ci,
      upper_ci_derivative = upper_ci,
      x_point = newdata[[variable_x]]
    ))
  }

  pool_derivatives_with_rubins_rules <- function(derivative_list, MI_method) {
    n_points <- nrow(derivative_list[[1]])
    pooled_derivatives <- numeric(n_points)
    pooled_var <- numeric(n_points)
    m <- length(derivative_list)

    for (i in 1:n_points) {
      derivs_i <- sapply(derivative_list, function(x) x$derivative[i])
      ses_i <- sapply(derivative_list, function(x) x$se_derivative[i])
      valid_idx <- !is.na(derivs_i) & !is.na(ses_i)

      if (sum(valid_idx) == 0) {
        pooled_derivatives[i] <- NA
        pooled_var[i] <- NA
        next
      }

      derivs_i <- derivs_i[valid_idx]
      ses_i <- ses_i[valid_idx]
      m_valid <- sum(valid_idx)

      pooled_derivatives[i] <- mean(derivs_i, na.rm = TRUE)

      if (MI_method == "Rubin") {
        W <- mean(ses_i^2, na.rm = TRUE)
        B <- var(derivs_i, na.rm = TRUE)
        if (m_valid > 1) {
          pooled_var[i] <- W + (1 + 1/m_valid) * B
        } else {
          pooled_var[i] <- W
        }
      } else if (MI_method == "average") {
        pooled_var[i] <- mean(ses_i^2, na.rm = TRUE)
      }
    }

    pooled_se <- sqrt(pooled_var)
    lower_ci <- pooled_derivatives - 1.96 * pooled_se
    upper_ci <- pooled_derivatives + 1.96 * pooled_se

    return(list(
      derivative = pooled_derivatives,
      se_derivative = pooled_se,
      lower_ci_derivative = lower_ci,
      upper_ci_derivative = upper_ci
    ))
  }

  # ====== Derivative calculation (optional) ======

  if (calculate_derivatives) {
    requireNamespace("rms", quietly = TRUE)

    if (is.null(derivative_points)) {
      derivative_points <- seq(
        from = quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
        to = quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
        length.out = 100
      )
    }

    if (is.null(subgroups)) {
      deriv_data <- data.frame(x_values = derivative_points)
      colnames(deriv_data) <- variable_x
    } else {
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

      if (is.factor(Data_Subset[[subgroups]])) {
        deriv_data[[subgroups]] <- factor(deriv_data[[subgroups]],
                                          levels = levels(Data_Subset[[subgroups]]))
      } else if (subgroup_as_factor) {
        deriv_data[[subgroups]] <- factor(deriv_data[[subgroups]])
        if (!is.null(subgroup_labels)) {
          levels(deriv_data[[subgroups]]) <- subgroup_labels
        }
      }
    }

    if (!is.null(expanded_covariables)) {
      covariate_vars_for_prediction <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
      covariate_vars_for_prediction <- trimws(covariate_vars_for_prediction)
      covariate_vars_for_prediction <- sapply(covariate_vars_for_prediction, function(term) {
        if (grepl("^rcs\\(", term) || grepl("^bs\\(", term) || grepl("^poly\\(", term)) {
          inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
          return(trimws(inside))
        } else {
          return(term)
        }
      })

      for (cov in covariate_vars_for_prediction) {
        if (cov %in% colnames(Data_Subset)) {
          if (is.factor(Data_Subset[[cov]])) {
            deriv_data[[cov]] <- as.factor(names(which.max(table(Data_Subset[[cov]]))))
          } else {
            deriv_data[[cov]] <- median(Data_Subset[[cov]], na.rm = TRUE)
          }
        }
      }
    }

    if (trial_factor == "Yes") {
      deriv_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
    }
    if (followup_offset == "Yes") {
      deriv_data[[followup_col]] <- 365
    }
    if (random_intercept == "Yes" || random_slope == "Yes") {
      deriv_data[[random_intercept_var]] <- names(which.max(table(Data_Subset[[random_intercept_var]])))
    }

    derivative_list <- list()
    derivative_data <- NULL

    if (MI_method == "first") {
      if (derivative_method == "basis") {
        derivative_list[[1]] <- calculate_basis_function_derivatives(
          model = model,
          newdata = deriv_data,
          variable_x = variable_x,
          knot_n = knot_n,
          model_type = model_type
        )
      } else if (derivative_method == "delta") {
        derivative_list[[1]] <- calculate_delta_method_derivatives(
          model = model,
          model_data = Data_Subset,
          newdata = deriv_data,
          variable_x = variable_x,
          model_type = model_type,
          knot_n = knot_n,
          random_intercept = random_intercept,
          formula_obj = formula_obj
        )
      } else if (derivative_method == "numeric") {
        derivative_list[[1]] <- calculate_numerical_derivatives(
          model = model,
          newdata = deriv_data,
          variable_x = variable_x,
          model_type = model_type
        )
      }
      derivative_data <- derivative_list[[1]]

    } else {
      imputation_ids <- unique(Data_Subset[[imp_col]])

      for (imp in imputation_ids) {
        Data_imp <- subset(Data_Subset, get(imp_col) == imp)

        if (random_intercept == "Yes" || random_slope == "Yes") {
          if (model_type == "nb" && requireNamespace("glmmTMB", quietly = TRUE)) {
            model_imp <- glmmTMB::glmmTMB(formula_obj, family = glmmTMB::nbinom2, data = Data_imp)
          } else if (model_type == "poisson") {
            model_imp <- lme4::glmer(formula_obj, family = poisson(link = "log"), data = Data_imp)
          } else if (model_type == "logistic") {
            model_imp <- lme4::glmer(formula_obj, family = binomial(link = "logit"), data = Data_imp)
          } else if (model_type == "lm") {
            model_imp <- lme4::lmer(formula_obj, data = Data_imp)
          } else if (model_type == "cox") {
            model_imp <- coxme::coxme(formula_obj, data = Data_imp)
          }
        } else {
          if (model_type == "nb") {
            model_imp <- MASS::glm.nb(formula_obj, data = Data_imp)
          } else if (model_type == "poisson") {
            model_imp <- glm(formula_obj, family = poisson(link = "log"), data = Data_imp)
          } else if (model_type == "logistic") {
            model_imp <- glm(formula_obj, family = binomial(link = "logit"), data = Data_imp)
          } else if (model_type == "lm") {
            model_imp <- glm(formula_obj, family = gaussian(), data = Data_imp)
          } else if (model_type == "cox") {
            model_imp <- survival::coxph(formula_obj, data = Data_imp)
          }
        }

        if (derivative_method == "basis") {
          derivative_list[[as.character(imp)]] <- calculate_basis_function_derivatives(
            model = model_imp,
            newdata = deriv_data,
            variable_x = variable_x,
            knot_n = knot_n,
            model_type = model_type
          )
        } else if (derivative_method == "delta") {
          derivative_list[[as.character(imp)]] <- calculate_delta_method_derivatives(
            model = model_imp,
            model_data = Data_imp,
            newdata = deriv_data,
            variable_x = variable_x,
            model_type = model_type,
            knot_n = knot_n,
            random_intercept = random_intercept,
            formula_obj = formula_obj
          )
        } else if (derivative_method == "numeric") {
          derivative_list[[as.character(imp)]] <- calculate_numerical_derivatives(
            model = model_imp,
            newdata = deriv_data,
            variable_x = variable_x,
            model_type = model_type
          )
        }
      }

      pooled_results <- pool_derivatives_with_rubins_rules(derivative_list, MI_method)

      derivative_data <- deriv_data
      derivative_data$derivative <- pooled_results$derivative
      derivative_data$se_derivative <- pooled_results$se_derivative
      derivative_data$lower_ci_derivative <- pooled_results$lower_ci_derivative
      derivative_data$upper_ci_derivative <- pooled_results$upper_ci_derivative
      derivative_data$x_point <- deriv_data[[variable_x]]
    }

    # Thresholds based on derivative CI
    if (!is.null(subgroups)) {
      derivative_data$significant_positive <- derivative_data$lower_ci_derivative > 0
      subgroup_levels <- unique(derivative_data[[subgroups]])
      threshold_values <- numeric(length(subgroup_levels))
      names(threshold_values) <- as.character(subgroup_levels)

      for (level in subgroup_levels) {
        subgroup_deriv <- derivative_data[derivative_data[[subgroups]] == level, ]
        subgroup_deriv <- subgroup_deriv[order(subgroup_deriv$x_point), ]
        sig_pos_idx <- which(subgroup_deriv$significant_positive)
        if (length(sig_pos_idx) > 0) {
          threshold_values[as.character(level)] <- subgroup_deriv$x_point[min(sig_pos_idx)]
        } else {
          threshold_values[as.character(level)] <- NA
        }
      }
    } else {
      derivative_data$significant_positive <- derivative_data$lower_ci_derivative > 0
      sorted_deriv <- derivative_data[order(derivative_data$x_point), ]
      sig_pos_idx <- which(sorted_deriv$significant_positive)
      if (length(sig_pos_idx) > 0) {
        threshold <- sorted_deriv$x_point[min(sig_pos_idx)]
      } else {
        threshold <- NA
      }
    }
  }

  # ====== Plotting (prediction curve) ======

  if (is.null(plot_options)) {
    plot_options <- list()
  }

  x_lab <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x

  if (!is.null(plot_options$use_superscript) && plot_options$use_superscript) {
    if (is.character(x_lab)) {
      x_lab <- gsub("\\^9", "\u2079", x_lab)
      x_lab <- gsub("\\^8", "\u2078", x_lab)
      x_lab <- gsub("\\^7", "\u2077", x_lab)
      x_lab <- gsub("\\^6", "\u2076", x_lab)
      x_lab <- gsub("\\^5", "\u2075", x_lab)
      x_lab <- gsub("\\^4", "\u2074", x_lab)
      x_lab <- gsub("\\^3", "\u00B3", x_lab)
      x_lab <- gsub("\\^2", "\u00B2", x_lab)
      x_lab <- gsub("\\^1", "\u00B9", x_lab)
      x_lab <- gsub("\\^0", "\u2070", x_lab)
      x_lab <- gsub("\\^-", "\u207B", x_lab)
    }
  }

  y_lab <- if (!is.null(plot_options$y_lab)) {
    plot_options$y_lab
  } else {
    if (model_type %in% c("nb", "poisson")) {
      "Predicted Rate"
    } else if (model_type == "logistic") {
      "Predicted Probability"
    } else if (model_type == "lm") {
      "Predicted Mean"
    } else {
      "Predicted Hazard Ratio"
    }
  }

  if (!is.null(plot_options$use_superscript) && plot_options$use_superscript) {
    if (is.character(y_lab)) {
      y_lab <- gsub("\\^9", "\u2079", y_lab)
      y_lab <- gsub("\\^8", "\u2078", y_lab)
      y_lab <- gsub("\\^7", "\u2077", y_lab)
      y_lab <- gsub("\\^6", "\u2076", y_lab)
      y_lab <- gsub("\\^5", "\u2075", y_lab)
      y_lab <- gsub("\\^4", "\u2074", y_lab)
      y_lab <- gsub("\\^3", "\u00B3", y_lab)
      y_lab <- gsub("\\^2", "\u00B2", y_lab)
      y_lab <- gsub("\\^1", "\u00B9", y_lab)
      y_lab <- gsub("\\^0", "\u2070", y_lab)
      y_lab <- gsub("\\^-", "\u207B", y_lab)
    }
  }

  if (!is.null(subgroups)) {
    if (MI_method == "average" || MI_method == "Rubin") {
      count_data <- subset(data, get(imp_col) == 1)
      if (!is.null(subgroup_mapping) && !is.null(original_subgroup_values)) {
        orig_subgroup_values_first_imp <- count_data[[subgroups]]
        subgroup_counts <- table(orig_subgroup_values_first_imp)
        names(subgroup_counts) <- subgroup_mapping[names(subgroup_counts)]
      } else {
        subgroup_counts <- table(count_data[[subgroups]])
      }
    } else {
      if (!is.null(subgroup_mapping) && !is.null(original_subgroup_values)) {
        subgroup_counts <- table(original_subgroup_values)
        names(subgroup_counts) <- subgroup_mapping[names(subgroup_counts)]
      } else {
        subgroup_counts <- table(Data_Subset[[subgroups]])
      }
    }

    if (!is.null(plot_options$legend_labels)) {
      legend_labels <- plot_options$legend_labels
      if (!is.null(plot_options$include_counts) && plot_options$include_counts) {
        plot_levels <- levels(factor(pred_data[[subgroups]]))
        if (!is.null(subgroup_mapping)) {
          legend_labels_with_counts <- character(length(plot_levels))
          for (i in seq_along(plot_levels)) {
            level <- plot_levels[i]
            label <- plot_options$legend_labels[i]
            count <- subgroup_counts[level]
            if (is.na(count)) count <- 0
            legend_labels_with_counts[i] <- paste0(label, " (N=", count, ")")
          }
          legend_labels <- legend_labels_with_counts
        } else {
          for (i in seq_along(legend_labels)) {
            count <- subgroup_counts[plot_levels[i]]
            if (is.na(count)) count <- 0
            legend_labels[i] <- paste0(legend_labels[i], " (N=", count, ")")
          }
        }
      }
    } else if (!is.null(plot_options$include_counts) && plot_options$include_counts) {
      legend_labels <- paste0(names(subgroup_counts), " (N=", subgroup_counts, ")")
    } else {
      legend_labels <- levels(factor(pred_data[[subgroups]]))
    }

    plot <- ggplot2::ggplot(pred_data,
                            ggplot2::aes_string(x = variable_x, y = "prediction",
                                                color = subgroups,
                                                fill = subgroups,
                                                linetype = subgroups))
  } else {
    plot <- ggplot2::ggplot(pred_data, ggplot2::aes_string(x = variable_x, y = "prediction"))
    line_color <- if (!is.null(plot_options$colors)) plot_options$colors[1] else "black"
    fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors[1] else "grey70"
  }

  line_size <- if (!is.null(plot_options$line_size)) plot_options$line_size else 1
  ribbon_alpha <- if (!is.null(plot_options$ribbon_alpha)) plot_options$ribbon_alpha else 0.3

  if (!is.null(subgroups)) {
    plot <- plot +
      ggplot2::geom_line(linewidth = line_size) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), alpha = ribbon_alpha)
  } else {
    plot <- plot +
      ggplot2::geom_line(linewidth = line_size, color = line_color) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
        fill = fill_colors,
        color = line_color,
        alpha = ribbon_alpha
      )
  }

  plot <- plot +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab)

  title_args <- list()
  if (!is.null(plot_options$title)) title_args$title <- plot_options$title
  if (!is.null(plot_options$subtitle)) title_args$subtitle <- plot_options$subtitle
  if (!is.null(plot_options$caption)) title_args$caption <- plot_options$caption
  if (length(title_args) > 0) {
    plot <- plot + do.call(ggplot2::labs, title_args)
  }

  if (!is.null(plot_options$use_log_x) && plot_options$use_log_x) {
    if (!is.null(plot_options$x_breaks)) {
      plot <- plot + ggplot2::scale_x_log10(breaks = plot_options$x_breaks)
    } else {
      plot <- plot + ggplot2::scale_x_log10()
    }
  } else if (!is.null(plot_options$x_breaks)) {
    plot <- plot + ggplot2::scale_x_continuous(breaks = plot_options$x_breaks)
  }

  if (!is.null(plot_options$y_breaks)) {
    plot <- plot + ggplot2::scale_y_continuous(
      limits = plot_options$y_limits,
      breaks = plot_options$y_breaks,
      labels = plot_options$y_labels
    )
  }

  coord_args <- list()
  if (!is.null(plot_options$y_limits)) coord_args$ylim <- plot_options$y_limits
  if (!is.null(plot_options$x_limits)) coord_args$xlim <- plot_options$x_limits
  if (length(coord_args) > 0) {
    plot <- plot + do.call(ggplot2::coord_cartesian, coord_args)
  }

  if (!is.null(plot_options$colors) && !is.null(subgroups)) {
    plot <- plot + ggplot2::scale_color_manual(
      values = plot_options$colors,
      labels = legend_labels,
      name = if (!is.null(plot_options$legend_title)) plot_options$legend_title else subgroups
    )
    fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors else plot_options$colors
    plot <- plot + ggplot2::scale_fill_manual(
      values = fill_colors,
      guide = "none"
    )
  }

  if (!is.null(plot_options$line_types) && !is.null(subgroups)) {
    plot <- plot + ggplot2::scale_linetype_manual(
      values = plot_options$line_types,
      labels = legend_labels,
      name = if (!is.null(plot_options$legend_title)) plot_options$legend_title else subgroups
    )
  }

  if (!is.null(plot_options$custom_guides)) {
    plot <- plot + plot_options$custom_guides
  } else if (!is.null(subgroups) && !is.null(plot_options$colors)) {
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

  if (!is.null(plot_options$facet_var)) {
    facet_scales <- if (!is.null(plot_options$facet_scales)) plot_options$facet_scales else "fixed"
    plot <- plot + ggplot2::facet_wrap(as.formula(paste("~", plot_options$facet_var)),
                                       scales = facet_scales)
  }

  theme_settings <- if (!is.null(plot_options$theme_settings)) plot_options$theme_settings else list()
  plot <- plot + ggplot2::theme_minimal(base_size = 10)

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

  if (!is.null(plot_options$legend_position)) {
    default_theme$legend.position <- plot_options$legend_position
  }

  possible_direct_theme_elements <- c(
    "legend.position", "legend.background", "legend.title", "legend.text",
    "legend.key.size", "legend.key.height", "legend.key.width",
    "legend.spacing", "legend.spacing.y", "legend.box.spacing",
    "legend.key.spacing", "legend.box.margin", "legend.margin",
    "axis.text.x", "axis.text.y", "axis.title.x", "axis.title.y", "axis.title",
    "plot.title", "plot.margin", "panel.grid.major", "panel.grid.minor",
    "axis.line", "axis.ticks", "axis.ticks.length", "panel.background", "plot.background"
  )

  for (element_name in possible_direct_theme_elements) {
    direct_element <- paste0("theme_", gsub("\\.", "_", element_name))
    if (!is.null(plot_options[[direct_element]])) {
      default_theme[[element_name]] <- plot_options[[direct_element]]
    }
  }

  for (name in names(theme_settings)) {
    default_theme[[name]] <- theme_settings[[name]]
  }

  plot <- plot + do.call(ggplot2::theme, default_theme)

  if (!is.null(plot_options$custom_scales)) {
    for (scale in plot_options$custom_scales) {
      plot <- plot + scale
    }
  }

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

  if (calculate_derivatives && !is.null(plot_options$show_threshold_lines) && plot_options$show_threshold_lines) {
    if (!is.null(subgroups)) {
      threshold_line_data <- data.frame()
      for (sg_level in names(threshold_values)) {
        threshold_val <- threshold_values[sg_level]
        if (!is.na(threshold_val)) {
          new_threshold_data <- data.frame(
            x = threshold_val,
            subgroup = sg_level,
            stringsAsFactors = FALSE
          )
          names(new_threshold_data)[2] <- subgroups
          threshold_line_data <- rbind(threshold_line_data, new_threshold_data)
        }
      }
      if (nrow(threshold_line_data) > 0) {
        if (is.factor(pred_data[[subgroups]])) {
          threshold_line_data[[subgroups]] <- factor(
            threshold_line_data[[subgroups]],
            levels = levels(pred_data[[subgroups]])
          )
        }
        plot <- plot + ggplot2::geom_vline(
          data = threshold_line_data,
          ggplot2::aes_string(xintercept = "x", color = subgroups),
          linetype = if (!is.null(plot_options$threshold_line_type)) plot_options$threshold_line_type else "dashed",
          linewidth = if (!is.null(plot_options$threshold_line_size)) plot_options$threshold_line_size else 0.5,
          alpha = if (!is.null(plot_options$threshold_line_alpha)) plot_options$threshold_line_alpha else 0.8,
          show.legend = FALSE
        )
      }
    } else {
      if (!is.na(threshold)) {
        plot <- plot + ggplot2::geom_vline(
          xintercept = threshold,
          color = if (!is.null(plot_options$threshold_line_color)) plot_options$threshold_line_color else "black",
          linetype = if (!is.null(plot_options$threshold_line_type)) plot_options$threshold_line_type else "dashed",
          linewidth = if (!is.null(plot_options$threshold_line_size)) plot_options$threshold_line_size else 0.5,
          alpha = if (!is.null(plot_options$threshold_line_alpha)) plot_options$threshold_line_alpha else 0.8
        )
      }
    }
  }

  # Derivative plot (optional)
  derivative_plot <- NULL

  if (calculate_derivatives &&
      !is.null(plot_options$create_derivative_plot) &&
      plot_options$create_derivative_plot) {

    deriv_y_lab <- if (!is.null(plot_options$derivative_y_lab)) {
      plot_options$derivative_y_lab
    } else {
      if (model_type %in% c("nb", "poisson")) {
        "Derivative of Predicted Rate"
      } else if (model_type == "logistic") {
        "Derivative of Predicted Probability"
      } else if (model_type == "lm") {
        "Derivative of Predicted Mean"
      } else {
        "Derivative of Predicted Hazard Ratio"
      }
    }

    if (!is.null(plot_options$use_superscript) && plot_options$use_superscript) {
      if (is.character(deriv_y_lab)) {
        deriv_y_lab <- gsub("\\^9", "\u2079", deriv_y_lab)
        deriv_y_lab <- gsub("\\^8", "\u2078", deriv_y_lab)
        deriv_y_lab <- gsub("\\^7", "\u2077", deriv_y_lab)
        deriv_y_lab <- gsub("\\^6", "\u2076", deriv_y_lab)
        deriv_y_lab <- gsub("\\^5", "\u2075", deriv_y_lab)
        deriv_y_lab <- gsub("\\^4", "\u2074", deriv_y_lab)
        deriv_y_lab <- gsub("\\^3", "\u00B3", deriv_y_lab)
        deriv_y_lab <- gsub("\\^2", "\u00B2", deriv_y_lab)
        deriv_y_lab <- gsub("\\^1", "\u00B9", deriv_y_lab)
        deriv_y_lab <- gsub("\\^0", "\u2070", deriv_y_lab)
        deriv_y_lab <- gsub("\\^-", "\u207B", deriv_y_lab)
      }
    }

    if (!is.null(subgroups)) {
      pred_levels <- levels(pred_data[[subgroups]])
      if (is.numeric(derivative_data[[subgroups]])) {
        numeric_values <- sort(unique(derivative_data[[subgroups]]))
        if (length(numeric_values) == length(pred_levels)) {
          derivative_data[[subgroups]] <- factor(derivative_data[[subgroups]],
                                                 levels = numeric_values,
                                                 labels = pred_levels)
        }
      } else {
        derivative_data[[subgroups]] <- factor(derivative_data[[subgroups]],
                                               levels = pred_levels)
      }
    }

    if (!is.null(subgroups)) {
      deriv_plot <- ggplot2::ggplot(derivative_data,
                                    ggplot2::aes_string(x = "x_point",
                                                        y = "derivative",
                                                        color = subgroups,
                                                        fill = subgroups,
                                                        linetype = subgroups))
    } else {
      deriv_plot <- ggplot2::ggplot(derivative_data,
                                    ggplot2::aes_string(x = "x_point",
                                                        y = "derivative"))
      deriv_line_color <- if (!is.null(plot_options$colors)) plot_options$colors[1] else "black"
      deriv_fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors[1] else "grey70"
    }

    line_size <- if (!is.null(plot_options$line_size)) plot_options$line_size else 1
    ribbon_alpha <- if (!is.null(plot_options$ribbon_alpha)) plot_options$ribbon_alpha else 0.3

    if (!is.null(subgroups)) {
      deriv_plot <- deriv_plot +
        ggplot2::geom_line(linewidth = line_size) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci_derivative,
                                          ymax = upper_ci_derivative),
                             alpha = ribbon_alpha)
    } else {
      deriv_plot <- deriv_plot +
        ggplot2::geom_line(linewidth = line_size, color = deriv_line_color) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci_derivative,
                                          ymax = upper_ci_derivative),
                             fill = deriv_fill_colors,
                             alpha = ribbon_alpha)
    }

    deriv_plot <- deriv_plot +
      ggplot2::xlab(x_lab) +
      ggplot2::ylab(deriv_y_lab)

    deriv_title_args <- list()
    if (!is.null(plot_options$derivative_title)) {
      deriv_title_args$title <- plot_options$derivative_title
    } else {
      method_name <- switch(derivative_method,
                            "basis" = "Basis Function",
                            "delta" = "Delta Method",
                            "numeric" = "Numerical",
                            "Unknown")
      deriv_title_args$title <- paste0(method_name, " Derivative of ", variable_x)
    }
    if (!is.null(plot_options$derivative_subtitle)) {
      deriv_title_args$subtitle <- plot_options$derivative_subtitle
    }
    if (!is.null(plot_options$derivative_caption)) {
      deriv_title_args$caption <- plot_options$derivative_caption
    }
    if (length(deriv_title_args) > 0) {
      deriv_plot <- deriv_plot + do.call(ggplot2::labs, deriv_title_args)
    }

    if (!is.null(plot_options$use_log_x) && plot_options$use_log_x) {
      if (!is.null(plot_options$x_breaks)) {
        deriv_plot <- deriv_plot + ggplot2::scale_x_log10(breaks = plot_options$x_breaks)
      } else {
        deriv_plot <- deriv_plot + ggplot2::scale_x_log10()
      }
    } else if (!is.null(plot_options$x_breaks)) {
      deriv_plot <- deriv_plot + ggplot2::scale_x_continuous(breaks = plot_options$x_breaks)
    }

    if (!is.null(plot_options$derivative_y_breaks)) {
      deriv_plot <- deriv_plot + ggplot2::scale_y_continuous(breaks = plot_options$derivative_y_breaks)
    }

    coord_args_deriv <- list()
    if (!is.null(plot_options$derivative_y_limits)) coord_args_deriv$ylim <- plot_options$derivative_y_limits
    if (!is.null(plot_options$x_limits)) coord_args_deriv$xlim <- plot_options$x_limits
    if (length(coord_args_deriv) > 0) {
      deriv_plot <- deriv_plot + do.call(ggplot2::coord_cartesian, coord_args_deriv)
    }

    if (!is.null(plot_options$colors) && !is.null(subgroups)) {
      if (!is.null(plot_options$legend_labels)) {
        legend_labels_deriv <- plot_options$legend_labels
      } else {
        legend_labels_deriv <- levels(factor(pred_data[[subgroups]]))
      }
      deriv_plot <- deriv_plot + ggplot2::scale_color_manual(
        values = plot_options$colors,
        labels = legend_labels_deriv,
        name = if (!is.null(plot_options$legend_title)) plot_options$legend_title else subgroups
      )
      fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors else plot_options$colors
      deriv_plot <- deriv_plot + ggplot2::scale_fill_manual(
        values = fill_colors,
        guide = "none"
      )
    }

    if (!is.null(plot_options$line_types) && !is.null(subgroups)) {
      deriv_plot <- deriv_plot + ggplot2::scale_linetype_manual(
        values = plot_options$line_types,
        labels = legend_labels,
        name = if (!is.null(plot_options$legend_title)) plot_options$legend_title else subgroups
      )
    }

    deriv_plot <- deriv_plot + ggplot2::theme_minimal(base_size = 10)
    deriv_plot <- deriv_plot + do.call(ggplot2::theme, default_theme)

    deriv_plot <- deriv_plot + ggplot2::geom_hline(
      yintercept = 0,
      color = "gray50",
      linetype = "solid",
      linewidth = 0.5
    )

    derivative_plot <- deriv_plot
  }

  # ====== Return ======

  return(list(
    predictions = pred_data,
    model = model,
    plot = plot,
    plot_data = pred_data,
    prediction_range_values = c(
      quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
      quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE)
    ),
    derivatives = if (calculate_derivatives) derivative_data else NULL,
    derivative_method = if (calculate_derivatives) derivative_method else NULL,
    threshold = if (calculate_derivatives && !is.null(subgroups)) {
      threshold_values
    } else if (calculate_derivatives) {
      threshold
    } else {
      NULL
    },
    derivative_plot = derivative_plot
  ))
}
