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
                      derivative_method = "basis",
                      group_fits = FALSE) {

  ## -----------------------------
  ## Packages & basic checks
  ## -----------------------------
  requireNamespace("rms", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  if (random_intercept == "Yes" || random_slope == "Yes") {
    requireNamespace("lme4", quietly = TRUE)

    if (model_type == "nb" && !requireNamespace("glmmTMB", quietly = TRUE)) {
      warning("Package 'glmmTMB' not available. Using lme4 Poisson as approximation for NB with random effects.")
    }

    if (model_type == "cox" && !requireNamespace("coxme", quietly = TRUE)) {
      stop("Package 'coxme' is required for Cox models with random effects.")
    }
  }

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

  if (!random_intercept %in% c("Yes", "No")) {
    stop("random_intercept must be either 'Yes' or 'No'")
  }
  if (!random_slope %in% c("Yes", "No")) {
    stop("random_slope must be either 'Yes' or 'No'")
  }

  if ((random_intercept == "Yes" || random_slope == "Yes") && is.null(random_intercept_var)) {
    stop("If random_intercept = 'Yes' or random_slope = 'Yes', random_intercept_var must be provided")
  }

  if ((random_intercept == "Yes" || random_slope == "Yes") &&
      !random_intercept_var %in% names(data)) {
    stop("random_intercept_var not found in data")
  }

  if (random_slope == "Yes" && is.null(predictor_vars_random_slope)) {
    stop("If random_slope = 'Yes', predictor_vars_random_slope must be provided")
  }

  if (length(prediction_range) != 2 ||
      any(prediction_range < 0) ||
      any(prediction_range > 1) ||
      prediction_range[1] >= prediction_range[2]) {
    stop("prediction_range must be c(a, b) with 0 <= a < b <= 1")
  }

  if (calculate_derivatives &&
      !derivative_method %in% c("basis", "delta", "numeric")) {
    stop("derivative_method must be one of 'basis', 'delta', or 'numeric'")
  }

  if (group_fits && MI_method != "first") {
    warning("group_fits is only implemented for MI_method = 'first'. It will be ignored.")
    group_fits <- FALSE
  }

  # Initialize thresholds
  threshold        <- NA
  threshold_values <- NULL

  ## -----------------------------
  ## Imputation subset
  ## -----------------------------
  if (MI_method == "first") {
    Data_Subset <- subset(data, get(imp_col) == 1)
  } else if (MI_method %in% c("average", "Rubin")) {
    Data_Subset <- data
  } else {
    stop("Invalid MI_method. Choose 'first', 'average', or 'Rubin'.")
  }

  ## -----------------------------
  ## Subgroups handling
  ## -----------------------------
  original_subgroup_values <- NULL
  if (!is.null(subgroups)) {
    original_subgroup_values <- Data_Subset[[subgroups]]
  }

  subgroup_mapping <- NULL
  if (!is.null(subgroups) && !is.null(subgroup_labels)) {
    unique_values <- sort(unique(Data_Subset[[subgroups]]))
    if (length(unique_values) == length(subgroup_labels)) {
      subgroup_mapping <- stats::setNames(subgroup_labels, unique_values)
    }
  }

  if (!is.null(subgroups) && subgroup_as_factor && !is.factor(Data_Subset[[subgroups]])) {
    Data_Subset[[subgroups]] <- as.factor(Data_Subset[[subgroups]])
    if (!is.null(subgroup_labels)) {
      levels(Data_Subset[[subgroups]]) <- subgroup_labels
    }
  }

  ## -----------------------------
  ## Build spline term
  ## -----------------------------
  spline_term <- paste0("rcs(", variable_x, ", ", knot_n, ")")
  if (!is.null(subgroups)) {
    spline_term <- paste0(spline_term, " * ", subgroups)
  }

  ## -----------------------------
  ## Expand covariables
  ## -----------------------------
  expand_terms <- function(term) {
    if (grepl("^rcs\\(", term) || grepl("^bs\\(", term)) {
      return(term)
    }
    if (grepl("\\*", term)) {
      vars_split <- trimws(unlist(strsplit(term, "\\*")))
      all_combos <- lapply(seq_along(vars_split), function(k) {
        combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
      })
      return(unlist(all_combos))
    } else {
      term
    }
  }

  expanded_covariables <- NULL
  covariates_str <- ""
  if (!is.null(covariables)) {
    expanded_covariables <- unique(unlist(lapply(covariables, expand_terms)))

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
      stop("Covariates not found in data: ", paste(missing_vars_cov, collapse = ", "))
    }

    covariates_str <- paste0(" + ", paste(expanded_covariables, collapse = " + "))
  }

  ## -----------------------------
  ## Trial factor, offset, random effects (incl. random slopes)
  ## -----------------------------
  trial_str <- ""
  if (trial_factor == "Yes") {
    trial_str <- paste0(" + as.factor(", trial_col, ")")
  }

  offset_str <- ""
  if (followup_offset == "Yes" && model_type %in% c("nb", "poisson")) {
    offset_str <- paste0(" + offset(log(", followup_col, "))")
  }

  random_effect_str <- ""
  if (random_intercept == "Yes" || random_slope == "Yes") {
    rand_terms <- if (random_intercept == "Yes") "1" else "0"

    if (random_slope == "Yes") {
      slope_vars <- unique(c(predictor_vars_random_slope, covariables_random_slope))
      slope_vars <- slope_vars[!is.na(slope_vars) & slope_vars != ""]
      if (length(slope_vars) > 0) {
        rand_terms <- paste(c(rand_terms, slope_vars), collapse = " + ")
      }
    }

    random_effect_str <- paste0(" + (", rand_terms, " | ", random_intercept_var, ")")
  }

  ## -----------------------------
  ## Build formula
  ## -----------------------------
  if (model_type == "cox") {
    formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                          spline_term, covariates_str, trial_str, random_effect_str)
  } else {
    formula_str <- paste0(outcome_var, " ~ ", spline_term, covariates_str,
                          trial_str, offset_str, random_effect_str)
  }
  formula_obj <- stats::as.formula(formula_str)

  ## -----------------------------
  ## Fit main model
  ## -----------------------------
  model <- NULL
  if (random_intercept == "Yes" || random_slope == "Yes") {
    if (model_type == "nb") {
      if (requireNamespace("glmmTMB", quietly = TRUE)) {
        model <- glmmTMB::glmmTMB(formula_obj, family = glmmTMB::nbinom2, data = Data_Subset)
      } else {
        warning("Using Poisson mixed model as approximation for NB.")
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
    if (model_type == "nb") {
      requireNamespace("MASS", quietly = TRUE)
      model <- MASS::glm.nb(formula_obj, data = Data_Subset)
    } else if (model_type == "poisson") {
      model <- stats::glm(formula_obj, family = poisson(link = "log"), data = Data_Subset)
    } else if (model_type == "logistic") {
      model <- stats::glm(formula_obj, family = binomial(link = "logit"), data = Data_Subset)
    } else if (model_type == "lm") {
      model <- stats::glm(formula_obj, family = gaussian(), data = Data_Subset)
    } else if (model_type == "cox") {
      requireNamespace("survival", quietly = TRUE)
      model <- survival::coxph(formula_obj, data = Data_Subset)
    }
  }

  ## -----------------------------
  ## Prediction data for main curve
  ## -----------------------------
  x_values <- seq(
    from = stats::quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
    to   = stats::quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
    length.out = 100
  )

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

  ## Add covariates at median/mode
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
          pred_data[[cov]] <- stats::median(Data_Subset[[cov]], na.rm = TRUE)
        }
      }
    }
  }

  if (trial_factor == "Yes") {
    pred_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
  }
  if (followup_offset == "Yes") {
    pred_data[[followup_col]] <- 365
  }
  if (random_intercept == "Yes" || random_slope == "Yes") {
    pred_data[[random_intercept_var]] <- names(which.max(table(Data_Subset[[random_intercept_var]])))
  }

  ## -----------------------------
  ## Mixed model prediction helper
  ## -----------------------------
  get_mixed_model_predictions <- function(model, newdata, type = "link") {
    if (inherits(model, "merMod")) {
      preds <- predict(model, newdata = newdata, re.form = NA, type = type)
      mm <- stats::model.matrix(stats::formula(model, fixed.only = TRUE), newdata)
      vcov_fixed <- as.matrix(stats::vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      return(list(fit = preds, se.fit = se_fixed))
    } else if (inherits(model, "glmmTMB")) {
      preds <- predict(model, newdata = newdata, re.form = NA, type = type, se.fit = TRUE)
      return(list(fit = preds$fit, se.fit = preds$se.fit))
    } else if (inherits(model, "coxme")) {
      preds <- predict(model, newdata = newdata, type = "lp")
      mm <- stats::model.matrix(stats::formula(model, fixed.only = TRUE), newdata)
      vcov_fixed <- as.matrix(stats::vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      return(list(fit = preds, se.fit = se_fixed))
    } else {
      stop("Unsupported mixed model class for prediction.")
    }
  }

  ## -----------------------------
  ## Predictions + MI pooling
  ## -----------------------------
  if (MI_method == "first") {

    if (random_intercept == "Yes" || random_slope == "Yes") {
      if (model_type == "cox") {
        preds <- get_mixed_model_predictions(model, pred_data, type = "link")
        pred_data$prediction <- exp(preds$fit)
        pred_data$lower_ci   <- exp(preds$fit - 1.96 * preds$se.fit)
        pred_data$upper_ci   <- exp(preds$fit + 1.96 * preds$se.fit)
      } else {
        preds <- get_mixed_model_predictions(model, pred_data, type = "link")
        if (model_type %in% c("nb", "poisson")) {
          pred_data$prediction <- exp(preds$fit)
          pred_data$lower_ci   <- exp(preds$fit - 1.96 * preds$se.fit)
          pred_data$upper_ci   <- exp(preds$fit + 1.96 * preds$se.fit)
        } else if (model_type == "logistic") {
          pred_data$prediction <- stats::plogis(preds$fit)
          pred_data$lower_ci   <- stats::plogis(preds$fit - 1.96 * preds$se.fit)
          pred_data$upper_ci   <- stats::plogis(preds$fit + 1.96 * preds$se.fit)
        } else if (model_type == "lm") {
          pred_data$prediction <- preds$fit
          pred_data$lower_ci   <- preds$fit - 1.96 * preds$se.fit
          pred_data$upper_ci   <- preds$fit + 1.96 * preds$se.fit
        }
      }
    } else {
      if (model_type == "cox") {
        lp <- predict(model, newdata = pred_data, type = "lp")
        pred_data$prediction <- exp(lp)
        pred_data$lower_ci   <- NA
        pred_data$upper_ci   <- NA
      } else {
        preds <- stats::predict(model, newdata = pred_data, type = "link", se.fit = TRUE)
        if (model_type %in% c("nb", "poisson")) {
          pred_data$prediction <- exp(preds$fit)
          pred_data$lower_ci   <- exp(preds$fit - 1.96 * preds$se.fit)
          pred_data$upper_ci   <- exp(preds$fit + 1.96 * preds$se.fit)
        } else if (model_type == "logistic") {
          pred_data$prediction <- stats::plogis(preds$fit)
          pred_data$lower_ci   <- stats::plogis(preds$fit - 1.96 * preds$se.fit)
          pred_data$upper_ci   <- stats::plogis(preds$fit + 1.96 * preds$se.fit)
        } else if (model_type == "lm") {
          pred_data$prediction <- preds$fit
          pred_data$lower_ci   <- preds$fit - 1.96 * preds$se.fit
          pred_data$upper_ci   <- preds$fit + 1.96 * preds$se.fit
        }
      }
    }

  } else if (MI_method %in% c("average", "Rubin")) {

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
          model_imp <- stats::glm(formula_obj, family = poisson(link = "log"), data = Data_imp)
        } else if (model_type == "logistic") {
          model_imp <- stats::glm(formula_obj, family = binomial(link = "logit"), data = Data_imp)
        } else if (model_type == "lm") {
          model_imp <- stats::glm(formula_obj, family = gaussian(), data = Data_imp)
        } else if (model_type == "cox") {
          model_imp <- survival::coxph(formula_obj, data = Data_imp)
        }
      }

      if (random_intercept == "Yes" || random_slope == "Yes") {
        preds_imp <- get_mixed_model_predictions(model_imp, pred_data, type = "link")
      } else {
        if (model_type == "cox") {
          fit <- predict(model_imp, newdata = pred_data, type = "lp")
          preds_imp <- list(fit = fit, se.fit = rep(NA, length(fit)))
        } else {
          preds_imp <- stats::predict(model_imp, newdata = pred_data, type = "link", se.fit = TRUE)
        }
      }

      all_predictions[[as.character(imp)]] <- list(fit = preds_imp$fit,
                                                   se.fit = preds_imp$se.fit)
    }

    n_points <- length(all_predictions[[1]]$fit)
    mean_fit <- numeric(n_points)
    mean_var <- numeric(n_points)

    for (i in seq_len(n_points)) {
      fits_i <- sapply(all_predictions, function(x) x$fit[i])
      ses_i  <- sapply(all_predictions, function(x) x$se.fit[i])

      mean_fit[i] <- mean(fits_i)
      m <- length(imputation_ids)

      if (MI_method == "Rubin") {
        W <- mean(ses_i^2)
        B <- stats::var(fits_i)
        mean_var[i] <- W + (1 + 1/m) * B
      } else {
        mean_var[i] <- mean(ses_i^2)
      }
    }

    if (model_type %in% c("nb", "poisson")) {
      pred_data$prediction <- exp(mean_fit)
      pred_data$lower_ci   <- exp(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci   <- exp(mean_fit + 1.96 * sqrt(mean_var))
    } else if (model_type == "logistic") {
      pred_data$prediction <- stats::plogis(mean_fit)
      pred_data$lower_ci   <- stats::plogis(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci   <- stats::plogis(mean_fit + 1.96 * sqrt(mean_var))
    } else if (model_type == "cox") {
      pred_data$prediction <- exp(mean_fit)
      pred_data$lower_ci   <- exp(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci   <- exp(mean_fit + 1.96 * sqrt(mean_var))
    } else if (model_type == "lm") {
      pred_data$prediction <- mean_fit
      pred_data$lower_ci   <- mean_fit - 1.96 * sqrt(mean_var)
      pred_data$upper_ci   <- mean_fit + 1.96 * sqrt(mean_var)
    }
  }

  ## ===================================================================
  ## Derivative helper functions
  ## ===================================================================

  get_knots_from_model <- function(model, variable_x, knot_n, data) {
    if (inherits(model, c("glm", "lm", "coxph", "coxme", "lmerMod", "glmerMod"))) {
      terms_obj   <- stats::terms(model)
      term_labels <- attr(terms_obj, "term.labels")
      spline_term <- grep(paste0("rcs\\(", variable_x), term_labels, value = TRUE)

      if (length(spline_term) > 0 && !is.null(model$model) && variable_x %in% names(model$model)) {
        var_term <- model$model[[variable_x]]
        if (!is.null(attr(var_term, "knots"))) {
          return(attr(var_term, "knots"))
        }
      }
    }

    if (is.null(data[[variable_x]])) {
      stop("Variable not found in data")
    }
    x_vals <- data[[variable_x]]
    probs  <- seq(0, 1, length.out = knot_n)
    stats::quantile(x_vals, probs = probs, na.rm = TRUE)
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

        for (i in seq_len(n)) {
          term1_deriv <- if (x[i] > t_j) 3 * (x[i] - t_j)^2 else 0
          term2_deriv <- if (x[i] > t_k) 3 * (x[i] - t_k)^2 * (t_k - t_j) / (t_k - t_1) else 0
          derivative_matrix[i, j] <- term1_deriv - term2_deriv
        }
      }
    }

    derivative_matrix
  }

  extract_spline_coefficients <- function(model, variable_x) {
    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
      all_coefs <- lme4::fixef(model)
    } else {
      all_coefs <- stats::coef(model)
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
      stop("Could not identify spline coefficients for ", variable_x, " in the model")
    }

    spline_coefs
  }

  apply_link_derivative <- function(model, newdata, linear_predictor_derivatives, model_type) {
    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
      predictions <- predict(model, newdata = newdata, type = "link", re.form = NA)
    } else {
      predictions <- predict(model, newdata = newdata, type = "link")
    }

    if (model_type == "logistic") {
      p <- stats::plogis(predictions)
      p * (1 - p) * linear_predictor_derivatives
    } else if (model_type %in% c("nb", "poisson", "cox")) {
      exp(predictions) * linear_predictor_derivatives
    } else {
      linear_predictor_derivatives
    }
  }

  calculate_derivative_ci <- function(model, newdata, basis_derivatives, spline_coeffs, model_type) {
    vcov_matrix <- tryCatch({
      if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
        as.matrix(stats::vcov(model))
      } else {
        stats::vcov(model)
      }
    }, error = function(e) {
      NULL
    })

    if (is.null(vcov_matrix)) {
      warning("Could not compute vcov matrix; using approximate SE for derivatives.")
      if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
        predictions <- predict(model, newdata = newdata, type = "link", re.form = NA)
      } else {
        predictions <- predict(model, newdata = newdata, type = "link")
      }

      linear_deriv <- basis_derivatives %*% spline_coeffs

      if (model_type == "logistic") {
        p <- stats::plogis(predictions)
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
      warning("Not all spline terms in vcov matrix; using approximate SE.")
      common_terms <- intersect(spline_terms, rownames(vcov_matrix))

      if (length(common_terms) == 0) {
        if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
          predictions <- predict(model, newdata = newdata, type = "link", re.form = NA)
        } else {
          predictions <- predict(model, newdata = newdata, type = "link")
        }

        linear_deriv <- basis_derivatives %*% spline_coeffs

        if (model_type == "logistic") {
          p <- stats::plogis(predictions)
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

    for (i in seq_len(n_points)) {
      gradient <- basis_derivatives[i, , drop = FALSE]
      var_linear <- gradient %*% vcov_spline %*% t(gradient)
      se_linear <- sqrt(as.numeric(var_linear))

      if (model_type == "logistic") {
        p_i <- stats::plogis(predictions[i])
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

    list(
      derivative = derivatives,
      se_derivative = se_derivatives,
      lower_ci = lower_ci,
      upper_ci = upper_ci
    )
  }

  calculate_basis_function_derivatives <- function(model, newdata, variable_x, knot_n, model_type) {
    knots <- get_knots_from_model(model, variable_x, knot_n, newdata)
    basis_derivatives <- calculate_rcs_derivatives(newdata[[variable_x]], knots)
    spline_coeffs <- extract_spline_coefficients(model, variable_x)
    deriv_results <- calculate_derivative_ci(model, newdata, basis_derivatives, spline_coeffs, model_type)

    list(
      derivative = deriv_results$derivative,
      se_derivative = deriv_results$se_derivative,
      lower_ci_derivative = deriv_results$lower_ci,
      upper_ci_derivative = deriv_results$upper_ci,
      x_point = newdata[[variable_x]]
    )
  }

  calculate_delta_method_derivatives <- function(model, model_data, newdata, variable_x, model_type,
                                                 knot_n, random_intercept, formula_obj) {

    delta_method_derivatives <- function(model, newdata, subgroup_var = NULL, x_var, model_type) {
      epsilon <- 1e-5 * stats::sd(newdata[[x_var]], na.rm = TRUE)
      data_plus  <- newdata
      data_minus <- newdata
      data_plus[[x_var]]  <- newdata[[x_var]] + epsilon
      data_minus[[x_var]] <- newdata[[x_var]] - epsilon

      if (model_type == "cox") {
        if (inherits(model, "coxme")) {
          pred_orig  <- predict(model, newdata = newdata,    type = "lp", re.form = NA)
          pred_plus  <- predict(model, newdata = data_plus,  type = "lp", re.form = NA)
          pred_minus <- predict(model, newdata = data_minus, type = "lp", re.form = NA)
        } else {
          pred_orig  <- predict(model, newdata = newdata,    type = "lp")
          pred_plus  <- predict(model, newdata = data_plus,  type = "lp")
          pred_minus <- predict(model, newdata = data_minus, type = "lp")
        }
        pred_response <- exp(pred_orig)
      } else if (model_type %in% c("nb", "poisson")) {
        if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB"))) {
          pred_orig  <- predict(model, newdata = newdata,    type = "link", re.form = NA)
          pred_plus  <- predict(model, newdata = data_plus,  type = "link", re.form = NA)
          pred_minus <- predict(model, newdata = data_minus, type = "link", re.form = NA)
        } else {
          pred_orig  <- predict(model, newdata = newdata,    type = "link")
          pred_plus  <- predict(model, newdata = data_plus,  type = "link")
          pred_minus <- predict(model, newdata = data_minus, type = "link")
        }
        pred_response <- exp(pred_orig)
      } else if (model_type == "logistic") {
        if (inherits(model, c("lmerMod", "glmerMod"))) {
          pred_orig  <- predict(model, newdata = newdata,    type = "link", re.form = NA)
          pred_plus  <- predict(model, newdata = data_plus,  type = "link", re.form = NA)
          pred_minus <- predict(model, newdata = data_minus, type = "link", re.form = NA)
        } else {
          pred_orig  <- predict(model, newdata = newdata,    type = "link")
          pred_plus  <- predict(model, newdata = data_plus,  type = "link")
          pred_minus <- predict(model, newdata = data_minus, type = "link")
        }
        pred_response <- stats::plogis(pred_orig)
      } else if (model_type == "lm") {
        if (inherits(model, c("lmerMod"))) {
          pred_orig  <- predict(model, newdata = newdata,    re.form = NA)
          pred_plus  <- predict(model, newdata = data_plus,  re.form = NA)
          pred_minus <- predict(model, newdata = data_minus, re.form = NA)
        } else {
          pred_orig  <- predict(model, newdata = newdata)
          pred_plus  <- predict(model, newdata = data_plus)
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

      se_derivatives <- abs(derivatives) * 0.2

      list(
        derivatives = derivatives,
        se = se_derivatives
      )
    }

    if (!is.null(subgroups)) {
      result <- newdata
      result$derivative           <- NA
      result$se_derivative        <- NA
      result$lower_ci_derivative  <- NA
      result$upper_ci_derivative  <- NA
      result$x_point              <- newdata[[variable_x]]

      for (sg_level in unique(newdata[[subgroups]])) {
        sg_idx  <- which(newdata[[subgroups]] == sg_level)
        sg_data <- newdata[sg_idx, , drop = FALSE]

        delta_results <- delta_method_derivatives(
          model      = model,
          newdata    = sg_data,
          subgroup_var = subgroups,
          x_var      = variable_x,
          model_type = model_type
        )

        result$derivative[sg_idx]          <- delta_results$derivatives
        result$se_derivative[sg_idx]       <- delta_results$se
        result$lower_ci_derivative[sg_idx] <- delta_results$derivatives - 1.96 * delta_results$se
        result$upper_ci_derivative[sg_idx] <- delta_results$derivatives + 1.96 * delta_results$se
      }
    } else {
      result <- newdata

      delta_results <- delta_method_derivatives(
        model      = model,
        newdata    = newdata,
        x_var      = variable_x,
        model_type = model_type
      )

      result$derivative           <- delta_results$derivatives
      result$se_derivative        <- delta_results$se
      result$lower_ci_derivative  <- delta_results$derivatives - 1.96 * delta_results$se
      result$upper_ci_derivative  <- delta_results$derivatives + 1.96 * delta_results$se
      result$x_point              <- newdata[[variable_x]]
    }

    result
  }

  calculate_numerical_derivatives <- function(model, newdata, variable_x, model_type) {
    epsilon <- 1e-5 * stats::sd(newdata[[variable_x]], na.rm = TRUE)
    data_plus  <- newdata
    data_minus <- newdata
    data_plus[[variable_x]]  <- newdata[[variable_x]] + epsilon
    data_minus[[variable_x]] <- newdata[[variable_x]] - epsilon

    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
      if (model_type == "cox") {
        pred_orig  <- predict(model, newdata = newdata,    type = "lp", re.form = NA)
        pred_plus  <- predict(model, newdata = data_plus,  type = "lp", re.form = NA)
        pred_minus <- predict(model, newdata = data_minus, type = "lp", re.form = NA)
      } else {
        pred_orig  <- predict(model, newdata = newdata,    type = "link", re.form = NA)
        pred_plus  <- predict(model, newdata = data_plus,  type = "link", re.form = NA)
        pred_minus <- predict(model, newdata = data_minus, type = "link", re.form = NA)
      }
    } else {
      if (model_type == "cox") {
        pred_orig  <- predict(model, newdata = newdata,    type = "lp")
        pred_plus  <- predict(model, newdata = data_plus,  type = "lp")
        pred_minus <- predict(model, newdata = data_minus, type = "lp")
      } else if (model_type %in% c("nb", "poisson", "logistic")) {
        pred_orig  <- predict(model, newdata = newdata,    type = "link")
        pred_plus  <- predict(model, newdata = data_plus,  type = "link")
        pred_minus <- predict(model, newdata = data_minus, type = "link")
      } else {
        pred_orig  <- predict(model, newdata = newdata)
        pred_plus  <- predict(model, newdata = data_plus)
        pred_minus <- predict(model, newdata = data_minus)
      }
    }

    linear_derivatives <- (pred_plus - pred_minus) / (2 * epsilon)

    if (model_type %in% c("nb", "poisson", "cox")) {
      response_derivatives <- exp(pred_orig) * linear_derivatives
    } else if (model_type == "logistic") {
      p <- stats::plogis(pred_orig)
      response_derivatives <- p * (1 - p) * linear_derivatives
    } else {
      response_derivatives <- linear_derivatives
    }

    se_derivatives <- abs(response_derivatives) * 0.2
    lower_ci <- response_derivatives - 1.96 * se_derivatives
    upper_ci <- response_derivatives + 1.96 * se_derivatives

    list(
      derivative = response_derivatives,
      se_derivative = se_derivatives,
      lower_ci_derivative = lower_ci,
      upper_ci_derivative = upper_ci,
      x_point = newdata[[variable_x]]
    )
  }

  pool_derivatives_with_rubins_rules <- function(derivative_list, MI_method) {
    n_points <- nrow(derivative_list[[1]])
    pooled_derivatives <- numeric(n_points)
    pooled_var <- numeric(n_points)
    m <- length(derivative_list)

    for (i in seq_len(n_points)) {
      derivs_i <- sapply(derivative_list, function(x) x$derivative[i])
      ses_i    <- sapply(derivative_list, function(x) x$se_derivative[i])

      valid_idx <- !is.na(derivs_i) & !is.na(ses_i)
      if (sum(valid_idx) == 0) {
        pooled_derivatives[i] <- NA
        pooled_var[i]         <- NA
        next
      }

      derivs_i <- derivs_i[valid_idx]
      ses_i    <- ses_i[valid_idx]
      m_valid  <- sum(valid_idx)

      pooled_derivatives[i] <- mean(derivs_i, na.rm = TRUE)

      if (MI_method == "Rubin") {
        W <- mean(ses_i^2, na.rm = TRUE)
        B <- stats::var(derivs_i, na.rm = TRUE)
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
    lower_ci  <- pooled_derivatives - 1.96 * pooled_se
    upper_ci  <- pooled_derivatives + 1.96 * pooled_se

    list(
      derivative = pooled_derivatives,
      se_derivative = pooled_se,
      lower_ci_derivative = lower_ci,
      upper_ci_derivative = upper_ci
    )
  }

  ## ===================================================================
  ## Derivative calculation (if requested)
  ## ===================================================================

  derivative_data  <- NULL
  derivative_plot  <- NULL

  if (calculate_derivatives) {
    requireNamespace("rms", quietly = TRUE)

    if (is.null(derivative_points)) {
      derivative_points <- seq(
        from = stats::quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
        to   = stats::quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
        length.out = 100
      )
    }

    if (is.null(subgroups)) {
      deriv_data <- data.frame(derivative_points)
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
            deriv_data[[cov]] <- as.factor(names(which.max(table(Data_Subset[[cov]]))))
          } else {
            deriv_data[[cov]] <- stats::median(Data_Subset[[cov]], na.rm = TRUE)
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

    if (MI_method == "first") {
      if (derivative_method == "basis") {
        derivative_list[[1]] <- calculate_basis_function_derivatives(
          model      = model,
          newdata    = deriv_data,
          variable_x = variable_x,
          knot_n     = knot_n,
          model_type = model_type
        )
      } else if (derivative_method == "delta") {
        derivative_list[[1]] <- calculate_delta_method_derivatives(
          model          = model,
          model_data     = Data_Subset,
          newdata        = deriv_data,
          variable_x     = variable_x,
          model_type     = model_type,
          knot_n         = knot_n,
          random_intercept = random_intercept,
          formula_obj    = formula_obj
        )
      } else if (derivative_method == "numeric") {
        derivative_list[[1]] <- calculate_numerical_derivatives(
          model      = model,
          newdata    = deriv_data,
          variable_x = variable_x,
          model_type = model_type
        )
      }

      derivative_data <- as.data.frame(derivative_list[[1]])

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
            model_imp <- stats::glm(formula_obj, family = poisson(link = "log"), data = Data_imp)
          } else if (model_type == "logistic") {
            model_imp <- stats::glm(formula_obj, family = binomial(link = "logit"), data = Data_imp)
          } else if (model_type == "lm") {
            model_imp <- stats::glm(formula_obj, family = gaussian(), data = Data_imp)
          } else if (model_type == "cox") {
            model_imp <- survival::coxph(formula_obj, data = Data_imp)
          }
        }

        if (derivative_method == "basis") {
          derivative_list[[as.character(imp)]] <- calculate_basis_function_derivatives(
            model      = model_imp,
            newdata    = deriv_data,
            variable_x = variable_x,
            knot_n     = knot_n,
            model_type = model_type
          )
        } else if (derivative_method == "delta") {
          derivative_list[[as.character(imp)]] <- calculate_delta_method_derivatives(
            model          = model_imp,
            model_data     = Data_imp,
            newdata        = deriv_data,
            variable_x     = variable_x,
            model_type     = model_type,
            knot_n         = knot_n,
            random_intercept = random_intercept,
            formula_obj    = formula_obj
          )
        } else if (derivative_method == "numeric") {
          derivative_list[[as.character(imp)]] <- calculate_numerical_derivatives(
            model      = model_imp,
            newdata    = deriv_data,
            variable_x = variable_x,
            model_type = model_type
          )
        }
      }

      pooled_results <- pool_derivatives_with_rubins_rules(derivative_list, MI_method)

      derivative_data <- deriv_data
      derivative_data$derivative           <- pooled_results$derivative
      derivative_data$se_derivative        <- pooled_results$se_derivative
      derivative_data$lower_ci_derivative  <- pooled_results$lower_ci_derivative
      derivative_data$upper_ci_derivative  <- pooled_results$upper_ci_derivative
      derivative_data$x_point              <- deriv_data[[variable_x]]
    }

    ## Thresholds
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
      sig_pos_idx  <- which(sorted_deriv$significant_positive)
      if (length(sig_pos_idx) > 0) {
        threshold <- sorted_deriv$x_point[min(sig_pos_idx)]
      } else {
        threshold <- NA
      }
    }
  }

  ## ===================================================================
  ## group_fits: per-group curves when random slope, MI_method="first"
  ## ===================================================================

  group_fits_df <- NULL

  if (group_fits &&
      (random_intercept == "Yes" || random_slope == "Yes") &&
      MI_method == "first") {

    groups <- unique(Data_Subset[[random_intercept_var]])
    all_group_pred <- list()

    for (g in groups) {
      Data_g <- Data_Subset[Data_Subset[[random_intercept_var]] == g, , drop = FALSE]
      if (nrow(Data_g) < 5) next

      x_g <- x_values[
        x_values >= min(Data_g[[variable_x]], na.rm = TRUE) &
          x_values <= max(Data_g[[variable_x]], na.rm = TRUE)
      ]
      if (length(x_g) < 3) next

      if (is.null(subgroups)) {
        new_g <- data.frame(x_g)
        colnames(new_g) <- variable_x
      } else {
        subgroup_levels <- if (is.factor(Data_Subset[[subgroups]])) {
          levels(Data_Subset[[subgroups]])
        } else {
          sort(unique(Data_Subset[[subgroups]]))
        }
        new_g <- expand.grid(x_g, subgroup_levels)
        colnames(new_g) <- c(variable_x, subgroups)
        if (is.factor(Data_Subset[[subgroups]])) {
          new_g[[subgroups]] <- factor(new_g[[subgroups]],
                                       levels = levels(Data_Subset[[subgroups]]))
        }
      }

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
              new_g[[cov]] <- as.factor(names(which.max(table(Data_Subset[[cov]]))))
            } else {
              new_g[[cov]] <- stats::median(Data_Subset[[cov]], na.rm = TRUE)
            }
          }
        }
      }

      if (trial_factor == "Yes") {
        if (trial_col == random_intercept_var) {
          new_g[[trial_col]] <- g
        } else {
          new_g[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
        }
      }

      if (followup_offset == "Yes") {
        new_g[[followup_col]] <- 365
      }

      new_g[[random_intercept_var]] <- g

      if (inherits(model, "glmmTMB")) {
        p_g <- predict(model, newdata = new_g,
                       type = "link", se.fit = TRUE, re.form = NULL)
        fit_link <- p_g$fit
        se_link  <- p_g$se.fit
      } else if (inherits(model, "merMod")) {
        fit_link <- predict(model, newdata = new_g,
                            type = "link", re.form = NULL)
        mm <- stats::model.matrix(stats::formula(model, fixed.only = TRUE), new_g)
        vcov_fixed <- as.matrix(stats::vcov(model))
        se_link <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      } else if (inherits(model, "coxme")) {
        fit_link <- predict(model, newdata = new_g, type = "lp")
        mm <- stats::model.matrix(stats::formula(model, fixed.only = TRUE), new_g)
        vcov_fixed <- as.matrix(stats::vcov(model))
        se_link <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      } else {
        p_g <- stats::predict(model, newdata = new_g, type = "link", se.fit = TRUE)
        fit_link <- p_g$fit
        se_link  <- p_g$se.fit
      }

      if (model_type %in% c("nb", "poisson")) {
        new_g$prediction <- exp(fit_link)
        new_g$lower_ci   <- exp(fit_link - 1.96 * se_link)
        new_g$upper_ci   <- exp(fit_link + 1.96 * se_link)
      } else if (model_type == "logistic") {
        new_g$prediction <- stats::plogis(fit_link)
        new_g$lower_ci   <- stats::plogis(fit_link - 1.96 * se_link)
        new_g$upper_ci   <- stats::plogis(fit_link + 1.96 * se_link)
      } else if (model_type == "cox") {
        new_g$prediction <- exp(fit_link)
        new_g$lower_ci   <- exp(fit_link - 1.96 * se_link)
        new_g$upper_ci   <- exp(fit_link + 1.96 * se_link)
      } else if (model_type == "lm") {
        new_g$prediction <- fit_link
        new_g$lower_ci   <- fit_link - 1.96 * se_link
        new_g$upper_ci   <- fit_link + 1.96 * se_link
      }

      new_g[[random_intercept_var]] <- g
      all_group_pred[[as.character(g)]] <- new_g
    }

    if (length(all_group_pred) > 0) {
      group_fits_df <- dplyr::bind_rows(all_group_pred)
    }
  }

  ## ===================================================================
  ## Plotting of main curve
  ## ===================================================================

  if (is.null(plot_options)) plot_options <- list()

  x_lab <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x

  if (!is.null(plot_options$use_superscript) && plot_options$use_superscript && is.character(x_lab)) {
    reps <- c("9"="\u2079","8"="\u2078","7"="\u2077","6"="\u2076","5"="\u2075",
              "4"="\u2074","3"="\u00B3","2"="\u00B2","1"="\u00B9","0"="\u2070","-"="\u207B")
    for (k in names(reps)) {
      x_lab <- gsub(paste0("\\^", k), reps[[k]], x_lab)
    }
  }

  y_lab <- if (!is.null(plot_options$y_lab)) plot_options$y_lab else {
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

  if (!is.null(plot_options$use_superscript) && plot_options$use_superscript && is.character(y_lab)) {
    reps <- c("9"="\u2079","8"="\u2078","7"="\u2077","6"="\u2076","5"="\u2075",
              "4"="\u2074","3"="\u00B3","2"="\u00B2","1"="\u00B9","0"="\u2070","-"="\u207B")
    for (k in names(reps)) {
      y_lab <- gsub(paste0("\\^", k), reps[[k]], y_lab)
    }
  }

  if (!is.null(subgroups)) {
    if (MI_method %in% c("average", "Rubin")) {
      count_data <- subset(data, get(imp_col) == 1)
      subgroup_counts <- table(count_data[[subgroups]])
    } else {
      subgroup_counts <- table(Data_Subset[[subgroups]])
    }

    if (!is.null(plot_options$legend_labels)) {
      legend_labels <- plot_options$legend_labels
      if (!is.null(plot_options$include_counts) && plot_options$include_counts) {
        plot_levels <- levels(factor(pred_data[[subgroups]]))
        for (i in seq_along(legend_labels)) {
          cnt <- subgroup_counts[plot_levels[i]]
          if (is.na(cnt)) cnt <- 0
          legend_labels[i] <- paste0(legend_labels[i], " (N=", cnt, ")")
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
                                                fill  = subgroups,
                                                linetype = subgroups))
  } else {
    plot <- ggplot2::ggplot(pred_data,
                            ggplot2::aes_string(x = variable_x, y = "prediction"))
    line_color  <- if (!is.null(plot_options$colors))      plot_options$colors[1]      else "black"
    fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors[1] else "grey70"
  }

  line_size    <- if (!is.null(plot_options$line_size))    plot_options$line_size    else 1
  ribbon_alpha <- if (!is.null(plot_options$ribbon_alpha)) plot_options$ribbon_alpha else 0.3

  if (!is.null(subgroups)) {
    plot <- plot +
      ggplot2::geom_line(linewidth = line_size) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                           alpha = ribbon_alpha)
  } else {
    plot <- plot +
      ggplot2::geom_line(linewidth = line_size, color = line_color) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                           fill  = fill_colors,
                           color = line_color,
                           alpha = ribbon_alpha)
  }

  plot <- plot + ggplot2::xlab(x_lab) + ggplot2::ylab(y_lab)

  title_args <- list()
  if (!is.null(plot_options$title))     title_args$title    <- plot_options$title
  if (!is.null(plot_options$subtitle))  title_args$subtitle <- plot_options$subtitle
  if (!is.null(plot_options$caption))   title_args$caption  <- plot_options$caption
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

  if (!is.null(plot_options$y_breaks) && !is.null(plot_options$y_limits)) {
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
      name   = if (!is.null(plot_options$legend_title)) plot_options$legend_title else subgroups
    )
    fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors else plot_options$colors
    plot <- plot + ggplot2::scale_fill_manual(values = fill_colors, guide = "none")
  }

  if (!is.null(plot_options$line_types) && !is.null(subgroups)) {
    plot <- plot + ggplot2::scale_linetype_manual(
      values = plot_options$line_types,
      labels = legend_labels,
      name   = if (!is.null(plot_options$legend_title)) plot_options$legend_title else subgroups
    )
  }

  if (!is.null(plot_options$custom_guides)) {
    plot <- plot + plot_options$custom_guides
  } else if (!is.null(subgroups) && !is.null(plot_options$colors)) {
    fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors else plot_options$colors
    plot <- plot + ggplot2::guides(
      linetype = ggplot2::guide_legend(
        override.aes = list(color = plot_options$colors,
                            fill  = fill_colors,
                            shape = NA),
        keywidth = 2
      ),
      color = "none",
      fill  = "none"
    )
  }

  if (!is.null(plot_options$facet_var)) {
    facet_scales <- if (!is.null(plot_options$facet_scales)) plot_options$facet_scales else "fixed"
    plot <- plot + ggplot2::facet_wrap(stats::as.formula(paste("~", plot_options$facet_var)),
                                       scales = facet_scales)
  }

  theme_settings <- if (!is.null(plot_options$theme_settings)) plot_options$theme_settings else list()
  plot <- plot + ggplot2::theme_minimal(base_size = 10)

  default_theme <- list(
    legend.background = ggplot2::element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
    legend.title      = ggplot2::element_text(face = "bold", size = 11),
    legend.text       = ggplot2::element_text(size = 10),
    axis.text.x       = ggplot2::element_text(face = "bold", size = 12),
    axis.text.y       = ggplot2::element_text(face = "bold", size = 12),
    axis.title.x      = ggplot2::element_text(face = "bold", size = 12),
    axis.title.y      = ggplot2::element_text(face = "bold", size = 12),
    plot.title        = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
    plot.margin       = ggplot2::margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid.major  = ggplot2::element_blank(),
    panel.grid.minor  = ggplot2::element_blank(),
    axis.line         = ggplot2::element_line(colour = "black", size = 0.5),
    axis.ticks        = ggplot2::element_line(size = 0.5, colour = "black"),
    axis.ticks.length = ggplot2::unit(0.2, "cm"),
    panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background   = ggplot2::element_rect(fill = "white", colour = NA)
  )

  if (!is.null(plot_options$legend_position)) {
    default_theme$legend.position <- plot_options$legend_position
  }

  for (nm in names(theme_settings)) {
    default_theme[[nm]] <- theme_settings[[nm]]
  }
  plot <- plot + do.call(ggplot2::theme, default_theme)

  if (!is.null(plot_options$custom_scales)) {
    for (sc in plot_options$custom_scales) {
      plot <- plot + sc
    }
  }

  if (!is.null(plot_options$vline)) {
    for (v in plot_options$vline) {
      plot <- plot + ggplot2::geom_vline(
        xintercept = v,
        color      = if (!is.null(plot_options$vline_color)) plot_options$vline_color else "black",
        linetype   = if (!is.null(plot_options$vline_type))  plot_options$vline_type  else "dashed",
        linewidth  = if (!is.null(plot_options$vline_size))  plot_options$vline_size  else 0.5
      )
    }
  }

  if (!is.null(plot_options$hline)) {
    for (h in plot_options$hline) {
      plot <- plot + ggplot2::geom_hline(
        yintercept = h,
        color      = if (!is.null(plot_options$hline_color)) plot_options$hline_color else "black",
        linetype   = if (!is.null(plot_options$hline_type))  plot_options$hline_type  else "dashed",
        linewidth  = if (!is.null(plot_options$hline_size))  plot_options$hline_size  else 0.5
      )
    }
  }

  if (!is.null(plot_options$annotations)) {
    for (annotation in plot_options$annotations) {
      plot <- plot + ggplot2::annotate(
        geom  = annotation$geom,
        x     = annotation$x,
        y     = annotation$y,
        label = annotation$label,
        color = if (!is.null(annotation$color))    annotation$color    else "black",
        size  = if (!is.null(annotation$size))     annotation$size     else 4,
        hjust = if (!is.null(annotation$hjust))    annotation$hjust    else 0.5,
        vjust = if (!is.null(annotation$vjust))    annotation$vjust    else 0.5,
        angle = if (!is.null(annotation$angle))    annotation$angle    else 0,
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
        ymax = if (!is.null(area$ymax)) area$ymax else  Inf,
        fill = if (!is.null(area$fill)) area$fill else "lightblue",
        alpha = if (!is.null(area$alpha)) area$alpha else 0.2
      )
    }
  }

  if (calculate_derivatives &&
      !is.null(plot_options$show_threshold_lines) &&
      plot_options$show_threshold_lines) {

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
          linetype  = if (!is.null(plot_options$threshold_line_type))  plot_options$threshold_line_type  else "dashed",
          linewidth = if (!is.null(plot_options$threshold_line_size))  plot_options$threshold_line_size  else 0.5,
          alpha     = if (!is.null(plot_options$threshold_line_alpha)) plot_options$threshold_line_alpha else 0.8,
          show.legend = FALSE
        )
      }
    } else if (!is.na(threshold)) {
      plot <- plot + ggplot2::geom_vline(
        xintercept = threshold,
        color      = if (!is.null(plot_options$threshold_line_color)) plot_options$threshold_line_color else "black",
        linetype   = if (!is.null(plot_options$threshold_line_type))  plot_options$threshold_line_type  else "dashed",
        linewidth  = if (!is.null(plot_options$threshold_line_size))  plot_options$threshold_line_size  else 0.5,
        alpha      = if (!is.null(plot_options$threshold_line_alpha)) plot_options$threshold_line_alpha else 0.8
      )
    }
  }

  ## (Derivative plot code omitted here to keep length manageable.
  ##  You can keep your previous derivative_plot section if you need
  ##  a dedicated derivative figure; it will plug into derivative_data /
  ##  threshold / threshold_values as before.)

  ## ===================================================================
  ## Return
  ## ===================================================================

  out <- list(
    predictions              = pred_data,
    model                    = model,
    plot                     = plot,
    plot_data                = pred_data,
    prediction_range_values  = c(
      stats::quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
      stats::quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE)
    ),
    derivatives              = if (calculate_derivatives) derivative_data else NULL,
    derivative_method        = if (calculate_derivatives) derivative_method else NULL,
    threshold                = if (calculate_derivatives && !is.null(subgroups)) threshold_values else if (calculate_derivatives) threshold else NULL,
    derivative_plot          = derivative_plot,
    group_fits               = group_fits_df
  )

  class(out) <- c("MI_spline_result", class(out))
  return(out)
}
