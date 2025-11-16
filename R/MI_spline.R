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
                      derivative_method = "numeric",
                      group_fits = FALSE) {

  #------------------------------
  # 0. Required packages
  #------------------------------
  requireNamespace("rms", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  if (random_intercept == "Yes" || random_slope == "Yes") {
    requireNamespace("lme4", quietly = TRUE)
    if (model_type == "nb" && !requireNamespace("glmmTMB", quietly = TRUE)) {
      warning("Package 'glmmTMB' not available. Using Poisson mixed model as approx NB.")
    }
    if (model_type == "cox" && !requireNamespace("coxme", quietly = TRUE)) {
      stop("Package 'coxme' is required for Cox mixed models.")
    }
  }

  #------------------------------
  # 1. Basic checks
  #------------------------------
  if (!model_type %in% c("nb", "poisson", "cox", "logistic", "lm")) {
    stop("model_type must be one of 'nb', 'poisson', 'cox', 'logistic', 'lm'")
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
    stop("random_intercept must be 'Yes' or 'No'")
  }
  if (!random_slope %in% c("Yes", "No")) {
    stop("random_slope must be 'Yes' or 'No'")
  }

  if ((random_intercept == "Yes" || random_slope == "Yes") && is.null(random_intercept_var)) {
    stop("If random_intercept = 'Yes' or random_slope = 'Yes', random_intercept_var must be provided")
  }
  if ((random_intercept == "Yes" || random_scept == "Yes") &&
      !random_intercept_var %in% names(data)) {
    stop("random_intercept_var not found in data")
  }

  if (random_slope == "Yes") {
    slope_vars <- unique(c(predictor_vars_random_slope, covariables_random_slope))
    slope_vars <- slope_vars[!is.na(slope_vars) & slope_vars != ""]
    if (length(slope_vars) == 0) {
      stop("random_slope = 'Yes' but no predictor_vars_random_slope/covariables_random_slope provided")
    }
    missing_slope <- slope_vars[!slope_vars %in% names(data)]
    if (length(missing_slope) > 0) {
      stop("Random slope variables not found in data: ",
           paste(missing_slope, collapse = ", "))
    }
  }

  if (length(prediction_range) != 2 ||
      any(prediction_range < 0) || any(prediction_range > 1) ||
      prediction_range[1] >= prediction_range[2]) {
    stop("prediction_range must be c(q_low, q_high) with 0 < q_low < q_high < 1")
  }

  if (calculate_derivatives && derivative_method != "numeric") {
    stop("In this version, derivative_method must be 'numeric' if calculate_derivatives = TRUE")
  }

  # thresholds (for derivative-based thresholds if used later)
  threshold <- NA
  threshold_values <- NULL

  #------------------------------
  # 2. Select imputations
  #------------------------------
  if (MI_method == "first") {
    Data_Subset <- subset(data, get(imp_col) == 1)
  } else if (MI_method %in% c("average", "Rubin")) {
    Data_Subset <- data
  } else {
    stop("MI_method must be 'first', 'average', or 'Rubin'")
  }

  # store original subgroup values
  original_subgroup_values <- NULL
  if (!is.null(subgroups)) {
    original_subgroup_values <- Data_Subset[[subgroups]]
  }

  # mapping of original subgroup values to labels
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

  #------------------------------
  # 3. Build spline term & covariables
  #------------------------------
  spline_term <- paste0("rcs(", variable_x, ", ", knot_n, ")")
  if (!is.null(subgroups)) {
    spline_term <- paste0(spline_term, " * ", subgroups)
  }

  expand_terms <- function(term) {
    if (grepl("^rcs\\(", term) || grepl("^bs\\(", term)) return(term)
    if (grepl("\\*", term)) {
      vars_split <- unlist(strsplit(term, "\\*"))
      vars_split <- trimws(vars_split)
      all_combinations <- lapply(1:length(vars_split), function(k) {
        combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
      })
      unlist(all_combinations)
    } else term
  }

  expanded_covariables <- NULL
  if (!is.null(covariables)) {
    expanded_covariables <- unique(unlist(lapply(covariables, expand_terms)))

    extract_variables_cov <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^bs\\(", term)) {
        inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
        return(trimws(inside))
      } else {
        term
      }
    }

    all_base_vars_cov <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
    all_base_vars_cov <- trimws(all_base_vars_cov)
    all_base_vars_cov <- sapply(all_base_vars_cov, extract_variables_cov)

    missing_vars_cov <- all_base_vars_cov[!all_base_vars_cov %in% names(Data_Subset)]
    if (length(missing_vars_cov) > 0) {
      stop("Covariates not found in data: ", paste(missing_vars_cov, collapse = ", "))
    }

    covariables_str <- paste0(" + ", paste(expanded_covariables, collapse = " + "))
  } else {
    covariables_str <- ""
  }

  # trial factor
  trial_str <- ""
  if (trial_factor == "Yes") {
    trial_str <- paste0(" + as.factor(", trial_col, ")")
  }

  # offset
  offset_str <- ""
  if (followup_offset == "Yes" && model_type %in% c("nb", "poisson")) {
    offset_str <- paste0(" + offset(log(", followup_col, "))")
  }

  #------------------------------
  # 4. Random effects structure
  #------------------------------
  random_effect_str <- ""
  slope_vars <- NULL
  if (random_slope == "Yes") {
    slope_vars <- unique(c(predictor_vars_random_slope, covariables_random_slope))
    slope_vars <- slope_vars[!is.na(slope_vars) & slope_vars != ""]
  }

  if (random_intercept == "Yes" && random_slope == "Yes") {
    slope_part <- paste(c("1", slope_vars), collapse = " + ")
    random_effect_str <- paste0(" + (", slope_part, " | ", random_intercept_var, ")")
  } else if (random_slope == "Yes") {
    slope_part <- paste(slope_vars, collapse = " + ")
    random_effect_str <- paste0(" + (0 + ", slope_part, " | ", random_intercept_var, ")")
  } else if (random_intercept == "Yes") {
    random_effect_str <- paste0(" + (1 | ", random_intercept_var, ")")
  }

  #------------------------------
  # 5. Build formula and fit model
  #------------------------------
  if (model_type == "cox") {
    formula_str <- paste0(
      "Surv(", time_col, ", ", event_col, ") ~ ",
      spline_term, covariables_str, trial_str, random_effect_str
    )
  } else {
    formula_str <- paste0(
      outcome_var, " ~ ",
      spline_term, covariables_str, trial_str, offset_str, random_effect_str
    )
  }
  formula_obj <- as.formula(formula_str)

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

  #------------------------------
  # 6. Prediction grid (population level)
  #------------------------------
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

    pred_data <- expand.grid(x_values, subgroup_levels)
    colnames(pred_data) <- c(variable_x, subgroups)

    if (is.factor(Data_Subset[[subgroups]])) {
      pred_data[[subgroups]] <- factor(pred_data[[subgroups]],
                                       levels = levels(Data_Subset[[subgroups]]))
    } else if (subgroup_as_factor) {
      pred_data[[subgroups]] <- factor(pred_data[[subgroups]])
      if (!is.null(subgroup_labels)) {
        levels(pred_data[[subgroups]]) <- subgroup_labels
      }
    }
  }

  # fill covariates in pred_data
  if (!is.null(expanded_covariables)) {
    extract_variables_cov <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^bs\\(", term)) {
        inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
        return(trimws(inside))
      } else {
        term
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

  # helper: mixed model predictions with SE on link scale
  get_mixed_model_predictions <- function(model, newdata, type = "link") {
    if (inherits(model, "merMod")) {
      preds <- predict(model, newdata = newdata, re.form = NA, type = type)
      mm <- model.matrix(stats::formula(model, fixed.only = TRUE), newdata)
      vcov_fixed <- as.matrix(stats::vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      list(fit = preds, se.fit = se_fixed)
    } else if (inherits(model, "glmmTMB")) {
      preds <- predict(model, newdata = newdata, re.form = NA, type = type, se.fit = TRUE)
      list(fit = preds$fit, se.fit = preds$se.fit)
    } else if (inherits(model, "coxme")) {
      preds <- predict(model, newdata = newdata, type = "lp")
      mm <- model.matrix(stats::formula(model, fixed.only = TRUE), newdata)
      vcov_fixed <- as.matrix(stats::vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      list(fit = preds, se.fit = se_fixed)
    } else {
      stop("Unsupported mixed model class for get_mixed_model_predictions")
    }
  }

  #------------------------------
  # 7. Predictions + MI pooling
  #------------------------------
  if (MI_method == "first") {
    if (random_intercept == "Yes" || random_slope == "Yes") {
      if (model_type == "cox") {
        pr <- get_mixed_model_predictions(model, pred_data)
        pred_data$prediction <- exp(pr$fit)
        pred_data$lower_ci   <- exp(pr$fit - 1.96*pr$se.fit)
        pred_data$upper_ci   <- exp(pr$fit + 1.96*pr$se.fit)
      } else {
        pr <- get_mixed_model_predictions(model, pred_data)
        if (model_type %in% c("nb", "poisson")) {
          pred_data$prediction <- exp(pr$fit)
          pred_data$lower_ci   <- exp(pr$fit - 1.96*pr$se.fit)
          pred_data$upper_ci   <- exp(pr$fit + 1.96*pr$se.fit)
        } else if (model_type == "logistic") {
          pred_data$prediction <- stats::plogis(pr$fit)
          pred_data$lower_ci   <- stats::plogis(pr$fit - 1.96*pr$se.fit)
          pred_data$upper_ci   <- stats::plogis(pr$fit + 1.96*pr$se.fit)
        } else if (model_type == "lm") {
          pred_data$prediction <- pr$fit
          pred_data$lower_ci   <- pr$fit - 1.96*pr$se.fit
          pred_data$upper_ci   <- pr$fit + 1.96*pr$se.fit
        }
      }
    } else {
      # standard models
      if (model_type == "cox") {
        fit_link <- predict(model, newdata = pred_data, type = "lp")
        pred_data$prediction <- exp(fit_link)
        pred_data$lower_ci   <- NA
        pred_data$upper_ci   <- NA
      } else {
        pr <- predict(model, newdata = pred_data, type = "link", se.fit = TRUE)
        if (model_type %in% c("nb", "poisson")) {
          pred_data$prediction <- exp(pr$fit)
          pred_data$lower_ci   <- exp(pr$fit - 1.96*pr$se.fit)
          pred_data$upper_ci   <- exp(pr$fit + 1.96*pr$se.fit)
        } else if (model_type == "logistic") {
          pred_data$prediction <- stats::plogis(pr$fit)
          pred_data$lower_ci   <- stats::plogis(pr$fit - 1.96*pr$se.fit)
          pred_data$upper_ci   <- stats::plogis(pr$fit + 1.96*pr$se.fit)
        } else if (model_type == "lm") {
          pred_data$prediction <- pr$fit
          pred_data$lower_ci   <- pr$fit - 1.96*pr$se.fit
          pred_data$upper_ci   <- pr$fit + 1.96*pr$se.fit
        }
      }
    }
  } else if (MI_method %in% c("average", "Rubin")) {
    imps <- unique(Data_Subset[[imp_col]])
    all_preds <- vector("list", length(imps))
    names(all_preds) <- as.character(imps)

    for (imp in imps) {
      Data_imp <- subset(Data_Subset, get(imp_col) == imp)

      # fit
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

      # predict on link scale
      if (random_intercept == "Yes" || random_slope == "Yes") {
        pr <- get_mixed_model_predictions(model_imp, pred_data)
      } else {
        if (model_type == "cox") {
          pr <- list(
            fit   = predict(model_imp, newdata = pred_data, type = "lp"),
            se.fit = rep(NA, nrow(pred_data))
          )
        } else {
          pr <- predict(model_imp, newdata = pred_data, type = "link", se.fit = TRUE)
        }
      }

      all_preds[[as.character(imp)]] <- list(fit = pr$fit, se.fit = pr$se.fit)
    }

    n_points <- length(all_preds[[1]]$fit)
    mean_fit <- numeric(n_points)
    mean_var <- numeric(n_points)

    for (i in seq_len(n_points)) {
      fi <- sapply(all_preds, function(z) z$fit[i])
      sei <- sapply(all_preds, function(z) z$se.fit[i])
      m <- length(imps)

      mean_fit[i] <- mean(fi, na.rm = TRUE)
      if (MI_method == "Rubin") {
        W <- mean(sei^2, na.rm = TRUE)
        B <- stats::var(fi, na.rm = TRUE)
        mean_var[i] <- if (m > 1) W + (1 + 1/m)*B else W
      } else {
        mean_var[i] <- mean(sei^2, na.rm = TRUE)
      }
    }

    if (model_type %in% c("nb", "poisson")) {
      pred_data$prediction <- exp(mean_fit)
      pred_data$lower_ci   <- exp(mean_fit - 1.96*sqrt(mean_var))
      pred_data$upper_ci   <- exp(mean_fit + 1.96*sqrt(mean_var))
    } else if (model_type == "logistic") {
      pred_data$prediction <- stats::plogis(mean_fit)
      pred_data$lower_ci   <- stats::plogis(mean_fit - 1.96*sqrt(mean_var))
      pred_data$upper_ci   <- stats::plogis(mean_fit + 1.96*sqrt(mean_var))
    } else if (model_type == "cox") {
      pred_data$prediction <- exp(mean_fit)
      pred_data$lower_ci   <- exp(mean_fit - 1.96*sqrt(mean_var))
      pred_data$upper_ci   <- exp(mean_fit + 1.96*sqrt(mean_var))
    } else if (model_type == "lm") {
      pred_data$prediction <- mean_fit
      pred_data$lower_ci   <- mean_fit - 1.96*sqrt(mean_var)
      pred_data$upper_ci   <- mean_fit + 1.96*sqrt(mean_var)
    }
  }

  #------------------------------
  # 8. Optional: per-group fits with random slope
  #------------------------------
  group_predictions <- NULL
  group_plot <- NULL

  if (group_fits &&
      MI_method == "first" &&
      random_slope == "Yes" &&
      !is.null(random_intercept_var) &&
      random_intercept_var %in% names(Data_Subset)) {

    group_levels <- sort(unique(Data_Subset[[random_intercept_var]]))

    group_pred_list <- lapply(group_levels, function(g) {
      newdat <- pred_data
      newdat[[random_intercept_var]] <- g

      if (inherits(model, "merMod")) {
        if (model_type == "cox") {
          lp <- predict(model, newdata = newdat, type = "lp", re.form = NULL)
          newdat$prediction_group <- exp(lp)
        } else if (model_type %in% c("nb", "poisson")) {
          lp <- predict(model, newdata = newdat, type = "link", re.form = NULL)
          newdat$prediction_group <- exp(lp)
        } else if (model_type == "logistic") {
          lp <- predict(model, newdata = newdat, type = "link", re.form = NULL)
          newdat$prediction_group <- stats::plogis(lp)
        } else if (model_type == "lm") {
          newdat$prediction_group <- predict(model, newdata = newdat, re.form = NULL)
        }
      } else if (inherits(model, "glmmTMB")) {
        lp <- predict(model, newdata = newdat, type = "link", re.form = NULL)
        if (model_type %in% c("nb", "poisson", "cox")) {
          newdat$prediction_group <- exp(lp)
        } else if (model_type == "logistic") {
          newdat$prediction_group <- stats::plogis(lp)
        } else {
          newdat$prediction_group <- lp
        }
      } else if (inherits(model, "coxme")) {
        lp <- predict(model, newdata = newdat, type = "lp")
        newdat$prediction_group <- exp(lp)
      } else {
        # fallback: no random effect
        if (model_type == "cox") {
          lp <- predict(model, newdata = newdat, type = "lp")
          newdat$prediction_group <- exp(lp)
        } else if (model_type %in% c("nb", "poisson", "logistic")) {
          lp <- predict(model, newdata = newdat, type = "link")
          if (model_type %in% c("nb", "poisson")) {
            newdat$prediction_group <- exp(lp)
          } else {
            newdat$prediction_group <- stats::plogis(lp)
          }
        } else if (model_type == "lm") {
          newdat$prediction_group <- predict(model, newdata = newdat)
        }
      }

      newdat[[random_intercept_var]] <- factor(newdat[[random_intercept_var]],
                                               levels = group_levels)
      newdat
    })

    group_predictions <- do.call(rbind, group_pred_list)

    x_lab_tmp <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x
    if (model_type %in% c("nb", "poisson")) {
      y_lab_tmp <- "Predicted Rate"
    } else if (model_type == "logistic") {
      y_lab_tmp <- "Predicted Probability"
    } else if (model_type == "lm") {
      y_lab_tmp <- "Predicted Mean"
    } else {
      y_lab_tmp <- "Predicted Hazard Ratio"
    }

    group_plot <- ggplot2::ggplot(
      group_predictions,
      ggplot2::aes_string(
        x = variable_x,
        y = "prediction_group",
        group = random_intercept_var,
        color = random_intercept_var
      )
    ) +
      ggplot2::geom_line(alpha = 0.4, linewidth = 0.6) +
      ggplot2::xlab(x_lab_tmp) +
      ggplot2::ylab(y_lab_tmp) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        legend.position  = "none",
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  #------------------------------
  # 9. Optional numeric derivatives
  #------------------------------
  derivative_data  <- NULL
  derivative_plot  <- NULL

  if (calculate_derivatives) {
    if (is.null(derivative_points)) {
      derivative_points <- seq(
        from = stats::quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
        to   = stats::quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
        length.out = 100
      )
    }

    if (is.null(subgroups)) {
      deriv_data <- data.frame(x = derivative_points)
      colnames(deriv_data) <- variable_x
    } else {
      subgroup_levels <- if (is.factor(Data_Subset[[subgroups]])) {
        levels(Data_Subset[[subgroups]])
      } else {
        sort(unique(Data_Subset[[subgroups]]))
      }
      deriv_data <- expand.grid(derivative_points, subgroup_levels)
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

    # covariates in deriv_data
    if (!is.null(expanded_covariables)) {
      covariate_vars_for_prediction <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
      covariate_vars_for_prediction <- trimws(covariate_vars_for_prediction)
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

    # numeric derivative: central difference on response scale
    epsilon <- 1e-5 * stats::sd(deriv_data[[variable_x]], na.rm = TRUE)
    data_plus  <- deriv_data
    data_minus <- deriv_data
    data_plus[[variable_x]]  <- deriv_data[[variable_x]] + epsilon
    data_minus[[variable_x]] <- deriv_data[[variable_x]] - epsilon

    predict_resp <- function(model, newdata) {
      if (inherits(model, "merMod")) {
        if (model_type == "cox") {
          lp <- predict(model, newdata = newdata, type = "lp", re.form = NA)
          exp(lp)
        } else if (model_type %in% c("nb","poisson")) {
          lp <- predict(model, newdata = newdata, type = "link", re.form = NA)
          exp(lp)
        } else if (model_type == "logistic") {
          lp <- predict(model, newdata = newdata, type = "link", re.form = NA)
          stats::plogis(lp)
        } else if (model_type == "lm") {
          predict(model, newdata = newdata, re.form = NA)
        }
      } else if (inherits(model, "glmmTMB")) {
        lp <- predict(model, newdata = newdata, type = "link", re.form = NA)
        if (model_type %in% c("nb","poisson","cox")) {
          exp(lp)
        } else if (model_type == "logistic") {
          stats::plogis(lp)
        } else lp
      } else if (inherits(model,"coxme")) {
        lp <- predict(model, newdata = newdata, type = "lp")
        exp(lp)
      } else {
        if (model_type == "cox") {
          lp <- predict(model, newdata = newdata, type = "lp")
          exp(lp)
        } else if (model_type %in% c("nb","poisson")) {
          lp <- predict(model, newdata = newdata, type = "link")
          exp(lp)
        } else if (model_type == "logistic") {
          lp <- predict(model, newdata = newdata, type = "link")
          stats::plogis(lp)
        } else if (model_type == "lm") {
          predict(model, newdata = newdata)
        }
      }
    }

    y_plus  <- predict_resp(model, data_plus)
    y_minus <- predict_resp(model, data_minus)
    deriv   <- (y_plus - y_minus) / (2 * epsilon)
    se_deriv <- abs(deriv) * 0.2  # simple conservative approx
    lower_d <- deriv - 1.96 * se_deriv
    upper_d <- deriv + 1.96 * se_deriv

    derivative_data <- deriv_data
    derivative_data$derivative          <- deriv
    derivative_data$se_derivative       <- se_deriv
    derivative_data$lower_ci_derivative <- lower_d
    derivative_data$upper_ci_derivative <- upper_d
    derivative_data$x_point             <- derivative_data[[variable_x]]

    # find threshold if desired (first point where lower_ci > 0)
    derivative_data$significant_positive <- derivative_data$lower_ci_derivative > 0
    if (is.null(subgroups)) {
      sorted <- derivative_data[order(derivative_data$x_point), ]
      idx <- which(sorted$significant_positive)
      threshold <- if (length(idx) > 0) sorted$x_point[min(idx)] else NA
    } else {
      threshold_values <- tapply(derivative_data$x_point,
                                 derivative_data[[subgroups]],
                                 function(x) NA_real_)
      for (lvl in unique(derivative_data[[subgroups]])) {
        tmp <- derivative_data[derivative_data[[subgroups]] == lvl, ]
        tmp <- tmp[order(tmp$x_point), ]
        idx <- which(tmp$significant_positive)
        if (length(idx) > 0) {
          threshold_values[as.character(lvl)] <- tmp$x_point[min(idx)]
        }
      }
    }

    # small derivative plot (you can extend if desired)
    derivative_plot <- ggplot2::ggplot(
      derivative_data,
      ggplot2::aes_string(x = "x_point", y = "derivative")
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::theme_minimal(base_size = 10)
  }

  #------------------------------
  # 10. Main plot (population)
  #------------------------------
  if (is.null(plot_options)) plot_options <- list()

  x_lab <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x
  y_lab <- if (!is.null(plot_options$y_lab)) {
    plot_options$y_lab
  } else {
    if (model_type %in% c("nb","poisson")) "Predicted Rate"
    else if (model_type == "logistic") "Predicted Probability"
    else if (model_type == "lm") "Predicted Mean"
    else "Predicted Hazard Ratio"
  }

  if (!is.null(subgroups)) {
    plot <- ggplot2::ggplot(
      pred_data,
      ggplot2::aes_string(
        x = variable_x,
        y = "prediction",
        color = subgroups,
        fill  = subgroups,
        linetype = subgroups
      )
    ) +
      ggplot2::geom_line(linewidth = if (!is.null(plot_options$line_size)) plot_options$line_size else 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                           alpha = if (!is.null(plot_options$ribbon_alpha)) plot_options$ribbon_alpha else 0.3)
  } else {
    line_color  <- if (!is.null(plot_options$colors)) plot_options$colors[1] else "black"
    fill_colors <- if (!is.null(plot_options$fill_colors)) plot_options$fill_colors[1] else "grey70"
    plot <- ggplot2::ggplot(
      pred_data,
      ggplot2::aes_string(x = variable_x, y = "prediction")
    ) +
      ggplot2::geom_line(
        linewidth = if (!is.null(plot_options$line_size)) plot_options$line_size else 1,
        color = line_color
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
        fill  = fill_colors,
        alpha = if (!is.null(plot_options$ribbon_alpha)) plot_options$ribbon_alpha else 0.3
      )
  }

  plot <- plot +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::theme_minimal(base_size = 10)

  #------------------------------
  # 11. Return
  #------------------------------
  list(
    predictions             = pred_data,
    model                   = model,
    plot                    = plot,
    plot_data               = pred_data,
    prediction_range_values = c(
      stats::quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
      stats::quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE)
    ),
    derivatives             = derivative_data,
    derivative_method       = if (calculate_derivatives) derivative_method else NULL,
    threshold               = if (calculate_derivatives && !is.null(subgroups)) {
      threshold_values
    } else if (calculate_derivatives) {
      threshold
    } else {
      NULL
    },
    derivative_plot         = derivative_plot,
    group_predictions       = group_predictions,
    group_plot              = group_plot
  )
}
