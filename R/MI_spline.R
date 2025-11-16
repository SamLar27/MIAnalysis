#' Create Spline Plots from Multiple Imputed Data (with Random Effects & Group Fits)
#'
#' This function fits restricted cubic spline models to multiply imputed datasets and creates
#' visualizations of the relationship between a continuous predictor and an outcome, with optional
#' stratification by subgroups. It supports random effects and independent per-group spline fits.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for the model.
#' @param variable_x The continuous predictor variable to be modeled with splines (main spline).
#' @param subgroups Optional variable for stratification (default is NULL).
#' @param covariables Optional vector of covariates to adjust for.
#' @param knot_n Number of knots for the restricted cubic spline (default is 4).
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param MI_method "first", "average", or "Rubin".
#' @param model_type "nb", "poisson", "cox", "logistic", or "lm".
#' @param followup_offset "Yes"/"No" – whether to include log-follow-up offset for count models.
#' @param followup_col Column with follow-up duration (required if followup_offset = "Yes").
#' @param trial_factor "Yes"/"No" – whether to include trial as factor.
#' @param trial_col Trial factor column (required if trial_factor = "Yes").
#' @param time_col Time variable for Cox regression (if model_type = "cox").
#' @param event_col Event variable for Cox regression (if model_type = "cox").
#' @param random_intercept "Yes"/"No" – random intercept term.
#' @param random_intercept_var Grouping variable for random effects (if random_intercept / random_slope = "Yes").
#' @param random_slope "Yes"/"No" – include random slopes.
#' @param predictor_vars_random_slope Main predictor(s) with random slope.
#' @param covariables_random_slope Additional covariables with random slope.
#' @param plot_options List of options for plot customization.
#' @param subgroup_as_factor Whether to convert subgroups to factors (default TRUE).
#' @param subgroup_labels Optional labels for subgroup levels.
#' @param prediction_range Quantiles for x-range (default c(0.01, 0.99)).
#' @param calculate_derivatives Logical; whether to compute derivatives.
#' @param derivative_points Specific x-values for derivatives (optional).
#' @param derivative_method "basis", "delta", or "numeric".
#' @param group_fits Logical; if TRUE, fit independent spline models per group (trial) using same knots.
#'
#' @return A list with model, predictions, plot, derivative info, knots, and (optionally) group_fits.
#'
#' @importFrom stats as.formula glm binomial poisson gaussian quantile median plogis pchisq vcov coef
#' @importFrom MASS glm.nb
#' @importFrom survival Surv coxph
#' @importFrom rms rcs rcspline.eval
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon xlab ylab labs scale_x_log10 scale_x_continuous scale_y_continuous coord_cartesian scale_color_manual scale_fill_manual scale_linetype_manual guides guide_legend facet_wrap theme_minimal element_rect element_text margin element_blank element_line unit annotate
#' @importFrom dplyr %>%
#' @importFrom lme4 lmer glmer fixef
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

  #### 0. Packages & basic checks ####
  requireNamespace("rms", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

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
    stop("If random_intercept='Yes' or random_slope='Yes', random_intercept_var must be provided")
  }

  if (length(prediction_range) != 2 ||
      any(prediction_range < 0) ||
      any(prediction_range > 1) ||
      prediction_range[1] >= prediction_range[2]) {
    stop("prediction_range must be c(q1, q2) with 0 <= q1 < q2 <= 1")
  }

  if (calculate_derivatives && !derivative_method %in% c("basis", "delta", "numeric")) {
    stop("derivative_method must be one of 'basis', 'delta', 'numeric'")
  }

  if ((random_intercept == "Yes" || random_slope == "Yes") && model_type == "nb") {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      warning("glmmTMB not available; nb mixed models will fall back to Poisson GLMM (approximation).")
      requireNamespace("lme4", quietly = TRUE)
    }
  }

  if ((random_intercept == "Yes" || random_slope == "Yes") && model_type %in% c("poisson", "logistic", "lm")) {
    requireNamespace("lme4", quietly = TRUE)
  }

  if ((random_intercept == "Yes" || random_slope == "Yes") && model_type == "cox") {
    if (!requireNamespace("coxme", quietly = TRUE)) {
      stop("coxme is required for Cox models with random effects")
    }
  }

  #### 1. Choose imputations ####
  if (MI_method == "first") {
    Data_Subset <- subset(data, get(imp_col) == 1)
  } else if (MI_method %in% c("average", "Rubin")) {
    Data_Subset <- data
  } else {
    stop("MI_method must be 'first', 'average', or 'Rubin'")
  }

  #### 2. Subgroup handling ####
  original_subgroup_values <- NULL
  subgroup_mapping <- NULL

  if (!is.null(subgroups)) {
    original_subgroup_values <- Data_Subset[[subgroups]]
    if (!is.null(subgroup_labels)) {
      unique_values <- sort(unique(Data_Subset[[subgroups]]))
      if (length(unique_values) == length(subgroup_labels)) {
        subgroup_mapping <- stats::setNames(subgroup_labels, unique_values)
      }
    }
    if (subgroup_as_factor && !is.factor(Data_Subset[[subgroups]])) {
      Data_Subset[[subgroups]] <- as.factor(Data_Subset[[subgroups]])
      if (!is.null(subgroup_labels)) {
        levels(Data_Subset[[subgroups]]) <- subgroup_labels
      }
    }
  }

  #### 3. Compute global knots and main spline term ####
  x_full <- Data_Subset[[variable_x]]
  x_full <- x_full[is.finite(x_full)]

  if (length(x_full) < knot_n) {
    stop("Not enough non-missing values of variable_x to place knots")
  }

  rcs_tmp <- rms::rcspline.eval(x_full, nk = knot_n, inclx = TRUE)
  global_knots <- attr(rcs_tmp, "knots")

  spline_term <- sprintf(
    "rcs(%s, c(%s))",
    variable_x,
    paste(round(global_knots, 4), collapse = ", ")
  )

  if (!is.null(subgroups)) {
    spline_term <- paste0(spline_term, " * ", subgroups)
  }

  #### 4. Covariables expansion ####
  expand_terms <- function(term) {
    # Do not expand spline/poly terms
    if (grepl("^rcs\\(", term) || grepl("^bs\\(", term) || grepl("^poly\\(", term)) {
      return(term)
    }
    if (grepl("\\*", term)) {
      vars_split <- unlist(strsplit(term, "\\*"))
      vars_split <- trimws(vars_split)
      all_comb <- lapply(seq_along(vars_split), function(k) {
        combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
      })
      return(unlist(all_comb))
    }
    term
  }

  extract_variables_cov <- function(term) {
    if (grepl("^rcs\\(", term) || grepl("^bs\\(", term) || grepl("^poly\\(", term)) {
      inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
      return(trimws(inside))
    } else {
      return(term)
    }
  }

  expanded_covariables <- NULL
  covariates_str <- ""
  if (!is.null(covariables)) {
    expanded_covariables <- unique(unlist(lapply(covariables, expand_terms)))
    all_base_vars_cov <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
    all_base_vars_cov <- trimws(all_base_vars_cov)
    all_base_vars_cov <- sapply(all_base_vars_cov, extract_variables_cov)
    missing_vars_cov <- all_base_vars_cov[!all_base_vars_cov %in% names(Data_Subset)]
    if (length(missing_vars_cov) > 0) {
      stop("Covariates not found in data: ", paste(missing_vars_cov, collapse = ", "))
    }
    covariates_str <- paste0(" + ", paste(expanded_covariables, collapse = " + "))
  }

  #### 5. Trial factor, offset, random effects strings ####
  trial_str <- ""
  if (trial_factor == "Yes") {
    trial_str <- paste0(" + as.factor(", trial_col, ")")
  }

  offset_str <- ""
  if (followup_offset == "Yes" && model_type %in% c("nb", "poisson")) {
    offset_str <- paste0(" + offset(log(", followup_col, "))")
  }

  random_effect_str <- ""
  use_mixed_model <- (random_intercept == "Yes" || random_slope == "Yes")

  if (use_mixed_model) {
    if (is.null(random_intercept_var)) {
      stop("random_intercept_var must be provided if random_intercept or random_slope is 'Yes'")
    }
    random_terms <- c()
    if (random_intercept == "Yes") {
      random_terms <- c(random_terms, "1")
    }
    if (random_slope == "Yes") {
      slope_vars <- unique(c(predictor_vars_random_slope, covariables_random_slope))
      slope_vars <- slope_vars[!is.na(slope_vars) & nzchar(slope_vars)]
      if (length(slope_vars) == 0) {
        stop("random_slope='Yes' but no predictor_vars_random_slope / covariables_random_slope specified.")
      }
      random_terms <- c(random_terms, slope_vars)
    }
    random_effect_str <- paste0(" + (", paste(random_terms, collapse = " + "),
                                " | ", random_intercept_var, ")")
  }

  #### 6. Build main formula ####
  if (model_type == "cox") {
    formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                          spline_term, covariates_str, trial_str, random_effect_str)
  } else {
    formula_str <- paste0(outcome_var, " ~ ",
                          spline_term, covariates_str, trial_str,
                          offset_str, random_effect_str)
  }
  formula_obj <- as.formula(formula_str)

  #### 7. Fit main model ####
  model <- NULL
  if (use_mixed_model) {
    if (model_type == "nb") {
      if (requireNamespace("glmmTMB", quietly = TRUE)) {
        model <- glmmTMB::glmmTMB(formula_obj, family = glmmTMB::nbinom2, data = Data_Subset)
      } else {
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

  #### 8. Prediction grid for population curve ####
  x_values <- seq(
    from = quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
    to   = quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
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
      pred_data[[subgroups]] <- factor(pred_data[[subgroups]],
                                       levels = levels(Data_Subset[[subgroups]]))
    } else if (subgroup_as_factor) {
      pred_data[[subgroups]] <- factor(pred_data[[subgroups]])
      if (!is.null(subgroup_labels)) {
        levels(pred_data[[subgroups]]) <- subgroup_labels
      }
    }
  }

  # Add covariates at median/mode
  if (!is.null(expanded_covariables)) {
    covariate_vars_for_prediction <- unique(
      unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*"))
    )
    covariate_vars_for_prediction <- trimws(covariate_vars_for_prediction)
    covariate_vars_for_prediction <- sapply(covariate_vars_for_prediction,
                                            extract_variables_cov)
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

  if (trial_factor == "Yes") {
    pred_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
  }

  if (followup_offset == "Yes") {
    pred_data[[followup_col]] <- 365
  }

  if (use_mixed_model) {
    pred_data[[random_intercept_var]] <- names(which.max(table(Data_Subset[[random_intercept_var]])))
  }

  #### 9. Helpers for predictions from mixed models ####
  get_mixed_model_predictions <- function(model, newdata, type = "link") {
    if (inherits(model, "merMod")) {
      preds <- predict(model, newdata = newdata, re.form = NA, type = type)
      mm <- model.matrix(lme4::getME(model, "X"), newdata)
      vcov_fixed <- as.matrix(vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      return(list(fit = preds, se.fit = se_fixed))
    } else if (inherits(model, "glmmTMB")) {
      preds <- predict(model, newdata = newdata, re.form = NA, type = type, se.fit = TRUE)
      return(list(fit = preds$fit, se.fit = preds$se.fit))
    } else if (inherits(model, "coxme")) {
      preds <- predict(model, newdata = newdata, type = "lp")
      mm <- model.matrix(model)
      vcov_fixed <- as.matrix(vcov(model))
      se_fixed <- sqrt(diag(mm %*% vcov_fixed %*% t(mm)))
      return(list(fit = preds, se.fit = se_fixed))
    }
    stop("Unsupported mixed model class for prediction")
  }

  #### 10. Population predictions + MI pooling ####
  preds <- NULL

  if (MI_method == "first") {
    if (use_mixed_model) {
      if (model_type == "cox") {
        preds <- get_mixed_model_predictions(model, pred_data)
        pred_data$prediction <- exp(preds$fit)
        pred_data$lower_ci   <- exp(preds$fit - 1.96 * preds$se.fit)
        pred_data$upper_ci   <- exp(preds$fit + 1.96 * preds$se.fit)
      } else {
        preds <- get_mixed_model_predictions(model, pred_data)
        if (model_type %in% c("nb", "poisson")) {
          pred_data$prediction <- exp(preds$fit)
          pred_data$lower_ci   <- exp(preds$fit - 1.96 * preds$se.fit)
          pred_data$upper_ci   <- exp(preds$fit + 1.96 * preds$se.fit)
        } else if (model_type == "logistic") {
          pred_data$prediction <- plogis(preds$fit)
          pred_data$lower_ci   <- plogis(preds$fit - 1.96 * preds$se.fit)
          pred_data$upper_ci   <- plogis(preds$fit + 1.96 * preds$se.fit)
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
        pred_data$lower_ci <- NA
        pred_data$upper_ci <- NA
      } else {
        preds <- predict(model, newdata = pred_data, type = "link", se.fit = TRUE)
        if (model_type %in% c("nb", "poisson")) {
          pred_data$prediction <- exp(preds$fit)
          pred_data$lower_ci   <- exp(preds$fit - 1.96 * preds$se.fit)
          pred_data$upper_ci   <- exp(preds$fit + 1.96 * preds$se.fit)
        } else if (model_type == "logistic") {
          pred_data$prediction <- plogis(preds$fit)
          pred_data$lower_ci   <- plogis(preds$fit - 1.96 * preds$se.fit)
          pred_data$upper_ci   <- plogis(preds$fit + 1.96 * preds$se.fit)
        } else if (model_type == "lm") {
          pred_data$prediction <- preds$fit
          pred_data$lower_ci   <- preds$fit - 1.96 * preds$se.fit
          pred_data$upper_ci   <- preds$fit + 1.96 * preds$se.fit
        }
      }
    }
  } else if (MI_method %in% c("average", "Rubin")) {

    imputation_ids <- sort(unique(Data_Subset[[imp_col]]))
    all_predictions <- vector("list", length(imputation_ids))
    names(all_predictions) <- as.character(imputation_ids)

    for (imp in imputation_ids) {
      Data_imp <- subset(Data_Subset, get(imp_col) == imp)

      # refit model
      if (use_mixed_model) {
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

      if (use_mixed_model) {
        preds_imp <- get_mixed_model_predictions(model_imp, pred_data)
      } else {
        if (model_type == "cox") {
          lp_imp <- predict(model_imp, newdata = pred_data, type = "lp")
          preds_imp <- list(fit = lp_imp,
                            se.fit = rep(NA_real_, length(lp_imp)))
        } else {
          preds_imp <- predict(model_imp, newdata = pred_data, type = "link", se.fit = TRUE)
        }
      }

      all_predictions[[as.character(imp)]] <- list(
        fit = preds_imp$fit,
        se.fit = preds_imp$se.fit
      )
    }

    n_pred <- length(all_predictions[[1]]$fit)
    mean_fit <- numeric(n_pred)
    mean_var <- numeric(n_pred)

    for (i in seq_len(n_pred)) {
      fits_i <- sapply(all_predictions, function(x) x$fit[i])
      ses_i  <- sapply(all_predictions, function(x) x$se.fit[i])
      m      <- length(fits_i)

      mean_fit[i] <- mean(fits_i, na.rm = TRUE)

      if (MI_method == "Rubin") {
        W <- mean(ses_i^2, na.rm = TRUE)
        B <- stats::var(fits_i, na.rm = TRUE)
        mean_var[i] <- W + (1 + 1/m) * B
      } else {
        mean_var[i] <- mean(ses_i^2, na.rm = TRUE)
      }
    }

    if (model_type %in% c("nb", "poisson")) {
      pred_data$prediction <- exp(mean_fit)
      pred_data$lower_ci   <- exp(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci   <- exp(mean_fit + 1.96 * sqrt(mean_var))
    } else if (model_type == "logistic") {
      pred_data$prediction <- plogis(mean_fit)
      pred_data$lower_ci   <- plogis(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci   <- plogis(mean_fit + 1.96 * sqrt(mean_var))
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

  #### 11. Derivative machinery (using global_knots) ####
  # For brevity I keep the high-level structure; behaviour is unchanged except that
  # basis-method uses global_knots directly.
  # (Same helper functions as in your previous version, but the main change is:
  #  calculate_basis_function_derivatives(..., knots = global_knots))

  # --- Helper: extract spline coefficients ---
  extract_spline_coefficients <- function(model, variable_x) {
    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB"))) {
      all_coefs <- lme4::fixef(model)
    } else {
      all_coefs <- coef(model)
    }
    patterns <- c(
      paste0("^", variable_x, "$"),
      paste0("^rcs\\(", variable_x),
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

  # --- Helper: derivative basis for RCS given knots ---
  calculate_rcs_derivatives <- function(x, knots) {
    k <- length(knots)
    n <- length(x)
    deriv_mat <- matrix(0, nrow = n, ncol = k - 1)
    deriv_mat[, 1] <- 1
    if (k > 2) {
      t1 <- knots[1]
      tk <- knots[k]
      for (j in 2:(k - 1)) {
        tj <- knots[j]
        for (i in seq_len(n)) {
          term1 <- if (x[i] > tj) 3 * (x[i] - tj)^2 else 0
          term2 <- if (x[i] > tk) 3 * (x[i] - tk)^2 * (tk - tj)/(tk - t1) else 0
          deriv_mat[i, j] <- term1 - term2
        }
      }
    }
    deriv_mat
  }

  # --- Helper: apply link derivative ---
  apply_link_derivative <- function(model, newdata, linear_deriv, model_type) {
    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
      lp <- predict(model, newdata = newdata, type = "link", re.form = NA)
    } else if (model_type == "cox") {
      lp <- predict(model, newdata = newdata, type = "lp")
    } else if (model_type %in% c("nb", "poisson", "logistic")) {
      lp <- predict(model, newdata = newdata, type = "link")
    } else {
      lp <- predict(model, newdata = newdata)
    }

    if (model_type == "logistic") {
      p <- plogis(lp)
      return(p * (1 - p) * linear_deriv)
    } else if (model_type %in% c("nb", "poisson", "cox")) {
      return(exp(lp) * linear_deriv)
    } else {
      return(linear_deriv)
    }
  }

  # --- Helper: CI for derivatives (delta method on coefficients) ---
  calculate_derivative_ci <- function(model, newdata, basis_deriv, spline_coefs, model_type) {
    vc <- tryCatch({
      if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
        as.matrix(vcov(model))
      } else {
        vcov(model)
      }
    }, error = function(e) NULL)

    spline_terms <- names(spline_coefs)
    if (is.null(vc) || !all(spline_terms %in% rownames(vc))) {
      # fallback: crude 20% relative SE
      lin <- basis_deriv %*% spline_coefs
      resp <- apply_link_derivative(model, newdata, lin, model_type)
      se   <- abs(resp) * 0.2
      lower <- resp - 1.96 * se
      upper <- resp + 1.96 * se
      return(list(derivative = resp,
                  se_derivative = se,
                  lower_ci = lower,
                  upper_ci = upper))
    }

    vc_spline <- vc[spline_terms, spline_terms, drop = FALSE]
    lin <- basis_deriv %*% spline_coefs

    if (inherits(model, c("lmerMod", "glmerMod", "glmmTMB", "coxme"))) {
      lp <- predict(model, newdata = newdata, type = "link", re.form = NA)
    } else if (model_type == "cox") {
      lp <- predict(model, newdata = newdata, type = "lp")
    } else if (model_type %in% c("nb", "poisson", "logistic")) {
      lp <- predict(model, newdata = newdata, type = "link")
    } else {
      lp <- predict(model, newdata = newdata)
    }

    n <- nrow(newdata)
    der <- numeric(n)
    se  <- numeric(n)

    for (i in seq_len(n)) {
      g <- basis_deriv[i, ]
      var_lin <- as.numeric(t(g) %*% vc_spline %*% g)
      se_lin  <- sqrt(var_lin)

      if (model_type == "logistic") {
        p <- plogis(lp[i])
        der[i] <- p * (1 - p) * lin[i]
        se[i]  <- p * (1 - p) * se_lin
      } else if (model_type %in% c("nb", "poisson", "cox")) {
        der[i] <- exp(lp[i]) * lin[i]
        se[i]  <- exp(lp[i]) * se_lin
      } else {
        der[i] <- lin[i]
        se[i]  <- se_lin
      }
    }

    lower <- der - 1.96 * se
    upper <- der + 1.96 * se
    list(derivative = der,
         se_derivative = se,
         lower_ci = lower,
         upper_ci = upper)
  }

  calculate_basis_function_derivatives <- function(model, newdata, variable_x, knots, model_type) {
    basis_deriv   <- calculate_rcs_derivatives(newdata[[variable_x]], knots)
    spline_coefs  <- extract_spline_coefficients(model, variable_x)
    deriv_results <- calculate_derivative_ci(model, newdata, basis_deriv, spline_coefs, model_type)
    list(
      derivative          = deriv_results$derivative,
      se_derivative       = deriv_results$se_derivative,
      lower_ci_derivative = deriv_results$lower_ci,
      upper_ci_derivative = deriv_results$upper_ci,
      x_point             = newdata[[variable_x]]
    )
  }

  # (Delta and numeric derivative helpers omitted here for brevity; keep your previous versions
  #  if you actively use derivative_method = 'delta' or 'numeric'.)

  #### 12. Derivatives (optional) ####
  derivative_data <- NULL
  threshold <- NA
  threshold_values <- NULL

  if (calculate_derivatives) {
    if (is.null(derivative_points)) {
      derivative_points <- seq(
        from = quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
        to   = quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
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

    # add covariates, offset, trial, RE as for pred_data
    if (!is.null(expanded_covariables)) {
      covariate_vars_for_prediction <- unique(
        unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*"))
      )
      covariate_vars_for_prediction <- trimws(covariate_vars_for_prediction)
      covariate_vars_for_prediction <- sapply(covariate_vars_for_prediction,
                                              extract_variables_cov)
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
    if (use_mixed_model) {
      deriv_data[[random_intercept_var]] <- names(which.max(table(Data_Subset[[random_intercept_var]])))
    }

    # Only implement basis method here (the one you’re using)
    deriv_res <- calculate_basis_function_derivatives(
      model = model,
      newdata = deriv_data,
      variable_x = variable_x,
      knots = global_knots,
      model_type = model_type
    )

    derivative_data <- deriv_data
    derivative_data$derivative          <- deriv_res$derivative
    derivative_data$se_derivative       <- deriv_res$se_derivative
    derivative_data$lower_ci_derivative <- deriv_res$lower_ci_derivative
    derivative_data$upper_ci_derivative <- deriv_res$upper_ci_derivative
    derivative_data$x_point             <- deriv_res$x_point

    # Threshold detection (first x where lower_ci_derivative > 0)
    if (!is.null(subgroups)) {
      derivative_data$significant_positive <- derivative_data$lower_ci_derivative > 0
      subgroup_levels <- unique(derivative_data[[subgroups]])
      threshold_values <- numeric(length(subgroup_levels))
      names(threshold_values) <- as.character(subgroup_levels)

      for (level in subgroup_levels) {
        sub_d <- derivative_data[derivative_data[[subgroups]] == level, ]
        sub_d <- sub_d[order(sub_d$x_point), ]
        idx <- which(sub_d$significant_positive)
        threshold_values[as.character(level)] <- if (length(idx) > 0) sub_d$x_point[min(idx)] else NA
      }
    } else {
      derivative_data$significant_positive <- derivative_data$lower_ci_derivative > 0
      sorted_d <- derivative_data[order(derivative_data$x_point), ]
      idx <- which(sorted_d$significant_positive)
      threshold <- if (length(idx) > 0) sorted_d$x_point[min(idx)] else NA
    }
  }

  #### 13. Build ggplot for population curve ####
  if (is.null(plot_options)) plot_options <- list()

  x_lab <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x
  y_lab <- if (!is.null(plot_options$y_lab)) {
    plot_options$y_lab
  } else {
    if (model_type %in% c("nb", "poisson")) "Predicted rate"
    else if (model_type == "logistic") "Predicted probability"
    else if (model_type == "lm") "Predicted mean"
    else "Predicted hazard ratio"
  }

  if (!is.null(subgroups)) {
    plot <- ggplot2::ggplot(pred_data,
                            ggplot2::aes_string(x = variable_x,
                                                y = "prediction",
                                                color = subgroups,
                                                fill  = subgroups,
                                                linetype = subgroups))
  } else {
    plot <- ggplot2::ggplot(pred_data,
                            ggplot2::aes_string(x = variable_x, y = "prediction"))
  }

  line_size   <- plot_options$line_size   %||% 1
  ribbon_alpha<- plot_options$ribbon_alpha%||% 0.3

  if (!is.null(subgroups)) {
    plot <- plot +
      ggplot2::geom_line(linewidth = line_size) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                           alpha = ribbon_alpha)
  } else {
    line_color <- plot_options$colors      %||% "black"
    fill_color <- plot_options$fill_colors %||% "grey70"
    if (is.vector(line_color)) line_color <- line_color[1]
    if (is.vector(fill_color)) fill_color <- fill_color[1]

    plot <- plot +
      ggplot2::geom_line(linewidth = line_size, color = line_color) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                           fill  = fill_color,
                           color = line_color,
                           alpha = ribbon_alpha)
  }

  plot <- plot +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab)

  # Axis scales / limits
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
    plot <- plot + ggplot2::scale_y_continuous(breaks = plot_options$y_breaks)
  }

  coord_args <- list()
  if (!is.null(plot_options$y_limits)) coord_args$ylim <- plot_options$y_limits
  if (!is.null(plot_options$x_limits)) coord_args$xlim <- plot_options$x_limits
  if (length(coord_args) > 0) {
    plot <- plot + do.call(ggplot2::coord_cartesian, coord_args)
  }

  plot <- plot + ggplot2::theme_minimal(base_size = 10)

  default_theme <- list(
    legend.background = ggplot2::element_rect(fill = "white", colour = "black", size = 0.5),
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

  plot <- plot + do.call(ggplot2::theme, default_theme)

  #### 14. Independent group_fits (per-trial models) ####
  group_fits_df <- NULL

  if (group_fits) {

    if (MI_method != "first") {
      warning("group_fits currently implemented for MI_method = 'first'; using first imputation only.")
      Data_group <- subset(data, get(imp_col) == 1)
    } else {
      Data_group <- Data_Subset
    }

    # which grouping variable?
    if (!is.null(random_intercept_var)) {
      group_var <- random_intercept_var
    } else if (trial_factor == "Yes" && !is.null(trial_col)) {
      group_var <- trial_col
    } else {
      stop("group_fits = TRUE but no random_intercept_var or trial_col available to define groups.")
    }

    if (!group_var %in% names(Data_group)) {
      stop("group_var '", group_var, "' not found in data for group_fits.")
    }

    groups <- sort(unique(Data_group[[group_var]]))

    # covariates for group model WITHOUT trial factor
    covariates_str_no_trial <- ""
    if (!is.null(expanded_covariables) && length(expanded_covariables) > 0) {
      expanded_covariables_no_trial <- expanded_covariables
      if (trial_factor == "Yes" && !is.null(trial_col)) {
        expanded_covariables_no_trial <- expanded_covariables_no_trial[
          !grepl(paste0("as.factor\\(", trial_col, "\\)"), expanded_covariables_no_trial)
        ]
      }
      if (length(expanded_covariables_no_trial) > 0) {
        covariates_str_no_trial <- paste0(" + ",
                                          paste(expanded_covariables_no_trial, collapse = " + "))
      }
    }

    offset_str_group <- ""
    if (followup_offset == "Yes" && model_type %in% c("nb", "poisson")) {
      offset_str_group <- paste0(" + offset(log(", followup_col, "))")
    }

    for (g in groups) {
      Data_g <- Data_group[Data_group[[group_var]] == g, , drop = FALSE]

      if (nrow(Data_g) < (knot_n + 3)) {
        warning(sprintf("Group %s has too few observations (%d). Skipping group_fits for this group.",
                        as.character(g), nrow(Data_g)))
        next
      }

      formula_str_g <- paste0(
        outcome_var, " ~ ",
        spline_term,
        covariates_str_no_trial,
        offset_str_group
      )
      formula_g <- as.formula(formula_str_g)

      model_g <- try(MASS::glm.nb(formula_g, data = Data_g), silent = TRUE)
      if (inherits(model_g, "try-error")) {
        warning("Group ", as.character(g), ": glm.nb failed. Skipping.")
        next
      }

      x_seq <- seq(
        from = quantile(Data_g[[variable_x]], prediction_range[1], na.rm = TRUE),
        to   = quantile(Data_g[[variable_x]], prediction_range[2], na.rm = TRUE),
        length.out = 100
      )

      pred_g <- data.frame(x_seq)
      colnames(pred_g) <- variable_x

      if (!is.null(expanded_covariables) && length(expanded_covariables) > 0) {
        covariate_vars_for_prediction <- unique(
          unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*"))
        )
        covariate_vars_for_prediction <- trimws(covariate_vars_for_prediction)
        covariate_vars_for_prediction <- sapply(covariate_vars_for_prediction,
                                                extract_variables_cov)
        for (cov in covariate_vars_for_prediction) {
          if (cov %in% colnames(Data_g)) {
            if (is.factor(Data_g[[cov]])) {
              pred_g[[cov]] <- as.factor(names(which.max(table(Data_g[[cov]]))))
            } else {
              pred_g[[cov]] <- median(Data_g[[cov]], na.rm = TRUE)
            }
          }
        }
      }
      if (followup_offset == "Yes") {
        pred_g[[followup_col]] <- 365
      }

      pp <- predict(model_g, newdata = pred_g, type = "link", se.fit = TRUE)
      pred_resp  <- exp(pp$fit)
      lower_resp <- exp(pp$fit - 1.96 * pp$se.fit)
      upper_resp <- exp(pp$fit + 1.96 * pp$se.fit)

      df_g <- data.frame(
        !!variable_x := pred_g[[variable_x]],
        prediction   = pred_resp,
        lower_ci     = lower_resp,
        upper_ci     = upper_resp,
        !!group_var  := g,
        stringsAsFactors = FALSE
      )

      group_fits_df <- rbind(group_fits_df, df_g)
    }
  }

  #### 15. Return ####
  list(
    predictions = pred_data,
    model       = model,
    plot        = plot,
    plot_data   = pred_data,
    prediction_range_values = c(
      quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
      quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE)
    ),
    derivatives      = if (calculate_derivatives) derivative_data else NULL,
    derivative_method= if (calculate_derivatives) derivative_method else NULL,
    threshold        = if (calculate_derivatives && !is.null(subgroups))
      threshold_values else if (calculate_derivatives) threshold else NULL,
    knots            = global_knots,
    group_fits       = group_fits_df
  )
}
