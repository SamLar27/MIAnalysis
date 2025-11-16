#' Create Spline Plots from Multiple Imputed Data (with Random Effects & Group Fits)
#'
#' Fits restricted cubic spline models to multiply imputed datasets and creates
#' visualizations of the relationship between a continuous predictor and the outcome,
#' with optional random effects and trial-specific "group_fits".
#'
#' @param data Data frame containing all imputations.
#' @param outcome_var Name of outcome variable (string).
#' @param variable_x Continuous predictor modelled with rcs().
#' @param subgroups Optional subgroup variable for interaction with spline.
#' @param covariables Optional vector of additional covariates (strings). Can contain
#'        interaction terms like "A*B".
#' @param knot_n Number of knots for rcs (default 4).
#' @param imp_col Column name with imputation index (default ".imp").
#' @param MI_method "first", "average", or "Rubin".
#' @param model_type "nb", "poisson", "logistic", "lm", or "cox".
#' @param followup_offset "Yes"/"No" for offset(log(followup_col)) in count models.
#' @param followup_col Name of follow-up time column (required if offset = "Yes").
#' @param trial_factor "Yes"/"No": add trial as fixed factor.
#' @param trial_col Trial factor name (required if trial_factor = "Yes").
#' @param time_col Time variable for Cox models.
#' @param event_col Event variable for Cox models.
#' @param random_intercept "Yes"/"No" for random intercept.
#' @param random_intercept_var Grouping variable for random effects.
#' @param random_slope "Yes"/"No" for random slopes.
#' @param predictor_vars_random_slope Variables with random slopes (character vector).
#' @param covariables_random_slope Additional covariates with random slopes.
#' @param plot_options List of ggplot options.
#' @param subgroup_as_factor Convert subgroups to factor (default TRUE).
#' @param subgroup_labels Optional custom labels for subgroup levels.
#' @param prediction_range Quantiles for x-range (default c(0.01, 0.99)).
#' @param calculate_derivatives Whether to compute derivatives (default FALSE).
#' @param derivative_points Optional explicit x-values for derivatives.
#' @param derivative_method Currently ignored; numeric finite-difference is used.
#' @param group_fits If TRUE, refit independent models per group (e.g. per trial)
#'        using the same global knots; returns spline per group.
#'
#' @return A list with:
#'   - predictions: prediction data for population curve
#'   - model: fitted global model
#'   - plot: ggplot object for population curve
#'   - plot_data: same as predictions
#'   - prediction_range_values: numeric range used for x
#'   - derivatives: derivative data (if requested)
#'   - derivative_plot: reserved (currently NULL)
#'   - threshold: approximate threshold where derivative CI>0 (if requested)
#'   - group_fits: data frame of trial-specific fits (if group_fits = TRUE)
#'
#' @importFrom stats as.formula glm binomial poisson gaussian predict quantile median plogis setNames vcov
#' @importFrom MASS glm.nb
#' @importFrom survival Surv coxph
#' @importFrom rms rcs
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon xlab ylab labs scale_x_log10 scale_x_continuous scale_y_continuous coord_cartesian scale_color_manual scale_fill_manual scale_linetype_manual guides guide_legend facet_wrap theme_minimal element_rect element_text margin element_blank element_line unit annotate
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
                      random_slope  = "No",
                      predictor_vars_random_slope = NULL,
                      covariables_random_slope    = NULL,
                      plot_options = NULL,
                      subgroup_as_factor = TRUE,
                      subgroup_labels = NULL,
                      prediction_range = c(0.01, 0.99),
                      calculate_derivatives = FALSE,
                      derivative_points = NULL,
                      derivative_method = "numeric",
                      group_fits = FALSE) {

  # small helper
  `%||%` <- function(x, y) if (!is.null(x)) x else y

  # ---- Basic checks ----
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
  if ((random_intercept == "Yes" || random_slope == "Yes") &&
      is.null(random_intercept_var)) {
    stop("If random_intercept='Yes' or random_slope='Yes', random_intercept_var must be provided")
  }

  # ---- MI subset ----
  if (MI_method == "first") {
    Data_Subset <- subset(data, get(imp_col) == 1)
  } else if (MI_method %in% c("average", "Rubin")) {
    Data_Subset <- data
  } else {
    stop("MI_method must be 'first', 'average', or 'Rubin'")
  }

  # ---- Subgroups handling ----
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

  # ---- Global spline knots (explicit) ----
  x_all <- Data_Subset[[variable_x]]
  x_all <- x_all[is.finite(x_all)]
  if (length(x_all) == 0) stop("No valid values for variable_x in data")

  probs_knots  <- seq(0, 1, length.out = knot_n)
  global_knots <- as.numeric(stats::quantile(x_all, probs = probs_knots, na.rm = TRUE))

  spline_term <- paste0(
    "rcs(", variable_x, ", c(",
    paste(round(global_knots, 4), collapse = ", "),
    "))"
  )
  if (!is.null(subgroups)) {
    spline_term <- paste0(spline_term, " * ", subgroups)
  }

  # ---- Expand covariables ----
  expand_terms <- function(term) {
    if (grepl("^rcs\\(", term) || grepl("^bs\\(", term) || grepl("^poly\\(", term)) {
      return(term)
    }
    if (grepl("\\*", term)) {
      vars_split <- trimws(unlist(strsplit(term, "\\*")))
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
  }

  covariates_str <- if (!is.null(expanded_covariables) && length(expanded_covariables) > 0) {
    paste("+", paste(expanded_covariables, collapse = " + "))
  } else {
    ""
  }

  # global covariate vars for prediction templates
  global_covariate_vars <- NULL
  if (!is.null(expanded_covariables)) {
    extract_variables_cov <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^bs\\(", term) || grepl("^poly\\(", term)) {
        inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
        return(trimws(inside))
      } else {
        return(term)
      }
    }
    global_covariate_vars <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
    global_covariate_vars <- trimws(global_covariate_vars)
    global_covariate_vars <- sapply(global_covariate_vars, extract_variables_cov)
  }

  # ---- Trial factor & offset ----
  trial_str <- ""
  if (trial_factor == "Yes") {
    trial_str <- paste0(" + as.factor(", trial_col, ")")
  }

  offset_str <- ""
  if (followup_offset == "Yes" && model_type %in% c("nb", "poisson")) {
    offset_str <- paste0(" + offset(log(", followup_col, "))")
  }

  # ---- Random effects syntax ----
  random_effect_str <- ""
  use_mixed <- (random_intercept == "Yes" || random_slope == "Yes")

  slope_terms <- character(0)
  if (!is.null(predictor_vars_random_slope)) {
    slope_terms <- c(slope_terms, predictor_vars_random_slope)
  }
  if (!is.null(covariables_random_slope)) {
    slope_terms <- c(slope_terms, covariables_random_slope)
  }
  slope_terms <- unique(slope_terms)

  if (use_mixed) {
    if (length(slope_terms) == 0 && random_slope == "Yes") {
      stop("random_slope='Yes' but no predictor_vars_random_slope/covariables_random_slope provided.")
    }

    if (random_intercept == "Yes" && random_slope == "Yes" && length(slope_terms) > 0) {
      random_effect_str <- paste0(" + (1 + ", paste(slope_terms, collapse = " + "),
                                  " | ", random_intercept_var, ")")
    } else if (random_intercept == "Yes") {
      random_effect_str <- paste0(" + (1 | ", random_intercept_var, ")")
    } else if (random_slope == "Yes" && length(slope_terms) > 0) {
      random_effect_str <- paste0(" + (0 + ", paste(slope_terms, collapse = " + "),
                                  " | ", random_intercept_var, ")")
    }
  }

  # ---- Build main formula ----
  if (model_type == "cox") {
    formula_str <- paste0(
      "survival::Surv(", time_col, ", ", event_col, ") ~ ",
      spline_term, " ", covariates_str, trial_str, random_effect_str
    )
  } else {
    formula_str <- paste0(
      outcome_var, " ~ ",
      spline_term, " ", covariates_str, trial_str, offset_str, random_effect_str
    )
  }
  formula_obj <- stats::as.formula(formula_str)

  # ---- Fit main model ----
  fit_single_model <- function(dat) {
    if (use_mixed) {
      if (model_type == "nb") {
        if (requireNamespace("glmmTMB", quietly = TRUE)) {
          glmmTMB::glmmTMB(formula_obj, family = glmmTMB::nbinom2, data = dat)
        } else {
          warning("glmmTMB not available, using Poisson glmer approximation for negative binomial.")
          lme4::glmer(formula_obj, family = stats::poisson(link = "log"), data = dat)
        }
      } else if (model_type == "poisson") {
        lme4::glmer(formula_obj, family = stats::poisson(link = "log"), data = dat)
      } else if (model_type == "logistic") {
        lme4::glmer(formula_obj, family = stats::binomial(link = "logit"), data = dat)
      } else if (model_type == "lm") {
        lme4::lmer(formula_obj, data = dat)
      } else if (model_type == "cox") {
        coxme::coxme(formula_obj, data = dat)
      }
    } else {
      if (model_type == "nb") {
        MASS::glm.nb(formula_obj, data = dat)
      } else if (model_type == "poisson") {
        stats::glm(formula_obj, family = stats::poisson(link = "log"), data = dat)
      } else if (model_type == "logistic") {
        stats::glm(formula_obj, family = stats::binomial(link = "logit"), data = dat)
      } else if (model_type == "lm") {
        stats::glm(formula_obj, family = stats::gaussian(), data = dat)
      } else if (model_type == "cox") {
        survival::coxph(formula_obj, data = dat)
      }
    }
  }

  if (MI_method == "first") {
    model <- fit_single_model(Data_Subset)
  } else {
    model <- NULL
  }

  # ---- Prediction data (population curve) ----
  x_pred <- seq(
    from = stats::quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
    to   = stats::quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
    length.out = 200
  )

  if (is.null(subgroups)) {
    pred_data <- data.frame(x_pred)
    colnames(pred_data) <- variable_x
  } else {
    subgroup_levels <- if (is.factor(Data_Subset[[subgroups]])) {
      levels(Data_Subset[[subgroups]])
    } else {
      sort(unique(Data_Subset[[subgroups]]))
    }
    pred_data <- expand.grid(x_pred, subgroup_levels)
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

  # add covariates at median/mode (using global_covariate_vars)
  if (!is.null(global_covariate_vars)) {
    for (cov in global_covariate_vars) {
      if (cov %in% colnames(Data_Subset)) {
        if (is.factor(Data_Subset[[cov]])) {
          pred_data[[cov]] <- as.factor(names(which.max(table(Data_Subset[[cov]]))))
        } else {
          pred_data[[cov]] <- stats::median(Data_Subset[[cov]], na.rm = TRUE)
        }
      } else {
        pred_data[[cov]] <- NA_real_
      }
    }
  }

  if (trial_factor == "Yes") {
    pred_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
  }
  if (followup_offset == "Yes") {
    pred_data[[followup_col]] <- 365
  }
  if (use_mixed) {
    pred_data[[random_intercept_var]] <- names(which.max(table(Data_Subset[[random_intercept_var]])))
  }

  # ---- Helper for mixed-model predictions ----
  get_mixed_model_predictions <- function(model, newdata, type = "link") {
    if (inherits(model, "merMod")) {
      fit <- predict(model, newdata = newdata, re.form = NA, type = type)
      mm  <- model.matrix(stats::formula(model, fixed.only = TRUE), newdata)
      vc  <- as.matrix(stats::vcov(model))
      se  <- sqrt(diag(mm %*% vc %*% t(mm)))
      list(fit = fit, se.fit = se)
    } else if (inherits(model, "glmmTMB")) {
      p <- predict(model, newdata = newdata, re.form = NA, type = type, se.fit = TRUE)
      list(fit = p$fit, se.fit = p$se.fit)
    } else if (inherits(model, "coxme")) {
      fit <- predict(model, newdata = newdata, type = "lp")
      mm  <- model.matrix(stats::formula(model, fixed.only = TRUE), newdata)
      vc  <- as.matrix(stats::vcov(model))
      se  <- sqrt(diag(mm %*% vc %*% t(mm)))
      list(fit = fit, se.fit = se)
    } else {
      stop("Unsupported mixed model class for prediction.")
    }
  }

  # ---- Predictions for population curve ----
  pred_on_scale <- function(fit_obj, newdata) {
    if (use_mixed) {
      preds <- get_mixed_model_predictions(fit_obj, newdata, type = "link")
      fit   <- preds$fit
      se    <- preds$se.fit
    } else {
      if (model_type == "cox") {
        fit <- predict(fit_obj, newdata = newdata, type = "lp")
        se  <- rep(NA_real_, length(fit))
      } else {
        p <- predict(fit_obj, newdata = newdata, type = "link", se.fit = TRUE)
        fit <- p$fit
        se  <- p$se.fit
      }
    }

    if (model_type %in% c("nb", "poisson", "cox")) {
      pred  <- exp(fit)
      lower <- exp(fit - 1.96 * se)
      upper <- exp(fit + 1.96 * se)
    } else if (model_type == "logistic") {
      pred  <- plogis(fit)
      lower <- plogis(fit - 1.96 * se)
      upper <- plogis(fit + 1.96 * se)
    } else { # lm
      pred  <- fit
      lower <- fit - 1.96 * se
      upper <- fit + 1.96 * se
    }

    list(prediction = pred, lower_ci = lower, upper_ci = upper)
  }

  if (MI_method == "first") {
    pred_res <- pred_on_scale(model, pred_data)
    pred_data$prediction <- pred_res$prediction
    pred_data$lower_ci   <- pred_res$lower_ci
    pred_data$upper_ci   <- pred_res$upper_ci
  } else {
    # Pool across imputations
    imps <- sort(unique(Data_Subset[[imp_col]]))
    all_fit  <- list()
    all_var  <- list()

    for (imp in imps) {
      dat_imp   <- subset(Data_Subset, get(imp_col) == imp)
      model_imp <- fit_single_model(dat_imp)
      pr        <- pred_on_scale(model_imp, pred_data)

      all_fit[[as.character(imp)]] <- log(pr$prediction)  # work on link scale
      all_var[[as.character(imp)]] <- ((log(pr$upper_ci) - log(pr$lower_ci)) / (2*1.96))^2
    }

    n_pt <- length(all_fit[[1]])
    mean_fit <- numeric(n_pt)
    mean_var <- numeric(n_pt)

    for (i in seq_len(n_pt)) {
      fits_i <- sapply(all_fit, `[`, i)
      vars_i <- sapply(all_var, `[`, i)
      m      <- length(fits_i)

      mean_fit[i] <- mean(fits_i)
      W <- mean(vars_i)
      if (MI_method == "Rubin") {
        B <- stats::var(fits_i)
        mean_var[i] <- W + (1 + 1/m) * B
      } else {
        mean_var[i] <- W
      }
    }

    if (model_type %in% c("nb", "poisson", "cox")) {
      pred_data$prediction <- exp(mean_fit)
      pred_data$lower_ci   <- exp(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci   <- exp(mean_fit + 1.96 * sqrt(mean_var))
    } else if (model_type == "logistic") {
      pred_data$prediction <- plogis(mean_fit)
      pred_data$lower_ci   <- plogis(mean_fit - 1.96 * sqrt(mean_var))
      pred_data$upper_ci   <- plogis(mean_fit + 1.96 * sqrt(mean_var))
    } else { # lm
      pred_data$prediction <- mean_fit
      pred_data$lower_ci   <- mean_fit - 1.96 * sqrt(mean_var)
      pred_data$upper_ci   <- mean_fit + 1.96 * sqrt(mean_var)
    }
  }

  # ---- Derivatives (numeric finite-difference, optional) ----
  derivative_data  <- NULL
  threshold        <- NA
  threshold_values <- NULL

  if (calculate_derivatives) {
    if (is.null(derivative_points)) {
      derivative_points <- seq(
        from = stats::quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
        to   = stats::quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE),
        length.out = 200
      )
    }

    if (is.null(subgroups)) {
      deriv_data <- data.frame(x_point = derivative_points)
      colnames(deriv_data)[1] <- variable_x
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

    # covariates template
    if (!is.null(global_covariate_vars)) {
      for (cov in global_covariate_vars) {
        if (cov %in% colnames(Data_Subset)) {
          if (is.factor(Data_Subset[[cov]])) {
            deriv_data[[cov]] <- as.factor(names(which.max(table(Data_Subset[[cov]]))))
          } else {
            deriv_data[[cov]] <- stats::median(Data_Subset[[cov]], na.rm = TRUE)
          }
        } else {
          deriv_data[[cov]] <- NA_real_
        }
      }
    }
    if (trial_factor == "Yes") {
      deriv_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
    }
    if (followup_offset == "Yes") {
      deriv_data[[followup_col]] <- 365
    }
    if (use_mixed) {
      deriv_data[[random_intercept_var]] <- names(which.max(table(Data_Subset[[random_intercept_var]])))
    }

    if (MI_method != "first") {
      warning("Derivatives currently implemented for MI_method='first' only. Returning NULL.")
      derivative_data <- NULL
    } else {
      eps <- 1e-5 * stats::sd(Data_Subset[[variable_x]], na.rm = TRUE)
      data_plus  <- deriv_data
      data_minus <- deriv_data
      data_plus[[variable_x]]  <- deriv_data[[variable_x]] + eps
      data_minus[[variable_x]] <- deriv_data[[variable_x]] - eps

      # link-scale predictions
      if (use_mixed) {
        lp_center <- get_mixed_model_predictions(model, deriv_data, type = "link")$fit
        lp_plus   <- get_mixed_model_predictions(model, data_plus,  type = "link")$fit
        lp_minus  <- get_mixed_model_predictions(model, data_minus, type = "link")$fit
      } else {
        if (model_type == "cox") {
          lp_center <- predict(model, newdata = deriv_data, type = "lp")
          lp_plus   <- predict(model, newdata = data_plus,  type = "lp")
          lp_minus  <- predict(model, newdata = data_minus, type = "lp")
        } else {
          lp_center <- predict(model, newdata = deriv_data, type = "link")
          lp_plus   <- predict(model, newdata = data_plus,  type = "link")
          lp_minus  <- predict(model, newdata = data_minus, type = "link")
        }
      }

      d_link <- (lp_plus - lp_minus) / (2 * eps)

      if (model_type %in% c("nb", "poisson", "cox")) {
        resp_center <- exp(lp_center)
        d_resp      <- resp_center * d_link
      } else if (model_type == "logistic") {
        p_center <- plogis(lp_center)
        d_resp   <- p_center * (1 - p_center) * d_link
      } else {
        d_resp <- d_link
      }

      se_der <- abs(d_resp) * 0.2
      lower  <- d_resp - 1.96 * se_der
      upper  <- d_resp + 1.96 * se_der

      derivative_data <- deriv_data
      derivative_data$derivative          <- d_resp
      derivative_data$se_derivative       <- se_der
      derivative_data$lower_ci_derivative <- lower
      derivative_data$upper_ci_derivative <- upper
      derivative_data$x_point             <- deriv_data[[variable_x]]

      if (is.null(subgroups)) {
        derivative_data <- derivative_data[order(derivative_data$x_point), ]
        sig_pos <- which(derivative_data$lower_ci_derivative > 0)
        if (length(sig_pos) > 0) {
          threshold <- derivative_data$x_point[min(sig_pos)]
        }
      } else {
        derivative_data$significant_positive <- derivative_data$lower_ci_derivative > 0
        lvl <- unique(derivative_data[[subgroups]])
        threshold_values <- rep(NA_real_, length(lvl))
        names(threshold_values) <- as.character(lvl)
        for (lv in lvl) {
          sub_d <- derivative_data[derivative_data[[subgroups]] == lv, ]
          sub_d <- sub_d[order(sub_d$x_point), ]
          sig   <- which(sub_d$significant_positive)
          if (length(sig) > 0) {
            threshold_values[as.character(lv)] <- sub_d$x_point[min(sig)]
          }
        }
      }
    }
  }

  # ---- Basic population plot ----
  if (is.null(plot_options)) plot_options <- list()

  x_lab <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x
  y_lab <- if (!is.null(plot_options$y_lab)) {
    plot_options$y_lab
  } else if (model_type %in% c("nb", "poisson")) {
    "Predicted rate"
  } else if (model_type == "logistic") {
    "Predicted probability"
  } else if (model_type == "lm") {
    "Predicted mean"
  } else {
    "Predicted hazard ratio"
  }

  if (!is.null(subgroups)) {
    plot <- ggplot2::ggplot(pred_data,
                            ggplot2::aes_string(x = variable_x,
                                                y = "prediction",
                                                color = subgroups,
                                                fill  = subgroups,
                                                linetype = subgroups)) +
      ggplot2::geom_line(linewidth = plot_options$line_size %||% 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                           alpha = plot_options$ribbon_alpha %||% 0.3)
  } else {
    line_color  <- (plot_options$colors     %||% "black")
    fill_colors <- (plot_options$fill_colors %||% "grey70")
    plot <- ggplot2::ggplot(pred_data,
                            ggplot2::aes_string(x = variable_x, y = "prediction")) +
      ggplot2::geom_line(linewidth = plot_options$line_size %||% 1,
                         color = line_color) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                           fill  = fill_colors,
                           alpha = plot_options$ribbon_alpha %||% 0.3)
  }

  plot <- plot +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::theme_minimal(base_size = 10)

  coord_args <- list()
  if (!is.null(plot_options$y_limits)) coord_args$ylim <- plot_options$y_limits
  if (!is.null(plot_options$x_limits)) coord_args$xlim <- plot_options$x_limits
  if (length(coord_args) > 0) {
    plot <- plot + do.call(ggplot2::coord_cartesian, coord_args)
  }

  # ---- group_fits: per-trial independent refits with same knots ----
  group_fits_df <- NULL

  if (group_fits) {
    if (MI_method != "first") {
      warning("group_fits currently implemented for MI_method='first' only. Skipping.")
    } else if (is.null(trial_col)) {
      warning("group_fits requires trial_col (e.g., 'Enrolled_Trial_name'). Skipping.")
    } else {
      groups <- sort(unique(Data_Subset[[trial_col]]))
      all_g  <- list()

      for (g in groups) {
        dat_g <- subset(Data_Subset, get(trial_col) == g)
        dat_g <- droplevels(dat_g)

        if (nrow(dat_g) < (knot_n + 2)) {
          warning(sprintf("Group '%s': too few observations, skipping in group_fits.", g))
          next
        }

        # find 1-level factors (for formula, not for prediction template)
        one_level_factors <- vapply(
          dat_g,
          function(col) is.factor(col) && nlevels(col) < 2,
          logical(1)
        )
        bad_vars <- names(dat_g)[one_level_factors]

        covariates_group <- expanded_covariables
        if (!is.null(covariates_group) && length(bad_vars) > 0) {
          covariates_group <- covariates_group[!sapply(covariates_group, function(term) {
            pieces <- trimws(unlist(strsplit(gsub(":", "*", term), "\\*")))
            any(pieces %in% bad_vars)
          })]
        }

        # group-specific spline term with SAME global knots
        spline_term_g <- paste0(
          "rcs(", variable_x, ", c(",
          paste(round(global_knots, 4), collapse = ", "),
          "))"
        )

        rhs_terms <- c(
          spline_term_g,
          if (!is.null(covariates_group) && length(covariates_group) > 0)
            paste(covariates_group, collapse = " + ")
          else
            NULL,
          if (followup_offset == "Yes") paste0("offset(log(", followup_col, "))") else NULL
        )
        rhs_g <- paste(rhs_terms, collapse = " + ")

        form_g_str <- if (model_type == "cox") {
          paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", rhs_g)
        } else {
          paste0(outcome_var, " ~ ", rhs_g)
        }
        form_g <- stats::as.formula(form_g_str)

        # fit per-group model (no random effects here)
        mod_g <- tryCatch({
          if (model_type == "nb") {
            MASS::glm.nb(form_g, data = dat_g)
          } else if (model_type == "poisson") {
            stats::glm(form_g, family = stats::poisson(link = "log"), data = dat_g)
          } else if (model_type == "logistic") {
            stats::glm(form_g, family = stats::binomial(link = "logit"), data = dat_g)
          } else if (model_type == "lm") {
            stats::glm(form_g, family = stats::gaussian(), data = dat_g)
          } else if (model_type == "cox") {
            survival::coxph(form_g, data = dat_g)
          }
        }, error = function(e) {
          warning(sprintf("Group '%s': model failed (%s). Skipping this group in group_fits.",
                          g, e$message))
          NULL
        })

        if (is.null(mod_g)) next

        # group-specific x-range
        x_g <- seq(
          min(dat_g[[variable_x]], na.rm = TRUE),
          max(dat_g[[variable_x]], na.rm = TRUE),
          length.out = 80
        )

        new_g <- data.frame(x_g)
        colnames(new_g)[1] <- variable_x

        # add global covariate columns for prediction template (same columns across groups)
        if (!is.null(global_covariate_vars)) {
          for (cv in global_covariate_vars) {
            if (cv %in% names(dat_g)) {
              if (is.factor(dat_g[[cv]])) {
                if (nlevels(dat_g[[cv]]) >= 2) {
                  new_g[[cv]] <- as.factor(names(which.max(table(dat_g[[cv]]))))
                } else {
                  new_g[[cv]] <- dat_g[[cv]][1]
                }
              } else {
                new_g[[cv]] <- stats::median(dat_g[[cv]], na.rm = TRUE)
              }
            } else {
              new_g[[cv]] <- NA_real_
            }
          }
        }

        if (followup_offset == "Yes") {
          new_g[[followup_col]] <- 365
        }

        # predictions with CI
        if (model_type == "cox") {
          lp <- predict(mod_g, newdata = new_g, type = "lp", se.fit = TRUE)
          fit <- lp$fit
          se  <- lp$se.fit
          new_g$prediction <- exp(fit)
          new_g$lower_ci   <- exp(fit - 1.96 * se)
          new_g$upper_ci   <- exp(fit + 1.96 * se)
        } else {
          p <- predict(mod_g, newdata = new_g, type = "link", se.fit = TRUE)
          fit <- p$fit
          se  <- p$se.fit
          if (model_type %in% c("nb", "poisson")) {
            new_g$prediction <- exp(fit)
            new_g$lower_ci   <- exp(fit - 1.96 * se)
            new_g$upper_ci   <- exp(fit + 1.96 * se)
          } else if (model_type == "logistic") {
            new_g$prediction <- plogis(fit)
            new_g$lower_ci   <- plogis(fit - 1.96 * se)
            new_g$upper_ci   <- plogis(fit + 1.96 * se)
          } else { # lm
            new_g$prediction <- fit
            new_g$lower_ci   <- fit - 1.96 * se
            new_g$upper_ci   <- fit + 1.96 * se
          }
        }

        new_g[[trial_col]] <- g
        all_g[[g]] <- new_g
      }

      if (length(all_g) > 0) {
        # Ensure identical column order/names before rbind
        all_names <- unique(unlist(lapply(all_g, names)))
        all_g_aligned <- lapply(all_g, function(df) {
          missing_cols <- setdiff(all_names, names(df))
          for (mc in missing_cols) df[[mc]] <- NA
          df[, all_names]
        })
        group_fits_df <- do.call(rbind, all_g_aligned)
        rownames(group_fits_df) <- NULL
      }
    }
  }

  # ---- Return ----
  list(
    predictions = pred_data,
    model       = model,
    plot        = plot,
    plot_data   = pred_data,
    prediction_range_values = c(
      stats::quantile(Data_Subset[[variable_x]], prediction_range[1], na.rm = TRUE),
      stats::quantile(Data_Subset[[variable_x]], prediction_range[2], na.rm = TRUE)
    ),
    derivatives       = derivative_data,
    derivative_method = derivative_method,
    threshold         = if (!calculate_derivatives) NULL else {
      if (!is.null(subgroups)) threshold_values else threshold
    },
    derivative_plot = NULL,
    group_fits      = group_fits_df
  )
}
