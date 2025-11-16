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
                      MI_method = "first",      # "first", "average", "Rubin"
                      model_type = "nb",        # "nb", "poisson", "logistic", "lm", "cox"
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
                      derivative_method  = "numeric",  # kept for compatibility
                      group_fits = FALSE) {

  ## ────────────────────────────────────────────────────────────────
  ## Packages
  ## ────────────────────────────────────────────────────────────────
  requireNamespace("rms",    quietly = TRUE)
  requireNamespace("Hmisc",  quietly = TRUE)
  requireNamespace("MASS",   quietly = TRUE)
  requireNamespace("ggplot2",quietly = TRUE)
  requireNamespace("dplyr",  quietly = TRUE)

  ## ────────────────────────────────────────────────────────────────
  ## Basic checks
  ## ────────────────────────────────────────────────────────────────
  if (!model_type %in% c("nb","poisson","logistic","lm","cox"))
    stop("model_type must be one of 'nb','poisson','logistic','lm','cox'")

  if (followup_offset == "Yes" && is.null(followup_col))
    stop("followup_col must be provided when followup_offset = 'Yes'")

  if (trial_factor == "Yes" && is.null(trial_col))
    stop("trial_col must be provided when trial_factor = 'Yes'")

  if (model_type == "cox" && (is.null(time_col) || is.null(event_col)))
    stop("time_col and event_col must be provided for Cox models")

  if (!random_intercept %in% c("Yes","No"))
    stop("random_intercept must be 'Yes' or 'No'")

  if (!random_slope %in% c("Yes","No"))
    stop("random_slope must be 'Yes' or 'No'")

  has_random_effects <- (random_intercept == "Yes" || random_slope == "Yes")

  if (has_random_effects) {
    requireNamespace("lme4", quietly = TRUE)
    if (model_type == "nb" && !requireNamespace("glmmTMB", quietly = TRUE)) {
      warning("Package 'glmmTMB' not available. Will approximate NB with Poisson mixed model.")
    }
    if (model_type == "cox" && !requireNamespace("coxme", quietly = TRUE)) {
      stop("Package 'coxme' is required for Cox models with random effects.")
    }
  }

  if (has_random_effects && is.null(random_intercept_var)) {
    stop("If random_intercept = 'Yes' or random_slope = 'Yes', random_intercept_var must be provided")
  }

  if (length(prediction_range) != 2 ||
      any(prediction_range < 0) || any(prediction_range > 1) ||
      prediction_range[1] >= prediction_range[2]) {
    stop("prediction_range must be c(p_low, p_high) with 0 <= p_low < p_high <= 1")
  }

  if (group_fits && MI_method != "first") {
    stop("group_fits = TRUE is currently implemented only for MI_method = 'first'")
  }

  ## ────────────────────────────────────────────────────────────────
  ## Select imputations
  ## ────────────────────────────────────────────────────────────────
  if (MI_method == "first") {
    Data_Subset <- subset(data, get(imp_col) == 1)
  } else if (MI_method %in% c("average","Rubin")) {
    Data_Subset <- data
  } else {
    stop("MI_method must be 'first', 'average', or 'Rubin'")
  }

  ## ────────────────────────────────────────────────────────────────
  ## Subgroups handling
  ## ────────────────────────────────────────────────────────────────
  if (!is.null(subgroups)) {
    if (subgroup_as_factor && !is.factor(Data_Subset[[subgroups]])) {
      Data_Subset[[subgroups]] <- factor(Data_Subset[[subgroups]])
    }
    if (!is.null(subgroup_labels)) {
      levels(Data_Subset[[subgroups]]) <- subgroup_labels
    }
  }

  ## ────────────────────────────────────────────────────────────────
  ## Build spline term – first compute global knots via Hmisc
  ## ────────────────────────────────────────────────────────────────
  x_full <- Data_Subset[[variable_x]]
  x_full <- x_full[is.finite(x_full)]

  if (length(x_full) < knot_n) {
    stop("Not enough non-missing values of variable_x to place knots")
  }

  # IMPORTANT: use Hmisc::rcspline.eval
  rcs_tmp     <- Hmisc::rcspline.eval(x_full, nk = knot_n, inclx = TRUE)
  global_knots <- attr(rcs_tmp, "knots")

  # rcs() with explicit knots
  spline_term_global <- paste0(
    "rcs(", variable_x, ", c(",
    paste(round(global_knots, 4), collapse = ", "),
    "))"
  )

  # If subgroups: interaction
  if (!is.null(subgroups)) {
    spline_term_global <- paste0(spline_term_global, " * ", subgroups)
  }

  ## ────────────────────────────────────────────────────────────────
  ## Covariables expansion / validation
  ## ────────────────────────────────────────────────────────────────
  expand_terms <- function(term) {
    # do not expand actual spline calls
    if (grepl("^rcs\\(", term)) return(term)

    if (grepl("\\*", term)) {
      parts <- trimws(unlist(strsplit(term, "\\*")))
      all_combos <- lapply(seq_along(parts), function(k) {
        combn(parts, k, FUN = function(x) paste(x, collapse = ":"))
      })
      unique(unlist(all_combos))
    } else {
      term
    }
  }

  expanded_covariables <- NULL
  covariates_str <- ""
  if (!is.null(covariables)) {
    expanded_covariables <- unique(unlist(lapply(covariables, expand_terms)))

    # validate that all base variables exist
    base_vars <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
    base_vars <- trimws(base_vars)
    missing_vars <- base_vars[!base_vars %in% names(Data_Subset)]
    if (length(missing_vars) > 0) {
      stop("Covariates not found in data: ", paste(missing_vars, collapse = ", "))
    }

    covariates_str <- paste(expanded_covariables, collapse = " + ")
  }

  ## ────────────────────────────────────────────────────────────────
  ## Trial factor, offset, random effects strings
  ## ────────────────────────────────────────────────────────────────
  trial_str <- ""
  if (trial_factor == "Yes") {
    trial_str <- paste0("as.factor(", trial_col, ")")
  }

  offset_str <- ""
  if (followup_offset == "Yes" && model_type %in% c("nb","poisson")) {
    offset_str <- paste0("offset(log(", followup_col, "))")
  }

  random_str <- ""
  if (has_random_effects) {
    random_parts <- c()
    if (random_intercept == "Yes") {
      random_parts <- c(random_parts, "1")
    }
    if (random_slope == "Yes") {
      slope_terms <- unique(na.omit(c(predictor_vars_random_slope, covariables_random_slope)))
      if (length(slope_terms) > 0) {
        random_parts <- c(random_parts, paste0("0 + ", slope_terms))
      }
    }
    if (length(random_parts) > 0) {
      random_str <- paste0(
        "(", paste(random_parts, collapse = " + "), " | ", random_intercept_var, ")"
      )
    }
  }

  ## ────────────────────────────────────────────────────────────────
  ## Build full formula for the GLOBAL model
  ## ────────────────────────────────────────────────────────────────
  rhs_terms <- c(
    spline_term_global,
    if (covariates_str != "") covariates_str else NULL,
    if (trial_str != "")      trial_str      else NULL,
    if (offset_str != "")     offset_str     else NULL,
    if (random_str != "")     random_str     else NULL
  )

  rhs <- paste(rhs_terms, collapse = " + ")

  if (model_type == "cox") {
    form_str <- paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", rhs)
  } else {
    form_str <- paste0(outcome_var, " ~ ", rhs)
  }

  formula_obj <- as.formula(form_str)

  ## ────────────────────────────────────────────────────────────────
  ## Fit GLOBAL model (one per imputation or one overall)
  ## ────────────────────────────────────────────────────────────────
  fit_one_model <- function(dat) {
    if (has_random_effects) {
      if (model_type == "nb") {
        if (requireNamespace("glmmTMB", quietly = TRUE)) {
          glmmTMB::glmmTMB(formula_obj, family = glmmTMB::nbinom2, data = dat)
        } else {
          lme4::glmer(formula_obj, family = poisson(link = "log"), data = dat)
        }
      } else if (model_type == "poisson") {
        lme4::glmer(formula_obj, family = poisson(link = "log"), data = dat)
      } else if (model_type == "logistic") {
        lme4::glmer(formula_obj, family = binomial(link = "logit"), data = dat)
      } else if (model_type == "lm") {
        lme4::lmer(formula_obj, data = dat)
      } else if (model_type == "cox") {
        coxme::coxme(formula_obj, data = dat)
      }
    } else {
      if (model_type == "nb") {
        MASS::glm.nb(formula_obj, data = dat)
      } else if (model_type == "poisson") {
        glm(formula_obj, family = poisson(link = "log"), data = dat)
      } else if (model_type == "logistic") {
        glm(formula_obj, family = binomial(link = "logit"), data = dat)
      } else if (model_type == "lm") {
        glm(formula_obj, family = gaussian(), data = dat)
      } else if (model_type == "cox") {
        survival::coxph(formula_obj, data = dat)
      }
    }
  }

  if (MI_method == "first") {
    model <- fit_one_model(Data_Subset)
    model_list <- list("1" = model)
  } else {
    imp_ids <- sort(unique(Data_Subset[[imp_col]]))
    model_list <- lapply(imp_ids, function(imp) {
      dat_imp <- subset(Data_Subset, get(imp_col) == imp)
      fit_one_model(dat_imp)
    })
    names(model_list) <- as.character(imp_ids)
    model <- model_list[[1]]  # for returning a representative model
  }

  ## ────────────────────────────────────────────────────────────────
  ## Create GLOBAL prediction grid
  ## ────────────────────────────────────────────────────────────────
  x_range <- quantile(Data_Subset[[variable_x]],
                      probs = prediction_range,
                      na.rm = TRUE)
  x_values <- seq(from = x_range[1], to = x_range[2], length.out = 100)

  if (is.null(subgroups)) {
    pred_data <- data.frame(x_values)
    colnames(pred_data) <- variable_x
  } else {
    sg_levels <- if (is.factor(Data_Subset[[subgroups]])) {
      levels(Data_Subset[[subgroups]])
    } else {
      sort(unique(Data_Subset[[subgroups]]))
    }
    pred_data <- expand.grid(
      x_values,
      sg_levels
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

  # fill covariates with median/mode
  if (!is.null(expanded_covariables)) {
    cov_vars <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
    cov_vars <- trimws(cov_vars)
    for (cv in cov_vars) {
      if (cv %in% colnames(Data_Subset)) {
        if (is.factor(Data_Subset[[cv]])) {
          pred_data[[cv]] <- factor(
            names(which.max(table(Data_Subset[[cv]]))),
            levels = levels(Data_Subset[[cv]])
          )
        } else {
          pred_data[[cv]] <- stats::median(Data_Subset[[cv]], na.rm = TRUE)
        }
      }
    }
  }

  # trial factor: most frequent trial
  if (trial_factor == "Yes") {
    pred_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
  }

  # follow-up offset: standardised to 365 days
  if (followup_offset == "Yes") {
    pred_data[[followup_col]] <- 365
  }

  # random_intercept_var for prediction (population-level, re.form = NA)
  if (has_random_effects) {
    pred_data[[random_intercept_var]] <-
      names(which.max(table(Data_Subset[[random_intercept_var]])))
  }

  ## ────────────────────────────────────────────────────────────────
  ## Prediction helpers
  ## ────────────────────────────────────────────────────────────────
  predict_one <- function(model, newdata, model_type, has_random) {
    if (model_type == "cox") {
      if (inherits(model, "coxme")) {
        lp <- predict(model, newdata = newdata, type = "lp", re.form = NA)
        # no straightforward se.fit
        list(fit = lp, se.fit = rep(NA_real_, length(lp)), link = "log")
      } else {
        lp <- predict(model, newdata = newdata, type = "lp", se.fit = FALSE)
        list(fit = lp, se.fit = rep(NA_real_, length(lp)), link = "log")
      }
    } else if (has_random) {
      if (inherits(model, "glmmTMB")) {
        p <- predict(model, newdata = newdata, type = "link",
                     se.fit = TRUE, re.form = NA)
        list(fit = p$fit, se.fit = p$se.fit, link = "logit_or_log")
      } else if (inherits(model, "merMod")) {
        # lme4 doesn't give se.fit; approximate from vcov of fixed effects
        mm <- model.matrix(stats::formula(model, fixed.only = TRUE), newdata)
        beta <- lme4::fixef(model)
        vc  <- stats::vcov(model)
        fit <- as.vector(mm %*% beta)
        se  <- sqrt(diag(mm %*% vc %*% t(mm)))
        list(fit = fit, se.fit = se, link = "logit_or_log")
      } else {
        stop("Unknown mixed model class for prediction")
      }
    } else {  # standard glm/glm.nb/coxph
      if (model_type == "cox") {
        lp <- predict(model, newdata = newdata, type = "lp", se.fit = FALSE)
        list(fit = lp, se.fit = rep(NA_real_, length(lp)), link = "log")
      } else {
        p <- stats::predict(model, newdata = newdata, type = "link", se.fit = TRUE)
        list(fit = p$fit, se.fit = p$se.fit, link = "logit_or_log")
      }
    }
  }

  transform_link <- function(fit, se, model_type) {
    if (model_type %in% c("nb","poisson","cox")) {
      mu <- exp(fit)
      if (all(is.na(se))) {
        lower <- upper <- NA_real_
      } else {
        lower <- exp(fit - 1.96 * se)
        upper <- exp(fit + 1.96 * se)
      }
    } else if (model_type == "logistic") {
      mu <- stats::plogis(fit)
      if (all(is.na(se))) {
        lower <- upper <- NA_real_
      } else {
        lower <- stats::plogis(fit - 1.96 * se)
        upper <- stats::plogis(fit + 1.96 * se)
      }
    } else { # lm
      mu <- fit
      if (all(is.na(se))) {
        lower <- upper <- NA_real_
      } else {
        lower <- fit - 1.96 * se
        upper <- fit + 1.96 * se
      }
    }
    list(mu = mu, lower = lower, upper = upper)
  }

  ## ────────────────────────────────────────────────────────────────
  ## GLOBAL predictions (MI pooled if needed)
  ## ────────────────────────────────────────────────────────────────
  if (MI_method == "first") {
    pr <- predict_one(model, pred_data, model_type, has_random_effects)
    tr <- transform_link(pr$fit, pr$se.fit, model_type)
    pred_data$prediction <- tr$mu
    pred_data$lower_ci   <- tr$lower
    pred_data$upper_ci   <- tr$upper
  } else {
    imp_ids <- names(model_list)
    all_fit <- lapply(model_list, function(m) {
      pr <- predict_one(m, pred_data, model_type, has_random_effects)
      pr
    })

    n_pred <- nrow(pred_data)
    m_imp  <- length(all_fit)
    fit_mat <- sapply(all_fit, `[[`, "fit")
    se_mat  <- sapply(all_fit, `[[`, "se.fit")

    mean_fit <- rowMeans(fit_mat)
    mean_var <- numeric(n_pred)

    for (i in seq_len(n_pred)) {
      fits_i <- fit_mat[i, ]
      ses_i  <- se_mat[i, ]
      W <- mean(ses_i^2, na.rm = TRUE)
      B <- stats::var(fits_i, na.rm = TRUE)
      if (MI_method == "Rubin") {
        mean_var[i] <- W + (1 + 1/m_imp) * B
      } else { # "average"
        mean_var[i] <- W
      }
    }
    se_pooled <- sqrt(mean_var)
    tr <- transform_link(mean_fit, se_pooled, model_type)
    pred_data$prediction <- tr$mu
    pred_data$lower_ci   <- tr$lower
    pred_data$upper_ci   <- tr$upper
  }

  ## ────────────────────────────────────────────────────────────────
  ## GROUP-specific fits (one model per trial / random_intercept_var)
  ## – each group refitted independently, using GLOBAL knots
  ## ────────────────────────────────────────────────────────────────
  group_fits_df <- NULL

  if (group_fits) {
    if (is.null(trial_col)) {
      stop("group_fits=TRUE requires trial_col (grouping variable for separate fits)")
    }

    group_var <- trial_col
    groups <- sort(unique(Data_Subset[[group_var]]))

    group_list <- vector("list", length(groups))
    names(group_list) <- groups

    for (g in groups) {
      dat_g <- subset(Data_Subset, get(group_var) == g)

      # if too few observations, skip
      if (nrow(dat_g) < knot_n + 3) next

      # local x-range for this group
      xg <- dat_g[[variable_x]]
      xg <- xg[is.finite(xg)]
      if (length(xg) < knot_n) next

      xg_range <- quantile(xg, probs = prediction_range, na.rm = TRUE)
      xg_seq   <- seq(from = xg_range[1], to = xg_range[2], length.out = 100)

      new_g <- data.frame(xg_seq)
      colnames(new_g) <- variable_x

      # fill covariates for this group
      if (!is.null(expanded_covariables)) {
        cov_vars <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
        cov_vars <- trimws(cov_vars)
        for (cv in cov_vars) {
          if (cv %in% colnames(dat_g)) {
            if (is.factor(dat_g[[cv]])) {
              new_g[[cv]] <- factor(
                names(which.max(table(dat_g[[cv]]))),
                levels = levels(dat_g[[cv]])
              )
            } else {
              new_g[[cv]] <- stats::median(dat_g[[cv]], na.rm = TRUE)
            }
          }
        }
      }

      # followup offset fixed to 365
      if (followup_offset == "Yes") {
        new_g[[followup_col]] <- 365
      }

      # build group formula WITHOUT random effects, but WITH same knots
      spline_term_g <- paste0(
        "rcs(", variable_x, ", c(",
        paste(round(global_knots, 4), collapse = ", "),
        "))"
      )

      rhs_g_terms <- c(
        spline_term_g,
        if (covariates_str != "") covariates_str else NULL
      )
      rhs_g <- paste(rhs_g_terms, collapse = " + ")

      form_g_str <- if (model_type == "cox") {
        paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", rhs_g)
      } else {
        paste0(outcome_var, " ~ ", rhs_g)
      }

      form_g <- as.formula(form_g_str)

      # fit single-group model
      if (model_type == "nb") {
        mod_g <- MASS::glm.nb(form_g, data = dat_g)
      } else if (model_type == "poisson") {
        mod_g <- glm(form_g, family = poisson(link = "log"), data = dat_g)
      } else if (model_type == "logistic") {
        mod_g <- glm(form_g, family = binomial(link = "logit"), data = dat_g)
      } else if (model_type == "lm") {
        mod_g <- glm(form_g, family = gaussian(), data = dat_g)
      } else if (model_type == "cox") {
        mod_g <- survival::coxph(form_g, data = dat_g)
      }

      pr_g <- stats::predict(mod_g, newdata = new_g, type = "link", se.fit = TRUE)

      # transform to response scale
      if (model_type %in% c("nb","poisson","cox")) {
        mu_g <- exp(pr_g$fit)
        lower_g <- exp(pr_g$fit - 1.96 * pr_g$se.fit)
        upper_g <- exp(pr_g$fit + 1.96 * pr_g$se.fit)
      } else if (model_type == "logistic") {
        mu_g <- stats::plogis(pr_g$fit)
        lower_g <- stats::plogis(pr_g$fit - 1.96 * pr_g$se.fit)
        upper_g <- stats::plogis(pr_g$fit + 1.96 * pr_g$se.fit)
      } else {
        mu_g <- pr_g$fit
        lower_g <- pr_g$fit - 1.96 * pr_g$se.fit
        upper_g <- pr_g$fit + 1.96 * pr_g$se.fit
      }

      new_g$prediction <- mu_g
      new_g$lower_ci   <- lower_g
      new_g$upper_ci   <- upper_g
      new_g[[group_var]] <- g

      group_list[[g]] <- new_g
    }

    group_fits_df <- dplyr::bind_rows(group_list)
  }

  ## ────────────────────────────────────────────────────────────────
  ## Derivative calculation (simplified: numeric, no random effects)
  ## ────────────────────────────────────────────────────────────────
  derivative_data <- NULL

  if (calculate_derivatives) {
    if (has_random_effects || MI_method != "first") {
      warning("Derivatives are currently implemented only for MI_method='first' and models without random effects. Returning NULL for derivatives.")
      derivative_data <- NULL
    } else {
      # numeric derivatives on response scale
      if (is.null(derivative_points)) {
        derivative_points <- x_values
      }

      deriv_data <- data.frame(x = derivative_points)
      colnames(deriv_data) <- variable_x

      # fill covariates as in pred_data
      if (!is.null(expanded_covariables)) {
        cov_vars <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
        cov_vars <- trimws(cov_vars)
        for (cv in cov_vars) {
          if (cv %in% colnames(Data_Subset)) {
            if (is.factor(Data_Subset[[cv]])) {
              deriv_data[[cv]] <- factor(
                names(which.max(table(Data_Subset[[cv]]))),
                levels = levels(Data_Subset[[cv]])
              )
            } else {
              deriv_data[[cv]] <- stats::median(Data_Subset[[cv]], na.rm = TRUE)
            }
          }
        }
      }

      if (followup_offset == "Yes") {
        deriv_data[[followup_col]] <- 365
      }
      if (trial_factor == "Yes") {
        deriv_data[[trial_col]] <- names(which.max(table(Data_Subset[[trial_col]])))
      }

      # step size
      eps <- 1e-5 * stats::sd(Data_Subset[[variable_x]], na.rm = TRUE)
      data_plus  <- deriv_data
      data_minus <- deriv_data
      data_plus[[variable_x]]  <- deriv_data[[variable_x]] + eps
      data_minus[[variable_x]] <- deriv_data[[variable_x]] - eps

      pr0 <- predict_one(model, deriv_data, model_type, FALSE)
      prp <- predict_one(model, data_plus,  model_type, FALSE)
      prm <- predict_one(model, data_minus, model_type, FALSE)

      tr0 <- transform_link(pr0$fit, pr0$se.fit, model_type)
      trp <- transform_link(prp$fit, prp$se.fit, model_type)
      trm <- transform_link(prm$fit, prm$se.fit, model_type)

      d_resp <- (trp$mu - trm$mu) / (2 * eps)
      se_d   <- abs(d_resp) * 0.2   # rough 20% relative error
      lower_d <- d_resp - 1.96 * se_d
      upper_d <- d_resp + 1.96 * se_d

      derivative_data <- data.frame(
        x_point             = deriv_data[[variable_x]],
        derivative          = d_resp,
        se_derivative       = se_d,
        lower_ci_derivative = lower_d,
        upper_ci_derivative = upper_d
      )
    }
  }

  ## ────────────────────────────────────────────────────────────────
  ## Build basic GLOBAL plot
  ## ────────────────────────────────────────────────────────────────
  if (is.null(plot_options)) plot_options <- list()

  x_lab <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x
  y_lab <- if (!is.null(plot_options$y_lab)) {
    plot_options$y_lab
  } else if (model_type %in% c("nb","poisson")) {
    "Predicted rate"
  } else if (model_type == "logistic") {
    "Predicted probability"
  } else if (model_type == "lm") {
    "Predicted mean"
  } else {
    "Predicted hazard ratio"
  }

  line_size    <- plot_options$line_size    %||% 1
  ribbon_alpha <- plot_options$ribbon_alpha %||% 0.3

  `%||%` <- function(a,b) if (!is.null(a)) a else b

  if (is.null(subgroups)) {
    line_color  <- plot_options$colors      %||% "black"
    fill_color  <- plot_options$fill_colors %||% "grey80"
    p <- ggplot2::ggplot(pred_data,
                         ggplot2::aes_string(x = variable_x, y = "prediction")) +
      ggplot2::geom_line(linewidth = line_size, color = line_color) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
        fill = fill_color, alpha = ribbon_alpha, color = NA
      )
  } else {
    p <- ggplot2::ggplot(pred_data,
                         ggplot2::aes_string(x = variable_x,
                                             y = "prediction",
                                             color = subgroups,
                                             fill  = subgroups,
                                             linetype = subgroups)) +
      ggplot2::geom_line(linewidth = line_size) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
        alpha = ribbon_alpha, color = NA
      )
  }

  p <- p +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::theme_minimal()

  if (!is.null(plot_options$y_limits) || !is.null(plot_options$x_limits)) {
    coord_args <- list()
    if (!is.null(plot_options$y_limits)) coord_args$ylim <- plot_options$y_limits
    if (!is.null(plot_options$x_limits)) coord_args$xlim <- plot_options$x_limits
    p <- p + do.call(ggplot2::coord_cartesian, coord_args)
  }

  ## ────────────────────────────────────────────────────────────────
  ## Return
  ## ────────────────────────────────────────────────────────────────
  list(
    predictions             = pred_data,
    model                   = model,
    model_list              = if (MI_method == "first") NULL else model_list,
    plot                    = p,
    plot_data               = pred_data,
    prediction_range_values = x_range,
    derivatives             = derivative_data,
    derivative_method       = if (calculate_derivatives) derivative_method else NULL,
    group_fits              = group_fits_df,
    global_knots            = global_knots,
    formula                 = formula_obj
  )
}
