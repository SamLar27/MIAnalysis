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

                      ## NEW spline arguments
                      spline_knots_n          = 4,
                      spline_knots_percentile = NULL,

                      imp_col     = ".imp",
                      MI_method   = "first",   # "first", "average", "Rubin"
                      model_type  = "nb",      # "nb", "poisson", "logistic", "lm", "cox"

                      followup_offset = "No",
                      followup_col    = NULL,

                      trial_factor    = "No",
                      trial_col       = NULL,

                      time_col  = NULL,        # for Cox
                      event_col = NULL,        # for Cox

                      random_intercept     = "No",
                      random_intercept_var = NULL,

                      random_slope  = "No",
                      predictor_vars_random_slope = NULL,
                      covariables_random_slope    = NULL,

                      prediction_range = c(0.01, 0.99),
                      calculate_derivatives = FALSE,
                      derivative_points     = NULL,

                      plot_options = NULL,
                      group_fits   = FALSE) {
  ## ------------------------------------------------------------------
  ## 0. BASIC CHECKS & SETUP
  ## ------------------------------------------------------------------
  requireNamespace("rms", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  if (!model_type %in% c("nb", "poisson", "logistic", "lm", "cox")) {
    stop("model_type must be one of 'nb', 'poisson', 'logistic', 'lm', 'cox'")
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
  if (random_intercept == "Yes" && is.null(random_intercept_var)) {
    stop("If random_intercept = 'Yes', random_intercept_var must be provided")
  }

  if (!random_slope %in% c("Yes", "No")) {
    stop("random_slope must be 'Yes' or 'No'")
  }
  if (random_slope == "Yes" && is.null(predictor_vars_random_slope)) {
    stop("If random_slope = 'Yes', predictor_vars_random_slope must be provided")
  }

  if (length(prediction_range) != 2 ||
      any(prediction_range < 0) ||
      any(prediction_range > 1) ||
      prediction_range[1] >= prediction_range[2]) {
    stop("prediction_range must be c(p1, p2) with 0 <= p1 < p2 <= 1")
  }

  if (calculate_derivatives && MI_method != "first") {
    stop("Derivatives currently implemented only for MI_method = 'first'")
  }

  if (is.null(plot_options)) plot_options <- list()

  ## ------------------------------------------------------------------
  ## 1. SELECT DATA (IMPUTATION SUBSET)
  ## ------------------------------------------------------------------
  if (!imp_col %in% names(data)) {
    stop("imp_col not found in data")
  }

  if (MI_method == "first") {
    Data_Subset <- subset(data, get(imp_col) == 1)
  } else if (MI_method %in% c("average", "Rubin")) {
    Data_Subset <- data
  } else {
    stop("MI_method must be 'first', 'average', or 'Rubin'")
  }

  if (!variable_x %in% names(Data_Subset)) {
    stop("variable_x not found in data")
  }
  if (!outcome_var %in% names(Data_Subset)) {
    stop("outcome_var not found in data")
  }

  ## ------------------------------------------------------------------
  ## 2. GLOBAL KNOTS FOR THE SPLINE
  ## ------------------------------------------------------------------
  x_all <- Data_Subset[[variable_x]]
  if (!is.numeric(x_all)) {
    stop("variable_x must be numeric for restricted cubic spline")
  }

  if (!is.null(spline_knots_percentile)) {
    if (length(spline_knots_percentile) != spline_knots_n) {
      stop("length(spline_knots_percentile) must equal spline_knots_n")
    }
    probs_knots <- spline_knots_percentile / 100
  } else {
    probs_knots <- seq(0, 1, length.out = spline_knots_n)
  }

  global_knots <- as.numeric(stats::quantile(x_all, probs = probs_knots, na.rm = TRUE))
  knots_str    <- paste0("c(", paste(signif(global_knots, 6), collapse = ", "), ")")

  ## ------------------------------------------------------------------
  ## 3. SPLINE TERM & COVARIABLES (population model)
  ## ------------------------------------------------------------------
  spline_term <- paste0("rcs(", variable_x, ", ", knots_str, ")")
  if (!is.null(subgroups)) {
    spline_term <- paste0(spline_term, " * ", subgroups)
  }

  ## expand covariables with "*" into main effects and interactions
  expand_terms <- function(term) {
    if (grepl("^rcs\\(", term)) return(term)
    if (grepl("\\*", term)) {
      vars_split <- unlist(strsplit(term, "\\*"))
      vars_split <- trimws(vars_split)
      all_combinations <- lapply(seq_along(vars_split), function(k) {
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

    # helper: extract base variable name from rcs()/poly() if needed
    extract_var_from_term <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^poly\\(", term)) {
        inside <- sub("^[^(]+\\(([^,]+),.*$", "\\1", term)
        trimws(inside)
      } else {
        term
      }
    }

    all_base_cov <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
    all_base_cov <- trimws(all_base_cov)
    all_base_cov <- vapply(all_base_cov, extract_var_from_term, character(1))

    miss_cov <- setdiff(all_base_cov, names(Data_Subset))
    if (length(miss_cov) > 0) {
      stop("Covariables not found in data: ", paste(miss_cov, collapse = ", "))
    }

    covariables_str <- paste(expanded_covariables, collapse = " + ")
    covariables_str <- paste0(" + ", covariables_str)
  } else {
    covariables_str <- ""
  }

  ## trial factor for POPULATION model
  trial_str <- ""
  if (trial_factor == "Yes") {
    trial_str <- paste0(" + as.factor(", trial_col, ")")
  }

  ## offsets
  offset_str <- ""
  if (followup_offset == "Yes" && model_type %in% c("nb", "poisson")) {
    offset_str <- paste0(" + offset(log(", followup_col, "))")
  }

  ## random effects (population model only)
  rand_eff_str <- ""
  if (random_intercept == "Yes") {
    rand_eff_str <- paste0(rand_eff_str, " + (1 | ", random_intercept_var, ")")
  }
  if (random_slope == "Yes" && !is.null(predictor_vars_random_slope)) {
    for (v in predictor_vars_random_slope) {
      rand_eff_str <- paste0(rand_eff_str,
                             " + (0 + ", v, " | ", random_intercept_var, ")")
    }
    if (!is.null(covariables_random_slope)) {
      for (v in covariables_random_slope) {
        rand_eff_str <- paste0(rand_eff_str,
                               " + (0 + ", v, " | ", random_intercept_var, ")")
      }
    }
  }

  ## ------------------------------------------------------------------
  ## 4. BUILD POPULATION MODEL FORMULA
  ## ------------------------------------------------------------------
  if (model_type == "cox") {
    rhs <- paste0(spline_term, covariables_str, trial_str)
    formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ ", rhs)
  } else {
    rhs <- paste0(spline_term, covariables_str, trial_str, offset_str)
    formula_str <- paste0(outcome_var, " ~ ", rhs)
  }

  ## random effects only added for mixed models (nb/poisson/logistic/lm, Cox via coxme)
  mixed_formula_str <- formula_str
  if (random_intercept == "Yes" || random_slope == "Yes") {
    if (model_type != "cox") {
      mixed_formula_str <- paste0(formula_str, rand_eff_str)
    }
  }

  formula_obj       <- stats::as.formula(formula_str)
  mixed_formula_obj <- stats::as.formula(mixed_formula_str)

  ## ------------------------------------------------------------------
  ## 5. FIT POPULATION MODEL (GLM/GLMM + MI pooling)
  ## ------------------------------------------------------------------
  fit_single_model <- function(dat, use_mixed = FALSE) {
    fm <- if (use_mixed) mixed_formula_obj else formula_obj

    if (!use_mixed) {
      if (model_type == "nb") {
        requireNamespace("MASS", quietly = TRUE)
        MASS::glm.nb(fm, data = dat)
      } else if (model_type == "poisson") {
        glm(fm, family = poisson(link = "log"), data = dat)
      } else if (model_type == "logistic") {
        glm(fm, family = binomial(link = "logit"), data = dat)
      } else if (model_type == "lm") {
        glm(fm, family = gaussian(), data = dat)
      } else if (model_type == "cox") {
        requireNamespace("survival", quietly = TRUE)
        survival::coxph(fm, data = dat)
      }
    } else {
      ## mixed models
      if (model_type %in% c("nb", "poisson", "logistic", "lm")) {
        requireNamespace("lme4", quietly = TRUE)
        if (model_type == "nb") {
          ## approximate: Poisson glmer if glmmTMB not used here
          lme4::glmer(mixed_formula_obj, family = poisson(link = "log"),
                      data = dat)
        } else if (model_type == "poisson") {
          lme4::glmer(mixed_formula_obj, family = poisson(link = "log"),
                      data = dat)
        } else if (model_type == "logistic") {
          lme4::glmer(mixed_formula_obj, family = binomial(link = "logit"),
                      data = dat)
        } else if (model_type == "lm") {
          lme4::lmer(mixed_formula_obj, data = dat)
        }
      } else if (model_type == "cox") {
        requireNamespace("coxme", quietly = TRUE)
        coxme::coxme(mixed_formula_obj, data = dat)
      }
    }
  }

  ## --- main fit across imputations ---
  population_model <- NULL
  pooled_preds     <- NULL

  ## build prediction frame template (used later)
  x_seq <- seq(
    from = stats::quantile(x_all, prediction_range[1], na.rm = TRUE),
    to   = stats::quantile(x_all, prediction_range[2], na.rm = TRUE),
    length.out = 100
  )

  build_pred_frame <- function(base_data, x_vals) {
    if (is.null(subgroups)) {
      pd <- data.frame(x_vals)
      colnames(pd) <- variable_x
    } else {
      sub_lvls <- if (is.factor(base_data[[subgroups]])) {
        levels(base_data[[subgroups]])
      } else {
        sort(unique(base_data[[subgroups]]))
      }
      pd <- expand.grid(x_vals, sub_lvls, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      colnames(pd) <- c(variable_x, subgroups)
      if (is.factor(base_data[[subgroups]])) {
        pd[[subgroups]] <- factor(pd[[subgroups]], levels = levels(base_data[[subgroups]]))
      } else {
        pd[[subgroups]] <- factor(pd[[subgroups]])
      }
    }

    ## fill covariates with median/mode
    if (!is.null(expanded_covariables)) {
      extract_var_from_term <- function(term) {
        if (grepl("^rcs\\(", term) || grepl("^poly\\(", term)) {
          inside <- sub("^[^(]+\\(([^,]+),.*$", "\\1", term)
          trimws(inside)
        } else {
          term
        }
      }
      cov_vars <- unique(unlist(strsplit(gsub(":", "*", expanded_covariables), "\\*")))
      cov_vars <- trimws(cov_vars)
      cov_vars <- vapply(cov_vars, extract_var_from_term, character(1))
      for (cv in cov_vars) {
        if (cv %in% names(base_data)) {
          if (is.factor(base_data[[cv]])) {
            pd[[cv]] <- factor(names(which.max(table(base_data[[cv]]))),
                               levels = levels(base_data[[cv]]))
          } else {
            pd[[cv]] <- stats::median(base_data[[cv]], na.rm = TRUE)
          }
        }
      }
    }

    if (trial_factor == "Yes") {
      pd[[trial_col]] <- names(which.max(table(base_data[[trial_col]])))
    }
    if (followup_offset == "Yes") {
      pd[[followup_col]] <- 365
    }
    if (random_intercept == "Yes") {
      pd[[random_intercept_var]] <- names(which.max(table(base_data[[random_intercept_var]])))
    }

    pd
  }

  ## predict helper on link/responsescale + se
  predict_with_se <- function(mod, newdata) {
    if (inherits(mod, c("lmerMod", "glmerMod"))) {
      lp <- predict(mod, newdata = newdata, type = "link", se.fit = TRUE, re.form = NA)
      fit <- lp$fit
      se  <- lp$se.fit
    } else if (inherits(mod, "coxph")) {
      lp <- predict(mod, newdata = newdata, type = "lp", se.fit = TRUE)
      fit <- lp$fit
      se  <- lp$se.fit
    } else {
      lp <- predict(mod, newdata = newdata, type = "link", se.fit = TRUE)
      fit <- lp$fit
      se  <- lp$se.fit
    }
    list(fit = fit, se.fit = se)
  }

  ## link -> response transform
  link_to_response <- function(fit, se, type) {
    if (type %in% c("nb", "poisson", "cox")) {
      pred <- exp(fit)
      lower <- exp(fit - 1.96 * se)
      upper <- exp(fit + 1.96 * se)
    } else if (type == "logistic") {
      pred <- stats::plogis(fit)
      lower <- stats::plogis(fit - 1.96 * se)
      upper <- stats::plogis(fit + 1.96 * se)
    } else if (type == "lm") {
      pred <- fit
      lower <- fit - 1.96 * se
      upper <- fit + 1.96 * se
    } else {
      pred <- exp(fit)
      lower <- exp(fit - 1.96 * se)
      upper <- exp(fit + 1.96 * se)
    }
    list(pred = pred, lower = lower, upper = upper)
  }

  ## --- MI_method implementation ---
  if (MI_method == "first") {
    use_mixed <- (random_intercept == "Yes" | random_slope == "Yes")
    population_model <- fit_single_model(Data_Subset, use_mixed = use_mixed)
    pred_data <- build_pred_frame(Data_Subset, x_seq)

    if (model_type == "cox") {
      ## coxph/coxme: predict on lp, then exp
      if (inherits(population_model, "coxph")) {
        pr <- predict(population_model, newdata = pred_data, type = "lp", se.fit = TRUE)
        resp <- link_to_response(pr$fit, pr$se.fit, "cox")
      } else {
        ## coxme: no se.fit, approximate (no CI)
        lp <- predict(population_model, newdata = pred_data, type = "lp")
        resp <- list(pred = exp(lp),
                     lower = NA_real_,
                     upper = NA_real_)
      }
    } else {
      pr <- predict_with_se(population_model, pred_data)
      resp <- link_to_response(pr$fit, pr$se.fit, model_type)
    }

    pred_data$prediction <- resp$pred
    pred_data$lower_ci   <- resp$lower
    pred_data$upper_ci   <- resp$upper

  } else {
    ## average / Rubin
    imps <- sort(unique(Data_Subset[[imp_col]]))
    pred_data <- build_pred_frame(Data_Subset, x_seq)

    all_fits <- matrix(NA_real_, nrow = length(x_seq) *
                         ifelse(is.null(subgroups), 1, length(unique(pred_data[[subgroups]]))),
                       ncol = length(imps))
    all_ses  <- all_fits

    row_template <- NULL

    for (j in seq_along(imps)) {
      imp <- imps[j]
      dat_imp <- subset(Data_Subset, get(imp_col) == imp)
      use_mixed <- (random_intercept == "Yes" | random_slope == "Yes")
      mod_imp <- fit_single_model(dat_imp, use_mixed = use_mixed)
      pd_imp  <- build_pred_frame(dat_imp, x_seq)
      if (is.null(row_template)) {
        row_template <- nrow(pd_imp)
      }
      pr_imp <- predict_with_se(mod_imp, pd_imp)
      all_fits[1:nrow(pd_imp), j] <- pr_imp$fit
      all_ses[1:nrow(pd_imp),  j] <- pr_imp$se.fit
    }

    mean_fit <- rowMeans(all_fits, na.rm = TRUE)
    m        <- length(imps)

    mean_var <- numeric(length(mean_fit))
    if (MI_method == "Rubin") {
      for (i in seq_along(mean_fit)) {
        fits_i <- all_fits[i, ]
        ses_i  <- all_ses[i, ]
        fits_i <- fits_i[is.finite(fits_i)]
        ses_i  <- ses_i[is.finite(ses_i)]
        if (length(fits_i) == 0) {
          mean_var[i] <- NA_real_
        } else {
          W <- mean(ses_i^2)
          B <- stats::var(fits_i)
          mean_var[i] <- W + (1 + 1/length(fits_i)) * B
        }
      }
    } else { ## "average"
      for (i in seq_along(mean_fit)) {
        ses_i <- all_ses[i, ]
        ses_i <- ses_i[is.finite(ses_i)]
        if (length(ses_i) == 0) {
          mean_var[i] <- NA_real_
        } else {
          mean_var[i] <- mean(ses_i^2)
        }
      }
    }

    se_pool <- sqrt(mean_var)
    resp <- link_to_response(mean_fit, se_pool, model_type)

    pred_data$prediction <- resp$pred
    pred_data$lower_ci   <- resp$lower
    pred_data$upper_ci   <- resp$upper
  }

  ## ------------------------------------------------------------------
  ## 6. OPTIONAL DERIVATIVES (simple numeric, MI_method = "first")
  ## ------------------------------------------------------------------
  derivative_data <- NULL
  if (calculate_derivatives) {
    if (is.null(derivative_points)) {
      derivative_points <- x_seq
    }

    deriv_data <- build_pred_frame(Data_Subset, derivative_points)

    eps <- 1e-5 * stats::sd(Data_Subset[[variable_x]], na.rm = TRUE)
    data_plus  <- deriv_data
    data_minus <- deriv_data
    data_plus[[variable_x]]  <- deriv_data[[variable_x]] + eps
    data_minus[[variable_x]] <- deriv_data[[variable_x]] - eps

    ## predictions on response scale
    pred_resp <- function(mod, newdata) {
      if (model_type == "cox") {
        lp <- predict(mod, newdata = newdata, type = "lp")
        exp(lp)
      } else if (model_type %in% c("nb", "poisson")) {
        lp <- predict(mod, newdata = newdata, type = "link")
        exp(lp)
      } else if (model_type == "logistic") {
        lp <- predict(mod, newdata = newdata, type = "link")
        stats::plogis(lp)
      } else {
        predict(mod, newdata = newdata)
      }
    }

    f_plus  <- pred_resp(population_model, data_plus)
    f_minus <- pred_resp(population_model, data_minus)

    deriv <- (f_plus - f_minus) / (2 * eps)
    se_approx <- abs(deriv) * 0.2  # crude approximation
    lower_d <- deriv - 1.96 * se_approx
    upper_d <- deriv + 1.96 * se_approx

    deriv_data$derivative          <- deriv
    deriv_data$se_derivative       <- se_approx
    deriv_data$lower_ci_derivative <- lower_d
    deriv_data$upper_ci_derivative <- upper_d
    deriv_data$x_point             <- deriv_data[[variable_x]]

    derivative_data <- deriv_data
  }

  ## ------------------------------------------------------------------
  ## 7. PLOT (population curve)
  ## ------------------------------------------------------------------
  x_lab <- if (!is.null(plot_options$x_lab)) plot_options$x_lab else variable_x
  if (model_type %in% c("nb", "poisson")) {
    y_lab_default <- "Predicted rate"
  } else if (model_type == "logistic") {
    y_lab_default <- "Predicted probability"
  } else if (model_type == "lm") {
    y_lab_default <- "Predicted mean"
  } else {
    y_lab_default <- "Predicted hazard ratio"
  }
  y_lab <- if (!is.null(plot_options$y_lab)) plot_options$y_lab else y_lab_default

  plt <- ggplot2::ggplot(pred_data,
                         ggplot2::aes_string(x = variable_x, y = "prediction")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                         alpha = ifelse(is.null(plot_options$ribbon_alpha),
                                        0.3, plot_options$ribbon_alpha),
                         fill = ifelse(is.null(plot_options$fill_colors),
                                       "grey80", plot_options$fill_colors[1])) +
    ggplot2::geom_line(linewidth = ifelse(is.null(plot_options$line_size),
                                          1, plot_options$line_size),
                       color = ifelse(is.null(plot_options$colors),
                                      "black", plot_options$colors[1])) +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::theme_minimal()

  ## ------------------------------------------------------------------
  ## 8. group_fits: independent per-trial spline fits (no random effects)
  ## ------------------------------------------------------------------
  group_fits_df <- NULL

  if (group_fits) {
    if (is.null(trial_col)) {
      stop("group_fits = TRUE requires trial_col to be specified")
    }
    if (!trial_col %in% names(Data_Subset)) {
      stop("trial_col not found in data")
    }

    ## helper: identify non-constant covariates in this group
    non_constant_covs <- function(dat, vars) {
      keep <- logical(length(vars))
      for (i in seq_along(vars)) {
        v <- vars[i]
        if (!v %in% names(dat)) {
          keep[i] <- FALSE
        } else {
          x <- dat[[v]]
          if (is.factor(x)) {
            keep[i] <- length(unique(stats::na.omit(x))) > 1
          } else {
            keep[i] <- stats::sd(x, na.rm = TRUE) > 0
          }
        }
      }
      vars[keep]
    }

    ## base covariates (no interactions) for group fits
    base_cov_vars <- NULL
    if (!is.null(covariables)) {
      base_cov_vars <- unique(trimws(unlist(strsplit(gsub("\\*", ":", covariables), ":"))))
      base_cov_vars <- setdiff(base_cov_vars, variable_x)
    }

    all_groups <- sort(unique(Data_Subset[[trial_col]]))
    group_list <- vector("list", length(all_groups))

    for (gi in seq_along(all_groups)) {
      g  <- all_groups[gi]
      dg <- droplevels(subset(Data_Subset, get(trial_col) == g))
      if (nrow(dg) < spline_knots_n) next

      ## x-range for this group
      xg <- dg[[variable_x]]
      if (all(is.na(xg))) next

      xg_seq <- seq(
        from = stats::quantile(xg, prediction_range[1], na.rm = TRUE),
        to   = stats::quantile(xg, prediction_range[2], na.rm = TRUE),
        length.out = 100
      )

      ## covariates to use (non-constant in this group)
      group_covs <- character(0)
      if (!is.null(base_cov_vars)) {
        group_covs <- non_constant_covs(dg, base_cov_vars)
      }

      ## build formula for this group: outcome ~ rcs(x, global_knots) + group_covs + offset
      knots_str_g  <- knots_str  # reuse global knots
      spline_term_g <- paste0("rcs(", variable_x, ", ", knots_str_g, ")")

      rhs_g <- spline_term_g
      if (length(group_covs) > 0) {
        rhs_g <- paste(rhs_g, paste(group_covs, collapse = " + "), sep = " + ")
      }
      if (followup_offset == "Yes" && model_type %in% c("nb", "poisson")) {
        rhs_g <- paste0(rhs_g, " + offset(log(", followup_col, "))")
      }

      form_g <- stats::as.formula(paste0(outcome_var, " ~ ", rhs_g))

      ## fit simple GLM/NB in this group
      mod_g <- switch(model_type,
                      "nb" = { requireNamespace("MASS", quietly = TRUE);
                        MASS::glm.nb(form_g, data = dg) },
                      "poisson" = glm(form_g, family = poisson(link = "log"), data = dg),
                      "logistic" = glm(form_g, family = binomial(link = "logit"), data = dg),
                      "lm" = glm(form_g, family = gaussian(), data = dg),
                      "cox" = {
                        requireNamespace("survival", quietly = TRUE)
                        survival::coxph(form_g, data = dg)
                      })

      ## prediction frame for this group
      pg <- data.frame(xg_seq)
      colnames(pg) <- variable_x

      ## fill covariates with group medians/modes
      if (length(group_covs) > 0) {
        for (cv in group_covs) {
          if (cv %in% names(dg)) {
            if (is.factor(dg[[cv]])) {
              pg[[cv]] <- factor(names(which.max(table(dg[[cv]]))),
                                 levels = levels(dg[[cv]]))
            } else {
              pg[[cv]] <- stats::median(dg[[cv]], na.rm = TRUE)
            }
          }
        }
      }
      if (followup_offset == "Yes") {
        pg[[followup_col]] <- 365
      }

      ## predict with se
      if (model_type == "cox") {
        prg <- predict(mod_g, newdata = pg, type = "lp", se.fit = TRUE)
        resp_g <- link_to_response(prg$fit, prg$se.fit, "cox")
      } else {
        prg <- predict_with_se(mod_g, pg)
        resp_g <- link_to_response(prg$fit, prg$se.fit, model_type)
      }

      pg$prediction <- resp_g$pred
      pg$lower_ci   <- resp_g$lower
      pg$upper_ci   <- resp_g$upper
      pg[[trial_col]] <- g

      group_list[[gi]] <- pg
    }

    group_list <- group_list[!vapply(group_list, is.null, logical(1))]
    if (length(group_list) > 0) {
      group_fits_df <- do.call(rbind, group_list)
      rownames(group_fits_df) <- NULL
      names(group_fits_df)[names(group_fits_df) == trial_col] <- trial_col
    }
  }

  ## ------------------------------------------------------------------
  ## 9. RETURN
  ## ------------------------------------------------------------------
  res <- list(
    predictions = pred_data,
    model       = population_model,
    plot        = plt,
    plot_data   = pred_data,
    knots       = global_knots,
    spline_knots_n          = spline_knots_n,
    spline_knots_percentile = spline_knots_percentile,
    derivatives = derivative_data,
    group_fits  = group_fits_df
  )

  attr(res, "formula") <- format(formula_obj)
  class(res) <- c("MI_spline_result", class(res))
  res
}
