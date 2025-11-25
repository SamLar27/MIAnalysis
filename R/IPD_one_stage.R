#' IPD One-Stage Modelling on Multiply Imputed Data (with RCS support)
#'
#' Fits models to multiply imputed datasets and pools results using Rubin's Rules.
#' Supports GLM, mixed models (via lme4 / glmmTMB / coxme), interactions, and
#' restricted cubic splines (rcs) via the `rms` package. Designed for IPD
#' one-stage analyses with optional random effects / stratified intercepts.
#'
#' Intercept / slope structure is inferred automatically:
#'   - If `stratified_intercept_var` is provided:
#'       * GLM: adds `+ as.factor(stratified_intercept_var)` (keeps regular intercept)
#'       * Cox: adds `+ strata(stratified_intercept_var)`
#'   - If `random_intercept_var` is provided:
#'       * Mixed model is used (glmmTMB/lme4/coxme)
#'       * Random intercept and/or random slopes are defined by
#'         `predictor_vars_random_slope` and `covariables_random_slope`.
#'         Special case when `random_intercept_var == stratified_intercept_var`:
#'           - and no random slopes      -> no random effect (pure stratified intercept)
#'           - and at least one slope    -> random slopes only:
#'                (0 + slopes | random_intercept_var)
#'
#' Restricted cubic splines:
#'   - Use `spline_terms` (character vector of variable names).
#'   - Knots are defined ONCE on the full stacked data via:
#'       * `spline_knots_percentile`: numeric vector of percentiles (e.g. c(5,50,95)),
#'         applied to the whole distribution of each variable, OR
#'       * `spline_knots_n`: number of knots; default percentiles are spread between 5 and 95.
#'   - Requires package `rms`; internally uses `rms::rcs(x, parms = knots)`.
#'   - Spline coefficients in the output are renamed as:
#'       `<var>_rcs_linear`, `<var>_rcs_nl1`, `<var>_rcs_nl2`, ...
#'
#' Random slopes + splines:
#'   - If a variable appears in `spline_terms` and also in
#'     `predictor_vars_random_slope` or `covariables_random_slope`, the raw
#'     variable name is automatically expanded to its spline basis terms for
#'     the random-effects part, e.g. random slopes on
#'     `FEV1_preBD_PCT_0W` become random slopes on:
#'       `FEV1_preBD_PCT_0W_rcs_linear` and `FEV1_preBD_PCT_0W_rcs_nl*`.
#'
#' @param data Data frame with all imputations stacked.
#' @param outcome_var Dependent variable (for non-Cox models).
#' @param predictor_vars Character vector of predictors to be included as fixed effects.
#' @param covariables Character vector of covariables (always included as fixed effects).
#' @param imp_col Column name indicating imputation index (default ".imp").
#' @param imp_n Number of imputations (if NULL, detected automatically).
#' @param model_type One of: "nb","lm","bin","poisson","gamma",
#'   "quasipoisson","quasibinomial","cox".
#' @param followup_offset "Yes"/"No" – whether to include offset(log(followup_col)).
#' @param followup_col Name of follow-up duration column if offset used.
#' @param trial_factor "Yes"/"No" – include trial as additional fixed factor.
#' @param trial_col Trial column (if trial_factor = "Yes").
#' @param time_col Time variable for Cox models.
#' @param event_col Event variable for Cox models.
#' @param formula_string Optional custom formula (overrides automatic building).
#' @param highlight_interactions Logical, add flags for interactions/splines.
#' @param spline_terms Character vector of variable names for rcs() terms.
#' @param spline_knots_n Number of knots (if percentiles are not given).
#' @param spline_knots_percentile Numeric vector of percentiles for knots (0–100).
#' @param random_intercept_var Grouping variable for random effects (NULL = no random effects).
#' @param predictor_vars_random_slope Character vector of predictors with random slope.
#' @param covariables_random_slope Character vector of covariables with random slope.
#' @param stratified_intercept_var Variable defining stratified intercept (GLM) / baseline hazard (Cox).
#' @param model_performance Logical; if TRUE, compute model performance via `performance::model_performance()`.
#' @param weighted_intercept Logical; if TRUE and `stratified_intercept_var` is not NULL,
#'   computes a single weighted intercept across strata.
#'
#' @return
#' A list of class "IPD_one_stage" with elements:
#'   - table              : pooled fixed-effect estimates (without stratified intercept rows)
#'   - term               : same as table$term
#'   - Intercept          : intercept + stratified intercept rows
#'   - Weighted_intercept : one-row data.frame with weighted intercept (if requested)
#' and attributes:
#'   - "fit_log", "models", "performance_per_imp", "model_type", "formula",
#'     "imputations", "n_imp", "has_random_effects", "random_intercept_var",
#'     "stratified_intercept_var", "spline_info".
#'
#' @export
IPD_one_stage <- function(data,
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
                          highlight_interactions = TRUE,
                          spline_terms = NULL,
                          spline_knots_n = NULL,
                          spline_knots_percentile = NULL,
                          random_intercept_var = NULL,
                          predictor_vars_random_slope = NULL,
                          covariables_random_slope = NULL,
                          stratified_intercept_var = NULL,
                          model_performance = FALSE,
                          weighted_intercept = FALSE) {

  ## ------------------------------------------------------------
  ## Basic checks
  ## ------------------------------------------------------------
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  if (!imp_col %in% names(data)) stop("imp_col not found in data.")
  if (model_type != "cox" && !outcome_var %in% names(data)) {
    stop("outcome_var not found in data.")
  }
  if (!followup_offset %in% c("Yes", "No")) {
    stop("followup_offset must be 'Yes' or 'No'.")
  }
  if (followup_offset == "Yes") {
    if (is.null(followup_col) || !followup_col %in% names(data)) {
      stop("followup_offset = 'Yes' but followup_col is missing or not in data.")
    }
    if (any(data[[followup_col]] <= 0, na.rm = TRUE)) {
      stop("followup_col must be strictly positive for offset(log(followup_col)).")
    }
  }
  if (!trial_factor %in% c("Yes", "No")) {
    stop("trial_factor must be 'Yes' or 'No'.")
  }
  if (trial_factor == "Yes" && (is.null(trial_col) || !trial_col %in% names(data))) {
    stop("trial_factor = 'Yes' but trial_col is missing or not in data.")
  }
  if (!is.null(random_intercept_var) && !random_intercept_var %in% names(data)) {
    stop("random_intercept_var not found in data.")
  }
  if (!is.null(stratified_intercept_var) && !stratified_intercept_var %in% names(data)) {
    stop("stratified_intercept_var not found in data.")
  }
  if (model_type == "cox") {
    if (is.null(time_col) || is.null(event_col)) {
      stop("For Cox models, time_col and event_col must be provided.")
    }
    if (!time_col %in% names(data))  stop("time_col not found in data.")
    if (!event_col %in% names(data)) stop("event_col not found in data.")
  }

  # Random slope checks: allow spline_terms to appear even if basis not yet built
  if (!is.null(predictor_vars_random_slope)) {
    miss <- setdiff(predictor_vars_random_slope,
                    c(names(data), spline_terms %||% character(0)))
    if (length(miss) > 0) {
      stop("predictor_vars_random_slope not found in data (and not declared as spline_terms): ",
           paste(miss, collapse = ", "))
    }
  }
  if (!is.null(covariables_random_slope)) {
    miss <- setdiff(covariables_random_slope, names(data))
    if (length(miss) > 0) {
      stop("covariables_random_slope not found in data: ",
           paste(miss, collapse = ", "))
    }
  }
  use_random_effects <- !is.null(random_intercept_var)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  ## ------------------------------------------------------------
  ## 1) Prepare imputation list
  ## ------------------------------------------------------------
  actual_imps <- sort(unique(data[[imp_col]]))
  if (is.null(imp_n)) {
    imp_n <- length(actual_imps)
  } else if (imp_n != length(actual_imps)) {
    warning("imp_n does not match number of unique imputations in data; using present imputations only.")
    imp_n <- length(actual_imps)
  }

  implist <- vector("list", length(actual_imps))
  for (i in seq_along(actual_imps)) {
    imp_val <- actual_imps[i]
    dat_i   <- data[data[[imp_col]] == imp_val, , drop = FALSE]
    if (nrow(dat_i) == 0) stop("No data for imputation ", imp_val)
    implist[[i]] <- dat_i
  }

  ## ------------------------------------------------------------
  ## 2) Global spline info (knots on full stacked data)
  ## ------------------------------------------------------------
  spline_info <- NULL
  if (!is.null(spline_terms)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("spline_terms specified, but package 'rms' is not installed.")
    }

    spline_terms <- intersect(spline_terms, names(data))
    if (length(spline_terms) == 0L) {
      stop("spline_terms not found in data.")
    }

    spline_info <- list()
    for (var in spline_terms) {
      x_all <- data[[var]]

      if (!is.null(spline_knots_percentile)) {
        probs  <- spline_knots_percentile / 100
        knots  <- as.numeric(stats::quantile(x_all, probs = probs, na.rm = TRUE))
        knots  <- sort(unique(knots))
      } else if (!is.null(spline_knots_n)) {
        k      <- spline_knots_n
        probs  <- seq(5, 95, length.out = k) / 100
        knots  <- as.numeric(stats::quantile(x_all, probs = probs, na.rm = TRUE))
        knots  <- sort(unique(knots))
      } else {
        # default: 4 knots at 5,35,65,95
        probs  <- c(5, 35, 65, 95) / 100
        knots  <- as.numeric(stats::quantile(x_all, probs = probs, na.rm = TRUE))
        knots  <- sort(unique(knots))
      }

      if (length(knots) < 3L) {
        stop("Not enough distinct knot locations for spline term '", var, "'.")
      }

      # For rcs with K knots, df = K-1 basis functions
      K          <- length(knots)
      n_basis    <- K - 1L
      basis_names <- c(
        paste0(var, "_rcs_linear"),
        if (n_basis > 1L) paste0(var, "_rcs_nl", seq_len(n_basis - 1L)) else character(0)
      )

      spline_info[[var]] <- list(
        var         = var,
        knots       = knots,
        basis_names = basis_names
      )
    }
  }

  ## ------------------------------------------------------------
  ## 3) Expand rhs_main_terms with spline basis names
  ## ------------------------------------------------------------
  rhs_main_terms <- c(predictor_vars, covariables)

  if (!is.null(spline_info)) {
    for (var in names(spline_info)) {
      if (var %in% rhs_main_terms) {
        rhs_main_terms <- setdiff(rhs_main_terms, var)
        rhs_main_terms <- c(rhs_main_terms, spline_info[[var]]$basis_names)
      }
    }
  }

  ## ------------------------------------------------------------
  ## 4) Expand random slopes w.r.t. spline basis
  ## ------------------------------------------------------------
  expand_for_spline <- function(vars_vec) {
    if (is.null(vars_vec)) return(NULL)
    v <- vars_vec
    if (!is.null(spline_info)) {
      for (var in names(spline_info)) {
        if (var %in% v) {
          v <- union(setdiff(v, var), spline_info[[var]]$basis_names)
        }
      }
    }
    unique(v)
  }

  predictor_vars_random_slope_exp   <- expand_for_spline(predictor_vars_random_slope)
  covariables_random_slope_exp      <- expand_for_spline(covariables_random_slope)

  ## ------------------------------------------------------------
  ## 5) Build formula (fixed + random parts)
  ## ------------------------------------------------------------
  trial_term <- if (trial_factor == "Yes") {
    paste0("+ as.factor(", trial_col, ")")
  } else {
    ""
  }

  strat_glm_term <- ""
  strat_cox_term <- ""
  if (!is.null(stratified_intercept_var)) {
    if (model_type == "cox") {
      strat_cox_term <- paste0("+ strata(", stratified_intercept_var, ")")
    } else {
      strat_glm_term <- paste0("+ as.factor(", stratified_intercept_var, ")")
    }
  }

  offset_term <- if (followup_offset == "Yes") {
    paste0("+ offset(log(", followup_col, "))")
  } else {
    ""
  }

  build_random_term <- function() {
    if (!use_random_effects) return("")

    rs_vars <- character(0)
    if (!is.null(predictor_vars_random_slope_exp)) {
      rs_vars <- c(rs_vars, predictor_vars_random_slope_exp)
    }
    if (!is.null(covariables_random_slope_exp)) {
      rs_vars <- c(rs_vars, covariables_random_slope_exp)
    }
    rs_vars <- unique(rs_vars)

    same_group_as_strata <- (!is.null(stratified_intercept_var) &&
                               stratified_intercept_var == random_intercept_var)

    if (length(rs_vars) == 0) {
      if (same_group_as_strata) {
        # Only stratified intercept, no random effect
        return("")
      } else {
        return(paste0("+ (1 | ", random_intercept_var, ")"))
      }
    } else {
      rs_str <- paste(rs_vars, collapse = " + ")
      if (same_group_as_strata) {
        # random slopes only, no random intercept (0 + ...)
        return(paste0("+ (0 + ", rs_str, " | ", random_intercept_var, ")"))
      } else {
        # random intercept + random slopes
        return(paste0("+ (1 + ", rs_str, " | ", random_intercept_var, ")"))
      }
    }
  }

  random_term <- build_random_term()

  if (is.null(formula_string)) {
    rhs_core <- paste(rhs_main_terms, collapse = " + ")
    if (rhs_core == "") rhs_core <- "1"

    if (model_type == "cox") {
      formula_str <- paste0(
        "survival::Surv(", time_col, ",", event_col, ") ~ ",
        rhs_core, " ",
        strat_cox_term, " ",
        trial_term, " ",
        random_term
      )
    } else {
      formula_str <- paste0(
        outcome_var, " ~ ",
        rhs_core, " ",
        offset_term, " ",
        strat_glm_term, " ",
        trial_term, " ",
        random_term
      )
    }
  } else {
    formula_str <- formula_string
    if (model_type == "cox") {
      if (!grepl("^Surv\\(", formula_str)) {
        formula_str <- paste0(
          "survival::Surv(", time_col, ",", event_col, ") ~ ",
          formula_str
        )
      }
      if (!identical(trial_term, "") &&
          !grepl(paste0("as\\.factor\\(", trial_col, "\\)"), formula_str)) {
        formula_str <- paste(formula_str, trial_term)
      }
      if (!identical(strat_cox_term, "") &&
          !grepl(paste0("strata\\(", stratified_intercept_var, "\\)"), formula_str)) {
        formula_str <- paste(formula_str, strat_cox_term)
      }
      if (use_random_effects && !grepl("\\|", formula_str)) {
        formula_str <- paste(formula_str, random_term)
      }
    } else {
      if (!grepl(paste0("^", outcome_var, "\\s*~"), formula_str)) {
        formula_str <- paste(outcome_var, "~", formula_str)
      }
      if (followup_offset == "Yes" &&
          !grepl("offset\\(log\\(", formula_str)) {
        formula_str <- paste(formula_str, offset_term)
      }
      if (!identical(strat_glm_term, "") &&
          !grepl(paste0("as\\.factor\\(", stratified_intercept_var, "\\)"), formula_str)) {
        formula_str <- paste(formula_str, strat_glm_term)
      }
      if (!identical(trial_term, "") &&
          !grepl(paste0("as\\.factor\\(", trial_col, "\\)"), formula_str)) {
        formula_str <- paste(formula_str, trial_term)
      }
      if (use_random_effects && !grepl("\\|", formula_str)) {
        formula_str <- paste(formula_str, random_term)
      }
    }
  }

  model_formula <- stats::as.formula(formula_str)

  ## ------------------------------------------------------------
  ## 6) Fit per imputation (adding spline basis columns inside each dat_i)
  ## ------------------------------------------------------------
  models_list <- vector("list", length(actual_imps))
  names(models_list) <- paste0("imp_", actual_imps)

  fit_log <- data.frame(
    imp   = actual_imps,
    ok    = FALSE,
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(actual_imps)) {
    dat_i   <- implist[[i]]
    imp_val <- actual_imps[i]

    # add spline basis columns (same knots for all imputations)
    if (!is.null(spline_info)) {
      for (var in names(spline_info)) {
        info  <- spline_info[[var]]
        x     <- dat_i[[var]]
        # rcs with fixed knots
        basis <- rms::rcs(x, parms = info$knots)
        basis <- as.matrix(basis)
        if (ncol(basis) != length(info$basis_names)) {
          stop("Mismatch between expected and actual number of spline basis columns for '", var, "'.")
        }
        colnames(basis) <- info$basis_names
        dat_i[, info$basis_names] <- basis
      }
    }

    fit <- tryCatch({
      if (use_random_effects) {
        if (model_type == "lm") {
          if (!requireNamespace("lme4", quietly = TRUE)) {
            stop("Package 'lme4' is required for lm with random effects.")
          }
          lme4::lmer(model_formula, data = dat_i)
        } else if (model_type %in% c("bin")) {
          if (!requireNamespace("lme4", quietly = TRUE)) {
            stop("Package 'lme4' is required for binomial mixed models.")
          }
          lme4::glmer(model_formula, family = stats::binomial(), data = dat_i)
        } else if (model_type %in% c("poisson", "quasipoisson")) {
          if (!requireNamespace("lme4", quietly = TRUE)) {
            stop("Package 'lme4' is required for Poisson mixed models.")
          }
          lme4::glmer(model_formula, family = stats::poisson(), data = dat_i)
        } else if (model_type == "nb") {
          if (requireNamespace("glmmTMB", quietly = TRUE)) {
            glmmTMB::glmmTMB(model_formula,
                             family = glmmTMB::nbinom2,
                             data   = dat_i)
          } else {
            warning("glmmTMB not installed; using Poisson GLMM (approximation to NB).")
            if (!requireNamespace("lme4", quietly = TRUE)) {
              stop("Package 'lme4' is required for Poisson mixed models.")
            }
            lme4::glmer(model_formula, family = stats::poisson(), data = dat_i)
          }
        } else if (model_type == "cox") {
          if (!requireNamespace("coxme", quietly = TRUE)) {
            stop("Package 'coxme' is required for Cox models with random effects.")
          }
          coxme::coxme(model_formula, data = dat_i)
        } else {
          stop("Unsupported model_type with random effects.")
        }
      } else {
        switch(model_type,
               "nb"            = MASS::glm.nb(model_formula, data = dat_i),
               "lm"            = stats::glm(model_formula, family = stats::gaussian(), data = dat_i),
               "bin"           = stats::glm(model_formula, family = stats::binomial(), data = dat_i),
               "poisson"       = stats::glm(model_formula, family = stats::poisson(), data = dat_i),
               "gamma"         = stats::glm(model_formula, family = stats::Gamma(), data = dat_i),
               "quasipoisson"  = stats::glm(model_formula, family = stats::quasipoisson(), data = dat_i),
               "quasibinomial" = stats::glm(model_formula, family = stats::quasibinomial(), data = dat_i),
               "cox"           = survival::coxph(model_formula, data = dat_i),
               stop("Unsupported model_type.")
        )
      }
    },
    error = function(e) {
      fit_log$error[fit_log$imp == imp_val] <<- e$message
      NULL
    })

    if (!is.null(fit)) {
      fit_log$ok[fit_log$imp == imp_val] <- TRUE
    }

    models_list[[i]] <- fit
  }

  ok_idx    <- which(fit_log$ok)
  ok_models <- models_list[ok_idx]

  if (length(ok_models) == 0) {
    stop("All models failed to fit in IPD_one_stage.")
  }

  ## ------------------------------------------------------------
  ## 7) Pool coefficients across imputations
  ## ------------------------------------------------------------
  if (use_random_effects) {

    coefs_list <- lapply(ok_models, function(m) {
      if (is.null(m)) return(NULL)

      if (model_type == "cox") {
        c_est  <- coxme::fixef(m)
        c_vcov <- stats::vcov(m)
        c_se   <- sqrt(diag(c_vcov))
      } else if (inherits(m, "glmmTMB")) {
        c_est <- glmmTMB::fixef(m)$cond
        v     <- stats::vcov(m)$cond
        c_se  <- sqrt(diag(as.matrix(v)))
      } else if (inherits(m, "lmerMod") || inherits(m, "glmerMod")) {
        c_est <- lme4::fixef(m)
        c_se  <- sqrt(diag(as.matrix(stats::vcov(m))))
      } else {
        c_est <- stats::coef(m)
        c_se  <- sqrt(diag(as.matrix(stats::vcov(m))))
      }

      data.frame(
        term      = names(c_est),
        estimate  = as.numeric(c_est),
        std.error = as.numeric(c_se),
        stringsAsFactors = FALSE
      )
    })

    coefs_list <- Filter(Negate(is.null), coefs_list)
    if (length(coefs_list) == 0) stop("No valid coefficients for pooling.")

    all_terms <- unique(unlist(lapply(coefs_list, function(df) df$term)))

    pooled_results <- data.frame(
      term      = all_terms,
      estimate  = NA_real_,
      std.error = NA_real_,
      stringsAsFactors = FALSE
    )

    for (tt in all_terms) {
      est_vec <- sapply(coefs_list, function(df) {
        if (tt %in% df$term) df$estimate[df$term == tt] else NA_real_
      })
      se_vec <- sapply(coefs_list, function(df) {
        if (tt %in% df$term) df$std.error[df$term == tt] else NA_real_
      })

      est_vec <- est_vec[!is.na(est_vec)]
      se_vec  <- se_vec[!is.na(se_vec)]

      if (length(est_vec) < 1 || length(se_vec) < 1) next

      Q_bar <- mean(est_vec)
      U_bar <- mean(se_vec^2)
      if (length(est_vec) > 1) {
        B <- stats::var(est_vec)
      } else {
        B <- 0
      }
      T_var <- U_bar + (1 + 1 / length(est_vec)) * B
      pooled_results$estimate[pooled_results$term == tt]  <- Q_bar
      pooled_results$std.error[pooled_results$term == tt] <- sqrt(T_var)
    }

    Results <- pooled_results
    Results$`2.5 %`  <- Results$estimate - 1.96 * Results$std.error
    Results$`97.5 %` <- Results$estimate + 1.96 * Results$std.error
    Results$p.value  <- 2 * stats::pnorm(-abs(Results$estimate / Results$std.error))

  } else {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package 'mice' is required to pool non-mixed models.")
    }
    ok_models_clean <- Filter(Negate(is.null), models_list)
    pooled <- mice::pool(ok_models_clean)
    Results <- summary(pooled, conf.int = TRUE, exponentiate = FALSE)
    if (!"term" %in% names(Results)) {
      Results$term <- rownames(Results)
    }
    Results <- Results[, c("term", "estimate", "std.error", "2.5 %", "97.5 %", "p.value")]
  }

  ## ------------------------------------------------------------
  ## 8) Exponentiate
  ## ------------------------------------------------------------
  Results$exp_estimate   <- exp(Results$estimate)
  Results$exp_CI95_lower <- exp(Results$`2.5 %`)
  Results$exp_CI95_upper <- exp(Results$`97.5 %`)

  ## ------------------------------------------------------------
  ## 9) Flags
  ## ------------------------------------------------------------
  if (highlight_interactions) {
    Results$is_interaction <- grepl(":", Results$term)
  } else {
    Results$is_interaction <- FALSE
  }
  # spline terms identified by suffix convention "_rcs_"
  Results$is_spline     <- grepl("_rcs_", Results$term, fixed = TRUE)
  Results$is_polynomial <- FALSE  # not implemented here

  ## ------------------------------------------------------------
  ## 10) Split intercept vs other terms (including stratified intercepts)
  ## ------------------------------------------------------------
  pattern_strata <- if (!is.null(stratified_intercept_var)) {
    paste0("^as\\.factor\\(", stratified_intercept_var, "\\)")
  } else {
    NULL
  }

  is_intercept_row <- Results$term == "(Intercept)"
  is_strata_row <- if (!is.null(pattern_strata)) {
    grepl(pattern_strata, Results$term)
  } else {
    rep(FALSE, nrow(Results))
  }

  Intercept_table <- Results[is_intercept_row | is_strata_row, , drop = FALSE]
  Table_no_int    <- Results[!(is_intercept_row | is_strata_row), , drop = FALSE]

  ## ------------------------------------------------------------
  ## 11) Weighted intercept (optional)
  ## ------------------------------------------------------------
  Weighted_intercept <- NULL

  if (isTRUE(weighted_intercept) && !is.null(stratified_intercept_var)) {
    if (nrow(Intercept_table) > 0 && nrow(data) > 0) {
      strat_fac <- as.factor(data[[stratified_intercept_var]])
      strat_fac <- droplevels(strat_fac)

      N_by_stratum <- as.data.frame(table(strat_fac), stringsAsFactors = FALSE)
      colnames(N_by_stratum) <- c("level", "N")
      N_by_stratum$w <- N_by_stratum$N / sum(N_by_stratum$N)

      beta0 <- Intercept_table$estimate[Intercept_table$term == "(Intercept)"]
      if (length(beta0) == 0L) beta0 <- 0

      if (!is.null(pattern_strata)) {
        strata_coefs <- Results[grepl(pattern_strata, Results$term),
                                c("term", "estimate"), drop = FALSE]
        if (nrow(strata_coefs) > 0L) {
          strata_coefs$level <- sub(pattern_strata, "", strata_coefs$term)
        } else {
          strata_coefs <- data.frame(level = character(0),
                                     estimate = numeric(0))
        }
      } else {
        strata_coefs <- data.frame(level = character(0),
                                   estimate = numeric(0))
      }

      merged <- merge(N_by_stratum, strata_coefs,
                      by = "level", all.x = TRUE)
      merged$estimate[is.na(merged$estimate)] <- 0
      merged$beta_level <- beta0 + merged$estimate

      beta_w <- sum(merged$w * merged$beta_level)

      Weighted_intercept <- data.frame(
        term           = paste0("(Intercept)_weighted_", stratified_intercept_var),
        estimate       = beta_w,
        std.error      = NA_real_,
        `2.5 %`        = NA_real_,
        `97.5 %`       = NA_real_,
        p.value        = NA_real_,
        exp_estimate   = exp(beta_w),
        exp_CI95_lower = NA_real_,
        exp_CI95_upper = NA_real_,
        is_interaction = FALSE,
        is_spline      = FALSE,
        is_polynomial  = FALSE,
        stringsAsFactors = FALSE,
        check.names      = FALSE
      )
    }
  }

  ## ------------------------------------------------------------
  ## 12) Performance metrics per imputation (optional)
  ## ------------------------------------------------------------
  performance_per_imp <- NULL
  if (isTRUE(model_performance)) {
    if (!requireNamespace("performance", quietly = TRUE)) {
      warning("model_performance = TRUE, but package 'performance' is not installed. Skipping performance metrics.")
    } else {
      performance_per_imp <- lapply(seq_along(models_list), function(i) {
        m <- models_list[[i]]
        if (is.null(m)) return(NULL)
        tryCatch(
          performance::model_performance(m),
          error = function(e) NULL
        )
      })
      names(performance_per_imp) <- paste0("imp_", actual_imps)
    }
  }

  ## ------------------------------------------------------------
  ## 13) Return object
  ## ------------------------------------------------------------
  out <- list(
    table              = Table_no_int,
    term               = Table_no_int$term,
    Intercept          = Intercept_table,
    Weighted_intercept = Weighted_intercept
  )

  class(out) <- c("IPD_one_stage", "list")

  attr(out, "fit_log")             <- fit_log
  attr(out, "models")              <- models_list
  attr(out, "performance_per_imp") <- performance_per_imp
  attr(out, "model_type")          <- model_type
  attr(out, "formula")             <- formula_str
  attr(out, "imputations")         <- actual_imps
  attr(out, "n_imp")               <- length(actual_imps)
  attr(out, "has_random_effects")  <- use_random_effects
  attr(out, "random_intercept_var")    <- random_intercept_var
  attr(out, "stratified_intercept_var")<- stratified_intercept_var
  attr(out, "spline_info")         <- spline_info  # contains knots + basis names

  out
}
