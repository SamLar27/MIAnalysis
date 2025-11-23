#' IPD One-Stage Modelling on Multiply Imputed Data
#'
#' Fits models to multiply imputed datasets and pools results using Rubin's Rules.
#' Supports GLM, mixed models (via lme4 / glmmTMB / coxme), interactions,
#' restricted cubic splines (rcs) and polynomial terms. Designed for IPD
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
#'   - Knots are defined via:
#'       * `spline_knots_percentile`: numeric vector of percentiles (e.g. c(10,35,65,90)),
#'         applied to the whole distribution of each variable, OR
#'       * `spline_knots_n`: number of knots; default percentiles are spread between 5 and 95.
#'   - If `rms` is installed, uses `rms::rcs()`; otherwise falls back to `splines::bs()`.
#'   - The spline coefficients in the output are renamed as:
#'       `<var>_rcs_linear`, `<var>_rcs_nl1`, `<var>_rcs_nl2`, ...
#'
#' Polynomial terms:
#'   - Use `poly_terms` and `poly_degree` (2 or 3).
#'   - Internally uses `poly(var, degree, raw = TRUE)`.
#'
#' Performance indices:
#'   - Controlled by `performance_index`, which can be:
#'       * NULL (no performance metrics computed)
#'       * Any subset of:
#'         `c("Log_lik","AIC","AICc","BIC","BICc","RMSE","C_index")`
#'   - Metrics are computed per imputation using base R:
#'       * `Log_lik` : `logLik()`
#'       * `AIC`     : `AIC()`
#'       * `AICc`    : small-sample correction from AIC, using `nobs()` and df
#'       * `BIC`     : `BIC()`
#'       * `BICc`    : simple corrected BIC (approximate)
#'       * `RMSE`    : sqrt(mean((y - fitted)^2)), for non-Cox models
#'       * `C_index` : `survival::concordance()` (C-index) for Cox models (when supported)
#'
#' @param data Data frame with all imputations stacked.
#' @param outcome_var Dependent variable (for non-Cox models).
#' @param predictor_vars Character vector of predictors to be tested one by one.
#' @param covariables Character vector of covariables (always included as fixed effects).
#' @param imp_col Column name indicating imputation index (default ".imp").
#' @param imp_n Number of imputations (if NULL, detected automatically).
#' @param model_type One of: "nb","lm","bin","poisson","gamma",
#'   "quasipoisson","quasibinomial","cox".
#' @param followup_offset "Yes"/"No" – whether to include offset(log(followup_col)).
#' @param followup_col Name of follow-up duration column if offset used.
#' @param trial_factor "Yes"/"No" – include trial as fixed factor.
#' @param trial_col Trial column (if trial_factor = "Yes").
#' @param time_col Time variable for Cox models.
#' @param event_col Event variable for Cox models.
#' @param formula_string Optional custom formula (overrides automatic building).
#' @param highlight_interactions Logical, add flags for interactions/splines/polynomials.
#' @param spline_terms Character vector of variable names for rcs() / bs() terms.
#' @param spline_knots_n Number of knots (if percentiles are not given).
#' @param spline_knots_percentile Numeric vector of percentiles for knots (0–100).
#' @param poly_terms Character vector of variables with polynomial terms.
#' @param poly_degree Scalar or vector (2 or 3) giving polynomial degree(s).
#' @param include_poly_terms Logical, keep polynomial basis terms (default TRUE).
#' @param random_intercept_var Grouping variable for random effects (NULL = no random effects).
#' @param predictor_vars_random_slope Character vector of predictors with random slope.
#' @param covariables_random_slope Character vector of covariables with random slope.
#' @param stratified_intercept_var Variable defining stratified intercept (GLM) / baseline hazard (Cox).
#' @param performance_index Either NULL (no performance metrics) or a character vector
#'   with any subset of:
#'   `c("Log_lik","AIC","AICc","BIC","BICc","RMSE","C_index")`.
#'
#' @return
#'  If one predictor: a data.frame of pooled coefficients with attributes:
#'    - "models"              : list of fitted models (one per imputation)
#'    - "performance_per_imp" : list of lists with requested performance indices
#'    - "fit_log"             : data.frame with columns `imp`, `ok`, `error` (one row per imputation)
#'    - plus other metadata flags.
#'
#'  If multiple predictors: list of such data.frames with class "IPD_one_stage_multi"
#'    and attribute "combined_results".
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

  if (!is.null(predictor_vars_random_slope)) {
    miss <- setdiff(predictor_vars_random_slope, names(data))
    if (length(miss) > 0) {
      stop("predictor_vars_random_slope not found in data: ",
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

  ## ------------------------------------------------------------
  ## Imputation list
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
  ## Build formula
  ## ------------------------------------------------------------
  rhs_main_terms <- c(predictor_vars, covariables)

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

    if (!is.null(predictor_vars_random_slope)) {
      rs_vars <- c(rs_vars, predictor_vars_random_slope)
    }
    if (!is.null(covariables_random_slope)) {
      rs_vars <- c(rs_vars, covariables_random_slope)
    }
    rs_vars <- unique(rs_vars)

    same_group_as_strata <- (!is.null(stratified_intercept_var) &&
                               stratified_intercept_var == random_intercept_var)

    if (length(rs_vars) == 0) {
      if (same_group_as_strata) {
        return("")
      } else {
        return(paste0("+ (1 | ", random_intercept_var, ")"))
      }
    } else {
      rs_str <- paste(rs_vars, collapse = " + ")
      if (same_group_as_strata) {
        return(paste0("+ (0 + ", rs_str, " | ", random_intercept_var, ")"))
      } else {
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
  ## Fit per imputation
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
    stop("All models failed to fit in IPD_one_stage (random-effects block).")
  }

  ## ------------------------------------------------------------
  ## Pool coefficients across imputations
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
        v     <- stats::vcov(m)$cond   # <- FIX HERE
        c_se  <- sqrt(diag(v))
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
  ## Exponentiate
  ## ------------------------------------------------------------
  Results$exp_estimate   <- exp(Results$estimate)
  Results$exp_CI95_lower <- exp(Results$`2.5 %`)
  Results$exp_CI95_upper <- exp(Results$`97.5 %`)

  ## ------------------------------------------------------------
  ## Flags
  ## ------------------------------------------------------------
  if (highlight_interactions) {
    Results$is_interaction <- grepl(":", Results$term)
  } else {
    Results$is_interaction <- FALSE
  }
  Results$is_spline     <- FALSE
  Results$is_polynomial <- grepl("poly\\(", Results$term)

  ## ------------------------------------------------------------
  ## Split intercept vs other terms (including stratified intercepts)
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
  ## Weighted intercept (optional)
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
  ## Performance metrics per imputation (optional)
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
  ## Return object
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

  out
}
