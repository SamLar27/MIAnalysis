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
                          # --- NEW: spline options ---
                          spline_terms = NULL,
                          spline_knots_n = NULL,
                          spline_knots_percentile = NULL,
                          # (poly_terms / poly_degree could be added later if needed)
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

  use_random_effects <- !is.null(random_intercept_var)

  ## ------------------------------------------------------------
  ## 0. PREP: Build rcs() spline basis on FULL DATA (all imputations)
  ##          with GLOBAL knots, then reuse in each imputation.
  ## ------------------------------------------------------------
  spline_info <- list()

  if (!is.null(spline_terms)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("Package 'rms' is required for spline_terms (restricted cubic splines).")
    }

    spline_terms <- unique(spline_terms)

    # Default: if only spline_knots_n is given, use equally spaced percentiles between 5 and 95
    if (is.null(spline_knots_percentile)) {
      if (is.null(spline_knots_n)) {
        # default to 5 knots (like rms)
        spline_knots_n <- 5
      }
      # spread between 5 and 95
      spline_knots_percentile <- seq(5, 95, length.out = spline_knots_n)
    } else {
      # if user supplies percentiles, they define the number of knots
      spline_knots_n <- length(spline_knots_percentile)
    }

    # Helper: build rcs basis on full data for one variable
    build_rcs_for_var <- function(var_name) {
      if (!var_name %in% names(data)) {
        warning("spline term '", var_name, "' not found in data, skipping.")
        return(NULL)
      }
      x_all <- data[[var_name]]
      x_non_na <- x_all[!is.na(x_all)]
      if (length(x_non_na) < spline_knots_n) {
        warning("Not enough non-missing values for spline variable '",
                var_name, "'; skipping spline.")
        return(NULL)
      }

      knots <- as.numeric(stats::quantile(
        x_non_na,
        probs = spline_knots_percentile / 100,
        na.rm = TRUE,
        type = 2
      ))

      # Use rms::rcs with explicit knots
      # Result is an n x (K-1) matrix, where K = number of knots
      basis_full <- rms::rcs(x_all, knots = knots)
      basis_full <- as.matrix(basis_full)

      # Name columns: linear + nl1 + nl2 + ...
      n_basis <- ncol(basis_full)
      if (n_basis < 1L) {
        warning("No spline basis columns created for '", var_name, "'. Skipping.")
        return(NULL)
      }

      basis_names <- character(n_basis)
      basis_names[1] <- paste0(var_name, "_rcs_linear")
      if (n_basis > 1) {
        basis_names[2:n_basis] <- paste0(var_name, "_rcs_nl", seq_len(n_basis - 1L))
      }
      colnames(basis_full) <- basis_names

      # Bind to data (keep original column too)
      for (j in seq_len(n_basis)) {
        nm <- basis_names[j]
        if (nm %in% names(data)) {
          warning("Overwriting existing column '", nm, "' when adding spline basis.")
        }
        data[[nm]] <- basis_full[, j]
      }

      # Store info for later (random slopes + diagnostics)
      spline_info[[var_name]] <<- list(
        knots        = knots,
        percentiles  = spline_knots_percentile,
        basis_names  = basis_names
      )

      invisible(NULL)
    }

    # Build splines for each requested term
    for (v in spline_terms) {
      build_rcs_for_var(v)
    }
  }

  ## ------------------------------------------------------------
  ## 1. Imputation list (after spline columns are added)
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
  ## 2. Build RHS terms, incorporating spline basis instead of raw vars
  ## ------------------------------------------------------------
  rhs_main_terms <- c(predictor_vars, covariables)
  rhs_main_terms <- rhs_main_terms[!is.na(rhs_main_terms) & rhs_main_terms != ""]

  # If spline_terms are used, remove the raw variable from RHS and add basis names instead
  if (!is.null(spline_terms) && length(spline_info) > 0L) {
    for (v in names(spline_info)) {
      basis_names <- spline_info[[v]]$basis_names
      # remove original var from RHS if present
      rhs_main_terms <- setdiff(rhs_main_terms, v)
      # add basis names
      rhs_main_terms <- unique(c(rhs_main_terms, basis_names))
    }
  }

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

  ## ------------------------------------------------------------
  ## 3. Expand random slope variables if they refer to spline_terms
  ## ------------------------------------------------------------
  expand_random_vars <- function(vars) {
    if (is.null(vars)) return(character(0))
    vars <- unique(vars)
    out <- unlist(lapply(vars, function(v) {
      if (!is.null(spline_info[[v]])) {
        # This is a spline base var: use its basis columns
        spline_info[[v]]$basis_names
      } else {
        v
      }
    }))
    unique(out)
  }

  rs_predictor_vars  <- expand_random_vars(predictor_vars_random_slope)
  rs_covariables     <- expand_random_vars(covariables_random_slope)

  ## ------------------------------------------------------------
  ## 4. Build random-effect term string
  ## ------------------------------------------------------------
  build_random_term <- function() {
    if (!use_random_effects) return("")

    rs_vars <- unique(c(rs_predictor_vars, rs_covariables))

    same_group_as_strata <- (!is.null(stratified_intercept_var) &&
                               identical(stratified_intercept_var, random_intercept_var))

    if (length(rs_vars) == 0) {
      if (same_group_as_strata) {
        # pure stratification, no random intercept
        return("")
      } else {
        return(paste0("+ (1 | ", random_intercept_var, ")"))
      }
    } else {
      rs_str <- paste(rs_vars, collapse = " + ")
      if (same_group_as_strata) {
        # random slopes only (no random intercept)
        return(paste0("+ (0 + ", rs_str, " | ", random_intercept_var, ")"))
      } else {
        return(paste0("+ (1 + ", rs_str, " | ", random_intercept_var, ")"))
      }
    }
  }

  random_term <- build_random_term()

  ## ------------------------------------------------------------
  ## 5. Build model formula (string → formula)
  ## ------------------------------------------------------------
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
  ## 6. Fit per imputation
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
    stop("All models failed to fit in IPD_one_stage.")
  }

  ## ------------------------------------------------------------
  ## 7. Pool coefficients across imputations (Rubin)
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
  ## 8. Exponentiate for NB / Poisson / Binomial
  ## ------------------------------------------------------------
  Results$exp_estimate   <- exp(Results$estimate)
  Results$exp_CI95_lower <- exp(Results$`2.5 %`)
  Results$exp_CI95_upper <- exp(Results$`97.5 %`)

  ## ------------------------------------------------------------
  ## 9. Flags (interaction / spline / polynomial)
  ## ------------------------------------------------------------
  if (highlight_interactions) {
    Results$is_interaction <- grepl(":", Results$term)
  } else {
    Results$is_interaction <- FALSE
  }

  # mark spline basis terms
  all_spline_basis <- if (length(spline_info) > 0L) {
    unique(unlist(lapply(spline_info, function(z) z$basis_names)))
  } else character(0)

  Results$is_spline <- Results$term %in% all_spline_basis
  Results$is_polynomial <- grepl("poly\\(", Results$term)

  ## ------------------------------------------------------------
  ## 10. Split intercept vs other terms (including stratified intercepts)
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
  ## 11. Weighted intercept across strata (optional)
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
  ## 12. Performance metrics per imputation (optional)
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
  ## 13. Return object
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
  attr(out, "spline_info")         <- spline_info  # 🔹 NEW: for knot inspection

  out
}
