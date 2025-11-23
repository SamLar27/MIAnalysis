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
                          random_intercept_var = NULL,
                          predictor_vars_random_slope = NULL,
                          covariables_random_slope = NULL,
                          stratified_intercept_var = NULL,
                          time_col = NULL,
                          event_col = NULL,
                          formula_string = NULL,
                          highlight_interactions = TRUE,
                          performance_index = NULL) {
  ## ------------------------------ Flags & packages ------------------------------
  use_random_effects <- !is.null(random_intercept_var)
  use_strata         <- !is.null(stratified_intercept_var)

  # we’ll use the attached packages (typical for analysis scripts);
  # for a CRAN package you’d replace with Imports + ::
  require(dplyr)
  require(MASS)
  require(mice)
  require(rlang)
  require(survival)

  if (use_random_effects) {
    require(lme4)

    if (model_type == "cox") {
      if (!requireNamespace("coxme", quietly = TRUE)) {
        stop("Package 'coxme' is required for Cox models with random effects.")
      }
      require(coxme)
    }

    if (model_type == "nb" || model_type == "gamma") {
      if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        warning("Package 'glmmTMB' not available. Using 'lme4' with Poisson/Gamma as an approximation.")
      } else {
        require(glmmTMB)
      }
    }
  }

  ## ------------------------------ Basic checks ------------------------------
  if (!imp_col %in% names(data)) stop("imp_col not found in data.")
  if (!outcome_var %in% names(data) && model_type != "cox") stop("outcome_var not found in data.")
  if (!followup_offset %in% c("Yes", "No")) stop("followup_offset must be either 'Yes' or 'No'.")
  if (followup_offset == "Yes" && is.null(followup_col)) stop("If followup_offset = 'Yes', followup_col must be provided.")
  if (!is.null(followup_col) && any(data[[followup_col]] <= 0, na.rm = TRUE)) {
    stop("followup_col must be strictly positive for offset.")
  }
  if (use_random_effects && !random_intercept_var %in% names(data)) {
    stop("random_intercept_var not found in data.")
  }
  if (use_strata && !stratified_intercept_var %in% names(data)) {
    stop("stratified_intercept_var not found in data.")
  }
  if (!is.null(followup_col) && !followup_col %in% names(data)) {
    stop("followup_col not found in data.")
  }
  if (!is.null(time_col) && !time_col %in% names(data)) stop("time_col not found in data.")
  if (!is.null(event_col) && !event_col %in% names(data)) stop("event_col not found in data.")

  if (!is.null(predictor_vars_random_slope)) {
    if (!all(predictor_vars_random_slope %in% names(data))) {
      missing_rs <- predictor_vars_random_slope[!predictor_vars_random_slope %in% names(data)]
      stop(paste("predictor_vars_random_slope not found in data:", paste(missing_rs, collapse = ", ")))
    }
  }
  if (!is.null(covariables_random_slope)) {
    if (!all(covariables_random_slope %in% names(data))) {
      missing_rs <- covariables_random_slope[!covariables_random_slope %in% names(data)]
      stop(paste("covariables_random_slope not found in data:", paste(missing_rs, collapse = ", ")))
    }
  }
  if (!use_random_effects &&
      (!is.null(predictor_vars_random_slope) || !is.null(covariables_random_slope))) {
    warning("Random slopes specified but random_intercept_var is NULL. Random slopes will be ignored (no mixed model).")
  }

  ## ------------------------------ Imputations ------------------------------
  actual_imps <- sort(unique(data[[imp_col]]))
  if (is.null(imp_n)) {
    imp_n <- length(actual_imps)
  } else if (imp_n != length(actual_imps)) {
    warning("imp_n does not match the number of unique imputations found in the data. Using only present imputations.")
    imp_n <- length(actual_imps)
  }

  implist <- vector("list", length(actual_imps))
  for (i in seq_along(actual_imps)) {
    imp_val <- actual_imps[i]
    data_i <- dplyr::filter(data, !!rlang::sym(imp_col) == imp_val)
    if (nrow(data_i) == 0) stop(paste("No data for imputation", imp_val))
    implist[[i]] <- data_i
  }

  ## ------------------------------ Helpers ------------------------------

  # Optional: extract variables from interaction terms like "A*B"
  extract_vars_from_terms <- function(terms) {
    all_vars <- character(0)
    for (term in terms) {
      if (grepl("\\*", term)) {
        vars_in_interaction <- unlist(strsplit(term, "\\*"))
        vars_in_interaction <- trimws(vars_in_interaction)
        all_vars <- c(all_vars, vars_in_interaction)
      } else {
        all_vars <- c(all_vars, term)
      }
    }
    unique(all_vars)
  }

  # Expand "*" into ":" terms for interactions if needed
  expand_terms <- function(term) {
    if (grepl("\\*", term)) {
      vars_split <- trimws(unlist(strsplit(term, "\\*")))
      all_combinations <- lapply(seq_along(vars_split), function(k) {
        combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
      })
      unlist(all_combinations)
    } else term
  }

  # Check predictor / covariable presence (when no custom formula)
  if (is.null(formula_string)) {
    non_empty_predictors <- predictor_vars[predictor_vars != "" &
                                             !is.na(predictor_vars) &
                                             !is.null(predictor_vars)]
    individual_predictor_vars <- extract_vars_from_terms(non_empty_predictors)
    if (length(individual_predictor_vars) > 0 &&
        !all(individual_predictor_vars %in% names(data))) {
      missing_vars <- individual_predictor_vars[!individual_predictor_vars %in% names(data)]
      stop(paste("Predictor variables not found in data:",
                 paste(missing_vars, collapse = ", ")))
    }

    if (!is.null(covariables)) {
      individual_covariable_vars <- extract_vars_from_terms(covariables)
      if (!all(individual_covariable_vars %in% names(data))) {
        missing_vars <- individual_covariable_vars[!individual_covariable_vars %in% names(data)]
        stop(paste("Covariables not found in data:",
                   paste(missing_vars, collapse = ", ")))
      }
    }
  }

  # random term builder
  build_random_term_core <- function(rs_vars) {
    if (!use_random_effects) return("")
    rs_vars <- unique(rs_vars)
    same_group_as_strata <- (use_strata &&
                               stratified_intercept_var == random_intercept_var)

    if (length(rs_vars) == 0) {
      if (same_group_as_strata) {
        # same grouping as strata, no random intercept (stratified-only)
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

  build_random_term <- function(current_predictor, covariables_in_model) {
    if (!use_random_effects) return("")
    rs_vars <- character(0)

    if (!is.null(predictor_vars_random_slope) &&
        !is.null(current_predictor) &&
        nzchar(current_predictor) &&
        current_predictor %in% predictor_vars_random_slope) {
      rs_vars <- c(rs_vars, current_predictor)
    }

    if (!is.null(covariables_random_slope) && !is.null(covariables_in_model)) {
      rs_vars <- c(rs_vars, intersect(covariables_random_slope, covariables_in_model))
    }

    build_random_term_core(rs_vars)
  }

  # Build formula for each predictor
  build_formula_for_predictor <- function(current_predictor) {
    offset_term <- if (followup_offset == "Yes") paste("+ offset(log(", followup_col, "))", sep = "") else ""
    strat_glm_piece <- ""
    strat_cox_piece <- ""

    if (use_strata) {
      if (model_type == "cox") {
        strat_cox_piece <- paste("+ strata(", stratified_intercept_var, ")", sep = "")
      } else {
        strat_glm_piece <- paste(" + as.factor(", stratified_intercept_var, ")", sep = "")
      }
    }

    if (is.null(formula_string)) {
      # automatic build
      if (!is.null(covariables)) {
        expanded_covariables <- unique(unlist(lapply(covariables, expand_terms)))
      } else {
        expanded_covariables <- character(0)
      }

      covariables_in_model <- covariables

      if (current_predictor == "" || is.null(current_predictor) || is.na(current_predictor)) {
        # null model
        all_terms <- expanded_covariables
      } else {
        expanded_predictor <- expand_terms(current_predictor)
        current_covariables <- if (is.null(covariables)) character(0) else covariables[covariables != current_predictor]
        covariables_in_model <- current_covariables
        expanded_covariables <- unique(unlist(lapply(current_covariables, expand_terms)))
        all_terms <- c(expanded_predictor, expanded_covariables)
      }

      random_term <- build_random_term(current_predictor, covariables_in_model)

      if (model_type == "cox") {
        if (is.null(time_col) || is.null(event_col)) {
          stop("For Cox regression, time_col and event_col must be provided.")
        }
        formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ ",
                             paste(all_terms, collapse = " + "),
                             strat_cox_piece, random_term)
      } else {
        fixed_rhs   <- paste(all_terms, collapse = " + ")
        formula_str <- paste(outcome_var, "~", fixed_rhs,
                             offset_term, strat_glm_piece, random_term)
      }
    } else {
      # user-supplied formula_string; we just wrap and add strata/random bits if needed
      formula_str <- formula_string

      all_rs_vars <- unique(c(
        if (!is.null(predictor_vars_random_slope)) predictor_vars_random_slope else character(0),
        if (!is.null(covariables_random_slope))   covariables_random_slope   else character(0)
      ))

      if (model_type == "cox") {
        if (is.null(time_col) || is.null(event_col)) {
          stop("For Cox regression, time_col and event_col must be provided.")
        }
        if (!grepl("^Surv\\(", formula_str)) {
          formula_str <- paste("Surv(", time_col, ",", event_col, ") ~", formula_str)
        }
        if (use_strata && !grepl(paste0("strata\\(", stratified_intercept_var, "\\)"), formula_str)) {
          formula_str <- paste(formula_str, paste0("+ strata(", stratified_intercept_var, ")"))
        }
        if (use_random_effects && !grepl("\\|", formula_str)) {
          random_term <- build_random_term_core(all_rs_vars)
          if (nzchar(random_term)) {
            formula_str <- paste(formula_str, random_term)
          }
        }
      } else {
        if (!grepl(paste0("^", outcome_var, "\\s*~"), formula_str)) {
          formula_str <- paste(outcome_var, "~", formula_str)
        }
        if (followup_offset == "Yes" && !grepl("offset\\(log\\(.*\\)\\)", formula_str)) {
          formula_str <- paste(formula_str, offset_term)
        }
        if (use_strata && !grepl(paste0("as\\.factor\\(", stratified_intercept_var, "\\)"), formula_str)) {
          formula_str <- paste(formula_str, paste0("+ as.factor(", stratified_intercept_var, ")"))
        }
        if (use_random_effects && !grepl("\\|", formula_str)) {
          random_term <- build_random_term_core(all_rs_vars)
          if (nzchar(random_term)) {
            formula_str <- paste(formula_str, random_term)
          }
        }
      }
    }
    formula_str
  }

  ## -------------------------- Performance helper --------------------------
  # performance_index can be NULL or subset of:
  # c("Log_lik","AIC","AICc","BIC","BICc","RMSE","C_index")
  if (!is.null(performance_index)) {
    allowed_perf <- c("Log_lik","AIC","AICc","BIC","BICc","RMSE","C_index")
    if (!all(performance_index %in% allowed_perf)) {
      stop("performance_index must be subset of: ",
           paste(allowed_perf, collapse = ", "))
    }
  }

  compute_performance_for_model <- function(m, perf_names) {
    if (is.null(perf_names) || length(perf_names) == 0) return(NULL)
    res <- setNames(as.list(rep(NA_real_, length(perf_names))), perf_names)
    if (is.null(m)) return(res)

    # sample size & number of parameters for AICc/BICc
    n <- tryCatch(stats::nobs(m), error = function(e) NA_integer_)
    k <- tryCatch(length(stats::coef(m)), error = function(e) NA_integer_)

    if ("Log_lik" %in% perf_names) {
      res[["Log_lik"]] <- tryCatch(as.numeric(stats::logLik(m)),
                                   error = function(e) NA_real_)
    }
    if ("AIC" %in% perf_names) {
      res[["AIC"]] <- tryCatch(stats::AIC(m), error = function(e) NA_real_)
    }
    if ("BIC" %in% perf_names) {
      res[["BIC"]] <- tryCatch(stats::BIC(m), error = function(e) NA_real_)
    }
    if ("AICc" %in% perf_names && !is.na(n) && !is.na(k) && n > (k + 1)) {
      aic_val <- if (!is.null(res[["AIC"]]) && !is.na(res[["AIC"]])) {
        res[["AIC"]]
      } else {
        tryCatch(stats::AIC(m), error = function(e) NA_real_)
      }
      if (!is.na(aic_val)) {
        res[["AICc"]] <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
      }
    }
    if ("BICc" %in% perf_names && !is.na(n) && !is.na(k) && n > (k + 1)) {
      bic_val <- if (!is.null(res[["BIC"]]) && !is.na(res[["BIC"]])) {
        res[["BIC"]]
      } else {
        tryCatch(stats::BIC(m), error = function(e) NA_real_)
      }
      if (!is.na(bic_val)) {
        # simple corrected BIC; not universally standard but reasonable
        res[["BICc"]] <- bic_val + (k * (k + 1)) / (n - k - 1)
      }
    }
    if ("RMSE" %in% perf_names) {
      if (!inherits(m, "coxph") && !inherits(m, "coxme")) {
        y <- tryCatch({
          mf <- stats::model.frame(m)
          as.numeric(stats::model.response(mf))
        }, error = function(e) NA_real_)
        mu <- tryCatch(as.numeric(stats::predict(m, type = "response")),
                       error = function(e) NA_real_)
        if (!any(is.na(y)) && length(y) == length(mu) && length(y) > 0) {
          res[["RMSE"]] <- sqrt(mean((y - mu)^2, na.rm = TRUE))
        } else {
          res[["RMSE"]] <- NA_real_
        }
      } else {
        res[["RMSE"]] <- NA_real_
      }
    }
    if ("C_index" %in% perf_names) {
      if (inherits(m, "coxph")) {
        res[["C_index"]] <- tryCatch({
          survival::concordance(m)$concordance
        }, error = function(e) NA_real_)
      } else {
        res[["C_index"]] <- NA_real_
      }
    }

    res
  }

  ## ------------------------ Fit one predictor ------------------------
  fit_model_for_predictor <- function(current_predictor) {

    formula_string_current <- build_formula_for_predictor(current_predictor)
    model_formula <- as.formula(formula_string_current)

    # pick up interaction terms just for flags
    interaction_terms <- character(0)
    if (grepl(":", formula_string_current)) {
      terms_part <- strsplit(formula_string_current, "~")[[1]][2]
      terms <- trimws(strsplit(terms_part, "\\+")[[1]])
      interaction_terms <- terms[grepl(":", terms)]
    }

    is_null_model <- (current_predictor == "" ||
                        is.null(current_predictor) ||
                        is.na(current_predictor))

    current_models_list <- vector("list", length(actual_imps))
    fit_log <- data.frame(
      imp   = actual_imps,
      ok    = NA,
      error = NA_character_,
      stringsAsFactors = FALSE
    )

    ## ------------------------ Random-effects block ------------------------
    if (use_random_effects) {

      for (i in seq_along(actual_imps)) {
        data_i <- implist[[i]]
        imp_val <- actual_imps[i]

        current_models_list[[i]] <- tryCatch({
          if (model_type == "lm") {
            lme4::lmer(model_formula, data = data_i)
          } else if (model_type %in% c("bin")) {
            lme4::glmer(model_formula, family = binomial(), data = data_i)
          } else if (model_type %in% c("poisson", "quasipoisson")) {
            lme4::glmer(model_formula, family = poisson(), data = data_i)
          } else if (model_type == "nb") {
            if (requireNamespace("glmmTMB", quietly = TRUE)) {
              glmmTMB::glmmTMB(model_formula, family = glmmTMB::nbinom2, data = data_i)
            } else {
              warning("Using Poisson instead of Negative Binomial. Install glmmTMB for proper NB mixed models.")
              lme4::glmer(model_formula, family = poisson(), data = data_i)
            }
          } else if (model_type == "cox") {
            coxme::coxme(model_formula, data = data_i)
          } else if (model_type == "gamma") {
            if (requireNamespace("glmmTMB", quietly = TRUE)) {
              glmmTMB::glmmTMB(model_formula, family = glmmTMB::Gamma, data = data_i)
            } else {
              warning("Gamma family with random effects may not be stable in lme4.")
              lme4::glmer(model_formula, family = Gamma, data = data_i)
            }
          } else stop("Unsupported model type for random effects.")
        }, error = function(e) {
          fit_log$ok[fit_log$imp == imp_val]    <<- FALSE
          fit_log$error[fit_log$imp == imp_val] <<- e$message
          warning(sprintf("Model failed to fit for imputation %s: %s", imp_val, e$message))
          NULL
        })

        if (!is.null(current_models_list[[i]])) {
          fit_log$ok[fit_log$imp == imp_val]    <- TRUE
          fit_log$error[fit_log$imp == imp_val] <- NA_character_
        }
      }

      if (all(vapply(current_models_list, is.null, logical(1)))) {
        stop("All models failed to fit in IPD_one_stage (random-effects block).")
      }

      # --- relaxed coefficient extraction ---
      coefs <- lapply(current_models_list, function(m) {
        if (is.null(m)) return(NULL)

        # initialize
        c_est <- NULL
        c_se  <- NULL

        if (model_type == "cox") {
          # coxme: fixed effects via fixef; vcov may fail
          c_est <- tryCatch(coef(m), error = function(e) NULL)
          if (is.null(c_est)) return(NULL)
          c_se <- tryCatch({
            v <- as.matrix(vcov(m))
            sqrt(diag(v))
          }, error = function(e) {
            # last-resort: NA SE
            rep(NA_real_, length(c_est))
          })

        } else if (inherits(m, "glmmTMB")) {
          c_est <- tryCatch(glmmTMB::fixef(m)$cond,
                            error = function(e) NULL)
          if (is.null(c_est)) return(NULL)
          c_se <- tryCatch({
            v <- glmmTMB::vcov(m)$cond
            sqrt(diag(v))
          }, error = function(e) {
            # fallback to summary
            s <- tryCatch(summary(m), error = function(e2) NULL)
            if (!is.null(s) && !is.null(s$coefficients$cond)) {
              se_vec <- s$coefficients$cond[, "Std. Error"]
              se_vec[names(c_est)]
            } else {
              rep(NA_real_, length(c_est))
            }
          })

        } else if (inherits(m, "merMod")) {
          c_est <- lme4::fixef(m)
          c_se  <- tryCatch({
            s <- summary(m)
            s$coefficients[, "Std. Error"]
          }, error = function(e) {
            # fallback to vcov
            tryCatch({
              v <- as.matrix(vcov(m))
              sqrt(diag(v))
            }, error = function(e2) {
              rep(NA_real_, length(c_est))
            })
          })

        } else {
          # just in case we get a plain glm here
          c_est <- tryCatch(stats::coef(m), error = function(e) NULL)
          if (is.null(c_est)) return(NULL)
          c_se  <- tryCatch({
            v <- as.matrix(stats::vcov(m))
            sqrt(diag(v))
          }, error = function(e) {
            rep(NA_real_, length(c_est))
          })
        }

        data.frame(
          term      = names(c_est),
          estimate  = as.numeric(c_est),
          std.error = as.numeric(c_se),
          stringsAsFactors = FALSE
        )
      })

      coefs <- Filter(function(x) !is.null(x), coefs)
      if (length(coefs) == 0) {
        stop("No usable coefficient table could be extracted for any imputation (random-effects block).")
      }

      # pooled results using relaxed Rubin (allows NA SE, uses subset)
      all_terms <- unique(unlist(lapply(coefs, function(df) df$term)))
      pooled_results <- data.frame(
        term      = all_terms,
        estimate  = NA_real_,
        std.error = NA_real_,
        stringsAsFactors = FALSE
      )

      for (term in all_terms) {
        # gather estimates and SEs (some SE may be NA)
        term_ests <- sapply(coefs, function(df) {
          if (term %in% df$term) df$estimate[df$term == term] else NA_real_
        })
        term_ses <- sapply(coefs, function(df) {
          if (term %in% df$term) df$std.error[df$term == term] else NA_real_
        })

        term_ests <- term_ests[!is.na(term_ests)]
        valid_se  <- !is.na(term_ses)
        term_ses_valid <- term_ses[valid_se]

        if (length(term_ests) == 0) next

        Q_bar <- mean(term_ests)

        if (length(term_ses_valid) >= 2) {
          U_bar <- mean(term_ses_valid^2)
          B     <- stats::var(term_ests)
          m     <- length(term_ests)
          T_var <- U_bar + (1 + 1/m) * B
          pooled_se <- sqrt(T_var)
        } else {
          # not enough SE information; keep estimate, SE unknown
          pooled_se <- NA_real_
        }

        pooled_results$estimate[pooled_results$term == term]  <- Q_bar
        pooled_results$std.error[pooled_results$term == term] <- pooled_se
      }

      Results_multivariate_analysis <- pooled_results %>%
        mutate(
          `2.5 %` = ifelse(!is.na(std.error),
                           estimate - 1.96 * std.error,
                           NA_real_),
          `97.5 %` = ifelse(!is.na(std.error),
                            estimate + 1.96 * std.error,
                            NA_real_),
          p.value = ifelse(!is.na(std.error),
                           2 * pnorm(-abs(estimate / std.error)),
                           NA_real_)
        )

    } else {
      ## ------------------------ Non-random effects (GLM / Cox) ------------------------
      res_comb <- vector("list", length(actual_imps))
      for (i in seq_along(actual_imps)) {
        data_subset <- implist[[i]]
        imp_val     <- actual_imps[i]

        res_comb[[i]] <- tryCatch({
          switch(model_type,
                 "nb"            = MASS::glm.nb(model_formula, data = data_subset),
                 "lm"            = glm(model_formula, family = gaussian(), data = data_subset),
                 "bin"           = glm(model_formula, family = binomial(), data = data_subset),
                 "poisson"       = glm(model_formula, family = poisson(), data = data_subset),
                 "gamma"         = glm(model_formula, family = Gamma(), data = data_subset),
                 "quasipoisson"  = glm(model_formula, family = quasipoisson(), data = data_subset),
                 "quasibinomial" = glm(model_formula, family = quasibinomial(), data = data_subset),
                 "cox"           = coxph(model_formula, data = data_subset),
                 stop("Unsupported model type.")
          )
        }, error = function(e) {
          fit_log$ok[fit_log$imp == imp_val]    <<- FALSE
          fit_log$error[fit_log$imp == imp_val] <<- e$message
          warning(sprintf("Model failed to fit for imputation %s: %s", imp_val, e$message))
          NULL
        })

        if (!is.null(res_comb[[i]])) {
          fit_log$ok[fit_log$imp == imp_val]    <- TRUE
          fit_log$error[fit_log$imp == imp_val] <- NA_character_
        }
      }

      ok_models <- Filter(function(x) !is.null(x), res_comb)
      if (length(ok_models) == 0) {
        stop("All models failed to fit in IPD_one_stage.")
      }

      pooled <- mice::pool(ok_models)
      Results_multivariate_analysis <- summary(pooled, conf.int = TRUE, exp = FALSE)
      current_models_list <- res_comb
    }

    ## ------------------------ Performance indices (relaxed) ------------------------
    performance_per_imp <- NULL
    if (!is.null(performance_index)) {
      performance_per_imp <- lapply(seq_along(current_models_list), function(i) {
        m <- current_models_list[[i]]
        out <- compute_performance_for_model(m, performance_index)
        out
      })
      names(performance_per_imp) <- paste0("imp_", actual_imps)
    }

    ## ------------------------ Post-processing results table ------------------------
    Results_multivariate_analysis <- Results_multivariate_analysis %>%
      mutate(
        exp_estimate   = exp(estimate),
        exp_CI95_lower = exp(`2.5 %`),
        exp_CI95_upper = exp(`97.5 %`)
      ) %>%
      dplyr::select(term, estimate, std.error, `2.5 %`, `97.5 %`,
                    exp_estimate, exp_CI95_lower, exp_CI95_upper, p.value)

    # Null model row
    if (is_null_model) {
      null_row <- data.frame(
        term = "No_predictor",
        estimate = NA, std.error = NA, `2.5 %` = NA, `97.5 %` = NA,
        exp_estimate = NA, exp_CI95_lower = NA, exp_CI95_upper = NA, p.value = NA,
        stringsAsFactors = FALSE, check.names = FALSE
      )
      Results_multivariate_analysis <- rbind(null_row, Results_multivariate_analysis)
    }

    # Interaction / spline / polynomial flags (minimal)
    if (highlight_interactions && length(interaction_terms) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_interaction = grepl(":", term))
    } else {
      Results_multivariate_analysis$is_interaction <- FALSE
    }
    Results_multivariate_analysis$is_spline      <- FALSE
    Results_multivariate_analysis$is_polynomial  <- grepl("poly\\(", Results_multivariate_analysis$term)

    ## ------------------------ Attributes ------------------------
    attr(Results_multivariate_analysis, "has_random_effects")    <- use_random_effects
    attr(Results_multivariate_analysis, "random_intercept_var")  <- random_intercept_var
    attr(Results_multivariate_analysis, "predictor_tested")      <- current_predictor
    attr(Results_multivariate_analysis, "model_type_used")       <- model_type
    attr(Results_multivariate_analysis, "formula")               <- formula_string_current
    attr(Results_multivariate_analysis, "has_interactions")      <- length(interaction_terms) > 0
    attr(Results_multivariate_analysis, "interaction_terms")     <- if (length(interaction_terms) > 0) interaction_terms else NULL
    attr(Results_multivariate_analysis, "imputations")           <- actual_imps
    attr(Results_multivariate_analysis, "n_imp")                 <- length(actual_imps)
    attr(Results_multivariate_analysis, "stratified_intercept_var") <- stratified_intercept_var

    # Store models, performance, and fit_log
    attr(Results_multivariate_analysis, "models")              <- current_models_list
    attr(Results_multivariate_analysis, "performance_per_imp") <- performance_per_imp
    attr(Results_multivariate_analysis, "fit_log")             <- fit_log

    Results_multivariate_analysis
  }

  ## ------------------------ Main loop over predictors ------------------------
  results_list <- vector("list", length(predictor_vars))
  names(results_list) <- predictor_vars

  for (i in seq_along(predictor_vars)) {
    current_predictor <- predictor_vars[i]
    if (current_predictor == "" || is.null(current_predictor) || is.na(current_predictor)) {
      message("Fitting null model (no predictor)")
      results_list[["No_predictor"]] <- fit_model_for_predictor(current_predictor)
    } else {
      message(paste("Fitting model for predictor:", current_predictor))
      results_list[[i]] <- fit_model_for_predictor(current_predictor)
    }
  }

  if (length(predictor_vars) == 1) {
    return(results_list[[1]])
  } else {
    # combined table with main effect rows for each predictor
    combined_results <- data.frame()
    for (pred_name in names(results_list)) {
      result <- results_list[[pred_name]]
      if (is.null(result) || nrow(result) == 0) next
      if (pred_name == "No_predictor") {
        predictor_rows <- result[result$term == "No_predictor", ]
      } else {
        predictor_rows <- result[grepl(paste0("^", pred_name), result$term) &
                                   !grepl("Intercept", result$term), ]
      }
      if (nrow(predictor_rows) > 0) {
        predictor_rows <- predictor_rows[, c("term","estimate","std.error","2.5 %","97.5 %",
                                             "exp_estimate","exp_CI95_lower","exp_CI95_upper","p.value")]
        combined_results <- rbind(combined_results, predictor_rows)
      }
    }
    rownames(combined_results) <- NULL

    class(results_list) <- c("IPD_one_stage_multi", "list")
    attr(results_list, "predictors")       <- predictor_vars
    attr(results_list, "model_type")       <- model_type
    attr(results_list, "n_predictors")     <- length(predictor_vars)
    attr(results_list, "combined_results") <- combined_results
    return(results_list)
  }
}
