#' Compute Estimates for Multiple Imputed Data
#'
#' Fits models to multiply imputed datasets and pools results using Rubin's Rules.
#' Supports GLM, mixed models (via lme4 / glmmTMB / coxme), interactions,
#' restricted cubic splines (rcs) and polynomial terms.
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
#' @param data Data frame with all imputations stacked.
#' @param outcome_var Dependent variable (for non-Cox models).
#' @param predictor_vars Character vector of predictors to be tested one by one.
#' @param covariables Character vector of covariables (always included as fixed effects).
#' @param imp_col Column name indicating imputation index (default ".imp").
#' @param imp_n Number of imputations (if NULL, detected automatically).
#' @param model_type "nb","lm","bin","poisson","gamma","quasipoisson","quasibinomial","cox".
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
#'
#' @return
#'  If one predictor: a data.frame of pooled coefficients with attributes:
#'    - $term, estimate, std.error, 2.5%, 97.5%, exp_estimate, exp_CI95_lower, exp_CI95_upper, p.value
#'    - attributes: "fit_stats_per_imp", "fit_stats_pooled", "D1_components", etc.
#'  If multiple predictors: list of such data.frames with class "MI_estimates_multi"
#'    and attribute "combined_results".
#' @export
MI_estimates <- function(data,
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
                         poly_terms = NULL,
                         poly_degree = NULL,
                         include_poly_terms = TRUE,
                         random_intercept_var = NULL,
                         predictor_vars_random_slope = NULL,
                         covariables_random_slope = NULL,
                         stratified_intercept_var = NULL) {

  ## ---- Flags inferred from arguments ----
  use_random_effects <- !is.null(random_intercept_var)
  use_strata         <- !is.null(stratified_intercept_var)

  ## ---- Packages ----
  require(dplyr)
  require(MASS)
  require(mice)
  require(rlang)
  require(survival)
  require(splines)

  if (use_random_effects) {
    require(lme4)
    require(mitools)

    if (model_type == "cox") {
      if (!requireNamespace("coxme", quietly = TRUE)) {
        stop("Package 'coxme' is required for Cox models with random effects.")
      }
      require(coxme)
    }

    if (model_type == "nb") {
      if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        warning("Package 'glmmTMB' not available. Using 'lme4' with Poisson family as approximation.")
      } else {
        require(glmmTMB)
      }
    }
  }

  # Spline support (rcs via rms if available)
  use_rms <- FALSE
  if (!is.null(spline_terms)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      message("Package 'rms' is not available. Using bs() from splines instead of rcs().")
      use_rms <- FALSE
    } else {
      require(rms)
      use_rms <- TRUE
    }
  }

  ## ---- Core checks ----
  if (!imp_col %in% names(data)) stop("imp_col not found in data.")
  if (!outcome_var %in% names(data) && model_type != "cox") stop("outcome_var not found in data.")
  if (!followup_offset %in% c("Yes", "No")) stop("followup_offset must be either 'Yes' or 'No'.")
  if (followup_offset == "Yes" && is.null(followup_col)) stop("If followup_offset = 'Yes', followup_col must be provided.")
  if (!is.null(followup_col) && any(data[[followup_col]] <= 0, na.rm = TRUE)) stop("followup_col must be strictly positive for offset.")
  if (!trial_factor %in% c("Yes", "No")) stop("trial_factor must be either 'Yes' or 'No'.")
  if (trial_factor == "Yes" && is.null(trial_col)) stop("If trial_factor = 'Yes', trial_col must be provided.")
  if (use_random_effects && !random_intercept_var %in% names(data)) stop("random_intercept_var not found in data.")
  if (use_strata && !stratified_intercept_var %in% names(data)) stop("stratified_intercept_var not found in data.")
  if (!is.null(followup_col) && !followup_col %in% names(data)) stop("followup_col not found in data.")
  if (!is.null(trial_col) && !trial_col %in% names(data)) stop("trial_col not found in data.")
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

  ## ---- Imputations ----
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

  ## ---- Spline validation and knot computation ----
  spline_map <- list()  # var -> formula piece (rcs(...) or bs(...))
  if (!is.null(spline_terms)) {
    if (!is.character(spline_terms)) stop("spline_terms must be a character vector of variable names.")
    for (v in spline_terms) {
      if (!v %in% names(data)) stop(paste("Spline variable", v, "not found in data"))
    }
    if (is.null(spline_knots_n)) stop("spline_knots_n must be provided when using spline_terms.")
    if (is.null(spline_knots_percentile)) stop("spline_knots_percentile must be provided when using spline_terms.")
    if (length(spline_knots_percentile) != spline_knots_n)
      stop("Length of spline_knots_percentile must match spline_knots_n.")

    for (v in spline_terms) {
      knots <- as.numeric(stats::quantile(data[[v]],
                                          probs = spline_knots_percentile / 100,
                                          na.rm = TRUE))
      if (use_rms) {
        spline_map[[v]] <- paste0("rcs(", v, ", c(", paste(knots, collapse = ", "), "))")
      } else {
        df <- spline_knots_n + 1  # approx df
        spline_map[[v]] <- paste0("bs(", v, ", df = ", df, ", degree = 3)")
      }
    }
  }

  ## ---- Polynomial validation and degree handling ----
  poly_map <- list()      # var -> formula piece poly(var, degree, raw=TRUE)
  poly_degree_vec <- NULL # vector of degrees aligned with poly_terms

  if (!is.null(poly_terms)) {
    if (!is.character(poly_terms)) stop("poly_terms must be a character vector of variable names.")
    for (v in poly_terms) {
      if (!v %in% names(data)) stop(paste("Polynomial variable", v, "not found in data"))
    }

    if (is.null(poly_degree)) {
      poly_degree_vec <- rep(2L, length(poly_terms))
    } else {
      if (!is.numeric(poly_degree)) stop("poly_degree must be numeric (2 or 3).")
      if (length(poly_degree) == 1) {
        poly_degree_vec <- rep(as.integer(poly_degree), length(poly_terms))
      } else if (length(poly_degree) == length(poly_terms)) {
        poly_degree_vec <- as.integer(poly_degree)
      } else {
        stop("poly_degree must be length 1 or the same length as poly_terms.")
      }
      if (!all(poly_degree_vec %in% c(2L, 3L))) stop("poly_degree must be 2 or 3.")
    }

    for (i in seq_along(poly_terms)) {
      v   <- poly_terms[i]
      deg <- poly_degree_vec[i]
      poly_map[[v]] <- paste0("poly(", v, ", degree = ", deg, ", raw = TRUE)")
    }
  }

  ## ---- Variable presence checks (if no custom formula) ----
  if (is.null(formula_string)) {
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

    non_empty_predictors <- predictor_vars[predictor_vars != "" & !is.na(predictor_vars) & !is.null(predictor_vars)]
    individual_predictor_vars <- extract_vars_from_terms(non_empty_predictors)
    if (length(individual_predictor_vars) > 0 && !all(individual_predictor_vars %in% names(data))) {
      missing_vars <- individual_predictor_vars[!individual_predictor_vars %in% names(data)]
      stop(paste("Predictor variables not found in data:", paste(missing_vars, collapse = ", ")))
    }

    if (!is.null(covariables)) {
      individual_covariable_vars <- extract_vars_from_terms(covariables)
      if (!all(individual_covariable_vars %in% names(data))) {
        missing_vars <- individual_covariable_vars[!individual_covariable_vars %in% names(data)]
        stop(paste("Covariables not found in data:", paste(missing_vars, collapse = ", ")))
      }
    }
  }

  ## ---- Helpers ----

  expand_terms <- function(term) {
    if (grepl("^rcs\\(|^bs\\(|^poly\\(", term)) return(term)
    if (grepl("\\*", term)) {
      vars_split <- trimws(unlist(strsplit(term, "\\*")))
      all_combinations <- lapply(seq_along(vars_split), function(k) {
        combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
      })
      unlist(all_combinations)
    } else term
  }

  ## Build vectors of special terms (strings used in formulas)
  process_special_terms <- function() {
    spline_formula_parts <- character(0)
    if (length(spline_map) > 0) {
      spline_formula_parts <- unlist(unname(spline_map))
    }

    poly_formula_parts <- character(0)
    if (length(poly_map) > 0) {
      poly_formula_parts <- unlist(unname(poly_map))
    }

    list(spline_parts = spline_formula_parts, poly_parts = poly_formula_parts)
  }

  special_terms <- process_special_terms()
  spline_formula_parts <- special_terms$spline_parts
  poly_formula_parts   <- special_terms$poly_parts

  ## ---- Random-effects core builder ----
  build_random_term_core <- function(rs_vars) {
    if (!use_random_effects) return("")
    rs_vars <- unique(rs_vars)
    same_group_as_strata <- (use_strata &&
                               stratified_intercept_var == random_intercept_var)

    if (length(rs_vars) == 0) {
      if (same_group_as_strata) {
        # pure stratified intercept, no random effects
        return("")
      } else {
        return(paste0("+ (1 | ", random_intercept_var, ")"))
      }
    } else {
      rs_str <- paste(rs_vars, collapse = " + ")
      if (same_group_as_strata) {
        # stratified intercept + random slopes only
        return(paste0("+ (0 + ", rs_str, " | ", random_intercept_var, ")"))
      } else {
        # random intercept + random slopes
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

  ## ---- Fit stats helper (for each imputed model) ----
  compute_fit_stats <- function(model, outcome_var, model_type) {

    # If model is NULL or not a fitted object, return a 1-row NA table
    if (is.null(model)) {
      return(data.frame(
        logLik   = NA_real_,
        df       = NA_real_,
        AIC      = NA_real_,
        AICc     = NA_real_,
        BIC      = NA_real_,
        BICc     = NA_real_,
        RMSE     = NA_real_,
        MAE      = NA_real_,
        pseudoR2 = NA_real_
      ))
    }

    out <- list(
      logLik   = NA_real_,
      df       = NA_real_,
      AIC      = NA_real_,
      AICc     = NA_real_,
      BIC      = NA_real_,
      BICc     = NA_real_,
      RMSE     = NA_real_,
      MAE      = NA_real_,
      pseudoR2 = NA_real_
    )

    # logLik and df
    ll <- tryCatch(stats::logLik(model), error = function(e) NULL)
    if (!is.null(ll)) {
      out$logLik <- as.numeric(ll)
      out$df     <- attr(ll, "df")
    }

    # AIC / BIC
    out$AIC <- tryCatch(stats::AIC(model), error = function(e) NA_real_)
    out$BIC <- tryCatch(stats::BIC(model), error = function(e) NA_real_)

    # AICc / BICc
    n <- tryCatch(stats::nobs(model), error = function(e) NA_integer_)
    k <- out$df
    if (!is.na(n) && !is.na(k) && n > (k + 1) && !is.na(out$AIC) && !is.na(out$BIC)) {
      correction <- 2 * k * (k + 1) / (n - k - 1)
      out$AICc <- out$AIC + correction
      out$BICc <- out$BIC + correction
    }

    # RMSE / MAE / pseudoR2 – not meaningful for Cox, so skip
    if (model_type != "cox") {
      y_hat <- tryCatch(stats::predict(model, type = "response"), error = function(e) NULL)
      y_obs <- tryCatch({
        mf <- stats::model.frame(model)
        if (outcome_var %in% names(mf)) {
          mf[[outcome_var]]
        } else {
          stats::model.response(mf)
        }
      }, error = function(e) NULL)

      if (!is.null(y_hat) && !is.null(y_obs)) {
        if (is.factor(y_obs)) {
          if (nlevels(y_obs) == 2) {
            y_obs <- as.numeric(y_obs) - 1
          } else {
            y_obs <- as.numeric(y_obs)
          }
        }
        if (length(y_hat) == length(y_obs)) {
          resid <- y_obs - y_hat
          out$RMSE <- sqrt(mean(resid^2, na.rm = TRUE))
          out$MAE  <- mean(abs(resid), na.rm = TRUE)
        }
      }

      # pseudo-R²
      if (inherits(model, "glm") || inherits(model, "negbin")) {
        dev      <- tryCatch(model$deviance,      error = function(e) NA_real_)
        null_dev <- tryCatch(model$null.deviance, error = function(e) NA_real_)
        if (!is.na(dev) && !is.na(null_dev) && null_dev != 0) {
          out$pseudoR2 <- 1 - dev / null_dev
        }
      } else if (inherits(model, "lm")) {
        out$pseudoR2 <- tryCatch(summary(model)$r.squared, error = function(e) NA_real_)
      } else if (inherits(model, "glmmTMB") || inherits(model, "merMod")) {
        if (requireNamespace("performance", quietly = TRUE)) {
          out$pseudoR2 <- tryCatch({
            r2 <- performance::r2_nakagawa(model)
            as.numeric(r2$R2m)
          }, error = function(e) NA_real_)
        }
      }
    }

    as.data.frame(out)
  }

  ## ---- Formula builder ----
  build_formula_for_predictor <- function(current_predictor) {
    trial_term  <- if (trial_factor == "Yes") paste("+ as.factor(", trial_col, ")", sep = "") else ""
    offset_term <- if (followup_offset == "Yes") paste("+ offset(log(", followup_col, "))", sep = "") else ""

    strat_glm_piece <- ""
    strat_cox_piece <- ""

    if (use_strata) {
      if (model_type == "cox") {
        strat_cox_piece <- paste("+ strata(", stratified_intercept_var, ")", sep = "")
      } else {
        # GLM: keep intercept + add as.factor(strata)
        strat_glm_piece <- paste(" + as.factor(", stratified_intercept_var, ")", sep = "")
      }
    }

    if (is.null(formula_string)) {
      # Automatic formula construction
      if (!is.null(covariables)) {
        expanded_covariables <- unique(unlist(lapply(covariables, expand_terms)))
      } else {
        expanded_covariables <- character(0)
      }

      covariables_in_model <- covariables

      # Special handling for current predictor if spline/poly is requested
      if (current_predictor == "" || is.null(current_predictor) || is.na(current_predictor)) {
        # Null model: only covariables and special terms that are NOT tied to a specific predictor
        all_terms <- c(expanded_covariables, spline_formula_parts, poly_formula_parts)
      } else {
        current_spline_parts <- spline_formula_parts
        current_poly_parts   <- poly_formula_parts

        if (!is.null(spline_terms) && current_predictor %in% spline_terms) {
          main_pred_term <- spline_map[[current_predictor]]
          current_spline_parts <- setdiff(current_spline_parts, main_pred_term)
          expanded_predictor <- main_pred_term
        } else if (!is.null(poly_terms) && current_predictor %in% poly_terms) {
          main_pred_term <- poly_map[[current_predictor]]
          current_poly_parts <- setdiff(current_poly_parts, main_pred_term)
          expanded_predictor <- main_pred_term
        } else {
          expanded_predictor <- expand_terms(current_predictor)
        }

        current_covariables <- if (is.null(covariables)) character(0) else covariables[covariables != current_predictor]
        covariables_in_model <- current_covariables

        expanded_covariables <- unique(unlist(lapply(current_covariables, expand_terms)))
        all_terms <- c(expanded_predictor, expanded_covariables, current_spline_parts, current_poly_parts)
      }

      random_term <- build_random_term(current_predictor, covariables_in_model)

      if (model_type == "cox") {
        if (is.null(time_col) || is.null(event_col)) stop("For Cox regression, time_col and event_col must be provided.")
        formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ ",
                             paste(all_terms, collapse = " + "), trial_term, strat_cox_piece, random_term, sep = "")
      } else {
        fixed_rhs <- paste(all_terms, collapse = " + ")
        formula_str <- paste(outcome_var, "~", fixed_rhs, offset_term, trial_term, strat_glm_piece, random_term)
      }

    } else {
      # Custom formula: add trial/offset/strata/random if missing
      formula_str <- formula_string

      all_rs_vars <- unique(c(
        if (!is.null(predictor_vars_random_slope)) predictor_vars_random_slope else character(0),
        if (!is.null(covariables_random_slope))   covariables_random_slope   else character(0)
      ))

      if (model_type == "cox") {
        if (is.null(time_col) || is.null(event_col)) stop("For Cox regression, time_col and event_col must be provided.")
        if (!grepl("^Surv\\(", formula_str)) {
          formula_str <- paste("Surv(", time_col, ",", event_col, ") ~", formula_str)
        }
        if (trial_factor == "Yes" && !grepl(paste0("as\\.factor\\(", trial_col, "\\)"), formula_str)) {
          formula_str <- paste(formula_str, trial_term)
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
        if (trial_factor == "Yes" && !grepl(paste0("as\\.factor\\(", trial_col, "\\)"), formula_str)) {
          formula_str <- paste(formula_str, trial_term)
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

  ## ---- Fit one predictor ----
  fit_model_for_predictor <- function(current_predictor) {

    formula_string_current <- build_formula_for_predictor(current_predictor)
    model_formula <- as.formula(formula_string_current)

    interaction_terms <- character(0)
    if (grepl(":", formula_string_current)) {
      terms_part <- strsplit(formula_string_current, "~")[[1]][2]
      terms <- trimws(strsplit(terms_part, "\\+")[[1]])
      interaction_terms <- terms[grepl(":", terms)]
    }

    spline_terms_detected <- character(0)
    if (grepl("rcs\\(|bs\\(", formula_string_current)) {
      terms_part <- strsplit(formula_string_current, "~")[[1]][2]
      terms <- trimws(strsplit(terms_part, "\\+")[[1]])
      spline_terms_detected <- terms[grepl("rcs\\(|bs\\(", terms)]
    }

    poly_terms_detected <- character(0)
    if (grepl("poly\\(", formula_string_current)) {
      terms_part <- strsplit(formula_string_current, "~")[[1]][2]
      terms <- trimws(strsplit(terms_part, "\\+")[[1]])
      poly_terms_detected <- terms[grepl("poly\\(", terms)]
    }

    is_null_model <- (current_predictor == "" || is.null(current_predictor) || is.na(current_predictor))

    current_models_list <- NULL
    fit_stats_per_imp   <- NULL
    fit_stats_pooled    <- NULL
    D1_components       <- NULL

    if (use_random_effects) {

      models_list <- vector("list", imp_n)
      for (i in seq_along(actual_imps)) {
        data_i <- implist[[i]]
        models_list[[i]] <- tryCatch({
          if (model_type == "lm") {
            lme4::lmer(model_formula, data = data_i)
          } else if (model_type %in% c("bin")) {
            lme4::glmer(model_formula, family = binomial(), data = data_i)
          } else if (model_type %in% c("poisson", "quasipoisson")) {
            lme4::glmer(model_formula, family = poisson(), data = data_i)
          } else if (model_type == "nb") {
            if (requireNamespace("glmmTMB", quietly = TRUE)) {
              glmmTMB::glmmTMB(model_formula, family = nbinom2, data = data_i)
            } else {
              warning("Using Poisson instead of Negative Binomial. Install glmmTMB for proper NB mixed models.")
              lme4::glmer(model_formula, family = poisson(), data = data_i)
            }
          } else if (model_type == "cox") {
            coxme::coxme(model_formula, data = data_i)
          } else if (model_type == "gamma") {
            if (requireNamespace("glmmTMB", quietly = TRUE)) {
              glmmTMB::glmmTMB(model_formula, family = Gamma, data = data_i)
            } else {
              warning("Gamma family with random effects may not be stable in lme4.")
              lme4::glmer(model_formula, family = Gamma, data = data_i)
            }
          } else stop("Unsupported model type for random effects.")
        }, error = function(e) {
          warning(sprintf("Model failed to fit for imputation %s: %s", actual_imps[i], e$message))
          NULL
        })
      }

      current_models_list <- models_list

      # --- Pooled coefficients via Rubin ---
      coefs <- lapply(models_list, function(m) {
        if (is.null(m)) return(NULL)
        if (model_type == "cox") {
          c_est <- fixef(m)
          c_se  <- sqrt(diag(vcov(m)))
        } else if (inherits(m, "glmmTMB")) {
          c_est <- fixef(m)$cond
          v     <- vcov(m)$cond
          c_se  <- sqrt(diag(v))
        } else {
          c_est <- fixef(m)
          if (inherits(m, "lmerMod")) {
            model_summary <- summary(m)
            c_se <- model_summary$coefficients[, "Std. Error"]
          } else {
            c_se <- tryCatch({
              sqrt(diag(vcov(m)))
            }, error = function(e) {
              summary(m)$coefficients[, "Std. Error"]
            })
          }
        }
        data.frame(term = names(c_est), estimate = c_est, std.error = c_se, stringsAsFactors = FALSE)
      })

      coefs <- Filter(function(x) !is.null(x), coefs)
      all_terms <- unique(unlist(lapply(coefs, function(df) df$term)))
      pooled_results <- data.frame(term = all_terms, estimate = NA_real_, std.error = NA_real_, stringsAsFactors = FALSE)

      for (term in all_terms) {
        term_ests <- sapply(coefs, function(df) if (term %in% df$term) df$estimate[df$term == term] else NA_real_)
        term_ests <- term_ests[!is.na(term_ests)]
        term_ses  <- sapply(coefs, function(df) if (term %in% df$term) df$std.error[df$term == term] else NA_real_)
        term_ses  <- term_ses[!is.na(term_ses)]
        if (length(term_ests) < 2) next
        pooled_est <- mean(term_ests)
        w_var <- mean(term_ses^2)
        b_var <- sum((term_ests - pooled_est)^2) / (length(term_ests) - 1)
        total_var <- w_var + (1 + 1/length(term_ests)) * b_var
        pooled_results$estimate[pooled_results$term == term] <- pooled_est
        pooled_results$std.error[pooled_results$term == term] <- sqrt(total_var)
      }

      Results_multivariate_analysis <- pooled_results %>%
        mutate(
          `2.5 %` = estimate - 1.96 * std.error,
          `97.5 %` = estimate + 1.96 * std.error,
          p.value = 2 * pnorm(-abs(estimate / std.error))
        )

    } else {
      # ---- Non-random effects: classical glm/glm.nb, pooled with mice::pool ----
      res_comb <- vector("list", length(actual_imps))
      for (i in seq_along(actual_imps)) {
        data_subset <- implist[[i]]
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
          warning(sprintf("Model failed to fit for imputation %s: %s", actual_imps[i], e$message))
          NULL
        })
      }

      current_models_list <- res_comb

      ok_models <- Filter(function(x) !is.null(x), res_comb)
      if (length(ok_models) == 0) {
        stop("All models failed to fit in MI_estimates.")
      }

      pooled <- mice::pool(ok_models)
      Results_multivariate_analysis <- summary(pooled, conf.int = TRUE, exp = FALSE)
    }

    # --- Fit stats & D1 components from current_models_list ---
    if (!is.null(current_models_list) && length(current_models_list) > 0) {

      fit_list <- lapply(current_models_list,
                         compute_fit_stats,
                         outcome_var = outcome_var,
                         model_type  = model_type)

      fit_list <- Filter(function(x) {
        !is.null(x) && is.data.frame(x) && nrow(x) > 0
      }, fit_list)

      if (length(fit_list) > 0) {
        fit_stats_per_imp <- dplyr::bind_rows(fit_list, .id = "imp")
        fit_stats_per_imp$imp <- actual_imps[as.integer(fit_stats_per_imp$imp)]

        pooled_vals <- lapply(
          fit_stats_per_imp[, setdiff(names(fit_stats_per_imp), "imp"), drop = FALSE],
          function(x) mean(x, na.rm = TRUE)
        )
        fit_stats_pooled <- data.frame(imp = "pooled", pooled_vals, check.names = FALSE)

        if (all(c("logLik", "df") %in% names(fit_stats_per_imp))) {
          D1_components <- fit_stats_per_imp[, c("imp", "logLik", "df")]
        }
      } else {
        fit_stats_per_imp <- NULL
        fit_stats_pooled  <- NULL
        D1_components     <- NULL
      }
    }

    # exponentiation
    Results_multivariate_analysis <- Results_multivariate_analysis %>%
      mutate(exp_estimate   = exp(estimate),
             exp_CI95_lower = exp(`2.5 %`),
             exp_CI95_upper = exp(`97.5 %`)) %>%
      dplyr::select(term, estimate, std.error, `2.5 %`, `97.5 %`,
                    exp_estimate, exp_CI95_lower, exp_CI95_upper, p.value)

    # Normalize poly() term name a bit
    Results_multivariate_analysis$term <- gsub(
      "poly\\(([^,]+), degree = ([0-9]+), raw = TRUE\\)",
      "poly(\\1, \\2, raw = TRUE)",
      Results_multivariate_analysis$term
    )

    # Filter trial and nuisance terms
    if (!is.null(trial_col)) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        filter(!grepl(trial_col, term))
    }
    if (use_random_effects) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        filter(!grepl(paste0("\\(Intercept\\)|SD\\(", random_intercept_var, "\\)"), term))
    }
    if (use_strata && model_type != "cox") {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        filter(!grepl(paste0("^as\\.factor\\(", stratified_intercept_var, "\\)"), term))
    }

    # --- Rename spline terms nicely ( *_rcs_linear / *_rcs_nl1 / *_rcs_nl2 ) ---
    if (!is.null(spline_terms) && length(spline_terms) > 0) {
      for (v in spline_terms) {
        pattern_base <- paste0("rcs\\(", v, "[^)]*\\)")
        idx <- grepl(pattern_base, Results_multivariate_analysis$term)
        if (any(idx)) {
          old_names <- Results_multivariate_analysis$term[idx]
          new_names <- old_names

          two_p <- grepl("''$", old_names)
          one_p <- grepl("'$", old_names) & !two_p

          new_names[two_p] <- paste0(v, "_rcs_nl2")
          new_names[one_p] <- paste0(v, "_rcs_nl1")
          new_names[!one_p & !two_p] <- paste0(v, "_rcs_linear")

          Results_multivariate_analysis$term[idx] <- new_names
        }
      }
    }

    # --- Remove polynomial terms if requested (include_poly_terms = FALSE) ---
    if (length(poly_terms_detected) > 0 && !include_poly_terms) {
      poly_patterns <- sapply(poly_terms_detected, function(x) {
        if (grepl("poly\\(", x)) {
          var_name <- trimws(gsub("poly\\(([^,]+),.*", "\\1", x))
          return(paste0("^poly\\(", var_name, ".*\\)"))
        } else ""
      })
      pattern <- paste(poly_patterns, collapse = "|")
      if (nzchar(pattern)) {
        Results_multivariate_analysis <- Results_multivariate_analysis %>%
          filter(!grepl(pattern, term))
      }
    }

    # null model row
    if (is_null_model) {
      null_row <- data.frame(
        term = "No_predictor",
        estimate = NA, std.error = NA, `2.5 %` = NA, `97.5 %` = NA,
        exp_estimate = NA, exp_CI95_lower = NA, exp_CI95_upper = NA, p.value = NA,
        stringsAsFactors = FALSE, check.names = FALSE
      )
      Results_multivariate_analysis <- rbind(null_row, Results_multivariate_analysis)
    }

    # flags
    if (highlight_interactions && length(interaction_terms) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_interaction = sapply(term, function(t) any(sapply(interaction_terms, function(i) grepl(i, t, fixed = TRUE)))))
    } else {
      Results_multivariate_analysis$is_interaction <- FALSE
    }

    # spline flag
    Results_multivariate_analysis$is_spline <- FALSE
    if (!is.null(spline_terms) && length(spline_terms) > 0) {
      Results_multivariate_analysis$is_spline[grepl("_rcs_", Results_multivariate_analysis$term)] <- TRUE
    }

    # polynomial flag
    Results_multivariate_analysis$is_polynomial <- grepl("poly\\(", Results_multivariate_analysis$term)

    # attributes
    attr(Results_multivariate_analysis, "has_random_effects")      <- use_random_effects
    attr(Results_multivariate_analysis, "random_intercept_var")    <- random_intercept_var
    attr(Results_multivariate_analysis, "predictor_tested")        <- current_predictor
    attr(Results_multivariate_analysis, "model_type_used")         <- model_type
    attr(Results_multivariate_analysis, "formula")                 <- formula_string_current
    attr(Results_multivariate_analysis, "has_interactions")        <- length(interaction_terms) > 0
    attr(Results_multivariate_analysis, "interaction_terms")       <- if (length(interaction_terms) > 0) interaction_terms else NULL
    attr(Results_multivariate_analysis, "has_splines")             <- length(spline_terms_detected) > 0
    attr(Results_multivariate_analysis, "spline_terms")            <- if (length(spline_terms_detected) > 0) spline_terms_detected else NULL
    attr(Results_multivariate_analysis, "has_polynomials")         <- length(poly_terms_detected) > 0
    attr(Results_multivariate_analysis, "polynomial_terms")        <- if (length(poly_terms_detected) > 0) poly_terms_detected else NULL
    attr(Results_multivariate_analysis, "imputations")             <- actual_imps
    attr(Results_multivariate_analysis, "n_imp")                   <- length(actual_imps)
    attr(Results_multivariate_analysis, "stratified_intercept_var")<- stratified_intercept_var

    # variance components / ICC
    if (use_random_effects && !is.null(current_models_list) && length(current_models_list) > 0) {
      model1 <- current_models_list[[1]]
      if (inherits(model1, "glmmTMB")) {
        re_var <- tryCatch({ as.data.frame(glmmTMB::VarCorr(model1)$cond)[1, "vcov"] }, error = function(e) NA)
        attr(Results_multivariate_analysis, "variance_components") <- re_var
      } else if (model_type == "cox" && !is.null(model1$vcoef)) {
        var_comp <- as.numeric(model1$vcoef)
        names(var_comp) <- "Var(Random Intercept)"
        attr(Results_multivariate_analysis, "variance_components") <- var_comp
      } else if (inherits(model1, "merMod")) {
        var_comp <- as.data.frame(lme4::VarCorr(model1))
        var_comp <- setNames(var_comp$vcov, paste0("Var(", var_comp$grp, ")"))
        attr(Results_multivariate_analysis, "variance_components") <- var_comp
        if (model_type == "lm") {
          tau2   <- var_comp[1]
          sigma2 <- attr(lme4::VarCorr(model1), "sc")^2
          ICC    <- tau2 / (tau2 + sigma2)
          attr(Results_multivariate_analysis, "ICC") <- ICC
        }
      }
    }

    # attach fit stats / D1 components
    attr(Results_multivariate_analysis, "fit_stats_per_imp") <- fit_stats_per_imp
    attr(Results_multivariate_analysis, "fit_stats_pooled")  <- fit_stats_pooled
    attr(Results_multivariate_analysis, "D1_components")     <- D1_components

    Results_multivariate_analysis
  }

  ## ---- Main loop ----
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
    create_combined_results <- function() {
      combined_results <- data.frame()
      for (pred_name in names(results_list)) {
        result <- results_list[[pred_name]]
        if (is.null(result) || nrow(result) == 0) next
        if (pred_name == "No_predictor") {
          predictor_rows <- result[result$term == "No_predictor", ]
        } else {
          predictor_rows <- result[grepl(paste0("^", pred_name), result$term) & !grepl("Intercept", result$term), ]
          if (nrow(predictor_rows) == 0 && !is.null(covariables) && length(covariables) > 0) {
            covariable_pattern <- paste(covariables, collapse = "|")
            predictor_rows <- result[grepl(pred_name, result$term) & !grepl("Intercept", result$term), ]
            predictor_rows <- predictor_rows[!grepl(covariable_pattern, predictor_rows$term), ]
          }
        }
        if (nrow(predictor_rows) > 0) {
          predictor_rows <- predictor_rows[, c("term","estimate","std.error","2.5 %","97.5 %","exp_estimate","exp_CI95_lower","exp_CI95_upper","p.value")]
          combined_results <- rbind(combined_results, predictor_rows)
        }
      }
      rownames(combined_results) <- NULL
      combined_results
    }
    combined_predictors_results <- create_combined_results()

    class(results_list) <- c("MI_estimates_multi", "list")
    attr(results_list, "covariables")      <- covariables
    attr(results_list, "predictors")       <- predictor_vars
    attr(results_list, "model_type")       <- model_type
    attr(results_list, "n_predictors")     <- length(predictor_vars)
    attr(results_list, "combined_results") <- combined_predictors_results
    return(results_list)
  }
}
