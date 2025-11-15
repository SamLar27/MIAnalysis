#' Compute Estimates for Multiple Imputed Data
#'
#' This function fits statistical models to multiply imputed datasets and pools the results using Rubin's Rules.
#' It supports various regression models, including negative binomial, logistic, Poisson, linear, and Cox regression.
#' The function also supports mixed models (random effects) using `glmmTMB`, `lme4`, and `coxme`, and can handle interaction terms,
#' restricted cubic splines (rcs), and polynomial terms.
#'
#' Intercept / slope structure is inferred automatically:
#'   - If `stratified_intercept_var` is provided:
#'       * GLM: fixed stratified intercept via `0 + as.factor(stratified_intercept_var)`
#'       * Cox: stratified baseline hazard via `strata(stratified_intercept_var)`
#'   - If `random_intercept_var` is provided:
#'       * Mixed model is used (glmmTMB/lme4/coxme)
#'       * Random intercept and/or random slopes are defined by `predictor_vars_random_slope`
#'         and `covariables_random_slope`, with the following special case:
#'           - If `random_intercept_var == stratified_intercept_var`:
#'               + and no random slopes -> no random effect at all (pure stratified intercept)
#'               + and random slopes -> random slopes only (no random intercept):
#'                     (0 + slopes | random_intercept_var)
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for GLM models.
#' @param predictor_vars A vector of predictor variables to be tested in the model. Each predictor will be included individually along with all covariables.
#' @param covariables A vector of covariable names to be included as fixed effects in all models. These are always included regardless of which predictor is being tested.
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param imp_n Number of imputations (if NULL, it is detected automatically).
#' @param model_type The model type to be used:
#'   \itemize{
#'     \item "nb" for Negative Binomial regression (using `MASS::glm.nb` or `glmmTMB` if random effects are used)
#'     \item "lm" for Linear regression
#'     \item "bin" for Logistic regression (binomial family)
#'     \item "poisson" for Poisson regression
#'     \item "gamma" for Gamma regression
#'     \item "quasipoisson" for Overdispersed Poisson regression
#'     \item "quasibinomial" for Overdispersed logistic regression
#'     \item "cox" for Cox Proportional Hazards regression (using `survival` or `coxme`)
#'   }
#' @param followup_offset Whether to include an offset for follow-up duration ("Yes" or "No").
#' @param followup_col The column representing follow-up duration (required if `followup_offset = "Yes"`).
#' @param trial_factor Whether to include trial as a fixed factor ("Yes" or "No").
#' @param trial_col The column used for trial factor (required if `trial_factor = "Yes"`).
#' @param time_col The time variable for Cox regression (only required if `model_type = "cox"`).
#' @param event_col The event variable for Cox regression (only required if `model_type = "cox"`).
#' @param formula_string Optional custom formula string. If provided, overrides the formula created from predictor_vars and covariables.
#' @param highlight_interactions Whether to flag interaction terms in the output (default: TRUE).
#' @param spline_terms A named list to specify variables to be modeled with restricted cubic splines. Each element should be a list with 'var' (variable name) and 'knots' (number of knots or explicit knot positions).
#' @param poly_terms A named list to specify variables to be modeled with polynomial terms. Each element should be a list with 'var' (variable name) and 'degree' (2 for quadratic, 3 for cubic).
#' @param include_spline_terms Whether to include individual spline components in the output (default: FALSE).
#' @param include_poly_terms Whether to include individual polynomial components in the output (default: TRUE).
#' @param random_intercept_var Grouping variable for random effects.
#'        If NULL, no random effects are included.
#' @param predictor_vars_random_slope Character vector of predictor variable names for which a random slope should be added.
#' @param covariables_random_slope Character vector of covariable names for which a random slope should be added.
#' @param stratified_intercept_var Variable whose levels define a stratified intercept (GLM) or baseline hazard (Cox).
#'        If NULL, no stratification is applied.
#'
#' @return A list (or data frame if a single predictor) containing model results.
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
                         # --- splines ---
                         spline_terms = NULL,              # character vector of variables
                         spline_knots_n = NULL,            # number of knots (optional)
                         spline_knots_percentile = NULL,   # percentiles (0–100) for knots (optional)
                         # --- polynomials (kept as before) ---
                         poly_terms = NULL,
                         include_poly_terms = TRUE,
                         # --- random effects ---
                         random_intercept_var = NULL,
                         predictor_vars_random_slope = NULL,
                         covariables_random_slope = NULL,
                         # --- stratified intercept / baseline ---
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

  ## ---- Spline support (rms::rcs) ----
  use_rms <- FALSE
  if (!is.null(spline_terms)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("Package 'rms' is required when 'spline_terms' is specified (for rcs()).")
    } else {
      use_rms <- TRUE
      require(rms)
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

  ## ---- Validate splines / poly ----
  if (!is.null(spline_terms)) {
    if (!is.character(spline_terms)) stop("spline_terms must be a character vector of variable names.")
    for (v in spline_terms) {
      if (!v %in% names(data)) stop(paste("Spline variable", v, "not found in data"))
    }
  }

  if (!is.null(poly_terms)) {
    if (!is.list(poly_terms)) stop("poly_terms must be a list")
    for (i in seq_along(poly_terms)) {
      if (!is.list(poly_terms[[i]]) || is.null(poly_terms[[i]]$var))
        stop("Each element in poly_terms must be a list with at least a 'var' element")
      if (!poly_terms[[i]]$var %in% names(data))
        stop(paste("Polynomial variable", poly_terms[[i]]$var, "not found in data"))
      if (is.null(poly_terms[[i]]$degree)) {
        poly_terms[[i]]$degree <- 2
      } else if (!is.numeric(poly_terms[[i]]$degree) || !(poly_terms[[i]]$degree %in% c(2, 3))) {
        stop("Degree in poly_terms must be either 2 (quadratic) or 3 (cubic)")
      }
    }
  }

  ## ---- Predictor / covariable presence checks (if no custom formula) ----
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
    if (grepl("^rcs\\(|^bs\\(", term)) return(term)
    if (grepl("\\*", term)) {
      vars_split <- trimws(unlist(strsplit(term, "\\*")))
      all_combinations <- lapply(seq_along(vars_split), function(k) {
        combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
      })
      unlist(all_combinations)
    } else term
  }

  ## ---- Build spline / poly maps ----
  process_special_terms <- function() {

    # map: variable -> "rcs(var, ...)"
    spline_map <- list()
    if (!is.null(spline_terms)) {
      for (var_name in spline_terms) {

        x <- data[[var_name]]

        if (!is.null(spline_knots_percentile)) {
          probs <- spline_knots_percentile / 100
          knots <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, type = 2))
          knots <- sort(unique(knots))
          knot_str <- paste(knots, collapse = ", ")
          spline_map[[var_name]] <- paste0("rcs(", var_name, ", c(", knot_str, "))")
        } else if (!is.null(spline_knots_n)) {
          spline_map[[var_name]] <- paste0("rcs(", var_name, ", ", spline_knots_n, ")")
        } else {
          spline_map[[var_name]] <- paste0("rcs(", var_name, ", 4)")
        }
      }
    }

    # map: var -> "poly(var, degree = d, raw = TRUE)"
    poly_map <- list()
    if (!is.null(poly_terms)) {
      for (i in seq_along(poly_terms)) {
        var_name <- poly_terms[[i]]$var
        degree   <- poly_terms[[i]]$degree
        poly_map[[var_name]] <- paste0("poly(", var_name, ", degree = ", degree, ", raw = TRUE)")
      }
    }

    list(spline_map = spline_map, poly_map = poly_map)
  }

  special_terms <- process_special_terms()
  spline_map <- special_terms$spline_map
  poly_map   <- special_terms$poly_map

  ## ---- Random-effects core builder ----
  build_random_term_core <- function(rs_vars) {
    if (!use_random_effects) return("")
    rs_vars <- unique(rs_vars)
    same_group_as_strata <- (use_strata &&
                               stratified_intercept_var == random_intercept_var)

    if (length(rs_vars) == 0) {
      if (same_group_as_strata) {
        return("")  # pure stratified intercept, no random effect
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

  ## ---- Build formula for each predictor ----
  build_formula_for_predictor <- function(current_predictor) {

    trial_term  <- if (trial_factor == "Yes") paste("+ as.factor(", trial_col, ")", sep = "") else ""
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

      # covariables for this model (remove the predictor itself if present)
      if (is.null(covariables)) {
        current_covariables <- character(0)
      } else {
        current_covariables <- covariables
        if (!is.null(current_predictor) && nzchar(current_predictor)) {
          current_covariables <- setdiff(current_covariables, current_predictor)
        }
      }

      # build term for predictor (if any)
      predictor_term <- NULL
      if (!is.null(current_predictor) && nzchar(current_predictor)) {
        if (current_predictor %in% names(spline_map)) {
          predictor_term <- spline_map[[current_predictor]]
        } else if (current_predictor %in% names(poly_map)) {
          predictor_term <- poly_map[[current_predictor]]
        } else {
          predictor_term <- current_predictor
        }
      }

      # build terms for covariables
      cov_terms <- character(0)
      if (length(current_covariables) > 0) {
        for (cv in current_covariables) {
          if (cv %in% names(spline_map)) {
            cov_terms <- c(cov_terms, spline_map[[cv]])
          } else if (cv %in% names(poly_map)) {
            cov_terms <- c(cov_terms, poly_map[[cv]])
          } else {
            cov_terms <- c(cov_terms, cv)
          }
        }
      }

      # interactions expansion (only for non-spline/poly terms)
      expand_vec <- function(v) {
        unlist(lapply(v, expand_terms))
      }

      fixed_terms <- c(
        if (!is.null(predictor_term)) expand_vec(predictor_term) else NULL,
        expand_vec(cov_terms)
      )
      fixed_terms <- fixed_terms[nzchar(fixed_terms)]

      covariables_in_model <- current_covariables
      random_term <- build_random_term(current_predictor, covariables_in_model)

      if (model_type == "cox") {
        if (is.null(time_col) || is.null(event_col)) stop("For Cox regression, time_col and event_col must be provided.")
        formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ ",
                             paste(fixed_terms, collapse = " + "), trial_term, strat_cox_piece, random_term, sep = "")
      } else {
        rhs <- paste(fixed_terms, collapse = " + ")
        formula_str <- paste(outcome_var, "~", rhs, offset_term, trial_term, strat_glm_piece, random_term)
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

  ## ---- Helper: rename rcs() terms to clean names ----
  rename_rcs_terms <- function(df) {
    pattern <- "rcs\\(([^,]+),[^)]*\\)(.*)$"
    is_rcs  <- grepl(pattern, df$term)

    if (!any(is_rcs)) return(df)

    var_name <- gsub(pattern, "\\1", df$term[is_rcs])
    suffix   <- gsub(pattern, "\\2", df$term[is_rcs])

    new_names <- mapply(function(v, s) {
      s_trim <- trimws(s)
      if (s_trim == "")  return(paste0(v, "_rcs_linear"))
      if (s_trim == "'") return(paste0(v, "_rcs_nl1"))
      if (s_trim == "''")return(paste0(v, "_rcs_nl2"))
      paste0(v, "_rcs", s_trim)
    }, var_name, suffix, USE.NAMES = FALSE)

    df$term[is_rcs] <- new_names
    df
  }

  ## ---- Fit one predictor ----
  fit_model_for_predictor <- function(current_predictor) {

    formula_string_current <- build_formula_for_predictor(current_predictor)
    model_formula <- as.formula(formula_string_current)

    # interaction / spline / poly detection on raw formula
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

    if (use_random_effects) {

      models_list <- vector("list", imp_n)
      for (i in seq_along(actual_imps)) {
        data_i <- implist[[i]]
        if (model_type == "lm") {
          models_list[[i]] <- lme4::lmer(model_formula, data = data_i)
        } else if (model_type %in% c("bin")) {
          models_list[[i]] <- lme4::glmer(model_formula, family = binomial(), data = data_i)
        } else if (model_type %in% c("poisson", "quasipoisson")) {
          models_list[[i]] <- lme4::glmer(model_formula, family = poisson(), data = data_i)
        } else if (model_type == "nb") {
          if (requireNamespace("glmmTMB", quietly = TRUE)) {
            models_list[[i]] <- glmmTMB::glmmTMB(model_formula, family = nbinom2, data = data_i)
          } else {
            warning("Using Poisson instead of Negative Binomial. Install glmmTMB for proper NB mixed models.")
            models_list[[i]] <- lme4::glmer(model_formula, family = poisson(), data = data_i)
          }
        } else if (model_type == "cox") {
          models_list[[i]] <- coxme::coxme(model_formula, data = data_i)
        } else if (model_type == "gamma") {
          if (requireNamespace("glmmTMB", quietly = TRUE)) {
            models_list[[i]] <- glmmTMB::glmmTMB(model_formula, family = Gamma, data = data_i)
          } else {
            warning("Gamma family with random effects may not be stable in lme4.")
            models_list[[i]] <- lme4::glmer(model_formula, family = Gamma, data = data_i)
          }
        } else stop("Unsupported model type for random effects.")
      }

      coefs <- lapply(models_list, function(m) {
        if (model_type == "cox") {
          c_est <- fixef(m)
          c_se  <- sqrt(diag(vcov(m)))
        } else if (inherits(m, "glmmTMB")) {
          c_est <- fixef(m)$cond
          v <- vcov(m)$cond
          c_se <- sqrt(diag(v))
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
      current_models_list <- models_list

    } else {

      res_comb <- vector("list", length(actual_imps))
      for (i in seq_along(actual_imps)) {
        data_subset <- implist[[i]]
        model <- switch(model_type,
                        "nb" = MASS::glm.nb(model_formula, data = data_subset),
                        "lm" = glm(model_formula, family = gaussian(), data = data_subset),
                        "bin" = glm(model_formula, family = binomial(), data = data_subset),
                        "poisson" = glm(model_formula, family = poisson(), data = data_subset),
                        "gamma" = glm(model_formula, family = Gamma(), data = data_subset),
                        "quasipoisson" = glm(model_formula, family = quasipoisson(), data = data_subset),
                        "quasibinomial" = glm(model_formula, family = quasibinomial(), data = data_subset),
                        "cox" = coxph(model_formula, data = data_subset),
                        stop("Unsupported model type."))
        res_comb[[i]] <- model
      }
      pooled <- mice::pool(res_comb)
      Results_multivariate_analysis <- summary(pooled, conf.int = TRUE, exp = FALSE)
      current_models_list <- res_comb
    }

    # save spline flag *before* renaming
    is_spline_logical <- if (length(spline_terms_detected) > 0) {
      grepl("rcs\\(|bs\\(", Results_multivariate_analysis$term)
    } else {
      rep(FALSE, nrow(Results_multivariate_analysis))
    }

    # exponentiation
    Results_multivariate_analysis <- Results_multivariate_analysis %>%
      mutate(exp_estimate   = exp(estimate),
             exp_CI95_lower = exp(`2.5 %`),
             exp_CI95_upper = exp(`97.5 %`)) %>%
      dplyr::select(term, estimate, std.error, `2.5 %`, `97.5 %`,
                    exp_estimate, exp_CI95_lower, exp_CI95_upper, p.value)

    Results_multivariate_analysis$term <- gsub("poly\\(([^,]+), degree = ([0-9]+), raw = TRUE\\)",
                                               "poly(\\1, \\2, raw = TRUE)", Results_multivariate_analysis$term)

    # Filter trial and nuisance random-effects terms
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

    # optionally filter polynomial inner components
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
        is_spline_logical <- is_spline_logical[match(Results_multivariate_analysis$term,
                                                     Results_multivariate_analysis$term)]
      }
    }

    # ---- NEW: rename rcs() terms to nice names ----
    Results_multivariate_analysis <- rename_rcs_terms(Results_multivariate_analysis)

    # null model row
    if (is_null_model) {
      null_row <- data.frame(
        term = "No_predictor",
        estimate = NA, std.error = NA, `2.5 %` = NA, `97.5 %` = NA,
        exp_estimate = NA, exp_CI95_lower = NA, exp_CI95_upper = NA, p.value = NA,
        stringsAsFactors = FALSE, check.names = FALSE
      )
      Results_multivariate_analysis <- rbind(null_row, Results_multivariate_analysis)
      is_spline_logical <- c(FALSE, is_spline_logical)
    }

    # flags
    if (highlight_interactions && length(interaction_terms) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_interaction = sapply(term, function(t) any(sapply(interaction_terms, function(i) grepl(i, t, fixed = TRUE)))))
    }
    if (length(spline_terms_detected) > 0) {
      Results_multivariate_analysis$is_spline <- is_spline_logical
    }
    if (highlight_interactions && length(poly_terms_detected) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_polynomial = grepl("^poly\\(", term))
    }

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

    if (use_random_effects && exists("current_models_list") && length(current_models_list) > 0) {
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

    Results_multivariate_analysis
  }

  ## ---- Main loop over predictor_vars ----
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
