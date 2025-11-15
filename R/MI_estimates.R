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
                         # --- NEW INTERFACE FOR SPLINES ---
                         # You can pass either:
                         #   - character vector of variable names (e.g. "Tif_preBD_PCT_0W")
                         #   - or the old list-of-lists structure
                         spline_terms = NULL,
                         # number of knots (used if spline_knots_percentile is NULL)
                         spline_knots_n = 4,
                         # percentiles (0-100) to define knot locations, applied on the whole data
                         spline_knots_percentile = NULL,
                         # ---------------------------------
                         poly_terms = NULL,
                         include_spline_terms = FALSE,
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

  # Spline support
  if (!is.null(spline_terms)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      message("Package 'rms' is not available. Using splines::bs() instead of rcs().")
      use_rms <- FALSE
    } else {
      require(rms)
      use_rms <- TRUE
    }
  } else {
    use_rms <- FALSE
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

  ## ---- NEW: user-friendly rcs specification -> internal structure ----
  # Old API expected: list(list(var="x", knots=...), list(var="y", knots=...))
  # New API: spline_terms = "x" or c("x","y"), plus spline_knots_n / spline_knots_percentile.
  if (!is.null(spline_terms) && is.character(spline_terms)) {
    spline_vars <- spline_terms
    spline_list <- vector("list", length(spline_vars))

    # Convert percentiles to probabilities if provided
    if (!is.null(spline_knots_percentile)) {
      probs <- spline_knots_percentile / 100
      if (any(probs <= 0 | probs >= 1)) {
        stop("spline_knots_percentile must be strictly between 0 and 100 (exclusive).")
      }
    } else {
      probs <- NULL
    }

    for (i in seq_along(spline_vars)) {
      v <- spline_vars[i]
      if (!v %in% names(data)) stop(paste("Spline variable", v, "not found in data."))

      x <- data[[v]]

      if (!is.null(probs)) {
        # Global percentile-based knots over the *entire* distribution
        knots <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE))
      } else {
        # Fall back to "number of knots" logic used previously
        # We just store "k" and let process_special_terms() translate to df/knots
        knots <- spline_knots_n
      }

      spline_list[[i]] <- list(var = v, knots = knots)
    }

    spline_terms <- spline_list
  }

  ## ---- Validate spline/poly (internal representation) ----
  if (!is.null(spline_terms)) {
    if (!is.list(spline_terms)) stop("spline_terms must be a list (or character vector of variable names).")
    for (i in seq_along(spline_terms)) {
      if (!is.list(spline_terms[[i]]) || is.null(spline_terms[[i]]$var))
        stop("Each element in spline_terms must be a list with at least a 'var' element.")
      if (!spline_terms[[i]]$var %in% names(data))
        stop(paste("Spline variable", spline_terms[[i]]$var, "not found in data."))
      if (is.null(spline_terms[[i]]$knots)) {
        spline_terms[[i]]$knots <- spline_knots_n
      }
      # knots can be numeric scalar (k) or numeric vector (explicit positions)
      if (!is.numeric(spline_terms[[i]]$knots)) {
        stop("Knots must be numeric (either number of knots or vector of knot positions).")
      }
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
        stop("Degree must be either 2 (quadratic) or 3 (cubic)")
      }
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
    if (grepl("^rcs\\(|^bs\\(", term)) return(term)
    if (grepl("\\*", term)) {
      vars_split <- trimws(unlist(strsplit(term, "\\*")))
      all_combinations <- lapply(seq_along(vars_split), function(k) {
        combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
      })
      unlist(all_combinations)
    } else term
  }

  process_special_terms <- function() {
    spline_formula_parts <- character(0)
    if (!is.null(spline_terms)) {
      for (i in seq_along(spline_terms)) {
        var_name <- spline_terms[[i]]$var
        knots    <- spline_terms[[i]]$knots

        # Two modes:
        #  - scalar "knots" = number of knots (k) -> rcs(x, k) or bs(x, df = k+2)
        #  - vector "knots" = explicit knot positions -> rcs(x, c(...)) or bs(x, knots = c(...))
        if (length(knots) == 1) {
          k <- knots[1]
          if (use_rms) {
            spline_formula_parts <- c(spline_formula_parts,
                                      paste0("rcs(", var_name, ", ", k, ")"))
          } else {
            df <- k + 2
            spline_formula_parts <- c(spline_formula_parts,
                                      paste0("bs(", var_name, ", df = ", df, ", degree = 3)"))
          }
        } else {
          # explicit knot locations (from percentiles or user-provided)
          if (use_rms) {
            spline_formula_parts <- c(
              spline_formula_parts,
              paste0("rcs(", var_name, ", c(", paste(knots, collapse = ", "), "))")
            )
          } else {
            spline_formula_parts <- c(
              spline_formula_parts,
              paste0("bs(", var_name, ", knots = c(", paste(knots, collapse = ", "), "), degree = 3)")
            )
          }
        }
      }
    }

    poly_formula_parts <- character(0)
    if (!is.null(poly_terms)) {
      for (i in seq_along(poly_terms)) {
        var_name <- poly_terms[[i]]$var
        degree   <- poly_terms[[i]]$degree
        poly_formula_parts <- c(poly_formula_parts,
                                paste0("poly(", var_name, ", degree = ", degree, ", raw = TRUE)"))
      }
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

  ## ---- Formula builder ----
  build_formula_for_predictor <- function(current_predictor) {
    trial_term  <- if (trial_factor == "Yes") paste("+ as.factor(", trial_col, ")", sep = "") else ""
    offset_term <- if (followup_offset == "Yes") paste("+ offset(log(", followup_col, "))", sep = "") else ""

    strat_glm_piece  <- ""
    strat_cox_piece  <- ""
    intercept_prefix <- ""   # keep standard intercept

    if (use_strata) {
      if (model_type == "cox") {
        strat_cox_piece <- paste("+ strata(", stratified_intercept_var, ")", sep = "")
      } else {
        # GLM: add as.factor(strata) but keep intercept
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

      if (current_predictor == "" || is.null(current_predictor) || is.na(current_predictor)) {
        all_terms <- c(expanded_covariables, spline_formula_parts, poly_formula_parts)
      } else {
        current_spline_parts <- spline_formula_parts
        current_poly_parts   <- poly_formula_parts

        # Remove own spline/poly if present (to avoid double inclusion)
        if (!is.null(spline_terms)) {
          for (i in seq_along(spline_terms)) {
            if (spline_terms[[i]]$var == current_predictor) {
              pattern <- if (use_rms) paste0("rcs\\(", current_predictor, ",") else paste0("bs\\(", current_predictor, ",")
              current_spline_parts <- current_spline_parts[!grepl(pattern, current_spline_parts)]
            }
          }
        }
        if (!is.null(poly_terms)) {
          for (i in seq_along(poly_terms)) {
            if (poly_terms[[i]]$var == current_predictor) {
              pattern <- paste0("poly\\(", current_predictor, ",")
              current_poly_parts <- current_poly_parts[!grepl(pattern, current_poly_parts)]
            }
          }
        }

        current_covariables   <- if (is.null(covariables)) character(0) else covariables[covariables != current_predictor]
        covariables_in_model  <- current_covariables
        expanded_covariables  <- unique(unlist(lapply(current_covariables, expand_terms)))
        expanded_predictor    <- expand_terms(current_predictor)
        all_terms <- c(expanded_predictor, expanded_covariables,
                       current_spline_parts, current_poly_parts)
      }

      random_term <- build_random_term(current_predictor, covariables_in_model)

      if (model_type == "cox") {
        if (is.null(time_col) || is.null(event_col)) stop("For Cox regression, time_col and event_col must be provided.")
        formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ ",
                             paste(all_terms, collapse = " + "), trial_term, strat_cox_piece, random_term, sep = "")
      } else {
        fixed_rhs <- paste(c(if (intercept_prefix != "") intercept_prefix else NULL,
                             paste(all_terms, collapse = " + ")), collapse = " ")
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

      all_terms <- unique(unlist(lapply(coefs, function(df) df$term)))
      pooled_results <- data.frame(term = all_terms, estimate = NA_real_, std.error = NA_real_, stringsAsFactors = FALSE)

      for (term in all_terms) {
        term_ests <- sapply(coefs, function(df) if (term %in% df$term) df$estimate[df$term == term] else NA_real_)
        term_ests <- term_ests[!is.na(term_ests)]
        term_ses  <- sapply(coefs, function(df) if (term %in% df$term) df$std.error[df$term == term] else NA_real_)
        term_ses  <- term_ses[!is.na(term_ses)]
        if (length(term_ests) < 2) next
        pooled_est <- mean(term_ests)
        w_var      <- mean(term_ses^2)
        b_var      <- sum((term_ests - pooled_est)^2) / (length(term_ests) - 1)
        total_var  <- w_var + (1 + 1/length(term_ests)) * b_var
        pooled_results$estimate[pooled_results$term == term] <- pooled_est
        pooled_results$std.error[pooled_results$term == term] <- sqrt(total_var)
      }

      Results_multivariate_analysis <- pooled_results %>%
        mutate(
          `2.5 %` = estimate - 1.96 * std.error,
          `97.5 %` = estimate + 1.96 * std.error,
          p.value  = 2 * pnorm(-abs(estimate / std.error))
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

    # exponentiation
    Results_multivariate_analysis <- Results_multivariate_analysis %>%
      mutate(exp_estimate   = exp(estimate),
             exp_CI95_lower = exp(`2.5 %`),
             exp_CI95_upper = exp(`97.5 %`)) %>%
      dplyr::select(term, estimate, std.error, `2.5 %`, `97.5 %`,
                    exp_estimate, exp_CI95_lower, exp_CI95_upper, p.value)

    Results_multivariate_analysis$term <- gsub("poly\\(([^,]+), degree = ([0-9]+), raw = TRUE\\)",
                                               "poly(\\1, \\2, raw = TRUE)", Results_multivariate_analysis$term)

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

    # filter spline / poly components if required
    if (length(spline_terms_detected) > 0 && !include_spline_terms) {
      spline_patterns <- sapply(spline_terms_detected, function(x) {
        if (grepl("rcs\\(", x)) {
          var_name <- gsub("rcs\\(([^,]+),.*", "\\1", x)
        } else if (grepl("bs\\(", x)) {
          var_name <- gsub("bs\\(([^,]+),.*", "\\1", x)
        } else return("")
        paste0("^", trimws(var_name), "'")
      })
      pattern <- paste(spline_patterns, collapse = "|")
      if (nzchar(pattern)) {
        Results_multivariate_analysis <- Results_multivariate_analysis %>%
          filter(!grepl(pattern, term))
      }
    }
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
    }
    if (highlight_interactions && length(spline_terms_detected) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_spline = grepl("rcs\\(|bs\\(", term))
    }
    if (highlight_interactions && length(poly_terms_detected) > 0) {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        mutate(is_polynomial = grepl("poly\\(", term))
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
          predictor_rows <- predictor_rows[, c("term","estimate","std.error","2.5 %","97.5 %",
                                               "exp_estimate","exp_CI95_lower","exp_CI95_upper","p.value")]
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
