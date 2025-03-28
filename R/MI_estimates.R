#' Compute Estimates for Multiple Imputed Data
#'
#' This function fits statistical models to multiply imputed datasets and pools the results using Rubin's Rules.
#' It supports various regression models, including negative binomial, logistic, Poisson, linear, and Cox regression.
#' Enhanced to handle interaction terms, restricted cubic splines (rcs), and polynomial terms.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for GLM models.
#' @param predictor_vars A vector of predictor variables to include in the model. Can include interaction terms with : notation.
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param imp_n Number of imputations (if NULL, it is detected automatically).
#' @param model_type The model type:
#'   \itemize{
#'     \item "nb" for Negative Binomial regression
#'     \item "lm" for Linear regression
#'     \item "bin" for Logistic regression (binomial family)
#'     \item "poisson" for Poisson regression
#'     \item "gamma" for Gamma regression
#'     \item "quasipoisson" for Overdispersed Poisson regression
#'     \item "quasibinomial" for Overdispersed logistic regression
#'     \item "cox" for Cox Proportional Hazards regression
#'   }
#' @param followup_offset Whether to include offset for follow-up duration ("Yes" or "No").
#' @param followup_col The column representing follow-up duration (optional, required if followup_offset = "Yes").
#' @param trial_factor Whether to include trial as a factor ("Yes" or "No").
#' @param trial_col The column for trial factor adjustment (optional).
#' @param time_col The time variable for Cox regression (only required if model_type = "cox").
#' @param event_col The event variable for Cox regression (only required if model_type = "cox").
#' @param formula_string Optional custom formula string. If provided, overrides the formula created from predictor_vars.
#' @param highlight_interactions Whether to flag interaction terms in the output (default: TRUE).
#' @param spline_terms A named list to specify variables to be modeled with restricted cubic splines. Each element should be a list with 'var' (variable name) and 'knots' (number of knots or explicit knot positions).
#' @param poly_terms A named list to specify variables to be modeled with polynomial terms. Each element should be a list with 'var' (variable name) and 'degree' (2 for quadratic, 3 for cubic).
#' @param include_spline_terms Whether to include the individual spline terms in the output (default: FALSE).
#' @param include_poly_terms Whether to include the individual polynomial terms in the output (default: TRUE).
#' @return A data frame with model estimates, confidence intervals, and p-values.
#'   If the model type supports exponentiation, exponentiated estimates (e.g., odds ratios, rate ratios, hazard ratios) are provided.
#' @import dplyr MASS mice survival splines
#' @export
MI_estimates <- function(data,
                         outcome_var,
                         predictor_vars,
                         imp_col = ".imp",
                         imp_n = NULL,
                         model_type = "nb",
                         followup_offset = "No",  # Default to "No"
                         followup_col = NULL,
                         trial_factor = "No",  # Default to "No"
                         trial_col = NULL,
                         time_col = NULL,
                         event_col = NULL,
                         formula_string = NULL,
                         highlight_interactions = TRUE,
                         spline_terms = NULL,
                         poly_terms = NULL,
                         include_spline_terms = FALSE,
                         include_poly_terms = TRUE) {

  # Load required packages
  require(dplyr)
  require(MASS)         # For Negative Binomial
  require(mice)         # For multiple imputations
  require(rlang)        # For formula manipulation
  require(survival)     # For Cox regression
  require(splines)      # For basic spline functions

  # Check if rms package is available if spline_terms is provided
  if (!is.null(spline_terms)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      message("Package 'rms' is not available. Using splines package instead.")
      use_rms <- FALSE
    } else {
      require(rms)
      use_rms <- TRUE
    }
  } else {
    use_rms <- FALSE
  }

  # Input validation
  if (!imp_col %in% names(data)) stop("imp_col not found in data.")
  if (!outcome_var %in% names(data) && model_type != "cox") stop("outcome_var not found in data.")

  # Validate spline_terms if provided
  if (!is.null(spline_terms)) {
    if (!is.list(spline_terms)) {
      stop("spline_terms must be a list")
    }

    for (i in seq_along(spline_terms)) {
      if (!is.list(spline_terms[[i]]) || is.null(spline_terms[[i]]$var)) {
        stop("Each element in spline_terms must be a list with at least a 'var' element")
      }

      if (!spline_terms[[i]]$var %in% names(data)) {
        stop(paste("Spline variable", spline_terms[[i]]$var, "not found in data"))
      }

      if (is.null(spline_terms[[i]]$knots)) {
        # Default to 3 knots if not specified
        spline_terms[[i]]$knots <- 3
      } else if (!is.numeric(spline_terms[[i]]$knots)) {
        stop("Knots must be numeric")
      }
    }
  }

  # Validate polynomial terms if provided
  if (!is.null(poly_terms)) {
    if (!is.list(poly_terms)) {
      stop("poly_terms must be a list")
    }

    for (i in seq_along(poly_terms)) {
      if (!is.list(poly_terms[[i]]) || is.null(poly_terms[[i]]$var)) {
        stop("Each element in poly_terms must be a list with at least a 'var' element")
      }

      if (!poly_terms[[i]]$var %in% names(data)) {
        stop(paste("Polynomial variable", poly_terms[[i]]$var, "not found in data"))
      }

      if (is.null(poly_terms[[i]]$degree)) {
        # Default to quadratic if not specified
        poly_terms[[i]]$degree <- 2
      } else if (!is.numeric(poly_terms[[i]]$degree) || !(poly_terms[[i]]$degree %in% c(2, 3))) {
        stop("Degree must be either 2 (quadratic) or 3 (cubic)")
      }
    }
  }

  # Check if formula_string is provided, otherwise validate predictor_vars
  if (is.null(formula_string)) {
    # Check for interaction terms in predictor_vars
    has_interactions <- any(grepl(":", predictor_vars))

    if (!has_interactions) {
      # Regular variable check if no interactions
      if (!all(predictor_vars %in% names(data))) {
        missing_vars <- predictor_vars[!predictor_vars %in% names(data)]
        stop(paste("Predictor variables not found in data:", paste(missing_vars, collapse = ", ")))
      }
    } else {
      # For interactions, we need to parse them more carefully
      # Extract all unique variable names from interaction terms
      all_vars <- character(0)
      for (term in predictor_vars) {
        if (grepl(":", term)) {
          # Split interaction terms
          interaction_vars <- unlist(strsplit(term, ":"))
          # Clean up any whitespace
          interaction_vars <- trimws(interaction_vars)
          all_vars <- c(all_vars, interaction_vars)
        } else {
          all_vars <- c(all_vars, term)
        }
      }
      all_vars <- unique(all_vars)

      # Check if all extracted variables exist
      if (!all(all_vars %in% names(data))) {
        missing_vars <- all_vars[!all_vars %in% names(data)]
        stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
      }
    }
  }

  if (!is.null(followup_col) && !followup_col %in% names(data)) stop("followup_col not found in data.")
  if (!is.null(trial_col) && !trial_col %in% names(data)) stop("trial_col not found in data.")
  if (!is.null(time_col) && !time_col %in% names(data)) stop("time_col not found in data.")
  if (!is.null(event_col) && !event_col %in% names(data)) stop("event_col not found in data.")

  # Validate followup_offset input
  if (!followup_offset %in% c("Yes", "No")) stop("followup_offset must be either 'Yes' or 'No'.")

  if (followup_offset == "Yes" && is.null(followup_col)) {
    stop("If followup_offset = 'Yes', followup_col must be provided.")
  }

  if (!is.null(followup_col) && any(data[[followup_col]] <= 0, na.rm = TRUE)) {
    stop("followup_col must be strictly positive for offset.")
  }

  # Validate trial_factor input
  if (!trial_factor %in% c("Yes", "No")) stop("trial_factor must be either 'Yes' or 'No'.")

  if (trial_factor == "Yes" && is.null(trial_col)) {
    stop("If trial_factor = 'Yes', trial_col must be provided.")
  }

  # Determine number of imputations if not provided
  if (is.null(imp_n)) {
    imp_n <- length(unique(data[[imp_col]]))
    if (imp_n < 2) stop("At least 2 imputations required.")
  }

  # Construct formula with optional terms
  trial_term <- if (trial_factor == "Yes")
    paste("+ as.factor(", trial_col, ")") else ""

  offset_term <- if (followup_offset == "Yes")
    paste("+ offset(log(", followup_col, "))") else ""

  # Process spline terms if provided
  spline_formula_parts <- character(0)
  if (!is.null(spline_terms)) {
    for (i in seq_along(spline_terms)) {
      var_name <- spline_terms[[i]]$var
      knots <- spline_terms[[i]]$knots

      # Remove the original variable from predictor_vars if it's there
      if (!is.null(predictor_vars)) {
        predictor_vars <- predictor_vars[predictor_vars != var_name]
      }

      # Create the spline term
      if (use_rms) {
        # Use rcs() from rms package
        if (length(knots) == 1) {
          # Number of knots specified
          spline_formula_parts <- c(spline_formula_parts,
                                    paste0("rcs(", var_name, ", ", knots, ")"))
        } else {
          # Explicit knot positions specified
          spline_formula_parts <- c(spline_formula_parts,
                                    paste0("rcs(", var_name, ", c(",
                                           paste(knots, collapse = ", "), "))"))
        }
      } else {
        # Use bs() from splines package as a fallback
        if (length(knots) == 1) {
          # For bs(), df = knots + degree - 1 (for cubic splines with degree=3)
          df <- knots + 2  # degree=3, so 3-1=2
          spline_formula_parts <- c(spline_formula_parts,
                                    paste0("bs(", var_name, ", df = ", df, ", degree = 3)"))
        } else {
          # Explicit knot positions
          spline_formula_parts <- c(spline_formula_parts,
                                    paste0("bs(", var_name, ", knots = c(",
                                           paste(knots, collapse = ", "), "), degree = 3)"))
        }
      }
    }
  }

  # Process polynomial terms if provided
  poly_formula_parts <- character(0)
  if (!is.null(poly_terms)) {
    for (i in seq_along(poly_terms)) {
      var_name <- poly_terms[[i]]$var
      degree <- poly_terms[[i]]$degree

      # Remove the original variable from predictor_vars if it's there
      if (!is.null(predictor_vars)) {
        predictor_vars <- predictor_vars[predictor_vars != var_name]
      }

      # Create the polynomial term using poly() function
      # The raw=TRUE option creates simple polynomials rather than orthogonal ones
      poly_formula_parts <- c(poly_formula_parts,
                              paste0("poly(", var_name, ", degree = ", degree, ", raw = TRUE)"))
    }
  }

  # Define model formula
  if (is.null(formula_string)) {
    # Combine regular predictors, spline terms, and polynomial terms
    all_terms <- c(predictor_vars, spline_formula_parts, poly_formula_parts)

    if (model_type == "cox") {
      if (is.null(time_col) || is.null(event_col)) {
        stop("For Cox regression, time_col and event_col must be provided.")
      }
      formula_string <- paste("Surv(", time_col, ",", event_col, ") ~",
                              paste(all_terms, collapse = " + "), trial_term)
    } else {
      formula_string <- paste(outcome_var, "~",
                              paste(all_terms, collapse = " + "), offset_term, trial_term)
    }
  } else {
    # If a custom formula is provided, warn about potential conflicts
    if (!is.null(spline_terms) || !is.null(poly_terms)) {
      warning("Using custom formula_string with spline_terms or poly_terms. Make sure your formula includes these terms correctly.")
    }

    # Use the provided formula_string but ensure it has correct components
    if (model_type == "cox") {
      if (is.null(time_col) || is.null(event_col)) {
        stop("For Cox regression, time_col and event_col must be provided.")
      }

      # Check if formula already has Surv() term
      if (!grepl("^Surv\\(", formula_string)) {
        formula_string <- paste("Surv(", time_col, ",", event_col, ") ~", formula_string)
      }
    } else {
      # Check if formula already has outcome variable
      if (!grepl(paste0("^", outcome_var, "\\s*~"), formula_string)) {
        formula_string <- paste(outcome_var, "~", formula_string)
      }
    }

    # Add trial term if needed and not already in formula
    if (trial_factor == "Yes" && !grepl(paste0("as\\.factor\\(", trial_col, "\\)"), formula_string)) {
      formula_string <- paste(formula_string, trial_term)
    }

    # Add offset if needed and not already in formula
    if (followup_offset == "Yes" && !grepl("offset\\(log\\(.*\\)\\)", formula_string)) {
      formula_string <- paste(formula_string, offset_term)
    }
  }

  # Detect interaction terms
  interaction_terms <- character(0)
  if (grepl(":", formula_string)) {
    # Extract terms after the "~"
    terms_part <- strsplit(formula_string, "~")[[1]][2]
    # Split by "+"
    terms <- strsplit(terms_part, "\\+")[[1]]
    # Trim whitespace
    terms <- trimws(terms)
    # Find terms containing ":"
    interaction_terms <- terms[grepl(":", terms)]
  }

  # Detect spline terms
  spline_terms_detected <- character(0)
  if (grepl("rcs\\(|bs\\(", formula_string)) {
    # Extract terms after the "~"
    terms_part <- strsplit(formula_string, "~")[[1]][2]
    # Split by "+"
    terms <- strsplit(terms_part, "\\+")[[1]]
    # Trim whitespace
    terms <- trimws(terms)
    # Find terms containing "rcs(" or "bs("
    spline_terms_detected <- terms[grepl("rcs\\(|bs\\(", terms)]
  }

  # Detect polynomial terms
  poly_terms_detected <- character(0)
  if (grepl("poly\\(", formula_string)) {
    # Extract terms after the "~"
    terms_part <- strsplit(formula_string, "~")[[1]][2]
    # Split by "+"
    terms <- strsplit(terms_part, "\\+")[[1]]
    # Trim whitespace
    terms <- trimws(terms)
    # Find terms containing "poly("
    poly_terms_detected <- terms[grepl("poly\\(", terms)]
  }

  model_formula <- as.formula(formula_string)

  # Fit models to each imputed dataset
  res_comb <- vector("list", imp_n)
  for (i in 1:imp_n) {
    data_subset <- dplyr::filter(data, !!rlang::sym(imp_col) == i)
    if (nrow(data_subset) == 0) stop(paste("No data for imputation", i))

    # Ensure required packages are loaded for the modeling
    if (length(spline_terms_detected) > 0) {
      if (any(grepl("rcs\\(", spline_terms_detected)) && use_rms) {
        # For rcs(), we need rms package
        suppressPackageStartupMessages(require(rms))
      }
      if (any(grepl("bs\\(", spline_terms_detected))) {
        # For bs(), we need splines package
        suppressPackageStartupMessages(require(splines))
      }
    }

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

  # Pool results using mice::pool()
  pooled <- mice::pool(res_comb)
  summary_pool <- summary(pooled, conf.int = TRUE, exp = FALSE)

  # Determine if exponentiation is needed
  exp_required <- model_type %in% c("nb", "bin", "poisson", "quasipoisson", "quasibinomial", "cox")

  # Construct final results table, dynamically removing terms based on trial_col and spline terms
  Results_multivariate_analysis <- summary_pool %>%
    mutate(exp_estimate = exp(estimate),
           exp_CI95_lower = exp(`2.5 %`),
           exp_CI95_upper = exp(`97.5 %`)) %>%
    select(term, estimate, `2.5 %`, `97.5 %`, exp_estimate, exp_CI95_lower, exp_CI95_upper, p.value)

  # Filter out trial terms if needed
  if (!is.null(trial_col)) {
    Results_multivariate_analysis <- Results_multivariate_analysis %>%
      filter(!grepl(trial_col, term))
  }

  # Filter out individual spline terms if requested
  if (!include_spline_terms && length(spline_terms_detected) > 0) {
    # Create patterns to match individual spline components
    spline_patterns <- sapply(spline_terms_detected, function(x) {
      if (grepl("rcs\\(", x)) {
        # Extract variable name from rcs() term
        var_name <- gsub("rcs\\(([^,]+),.*", "\\1", x)
      } else if (grepl("bs\\(", x)) {
        # Extract variable name from bs() term
        var_name <- gsub("bs\\(([^,]+),.*", "\\1", x)
      } else {
        return("")
      }
      var_name <- trimws(var_name)
      # Pattern to match the variable's spline components
      paste0("^", var_name, "'")
    })

    # Combine patterns
    pattern <- paste(spline_patterns, collapse = "|")

    # Keep the main spline terms but filter out the individual components
    if (pattern != "") {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        filter(!grepl(pattern, term))
    }
  }

  # Filter out individual polynomial terms if requested
  if (!include_poly_terms && length(poly_terms_detected) > 0) {
    # Create patterns to match individual polynomial components
    poly_patterns <- sapply(poly_terms_detected, function(x) {
      if (grepl("poly\\(", x)) {
        # Extract variable name from poly() term
        var_name <- gsub("poly\\(([^,]+),.*", "\\1", x)
        var_name <- trimws(var_name)
        # Pattern to match the variable's polynomial components
        return(paste0("^poly\\(", var_name, ".*\\)"))
      } else {
        return("")
      }
    })

    # Combine patterns
    pattern <- paste(poly_patterns, collapse = "|")

    # Remove polynomial terms according to pattern
    if (pattern != "") {
      Results_multivariate_analysis <- Results_multivariate_analysis %>%
        filter(!grepl(pattern, term))
    }
  }

  # Add indicator for interaction terms if requested
  if (highlight_interactions && length(interaction_terms) > 0) {
    Results_multivariate_analysis <- Results_multivariate_analysis %>%
      mutate(is_interaction = sapply(term, function(t) any(sapply(interaction_terms, function(i) grepl(i, t, fixed = TRUE)))))
  }

  # Add indicator for spline terms if requested
  if (highlight_interactions && length(spline_terms_detected) > 0) {
    Results_multivariate_analysis <- Results_multivariate_analysis %>%
      mutate(is_spline = sapply(term, function(t) {
        is_spline_term <- grepl("rcs\\(|bs\\(", t)
        if (is_spline_term) return(TRUE)

        # Check if it's a component of a spline
        for (spline_term in spline_terms_detected) {
          if (grepl("rcs\\(", spline_term)) {
            var_name <- gsub("rcs\\(([^,]+),.*", "\\1", spline_term)
          } else if (grepl("bs\\(", spline_term)) {
            var_name <- gsub("bs\\(([^,]+),.*", "\\1", spline_term)
          } else {
            next
          }
          var_name <- trimws(var_name)
          if (grepl(paste0("^", var_name, "'"), t)) return(TRUE)
        }
        return(FALSE)
      }))
  }

  # Add indicator for polynomial terms if requested
  if (highlight_interactions && length(poly_terms_detected) > 0) {
    Results_multivariate_analysis <- Results_multivariate_analysis %>%
      mutate(is_polynomial = sapply(term, function(t) {
        # Check for direct polynomial term
        is_poly_term <- grepl("poly\\(", t)
        if (is_poly_term) return(TRUE)

        # Check for individual polynomial components (like .1, .2, .3)
        for (poly_term in poly_terms_detected) {
          if (grepl("poly\\(", poly_term)) {
            var_name <- gsub("poly\\(([^,]+),.*", "\\1", poly_term)
            var_name <- trimws(var_name)
            degree <- as.numeric(gsub(".*degree = ([0-9]+).*", "\\1", poly_term))
            # Check for each possible polynomial component
            for (d in 1:degree) {
              if (grepl(paste0("poly\\(", var_name, ".*\\)", d), t)) return(TRUE)
            }
          }
        }
        return(FALSE)
      }))
  }

  # Add formula and term information as attributes
  attr(Results_multivariate_analysis, "formula") <- formula_string
  attr(Results_multivariate_analysis, "has_interactions") <- length(interaction_terms) > 0
  attr(Results_multivariate_analysis, "interaction_terms") <- if(length(interaction_terms) > 0) interaction_terms else NULL
  attr(Results_multivariate_analysis, "has_splines") <- length(spline_terms_detected) > 0
  attr(Results_multivariate_analysis, "spline_terms") <- if(length(spline_terms_detected) > 0) spline_terms_detected else NULL
  attr(Results_multivariate_analysis, "spline_method") <- if(use_rms) "rcs" else "bs"
  attr(Results_multivariate_analysis, "has_polynomials") <- length(poly_terms_detected) > 0
  attr(Results_multivariate_analysis, "polynomial_terms") <- if(length(poly_terms_detected) > 0) poly_terms_detected else NULL

  return(Results_multivariate_analysis)
}
