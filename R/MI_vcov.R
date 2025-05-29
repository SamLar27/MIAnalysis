#' Calculate Variance-Covariance Matrix from Multiple Imputation Results
#'
#' @description
#' This function calculates the pooled variance-covariance matrix from multiple
#' imputation results following Rubin's rules. It supports a wide range of model types
#' (linear, binomial, Poisson, negative binomial, etc.) and special terms including
#' restricted cubic splines, interactions, polynomial terms, and random effects.
#'
#' @param data A data frame containing the imputed dataset. Each imputation should be
#'   identified by the column specified in \code{imp_col}.
#' @param outcome_var Character. The name of the dependent variable for the model.
#' @param predictor_vars Character vector. Names of predictor variables. Can include
#'   specially formatted terms like \code{"rcs(x1, 4)"} for restricted cubic splines,
#'   \code{"x1:x2"} or \code{"x1*x2"} for interactions.
#' @param imp_col Character. The name of the column identifying imputation numbers
#'   (default is ".imp").
#' @param imp_n Numeric. Number of imputations to use. If NULL (default), all
#'   unique values in \code{imp_col} are used.
#' @param model_type Character. The type of model to fit. Options are:
#'   \itemize{
#'     \item \code{"lm"}: Linear model
#'     \item \code{"bin"}: Binary/Logistic model
#'     \item \code{"nb"}: Negative binomial model (default)
#'     \item \code{"poisson"}: Poisson model
#'     \item \code{"gamma"}: Gamma model
#'     \item \code{"quasipoisson"}: Quasi-Poisson model
#'     \item \code{"quasibinomial"}: Quasi-binomial model
#'     \item \code{"cox"}: Cox proportional hazards model
#'   }
#' @param followup_offset Character. Whether to include offset for follow-up duration
#'   (\code{"Yes"} or \code{"No"}). Default is \code{"No"}.
#' @param followup_col Character. The name of the column representing follow-up duration
#'   (required if \code{followup_offset = "Yes"}).
#' @param trial_factor Character. Whether to include a trial factor (\code{"Yes"} or \code{"No"}).
#'   Default is \code{"No"}.
#' @param trial_col Character. The name of the column identifying trials
#'   (required if \code{trial_factor = "Yes"}).
#' @param time_col Character. The name of the time variable for Cox regression
#'   (only required if \code{model_type = "cox"}).
#' @param event_col Character. The name of the event variable for Cox regression
#'   (only required if \code{model_type = "cox"}).
#' @param formula_string Character. A custom formula string. If provided, it overrides
#'   the formula constructed from \code{predictor_vars} and other arguments.
#' @param spline_terms Named list. Specifications for restricted cubic splines
#'   not already included in \code{predictor_vars}. Names are variable names,
#'   values are vectors of knot positions or a single number indicating the
#'   number of knots.
#' @param poly_terms Named list. Specifications for polynomial terms.
#'   Names are variable names, values are polynomial degrees.
#' @param interaction_terms List. Specifications for interaction terms
#'   not already included in \code{predictor_vars}. Multiple formats supported:
#'   \itemize{
#'     \item Simple pairs: \code{c("x1", "sex")}
#'     \item Pre-formatted strings: \code{"x1:sex"}
#'     \item Specific variable mappings: \code{list(sex = c("x1", "x2"))}
#'     \item All pairwise interactions: Vector of variables
#'   }
#' @param random_intercept Character. Whether to include a random intercept in the model
#'   (\code{"Yes"} or \code{"No"}). Default is \code{"No"}.
#' @param random_intercept_var Character. The name of the grouping variable for the
#'   random intercept (required if \code{random_intercept = "Yes"}).
#' @param return_full Logical. Whether to return additional attributes with the
#'   variance-covariance matrix (default is FALSE).
#'
#' @return A variance-covariance matrix for the coefficients, pooled according to
#'   Rubin's rules. If \code{return_full = TRUE}, additional attributes are included:
#'   \itemize{
#'     \item \code{mean_coefficients}: The pooled coefficient estimates
#'     \item \code{W_matrices}: List of within-imputation variance-covariance matrices
#'     \item \code{W_pooled}: Pooled within-imputation variance-covariance matrix
#'     \item \code{B_matrix}: Between-imputation variance-covariance matrix
#'     \item \code{formula}: The formula used for modeling
#'     \item \code{model_type}: The model type used
#'   }
#'
#' @details
#' This function implements Rubin's rules for pooling variance-covariance matrices
#' from multiple imputation. It handles a wide range of model types and special terms.
#'
#' For restricted cubic splines, you can either:
#' 1. Include them directly in \code{predictor_vars}, e.g., \code{"rcs(x1, 4)"}
#' 2. Specify them via the \code{spline_terms} parameter
#'
#' For interactions, you can either:
#' 1. Include them directly in \code{predictor_vars}, e.g., \code{"x1:sex"} or \code{"x1*sex"}
#' 2. Specify them via the \code{interaction_terms} parameter
#'
#' The function tries to handle model fitting errors gracefully, omitting failed
#' imputations and providing warnings.
#'
#' @seealso \code{\link{MI_estimates}} for obtaining the pooled coefficient estimates.
#' @importFrom mice pool
#' @importFrom stats coef vcov as.formula
#' @importFrom utils installed.packages
#'
#' @export


MI_vcov <- function(data,
                    outcome_var,
                    predictor_vars,
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
                    spline_terms = NULL,
                    poly_terms = NULL,
                    interaction_terms = NULL,
                    random_intercept = "No",
                    random_intercept_var = NULL,
                    return_full = FALSE) {
  requireNamespace("mice", quietly = TRUE)
  requireNamespace("stats", quietly = TRUE)

  # Validate imputations
  actual_imps <- sort(unique(data[[imp_col]]))
  if (is.null(imp_n)) {
    imp_n <- length(actual_imps)
  }

  # Check for categorical variables with no variation in some imputations
  cat_vars <- predictor_vars[sapply(data[predictor_vars], function(x) is.factor(x) || is.character(x))]
  if (length(cat_vars) > 0) {
    problem_vars <- c()
    for (var in cat_vars) {
      # Check each imputation
      for (imp in actual_imps) {
        imp_data <- subset(data, data[[imp_col]] == imp)
        unique_values <- unique(imp_data[[var]])

        # If only one unique value exists in this imputation
        if (length(unique_values) == 1) {
          problem_vars <- c(problem_vars, var)
          warning(paste0("Variable '", var, "' has only one unique value in imputation ",
                         imp, ": '", unique_values, "'. This may cause issues with model matrix dimensions."))
          break  # No need to check other imputations for this variable
        }
      }
    }

    if (length(problem_vars) > 0) {
      warning(paste0("The following categorical variables have no variation in some imputations: ",
                     paste(unique(problem_vars), collapse=", "),
                     ". Consider removing them or combining levels to avoid dimension errors."))
    }
  }

  # Prepare formula if not provided
  if (is.null(formula_string)) {
    # Process predictor variables for interaction terms and spline terms
    standard_vars <- character(0)
    interaction_vars <- character(0)
    spline_vars <- character(0)

    for (term in predictor_vars) {
      if (grepl("^rcs\\(", term)) {
        # This is a restricted cubic spline term
        spline_vars <- c(spline_vars, term)
      } else if (grepl("[\\*\\:]", term)) {
        # This is an interaction term
        interaction_vars <- c(interaction_vars, term)
      } else {
        # This is a standard variable
        standard_vars <- c(standard_vars, term)
      }
    }

    # Start with basic predictor terms (non-interaction)
    if (length(standard_vars) > 0) {
      formula_terms <- paste(standard_vars, collapse = " + ")
    } else {
      formula_terms <- "1"  # Intercept only if no standard variables
    }

    # Add interaction terms from predictor_vars
    if (length(interaction_vars) > 0) {
      formula_terms <- paste(formula_terms, "+", paste(interaction_vars, collapse = " + "))
    }

    # Add spline terms from predictor_vars
    if (length(spline_vars) > 0) {
      formula_terms <- paste(formula_terms, "+", paste(spline_vars, collapse = " + "))
    }

    # Add trial factor if requested
    if (trial_factor == "Yes") {
      formula_terms <- paste0(formula_terms, " + as.factor(", trial_col, ")")
    }

    # Add follow-up offset if requested
    if (followup_offset == "Yes") {
      formula_terms <- paste0(formula_terms, " + offset(log(", followup_col, "))")
    }

    # Add spline terms if provided
    if (!is.null(spline_terms) && length(spline_terms) > 0) {
      for (term in names(spline_terms)) {
        # Get knots for this variable
        knots <- spline_terms[[term]]
        # Add rcs term to formula
        spline_term <- paste0("rcs(", term, ", ",
                              paste(knots, collapse = ", "), ")")
        formula_terms <- paste0(formula_terms, " + ", spline_term)
      }
    }

    # Add polynomial terms if provided
    if (!is.null(poly_terms) && length(poly_terms) > 0) {
      for (term in names(poly_terms)) {
        # Get degree for this variable
        degree <- poly_terms[[term]]
        # Add polynomial term to formula
        poly_term <- paste0("poly(", term, ", degree = ", degree, ", raw = TRUE)")
        formula_terms <- paste0(formula_terms, " + ", poly_term)
      }
    }

    # Add interaction terms if provided
    if (!is.null(interaction_terms) && length(interaction_terms) > 0) {
      for (interaction in interaction_terms) {
        # Check if interaction is a list with specific variables
        if (is.list(interaction)) {
          # For detailed interactions like list(var1 = c("x1", "x2"), var2 = "group")
          # This allows specific variables to interact with specific other variables
          for (var1 in names(interaction)) {
            for (var2 in interaction[[var1]]) {
              interaction_term <- paste0(var1, ":", var2)
              formula_terms <- paste0(formula_terms, " + ", interaction_term)
            }
          }
        } else if (length(interaction) == 2) {
          # For simple pairs like c("x1", "group")
          interaction_term <- paste0(interaction[1], ":", interaction[2])
          formula_terms <- paste0(formula_terms, " + ", interaction_term)
        } else if (is.character(interaction) && grepl(":", interaction)) {
          # For pre-formatted interactions like "x1:group"
          formula_terms <- paste0(formula_terms, " + ", interaction)
        } else if (is.character(interaction) && length(interaction) >= 2) {
          # For a vector of variables to create all possible two-way interactions
          for (i in 1:(length(interaction)-1)) {
            for (j in (i+1):length(interaction)) {
              interaction_term <- paste0(interaction[i], ":", interaction[j])
              formula_terms <- paste0(formula_terms, " + ", interaction_term)
            }
          }
        }
      }
    }

    # Add random intercept if requested
    if (random_intercept == "Yes" && !is.null(random_intercept_var)) {
      formula_terms <- paste0(formula_terms, " + (1|", random_intercept_var, ")")
    }

    # Construct the full formula
    formula_string <- paste(outcome_var, "~", formula_terms)
  }

  # Convert to formula object
  model_formula <- as.formula(formula_string)

  # Function to fit model per imputation
  fit_model <- function(data_sub) {
    # Check if we need to handle spline terms
    has_spline_terms <- any(grepl("^rcs\\(", predictor_vars))

    if (has_spline_terms && !requireNamespace("rms", quietly = TRUE)) {
      stop("Package 'rms' is needed for restricted cubic splines. Please install it.")
    }

    # If using spline terms, we need to use rms functions
    if (has_spline_terms) {
      # Setup rms design
      rms::dd <- rms::datadist(data_sub)
      options(datadist = "dd")

      # Use appropriate rms modeling function
      switch(model_type,
             "lm" = rms::ols(model_formula, data = data_sub),
             "bin" = rms::lrm(model_formula, data = data_sub),
             "nb" = {
               # rms doesn't have native nb model, use MASS with mgcv::gam wrapper for splines
               if (!requireNamespace("mgcv", quietly = TRUE)) {
                 stop("Package 'mgcv' is needed for negative binomial models with splines. Please install it.")
               }
               mgcv::gam(model_formula, family = MASS::negative.binomial(theta = 1), data = data_sub)
             },
             "poisson" = rms::poisson(model_formula, data = data_sub),
             "cox" = rms::cph(model_formula, data = data_sub),
             stop("Unsupported model type for splines."))
    } else if (random_intercept == "Yes") {
      # Mixed models (same as before)
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' is needed for mixed models. Please install it.")
      }

      # Fit appropriate mixed model
      switch(model_type,
             "lm" = lme4::lmer(model_formula, data = data_sub),
             "poisson" = lme4::glmer(model_formula, family = poisson(), data = data_sub),
             "bin" = lme4::glmer(model_formula, family = binomial(), data = data_sub),
             "nb" = {
               if (!requireNamespace("glmmTMB", quietly = TRUE)) {
                 stop("Package 'glmmTMB' is needed for negative binomial mixed models. Please install it.")
               }
               glmmTMB::glmmTMB(model_formula, family = nbinom2, data = data_sub)
             },
             stop("Unsupported model type for mixed models."))
    } else {
      # Regular GLM (same as before)
      switch(model_type,
             "nb" = MASS::glm.nb(model_formula, data = data_sub),
             "lm" = stats::glm(model_formula, family = gaussian(), data = data_sub),
             "bin" = stats::glm(model_formula, family = binomial(), data = data_sub),
             "poisson" = stats::glm(model_formula, family = poisson(), data = data_sub),
             "gamma" = stats::glm(model_formula, family = Gamma(), data = data_sub),
             "quasipoisson" = stats::glm(model_formula, family = quasipoisson(), data = data_sub),
             "quasibinomial" = stats::glm(model_formula, family = quasibinomial(), data = data_sub),
             stop("Unsupported model type."))
    }
  }

  # Loop over imputations
  models_list <- lapply(actual_imps, function(i) {
    tryCatch({
      data_i <- subset(data, data[[imp_col]] == i)
      fit_model(data_i)
    }, error = function(e) {
      warning(paste("Error in imputation", i, ":", e$message))
      return(NULL)
    })
  })

  # Remove any NULL models (failed fits)
  models_list <- models_list[!sapply(models_list, is.null)]

  if (length(models_list) == 0) {
    stop("All models failed to fit. Check your model specification.")
  }

  # Extract variance-covariance matrices and coefficients with appropriate handling for different model types
  W_matrices <- lapply(models_list, function(model) {
    if (inherits(model, "glmmTMB") || inherits(model, "merMod")) {
      # For mixed models
      if (inherits(model, "glmmTMB")) {
        vcov(model)$cond  # For glmmTMB, extract conditional component
      } else {
        vcov(model)  # For lme4 models
      }
    } else {
      # For regular GLMs
      stats::vcov(model)
    }
  })

  coef_list <- lapply(models_list, function(model) {
    if (inherits(model, "glmmTMB") || inherits(model, "merMod")) {
      # For mixed models
      if (inherits(model, "glmmTMB")) {
        fixef(model)$cond  # For glmmTMB, extract conditional component
      } else {
        fixef(model)  # For lme4 models
      }
    } else {
      # For regular GLMs
      stats::coef(model)
    }
  })

  # Check if all models have the same coefficient names
  coef_names_list <- lapply(coef_list, names)
  all_same_names <- all(sapply(coef_names_list, function(x) identical(x, coef_names_list[[1]])))

  if (!all_same_names) {
    warning("Not all models have the same coefficient names. This may cause issues.")

    # Get union of all coefficient names
    all_coef_names <- unique(unlist(coef_names_list))

    # Create new coefficient list with all names
    coef_list_complete <- lapply(coef_list, function(coefs) {
      result <- rep(NA, length(all_coef_names))
      names(result) <- all_coef_names
      result[names(coefs)] <- coefs
      return(result)
    })

    # Update coefficient list
    coef_list <- coef_list_complete

    # Similarly, expand W_matrices to include all coefficients
    W_matrices_complete <- lapply(seq_along(W_matrices), function(i) {
      W <- W_matrices[[i]]
      coef_names <- coef_names_list[[i]]

      # Create empty matrix with all coefficient names
      result <- matrix(0, length(all_coef_names), length(all_coef_names))
      rownames(result) <- colnames(result) <- all_coef_names

      # Fill in the values we have
      for (r in coef_names) {
        for (c in coef_names) {
          result[r, c] <- W[r, c]
        }
      }

      return(result)
    })

    # Update W_matrices
    W_matrices <- W_matrices_complete
  }

  # Convert to matrices
  W_array <- simplify2array(W_matrices)
  coef_mat <- do.call(cbind, coef_list)

  # Handle NA values in coefficient matrix
  if (any(is.na(coef_mat))) {
    warning("NA values found in coefficients. These will be omitted when calculating means.")
  }

  # Pooled mean of coefficients
  mean_coef <- rowMeans(coef_mat, na.rm = TRUE)

  # Within-imputation variance (mean of covariances)
  W_pooled <- apply(W_array, 1:2, mean, na.rm = TRUE)

  # Between-imputation variance
  # Handle NA values in coefficient matrix for centering
  coef_centered <- matrix(NA, nrow = nrow(coef_mat), ncol = ncol(coef_mat))
  for (i in 1:nrow(coef_mat)) {
    for (j in 1:ncol(coef_mat)) {
      if (!is.na(coef_mat[i, j])) {
        coef_centered[i, j] <- coef_mat[i, j] - mean_coef[i]
      }
    }
  }

  # Calculate B matrix using available data
  B_matrix <- matrix(0, nrow = length(mean_coef), ncol = length(mean_coef))
  rownames(B_matrix) <- colnames(B_matrix) <- names(mean_coef)

  # Count valid observations for each coefficient pair
  valid_count <- matrix(0, nrow = length(mean_coef), ncol = length(mean_coef))

  # Calculate sum of crossproducts
  for (i in 1:ncol(coef_centered)) {
    # Skip imputations with NA values
    if (any(is.na(coef_centered[, i]))) {
      next
    }

    # Outer product of coefficient differences
    outer_prod <- outer(coef_centered[, i], coef_centered[, i])
    B_matrix <- B_matrix + outer_prod
    valid_count <- valid_count + 1
  }

  # Divide by number of valid observations - 1
  for (i in 1:nrow(B_matrix)) {
    for (j in 1:ncol(B_matrix)) {
      if (valid_count[i, j] > 1) {
        B_matrix[i, j] <- B_matrix[i, j] / (valid_count[i, j] - 1)
      } else {
        # Not enough data for this pair
        B_matrix[i, j] <- 0
        warning(paste("Not enough valid imputations for coefficient pair:",
                      rownames(B_matrix)[i], "and", colnames(B_matrix)[j]))
      }
    }
  }

  # Total pooled variance
  Total_vcov <- W_pooled + (1 + 1/imp_n) * B_matrix

  # Output
  if (return_full) {
    attr(Total_vcov, "mean_coefficients") <- mean_coef
    attr(Total_vcov, "W_matrices") <- W_matrices
    attr(Total_vcov, "W_pooled") <- W_pooled
    attr(Total_vcov, "B_matrix") <- B_matrix
    attr(Total_vcov, "formula") <- formula_string
    attr(Total_vcov, "model_type") <- model_type
  }

  return(Total_vcov)
}
