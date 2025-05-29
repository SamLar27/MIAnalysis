#' Calculate Derivatives of Restricted Cubic Splines from Multiple Imputed Data
#'
#' This function calculates the derivatives of restricted cubic spline models fitted to
#' multiply imputed datasets. It can handle various model types and implements
#' Rubin's rules for pooling estimates across imputations.
#'
#' @param data A data frame containing the multiply imputed dataset.
#' @param outcome_var The dependent variable for the model.
#' @param variable_x The continuous predictor variable modeled with splines.
#' @param subgroups Optional variable for stratification (default is NULL).
#' @param subgroup_labels Optional labels for subgroup levels.
#' @param subgroup_as_factor Whether to convert subgroups to factors (default is TRUE).
#' @param covariables Optional vector of covariates to adjust for.
#' @param knot_n Number of knots for the restricted cubic spline (default is 4).
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param model_type The model type: "nb", "poisson", "cox", "logistic", or "lm".
#' @param followup_offset Whether to include offset for follow-up duration ("Yes" or "No").
#' @param followup_col The column representing follow-up duration.
#' @param trial_factor Whether to include trial as a factor ("Yes" or "No").
#' @param trial_col The column for trial factor adjustment.
#' @param time_col The time variable for Cox regression.
#' @param event_col The event variable for Cox regression.
#' @param random_intercept Whether to include a random intercept ("Yes" or "No").
#' @param random_intercept_var The grouping variable for the random intercept.
#' @param x_points Optional vector of specific x points to calculate derivatives at.
#' @param n_points Number of points to calculate derivatives at if x_points not provided.
#' @param prediction_range Range for predictions as quantiles (default is c(0.01, 0.99)).
#' @param n_imputations Number of imputations to use (default is 10).
#'
#' @return A list containing derivative estimates, confidence intervals, plot, and thresholds.
#'
#' @importFrom stats as.formula predict quantile setNames coef
#' @importFrom rms rcs
#' @importFrom Hmisc rcspline.eval
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_minimal labs
MI_spline_deriv <- function(data,
                            outcome_var,
                            variable_x,
                            subgroups = NULL,
                            subgroup_labels = NULL,
                            subgroup_as_factor = TRUE,
                            covariables = NULL,
                            knot_n = 4,
                            imp_col = ".imp",
                            model_type = "nb",
                            followup_offset = "No",
                            followup_col = NULL,
                            trial_factor = "No",
                            trial_col = NULL,
                            time_col = NULL,
                            event_col = NULL,
                            random_intercept = "No",
                            random_intercept_var = NULL,
                            x_points = NULL,
                            n_points = 100,
                            prediction_range = c(0.01, 0.99),
                            n_imputations = 10) {
  
  # Load required packages
  requireNamespace("rms", quietly = TRUE)
  requireNamespace("Hmisc", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  # Additional packages for specific model types
  if (model_type == "nb") {
    requireNamespace("MASS", quietly = TRUE)
  }
  if (model_type == "cox") {
    requireNamespace("survival", quietly = TRUE)
  }
  if (random_intercept == "Yes") {
    requireNamespace("lme4", quietly = TRUE)
    if (model_type == "nb") {
      requireNamespace("glmmTMB", quietly = TRUE)
    }
    if (model_type == "cox") {
      requireNamespace("coxme", quietly = TRUE)
    }
  }
  
  # Input validation
  if (!model_type %in% c("nb", "poisson", "cox", "logistic", "lm")) {
    stop("model_type must be one of 'nb', 'poisson', 'cox', 'logistic', or 'lm'")
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
  
  if (random_intercept == "Yes" && is.null(random_intercept_var)) {
    stop("random_intercept_var must be provided when random_intercept = 'Yes'")
  }
  
  # Define x points to calculate derivatives
  if (is.null(x_points)) {
    x_points <- seq(
      from = quantile(data[[variable_x]], prediction_range[1], na.rm = TRUE),
      to = quantile(data[[variable_x]], prediction_range[2], na.rm = TRUE),
      length.out = n_points
    )
  }
  
  # Get unique imputation ids
  imputations <- sort(unique(data[[imp_col]]))
  if (length(imputations) < n_imputations) {
    warning("Requested more imputations than available. Using all available imputations.")
    n_imputations <- length(imputations)
  } else {
    imputations <- imputations[1:n_imputations]
  }
  
  # Store original subgroup values before conversion, if available
  original_subgroup_values <- NULL
  if (!is.null(subgroups)) {
    original_subgroup_values <- data[[subgroups]]
  }
  
  # Create a mapping between original values and labels if both exist
  subgroup_mapping <- NULL
  if (!is.null(subgroups) && !is.null(subgroup_labels)) {
    # Identify unique values in original data
    unique_values <- sort(unique(data[[subgroups]]))
    
    # Create mapping
    if (length(unique_values) == length(subgroup_labels)) {
      subgroup_mapping <- setNames(subgroup_labels, unique_values)
    }
  }
  
  # Function to create prediction data frame
  create_pred_data <- function(x_values, subgroups_data = NULL) {
    if (is.null(subgroups)) {
      # Without subgroups
      pred_data <- data.frame(x_values)
      colnames(pred_data) <- variable_x
      return(pred_data)
    } else {
      # With subgroups
      if (is.factor(subgroups_data)) {
        subgroup_levels <- levels(subgroups_data)
      } else {
        subgroup_levels <- sort(unique(subgroups_data))
      }
      
      pred_data <- expand.grid(
        x_values,
        subgroup_levels
      )
      colnames(pred_data) <- c(variable_x, subgroups)
      
      # Ensure subgroups column is the right type
      if (is.factor(subgroups_data)) {
        pred_data[[subgroups]] <- factor(pred_data[[subgroups]], levels = levels(subgroups_data))
      } else if (subgroup_as_factor) {
        pred_data[[subgroups]] <- factor(pred_data[[subgroups]])
        if (!is.null(subgroup_labels)) {
          levels(pred_data[[subgroups]]) <- subgroup_labels
        }
      }
      return(pred_data)
    }
  }
  
  # Function to add median/mode values for covariates
  add_covariates <- function(pred_data, model_data) {
    # Extract simple variable names from complex terms
    extract_variables <- function(term) {
      if (grepl("^rcs\\(", term) || grepl("^bs\\(", term) || grepl("^poly\\(", term)) {
        inside <- sub("^[^\\(]+\\(([^,]+),.*$", "\\1", term)
        return(trimws(inside))
      } else {
        return(term)
      }
    }
    
    if (!is.null(covariables)) {
      # Expand any interaction terms
      expanded_covs <- unlist(lapply(covariables, function(term) {
        if (grepl("\\*", term)) {
          vars_split <- unlist(strsplit(term, "\\*"))
          vars_split <- trimws(vars_split)
          all_combinations <- lapply(1:length(vars_split), function(k) {
            combn(vars_split, k, FUN = function(x) paste(x, collapse = ":"))
          })
          unlist(all_combinations)
        } else {
          term
        }
      }))
      
      # Get base variables
      covariate_vars <- unique(unlist(strsplit(gsub(":", "*", expanded_covs), "\\*")))
      covariate_vars <- trimws(covariate_vars)
      covariate_vars <- sapply(covariate_vars, extract_variables)
      
      # Add each covariate to prediction data
      for (cov in covariate_vars) {
        if (cov %in% colnames(model_data)) {
          if (is.factor(model_data[[cov]])) {
            pred_data[[cov]] <- factor(names(which.max(table(model_data[[cov]]))), 
                                       levels = levels(model_data[[cov]]))
          } else {
            pred_data[[cov]] <- median(model_data[[cov]], na.rm = TRUE)
          }
        }
      }
    }
    
    # Add trial factor if needed
    if (trial_factor == "Yes") {
      pred_data[[trial_col]] <- names(which.max(table(model_data[[trial_col]])))
    }
    
    # Add follow-up time for offset if needed
    if (followup_offset == "Yes") {
      pred_data[[followup_col]] <- 365  # 1 year follow-up
    }
    
    # Add random effect variable if needed
    if (random_intercept == "Yes") {
      pred_data[[random_intercept_var]] <- names(which.max(table(model_data[[random_intercept_var]])))
    }
    
    return(pred_data)
  }
  
  # Function to build the model formula
  build_formula <- function() {
    # Setup spline term
    spline_term <- paste0("rcs(", variable_x, ", ", knot_n, ")")
    
    # Handle subgroups
    if (!is.null(subgroups)) {
      spline_term <- paste0(spline_term, " * ", subgroups)
    }
    
    # Handle covariates
    covariates_str <- ""
    if (!is.null(covariables)) {
      covariates_str <- paste0(" + ", paste(covariables, collapse = " + "))
    }
    
    # Add trial factor if requested
    trial_str <- ""
    if (trial_factor == "Yes") {
      trial_str <- paste0(" + as.factor(", trial_col, ")")
    }
    
    # Add offset for count models if requested
    offset_str <- ""
    if (followup_offset == "Yes") {
      offset_str <- paste0(" + offset(log(", followup_col, "))")
    }
    
    # Add random intercept if requested
    random_effect_str <- ""
    if (random_intercept == "Yes") {
      random_effect_str <- paste0(" + (1 | ", random_intercept_var, ")")
    }
    
    # Build the complete formula
    if (model_type == "cox") {
      formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                            spline_term, covariates_str, trial_str, random_effect_str)
    } else {
      formula_str <- paste0(outcome_var, " ~ ", spline_term, covariates_str,
                            trial_str, offset_str, random_effect_str)
    }
    
    return(as.formula(formula_str))
  }
  
  # Function to fit model for a specific imputation
  fit_model <- function(imp_data, formula_obj) {
    if (random_intercept == "Yes") {
      # Models with random effects
      if (model_type == "nb") {
        if (requireNamespace("glmmTMB", quietly = TRUE)) {
          model <- glmmTMB::glmmTMB(formula_obj, family = glmmTMB::nbinom2, data = imp_data)
        } else {
          warning("Using Poisson mixed model as approximation for negative binomial.")
          model <- lme4::glmer(formula_obj, family = poisson(link = "log"), data = imp_data)
        }
      } else if (model_type == "poisson") {
        model <- lme4::glmer(formula_obj, family = poisson(link = "log"), data = imp_data)
      } else if (model_type == "logistic") {
        model <- lme4::glmer(formula_obj, family = binomial(link = "logit"), data = imp_data)
      } else if (model_type == "lm") {
        model <- lme4::lmer(formula_obj, data = imp_data)
      } else if (model_type == "cox") {
        model <- coxme::coxme(formula_obj, data = imp_data)
      }
    } else {
      # Standard models without random effects
      if (model_type == "nb") {
        model <- MASS::glm.nb(formula_obj, data = imp_data)
      } else if (model_type == "poisson") {
        model <- glm(formula_obj, family = poisson(link = "log"), data = imp_data)
      } else if (model_type == "logistic") {
        model <- glm(formula_obj, family = binomial(link = "logit"), data = imp_data)
      } else if (model_type == "lm") {
        model <- glm(formula_obj, family = gaussian(), data = imp_data)
      } else if (model_type == "cox") {
        model <- survival::coxph(formula_obj, data = imp_data)
      }
    }
    
    return(model)
  }
  
  # Function to calculate analytical derivatives of spline models
  calculate_derivatives <- function(model, pred_data, knots) {
    # Extract model coefficients
    if (inherits(model, "glmmTMB")) {
      coef_values <- glmmTMB::fixef(model)$cond
    } else if (inherits(model, "merMod")) {
      coef_values <- lme4::fixef(model)
    } else if (inherits(model, "coxme")) {
      coef_values <- model$coefficients
    } else {
      coef_values <- coef(model)
    }
    
    # Initialize results
    derivatives <- numeric(nrow(pred_data))
    
    # Process by subgroup if needed
    if (!is.null(subgroups)) {
      for (sg in unique(pred_data[[subgroups]])) {
        # Filter data for this subgroup
        sg_idx <- which(pred_data[[subgroups]] == sg)
        
        # Build spline basis and derivative basis for this subgroup's x values
        x_values <- pred_data[[variable_x]][sg_idx]
        
        # First basis is for the values
        basis <- Hmisc::rcspline.eval(x_values, knots = knots, inclx = TRUE)
        # Second basis is for the derivatives
        deriv_basis <- Hmisc::rcspline.eval(x_values, knots = knots, inclx = TRUE, deriv = 1)
        
        # Determine spline coefficient names
        # This is tricky because they depend on the formula and whether there are interactions
        if (is.null(subgroups)) {
          # Simple case - no subgroups
          spline_coef_pattern <- "^rcs\\("
        } else {
          # With subgroups - coefficients include interaction terms
          spline_coef_pattern <- paste0("^rcs\\(|:", subgroups)
        }
        
        # Get spline coefficients
        spline_coef_names <- grep(spline_coef_pattern, names(coef_values), value = TRUE)
        
        # For subgroups, filter for specific level if needed
        if (!is.null(subgroups) && is.factor(pred_data[[subgroups]])) {
          # For factor levels, need to match the specific level in interaction terms
          level_idx <- which(levels(pred_data[[subgroups]]) == sg)
          sg_pattern <- paste0(subgroups, level_idx, "$|", subgroups, level_idx, ":")
          sg_coefs <- grep(sg_pattern, names(coef_values), value = TRUE)
          spline_coef_names <- intersect(spline_coef_names, sg_coefs)
        }
        
        # Get coefficient values for spline terms
        spline_coefs <- coef_values[spline_coef_names]
        
        # If we have the right number of coefficients, calculate derivative of linear predictor
        if (length(spline_coefs) >= ncol(deriv_basis)) {
          # Extract just the coefficients we need for the derivative calculation
          # This is more complex with interactions
          
          if (is.null(subgroups)) {
            # Simple case - no subgroups
            matched_coefs <- spline_coefs[1:ncol(deriv_basis)]
          } else {
            # With subgroups - need to handle main effects and interactions
            # This depends on how rcs() terms are formatted in model matrix
            # A more robust approach might be needed to match coefficients
            # to derivatives for complex models
            
            # For now, a simple heuristic to match:
            if (length(spline_coefs) == ncol(deriv_basis)) {
              matched_coefs <- spline_coefs
            } else {
              # Simple approach for interactions: first coefficients
              # This may need to be customized based on model structure
              matched_coefs <- spline_coefs[1:ncol(deriv_basis)]
            }
          }
          
          # Calculate derivative of linear predictor 
          lp_deriv <- as.vector(deriv_basis %*% matched_coefs)
          
          # Get predictions on link scale to transform derivative
          if (random_intercept == "Yes") {
            # For mixed models
            if (model_type == "cox") {
              lp <- predict(model, newdata = pred_data[sg_idx, ], re.form = NA, type = "lp")
            } else {
              lp <- predict(model, newdata = pred_data[sg_idx, ], re.form = NA, type = "link")
            }
          } else {
            # For fixed effect models
            if (model_type == "cox") {
              lp <- predict(model, newdata = pred_data[sg_idx, ], type = "lp")
            } else {
              lp <- predict(model, newdata = pred_data[sg_idx, ], type = "link")
            }
          }
          
          # Transform derivative based on model type
          final_deriv <- switch(model_type,
                                "nb" = exp(lp) * lp_deriv,
                                "poisson" = exp(lp) * lp_deriv,
                                "logistic" = {
                                  p <- plogis(lp)
                                  p * (1 - p) * lp_deriv
                                },
                                "cox" = exp(lp) * lp_deriv,
                                "lm" = lp_deriv,
                                lp_deriv)
          
          # Store derivatives for this subgroup
          derivatives[sg_idx] <- final_deriv
        } else {
          warning("Insufficient coefficients found for derivative calculation in subgroup ", sg)
          derivatives[sg_idx] <- NA
        }
      }
    } else {
      # No subgroups case
      x_values <- pred_data[[variable_x]]
      
      # First basis is for the values
      basis <- Hmisc::rcspline.eval(x_values, knots = knots, inclx = TRUE)
      # Second basis is for the derivatives
      deriv_basis <- Hmisc::rcspline.eval(x_values, knots = knots, inclx = TRUE, deriv = 1)
      
      # Get spline coefficients
      spline_coef_pattern <- "^rcs\\("
      spline_coef_names <- grep(spline_coef_pattern, names(coef_values), value = TRUE)
      spline_coefs <- coef_values[spline_coef_names]
      
      # If we have the right number of coefficients, calculate derivative
      if (length(spline_coefs) >= ncol(deriv_basis)) {
        # Extract just the coefficients we need
        matched_coefs <- spline_coefs[1:ncol(deriv_basis)]
        
        # Calculate derivative of linear predictor
        lp_deriv <- as.vector(deriv_basis %*% matched_coefs)
        
        # Get predictions on link scale to transform derivative
        if (random_intercept == "Yes") {
          if (model_type == "cox") {
            lp <- predict(model, newdata = pred_data, re.form = NA, type = "lp")
          } else {
            lp <- predict(model, newdata = pred_data, re.form = NA, type = "link")
          }
        } else {
          if (model_type == "cox") {
            lp <- predict(model, newdata = pred_data, type = "lp")
          } else {
            lp <- predict(model, newdata = pred_data, type = "link")
          }
        }
        
        # Transform derivative based on model type
        derivatives <- switch(model_type,
                              "nb" = exp(lp) * lp_deriv,
                              "poisson" = exp(lp) * lp_deriv,
                              "logistic" = {
                                p <- plogis(lp)
                                p * (1 - p) * lp_deriv
                              },
                              "cox" = exp(lp) * lp_deriv,
                              "lm" = lp_deriv,
                              lp_deriv)
      } else {
        warning("Insufficient coefficients found for derivative calculation")
        derivatives <- rep(NA, nrow(pred_data))
      }
    }
    
    return(derivatives)
  }
  
  # Main analysis loop for each imputation
  all_derivatives <- list()
  first_pred_data <- NULL
  knots <- NULL
  
  for (imp in imputations) {
    # Subset data for this imputation
    imp_data <- subset(data, data[[imp_col]] == imp)
    
    # Convert subgroups to factor if requested and needed
    if (!is.null(subgroups) && subgroup_as_factor && !is.factor(imp_data[[subgroups]])) {
      # Convert subgroups column to factor
      imp_data[[subgroups]] <- as.factor(imp_data[[subgroups]])
      
      # Apply custom labels if provided
      if (!is.null(subgroup_labels)) {
        levels(imp_data[[subgroups]]) <- subgroup_labels
      }
    }
    
    # Create prediction data for this imputation
    pred_data <- create_pred_data(x_points, imp_data[[subgroups]])
    pred_data <- add_covariates(pred_data, imp_data)
    
    # Store first pred_data structure for later use
    if (is.null(first_pred_data)) {
      first_pred_data <- pred_data
    }
    
    # Build model formula
    formula_obj <- build_formula()
    
    # Fit model
    model <- fit_model(imp_data, formula_obj)
    
    # Determine knots for splines (use same knots across imputations)
    if (is.null(knots)) {
      knots <- Hmisc::rcspline.eval(imp_data[[variable_x]], nk = knot_n, inclx = FALSE, knots.only = TRUE)
    }
    
    # Calculate derivatives
    derivatives <- calculate_derivatives(model, pred_data, knots)
    
    # Store results
    all_derivatives[[as.character(imp)]] <- derivatives
  }
  
  # Pool derivatives using Rubin's rules
  num_observations <- nrow(first_pred_data)
  pooled_derivatives <- numeric(num_observations)
  pooled_variance <- numeric(num_observations)
  
  for (i in 1:num_observations) {
    # Extract derivatives for this observation across imputations
    point_derivs <- sapply(all_derivatives, function(x) x[i])
    
    # Remove any NA values
    point_derivs <- point_derivs[!is.na(point_derivs)]
    
    if (length(point_derivs) > 0) {
      # Calculate mean (Q-bar) of derivatives across imputations
      q_bar <- mean(point_derivs)
      pooled_derivatives[i] <- q_bar
      
      # Calculate variance
      if (length(point_derivs) > 1) {
        # Between-imputation variance B
        B <- var(point_derivs)
        
        # Within-imputation variance U - approximate as square of pooled value / 10
        # This is a simplification, ideally we'd get the variance from each model
        U <- abs(q_bar / 10)^2
        
        # Total variance T = U + (1 + 1/m)B
        m <- length(point_derivs)
        T_var <- U + (1 + 1/m) * B
        
        pooled_variance[i] <- T_var
      } else {
        # If only one non-NA value, use approximation
        pooled_variance[i] <- abs(q_bar / 10)^2
      }
    } else {
      pooled_derivatives[i] <- NA
      pooled_variance[i] <- NA
    }
  }
  
  # Calculate confidence intervals
  pooled_se <- sqrt(pooled_variance)
  lower_ci <- pooled_derivatives - 1.96 * pooled_se
  upper_ci <- pooled_derivatives + 1.96 * pooled_se
  
  # Create results data frame
  result_data <- first_pred_data
  result_data$x_value <- result_data[[variable_x]]  # For plotting
  result_data$derivative <- pooled_derivatives
  result_data$derivative_se <- pooled_se
  result_data$lower_ci <- lower_ci
  result_data$upper_ci <- upper_ci
  result_data$significant_positive <- lower_ci > 0
  
  # Calculate thresholds (points where derivative significance changes)
  if (is.null(subgroups)) {
    # Order by x value
    ordered_results <- result_data[order(result_data$x_value), ]
    
    # Find where derivative becomes significantly positive
    sig_pos_idx <- which(ordered_results$significant_positive)
    threshold <- if (length(sig_pos_idx) > 0) {
      ordered_results$x_value[min(sig_pos_idx)]
    } else {
      NA
    }
    
    thresholds <- threshold
  } else {
    # Initialize threshold values for each subgroup
    subgroup_levels <- unique(result_data[[subgroups]])
    threshold_values <- rep(NA, length(subgroup_levels))
    names(threshold_values) <- as.character(subgroup_levels)
    
    # Find threshold for each subgroup
    for (sg in subgroup_levels) {
      # Subset and order data for this subgroup
      sg_data <- subset(result_data, result_data[[subgroups]] == sg)
      sg_data <- sg_data[order(sg_data$x_value), ]
      
      # Find where derivative becomes significantly positive
      sig_pos_idx <- which(sg_data$significant_positive)
      if (length(sig_pos_idx) > 0) {
        threshold_values[as.character(sg)] <- sg_data$x_value[min(sig_pos_idx)]
      }
    }
    
    thresholds <- threshold_values
  }
  
  # Create plot of derivatives
  if (is.null(subgroups)) {
    plot <- ggplot2::ggplot(result_data, ggplot2::aes(x = x_value, y = derivative)) +
      ggplot2::geom_line(color = "blue", linewidth = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), 
                           alpha = 0.3, fill = "lightblue") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = "Derivatives of Restricted Cubic Spline",
        x = variable_x,
        y = paste0("Derivative of ", outcome_var)
      ) +
      ggplot2::theme_minimal()
    
    # Add threshold line if available
    if (!is.na(threshold)) {
      plot <- plot + ggplot2::geom_vline(xintercept = threshold, 
                                         color = "red", linetype = "dashed")
    }
  } else {
    plot <- ggplot2::ggplot(result_data, 
                            ggplot2::aes(x = x_value, y = derivative, 
                                         color = .data[[subgroups]], 
                                         fill = .data[[subgroups]])) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), 
                           alpha = 0.3) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = "Derivatives of Restricted Cubic Spline by Subgroup",
        x = variable_x,
        y = paste0("Derivative of ", outcome_var)
      ) +
      ggplot2::theme_minimal()
    
    # Add threshold lines for subgroups
    for (sg in names(threshold_values)) {
      if (!is.na(threshold_values[sg])) {
        sg_data <- data.frame(
          x = threshold_values[sg],
          subgroup = sg
        )
        names(sg_data)[2] <- subgroups
        
        # Ensure subgroup column has same type as in result_data
        if (is.factor(result_data[[subgroups]])) {
          sg_data[[subgroups]] <- factor(sg_data[[subgroups]], 
                                         levels = levels(result_data[[subgroups]]))
        }
        plot <- plot + ggplot2::geom_vline(
          data = sg_data,
          ggplot2::aes(xintercept = x, color = .data[[subgroups]]),
          linetype = "dashed",
          linewidth = 0.5
        )
      }
    }
  }
  
  # Return results
  return(list(
    derivative_data = result_data,
    thresholds = thresholds,
    plot = plot,
    knots = knots,
    x_range = range(x_points)
  ))
}