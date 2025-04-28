#' Make predictions from multiple imputation models
#'
#' This function makes prediction based on the estimation of the coefficient and standard error obtained with the Rubin rule.
#' When subgroups are specified, separate models are fit for each subgroup.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for the model.
#' @param variable_x The continuous predictor variable to be modeled.
#' @param subgroups Optional variable for stratification (default is NULL).
#' @param covariables Optional vector of covariates to adjust for.
#' @param rcs_terms Optional list of restricted cubic spline specifications (default is NULL).
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param model_type The model type: "nb" for Negative Binomial, "poisson" for Poisson, "cox" for Cox, "bin"/"logistic" for Logistic, or "lm" for Linear.
#' @param followup_offset Whether to include offset for follow-up duration ("Yes" or "No").
#' @param followup_col The column representing follow-up duration (required if followup_offset = "Yes").
#' @param time_col The time variable for Cox regression (only required if model_type = "cox").
#' @param event_col The event variable for Cox regression (only required if model_type = "cox").
#' @param random_intercept Whether to include a random intercept in the model ("Yes" or "No").
#' @param random_intercept_var The grouping variable for the random intercept (required if random_intercept = "Yes").
#' @param plot_options List of options for plot customization.
#' @param prediction_range Range for predictions as quantiles (default is c(0.01, 0.99)).
#' @param n_points Number of points to use in the prediction sequence (default is 100).
#' @return A list containing predictions, plots, and model information for each subgroup.
#' @export

MI_predictions <- function(data,
                           outcome_var,
                           variable_x,
                           subgroups = NULL,
                           covariables = NULL,
                           rcs_terms = NULL,
                           imp_col = ".imp",
                           model_type = "lm",
                           followup_offset = "No",
                           followup_col = NULL,
                           time_col = NULL,
                           event_col = NULL,
                           random_intercept = "No",
                           random_intercept_var = NULL,
                           plot_options = list(
                             title = "Predicted Values with 95% Confidence Intervals",
                             x_label = NULL,
                             y_label = "Predicted Value",
                             line_color = "blue",
                             ribbon_fill = "blue",
                             ribbon_alpha = 0.2
                           ),
                           prediction_range = c(0.01, 0.99),
                           n_points = 100) {

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Data must be a data frame.")
  }

  if (!outcome_var %in% colnames(data)) {
    stop(paste("Outcome variable", outcome_var, "not found in data."))
  }

  if (!variable_x %in% colnames(data)) {
    stop(paste("Predictor variable", variable_x, "not found in data."))
  }

  if (!is.null(subgroups) && !subgroups %in% colnames(data)) {
    stop(paste("Subgroup variable", subgroups, "not found in data."))
  }

  if (!is.null(covariables)) {
    missing_covars <- covariables[!covariables %in% colnames(data)]
    if (length(missing_covars) > 0) {
      stop(paste("Covariates not found in data:", paste(missing_covars, collapse = ", ")))
    }
  }

  # Validate RCS terms if provided
  if (!is.null(rcs_terms)) {
    if (!is.list(rcs_terms)) {
      stop("rcs_terms must be a list where names are variable names and values are knot positions")
    }

    # Check if rms package is available
    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("The 'rms' package is needed for RCS terms. Please install it.")
    }

    # Check if variables in rcs_terms exist in data
    missing_vars <- names(rcs_terms)[!names(rcs_terms) %in% colnames(data)]
    if (length(missing_vars) > 0) {
      stop(paste("Variables specified in rcs_terms not found in data:",
                 paste(missing_vars, collapse = ", ")))
    }
  }

  # Adjust model_type if "logistic" is specified
  if (model_type == "logistic") model_type <- "bin"

  # Check for required columns based on model type
  if (model_type == "cox" && (is.null(time_col) || is.null(event_col))) {
    stop("For Cox models, time_col and event_col must be specified.")
  }

  if (followup_offset == "Yes" && is.null(followup_col)) {
    stop("followup_col must be specified when followup_offset is 'Yes'.")
  }

  if (random_intercept == "Yes" && is.null(random_intercept_var)) {
    stop("random_intercept_var must be specified when random_intercept is 'Yes'.")
  }

  # Importing required functions
  if (!exists("MI_estimates") || !exists("MI_vcov")) {
    stop("MI_estimates and MI_vcov functions must be available in the environment.")
  }

  # Create a function to process predictions for a given dataset and models
  process_predictions <- function(data, model_estimates, model_vcov, subgroup_value = NULL) {
    # Step 2: Create prediction data with sequence of x values
    x_seq <- seq(
      from = quantile(data[[variable_x]], prediction_range[1], na.rm = TRUE),
      to = quantile(data[[variable_x]], prediction_range[2], na.rm = TRUE),
      length.out = n_points
    )

    # Step 3: Create prediction dataframe with fixed values for other variables
    pred_data <- data.frame(x = x_seq)
    names(pred_data)[1] <- variable_x

    # Add median/mode values for covariables
    if (!is.null(covariables)) {
      for (cov in covariables) {
        if (is.numeric(data[[cov]])) {
          # For numeric variables, use median
          pred_data[[cov]] <- median(data[[cov]], na.rm = TRUE)
        } else if (is.factor(data[[cov]]) || is.character(data[[cov]])) {
          # For categorical variables, use most frequent value
          if (is.character(data[[cov]])) {
            data[[cov]] <- as.factor(data[[cov]])
          }
          majority_level <- names(which.max(table(data[[cov]])))
          pred_data[[cov]] <- factor(majority_level, levels = levels(data[[cov]]))
        } else {
          # For other types, use first value (with a warning)
          warning(paste("Variable", cov, "is of unsupported type. Using first value."))
          pred_data[[cov]] <- data[[cov]][1]
        }
      }
    }

    # Step 4: Handle the design matrix correctly
    # Get the coefficient names from the model
    coef_names <- model_estimates$term

    # Create a design matrix with the correct dimensions
    X_pred <- matrix(0, nrow = nrow(pred_data), ncol = length(coef_names))
    colnames(X_pred) <- coef_names

    # Fill in the intercept
    if("(Intercept)" %in% coef_names) {
      X_pred[, "(Intercept)"] <- 1
    }

    # Fill in the main predictor variable
    if(variable_x %in% coef_names) {
      X_pred[, variable_x] <- pred_data[[variable_x]]
    }

    # Handle RCS terms if present
    if (!is.null(rcs_terms)) {
      # Check for RCS terms in coefficient names
      for (var in names(rcs_terms)) {
        # Look for RCS terms like "rcs(var, knots)" or "var'"
        rcs_pattern <- paste0("^", var, "'")
        rcs_cols <- grep(rcs_pattern, coef_names, value = TRUE)

        if (length(rcs_cols) > 0) {
          # Require the rms package for RCS handling
          if (!requireNamespace("rms", quietly = TRUE)) {
            stop("The 'rms' package is needed for RCS terms. Please install it.")
          }

          # Get knots from rcs_terms
          knots <- rcs_terms[[var]]

          # Create rcs basis functions for prediction
          rcs_basis <- rms::rcspline.eval(pred_data[[var]], knots = knots, inclx = TRUE)
          colnames(rcs_basis) <- paste0(var, "'", seq_len(ncol(rcs_basis)))

          # Match column names with coefficient names
          for (i in seq_len(ncol(rcs_basis))) {
            col_name <- colnames(rcs_basis)[i]
            if (col_name %in% coef_names) {
              X_pred[, col_name] <- rcs_basis[, i]
            }
          }
        }
      }
    }

    # Fill in covariates
    for (var in covariables) {
      if (!is.null(var)) {
        if (var %in% coef_names) {
          X_pred[, var] <- pred_data[[var]]
        } else {
          # For categorical variables, we need to handle their dummy encoding
          if (is.factor(pred_data[[var]]) || is.character(pred_data[[var]])) {
            # Get the reference level (first level)
            ref_level <- if(is.factor(pred_data[[var]])) levels(pred_data[[var]])[1] else levels(factor(pred_data[[var]]))[1]

            # Get the value we're setting
            current_value <- as.character(pred_data[[var]][1])

            # Check for each possible column that could be from this variable
            for (coef in coef_names) {
              # Match pattern like "varNameLevelName"
              if (grepl(paste0("^", var), coef)) {
                level_name <- sub(paste0("^", var), "", coef)
                if (!is.na(level_name) && level_name != "") {
                  X_pred[, coef] <- ifelse(current_value == level_name, 1, 0)
                }
              }
            }
          }
        }
      }
    }

    # Step 5: Calculate predicted values
    beta_hat <- model_estimates$estimate

    # Calculate linear predictor
    linear_pred <- as.numeric(X_pred %*% beta_hat)

    # Transform predictions based on model type
    if (model_type %in% c("bin", "logistic")) {
      # Logistic regression - transform to probabilities
      pred_data$y_pred <- plogis(linear_pred)
    } else if (model_type %in% c("nb", "poisson")) {
      # Count models - transform to counts
      pred_data$y_pred <- exp(linear_pred)

      # Apply offset if specified
      if (followup_offset == "Yes" && !is.null(followup_col)) {
        # Use median follow-up time to adjust predictions
        median_followup <- median(data[[followup_col]], na.rm = TRUE)
        pred_data$y_pred <- pred_data$y_pred * median_followup
      }
    } else {
      # Linear model - no transformation needed
      pred_data$y_pred <- linear_pred
    }

    # Step 6: Calculate confidence intervals
    conf_level <- 0.95
    alpha <- 1 - conf_level
    z_crit <- qnorm(1 - alpha/2)

    # Calculate standard errors for the predictions
    pred_se <- numeric(nrow(X_pred))
    for (i in 1:nrow(X_pred)) {
      # Extract the row of the design matrix for this prediction
      x_i <- X_pred[i,]

      # Calculate the variance of the prediction
      pred_var <- t(x_i) %*% model_vcov %*% x_i

      # Standard error is the square root of the variance
      pred_se[i] <- sqrt(as.numeric(pred_var))
    }

    # Calculate the confidence intervals on the linear predictor scale
    lower_linear <- linear_pred - z_crit * pred_se
    upper_linear <- linear_pred + z_crit * pred_se

    # Transform confidence intervals based on model type
    if (model_type %in% c("bin", "logistic")) {
      pred_data$lower_ci <- plogis(lower_linear)
      pred_data$upper_ci <- plogis(upper_linear)
    } else if (model_type %in% c("nb", "poisson")) {
      pred_data$lower_ci <- exp(lower_linear)
      pred_data$upper_ci <- exp(upper_linear)

      # Apply offset if specified
      if (followup_offset == "Yes" && !is.null(followup_col)) {
        median_followup <- median(data[[followup_col]], na.rm = TRUE)
        pred_data$lower_ci <- pred_data$lower_ci * median_followup
        pred_data$upper_ci <- pred_data$upper_ci * median_followup
      }
    } else {
      # Linear model - no transformation needed
      pred_data$lower_ci <- lower_linear
      pred_data$upper_ci <- upper_linear
    }

    # Step 7: Create plot
    library(ggplot2)

    # Set default labels if not provided
    if (is.null(plot_options$x_label)) {
      plot_options$x_label <- variable_x
    }

    # Create caption with fixed values and subgroup information
    caption_parts <- c()

    # Add subgroup information if available
    if (!is.null(subgroup_value)) {
      caption_parts <- c(caption_parts, paste0(subgroups, " = ", subgroup_value))
    }

    # Add covariate information
    for (var in covariables) {
      if (!is.null(var)) {
        value <- if (is.numeric(data[[var]])) {
          round(median(data[[var]], na.rm = TRUE), 2)
        } else {
          as.character(pred_data[[var]][1])
        }
        caption_parts <- c(caption_parts, paste0(var, " = ", value))
      }
    }

    # Format caption
    if (length(caption_parts) > 0) {
      caption_text <- paste("Variables fixed at:", paste(caption_parts, collapse = ", "))
    } else {
      caption_text <- ""
    }

    # Create plot title with subgroup information
    if (!is.null(subgroup_value)) {
      main_title <- paste0(plot_options$title, " (", subgroups, " = ", subgroup_value, ")")
    } else {
      main_title <- plot_options$title
    }

    # Create the plot
    prediction_plot <- ggplot(pred_data, aes_string(x = variable_x, y = "y_pred")) +
      geom_line(color = plot_options$line_color, size = 1) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
                  alpha = plot_options$ribbon_alpha,
                  fill = plot_options$ribbon_fill) +
      labs(
        title = main_title,
        subtitle = paste0("Based on Multiple Imputation (", n_points, " points)"),
        x = plot_options$x_label,
        y = plot_options$y_label,
        caption = caption_text
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "italic"),
        plot.caption = element_text(hjust = 0)
      )

    # Return results
    return(list(
      predictions = pred_data,
      plot = prediction_plot,
      design_matrix = X_pred,
      coefficients = beta_hat,
      model_type = model_type,
      subgroup_value = subgroup_value
    ))
  }

  # Create formula arguments as a function to reuse for each subgroup
  create_formula_args <- function(data_subset, include_subgroups = TRUE) {
    # Determine which predictors to include
    subset_predictors <- c(variable_x, covariables)
    if (include_subgroups && !is.null(subgroups)) {
      subset_predictors <- c(subset_predictors, subgroups)
    }

    # Create the formula arguments
    args <- list(
      data = data_subset,
      outcome_var = outcome_var,
      predictor_vars = subset_predictors,
      imp_col = imp_col,
      model_type = model_type
    )

    # Handle RCS terms if provided
    if (!is.null(rcs_terms)) {
      args$spline_terms <- rcs_terms
    }

    # Add optional arguments
    if (followup_offset == "Yes") {
      args$followup_offset <- "Yes"
      args$followup_col <- followup_col
    }

    if (model_type == "cox") {
      args$time_col <- time_col
      args$event_col <- event_col
    }

    if (random_intercept == "Yes") {
      args$random_intercept <- "Yes"
      args$random_intercept_var <- random_intercept_var
    }

    return(args)
  }

  # Main processing based on whether subgroups are specified
  if (!is.null(subgroups)) {
    # Get unique subgroup values
    subgroup_values <- unique(data[[subgroups]])

    # First, fit an overall model including subgroups for comparison
    overall_formula_args <- create_formula_args(data, include_subgroups = TRUE)
    overall_MI_estimates <- do.call(MI_estimates, overall_formula_args)
    overall_MI_vcov <- do.call(MI_vcov, overall_formula_args)

    # Process each subgroup with a separate model
    results <- lapply(subgroup_values, function(val) {
      # Extract data for this subgroup
      subgroup_data <- data[data[[subgroups]] == val, ]

      # Create formula arguments for this subgroup
      # (subgroups parameter is excluded since we're already filtering)
      subgroup_formula_args <- create_formula_args(subgroup_data, include_subgroups = FALSE)

      # Fit models specific to this subgroup
      subgroup_MI_estimates <- do.call(MI_estimates, subgroup_formula_args)
      subgroup_MI_vcov <- do.call(MI_vcov, subgroup_formula_args)

      # Process predictions for this subgroup
      subgroup_results <- process_predictions(
        subgroup_data,
        subgroup_MI_estimates,
        subgroup_MI_vcov,
        val
      )

      # Add the model information
      subgroup_results$model_estimates <- subgroup_MI_estimates
      subgroup_results$model_vcov <- subgroup_MI_vcov

      return(subgroup_results)
    })
    names(results) <- subgroup_values

    # Also include the overall model
    results$overall <- process_predictions(data, overall_MI_estimates, overall_MI_vcov)
    results$overall$model_estimates <- overall_MI_estimates
    results$overall$model_vcov <- overall_MI_vcov

    # Add model comparison information
    model_summary <- list(
      overall_model = list(
        formula = paste(outcome_var, "~", paste(c(variable_x, covariables, subgroups), collapse=" + ")),
        coefficients = overall_MI_estimates
      ),
      subgroup_models = lapply(subgroup_values, function(val) {
        list(
          subgroup = val,
          formula = paste(outcome_var, "~", paste(c(variable_x, covariables), collapse=" + ")),
          coefficients = results[[val]]$model_estimates
        )
      })
    )
    results$model_summary <- model_summary

    return(results)

  } else {
    # No subgroups - just fit one model for all data
    formula_args <- create_formula_args(data)
    model_MI_estimates <- do.call(MI_estimates, formula_args)
    model_MI_vcov <- do.call(MI_vcov, formula_args)

    # Process predictions for entire dataset
    results <- process_predictions(data, model_MI_estimates, model_MI_vcov)
    results$model_estimates <- model_MI_estimates
    results$model_vcov <- model_MI_vcov

    return(results)
  }
}
