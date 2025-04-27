#' Make predictions from multiple imputation models
#'
#' This function make prediction based on the estimation of the coefficient and standard error obtained with the Rubin rule.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for the model.
#' @param variable_x The continuous predictor variable to be modeled with splines.
#' @param subgroups Optional variable for stratification (default is NULL).
#' @param covariables Optional vector of covariates to adjust for.
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
#' @return A list containing predictions, plot, and model information.
#' @export

MI_predictions <- function(data,
                           outcome_var,
                           variable_x,
                           subgroups = NULL,
                           covariables = NULL,
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

  # Construct predictor variables vector
  predictor_vars <- c(variable_x)
  if (!is.null(covariables)) {
    predictor_vars <- c(predictor_vars, covariables)
  }
  if (!is.null(subgroups)) {
    predictor_vars <- c(predictor_vars, subgroups)
  }

  # Step 1: Calculate the pooled estimation and variance-covariance matrix
  # Importing required functions (assuming MI_estimates and MI_vcov are available)
  if (!exists("MI_estimates") || !exists("MI_vcov")) {
    stop("MI_estimates and MI_vcov functions must be available in the environment.")
  }

  # Creating model formula
  formula_args <- list(
    outcome_var = outcome_var,
    predictor_vars = predictor_vars,
    imp_col = imp_col,
    model_type = model_type
  )

  # Add optional arguments based on user input
  if (followup_offset == "Yes") {
    formula_args$followup_offset <- "Yes"
    formula_args$followup_col <- followup_col
  }

  if (model_type == "cox") {
    formula_args$time_col <- time_col
    formula_args$event_col <- event_col
  }

  if (random_intercept == "Yes") {
    formula_args$random_intercept <- "Yes"
    formula_args$random_intercept_var <- random_intercept_var
  }

  # Make sure data is included in the arguments
  formula_args$data <- data

  # Call MI_estimates and MI_vcov with the constructed arguments
  model_MI_estimates <- do.call(MI_estimates, formula_args)
  model_MI_vcov <- do.call(MI_vcov, formula_args)

  # Create a function to process each subgroup
  process_subgroup <- function(data, subgroup_value = NULL) {
    # Filter data if needed
    if (!is.null(subgroups) && !is.null(subgroup_value)) {
      data <- data[data[[subgroups]] == subgroup_value, ]
    }

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

    # Add subgroup if we're processing all data together
    if (!is.null(subgroups) && is.null(subgroup_value)) {
      majority_level <- names(which.max(table(data[[subgroups]])))
      pred_data[[subgroups]] <- factor(majority_level, levels = levels(data[[subgroups]]))
    }

    # Step 4: Handle the design matrix correctly
    # Get the coefficient names from the model
    coef_names <- model_MI_estimates$term

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

    # Fill in covariates
    for (var in c(covariables, subgroups)) {
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
    beta_hat <- model_MI_estimates$estimate

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
      pred_var <- t(x_i) %*% model_MI_vcov %*% x_i

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

    # Create caption with fixed values
    caption_parts <- c()
    for (var in c(covariables, subgroups)) {
      if (!is.null(var)) {
        value <- if (is.numeric(data[[var]])) {
          round(median(data[[var]], na.rm = TRUE), 2)
        } else {
          as.character(pred_data[[var]][1])
        }
        caption_parts <- c(caption_parts, paste0(var, " = ", value))
      }
    }
    caption_text <- paste("Other variables fixed at:", paste(caption_parts, collapse = ", "))

    # Create the plot
    prediction_plot <- ggplot(pred_data, aes_string(x = variable_x, y = "y_pred")) +
      geom_line(color = plot_options$line_color, size = 1) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
                  alpha = plot_options$ribbon_alpha,
                  fill = plot_options$ribbon_fill) +
      labs(
        title = plot_options$title,
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
      model_type = model_type
    ))
  }

  # Process data based on whether subgroups are specified
  if (!is.null(subgroups)) {
    # Get unique subgroup values
    subgroup_values <- unique(data[[subgroups]])

    # Process each subgroup
    results <- lapply(subgroup_values, function(val) {
      process_subgroup(data, val)
    })
    names(results) <- subgroup_values

    # Also process the whole dataset with a representative subgroup
    results$overall <- process_subgroup(data)

    # Return results
    results$model_estimates <- model_MI_estimates
    results$model_vcov <- model_MI_vcov

    return(results)
  } else {
    # Process the whole dataset
    results <- process_subgroup(data)
    results$model_estimates <- model_MI_estimates
    results$model_vcov <- model_MI_vcov

    return(results)
  }
}
