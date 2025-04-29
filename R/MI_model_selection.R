#' Model Selection for Multiple Imputed Data
#'
#' This function performs systematic model selection on multiply imputed data.
#' It handles both core model selection and additional covariable selection
#' through forward, backward, or brute force approaches.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for the model.
#' @param core_vars_list List of core variables or combinations to test for the base model.
#'        Can be a list of character vectors, each representing a possible core model.
#' @param potential_vars_list Character vector of additional covariables to test.
#' @param imp_col The column name indicating imputation indices (default is ".imp").
#' @param followup_offset Whether to include an offset for follow-up duration ("Yes" or "No").
#' @param followup_col The column representing follow-up duration (required if followup_offset = "Yes").
#' @param trial_factor Whether to include trial as a factor ("Yes" or "No").
#' @param trial_col The column for trial factor (required if trial_factor = "Yes").
#' @param imp_n Number of imputations (if NULL, detected automatically).
#' @param model_type The model type: "nb", "poisson", "lm", or "bin"/"logistic".
#' @param random_intercept Whether to include a random intercept in the model ("Yes" or "No").
#' @param random_intercept_var The grouping variable for the random intercept (required if random_intercept = "Yes").
#' @param core_selection_criteria Criteria for selecting core model: "AIC", "AICc", "BIC", "BICc", or "LRT".
#' @param potential_selection_strategy Strategy for additional covariable selection: "forward", "backward", or "brute_force".
#' @param potential_selection_term Term to evaluate for covariable selection: "p_value", "AIC", "AICc", "BIC", or "BICc".
#' @param potential_selection_threshold Threshold for covariable selection (e.g., 0.05 for p-value).
#' @param max_models Maximum number of models to evaluate (to prevent computational explosion in brute force).
#' @param verbose Whether to display progress information (default: TRUE).
#'
#' @return A list containing:
#'   \item{best_model}{The final selected model output}
#'   \item{core_model_selection}{Results from core model selection}
#'   \item{covariable_selection}{Results from additional covariable selection}
#'   \item{model_comparison}{Full comparison of tested models}
#'   \item{final_formula}{Formula of the final selected model}
#'
#' @export
#' @importFrom stats as.formula
MI_model_selection <- function(
    data,
    outcome_var,
    core_vars_list,
    potential_vars_list = NULL,
    imp_col = ".imp",
    followup_offset = "No",
    followup_col = NULL,
    trial_factor = "No",
    trial_col = NULL,
    imp_n = NULL,
    model_type = "nb",
    random_intercept = "No",
    random_intercept_var = NULL,
    core_selection_criteria = "AICc",
    potential_selection_strategy = "forward",
    potential_selection_term = "p_value",
    potential_selection_threshold = 0.05,
    max_models = 1000,
    verbose = TRUE
) {
  # Validate inputs
  if (!exists("MI_model_performance")) {
    stop("The MI_model_performance function is required but not found.")
  }

  if (!exists("MI_models_comparison")) {
    stop("The MI_models_comparison function is required but not found.")
  }

  if (!core_selection_criteria %in% c("AIC", "AICc", "BIC", "BICc", "LRT")) {
    stop("core_selection_criteria must be one of: 'AIC', 'AICc', 'BIC', 'BICc', or 'LRT'")
  }

  if (!is.null(potential_vars_list)) {
    if (!potential_selection_strategy %in% c("forward", "backward", "brute_force")) {
      stop("potential_selection_strategy must be one of: 'forward', 'backward', or 'brute_force'")
    }

    if (!potential_selection_term %in% c("p_value", "AIC", "AICc", "BIC", "BICc")) {
      stop("potential_selection_term must be one of: 'p_value', 'AIC', 'AICc', 'BIC', or 'BICc'")
    }
  }

  # Generate all combinations of core variables
  # For each named group in core_vars_list, select one option
  generate_core_combinations <- function(vars_list) {
    # If it's empty, return empty list
    if (length(vars_list) == 0) {
      return(list(character(0)))
    }

    # Get the names of the groups
    group_names <- names(vars_list)

    # If no names provided, create default names
    if (is.null(group_names) || any(group_names == "")) {
      group_names <- paste0("Group", seq_along(vars_list))
      names(vars_list) <- group_names
    }

    # Generate all combinations using expand.grid
    # First, create a list where each element is a factor of options for a group
    group_options <- lapply(vars_list, function(options) {
      # Ensure each option within a group is a character
      if (is.list(options)) {
        options <- unlist(options)
      }
      factor(seq_along(options))
    })

    # Generate all combinations using expand.grid
    combinations <- expand.grid(group_options, stringsAsFactors = FALSE)

    # Convert these indices back to actual variable names
    result_combinations <- list()
    for (i in 1:nrow(combinations)) {
      row_vars <- c()
      for (j in 1:ncol(combinations)) {
        option_idx <- as.numeric(combinations[i, j])
        group_vars <- vars_list[[j]]

        # Handle nested lists if present
        if (is.list(group_vars)) {
          group_vars <- unlist(group_vars)
        }

        selected_var <- group_vars[option_idx]
        row_vars <- c(row_vars, selected_var)
      }
      result_combinations[[i]] <- row_vars
    }

    return(result_combinations)
  }

  # Process potential_vars_list to flatten if needed
  if (!is.null(potential_vars_list)) {
    if (!is.list(potential_vars_list)) {
      # Single vector provided directly
      potential_vars_list <- list(potential_vars_list)
    }

    # Check if potential_vars_list contains nested lists
    if (any(sapply(potential_vars_list, is.list))) {
      # Flatten nested lists
      potential_vars_flat <- unlist(potential_vars_list, recursive = TRUE)
      potential_vars_list <- potential_vars_flat
    } else {
      # Just unlist to get a character vector
      potential_vars_list <- unlist(potential_vars_list)
    }

    # Ensure it's a character vector
    if (!is.character(potential_vars_list)) {
      stop("Potential variables list must contain character vectors.")
    }
  }

  # Generate all possible combinations of core variable expressions
  core_combinations <- generate_core_combinations(core_vars_list)

  # Function to build a model and evaluate performance
  build_model <- function(variables, description = NULL) {

    # Add offset parameter for count models
    follow_offset_param <- followup_offset
    if (model_type %in% c("nb", "poisson") && !is.null(followup_col)) {
      follow_offset_param <- "Yes"
    }

    result <- MI_model_performance(
      data = data,
      outcome_var = outcome_var,
      predictor_vars = variables,
      imp_col = imp_col,
      followup_offset = follow_offset_param,
      followup_col = followup_col,
      trial_factor = trial_factor,
      trial_col = trial_col,
      imp_n = imp_n,
      model_type = model_type,
      random_intercept = random_intercept,
      random_intercept_var = random_intercept_var,
      verbose = FALSE
    )

    return(result)
  }

  #####################################################################
  # STEP 1: Core Model Selection
  #####################################################################
  if (verbose) {
    cat("\n=== STEP 1: CORE MODEL SELECTION ===\n")
    cat(sprintf("Generated %d core model combinations\n", length(core_combinations)))
    cat(sprintf("Evaluating core models using %s criterion\n", core_selection_criteria))
  }

  # Evaluate all core models
  core_models <- list()
  model_names <- character(length(core_combinations))

  for (i in seq_along(core_combinations)) {
    core_vars <- core_combinations[[i]]
    model_name <- paste("Core", i)
    model_names[i] <- model_name

    if (verbose) {
      cat(sprintf("  Model %s: %s\n", model_name, paste(core_vars, collapse = ", ")))
    }

    core_models[[model_name]] <- build_model(core_vars)

    # Print selection criteria value for this model
    if (verbose && core_selection_criteria != "LRT") {
      cat(sprintf("    %s: %.2f\n", core_selection_criteria,
                  core_models[[model_name]][[core_selection_criteria]]))
    }
  }

  # Store model variables for each core model to track them properly
  core_model_variables <- list()
  for (i in seq_along(core_combinations)) {
    model_name <- paste("Core", i)
    core_model_variables[[model_name]] <- core_combinations[[i]]
  }

  # Compare core models
  core_models_comparison <- do.call(
    MI_models_comparison,
    c(core_models, list(model_names = model_names, sort_by = core_selection_criteria))
  )

  # Select best core model based on specified criterion
  if (core_selection_criteria %in% c("AIC", "AICc", "BIC", "BICc")) {
    criterion_col <- paste0("Delta_", core_selection_criteria)
    best_core_idx <- which.min(core_models_comparison[[criterion_col]])
  } else if (core_selection_criteria == "LRT") {
    # For LRT, more complex models with non-significant p-values should be preferred
    # We select the model with the most variables that has LRT p > 0.05 compared to reference
    significant_models <- which(as.numeric(core_models_comparison$LRT) >= 0.05)
    if (length(significant_models) > 0) {
      # Find model with most variables among significant ones
      var_counts <- sapply(significant_models, function(idx) {
        model_name <- core_models_comparison$Model_Name[idx]
        return(length(core_model_variables[[model_name]]))
      })
      best_core_idx <- significant_models[which.max(var_counts)]
    } else {
      # If all are significant, use the most complex model
      var_counts <- sapply(seq_along(core_models_comparison$Model_Name), function(idx) {
        model_name <- core_models_comparison$Model_Name[idx]
        return(length(core_model_variables[[model_name]]))
      })
      best_core_idx <- which.max(var_counts)
    }
  }

  best_core_model_name <- core_models_comparison$Model_Name[best_core_idx]
  best_core_model <- core_models[[best_core_model_name]]
  best_core_vars <- core_model_variables[[best_core_model_name]]

  if (verbose) {
    cat(sprintf("\nBest core model: %s\n", best_core_model_name))
    cat(sprintf("Variables: %s\n", paste(best_core_vars, collapse = ", ")))
    cat(sprintf("%s: %.2f\n", core_selection_criteria,
                if (core_selection_criteria == "LRT")
                  as.numeric(core_models_comparison$LRT[best_core_idx])
                else best_core_model[[core_selection_criteria]]))

    # Print Delta values
    if (core_selection_criteria %in% c("AIC", "AICc", "BIC", "BICc")) {
      cat("\nComparison of all core models:\n")
      delta_col <- paste0("Delta_", core_selection_criteria)
      print_df <- data.frame(
        Model = core_models_comparison$Model_Name,
        Criterion = round(core_models_comparison[[core_selection_criteria]], 2),
        Delta = round(core_models_comparison[[delta_col]], 2)
      )
      print(print_df, row.names = FALSE)
    }
  }

  # If no potential variables provided, return best core model as final result
  if (is.null(potential_vars_list) || length(potential_vars_list) == 0) {
    if (verbose) {
      cat("\nNo potential covariables provided. Returning best core model as final result.\n")
    }

    return(list(
      best_model = best_core_model,
      core_model_selection = core_models_comparison,
      covariable_selection = NULL,
      model_comparison = core_models_comparison,
      final_formula = best_core_model$Model_Formula
    ))
  }

  #####################################################################
  # STEP 2: Additional Covariable Selection
  #####################################################################
  if (verbose) {
    cat("\n=== STEP 2: ADDITIONAL COVARIABLE SELECTION ===\n")
    cat(sprintf("Strategy: %s, Criterion: %s, Threshold: %.5f\n",
                potential_selection_strategy, potential_selection_term,
                potential_selection_threshold))
  }

  # Remove any potential variables that are already in the core model
  potential_vars_clean <- setdiff(potential_vars_list, best_core_vars)

  if (length(potential_vars_clean) == 0) {
    if (verbose) {
      cat("No additional potential covariables available (all already in core model).\n")
    }

    return(list(
      best_model = best_core_model,
      core_model_selection = core_models_comparison,
      covariable_selection = NULL,
      model_comparison = core_models_comparison,
      final_formula = best_core_model$Model_Formula
    ))
  }

  if (verbose) {
    cat(sprintf("Evaluating %d potential covariables: %s\n",
                length(potential_vars_clean),
                paste(potential_vars_clean, collapse = ", ")))
  }

  # Strategy-specific implementation
  covariable_models <- list()
  covariable_selection_steps <- list()

  if (potential_selection_strategy == "forward") {
    #####################################################################
    # Forward Selection
    #####################################################################
    current_model <- best_core_model
    current_vars <- best_core_vars
    remaining_vars <- potential_vars_clean
    step_counter <- 1
    improvement_found <- TRUE

    covariable_models[["Base"]] <- current_model

    while (improvement_found && length(remaining_vars) > 0) {
      improvement_found <- FALSE
      step_models <- list()
      step_vars <- list()

      # Add current model as baseline for comparison
      base_model_name <- paste("Step", step_counter, "base")
      step_models[[base_model_name]] <- current_model
      step_vars[[base_model_name]] <- current_vars

      # Try adding each remaining variable
      for (var in remaining_vars) {
        test_vars <- c(current_vars, var)
        model_name <- paste("Step", step_counter, "+", var)

        if (verbose) {
          cat(sprintf("Evaluating model with: Forward %s\n", model_name))
        }

        test_model <- build_model(test_vars,
                                  description = paste("Forward", model_name))

        step_models[[model_name]] <- test_model
        step_vars[[model_name]] <- test_vars
      }

      # Directly compare all models to base model
      base_model <- step_models[[base_model_name]]

      # Create comparison dataframe
      model_comparison <- data.frame(
        Model = names(step_models),
        stringsAsFactors = FALSE
      )

      # Add criterion values (AIC, BIC, etc)
      for (criterion in c("AIC", "AICc", "BIC", "BICc")) {
        model_comparison[[criterion]] <- sapply(step_models, function(m) m[[criterion]])
        model_comparison[[paste0("Delta_", criterion)]] <-
          model_comparison[[criterion]] - model_comparison[[criterion]][model_comparison$Model == base_model_name]
      }

      # Calculate p-values from LRT if available
      if ("logL" %in% names(base_model)) {
        model_comparison$LRT <- NA
        base_logL <- base_model$logL
        base_df <- base_model$df

        for (i in seq_along(step_models)) {
          if (names(step_models)[i] != base_model_name) {
            test_model <- step_models[[i]]
            # LRT statistic
            test_logL <- test_model$logL
            test_df <- test_model$df
            df_diff <- abs(test_df - base_df)
            if (df_diff == 0) df_diff <- 1

            # Calculate p-value
            lrt_stat <- 2 * abs(test_logL - base_logL)
            p_value <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
            model_comparison$LRT[i] <- p_value
          } else {
            model_comparison$LRT[i] <- 1.0  # Base model compared to itself
          }
        }
        # Format p-values
        model_comparison$LRT <- format.pval(model_comparison$LRT, digits = 5)
      }

      # Display comparison results
      if (verbose) {
        cat("  Comparison of all models from this step:\n")
        display_df <- data.frame(
          Model = model_comparison$Model,
          Criterion = round(model_comparison[[potential_selection_term]], 2),
          Delta = round(model_comparison[[paste0("Delta_", potential_selection_term)]], 2)
        )
        colnames(display_df)[2] <- potential_selection_term
        colnames(display_df)[3] <- paste0("Δ", potential_selection_term)

        # Highlight the base model
        base_idx <- which(display_df$Model == base_model_name)
        display_df$Model[base_idx] <- paste0(display_df$Model[base_idx], " (reference)")

        print(display_df, row.names = FALSE)
        cat("\n")
      }

      # Determine if any model is an improvement
      if (potential_selection_term == "p_value") {
        # Skip the base model for comparison
        p_idx <- which(model_comparison$Model == base_model_name)
        p_values <- as.numeric(model_comparison$LRT[-p_idx])

        if (any(p_values < potential_selection_threshold, na.rm = TRUE)) {
          # Find the best model by p-value
          best_p_idx <- which.min(p_values)
          best_idx <- setdiff(seq_along(model_comparison$Model), p_idx)[best_p_idx]
          improvement_found <- TRUE
        }
      } else {
        # For information criteria, check for improvement
        term_col <- paste0("Delta_", potential_selection_term)
        delta_idx <- which(model_comparison$Model == base_model_name)
        delta_values <- model_comparison[[term_col]][-delta_idx]

        if (any(delta_values < -potential_selection_threshold, na.rm = TRUE)) {
          # Find the best model by delta
          best_delta_idx <- which.min(delta_values)
          best_idx <- setdiff(seq_along(model_comparison$Model), delta_idx)[best_delta_idx]
          improvement_found <- TRUE
        }
      }

      # If an improvement was found, update the model
      if (improvement_found) {
        best_model_name <- model_comparison$Model[best_idx]

        # Extract the variable name from the model name
        added_var <- sub(paste0("Step ", step_counter, " \\+ "), "", best_model_name)

        if (verbose) {
          cat(sprintf("  Step %d: Added variable '%s'", step_counter, added_var))

          if (potential_selection_term == "p_value") {
            cat(sprintf(" (p-value: %s)\n", model_comparison$LRT[best_idx]))
          } else {
            term_col <- paste0("Delta_", potential_selection_term)
            cat(sprintf(" (%s: %.2f, Δ%s: %.2f)\n",
                        potential_selection_term,
                        model_comparison[[potential_selection_term]][best_idx],
                        potential_selection_term,
                        model_comparison[[term_col]][best_idx]))
          }
        }

        # Update current model and variables
        current_model <- step_models[[best_model_name]]
        current_vars <- step_vars[[best_model_name]]
        remaining_vars <- setdiff(remaining_vars, added_var)

        # Save model and step information
        covariable_models[[paste("Step", step_counter)]] <- current_model
        covariable_selection_steps[[step_counter]] <- list(
          step = step_counter,
          added_var = added_var,
          model_comparison = model_comparison,
          best_model_name = best_model_name
        )

        step_counter <- step_counter + 1
      } else {
        if (verbose) {
          cat("  No further improvements found.\n")
        }
      }
    }

  } else if (potential_selection_strategy == "backward") {
    #####################################################################
    # Backward Selection
    #####################################################################
    # Start with full model (core + all potential vars)
    full_vars <- c(best_core_vars, potential_vars_clean)

    # Check if the number of variables is reasonable
    if (length(full_vars) > 20) {
      warning("Backward selection with a large number of variables may be computationally intensive.")
    }

    full_model <- build_model(full_vars, description = "Full model (backward start)")
    current_model <- full_model
    current_vars <- full_vars
    step_counter <- 1

    covariable_models[["Full"]] <- current_model

    while (length(current_vars) > length(best_core_vars)) {
      step_models <- list()
      step_vars <- list()

      # Add current model as baseline for comparison
      base_model_name <- paste("Step", step_counter, "current")
      step_models[[base_model_name]] <- current_model
      step_vars[[base_model_name]] <- current_vars

      # Try removing each potential variable one at a time
      removable_vars <- setdiff(current_vars, best_core_vars)

      for (var in removable_vars) {
        test_vars <- setdiff(current_vars, var)
        model_name <- paste("Step", step_counter, "-", var)
        test_model <- build_model(test_vars,
                                  description = paste("Backward", model_name))
        step_models[[model_name]] <- test_model
        step_vars[[model_name]] <- test_vars
      }

      # Directly compare all models to base model
      base_model <- step_models[[base_model_name]]

      # Create comparison dataframe
      model_comparison <- data.frame(
        Model = names(step_models),
        stringsAsFactors = FALSE
      )

      # Add criterion values (AIC, BIC, etc)
      for (criterion in c("AIC", "AICc", "BIC", "BICc")) {
        model_comparison[[criterion]] <- sapply(step_models, function(m) m[[criterion]])
        model_comparison[[paste0("Delta_", criterion)]] <-
          model_comparison[[criterion]] - model_comparison[[criterion]][model_comparison$Model == base_model_name]
      }

      # Calculate p-values from LRT if available
      if ("logL" %in% names(base_model)) {
        model_comparison$LRT <- NA
        base_logL <- base_model$logL
        base_df <- base_model$df

        for (i in seq_along(step_models)) {
          if (names(step_models)[i] != base_model_name) {
            test_model <- step_models[[i]]
            # LRT statistic
            test_logL <- test_model$logL
            test_df <- test_model$df
            df_diff <- abs(test_df - base_df)
            if (df_diff == 0) df_diff <- 1

            # Calculate p-value
            lrt_stat <- 2 * abs(test_logL - base_logL)
            p_value <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
            model_comparison$LRT[i] <- p_value
          } else {
            model_comparison$LRT[i] <- 1.0  # Base model compared to itself
          }
        }
        # Format p-values
        model_comparison$LRT <- format.pval(model_comparison$LRT, digits = 5)
      }

      # Display comparison results
      if (verbose) {
        cat("  Comparison of all models from this step:\n")
        display_df <- data.frame(
          Model = model_comparison$Model,
          Criterion = round(model_comparison[[potential_selection_term]], 2),
          Delta = round(model_comparison[[paste0("Delta_", potential_selection_term)]], 2)
        )
        colnames(display_df)[2] <- potential_selection_term
        colnames(display_df)[3] <- paste0("Δ", potential_selection_term)

        # Highlight the base model
        base_idx <- which(display_df$Model == base_model_name)
        display_df$Model[base_idx] <- paste0(display_df$Model[base_idx], " (reference)")

        print(display_df, row.names = FALSE)
        cat("\n")
      }

      # Determine if we can remove a variable
      can_remove <- FALSE
      best_idx <- NULL

      if (potential_selection_term == "p_value") {
        # For p-value, we want variables that can be removed without significantly affecting the model
        p_idx <- which(model_comparison$Model == base_model_name)
        p_values <- as.numeric(model_comparison$LRT[-p_idx])

        if (any(p_values > potential_selection_threshold, na.rm = TRUE)) {
          # Find the variable with the highest p-value (least significant)
          best_p_idx <- which.max(p_values)
          best_idx <- setdiff(seq_along(model_comparison$Model), p_idx)[best_p_idx]
          can_remove <- TRUE
        }
      } else {
        # For information criteria, we want variables that, when removed, don't worsen the model too much
        term_col <- paste0("Delta_", potential_selection_term)
        delta_idx <- which(model_comparison$Model == base_model_name)
        delta_values <- model_comparison[[term_col]][-delta_idx]

        if (any(delta_values < potential_selection_threshold, na.rm = TRUE)) {
          # Find the model with the best (most negative or least positive) delta
          best_delta_idx <- which.min(delta_values)
          best_idx <- setdiff(seq_along(model_comparison$Model), delta_idx)[best_delta_idx]
          can_remove <- TRUE
        }
      }

      # If we can remove a variable, update the model
      if (can_remove) {
        best_model_name <- model_comparison$Model[best_idx]

        # Extract the variable name from the model name
        removed_var <- sub(paste0("Step ", step_counter, " - "), "", best_model_name)

        if (verbose) {
          cat(sprintf("  Step %d: Removed variable '%s'", step_counter, removed_var))

          if (potential_selection_term == "p_value") {
            cat(sprintf(" (p-value: %s)\n", model_comparison$LRT[best_idx]))
          } else {
            term_col <- paste0("Delta_", potential_selection_term)
            cat(sprintf(" (%s: %.2f, Δ%s: %.2f)\n",
                        potential_selection_term,
                        model_comparison[[potential_selection_term]][best_idx],
                        potential_selection_term,
                        model_comparison[[term_col]][best_idx]))
          }
        }

        # Update current model and variables
        current_model <- step_models[[best_model_name]]
        current_vars <- step_vars[[best_model_name]]

        # Save model and step information
        covariable_models[[paste("Step", step_counter)]] <- current_model
        covariable_selection_steps[[step_counter]] <- list(
          step = step_counter,
          removed_var = removed_var,
          model_comparison = model_comparison,
          best_model_name = best_model_name
        )

        step_counter <- step_counter + 1
      } else {
        if (verbose) {
          cat("  No variables can be removed without significant impact.\n")
        }
        break
      }
    }

  } else if (potential_selection_strategy == "brute_force") {
    #####################################################################
    # Brute Force Selection
    #####################################################################
    # Generate all possible combinations of potential variables
    n_potential <- length(potential_vars_clean)
    max_combinations <- 2^n_potential

    if (max_combinations > max_models) {
      warning(sprintf("The number of potential models (%d) exceeds max_models (%d). ",
                      max_combinations, max_models),
              "Using a probabilistic sampling of model space instead.")

      # Randomly sample from model space if too many combinations
      set.seed(42)  # For reproducibility
      combinations <- sample(0:(max_combinations-1), min(max_models, max_combinations))
    } else {
      combinations <- 0:(max_combinations-1)
    }

    # Initialize list for all models
    all_models <- list()
    all_models[["Core"]] <- best_core_model

    if (verbose) {
      cat(sprintf("Evaluating %d model combinations\n", length(combinations)))
      if (length(combinations) > 10) {
        pb <- txtProgressBar(min = 0, max = length(combinations), style = 3)
      }
    }

    # Evaluate all selected combinations
    for (i in seq_along(combinations)) {
      combo_idx <- combinations[i]

      # Convert index to binary representation for variable selection
      binary_rep <- integer(n_potential)
      temp <- combo_idx
      for (j in 1:n_potential) {
        binary_rep[j] <- temp %% 2
        temp <- temp %/% 2
        if (temp == 0) break
      }
      selection <- as.logical(binary_rep)

      if (sum(selection) == 0) {
        # Skip - this is just the core model
        next
      }

      # Get selected variables for this combination
      selected_vars <- potential_vars_clean[selection]
      all_vars <- c(best_core_vars, selected_vars)

      # Build the model
      model_name <- paste("Model", i)

      if (verbose) {
        description <- paste(model_name, ":", paste(selected_vars, collapse = ", "))
        all_models[[model_name]] <- build_model(all_vars, description = description)
      } else {
        all_models[[model_name]] <- build_model(all_vars)
      }

      if (verbose && length(combinations) > 10) {
        setTxtProgressBar(pb, i)
      }
    }

    if (verbose && length(combinations) > 10) {
      close(pb)
    }

    # Create comparison dataframe
    model_names <- names(all_models)
    model_comparison <- data.frame(
      Model = model_names,
      stringsAsFactors = FALSE
    )

    # Add criterion values (AIC, BIC, etc)
    for (criterion in c("AIC", "AICc", "BIC", "BICc")) {
      model_comparison[[criterion]] <- sapply(all_models, function(m) m[[criterion]])
      model_comparison[[paste0("Delta_", criterion)]] <-
        model_comparison[[criterion]] - model_comparison[[criterion]][model_comparison$Model == "Core"]
    }

    # Calculate p-values from LRT if available
    base_model <- all_models[["Core"]]
    if ("logL" %in% names(base_model)) {
      model_comparison$LRT <- NA
      base_logL <- base_model$logL
      base_df <- base_model$df

      for (i in seq_along(all_models)) {
        if (model_names[i] != "Core") {
          test_model <- all_models[[i]]
          # LRT statistic
          test_logL <- test_model$logL
          test_df <- test_model$df
          df_diff <- abs(test_df - base_df)
          if (df_diff == 0) df_diff <- 1

          # Calculate p-value
          lrt_stat <- 2 * abs(test_logL - base_logL)
          p_value <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
          model_comparison$LRT[i] <- p_value
        } else {
          model_comparison$LRT[i] <- 1.0  # Base model compared to itself
        }
      }
      # Format p-values
      model_comparison$LRT <- format.pval(model_comparison$LRT, digits = 5)
    }

    # Select best model based on criterion
    if (potential_selection_term == "p_value") {
      # For p-value, we select most complex model that's significantly better than core
      p_idx <- which(model_comparison$Model == "Core")
      p_values <- as.numeric(model_comparison$LRT[-p_idx])

      significant_models <- which(p_values < potential_selection_threshold)

      if (length(significant_models) > 0) {
        # Among significant models, choose one with best information criterion
        sig_model_indices <- setdiff(seq_along(model_comparison$Model), p_idx)[significant_models]
        aic_values <- model_comparison$AICc[sig_model_indices]
        best_idx <- sig_model_indices[which.min(aic_values)]
      } else {
        # If none are significant, stick with core model
        best_idx <- p_idx
      }
    } else {
      # For information criteria, choose model with lowest value
      term_col <- potential_selection_term
      best_idx <- which.min(model_comparison[[term_col]])
    }

    best_model_name <- model_comparison$Model[best_idx]
    final_model <- all_models[[best_model_name]]

    # Display results
    if (verbose) {
      if (best_model_name == "Core") {
        cat("\nBest model: Core model (no additional covariables selected)\n")
      } else {
        # Extract variables from the model formula
        model_formula <- final_model$Model_Formula
        added_vars_str <- setdiff(
          unlist(strsplit(model_formula, "\\s+|\\+|~")),
          c("", " ", unlist(strsplit(best_core_model$Model_Formula, "\\s+|\\+|~")))
        )
        added_vars <- added_vars_str[!grepl("offset|\\(|\\)", added_vars_str)]
        added_vars <- setdiff(added_vars, c("", " "))

        cat(sprintf("\nBest model: %s\n", best_model_name))
        if (length(added_vars) > 0) {
          cat(sprintf("Added variables: %s\n", paste(added_vars, collapse = ", ")))
        }
        cat(sprintf("%s: %.2f\n", potential_selection_term,
                    model_comparison[[potential_selection_term]][best_idx]))
      }

      # Print model comparison
      cat("\nComparison of top models:\n")

      # Sort by criterion and select top models to display
      sorted_idx <- order(model_comparison[[potential_selection_term]])
      top_n <- min(10, nrow(model_comparison))

      display_df <- data.frame(
        Model = model_comparison$Model[sorted_idx[1:top_n]],
        Criterion = round(model_comparison[[potential_selection_term]][sorted_idx[1:top_n]], 2),
        Delta = round(model_comparison[[paste0("Delta_", potential_selection_term)]][sorted_idx[1:top_n]], 2)
      )
      colnames(display_df)[2] <- potential_selection_term
      colnames(display_df)[3] <- paste0("Δ", potential_selection_term)

      print(display_df, row.names = FALSE)
      cat("\n")
    }

    # Set as selected model for output
    covariable_models <- all_models
    covariable_selection_steps <- model_comparison
  }

  #####################################################################
  # FINAL MODEL SELECTION
  #####################################################################
  # Determine the final best model
  if (potential_selection_strategy %in% c("forward", "backward")) {
    if (length(covariable_selection_steps) > 0) {
      # The last model in the step-wise process is the final selection
      final_model <- covariable_models[[paste("Step", length(covariable_selection_steps))]]
    } else {
      # No steps were taken, so the best core model is our final selection
      final_model <- best_core_model
    }
  } else {
    # For brute force, we already determined the final model
    # It's in 'final_model' from above
  }

  # Compare all models (core and covariable models)
  all_models <- c(core_models, covariable_models)

  # Make sure we don't have duplicate names
  if (any(duplicated(names(all_models)))) {
    # If there are duplicates, make names unique
    names(all_models) <- make.unique(names(all_models))
  }

  # Check if we have any models to compare
  if (length(all_models) > 0) {
    # Get names of models after making unique
    model_names <- names(all_models)

    # Use MI_models_comparison for final comparison
    final_model_comparison <- do.call(
      MI_models_comparison,
      c(all_models, list(model_names = model_names,
                         sort_by = core_selection_criteria))
    )
  } else {
    final_model_comparison <- NULL
    warning("No models to compare in final comparison.")
  }

  if (verbose) {
    cat("\n=== FINAL MODEL ===\n")
    cat(sprintf("Formula: %s\n", final_model$Model_Formula))
    cat(sprintf("%s: %.2f\n", core_selection_criteria, final_model[[core_selection_criteria]]))
  }

  return(list(
    best_model = final_model,
    core_model_selection = core_models_comparison,
    covariable_selection = covariable_selection_steps,
    model_comparison = final_model_comparison,
    final_formula = final_model$Model_Formula
  ))
}





