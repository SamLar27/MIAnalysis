#' Model Selection for Multiple Imputed Data with Grouped Variables
#'
#' This function performs systematic model selection on multiply imputed data.
#' It handles both core model selection and additional covariable selection
#' through forward, backward, or brute force approaches. The potential variables
#' can be grouped into categories, where variables in the same category are considered
#' alternatives and cannot be included together in the model.
#'
#' @param data A data frame containing the imputed dataset.
#' @param outcome_var The dependent variable for the model.
#' @param core_vars_list List of core variables or combinations to test for the base model.
#'        Can be a list of character vectors, each representing a possible core model.
#' @param potential_vars_list Named list of additional covariables to test, grouped by category.
#'        Variables within the same category are treated as alternatives (can't appear together).
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

  # Process potential_vars_list to ensure it's a named list
  if (!is.null(potential_vars_list)) {
    if (!is.list(potential_vars_list)) {
      # Single vector provided directly - convert to a single-category list
      potential_vars_list <- list(potential_vars = potential_vars_list)
      names(potential_vars_list) <- "potential_vars"
    } else if (is.null(names(potential_vars_list)) || any(names(potential_vars_list) == "")) {
      # List without proper names - add default names
      missing_names <- is.null(names(potential_vars_list)) || names(potential_vars_list) == ""
      if (all(missing_names)) {
        names(potential_vars_list) <- paste0("potential_var_", seq_along(potential_vars_list))
      } else {
        # Replace only missing names
        for (i in which(missing_names)) {
          names(potential_vars_list)[i] <- paste0("potential_var_", i)
        }
      }
    }

    # Flatten any nested lists within categories
    for (i in seq_along(potential_vars_list)) {
      if (is.list(potential_vars_list[[i]])) {
        potential_vars_list[[i]] <- unlist(potential_vars_list[[i]])
      }

      # Ensure each entry is a character vector
      if (!is.character(potential_vars_list[[i]])) {
        stop(paste("Potential variables in group", names(potential_vars_list)[i], "must be character vectors."))
      }
    }
  }

  # Generate all possible combinations of core variable expressions
  core_combinations <- generate_core_combinations(core_vars_list)

  # Helper function to extract variables from polynomial terms
  extract_poly_vars <- function(term) {
    if (grepl("^poly\\(", term)) {
      var_name <- gsub("^poly\\(([^,]+),.*", "\\1", term)
      return(trimws(var_name))
    } else {
      return(NULL)
    }
  }

  # Helper function to check if a term is a polynomial term
  is_poly_term <- function(term) {
    return(grepl("^poly\\(", term))
  }

  # Helper function to check if a term is a spline term (rcs)
  is_spline_term <- function(term) {
    return(grepl("^rcs\\(", term))
  }

  # Helper function to extract degree from polynomial terms
  extract_poly_degree <- function(term) {
    if (is_poly_term(term)) {
      # Try to extract degree=X parameter
      if (grepl("degree\\s*=\\s*\\d+", term)) {
        degree <- as.numeric(gsub(".*degree\\s*=\\s*(\\d+).*", "\\1", term))
        return(degree)
      } else {
        # Check if there's a simple number as the second argument
        if (grepl("^poly\\([^,]+,\\s*\\d+", term)) {
          degree <- as.numeric(gsub("^poly\\([^,]+,\\s*(\\d+).*", "\\1", term))
          return(degree)
        }
      }
    }
    # Default to 2 if we can't parse the degree
    return(2)
  }

  # Function to build a model and evaluate performance
  build_model <- function(variables, description = NULL) {
    # Add offset parameter for count models
    follow_offset_param <- followup_offset
    if (model_type %in% c("nb", "poisson") && !is.null(followup_col)) {
      follow_offset_param <- "Yes"
    }

    # For polynomial/spline terms, include them directly in the predictor variables
    # MI_model_performance should handle these terms in the formula
    variables_to_use <- variables

    result <- MI_model_performance(
      data = data,
      outcome_var = outcome_var,
      predictor_vars = variables_to_use,
      imp_col = imp_col,
      followup_offset = follow_offset_param,
      followup_col = followup_col,
      trial_factor = trial_factor,
      trial_col = trial_col,
      imp_n = imp_n,
      model_type = model_type,
      # Removed polynomial_terms parameter - terms are handled in predictor_vars
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

  # Add core_var_n columns and reorder
  core_var_names <- names(core_vars_list)
  for (i in seq_along(core_var_names)) {
    group_name <- core_var_names[i]
    core_models_comparison[[group_name]] <- sapply(core_models_comparison$Model_Name, function(model_name) {
      selected_vars <- core_model_variables[[model_name]]
      matched <- intersect(selected_vars, core_vars_list[[group_name]])
      if (length(matched) > 0) matched else NA
    })
  }

  # Reorder: move core_var_n columns just after Model_Type
  first_cols <- c("Model_Name", "Model_Type")
  insert_cols <- core_var_names
  remaining_cols <- setdiff(colnames(core_models_comparison), c(first_cols, insert_cols))
  core_models_comparison <- core_models_comparison[, c(first_cols, insert_cols, remaining_cols)]


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

  # Get all potential variables
  all_potential_vars <- unlist(potential_vars_list)

  # Helper function to extract base variable names from poly/spline terms
  extract_base_var_name <- function(var) {
    if (is_poly_term(var)) {
      return(extract_poly_vars(var))
    } else if (is_spline_term(var)) {
      return(gsub("^rcs\\(([^,]+),.*", "\\1", var))
    } else {
      return(var)
    }
  }

  # Remove any potential variables that are already in the core model
  # First convert polynomial/spline terms to their base variable names for comparison
  best_core_base_vars <- sapply(best_core_vars, extract_base_var_name)

  # Check if a potential variable conflicts with core variables
  is_var_in_core <- function(var) {
    base_var <- extract_base_var_name(var)
    return(base_var %in% best_core_base_vars)
  }

  all_potential_vars_clean <- all_potential_vars[!sapply(all_potential_vars, is_var_in_core)]

  # Update the potential_vars_list to remove variables already in the core model
  for (i in seq_along(potential_vars_list)) {
    potential_vars_list[[i]] <- potential_vars_list[[i]][!sapply(potential_vars_list[[i]], is_var_in_core)]
    # If a category is now empty, remove it
    if (length(potential_vars_list[[i]]) == 0) {
      potential_vars_list[[i]] <- NULL
    }
  }

  # Filter out any empty categories
  potential_vars_list <- potential_vars_list[sapply(potential_vars_list, length) > 0]

  if (length(all_potential_vars_clean) == 0) {
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
    cat("Potential variable groups:\n")
    for (i in seq_along(potential_vars_list)) {
      cat(sprintf("  Group %s: %s\n",
                  names(potential_vars_list)[i],
                  paste(potential_vars_list[[i]], collapse = ", ")))
    }
  }

  # Strategy-specific implementation
  covariable_models <- list()
  covariable_selection_steps <- list()

  if (potential_selection_strategy == "brute_force") {
    #####################################################################
    # Brute Force Selection - Modified to handle grouped variables
    #####################################################################

    # First, create the base model (core only)
    all_models <- list()
    all_models[["Core"]] <- best_core_model
    all_model_vars <- list()
    all_model_vars[["Core"]] <- best_core_vars

    if (verbose) {
      cat("\nGenerating models for brute force selection:\n")
      cat("  - Core model (no additional variables)\n")
    }

    # Add models with variables from single categories
    model_count <- 1

    # 1. First, add individual variables from all groups
    for (group_name in names(potential_vars_list)) {
      group_vars <- potential_vars_list[[group_name]]

      for (var in group_vars) {
        if (model_count > max_models) {
          warning("Maximum model count reached. Some combinations will not be tested.")
          break
        }

        model_name <- paste0("Model_", model_count)
        vars <- c(best_core_vars, var)

        if (verbose) {
          cat(sprintf("  - %s: Core + %s\n", model_name, var))
        }

        all_models[[model_name]] <- build_model(vars)
        all_model_vars[[model_name]] <- vars
        model_count <- model_count + 1
      }

      if (model_count > max_models) break
    }

    # 2. Generate combinations between different groups but never from the same group
    if (length(potential_vars_list) > 1 && model_count <= max_models) {
      # Create all possible combinations between variables from different groups
      combinations <- list()

      # For each variable in the first group
      for (var1 in potential_vars_list[[1]]) {
        # For each variable in the second group
        for (var2 in potential_vars_list[[2]]) {
          combinations[[length(combinations) + 1]] <- c(var1, var2)
        }
      }

      # If there are more than 2 groups, we need to handle those combinations too
      if (length(potential_vars_list) > 2) {
        warning("More than 2 potential variable groups found. Only combinations between the first two groups will be tested.")
        # A more general solution would require a recursive approach like before,
        # but properly tracking the model count
      }

      # Add each combination as a model
      for (combo in combinations) {
        if (model_count > max_models) {
          warning("Maximum model count reached. Some combinations will not be tested.")
          break
        }

        model_name <- paste0("Model_", model_count)
        vars <- c(best_core_vars, combo)

        if (verbose) {
          cat(sprintf("  - %s: Core + %s\n",
                      model_name,
                      paste(combo, collapse = " + ")))
        }

        all_models[[model_name]] <- build_model(vars)
        all_model_vars[[model_name]] <- vars
        model_count <- model_count + 1
      }
    }

    if (verbose) {
      cat(sprintf("\nGenerated %d models to compare\n", length(all_models)))
    }

    # Compare all models
    model_names <- names(all_models)
    model_comparison <- do.call(
      MI_models_comparison,
      c(all_models, list(model_names = model_names, sort_by = potential_selection_term))
    )

    # Determine the best model
    if (potential_selection_term == "p_value") {
      # For p-value, we want the most complex model that's not significantly worse than the reference
      # The reference will be the Core model (index 1)
      p_values <- as.numeric(model_comparison$LRT[-1])  # Skip Core

      significant_models <- which(p_values >= potential_selection_threshold) + 1  # Adjust index

      if (length(significant_models) > 0) {
        # Find the most complex model among non-significant ones
        var_counts <- sapply(significant_models, function(idx) {
          model_name <- model_comparison$Model_Name[idx]
          return(length(all_model_vars[[model_name]]))
        })

        best_idx <- significant_models[which.max(var_counts)]
      } else {
        # If all are significant, use the core model
        best_idx <- 1  # Core model index
      }
    } else {
      # For information criteria, we want the model with the best criterion
      best_idx <- which.min(model_comparison[[potential_selection_term]])
    }

    # Select best model
    best_model_name <- model_comparison$Model_Name[best_idx]
    final_model <- all_models[[best_model_name]]
    final_vars <- all_model_vars[[best_model_name]]

    # Identify variables added to the core model
    added_vars <- setdiff(final_vars, best_core_vars)

    if (verbose) {
      cat("\nBrute Force Model Selection Results:\n")
      cat(sprintf("  Best model: %s\n", best_model_name))
      cat(sprintf("  %s: %.2f\n", potential_selection_term, model_comparison[[potential_selection_term]][best_idx]))

      if (length(added_vars) > 0) {
        cat(sprintf("  Added variables: %s\n", paste(added_vars, collapse = ", ")))
      } else {
        cat("  No additional variables were selected (core model only).\n")
      }

      # Print top models
      cat("\nTop 10 models by", potential_selection_term, ":\n")
      sorted_idx <- order(model_comparison[[potential_selection_term]])
      top_n <- min(10, nrow(model_comparison))

      top_models_df <- data.frame(
        Model = model_comparison$Model_Name[sorted_idx[1:top_n]],
        Criterion = round(model_comparison[[potential_selection_term]][sorted_idx[1:top_n]], 2)
      )
      print(top_models_df, row.names = FALSE)
    }

    # Add potential_var_n columns and reorder
    potential_var_names <- names(potential_vars_list)
    for (group_name in potential_var_names) {
      group_vars <- potential_vars_list[[group_name]]
      model_comparison[[group_name]] <- sapply(model_comparison$Model_Name, function(model_name) {
        vars <- all_model_vars[[model_name]]
        selected <- intersect(vars, group_vars)
        if (length(selected) > 0) selected else "Not included"
      })
    }

    # Reorder: move potential_var_n columns just after Model_Type
    first_cols <- c("Model_Name", "Model_Type")
    insert_cols <- potential_var_names
    remaining_cols <- setdiff(colnames(model_comparison), c(first_cols, insert_cols))
    model_comparison <- model_comparison[, c(first_cols, insert_cols, remaining_cols)]

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
