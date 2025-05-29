#' Test multiple predictor combinations with multiple imputation models
#'
#' @description
#' This function systematically tests different combinations of predictors in models
#' with multiply imputed data. It explores all combinations of core variables and
#' optionally tests each combination with additional variables from a secondary list.
#' It's designed for model selection and comparison across multiple imputation sets.
#'
#' @param data A data frame containing multiply imputed datasets (in long format).
#' @param outcome_var Character string specifying the outcome variable name.
#' @param core_vars_list List of character vectors, each vector containing variable names
#'   to be included as a core set in the model.
#' @param others_vars_list Optional list of character vectors, each containing additional
#'   variable names to be tested alongside each core combination. Default is NULL.
#' @param imp_col Character string specifying the column name that identifies the
#'   imputation number. Default is ".imp".
#' @param followup_offset Character string, either "Yes" or "No", indicating whether to use
#'   a follow-up variable as an offset. Default is "No".
#' @param followup_col Character string specifying the follow-up time variable name if
#'   followup_offset is "Yes". Default is NULL.
#' @param trial_factor Character string, either "Yes" or "No", indicating whether to include
#'   a trial indicator as a factor. Default is "No".
#' @param trial_col Character string specifying the trial indicator variable name if
#'   trial_factor is "Yes". Default is NULL.
#' @param imp_n Integer specifying the number of imputations to use. If NULL, all
#'   imputations are used. Default is NULL.
#' @param model_type Character string specifying the type of model to fit. Options include
#'   "nb" (negative binomial), "poisson", "linear", "logistic", etc. Default is "nb".
#' @param random_intercept Character string, either "Yes" or "No", indicating whether to
#'   include a random intercept in the model. Default is "No".
#' @param random_intercept_var Character string specifying the variable to use for the
#'   random intercept if random_intercept is "Yes". Default is NULL.
#' @param aic_method Character string specifying the method for calculating AIC in
#'   multiply imputed datasets. Default is "Rubin".
#' @param r2_method Character string specifying the method for calculating R-squared.
#'   Default is "pseudo".
#'
#' @return A data frame with rows for each tested model and columns containing:
#'   \item{Core_*}{The core variables included in the model}
#'   \item{Other_Variable}{Additional variables included beyond the core set}
#'   \item{AIC}{Akaike's Information Criterion}
#'   \item{AICc}{AIC with small sample size correction}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{BICc}{BIC with small sample size correction}
#'   \item{R2}{R-squared value (method specified by r2_method)}
#'   \item{C_Index}{Concordance index (for applicable models)}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{MAE}{Mean Absolute Error}
#'
#' @examples
#' \dontrun{
#' # Using with a multiply imputed dataset
#' results <- MI_models_optimization(
#'   data = imputed_data,
#'   outcome_var = "count_outcome",
#'   core_vars_list = list(c("age", "sex"), c("bmi", "smoking")),
#'   others_vars_list = list("comorbidity", "treatment"),
#'   model_type = "nb"
#' )
#'
#' # View the results sorted by AIC
#' results_sorted <- results[order(results$AIC), ]
#' head(results_sorted)
#' }
#'
#' @details
#' The function systematically tests multiple combinations of predictor variables in models
#' fitted to multiply imputed data. It first tests all specified combinations of core variables,
#' then tests each core combination with each additional variable from the others_vars_list.
#'
#' For each model, it calls the `MI_model_performance()` function which runs the model and
#' calculates performance metrics. Progress messages are displayed during execution.
#'
#' If a model fails to converge or encounters other errors, a warning message is displayed
#' and that combination is skipped.
#'
#' @seealso \code{\link{MI_model_performance}}
#'
#' @export
MI_models_optimization <- function(
    data,
    outcome_var,
    core_vars_list,
    others_vars_list = NULL, # New argument
    imp_col = ".imp",
    followup_offset = "No",
    followup_col = NULL,
    trial_factor = "No",
    trial_col = NULL,
    imp_n = NULL,
    model_type = "nb",
    random_intercept = "No",
    random_intercept_var = NULL,
    aic_method = "Rubin",
    r2_method = "pseudo"
) {
  results_list <- list()
  counter <- 1
  # Calculate total number of core combinations
  grid_core <- expand.grid(core_vars_list, stringsAsFactors = FALSE)
  total_core_combinations <- nrow(grid_core)
  # Give readable names
  colnames(grid_core) <- paste0("Core_", seq_along(core_vars_list))
  # Loop over all combinations of core variables
  for (i in seq_len(nrow(grid_core))) {
    current_combination <- grid_core[i, ]
    predictor_set_core <- unlist(current_combination)
    # Define sets to test: core model + (optionally) with each others_var
    sets_to_test <- list(
      predictor_set_core  # Core model alone
    )
    # If others_vars_list is provided, add one set per additional variable
    if (!is.null(others_vars_list)) {
      for (j in seq_along(others_vars_list)) {
        sets_to_test[[length(sets_to_test) + 1]] <- c(predictor_set_core, others_vars_list[[j]])
      }
    }
    # Now, for each set of predictors (core alone + core + each others_var)
    for (k in seq_along(sets_to_test)) {
      current_predictors <- sets_to_test[[k]]
      # Progression message
      if (k == 1) {
        message(glue::glue("⏳ Running core model {i}/{total_core_combinations}: {paste(current_combination, collapse = ', ')}"))
      } else {
        message(glue::glue("➕ Running core model {i}/{total_core_combinations} + other var: {others_vars_list[[k-1]]}"))
      }
      model_result <- tryCatch({
        MI_model_performance(
          data = data,
          outcome_var = outcome_var,
          predictor_vars = current_predictors,
          imp_col = imp_col,
          followup_offset = followup_offset,
          followup_col = followup_col,
          trial_factor = trial_factor,
          trial_col = trial_col,
          imp_n = imp_n,
          model_type = model_type,
          random_intercept = random_intercept,
          random_intercept_var = random_intercept_var,
          aic_method = aic_method,
          r2_method = r2_method
        )
      }, error = function(e) {
        message(glue::glue("⚠️ Model failed: {e$message}"))
        NULL
      })
      if (!is.null(model_result)) {
        results_list[[counter]] <- data.frame(
          current_combination,
          Other_Variable = ifelse(k == 1, "None", paste(others_vars_list[[k-1]], collapse = "+")),
          AIC = model_result$AIC,
          AICc = model_result$AICc,
          BIC = model_result$BIC,
          BICc = model_result$BICc,
          R2 = model_result$R2,
          C_Index = model_result$C_Index,
          RMSE = model_result$RMSE,
          MAE = model_result$MAE
        )
        counter <- counter + 1
      }
    }
  }
  final_results <- do.call(rbind, results_list)
  return(final_results)
}
