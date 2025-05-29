#' Compare Multiple Imputation Models with Simplified Output
#'
#' This function compares different models created with MI_model_performance
#' and returns a simplified dataframe with key metrics. It automatically uses
#' input variable names when model_names are not specified.
#'
#' @param ... Model outputs from MI_model_performance function
#' @param model_names Character vector of model names (if NULL, uses input variable names)
#' @param reference_model Index of the reference model for LRT calculation (default: 1)
#' @param reference_name Name of the reference model (if NULL, uses the name from model_names)
#' @param sort_by Column to sort by (default: NULL for no sorting). Set to a column name (e.g., "AIC") to sort.
#' @param decreasing Logical. Should sorting be in decreasing order? Default is FALSE.
#' @param reference_first Logical. Should reference model be placed first in results? Default is TRUE.
#'
#' @return A data frame with specified model comparison metrics
#' @export
MI_models_comparison <- function(...,
                                 model_names = NULL,
                                 reference_model = 1,
                                 reference_name = NULL,
                                 sort_by = NULL,
                                 decreasing = FALSE,
                                 reference_first = TRUE) {

  # Capture all model outputs and their names
  models <- list(...)

  # Get the actual variable names used in the function call
  call_expr <- match.call()
  arg_names <- as.character(call_expr)[-1]  # Remove function name

  # Extract only the argument names that correspond to models
  model_arg_names <- arg_names[1:length(models)]

  # Remove any named arguments
  named_args <- names(call_expr)[-1]
  if (length(named_args) > 0) {
    model_arg_names <- model_arg_names[!model_arg_names %in% named_args]
  }

  # Check if any models were provided
  if (length(models) == 0) {
    stop("No models provided. Please supply at least one model output from MI_model_performance.")
  }

  # Verify each model is from MI_model_performance
  for (i in seq_along(models)) {
    if (!is.list(models[[i]]) || !all(c("Model_Type", "R2", "AIC") %in% names(models[[i]]))) {
      stop(paste("Model", i, "does not appear to be an output from MI_model_performance."))
    }
  }

  # Validate reference model
  if (reference_model < 1 || reference_model > length(models)) {
    warning("Invalid reference_model. Using first model as reference.")
    reference_model <- 1
  }

  # Set model names if not provided, using variable names if available
  if (is.null(model_names)) {
    if (length(model_arg_names) == length(models)) {
      model_names <- model_arg_names
    } else {
      model_names <- paste("Model", seq_along(models))
    }
  } else if (length(model_names) != length(models)) {
    warning("Number of model names does not match number of models. Using variable names or defaults.")
    if (length(model_arg_names) == length(models)) {
      model_names <- model_arg_names
    } else {
      model_names <- paste("Model", seq_along(models))
    }
  }

  # Set reference name if provided
  if (!is.null(reference_name)) {
    model_names[reference_model] <- reference_name
  }

  # Initialize the comparison data frame with specified columns
  comparison <- data.frame(
    Model_Name = model_names,
    Model_Type = sapply(models, function(x) x$Model_Type),
    df = sapply(models, function(x) x$df),
    logL = sapply(models, function(x) x$logL),
    LRT = NA,  # Will fill this in later
    R2 = sapply(models, function(x) x$R2),
    AIC = sapply(models, function(x) x$AIC),
    Delta_AIC = NA,  # Will fill this in later
    AICc = sapply(models, function(x) x$AICc),
    Delta_AICc = NA,  # Will fill this in later
    BIC = sapply(models, function(x) x$BIC),
    Delta_BIC = NA,  # Will fill this in later
    BICc = sapply(models, function(x) x$BICc),
    Delta_BICc = NA,  # Will fill this in later
    C_Index = sapply(models, function(x) x$C_Index),
    RMSE = sapply(models, function(x) x$RMSE),
    MAE = sapply(models, function(x) x$MAE)
  )

  # Calculate Delta AIC, Delta AICc, Delta BIC and Delta BICc
  min_aic <- min(comparison$AIC)
  comparison$Delta_AIC <- comparison$AIC - min_aic

  min_aicc <- min(comparison$AICc)
  comparison$Delta_AICc <- comparison$AICc - min_aicc

  min_bic <- min(comparison$BIC)
  comparison$Delta_BIC <- comparison$BIC - min_bic

  min_bicc <- min(comparison$BICc)
  comparison$Delta_BICc <- comparison$BICc - min_bicc


  # Calculate LRT p-values relative to the reference model
  if (length(models) > 1) {
    ref_logL <- models[[reference_model]]$logL
    ref_df <- models[[reference_model]]$df

    for (i in seq_along(models)) {
      if (i != reference_model) {
        # Calculate LRT statistic - always compute the p-value
        lrt_stat <- 2 * abs(models[[i]]$logL - ref_logL)
        df_diff <- abs(models[[i]]$df - ref_df)

        # If models have same degrees of freedom, use 1 df for chi-square test
        if (df_diff == 0) df_diff <- 1

        # Calculate p-value from chi-square distribution
        p_value <- 1 - pchisq(lrt_stat, df = df_diff)
        comparison$LRT[i] <- p_value
      } else {
        comparison$LRT[i] <- 1.0  # Reference model compared to itself - p-value is 1
      }
    }
  }

  # Add indication of reference model
  comparison$Reference <- (1:length(models) == reference_model)

  # Only sort if explicitly requested
  if (!is.null(sort_by) && sort_by %in% names(comparison)) {
    comparison <- comparison[order(comparison[[sort_by]], decreasing = decreasing), ]

    # If requested, move reference model to the top after sorting
    if (reference_first && !comparison$Reference[1]) {
      ref_row <- comparison[comparison$Reference, ]
      other_rows <- comparison[!comparison$Reference, ]
      comparison <- rbind(ref_row, other_rows)
    }
  } else if (reference_first) {
    # If no sorting but reference_first is TRUE, move reference to top
    if (!comparison$Reference[1]) {
      ref_row <- comparison[comparison$Reference, ]
      other_rows <- comparison[!comparison$Reference, ]
      comparison <- rbind(ref_row, other_rows)
    }
  }

  # Format LRT p-values
  comparison$LRT <- format.pval(comparison$LRT, digits = 3)

  # Return the data frame
  return(comparison)
}
