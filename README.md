# MIAnalysis

## Overview

MIAnalysis is an R package designed to simplify the analysis of multiply imputed datasets. It provides tools for model estimation, performance evaluation, and comparison of models across multiple imputations using Rubin's rules.

## Features

- Fit various regression models to multiply imputed data
- Automatically pool results across imputations
- Support for multiple regression types:
  - Negative Binomial
  - Linear (Gaussian)
  - Logistic (Binomial)
  - Poisson
  - Gamma
  - Quasi-Poisson
  - Quasi-Binomial
  - Cox Proportional Hazards
- Support for complex model terms:
  - Interaction terms
  - Restricted cubic splines
  - Polynomial terms
- Easy-to-use model performance assessment
- Simple model comparison framework
- Trial adjustment and follow-up offset options

## Installation

You can install the development version of MIAnalysis from GitHub with:

# Install devtools if not already installed
# install.packages("devtools")

# Install MIAnalysis from GitHub
devtools::install_github("SamLar27/MIAnalysis")

## Main Functions

### MI_estimates()

Computes estimates for statistical models using multiply imputed data.

library(MIAnalysis)

# Example usage
results <- MI_estimates(
  data = imputed_data,
  outcome_var = "outcome",
  predictor_vars = c("age", "sex", "treatment"),
  model_type = "nb"
)

### MI_model_performance()

Evaluates model performance across multiple imputations.

# Calculate performance metrics
performance <- MI_model_performance(
  data = imputed_data,
  outcome_var = "outcome",
  predictor_vars = c("age", "sex", "treatment"),
  model_type = "nb"
)

# View key metrics
performance$AIC
performance$R2
performance$C_Index

### MI_models_comparison()

Compares different models across multiple imputations.

# Create different models
model1 <- MI_model_performance(data, "outcome", c("age", "sex"))
model2 <- MI_model_performance(data, "outcome", c("age", "sex", "treatment"))
model3 <- MI_model_performance(data, "outcome", c("age", "sex", "treatment", "comorbidity"))

# Compare models
comparison <- MI_models_comparison(model1, model2, model3,
                                   model_names = c("Base", "Treatment", "Full"),
                                   sort_by = "AIC")

## Advanced Features

### Working with Splines

# Define spline terms
spline_terms <- list(
  list(var = "age", knots = 3),
  list(var = "bmi", knots = c(20, 25, 30))
)

# Fit model with splines
results <- MI_estimates(
  data = imputed_data,
  outcome_var = "outcome",
  predictor_vars = c("sex", "treatment"),
  spline_terms = spline_terms,
  model_type = "bin"
)

### Using Polynomial Terms

# Define polynomial terms
poly_terms <- list(
  list(var = "age", degree = 2),  # Quadratic
  list(var = "dose", degree = 3)  # Cubic
)

# Fit model with polynomials
results <- MI_estimates(
  data = imputed_data,
  outcome_var = "outcome",
  predictor_vars = c("sex", "treatment"),
  polynomial_terms = poly_terms,
  model_type = "lm"
)

## License

This package is licensed under the MIT License.
