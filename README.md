# Bayesian Survival Analysis of Breast Cancer Patients

## Overview

This project applies Bayesian regression and survival analysis to investigate associations between clinical and genomic characteristics and survival outcomes among breast cancer patients.

Using data from the MSK-IMPACT 2017 breast cancer cohort, I developed Bayesian multiple linear regression and Weibull survival models in R and JAGS. The analysis focuses not only on posterior estimation, but also on Bayesian inference, MCMC convergence, and evaluation of model fit.

## Objectives

The primary goals of this project were to:

- Identify clinical and genomic characteristics associated with differences in patient survival.
- Apply Bayesian regression methods to survival data.
- Estimate posterior distributions using Markov Chain Monte Carlo (MCMC).
- Evaluate the direction and strength of associations between predictors and survival.
- Assess MCMC convergence and model fit.
- Examine the suitability and limitations of the selected survival model.

## Data

The analysis uses clinical and genomic data from the **MSK-IMPACT 2017 breast cancer cohort**.

The dataset contains patient characteristics and genomic information that can be used to investigate variation in survival outcomes. Prior to modeling, the data were processed by:

- Selecting relevant clinical and genomic predictors.
- Encoding categorical variables for regression analysis.
- Standardizing continuous predictors.
- Constructing model matrices for Bayesian estimation.

## Methods

### Bayesian Multiple Linear Regression

An initial Bayesian multiple linear regression model was developed using log-transformed survival as the outcome.

The model was estimated using JAGS, allowing posterior distributions to be obtained for the regression coefficients and other model parameters.

### Bayesian Weibull Survival Model

A Bayesian Weibull model was subsequently used to model patient survival more directly.

The Weibull distribution provides a parametric framework for survival analysis while allowing the hazard function to change over time. Clinical and genomic predictors were incorporated into the model to estimate their associations with survival.

### MCMC Estimation

Posterior distributions were estimated using Markov Chain Monte Carlo sampling in JAGS.

Multiple chains were used to evaluate whether the sampling procedure converged to a stable posterior distribution. Posterior summaries and diagnostic plots were then examined to assess convergence and characterize uncertainty in the estimated parameters.

### Bayesian Inference

Posterior distributions were used to evaluate the direction and strength of associations between clinical and genomic predictors and patient survival.

Rather than relying solely on point estimates, the Bayesian framework allowed uncertainty in each parameter to be represented directly through its posterior distribution.

### Model Evaluation

Model performance was evaluated using posterior diagnostics and posterior predictive analysis.

The posterior predictive analysis indicated that the Weibull model tended to overestimate patient survival. This provided an opportunity to evaluate the assumptions and limitations of the model rather than relying solely on its estimated coefficients.

## Tools

- **R** — data preprocessing, analysis, and visualization
- **JAGS** — Bayesian model specification and MCMC estimation
- **MCMC** — posterior estimation and inference

## Key Takeaways

This project provided experience with the full Bayesian modeling workflow, from preprocessing clinical and genomic data through model specification, posterior estimation, inference, and model evaluation.

A particularly important result was that fitting a statistically sophisticated model does not necessarily mean that the model adequately represents the underlying data. Posterior predictive evaluation revealed systematic discrepancies between predicted and observed survival, highlighting the importance of model criticism and careful evaluation of modeling assumptions.

## Repository Structure

The repository contains the R code used for data preprocessing, Bayesian model specification, MCMC estimation, posterior analysis, and visualization.

## Limitations

The analysis was conducted as an academic statistical modeling project and should not be interpreted as a clinical prediction tool.

The Weibull model imposes a particular parametric structure on survival times and may not fully capture the complexity of the observed survival process. Results should therefore be interpreted in the context of the model assumptions and available clinical and genomic predictors.
