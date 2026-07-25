library(runjags)
library(ggplot2)
library(dplyr)

msk_data <- read.delim("~/Desktop/MATH-488S/Project/msk_impact_2017_clinical_data.tsv")

##### Data Cleaning

msk_data$Study.ID <- NULL
msk_data$Sample.ID <- NULL
msk_data$Metastatic.Site <- NULL
msk_data$Patient.ID <- NULL
msk_data$Patient.s.Vital.Status <- NULL

msk_data <- na.omit(msk_data)

breast <- msk_data %>%
  filter(Cancer.Type == "Breast Cancer")

breast$Cancer.Type <- NULL
breast$Primary.Tumor.Site <- NULL

# Remove survival status as a predictor because it is closely tied to the outcome
breast$Overall.Survival.Status <- NULL

# Keep only rows where survival time is positive before taking log
breast <- breast %>%
  filter(Overall.Survival..Months. > 0)

# Outcome variable
breast$log_survival_months <- log(breast$Overall.Survival..Months.)

# Remove original survival months from predictors
breast$Overall.Survival..Months. <- NULL


##### Set up predictors correctly

# Categorical predictors
categorical_vars <- c(
  "Cancer.Type.Detailed",
  "Matched.Status",
  "Oncotree.Code",
  "Sample.Class",
  "Sample.Collection.Source",
  "Sample.Type",
  "Sex",
  "Smoking.History",
  "Somatic.Status",
  "Specimen.Preservation.Type",
  "Specimen.Type"
)

# Only keep categorical variables that still exist
categorical_vars <- categorical_vars[categorical_vars %in% names(breast)]

# Convert to factors
for (v in categorical_vars) {
  breast[[v]] <- as.factor(breast[[v]])
}

# Drop columns with only one unique value
breast <- breast[, sapply(breast, function(x) length(unique(x)) > 1)]

# Create dummy variables
X <- model.matrix(log_survival_months ~ ., data = breast)

# Standardize numeric columns in X
# Dummy columns stay 0/1, so we only scale columns with more than 2 unique values
is_continuous <- apply(X, 2, function(col) length(unique(col)) > 2)

X_scaled <- X
X_scaled[, is_continuous] <- scale(X_scaled[, is_continuous])

# Outcome
y <- breast$log_survival_months

# Dimensions
N <- length(y)
P <- ncol(X_scaled)

##### Bayesian Multiple Linear Regression in JAGS

modelString <- "
model {

  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], invsigma2)

    mu[i] <- beta0 + inprod(beta[], X[i,])
  }

  beta0 ~ dnorm(0, 0.001)

  for (j in 1:P) {
    beta[j] ~ dnorm(0, 0.25)
  }

  invsigma2 ~ dgamma(1, 1)
  sigma <- sqrt(1 / invsigma2)
}
"

the_data <- list(
  y = y,
  X = X_scaled,
  N = N,
  P = P
)

initsfunction <- function(chain) {
  .RNG.seed <- c(1, 2)[chain]
  .RNG.name <- c("base::Super-Duper",
                 "base::Wichmann-Hill")[chain]
  
  return(list(
    .RNG.seed = .RNG.seed,
    .RNG.name = .RNG.name
  ))
}

posterior_MLR <- run.jags(
  modelString,
  n.chains = 2,
  data = the_data,
  monitor = c("beta0", "beta", "sigma"),
  adapt = 6000,
  burnin = 15000,
  sample =  15000,
  thin = 50,
  inits = initsfunction
)

summary(posterior_MLR)

plot(posterior_MLR, vars = "beta0")







# Weibull

library(runjags)
library(ggplot2)
library(dplyr)

msk_data <- read.delim("~/Desktop/MATH-488S/Project/msk_impact_2017_clinical_data.tsv")

##### Data Cleaning

msk_data$Study.ID <- NULL
msk_data$Sample.ID <- NULL
msk_data$Metastatic.Site <- NULL
msk_data$Patient.ID <- NULL
msk_data$Patient.s.Vital.Status <- NULL

msk_data <- na.omit(msk_data)

breast <- msk_data %>%
  filter(Cancer.Type == "Breast Cancer")

breast$Cancer.Type <- NULL
breast$Primary.Tumor.Site <- NULL

# Remove survival status as a predictor because it is directly related to survival
breast$Overall.Survival.Status <- NULL

# Keep only rows where survival time is positive
breast <- breast %>%
  filter(Overall.Survival..Months. > 0)

# Outcome variable for Weibull model
survival_months <- breast$Overall.Survival..Months.

# Remove outcome from predictors
breast$Overall.Survival..Months. <- NULL

##### Set up predictors correctly

categorical_vars <- c(
  "Cancer.Type.Detailed",
  "Matched.Status",
  "Oncotree.Code",
  "Sample.Class",
  "Sample.Collection.Source",
  "Sample.Type",
  "Sex",
  "Smoking.History",
  "Somatic.Status",
  "Specimen.Preservation.Type",
  "Specimen.Type"
)

# Only keep categorical variables that still exist
categorical_vars <- categorical_vars[categorical_vars %in% names(breast)]

# Convert categorical variables to factors
for (v in categorical_vars) {
  breast[[v]] <- as.factor(breast[[v]])
}

# Drop columns with only one unique value
breast <- breast[, sapply(breast, function(x) length(unique(x)) > 1)]

# Create dummy variables
X <- model.matrix(~ ., data = breast)

# Remove intercept because beta0 is already in JAGS
X <- X[, colnames(X) != "(Intercept)"]

# Standardize continuous predictors
is_continuous <- apply(X, 2, function(col) length(unique(col)) > 2)

X_scaled <- X

if (sum(is_continuous) > 0) {
  X_scaled[, is_continuous] <- scale(X_scaled[, is_continuous])
}

# Outcome
y <- survival_months

# Dimensions
N <- length(y)
P <- ncol(X_scaled)

##### Bayesian Weibull Regression in JAGS

modelString <- "
model {

  for (i in 1:N) {

    y[i] ~ dweib(shape, lambda[i])

    eta[i] <- beta0 + inprod(beta[], X[i,])

    # AFT-style Weibull model:
    # positive beta means longer expected survival time
    lambda[i] <- exp(-shape * eta[i])
  }

  beta0 ~ dnorm(0, 0.001)

  for (j in 1:P) {
    beta[j] ~ dnorm(0, 0.25)
  }

  shape ~ dgamma(1, 1)
}
"

the_data <- list(
  y = y,
  X = X_scaled,
  N = N,
  P = P
)

initsfunction <- function(chain) {
  .RNG.seed <- c(1, 2)[chain]
  .RNG.name <- c("base::Super-Duper",
                 "base::Wichmann-Hill")[chain]
  
  return(list(
    .RNG.seed = .RNG.seed,
    .RNG.name = .RNG.name,
    shape = 1
  ))
}

posterior_weibull <- run.jags(
  modelString,
  n.chains = 2,
  data = the_data,
  monitor = c("beta0", "beta", "shape"),
  adapt = 6000,
  burnin = 15000,
  sample = 15000,
  thin = 50,
  inits = initsfunction
)

summary(posterior_weibull)

plot(posterior_weibull, vars = "beta0")
