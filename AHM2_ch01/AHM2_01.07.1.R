#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

library(AHMbook)
library(jagsUI)

# 1.7 “Demographic” state-space models for inference about relative abundance
# ===========================================================================

# 1.7.1 Simulation assessment of a demographic state-space model
# --------------------------------------------------------------

set.seed(1)
str(tmp <- simPOP(M = 100, T = 10, mean.lam = 3, beta.lam = 0, sd.log.lam = 0,
  mean.gamma = 1.0, beta.gamma = 0, sd.log.gamma.site = 0,
  sd.log.gamma.time = 0, sd.log.gamma.survey = 0, sd.rho = 0, mean.p = 0.6,
  beta.p = 0, sd.logit.p.site = 0, sd.logit.p.time = 0, sd.logit.p.survey = 0,
  show.plot = TRUE))

# 1. Ideal case
betas <- 0.5 # Value of all covariate coefficients
str(data <- simPOP(mean.lam = 3, beta.lam = betas, mean.gamma = 1.0,
  beta.gamma = betas, sd.rho = 0, mean.p = 0.6, beta.p = betas))


# 2. Effects of “rescue”/immigration
sd.rho <- 0.2 # Value of random immigration parameter
str(data <- simPOP(mean.lam = 3, beta.lam = betas, mean.gamma = 1.0,
  beta.gamma = betas, sd.rho = sd.rho, mean.p = 0.6, beta.p = betas))

# 3. Effect of heterogeneity in lambda (in extended Markov model)
sd.log.lam <- 1 # Value of overdispersion in lambda
str(data <- simPOP(mean.lam = 3, beta.lam = betas, sd.log.lam = sd.log.lam,
  mean.gamma = 1.0, beta.gamma = betas, sd.rho = sd.rho, mean.p = 0.6,
  beta.p = betas))

# 4. Effect of heterogeneity in gamma (in extended Markov model)
sd.log.gamma.survey <- 0.5 # Value of overdispersion in gamma
str(data <- simPOP(mean.lam = 3, beta.lam = betas, mean.gamma = 1.0,
  beta.gamma = betas, sd.log.gamma.survey = sd.log.gamma.survey, sd.rho = 0.2,
  mean.p = 0.6, beta.p = betas))

# 5. Effect of heterogeneity in p (in extended Markov model)
sd.logit.p.survey <- 1 # Value of overdispersion in p
str(data <- simPOP(mean.lam = 3, beta.lam = betas, mean.gamma = 1.0,
  beta.gamma = betas, sd.rho = 0.2, mean.p = 0.6, beta.p = betas,
  sd.logit.p.survey = sd.logit.p.survey))

# 6. Effects of heterogeneity in lambda, gamma, and p simultaneously (in extended Markov model)
sd.log.lam <- 1 # Value of overdispersion in lambda
sd.log.gamma.survey <- 0.5 # Value of overdispersion in gamma
sd.logit.p.survey <- 1 # Value of overdispersion in p
str(data <- simPOP(mean.lam = 3, beta.lam = betas, sd.log.lam = sd.log.lam,
  mean.gamma = 1.0, beta.gamma = betas,
  sd.log.gamma.survey = sd.log.gamma.survey,
  sd.rho = 0.2, mean.p = 0.6, beta.p = betas,
  sd.logit.p.survey = sd.logit.p.survey))
