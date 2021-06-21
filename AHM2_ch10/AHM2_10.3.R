#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10 : INTEGRATED MODELS FOR MULTIPLE TYPES OF DATA
# =========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 1 min

library(AHMbook)
library(jagsUI)

# 10.3 Example 1: Combination of “raw” and aggregated occupancy data
# ==================================================================

nsites1 <- 267        # Sample size in full data set
nsites2 <- 2000       # Sample size in aggregated data set
nsurveys <- 3         # Number of surveys (= occasions)
mean.occ <- 0.4       # Average occupancy at mean covariate value
beta1 <- -3           # Coefficient of 'elevation' on occupancy
mean.det <- 0.4       # Average per-survey detection probability
alpha2 <- -3          # Effect of 'wind speed' on detection

# Create and summarize data set 1
library(AHMbook)
set.seed(1)
str(data1 <- simOcc(M = nsites1, J = nsurveys, mean.occ = mean.occ,
    beta1 = beta1, beta2 = 0, beta3 = 0, mean.det = mean.det,
    time.effects = c(0, 0), alpha1 = 0, alpha2 = alpha2, alpha3 = 0,
    sd.lp = 0, b = 0) )

# Create and summarize data set 2
# Note: assume everything as in data set 1 except for sample size
set.seed(24)
str(data2 <- simOcc(M = nsites2, J = nsurveys, mean.occ = mean.occ,
    beta1 = beta1, beta2 = 0, beta3 = 0, mean.det = mean.det,
    time.effects = c(0, 0), alpha1 = 0, alpha2 = alpha2, alpha3 = 0,
    sd.lp = 0, b = 0) )

# Pull out and inspect top rows of data set 1
head(y1 <- data1$y)

# Aggregate detection histories in data set 2 and first elements
head(y2agg <- apply(data2$y, 1, max))

# Pull out covariates necessary in analysis
elev1 <- data1$elev
elev2 <- data2$elev
wind1 <- data1$wind

# 10.3.1 Standard occupancy model fit to “raw” occupancy data
# -----------------------------------------------------------

# Bundle and summarize data set
str(bdata <- list(y1 = y1, nsite1 = nrow(y1), nrep1 = ncol(y1), elev1 = elev1,
    wind1 = wind1))

# Specify model in BUGS language
cat(file = "model1.txt", "
model {

  # Priors
  alpha.lpsi <- logit(mean.psi)         # Occupancy intercept
  mean.psi ~ dunif(0, 1)
  beta.lpsi ~ dnorm(0, 0.01)            # Coefficient occ. covariate
  alpha.lp <- logit(mean.p)             # Detection intercept
  mean.p ~ dunif(0, 1)
  beta.lp ~ dnorm(0, 0.01)              # Coefficient det. covariate

  # Likelihood data set 1
  for (i in 1:nsite1) {
    z1[i] ~ dbern(psi1[i])
    logit(psi1[i]) <- alpha.lpsi + beta.lpsi * elev1[i]
    for (j in 1:nrep1) {
      y1[i,j] ~ dbern(z1[i] * p1[i,j])
      logit(p1[i,j]) <- alpha.lp + beta.lp * wind1[i,j]
    }
  }

  # Derived quantity data set 1
  Nocc1 <- sum(z1[]) # Number of occupied sites in Data set 1
}
")

# Initial values
zst1 <- apply(y1, 1, max) # Avoid data/model/inits conflict
inits <- function(){list(z1 = zst1)}

# Parameters monitored
params <- c("mean.psi", "alpha.lpsi", "beta.lpsi", "mean.p", "alpha.lp",
    "beta.lp", "Nocc1")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nt <- 2 ; nb <- 1000 ; nc <- 3

# Call JAGS (ART 0.3 min), gauge convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "model1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  # ~~~ replaced with 'layout' argument
traceplot(out1, layout=c(2,2))
print(out1, 2) # not shown


# 10.3.2 Fitting the integrated model
# -----------------------------------

# Bundle and summarize data set
str(bdata <- list(y1 = y1, nsite1 = nrow(y1), nrep1 = ncol(y1), elev1 = elev1,
    wind1 = wind1, y2agg = y2agg, nsite2 = length(y2agg), elev2 = elev2))
# List of 8
# $ y1    : int [1:267, 1:3] 1 0 0 1 1 0 0 0 0 1 ...
# $ nsite1: int 267
# $ nrep1 : int 3
# $ elev1 : num [1:267] -0.469 -0.256 0.146 0.816 -0.597 ...
# $ wind1 : num [1:267, 1:3] -0.766 -0.914 -0.259 -0.326 -0.653 ...
# $ y2agg : int [1:2000] 1 1 0 0 1 0 0 0 0 0 ...
# $ nsite2: int 2000
# $ elev2 : num [1:2000] -0.4149 -0.5502 0.4084 0.0378 0.3252 ...

# Specify model in BUGS language
cat(file = "model2.txt", "
model {

  # Priors
  alpha.lpsi <- logit(mean.psi)           # Occupancy intercept
  mean.psi ~ dunif(0, 1)
  beta.lpsi ~ dnorm(0, 0.01)              # Coefficient occ. covariate
  alpha.lp <- logit(mean.p)               # Detection intercept
  mean.p ~ dunif(0, 1)
  beta.lp ~ dnorm(0, 0.01)                # Coefficient det. covariate

  # Likelihood data set 1 (raw data)
  for (i in 1:nsite1) {
    z1[i] ~ dbern(psi1[i])
    # Note identical parameters in state model for data set 1 ...
    logit(psi1[i]) <- alpha.lpsi + beta.lpsi * elev1[i]
    for (j in 1:nrep1) {
      y1[i,j] ~ dbern(z1[i] * p1[i,j])
      logit(p1[i,j]) <- alpha.lp + beta.lp * wind1[i,j]
    }
  }

  # Likelihood data set 2 (aggregated data)
  for (i in 1:nsite2) {
    z2[i] ~ dbern(psi2[i])
    # ... and for data set 2.
    logit(psi2[i]) <- alpha.lpsi + beta.lpsi * elev2[i]
    y2agg[i] ~ dbern(z2[i] * Pstar2)
  }
  Pstar2 <- 1 - pow((1 - mean.p), 3)      # = Prob detected at least once

  # Derived quantities
  Nocc1 <- sum(z1[])              # Number of occupied sites (data set 1)
  Nocc2 <- sum(z2[])              # Number of occupied sites (data set 2)
}
")

# Initial values
zst1 <- apply(y1, 1, max)         # Avoid data/model/inits conflict
zst2 <- rep(1, length(y2agg))
inits <- function(){list(z1 = zst1, z2 = zst2)}

# Parameters monitored
params <- c("mean.psi", "alpha.lpsi", "beta.lpsi", "mean.p", "alpha.lp",
    "beta.lp", "Pstar2", "Nocc1", "Nocc2")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nt <- 2 ; nb <- 1000 ; nc <- 3

# Call JAGS (ART 1.6 min), gauge convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "model2.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  # ~~~ replaced with 'layout' argument
traceplot(out2, layout=c(2,2))
print(out2, 2) # not shown

# ~~~  code to produce the first table ~~~
nocc.table <- rbind(c('truth' = data1$sumZ,
    'observed' = data1$sumZ.obs, out2$summary[8,c(1,3,7)]),
    c('true data set 2' = data2$sumZ, 'observed' = data2$sumZ.obs,
    out2$summary[9,c(1,3,7)]))
rownames(nocc.table) <- c('Data set 1', 'Data set 2')
nocc.table
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            truth observed    mean 2.5%   97.5%
# Data set 1   118       99 115.061  108 123.000
# Data set 2   892      732 869.712  819 923.025

# Comparison of truth and estimates from simple model and IM
# ~~~ code for second table ~~~~~~~~~~~~~~~~~~~~
truth <- c('Occupancy intercept' = mean.occ, 'Occupancy slope' = beta1,
    'Detection intercept' = mean.det, 'Detection slope' = alpha2)
esti.simple <- out1$summary[c(1,3,4,6), c(1,3,7)]
esti.IM <- out2$summary[c(1,3,4,6), c(1,3,7)]
print(cbind(truth, esti.simple, esti.IM), 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     truth   mean   2.5%  97.5%   mean   2.5%  97.5%
# Occupancy intercept   0.4  0.415  0.335  0.499  0.394  0.352  0.438
# Occupancy slope      -3.0 -2.617 -3.425 -1.892 -2.966 -3.340 -2.639
# Detection intercept   0.4  0.450  0.376  0.526  0.460  0.406  0.524
# Detection slope      -3.0 -2.885 -3.513 -2.281 -2.863 -3.491 -2.281
