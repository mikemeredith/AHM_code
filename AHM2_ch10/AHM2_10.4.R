#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10 : INTEGRATED MODELS FOR MULTIPLE TYPES OF DATA
# =========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 4 mins

library(AHMbook)
library(jagsUI)

# 10.4 Example 2: Combination of counts plus detection/nondetection data
# ======================================================================

library(AHMbook)

# Choose sample size and parameter values for both data sets
nsites1 <- 267        # Sample size for count data
nsites2 <- 520        # Sample size for detection/nondetection data
nsurveys <- 3         # Number of surveys in both data sets
mean.lam <- 3         # Average expected abundance (lambda) per site
beta3.lam <- -1       # Coefficient of site covariate on lambda
mean.p <- 0.4         # Average per-individual detection probability (p)
beta5.p <- 1          # Effect of a site covariate on p
beta.p.survey <- -2   # Effect of observational covariate on p

# Create and summarize data set 1
set.seed(1)
str(data1 <- simNmix(nsite = nsites1, nvisit = nsurveys,
    mean.lam = mean.lam, beta3.lam = beta3.lam, mean.p = mean.p,
    beta5.p = beta5.p, beta.p.survey = beta.p.survey) )

# Create and summarize data set 2
set.seed(24)
str(data2 <- simNmix(nsite = nsites2, nvisit = nsurveys,
    mean.lam = mean.lam, beta3.lam = beta3.lam, mean.p = mean.p,
    beta5.p = beta5.p, beta.p.survey = beta.p.survey) )

# Pull out and inspect top of data set 1
head(y1 <- data1$C)

# Turn count data set 2 into detection/nondetection data
head(y2 <- data2$C)
y2[y2 > 1] <- 1               # Reduce counts >1 to 1
head(y2)

# Pull out covariates
elev1 <- data1$site.cov[,3]   # Imagine site cov 3 to be site elevation
elev2 <- data2$site.cov[,3]
hcov1 <- data1$site.cov[,5]   # Imagine site cov 5 to be habitat cover
hcov2 <- data2$site.cov[,5]
wind1 <- data1$survey.cov     # Imagine this to be wind speed
wind2 <- data2$survey.cov

# 10.4.1 Fitting the binomial N-mixture model to the count data alone
# -------------------------------------------------------------------

# Bundle and summarize data set
str(bdata <- list(y1 = data1$C, y2 = y2, nsites1 = nsites1,
    nsites2 = nsites2, nsurveys = nsurveys, elev1 = elev1, elev2 = elev2,
    hcov1 = hcov1, hcov2 = hcov2, wind1 = wind1, wind2 = wind2))
# List of 11
# $ y1      : int [1:267, 1:3] 12 1 0 1 4 1 8 2 0 0 ...
# $ y2      : num [1:520, 1:3] 1 0 0 1 1 0 0 0 1 0 ...
# $ nsites1 : num 267
# $ nsites2 : num 520
# $ nsurveys: num 3
# $ elev1   : num [1:267] -1.531 -1.827 -0.519 -0.652 -1.305 ...
# $ elev2   : num [1:520] 0.717 1.305 1.579 -0.26 -1.236 ...
# $ hcov1   : num [1:267] 1.799 -0.151 -1.117 -1.992 0.477 ...
# $ hcov2   : num [1:520] 0.783 -0.0222 -1.1118 0.1099 1.572 ...
# $ wind1   : num [1:267, 1:3] -1.8643 1.6644 1.3608 -1.2845 -0.0189 ...
# $ wind2   : num [1:520, 1:3] -0.299 -0.318 1.319 0.136 1.72 ...

# ~~~~ extra code for model 1 ~~~~~~~
# Specify binomial Nmix model in BUGS language
cat(file = "model1.txt", "
model {

  # Priors
  alpha.lam ~ dunif(-10, 10)      # Abundance intercept
  mean.lam <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 0.01)
  alpha.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)            # Detection intercept
  beta.lp1 ~ dnorm(0, 0.01)
  beta.lp2 ~ dnorm(0, 0.01)

  # Process model for data set 1
  for (i in 1:nsites1){
    N1[i] ~ dpois(lambda1[i])
    lambda1[i] <- exp(alpha.lam + beta.lam * elev1[i])
    # Observation process in data set
    for (j in 1:nsurveys){
      y1[i,j] ~ dbinom(p1[i,j], N1[i])
      logit(p1[i,j]) <- alpha.lp + beta.lp1 * hcov1[i] + beta.lp2 * wind1[i,j]
    }
  }
  # Derived quantities
  Ntotal1 <- sum(N1[])            # Total population size in sample 1
}
")

# Initial values
Nst1 <- apply(bdata$y1, 1, max)      # Avoid data/model/inits conflict
inits <- function(){list(N1 = Nst1)}

# Parameters monitored
params <- c("mean.lam", "alpha.lam", "beta.lam", "mean.p", "alpha.lp",
    "beta.lp1", "beta.lp2", "Ntotal1")

# MCMC settings
na <- 1000  ;  ni <- 10000  ;  nt <- 8  ;  nb <- 2000  ;  nc <- 3

# Call JAGS (ART 1.1 min), gauge convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "model1.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  # ~~~ replaced with 'layout' argument
traceplot(out1, layout=c(2,2))
print(out1, dig = 2)
             # mean    sd    2.5%     50%   97.5% overlap0 f Rhat n.eff
# mean.lam     3.15  0.15    2.86    3.15    3.46    FALSE 1 1.00   430
# alpha.lam    1.15  0.05    1.05    1.15    1.24    FALSE 1 1.00   430
# beta.lam    -0.99  0.03   -1.05   -0.99   -0.92    FALSE 1 1.00  3000
# mean.p       0.38  0.02    0.34    0.38    0.41    FALSE 1 1.00   977
# alpha.lp    -0.51  0.07   -0.65   -0.51   -0.37    FALSE 1 1.00   977
# beta.lp1     1.03  0.06    0.91    1.02    1.15    FALSE 1 1.02    87
# beta.lp2    -1.98  0.07   -2.11   -1.98   -1.85    FALSE 1 1.09    26
# Ntotal1   1436.64 33.88 1370.98 1436.00 1504.00    FALSE 1 1.01   382
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 10.4.2 Fitting the integrated model
# -----------------------------------

# Specify model in BUGS language
cat(file = "model2.txt", "
model {

  # Priors
  alpha.lam ~ dunif(-10, 10)        # Abundance intercept
  mean.lam <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 0.01)
  alpha.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)              # Detection intercept
  beta.lp1 ~ dnorm(0, 0.01)
  beta.lp2 ~ dnorm(0, 0.01)

  # Process model with shared parameters for data sets 1 and 2
  # Note identical alpha.lam and beta.lam for both data sets
  for (i in 1:nsites1){
    N1[i] ~ dpois(lambda1[i])
    lambda1[i] <- exp(alpha.lam + beta.lam * elev1[i])
  }
  for (i in 1:nsites2){
    N2[i] ~ dpois(lambda2[i])
    lambda2[i] <- exp(alpha.lam + beta.lam * elev2[i])
  }

  # Observation process in data set 1: binomial of an Nmix model
  for (i in 1:nsites1){
    for (j in 1:nsurveys){
      y1[i,j] ~ dbinom(p1[i,j], N1[i])
      logit(p1[i,j]) <- alpha.lp + beta.lp1 * hcov1[i] + beta.lp2 * wind1[i,j]
    }
  }
  # Observation process in data set 2: Bernoulli of a RN model
  for (i in 1:nsites2){
    for (j in 1:nsurveys){
      y2[i,j] ~ dbern(Pstar2[i,j])
      Pstar2[i,j] <- 1 - pow((1 - p2[i,j]), N2[i])
      logit(p2[i,j]) <- alpha.lp + beta.lp1 * hcov2[i] + beta.lp2 * wind2[i,j]
    }
  }

  # Derived quantities
  Ntotal1 <- sum(N1[])          # Total population size in data set 1
  Ntotal2 <- sum(N2[])          # Total population size in data set 2
}
")

# Initial values
Nst1 <- apply(bdata$y1, 1, max) # Avoid data/model/inits conflict
Nst2 <- rep(round(mean(bdata$y1)), nsites2)
inits <- function(){list(N1 = Nst1, N2 = Nst2)}

# Parameters monitored
params <- c("mean.lam", "alpha.lam", "beta.lam", "mean.p", "alpha.lp",
    "beta.lp1", "beta.lp2", "Ntotal1", "Ntotal2")

# MCMC settings
na <- 1000 ; ni <- 10000 ; nt <- 8 ; nb <- 2000 ; nc <- 3

# Call JAGS (ART 9 min), gauge convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "model2.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  # ~~~ replaced with 'layout' argument
traceplot(out2, layout=c(2,2))
print(out2, 2) # not shown


# Comparison of truth and estimates from simple model and IM
# ~~~~ extra code for the table ~~~~~~~~~~~~~
truth <- c('mean.lam' = mean.lam, 'beta.lam' = beta3.lam, 'mean.p' = mean.p,
    'beta.lp1' = beta5.p, 'beta.lp2' = beta.p.survey, 'Ntotal1' = data1$Ntotal,
    'Ntotal2' = data2$Ntotal)
esti.Nmix <- rbind(cbind(out1$summary[c(1,3,4,6:8),c(1:2)],
    '%CV' = abs(100*( out1$summary[c(1,3,4,6:8),2] /
        out1$summary[c(1,3,4,6:8),1]))), rep(NA, 3))
esti.IM <- cbind(out2$summary[c(1,3,4,6:9),c(1:2)],
    '%CV' = abs(100*( out2$summary[c(1,3,4,6:9),2] /
        out2$summary[c(1,3,4,6:9),1])))
print(cbind(truth, esti.Nmix, esti.IM), 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           truth    mean     sd %CV    mean      sd %CV
# mean.lam    3.0    3.15  0.152 4.8    3.09   0.123 4.0
# beta.lam   -1.0   -0.99  0.034 3.5   -0.99   0.028 2.9
# mean.p      0.4    0.38  0.017 4.5    0.38   0.015 4.0
# beta.lp1    1.0    1.03  0.061 6.0    0.97   0.049 5.0
# beta.lp2   -2.0   -1.98  0.073 3.7   -2.00   0.061 3.0
# Ntotal1  1395.0 1436.66 33.484 2.3 1414.35  30.137 2.1
# Ntotal2  2873.0      NA     NA  NA 2982.09 106.994 3.6
