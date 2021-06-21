#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# Expected run time for this script: 2.2 hrs
# Run time with the full number of iterations: 15.5 hrs

library(AHMbook)
library(jagsUI)

# ~~~~~~~ change to RNG default in R 3.6.0 ~~~~~~~~~~~~~~~
# The values in the book were generated with R 3.5; to get the same
#   values in later versions you need the old buggy RNG for sample:
RNGversion("3.5.0")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 4.10 Analysis of citizen-science data using occupancy models
# ============================================================

# 4.10.2 Analysis of citizen-science data on Swiss middle-spotted woodpeckers
# ---------------------------------------------------------------------------

# Read in the data set from AHMbook package
data(spottedWoodpecker)
str(dat <- spottedWoodpecker)        # Look at overview of data set

# Add to data scaled date (such that 0 is 1 May and 1 unit is 1 month)
dat$date <- (dat$jdate - 121) / 30

# Check sample sizes in original data set
nsites <- length(unique(dat$site))   # 1545 sites
nyears <- length(unique(dat$year))   # 26 years
ndays <- length(unique(dat$jdate))   # 162 days in breeding season

# Randomly thin out the data set by subsampling 30%
dat.full <- dat                      # Make a copy of full data set
prop.data <- 0.3                     # Proportion of data to be used
ncase <- nrow(dat)                   # 116204
set.seed(1)                          # Ensures you get the same subset
sel.cases <- sort(sample(1:ncase, ncase * prop.data))
dat <- dat[sel.cases,]               # Smaller data set

# Look at subsampled data set
str(dat)
# 'data.frame' : 34861 obs. of 8 variables:
# $ site    : int 6 1155 1261 1262 608 741 821 907 1076 1262 ...
# $ coordx  : num 907942 1068942 1095942 1095942 1025942 ...
# $ coordy  : num 55276 169276 186276 187276 171276 ...
# $ year    : int 1990 1990 1990 1990 1992 1992 1992 1992 1992 1992 ...
# $ jdate   : int 51 51 51 51 51 51 51 51 51 51 ...
# $ y       : int 0 0 0 0 0 0 0 0 0 0 ...
# $ nsurveys: int 1 1 1 1 1 1 1 1 1 3 ...
# $ date    : num -2.33 -2.33 -2.33 -2.33 -2.33 ...

# Have to renumber the sites (since lost some in subsampling)
dat$site <- as.numeric(as.factor(dat[,"site"]))

# Sample sizes in new (subsampled) data set
nsites <- length(unique(dat$site))   # 1433 sites
nyears <- length(unique(dat$year))   # 26 years
ndays <- length(unique(dat$jdate))   # 162 days

# Plot annual total N of records and of middle spotted records (Fig. 4.24)
op <- par(mfrow = c(1,2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot(1990:2015, tapply(dat$nsurvey, list(dat$year), sum, na.rm = TRUE),
  cex = 2, type = 'b', pch = 16, ylab = 'Number of surveys', xlab = 'Year',
  frame = FALSE)
plot(1990:2015, tapply(dat$y, list(dat$year), sum, na.rm = TRUE),
  cex = 2, type = 'b', pch = 16, ylab = 'Number of middle spotted records',
  xlab = 'Year', frame = FALSE)
par(op)

# Compute number of middle spotted records per site/year (det. frequency)
table(df <- tapply(dat$y, list(dat$site, dat$year), sum, na.rm = TRUE))
#    0    1    2   3   4   5   6   7   8  10
# 8582 1063  178  40  21   6   9   8   2   3

# Proportion of missing values per site x year combo ?
(prop.NA <- sum(is.na(df)) / prod(dim(df)))
# [1] 0.7339632

# Proportion of missing value per site x year X day combo ?
(prop.NA <- 1 - (nrow(dat) / (nsites * nyears * ndays)))
# [1] 0.9942243                      # This is HUGE !!!

# Compute observed occupancy
zobs <- tapply(dat$y, list(dat$site, dat$year), max, na.rm = TRUE)
zobs[zobs>1] <- 1
psiobs <- apply(zobs, 2, mean, na.rm = TRUE)

# Bundle data (same for all models)
str(bdata <- list(y = dat[,'y'], nsurveys = dat[,'nsurvey'],
    site = dat[,'site'], year = dat[,'year']-1989, date = dat$date,
    nsites = nsites, nyears = nyears, nobs = nrow(dat)) )
# List of 8
# $ y       : int [1:34861] 0 0 0 0 0 0 0 0 0 0 ...
# $ nsurveys: int [1:34861] 1 1 1 1 1 1 1 1 1 3 ...
# $ site    : num [1:34861] 6 1092 1192 1193 573 ...
# $ year    : num [1:34861] 1 1 1 1 3 3 3 3 3 3 ...
# $ date    : num [1:34861] -2.33 -2.33 -2.33 -2.33 -2.33 ...
# $ nsites  : int 1433
# $ nyears  : int 26
# $ nobs    : int 34861

# Initial values
zst <- zobs ; zst[is.na(zst)] <- 1
inits <- function(){list(z = zst)}

# ~~~~~~ extra code for models 1 to 6 (only model 7 shown in the book) ~~~~~~
# Model 1: Static model with years as blocks (treated as fixed effects)
# --------------------------------------------------------------------

# In this model, all annual parameters (for occupancy and detection) are treated as fixed effects.

# Specify model in BUGS language for vertical data format
cat(file = "occmodel1.txt","
model {
  # Specify priors
  for (t in 1:nyears){
    psi[t] ~ dunif(0, 1)            # Occupancy
    alpha.lp[t] <- logit(mean.p[t]) # Detection parameters
    mean.p[t] ~ dunif(0, 1)
    beta.lp.1[t] ~ dnorm(0, 0.001)
    beta.lp.2[t] ~ dnorm(0, 0.001)
  }

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){
    for (t in 1:nyears){
      z[i,t] ~ dbern(psi[t])
    }
  }

  # Observation model
  for (i in 1:nobs){
    logit(p[i]) <- alpha.lp[year[i]] + beta.lp.1[year[i]] * date[i] + beta.lp.2[year[i]] * pow(date[i],2)
    y[i] ~ dbin(z[site[i],year[i]]*p[i], nsurveys[i])
  }

  # Derived parameters
  n.occ[1] <- sum(z[1:nsites,1])   # Number of occupied sites in sample
  for (t in 2:nyears){
    n.occ[t] <- sum(z[1:nsites,t])
  }
}
")

# Parameters monitored
params <- c("psi", "n.occ", "alpha.lp", "beta.lp1", "beta.lp2")

# MCMC settings
# na <- 1000  ;  ni <- 5000  ;  nb <- 1000  ;  nt <- 4  ;  nc <- 3 # 26 mins
na <- 1000  ;  ni <- 500  ;  nb <- 100  ;  nt <- 1  ;  nc <- 3 # ~~~ for testing

# Call JAGS from R, check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "occmodel1.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
print(out1, dig = 2)

# Model 2: Static model with years as blocks (treated as random effects)
# ---------------------------------------------------------------------

# This is the same as model 1, but where the annual occupancy and detection parameters are treated as random effects.

# Specify model in BUGS language for vertical data format
cat(file = "occmodel2.txt","
model {
  # Specify priors
  for (t in 1:nyears){
    logit(psi[t]) <- lpsi[t]
    lpsi[t] ~ dnorm(mu.lpsi, tau.lpsi)
    alpha.lp[t] ~ dnorm(mu.lp, tau.lp)
    beta.lp.1[t] ~ dnorm(mu.beta1.lp, tau.beta1.lp)
    beta.lp.2[t] ~ dnorm(mu.beta2.lp, tau.beta2.lp)
  }
  mu.lpsi <- logit(mean.psi)
  mean.psi ~ dunif(0, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  mu.beta1.lp ~ dnorm(0, 0.01)
  mu.beta2.lp ~ dnorm(0, 0.01)

  tau.lpsi <- pow(sd.lpsi, -2)
  sd.lpsi ~ dunif(0, 3)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 3)
  tau.beta1.lp <- pow(sd.beta1.lp, -2)
  sd.beta1.lp ~ dunif(0, 3)
  tau.beta2.lp <- pow(sd.beta2.lp, -2)
  sd.beta2.lp ~ dunif(0, 3)

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){
    for (t in 1:nyears){
      z[i,t] ~ dbern(psi[t])
    }
  }

  # Observation model
  for (i in 1:nobs){
    logit(p[i]) <- alpha.lp[year[i]] + beta.lp.1[year[i]] * date[i] + beta.lp.2[year[i]] * pow(date[i],2)
    y[i] ~ dbin(z[site[i],year[i]]*p[i], nsurveys[i])
  }

  # Derived parameters
  n.occ[1] <- sum(z[1:nsites,1])   # Number of occupied sites in sample
  for (t in 2:nyears){
    n.occ[t] <- sum(z[1:nsites,t])
  }
}
")

# Parameters monitored
params <- c("psi", "n.occ", "mean.psi", "mean.p", "mu.lpsi", "mu.lp",
    "mu.beta1.lp", "mu.beta2.lp", "sd.lpsi", "sd.lp", "sd.beta1.lp",
    "sd.beta2.lp", "alpha.lp", "beta.lp1", "beta.lp2")

# MCMC settings
# na <- 1000  ;  ni <- 5000  ;  nb <- 1000  ;  nt <- 4  ;  nc <- 3 # 27 mins
na <- 1000  ;  ni <- 500  ;  nb <- 100  ;  nt <- 1  ;  nc <- 3 # ~~~ for testing, 4 mins

# Call JAGS from R, check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "occmodel2.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out2)
print(out2, dig = 2)

# Model 3: Static model with a linear trend and random yearly deviations
# --------------------------------------------------------------------

# This model is similar to model 2, but it has a logit-linear trend in occupancy, with annual deviations around the trend which are treated as random effects.

# Specify model in BUGS language for vertical data format
cat(file = "occmodel3.txt","
model {
  # Specify priors and linear models
  # Occupancy
  alpha.lpsi <- logit(mean.psi)
  mean.psi ~ dunif(0, 1)
  beta.lpsi ~ dnorm(0, 0.01)
  for (t in 1:nyears){
    eps.lpsi[t] ~ dnorm(0, tau.lpsi)
  }
  tau.lpsi <- pow(sd.lpsi, -2)
  sd.lpsi ~ dunif(0, 3)

  # Detection
  for (t in 1:nyears){
    alpha.lp[t] ~ dnorm(mu.lp, tau.lp)
    beta.lp.1[t] ~ dnorm(mu.beta1.lp, tau.beta1.lp)
    beta.lp.2[t] ~ dnorm(mu.beta2.lp, tau.beta2.lp)
  }
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  mu.beta1.lp ~ dnorm(0, 0.01)
  mu.beta2.lp ~ dnorm(0, 0.01)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 3)
  tau.beta1.lp <- pow(sd.beta1.lp, -2)
  sd.beta1.lp ~ dunif(0, 3)
  tau.beta2.lp <- pow(sd.beta2.lp, -2)
  sd.beta2.lp ~ dunif(0, 3)

  # Ecological submodel: Define state conditional on parameters
  # Make t the outer loop
  for (t in 1:nyears){
    logit(psi[t]) <- alpha.lpsi + beta.lpsi * (t-13.5) + eps.lpsi[t]
    for (i in 1:nsites){
      z[i,t] ~ dbern(psi[t])
    }
  }

  # Observation model
  for (i in 1:nobs){
    logit(p[i]) <- alpha.lp[year[i]] + beta.lp.1[year[i]] * date[i] + beta.lp.2[year[i]] * pow(date[i],2)
    y[i] ~ dbin(z[site[i],year[i]]*p[i], nsurveys[i])
  }

  # Derived parameters: Number of occupied sites in sample and trendline
  n.occ[1] <- sum(z[1:nsites,1])
  logit(trendline[1]) <- alpha.lpsi + beta.lpsi * (1-13.5)
  for (t in 2:nyears){
    n.occ[t] <- sum(z[1:nsites,t])
    logit(trendline[t]) <- alpha.lpsi + beta.lpsi * (t-13.5)
  }
}
")

# Parameters monitored
params <- c("psi", "trendline", "n.occ", "mean.psi", "alpha.lpsi", "beta.lpsi", "sd.lpsi",
    "mean.p", "mu.lp", "mu.beta1.lp", "mu.beta2.lp", "sd.lp", "sd.beta1.lp", "sd.beta2.lp",
    "alpha.lp", "beta.lp1", "beta.lp2")

# MCMC settings
# na <- 1000  ;  ni <- 5000  ;  nb <- 1000  ;  nt <- 4  ;  nc <- 3  # 28 mins
na <- 1000  ;  ni <- 500  ;  nb <- 100  ;  nt <- 1  ;  nc <- 3  # ~~~ for testing, 4 mins

# Call JAGS from R, check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "occmodel3.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out3)
print(out3, dig = 2)

# Model 4: Static model with a quad. trend and random yearly deviations
# ---------------------------------------------------------------------

# This model has a quadratic time trend in occupancy, around which there
#   are random yearly deviations.

# Specify model in BUGS language for vertical data format
cat(file = "occmodel4.txt","
model {
  # Specify priors and linear models
  # Occupancy
  alpha.lpsi <- logit(mean.psi)
  mean.psi ~ dunif(0, 1)
  beta1.lpsi ~ dnorm(0, 0.01)
  beta2.lpsi ~ dnorm(0, 0.01)
  for (t in 1:nyears){
    eps.lpsi[t] ~ dnorm(0, tau.lpsi)
  }
  tau.lpsi <- pow(sd.lpsi, -2)
  sd.lpsi ~ dunif(0, 3)

  # Detection
  for (t in 1:nyears){
    alpha.lp[t] ~ dnorm(mu.lp, tau.lp)
    beta.lp.1[t] ~ dnorm(mu.beta1.lp, tau.beta1.lp)
    beta.lp.2[t] ~ dnorm(mu.beta2.lp, tau.beta2.lp)
  }
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  mu.beta1.lp ~ dnorm(0, 0.01)
  mu.beta2.lp ~ dnorm(0, 0.01)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 3)
  tau.beta1.lp <- pow(sd.beta1.lp, -2)
  sd.beta1.lp ~ dunif(0, 3)
  tau.beta2.lp <- pow(sd.beta2.lp, -2)
  sd.beta2.lp ~ dunif(0, 3)

  # Ecological submodel: Define state conditional on parameters
  # Make t the outer loop
  for (t in 1:nyears){
    logit(psi[t]) <- alpha.lpsi + beta1.lpsi * (t-13.5) + beta2.lpsi * pow((t-13.5), 2) + eps.lpsi[t]
    for (i in 1:nsites){
      z[i,t] ~ dbern(psi[t])
    }
  }

  # Observation model
  for (i in 1:nobs){
    logit(p[i]) <- alpha.lp[year[i]] + beta.lp.1[year[i]] * date[i] + beta.lp.2[year[i]] * pow(date[i],2)
    y[i] ~ dbin(z[site[i],year[i]]*p[i], nsurveys[i])
  }

  # Derived parameters: Number of occupied sites in sample and trendline
  n.occ[1]<-sum(z[1:nsites,1])
  logit(trendline[1]) <- alpha.lpsi + beta1.lpsi * (1-13.5) + beta2.lpsi * pow((1-13.5), 2)
  for (t in 2:nyears){
    n.occ[t] <- sum(z[1:nsites,t])
    logit(trendline[t]) <- alpha.lpsi + beta1.lpsi * (t-13.5) + beta2.lpsi * pow((t-13.5), 2)
  }
}
")

# Parameters monitored
params <- c("psi", "trendline", "n.occ", "mean.psi", "alpha.lpsi",
    "beta1.lpsi", "beta2.lpsi", "sd.lpsi",  "mean.p", "mu.lp", "mu.beta1.lp",
    "mu.beta2.lp", "sd.lp", "sd.beta1.lp", "sd.beta2.lp", "alpha.lp",
    "beta.lp1", "beta.lp2")

# MCMC settings
# na <- 1000  ;  ni <- 5000  ;  nb <- 1000  ;  nt <- 4  ;  nc <- 3  # 27 mins
na <- 1000  ;  ni <- 500  ;  nb <- 100  ;  nt <- 1  ;  nc <- 3  # ~~~ for testing, 4 mins

# Call JAGS from R, check convergence and summarize posteriors
out4 <- jags(bdata, inits, params, "occmodel4.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out4)
print(out4, dig = 2)

# Model 5: Dynamic model with fixed effects of year in phi, gamma, and p
# ----------------------------------------------------------------------

# Next, we fit one of the simplest dynamic occupancy model with fully time-dependent parameters in colonization, persistence and detection probability (with quadratic seasonal effects in the latter). All annual parameters are estimated as fixed effects.

# Specify model in BUGS language for vertical data format
cat(file = "occmodel5.txt","
model {
  # Specify priors
  psi1 ~ dunif(0, 1)            # Initial occupancy
  for (t in 1:(nyears-1)){       # For survival and persistence
    phi[t] ~ dunif(0, 1)
    gamma[t] ~ dunif(0, 1)
  }
  for (t in 1:nyears){           # For detection parameters
    alpha.lp[t] <- logit(mean.p[t])
    mean.p[t] ~ dunif(0, 1)
    beta.lp.1[t] ~ dnorm(0, 0.001)
    beta.lp.2[t] ~ dnorm(0, 0.001)
  }

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      z[i,t] ~ dbern(z[i,t-1]*phi[t-1] + (1-z[i,t-1])*gamma[t-1])
    }
  }

  # Observation model
  for (i in 1:nobs){
    logit(p[i]) <- alpha.lp[year[i]] + beta.lp.1[year[i]] * date[i] + beta.lp.2[year[i]] * pow(date[i],2)
    y[i] ~ dbin(z[site[i],year[i]]*p[i], nsurveys[i])
  }

  # Derived parameters
  psi[1] <- psi1                # Population occupancy
  n.occ[1] <- sum(z[1:nsites,1]) # Number of occupied sites in sample
  for (t in 2:nyears){
    psi[t] <- psi[t-1]*phi[t-1] + (1-psi[t-1])*gamma[t-1]
    n.occ[t] <- sum(z[1:nsites,t])
  }
}
")

# Parameters monitored
params <- c("psi", "phi", "gamma", "mean.p", "n.occ", "alpha.lp", "beta.lp.1", "beta.lp.2")

# MCMC settings
# na <- 1000  ;  ni <- 25000   ;   nb <- 5000   ;    nt <- 20   ;   nc <- 3  # 4.2 hrs
na <- 1000  ;  ni <- 2500   ;   nb <- 500   ;    nt <- 2   ;   nc <- 3  # ~~~ for testing, 20 mins

# Call JAGS from R, check convergence and summarize posteriors
out5 <- jags(bdata, inits, params, "occmodel5.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out5)
print(out5, dig = 3)

# Model 6: Dynamic model with random effects of year in phi, gamma, and p
# -----------------------------------------------------------------------

# Specify model in BUGS language for vertical data format
cat(file = "occmodel6.txt", "
model {
  # Specify priors
  psi1 ~ dunif(0, 1)            # Initial occupancy
  for (t in 1:(nyears-1)){       # For survival and persistence
    logit(phi[t]) <- lphi[t]
    lphi[t] ~ dnorm(mu.lphi, tau.lphi)
    logit(gamma[t]) <- lgamma[t]
    lgamma[t] ~ dnorm(mu.lgamma, tau.lgamma)
  }
  for (t in 1:nyears){           # For detection parameters
    alpha.lp[t] ~ dnorm(mu.alpha.lp, tau.alpha.lp)
    beta.lp.1[t] ~ dnorm(mu.beta.lp1, tau.beta.lp1)
    beta.lp.2[t] ~ dnorm(mu.beta.lp2, tau.beta.lp2)
  }

  # Hyperpriors for hyperparameters
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0, 1)
  mu.lgamma <- logit(mean.gamma)
  mean.gamma ~ dunif(0, 1)
  mu.alpha.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  mu.beta.lp1 ~ dnorm(0, 0.01)
  mu.beta.lp2 ~ dnorm(0, 0.01)
  tau.lphi <- pow(sd.lphi, -2)
  tau.lgamma <- pow(sd.lgamma, -2)
  tau.alpha.lp <- pow(sd.lp, -2)
  tau.beta.lp1 <- pow(sd.beta.lp1, -2)
  tau.beta.lp2 <- pow(sd.beta.lp2, -2)
  sd.lphi ~ dunif(0, 3)
  sd.lgamma  ~ dunif(0, 10)
  sd.lp  ~ dunif(0, 1)
  sd.beta.lp1  ~ dunif(0, 1)
  sd.beta.lp2  ~ dunif(0, 1)

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      z[i,t] ~ dbern(z[i,t-1]*phi[t-1] + (1-z[i,t-1])*gamma[t-1])
    }
  }

  # Observation model
  for (i in 1:nobs){
    logit(p[i]) <- alpha.lp[year[i]] + beta.lp.1[year[i]] * date[i] + beta.lp.2[year[i]] * pow(date[i],2)
    y[i] ~ dbin(z[site[i],year[i]]*p[i], nsurveys[i])
  }

  # Derived parameters
  psi[1] <- psi1                # Population occupancy
  n.occ[1] <- sum(z[1:nsites,1]) # Number of occupied sites in sample
  for (t in 2:nyears){
    psi[t] <- psi[t-1]*phi[t-1] + (1-psi[t-1])*gamma[t-1]
    n.occ[t] <- sum(z[1:nsites,t])
  }
}
")

# Parameters monitored
params <- c("psi", "phi", "gamma", "n.occ", "mean.phi", "mu.lphi",
    "sd.lphi", "mean.gamma", "mu.lgamma", "sd.lgamma", "mean.p", "mu.alpha.lp",
    "mu.beta.lp1", "mu.beta.lp2", "sd.lp", "sd.beta.lp1", "sd.beta.lp2",
    "alpha.lp", "beta.lp.1", "beta.lp.2")

# MCMC settings
# na <- 5000  ;  ni <- 50000  ;  nb <- 25000  ;   nt <- 25  ;  nc <- 3
na <- 5000  ;  ni <- 5000  ;  nb <- 2500  ;   nt <- 2  ;  nc <- 3  # ~~~~ for testing, 55 mins

# Call JAGS from R, check convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "occmodel6.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out6)
print(out6, dig = 2)
# ~~~~~~~~~~ end of extra code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Model 7: Same as 6, with annual heterogeneity in
# site-level extra dispersion in p
# ------------------------------------------------
# Specify model in BUGS language for vertical data format
cat(file = "occmodel7.txt", "
model {

  # Specify priors
  psi1 ~ dunif(0, 1) # Initial occupancy
  for (t in 1:(nyears-1)){ # For survival and persistence
    logit(phi[t]) <- lphi[t]
    lphi[t] ~ dnorm(mu.lphi, tau.lphi)
    logit(gamma[t]) <- lgamma[t]
    lgamma[t] ~ dnorm(mu.lgamma, tau.lgamma)
  }
  for (t in 1:nyears){ # For detection parameters
    alpha.lp[t] ~ dnorm(mu.alpha.lp, tau.alpha.lp)
    beta.lp.1[t] ~ dnorm(mu.beta.lp1, tau.beta.lp1)
    beta.lp.2[t] ~ dnorm(mu.beta.lp2, tau.beta.lp2)
  }

  # Hyperpriors for hyperparameters
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0, 1)
  mu.lgamma <- logit(mean.gamma)
  mean.gamma ~ dunif(0, 1)
  mu.alpha.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  mu.beta.lp1 ~ dnorm(0, 0.1)
  mu.beta.lp2 ~ dnorm(0, 0.1)
  tau.lphi <- pow(sd.lphi, -2)
  tau.lgamma <- pow(sd.lgamma, -2)
  tau.alpha.lp <- pow(sd.lp, -2)
  tau.beta.lp1 <- pow(sd.beta.lp1, -2)
  tau.beta.lp2 <- pow(sd.beta.lp2, -2)
  sd.lphi ~ dunif(0, 3)
  sd.lgamma ~ dunif(0, 10)
  sd.lp ~ dunif(0, 1)
  sd.beta.lp1 ~ dunif(0, 1)
  sd.beta.lp2 ~ dunif(0, 1)

  # Annually varying site random effects in detection
  for(t in 1:nyears){
    for(i in 1:nsites){
      eps.site[i,t] ~ dnorm(0, tau.p.site[year[i]])
    }
  tau.p.site[t] <- pow(sd.p.site[t],-2)
  sd.p.site[t] ~ dunif(0.001, 10)          # SD's estimated as fixed effects
  }

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      z[i,t] ~ dbern(z[i,t-1]*phi[t-1] + (1-z[i,t-1])*gamma[t-1])
    }
  }

  # Observation model
  for (i in 1:nobs){
    logit(p[i]) <- alpha.lp[year[i]] + beta.lp.1[year[i]] * date[i] +
        beta.lp.2[year[i]] * pow(date[i],2) + eps.site[site[i], year[i]]
    y[i] ~ dbin(z[site[i],year[i]]*p[i], nsurveys[i])
  }

  # Derived parameters
  psi[1] <- psi1                     # Population occupancy
  n.occ[1] <- sum(z[1:nsites,1])     # Number of occupied sites in sample
  for (t in 2:nyears){
    psi[t] <- psi[t-1]*phi[t-1] + (1-psi[t-1])*gamma[t-1]
    n.occ[t] <- sum(z[1:nsites,t])
  }
}
")

# Parameters monitored
params <- c("psi", "phi", "gamma", "n.occ", "mean.phi", "mu.lphi",
    "sd.lphi", "mean.gamma", "mu.lgamma", "sd.lgamma", "mean.p",
    "mu.alpha.lp", "mu.beta.lp1", "mu.beta.lp2", "sd.lp", "sd.beta.lp1",
    "sd.beta.lp2", "alpha.lp", "beta.lp.1", "beta.lp.2", "sd.p.site")

# MCMC settings
# na <- 5000 ; ni <- 50000 ; nb <- 25000 ; nt <- 25 ; nc <- 3
na <- 5000 ; ni <- 500 ; nb <- 250 ; nt <- 1 ; nc <- 3  # ~~~ testing, 50 mins

# Call JAGS (ART 832 min), check convergence and summarize posteriors
out7 <- jags(bdata, inits, params, "occmodel7.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out7)
print(out7, dig = 2)

# ~~~~ save output for use later ~~~~~~~~~~~~~~~~~
save(out1, out2, out3, out4, out5, out6, out7, file="AHM2_10.2_output.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~ code for Figures 4.25 - 4.27 ~~~~~~~~~~~~~~~
# Fig. 4.25
# Plot estimated occupancy trajectories for all models
plot(1990:2015, psiobs, type = "h", xlab = "Year", ylab = "Occupancy", main = '',
    ylim = c(0, 0.7), lwd = 20, lend = 'butt', frame = FALSE, col = 'grey', las = 1)
abline(h = seq(0, 0.75, 0.25), col = 'grey', lty = 2)
points(1990:2015, out1$mean$psi, pch = '1', type = 'b', col = 'blue', cex = 1.3)
points(1990:2015, out2$mean$psi, pch = '2', type = 'b', col = 'blue', cex = 1.3)
points(1990:2015, out3$mean$psi, pch = '3', type = 'b', col = 'blue', cex = 1.3)
points(1990:2015, out4$mean$psi, pch = '4', type = 'b', col = 'blue', cex = 1.3)
# lines(1990:2015, out3$mean$trendline, lty = 1, lwd = 2)
# lines(1990:2015, out4$mean$trendline, lty = 1, lwd = 2)

points(1990:2015, out5$mean$psi, pch = '5', type = 'b', col = 'red', cex = 1.3)
points(1990:2015, out6$mean$psi, pch = '6', type = 'b', col = 'red', cex = 1.3)
points(1990:2015, out7$mean$psi, pch = '7', type = 'b', col = 'red', cex = 1.3)

# Fig. 4.26
# Plot of yearly estimates of phi, gamma and p and site-level detection
#   heterogeneity
# For former three compare model 5 (fixed-effects) and 7 (random-effects)
op <- par(mfrow = c(2,2))
off <- 0.2
plot(1990:2014-off, out5$mean$gamma, xlab = "Year", ylab = "gamma", las = 1,
    frame = FALSE, pch = 16, ylim = c(0, 0.4) , col = 'blue')
segments(1990:2014-off, out5$q2.5$gamma, 1990:2014-off, out5$q97.5$gamma,
    col = 'blue')
points(1990:2014+off, out7$mean$gamma, pch = 16, col = 'red')
segments(1990:2014+off, out7$q2.5$gamma, 1990:2014+off, out7$q97.5$gamma,
    col = 'red')

off <- 0.2
plot(1990:2014-off, out5$mean$phi, xlab = "Year", ylab = "phi", las = 1,
    frame = FALSE, pch = 16, ylim = c(0, 1) , col = 'blue')
segments(1990:2014-off, out5$q2.5$phi, 1990:2014-off, out5$q97.5$phi,
    col = 'blue')
points(1990:2014+off, out7$mean$phi, pch = 16, col = 'red')
segments(1990:2014+off, out7$q2.5$phi, 1990:2014+off, out7$q97.5$phi,
    col = 'red')

off <- 0.2
plot(1990:2015-off, plogis(out5$mean$alpha.lp), xlab = "Year", ylab = "p",
    las = 1, frame = FALSE, pch = 16, ylim = c(0, 1) , col = 'blue')
segments(1990:2015-off, plogis(out5$q2.5$alpha.lp), 1990:2015-off,
    plogis(out5$q97.5$alpha.lp), col = 'blue')
points(1990:2015+off, plogis(out7$mean$alpha.lp), pch = 16,
    ylim = c(0.05, 0.2) , col = 'red')
segments(1990:2015+off, plogis(out7$q2.5$alpha.lp), 1990:2015+off,
    plogis(out7$q97.5$alpha.lp), col = 'red')
legend('topright', legend=c("Fixed effects (Model 5)",
    "Random effects (Model 7"), pch=16,
    col=c('blue', 'red'), bty='n')

plot(1990:2015, out7$mean$sd.p.site, xlab = "Year", ylab = "sd.p.site",
    las = 1, frame = FALSE, pch = 16, ylim = c(0, 10))
segments(1990:2015, out7$q2.5$sd.p.site, 1990:2015, out7$q97.5$sd.p.site)
par(op)

# Figure 4.27
# Seasonal pattern of detection probability
# Get date covariate for prediction and standardize in same way as in analysis
range(dat$jdat)
pred.date.original <- 51:212
pred.date <- (pred.date.original - 121) / 30

# Extract regression estimates
coef.alpha <- out7$mean$alpha.lp
coef.beta1 <- out7$mean$beta.lp.1
coef.beta2 <- out7$mean$beta.lp.2

# Compute predictions of p for each date and year
nyears <- length(coef.alpha)
pred.p <- array(NA, dim = c(length(pred.date), nyears))
for(k in 1:nyears){
  pred.p[,k] <- plogis(coef.alpha[k] + coef.beta1[k] * pred.date +
      coef.beta2[k] * pred.date^2)
}
# Average over all years
mean.pred <- apply(pred.p, 1, mean)

op <- par(mfrow = c(1, 2))
plot(pred.date.original, mean.pred, type = "l", lwd = 3, main = "",
    xlim = c(51, 210), ylim = c(0, 0.5), ylab = "Detection probability",
    xlab = "Julian Date", frame = FALSE)
for(t in 1:nyears){
  lines(pred.date.original, pred.p[,t], type = "l", col = "grey")
}
lines(pred.date.original, mean.pred, type = "l", lwd = 2)

# Model average annual occupancy with equal weights
str(out1$sims.list$psi)   # Look at format of data
str(allsamps <- rbind(out1$sims.list$psi, out2$sims.list$psi,
    out3$sims.list$psi, out4$sims.list$psi, out5$sims.list$psi,
    out6$sims.list$psi, out7$sims.list$psi))
pm <- apply(allsamps, 2, mean)
psd <- apply(allsamps, 2, sd)
CRI <- apply(allsamps, 2, quantile, probs=c(0.025, 0.975))

# Plot model-averaged occupancy trajectory
plot(1990:2015, psiobs, type = "h", xlab = "Year", ylab = "Occupancy",
    ylim = c(0, 0.7), lwd = 10, lend = 'butt', frame = FALSE,
    col = 'grey', las = 1)
abline(h = seq(0, 0.75, 0.25), col = 'grey', lty = 2)
points(1990:2015, pm, type = 'o', cex = 1.5, pch = 16)
segments(1990:2015, CRI[1,], 1990:2015, CRI[2,])
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
