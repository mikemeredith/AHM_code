#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6 : MULTISTATE OCCUPANCY MODELS
# =======================================
# Code from proofs dated 2020-08-19

# Approximate code execution time for this script: 2 hrs
# With full number of iterations: 19 hrs

library(jagsUI)

# ~~~~~~~~~~~ need data preparation from 6.4.1 ~~~~~~~~~~~~~~~~~~
source("AHM2_06.4.1.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 6.4 Case study: Swiss eagle owls
# ================================

# 6.4.3 Fitting a big dynamic multiseason model with covariates and
#       random effects
# ------------------------------------------------------------------

# This makes for much easier interpretation than use of scale()
date.scaled <- (date - 90) / 30.5
summary(date.scaled)

# Zero-impute missing dates
date.scaled[is.na(date.scaled)] <- 0
sum(is.na(date.scaled))

# Create region factor
region <- as.numeric(dat$sites$region)

# Create elevation and forest, both scaled
summary(elev <- dat$sites$elev)
summary(forest <- dat$sites$forest)
elev.scaled <- standardize(elev)
forest.scaled <- standardize(forest)

# Data bundle
str(bdata <- list(y = yms, nsites = dim(yms)[1], nsurveys = nsurveys,
    nyears = dim(yms)[3], region = region, elev = elev.scaled,
    forest = forest.scaled, date = date.scaled))
# List of 8
# $ y       : num [1:274, 1:20, 1:10] NA NA 1 NA NA NA NA 2 2 2 ...
# $ nsites  : int 274
# $ nsurveys: num [1:274, 1:10] 1 1 3 1 1 1 1 20 12 20 ...
# $ nyears  : int 10
# $ region  : num [1:274] 2 5 6 1 1 2 3 5 3 2 ...
# $ elev    : num [1:274] -0.84 0.115 -0.554 -0.267 1.07 ...
# $ forest  : num [1:274] -0.647 0.862 -0.459 -0.459 -1.213 ...
# $ date    : num [1:274, 1:20, 1:10] 0 0 -1.9 0 0 ...

# Specify model in BUGS language
cat(file = "dynMS2.txt", "
model {

  ### (1) Linear models and priors
  # ------------------------------
  ## (a) State process
  # Initial state vector Omega
  # Linear models
  for (i in 1:nsites){
    logit(psi[i]) <- alpha.lpsi[region[i]] + beta.lpsi[1] * elev[i] +
        beta.lpsi[2] * forest[i]
    logit(r[i]) <- alpha.lr[region[i]] + beta.lr[1] * elev[i] +
        beta.lr[2] * forest[i]
  }

  # Priors for parameters in the linear models of psi and r (Omega)
  # Region-specific intercepts (fixed effects)
  for(k in 1:6){
    alpha.lpsi[k] <- logit(mean.psi[k])
    mean.psi[k] ~ dunif(0, 1)
    alpha.lr[k] <- logit(mean.r[k])
    mean.r[k] ~ dunif(0, 1)
  }
  # Coefficients of 2 covariates in parameters in Omega
  for(k in 1:2){
    beta.lpsi[k] ~ dnorm(0, 0.1)
    beta.lr[k] ~ dnorm(0, 0.1)
  }

  # Transition probability matrix (PhiMat)
  # Linear models
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      logit(phi1[i,t]) <- alpha.lphi1[t] + beta.region.lphi1[region[i]] +
          beta.lphi1[1] * elev[i] + beta.lphi1[2] * forest[i]
      logit(rho1[i,t]) <- alpha.lrho1[t] + beta.region.lrho1[region[i]] +
          beta.lrho1[1] * elev[i] + beta.lrho1[2] * forest[i]
      logit(phi2[i,t]) <- alpha.lphi2[t] + beta.region.lphi2[region[i]] +
          beta.lphi2[1] * elev[i] + beta.lphi2[2] * forest[i]
      logit(rho2[i,t]) <- alpha.lrho2[t] + beta.region.lrho2[region[i]] +
          beta.lrho2[1] * elev[i] + beta.lrho2[2] * forest[i]
      logit(phi3[i,t]) <- alpha.lphi3[t] + beta.region.lphi3[region[i]] +
          beta.lphi3[1] * elev[i] + beta.lphi3[2] * forest[i]
      logit(rho3[i,t]) <- alpha.lrho3[t] + beta.region.lrho3[region[i]] +
          beta.lrho3[1] * elev[i] + beta.lrho3[2] * forest[i]
    }
  }
  # Priors for parameters in the linear models in PhiMat
  # Year-specific intercepts (random effects)
  for(t in 1:(nyears-1)){
    alpha.lphi1[t] ~ dnorm(mu.alpha.lphi1, tau.alpha.lphi1)
    alpha.lrho1[t] ~ dnorm(mu.alpha.lrho1, tau.alpha.lrho1)
    alpha.lphi2[t] ~ dnorm(mu.alpha.lphi2, tau.alpha.lphi2)
    alpha.lrho2[t] ~ dnorm(mu.alpha.lrho2, tau.alpha.lrho2)
    alpha.lphi3[t] ~ dnorm(mu.alpha.lphi3, tau.alpha.lphi3)
    alpha.lrho3[t] ~ dnorm(mu.alpha.lrho3, tau.alpha.lrho3)
  }

  # Hyperpriors for hyperparameters governing the random year effects
  mu.alpha.lphi1 <- logit(mean.phi1)
  mean.phi1 ~ dunif(0, 1)
  tau.alpha.lphi1 <- pow(sd.alpha.lphi1, -2)
  sd.alpha.lphi1 ~ dnorm(0, 0.5)I(0,)
  mu.alpha.lrho1 <- logit(mean.rho1)
  mean.rho1 ~ dunif(0, 1)
  tau.alpha.lrho1 <- pow(sd.alpha.lrho1, -2)
  sd.alpha.lrho1 ~ dnorm(0, 0.5)I(0,)
  # -------------------------------
  mu.alpha.lphi2 <- logit(mean.phi2)
  mean.phi2 ~ dunif(0, 1)
  tau.alpha.lphi2 <- pow(sd.alpha.lphi2, -2)
  sd.alpha.lphi2 ~ dnorm(0, 0.5)I(0,)
  mu.alpha.lrho2 <- logit(mean.rho2)
  mean.rho2 ~ dunif(0, 1)
  tau.alpha.lrho2 <- pow(sd.alpha.lrho2, -2)
  sd.alpha.lrho2 ~ dnorm(0, 0.5)I(0,)
  # -------------------------------
  mu.alpha.lphi3 <- logit(mean.phi3)
  mean.phi3 ~ dunif(0, 1)
  tau.alpha.lphi3 <- pow(sd.alpha.lphi3, -2)
  sd.alpha.lphi3 ~ dnorm(0, 0.5)I(0,)
  mu.alpha.lrho3 <- logit(mean.rho3)
  mean.rho3 ~ dunif(0, 1)
  tau.alpha.lrho3 <- pow(sd.alpha.lrho3, -2)
  sd.alpha.lrho3 ~ dnorm(0, 0.5)I(0,)

  # Fixed effects of region on the parameters in PhiMat
  beta.region.lphi1[1] <- 0 # Avoid overparameterization
  beta.region.lrho1[1] <- 0
  beta.region.lphi2[1] <- 0
  beta.region.lrho2[1] <- 0
  beta.region.lphi3[1] <- 0
  beta.region.lrho3[1] <- 0
  for(k in 2:6){
    beta.region.lphi1[k] ~ dnorm(0, 0.1)
    beta.region.lphi2[k] ~ dnorm(0, 0.1)
    beta.region.lphi3[k] ~ dnorm(0, 0.1)
    beta.region.lrho1[k] ~ dnorm(0, 0.1)
    beta.region.lrho2[k] ~ dnorm(0, 0.1)
    beta.region.lrho3[k] ~ dnorm(0, 0.1)
  }
  # Coefficients of 2 covariates in parameters in PhiMat
  for(k in 1:2){
    beta.lphi1[k] ~ dnorm(0, 0.1)
    beta.lphi2[k] ~ dnorm(0, 0.1)
    beta.lphi3[k] ~ dnorm(0, 0.1)
    beta.lrho1[k] ~ dnorm(0, 0.1)
    beta.lrho2[k] ~ dnorm(0, 0.1)
    beta.lrho3[k] ~ dnorm(0, 0.1)
  }

  # (b) Observation process
  # Observation matrix (Theta)
  # Linear models
  for (i in 1:nsites){
    for(t in 1:nyears){
      for(j in 1:nsurveys[i,t]){ # nsurveys is function of site and year
        # Observation model for sites in occupied state 1 (= single bird)
        logit(p2[i,j,t]) <- alpha.lp2 + beta.region.lp2[region[i]] +
            beta.lp2[1] * date[i,j,t] + beta.lp2[2] * pow(date[i,j,t],2)
        # Observation model for sites in occupied state 2 (= pairs)
        mlogit.p3[2,i,j,t] <- alpha.lp32 + beta.region.lp32[region[i]] +
            beta.lp32[1] * date[i,j,t] + beta.lp32[2] * pow(date[i,j,t],2)
        mlogit.p3[3,i,j,t] <- alpha.lp33 + beta.region.lp33[region[i]] +
            beta.lp33[1] * date[i,j,t] + beta.lp33[2] * pow(date[i,j,t],2)
      }
    }
  }

  # Priors for parameters in the linear models in Theta
  # Intercepts
  alpha.lp2 <- logit(mean.alpha.p2)
  mean.alpha.p2 ~ dunif(0, 1)
  alpha.lp32 ~ dnorm(0, 0.01)        # Must be normal for multinomial logit
  alpha.lp33 ~ dnorm(0, 0.01)
  # Fixed effects of region in parameters in Theta
  beta.region.lp2[1] <- 0            # Avoid overparameterization
  beta.region.lp32[1] <- 0
  beta.region.lp33[1] <- 0
  for(k in 2:6){
    beta.region.lp2[k] ~ dnorm(0, 0.1)
    beta.region.lp32[k] ~ dnorm(0, 0.1)
    beta.region.lp33[k] ~ dnorm(0, 0.1)
  }
  # Coefficients of survey date (linear and squared) in params in Theta
  for(k in 1:2){
    beta.lp2[k] ~ dnorm(0, 0.1)
    beta.lp32[k] ~ dnorm(0, 0.1)
    beta.lp33[k] ~ dnorm(0, 0.1)
  }

  # Implement Multinomial logit link for p3[2:3] in Theta
  for (i in 1:nsites){
    for (t in 1:nyears){
      for(j in 1:nsurveys[i,t]){
        p3[2,i,j,t] <- exp(mlogit.p3[2,i,j,t]) / (1 + exp(mlogit.p3[2,i,j,t]) +
            exp(mlogit.p3[3,i,j,t]))
        p3[3,i,j,t] <- exp(mlogit.p3[3,i,j,t]) / (1 + exp(mlogit.p3[2,i,j,t]) +
            exp(mlogit.p3[3,i,j,t]))
      }
    }
  }

  ### (2) Define Initial state vector (Omega), state transition matrix
  #(PhiMat) and observation matrix (Theta)
  # Initial state vector (Omega): Year 1
  for (i in 1:nsites){
    Omega[i,1] <- 1 - psi[i] # Prob. of non-occupation
    Omega[i,2] <- psi[i] * (1-r[i]) # Prob. of occ. by single bird
    Omega[i,3] <- psi[i] * r[i] # Prob. of occ. by pair
  }

  # State transition probability matrix (PhiMat): years 2:nyears
  # Note conditional Bernoulli parameterization of multinomial
  # Order of indices: Departing state, arrival state, site, year
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      PhiMat[1,1,i,t] <- 1 - phi1[i,t]
      PhiMat[1,2,i,t] <- phi1[i,t] * (1 - rho1[i,t])
      PhiMat[1,3,i,t] <- phi1[i,t] * rho1[i,t]
      PhiMat[2,1,i,t] <- 1 - phi2[i,t]
      PhiMat[2,2,i,t] <- phi2[i,t] * (1 - rho2[i,t])
      PhiMat[2,3,i,t] <- phi2[i,t] * rho2[i,t]
      PhiMat[3,1,i,t] <- 1 - phi3[i,t]
      PhiMat[3,2,i,t] <- phi3[i,t] * (1 - rho3[i,t])
      PhiMat[3,3,i,t] <- phi3[i,t] * rho3[i,t]
    }
  }

  # Observation probability matrix (Theta): years 1:nyears
  # Order of indices: true state, observed state, site, occasion, year
  # (No conditional Bernoulli reparameterization here)
  for(i in 1:nsites){
    for(t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        Theta[1,1,i,j,t] <- 1
        Theta[1,2,i,j,t] <- 0
        Theta[1,3,i,j,t] <- 0
        Theta[2,1,i,j,t] <- 1-p2[i,j,t]
        Theta[2,2,i,j,t] <- p2[i,j,t]
        Theta[2,3,i,j,t] <- 0
        Theta[3,1,i,j,t] <- 1-p3[2,i,j,t]-p3[3,i,j,t]
        Theta[3,2,i,j,t] <- p3[2,i,j,t]
        Theta[3,3,i,j,t] <- p3[3,i,j,t]
      }
    }
  }

  ### (3) Likelihood
  # Initial state: year 1
  for (i in 1:nsites){
    z[i,1] ~ dcat(Omega[i,])
  }

  # State transitions from yearly interval 1:(nyears-1)
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      z[i,t+1] ~ dcat(PhiMat[z[i,t],,i,t])
    }
  }

  # Observation equation
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        y[i,j,t] ~ dcat(Theta[z[i, t],,i,j,t])
      }
    }
  }

  ### (4) Derived quantities
  # Annual average PhiMat (averaged over sites)
  for(t in 1:(nyears-1)){
    PhiMat.annual[1,1,t] <- mean(PhiMat[1,1,,t])
    PhiMat.annual[1,2,t] <- mean(PhiMat[1,2,,t])
    PhiMat.annual[1,3,t] <- mean(PhiMat[1,3,,t])
    PhiMat.annual[2,1,t] <- mean(PhiMat[2,1,,t])
    PhiMat.annual[2,2,t] <- mean(PhiMat[2,2,,t])
    PhiMat.annual[2,3,t] <- mean(PhiMat[2,3,,t])
    PhiMat.annual[3,1,t] <- mean(PhiMat[3,1,,t])
    PhiMat.annual[3,2,t] <- mean(PhiMat[3,2,,t])
    PhiMat.annual[3,3,t] <- mean(PhiMat[3,3,,t])
  }

  # Grand average PhiMat (averaged over sites and years)
  PhiMat.avg[1,1] <- mean(PhiMat.annual[1,1,])
  PhiMat.avg[1,2] <- mean(PhiMat.annual[1,2,])
  PhiMat.avg[1,3] <- mean(PhiMat.annual[1,3,])
  PhiMat.avg[2,1] <- mean(PhiMat.annual[2,1,])
  PhiMat.avg[2,2] <- mean(PhiMat.annual[2,2,])
  PhiMat.avg[2,3] <- mean(PhiMat.annual[2,3,])
  PhiMat.avg[3,1] <- mean(PhiMat.annual[3,1,])
  PhiMat.avg[3,2] <- mean(PhiMat.annual[3,2,])
  PhiMat.avg[3,3] <- mean(PhiMat.annual[3,3,])

  # Number of sites in each state and occupied by single or pair
  for (t in 1:nyears){
    for (i in 1:nsites){
      state1[i,t] <- equals(z[i,t], 1)
      state2[i,t] <- equals(z[i,t], 2)
      state3[i,t] <- equals(z[i,t], 3)
      stateOcc[i,t] <- max(state2[i,t], state3[i,t])
    }
    n.occ[t,1] <- sum(state1[,t]) # Sites in state 1
    n.occ[t,2] <- sum(state2[,t]) # Sites in state 2
    n.occ[t,3] <- sum(state3[,t]) # Sites in state 3
    n.occ.total[t] <- sum(stateOcc[,t]) # All occupied
  }
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(bdata$nsites, bdata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("alpha.lpsi", "mean.psi", "beta.lpsi", "alpha.lr", "mean.r",
    "beta.lr", "alpha.lphi1", "alpha.lphi2", "alpha.lphi3", "alpha.lrho1",
    "alpha.lrho2", "alpha.lrho3", "mu.alpha.lphi1", "mean.phi1",
    "sd.alpha.lphi1", "mu.alpha.lphi2", "mean.phi2", "sd.alpha.lphi2",
    "mu.alpha.lphi3", "mean.phi3", "sd.alpha.lphi3", "mu.alpha.lrho1",
    "mean.rho1", "sd.alpha.lrho1", "mu.alpha.lrho2", "mean.rho2",
    "sd.alpha.lrho2", "mu.alpha.lrho3", "mean.rho3", "sd.alpha.lrho3",
    "beta.region.lphi1", "beta.region.lphi2", "beta.region.lphi3",
    "beta.region.lrho1", "beta.region.lrho2", "beta.region.lrho3",
    "beta.lphi1", "beta.lphi2", "beta.lphi3", "beta.lrho1", "beta.lrho2",
    "beta.lrho3", "alpha.lp2", "mean.alpha.p2", "alpha.lp32", "alpha.lp33",
    "beta.region.lp2", "beta.region.lp32", "beta.region.lp33", "beta.lp2",
    "beta.lp32", "beta.lp33", "PhiMat.avg", "n.occ", "n.occ.total", "z")

# MCMC settings
# na <- 1000 ; ni <- 60000 ; nt <- 20 ; nb <- 40000 ; nc <- 3
na <- 1000 ; ni <- 6000 ; nt <- 2 ; nb <- 4000 ; nc <- 3  # ~~~ for testing, 2 hrs

# Call JAGS (ART 22 h), check convergence and summarize posteriors
odms2 <- jags(bdata, inits, params, "dynMS2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(odms2)
summary(odms2) ; jags.View(odms2) # not shown

# ~~~ save output to use in next subsection ~~~
save(odms2, file="AHM2_06.4.3_odms2.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
