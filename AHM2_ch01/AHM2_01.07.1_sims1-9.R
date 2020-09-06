# Applied hierarchical modeling in ecology - Vol 2 - Marc Kéry & Andy Royle

# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================

# Code based on MS "Report on 2nd set of sims: focused comparisons just for the
#  DM model", Marc Kery, 1 August 2018 and later, largely reorganised by Mike.

# 1.7 “Demographic” state-space models for inference about relative abundance
# 1.7.1 Simulation assessment of a demographic state-space model
# --------------------------------------------------------------

# This file has code to run the simulations for cases 1 to 9 described on
#  pp.33-34 and to plot figures 1.8 and 1.9.

library(AHMbook)
library(jagsUI)

# Specify the 5 models in BUGS language
# =====================================

# MODEL 1 : all fine
# ------------------
cat(file = "simModel1.txt","
model {
  # Priors
  alpha.lam ~ dnorm(0, 0.1)    # Expected initial abundance with X = 0
  mean.lambda <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 0.1)     # Coefficient of covariate on lambda
  alpha.gam ~ dunif(-0.2, 0.2) # Population growth rate with X = 0
  mean.gamma <- exp(alpha.gam)
  beta.gam ~ dnorm(0, 0.1)     # Coefficient of covariate on gamma
  alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1)         # Detection probability
  beta.p ~ dnorm(0, 0.1)       # Coefficient of covariate on p

  # Likelihood
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- loglam[i]
    loglam[i] <- alpha.lam + beta.lam * Xlam[i]

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1])
      log(gamma[i, t-1]) <- loggam[i, t-1]
      loggam[i, t-1] <- alpha.gam + beta.gam * Xgam[i, t-1]
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p[i,t], N[i,t])
      logit(p[i,t]) <- lp[i,t]
      lp[i,t] <- alpha.p + beta.p * Xp[i, t]
    }
  }

  # Derived quantities: total realized abundance per year (sumN) and annual growth rate based on sumN (gammaX)
  for(t in 1:T){
    sumN[t] <- sum(N[,t])   # Total realized abundance per year
  }
  for(t in 2:T){
    gammaX[t-1] <- sumN[t] / sumN[t-1] # Realized average population growth rate
  }
}
")

# Parameters to monitor
params1 <- c("mean.lambda", "mean.gamma", "mean.p", "alpha.lam", "alpha.gam",
    "alpha.p", "beta.lam", "beta.gam", "beta.p", "sumN", "gammaX", "N")

# MODEL 2 : Extended Markov dynamics (rescue effect or random immigration)
# ------------------------------------------------------------------------
cat(file = "simModel2.txt","
model {
  # Priors
  alpha.lam ~ dnorm(0, 0.1)    # Expected initial abundance with X = 0
  mean.lambda <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 0.1)     # Coefficient of covariate on lambda
  alpha.gam ~ dunif(-0.2, 0.2) # Population growth rate with X = 0
  mean.gamma <- exp(alpha.gam)
  beta.gam ~ dnorm(0, 0.1)     # Coefficient of covariate on gamma
  alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1)         # Detection probability
  beta.p ~ dnorm(0, 0.1)       # Coefficient of covariate on p

  for(t in 1:(T-1)){           # half-normal prior for random immigration
    log(eps[t]) <- logeps[t]
    logeps[t] ~ dnorm(0, tau.eps)
  }
  tau.eps <- pow(sd.eps, -2)
  sd.eps ~ dnorm(0, 5)I(0,)
  # curve(dnorm(x, 0, sqrt(1 / 5)), 0, 10)

  # Likelihood
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- loglam[i]
    loglam[i] <- alpha.lam + beta.lam * Xlam[i]

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1] + eps[t-1])
      log(gamma[i, t-1]) <- loggam[i, t-1]
      loggam[i, t-1] <- alpha.gam + beta.gam * Xgam[i, t-1]
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p[i,t], N[i,t])
      logit(p[i,t]) <- lp[i,t]
      lp[i,t] <- alpha.p + beta.p * Xp[i, t]
    }
  }

  # Derived quantities: total realized abundance per year (sumN) and annual growth rate based on sumN (gammaX)
  for(t in 1:T){
    sumN[t] <- sum(N[,t])   # Total realized abundance per year
  }
  for(t in 2:T){
    gammaX[t-1] <- sumN[t] / sumN[t-1] # Realized average population growth rate
  }
}
")

# Parameters to monitor
params2 <- c("mean.lambda", "mean.gamma", "mean.p", "alpha.lam", "alpha.gam",
    "alpha.p", "beta.lam", "beta.gam", "sd.eps", "beta.p", "sumN", "gammaX",
    "eps", "N")

# MODEL 3 : Extended Markov model with OD in lambda
# -------------------------------------------------
cat(file = "simModel3.txt","
model {
# Priors
for(i in 1:M){
  alpha.lam[i] ~ dnorm(mu.alpha.lam, tau.alpha.lam)
}
mu.alpha.lam ~ dnorm(0, 0.1) # Expected initial abundance with X = 0
mean.lambda <- exp(mu.alpha.lam)
beta.lam ~ dnorm(0, 0.1)     # Coefficient of covariate on lambda
tau.alpha.lam <- pow(sd.alpha.lam, -2)
sd.alpha.lam ~ dunif(0, 3)
alpha.gam ~ dunif(-0.2, 0.2) # Population growth rate with X = 0
mean.gamma <- exp(alpha.gam)
beta.gam ~ dnorm(0, 0.1)     # Coefficient of covariate on gamma
alpha.p <- logit(mean.p)
mean.p ~ dunif(0, 1)         # Detection probability
beta.p ~ dnorm(0, 0.1)       # Coefficient of covariate on p

for(t in 1:(T-1)){           # half-normal prior for random immigration
  log(eps[t]) <- logeps[t]
  logeps[t] ~ dnorm(0, tau.eps)
}
tau.eps <- pow(sd.eps, -2)
sd.eps ~ dnorm(0, 5)I(0,)
# curve(dnorm(x, 0, sqrt(1 / 5)), 0, 10)

# Likelihood
# State process
for(i in 1:M){
  # Initial conditions
  N[i,1] ~ dpois(lambda[i])
  log(lambda[i]) <- loglam[i]
  loglam[i] <- alpha.lam[i] + beta.lam * Xlam[i]

  # Transition model
  for(t in 2:T){
    N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1] + eps[t-1])
    log(gamma[i, t-1]) <- loggam[i, t-1]
    loggam[i, t-1] <- alpha.gam + beta.gam * Xgam[i, t-1]
  }

  # Observation process
  for(t in 1:T){
    C[i,t] ~ dbin(p[i,t], N[i,t])
    logit(p[i,t]) <- lp[i,t]
    lp[i,t] <- alpha.p + beta.p * Xp[i, t]
  }
}

# Derived quantities: total realized abundance per year (sumN) and annual growth rate based on sumN (gammaX)
for(t in 1:T){
  sumN[t] <- sum(N[,t])   # Total realized abundance per year
}
for(t in 2:T){
  gammaX[t-1] <- sumN[t] / sumN[t-1] # Realized average population growth rate
}
}
")

# Parameters to monitor
params3 <- c("mean.lambda", "mean.gamma", "mean.p", "mu.alpha.lam",
    "sd.alpha.lam", "alpha.gam", "alpha.p", "beta.lam", "beta.gam",
    "sd.eps", "beta.p", "sumN", "gammaX", "eps", "alpha.lam", "N")

# MODEL 4 : Extended Markov model with site-by-interval-specific OD in gamma
# --------------------------------------------------------------------------
cat(file = "simModel4.txt","
model {
# Priors
alpha.lam ~ dnorm(0, 0.1)    # Expected initial abundance with X = 0
mean.lambda <- exp(alpha.lam)
beta.lam ~ dnorm(0, 0.1)     # Coefficient of covariate on lambda

for(i in 1:M){
  for(t in 1:(T-1)){
    alpha.gam[i,t] ~ dnorm(mu.alpha.gam, tau.alpha.gam)
  }
}
mu.alpha.gam ~ dunif(-0.2, 0.2) # Population growth rate with X = 0
mean.gamma <- exp(mu.alpha.gam)
tau.alpha.gam <- pow(sd.alpha.gam, -2)
sd.alpha.gam ~ dunif(0, 3)


beta.gam ~ dnorm(0, 0.1)     # Coefficient of covariate on gamma
alpha.p <- logit(mean.p)
mean.p ~ dunif(0, 1)         # Detection probability
beta.p ~ dnorm(0, 0.1)       # Coefficient of covariate on p

for(t in 1:(T-1)){           # half-normal prior for random immigration
  log(eps[t]) <- logeps[t]
  logeps[t] ~ dnorm(0, tau.eps)
}
tau.eps <- pow(sd.eps, -2)
sd.eps ~ dnorm(0, 5)I(0,)
# curve(dnorm(x, 0, sqrt(1 / 5)), 0, 10)


# Likelihood
# State process
for(i in 1:M){
  # Initial conditions
  N[i,1] ~ dpois(lambda[i])
  log(lambda[i]) <- loglam[i]
  loglam[i] <- alpha.lam + beta.lam * Xlam[i]

  # Transition model
  for(t in 2:T){
    N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1] + eps[t-1])
    log(gamma[i, t-1]) <- loggam[i, t-1]
    loggam[i, t-1] <- alpha.gam[i,t-1] + beta.gam * Xgam[i, t-1]
  }

  # Observation process
  for(t in 1:T){
    C[i,t] ~ dbin(p[i,t], N[i,t])
    logit(p[i,t]) <- lp[i,t]
    lp[i,t] <- alpha.p + beta.p * Xp[i, t]
  }
}

# Derived quantities: total realized abundance per year (sumN) and annual growth rate based on sumN (gammaX)
for(t in 1:T){
  sumN[t] <- sum(N[,t])   # Total realized abundance per year
}
for(t in 2:T){
  gammaX[t-1] <- sumN[t] / sumN[t-1] # Realized average population growth rate
}
}
")

# Parameters to monitor
params4 <- c("mean.lambda", "mean.gamma", "mean.p",  "alpha.lam",
      "mu.alpha.gam", "sd.alpha.gam", "alpha.p", "sd.eps", "beta.lam",
      "beta.gam",  "beta.p", "sumN", "gammaX", "eps", "alpha.gam", "N")

# MODEL 5 : Extended Markov model with site-by-interval-specific OD in p
# ----------------------------------------------------------------------
cat(file = "simModel5.txt","
model {
  # Priors
  alpha.lam ~ dnorm(0, 0.1)    # Expected initial abundance with X = 0
  mean.lambda <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 0.1)     # Coefficient of covariate on lambda

  alpha.gam ~ dunif(-0.2, 0.2) # Population growth rate with X = 0
  mean.gamma <- exp(alpha.gam)
  beta.gam ~ dnorm(0, 0.1)     # Coefficient of covariate on gamma

  for(i in 1:M){
    for(t in 1:T){
      alpha.p[i,t] ~ dnorm(mu.alpha.p, tau.alpha.p)
    }
  }
  mu.alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1)         # Detection probability
  tau.alpha.p <- pow(sd.alpha.p, -2)
  sd.alpha.p ~ dunif(0, 3)
  beta.p ~ dnorm(0, 0.1)       # Coefficient of covariate on p

  for(t in 1:(T-1)){           # half-normal prior for random immigration
    log(eps[t]) <- logeps[t]
    logeps[t] ~ dnorm(0, tau.eps)
  }
  tau.eps <- pow(sd.eps, -2)
  sd.eps ~ dnorm(0, 5)I(0,)
  # curve(dnorm(x, 0, sqrt(1 / 5)), 0, 10)


  # Likelihood
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- loglam[i]
    loglam[i] <- alpha.lam + beta.lam * Xlam[i]

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1] + eps[t-1])
      log(gamma[i, t-1]) <- loggam[i, t-1]
      loggam[i, t-1] <- alpha.gam + beta.gam * Xgam[i, t-1]
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p[i,t], N[i,t])
      logit(p[i,t]) <- lp[i,t]
      lp[i,t] <- alpha.p[i,t] + beta.p * Xp[i, t]
    }
  }

  # Derived quantities: total realized abundance per year (sumN) and annual growth rate based on sumN (gammaX)
  for(t in 1:T){
    sumN[t] <- sum(N[,t])   # Total realized abundance per year
  }
  for(t in 2:T){
    gammaX[t-1] <- sumN[t] / sumN[t-1] # Realized average population growth rate
  }
}
")

# Parameters to monitor
params5 <- c("mean.lambda", "mean.gamma", "mean.p", "alpha.lam", "alpha.gam",
      "mu.alpha.p", "sd.alpha.p", "sd.eps", "beta.lam", "beta.gam",  "beta.p",
      "sumN", "gammaX", "eps", "alpha.p", "N")


# Prepare stuff needed
# ====================

# Choose number of simreps and create structures to save data
#  (true values) and estimates

# simrep <- 100
simrep <- 3  # ~~~ for testing
# Each rep takes about 1 hr

# For truth, we need the N's, sumN's and gammaX's for each of the
#  5 scenarios and each simrep
trueN <- array(NA, dim=c(100, 10, 5, simrep))
truesumN <- array(NA, dim=c(10, 5, simrep))
truegammaX <- array(NA, dim=c(9, 5, simrep))

# We need the estimates of those for each of the 9 cases
estiN <- array(NA, dim=c(100, 10, 9, simrep))
estisumN <- array(NA, dim=c(10, 9, simrep))
estigammaX <- array(NA, dim=c(9, 9, simrep))

# and some other stuff
estiEtc <- array(NA, dim=c(5, 9, simrep))
dimnames(estiEtc)[[1]] <- c("beta.lam", "beta.gam", "beta.p",
    "sd.eps", "sd.alpha.x")

seeds <- c(24, 26, 104)  # ~~~ for testing, these converge, which we need for the plots
# Other seeds are possible, provided
length(seeds) >= simrep


# Run the simulations for scenarios 1 to 5
# ========================================

system.time(
for(i in 1:simrep){
  cat(paste("\n***** Simrep", i, "of", simrep, "*****\n"))

  # Scenario 1: Markov simulation model, Markov analysis model,
  # no heterogeneity and all assumptions met (case 1)
  # -------------------------------------------------------------
  set.seed(seeds[i])
  betas <- 0.5                    # Value of all covariate coefficients
  data <- simPOP(mean.lam = 3, beta.lam = betas, mean.gamma = 1.0,
      beta.gamma = betas, sd.rho = 0, mean.p = 0.6, beta.p = betas,
      show.plot=FALSE)
  trueN[, , 1, i] <- data$N
  truesumN[, 1, i] <- data$sumN
  truegammaX[, 1, i] <- data$gammaX

  # Case 1 : Fit dynamic Nmix model
  cat("Case 1") ; flush.console()
  # Bundle data
  bdata <- list(C = data$C, M = nrow(data$C),  T = ncol(data$C),
      Xlam = data$Xsite1, Xgam = data$ Xsiteyear1, Xp = data$ Xsiteyear2)

  # Initial values
  Nst <- data$C + 10
  inits <- function(){list(N = Nst)}

  # MCMC settings
  na <- 1000  ;  ni <- 5000  ;  nt <- 5   ;  nb <- 0  ;  nc <- 3 # 2.2 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out1 <- try(jags(bdata, inits, params1, "simModel1.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence (of first 20 nodes), save results
  if(!inherits(out1, "try-error") &&
      !any(out1$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 1, i] <- out1$mean$N
    estisumN[, 1, i] <- out1$mean$sumN
    estigammaX[, 1, i] <- out1$mean$gammaX
    estiEtc[, 1, i] <- with(out1$mean, c(beta.lam, beta.gam, beta.p, NA, NA))
  }

  # Scenario 2: Extended Markov simulation model, but no heterogeneity
  #             (cases 2 and 3)
  # ---------------------------------------------------------------------
  set.seed(seeds[i])
  betas <- 0.5 # Value of all covariate coefficients
  sd.rho <- 0.2 # Value of random immigration parameter
  data <- simPOP(mean.lam = 3, beta.lam = betas, mean.gamma = 1.0,
      beta.gamma = betas, sd.rho = sd.rho, mean.p = 0.6, beta.p = betas,
      show.plot=FALSE)
  trueN[, , 2, i] <- data$N
  truesumN[, 2, i] <- data$sumN
  truegammaX[, 2, i] <- data$gammaX

  # Bundle data
  bdata <- list(C = data$C, M = nrow(data$C),  T = ncol(data$C),
      Xlam = data$Xsite1, Xgam = data$Xsiteyear1, Xp = data$Xsiteyear2)

  # Initial values
  Nst <- data$C + 10
  inits <- function(){list(N = Nst)}

  # Case 2 : Model 1 (Markov)
  cat(" - 2") ; flush.console()
  # MCMC settings
  na <- 5000  ;  ni <- 10000  ;  nt <- 10   ;  nb <- 0  ;  nc <- 3 # 5 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out2 <- try(jags(bdata, inits, params1, "simModel1.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence, save results
  if(!inherits(out2, "try-error") &&
      !any(out2$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 2, i] <- out2$mean$N
    estisumN[, 2, i] <- out2$mean$sumN
    estigammaX[, 2, i] <- out2$mean$gammaX
    estiEtc[, 2, i] <- with(out2$mean, c(beta.lam, beta.gam, beta.p, NA, NA))
  }

  # Case 3 : Model 2 (extended Markov)
  cat(" - 3") ; flush.console()
  # MCMC settings
  na <- 1000  ;  ni <- 10000  ;  nt <- 10   ;  nb <- 0  ;  nc <- 3 # 5 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out3 <- try(jags(bdata, inits, params2, "simModel2.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence, save results
  if(!inherits(out3, "try-error") &&
      !any(out3$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 3, i] <- out3$mean$N
    estisumN[, 3, i] <- out3$mean$sumN
    estigammaX[, 3, i] <- out3$mean$gammaX
    estiEtc[, 3, i] <- with(out3$mean, c(beta.lam, beta.gam, beta.p, sd.eps, NA))
  }

  # Scenario 3: Extended Markov simulation model with heterogeneity in lambda
  #             (cases 4 and 5)
  # -------------------------------------------------------------------------
  # Scenario 3: Extended Markov dynamics (rescue effect or random immigration) and heterogeneity in lambda
  set.seed(seeds[i])
  betas <- 0.5 # Value of all covariate coefficients
  sd.rho <- 0.2 # Value of random immigration parameter
  sd.log.lam <- 1 # Value of overdispersion in lambda
  data <- simPOP(mean.lam = 3, beta.lam = betas, sd.log.lam = sd.log.lam,
      mean.gamma = 1.0, beta.gamma = betas, sd.rho = sd.rho, mean.p = 0.6,
      beta.p = betas, show.plot=FALSE)
  trueN[, , 3, i] <- data$N
  truesumN[, 3, i] <- data$sumN
  truegammaX[, 3, i] <- data$gammaX


  # Bundle data
  bdata <- list(C = data$C, M = nrow(data$C),  T = ncol(data$C),
      Xlam = data$Xsite1, Xgam = data$Xsiteyear1, Xp = data$Xsiteyear2)

  # Gets inits for N
  Nst <- data$C + 10
  inits <- function(){list(N = Nst)}

  # Case 4 = Model 2 (extended Markov without heterogeneity)
  cat(" - 4") ; flush.console()
  # MCMC settings
  na <- 1000  ;  ni <- 10000  ;  nt <- 10   ;  nb <- 0  ;  nc <- 3  # 5 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out4 <- try(jags(bdata, inits, params2, "simModel2.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence, save results
  if(!inherits(out4, "try-error") &&
      !any(out4$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 4, i] <- out4$mean$N
    estisumN[, 4, i] <- out4$mean$sumN
    estigammaX[, 4, i] <- out4$mean$gammaX
    estiEtc[, 4, i] <- with(out4$mean, c(beta.lam, beta.gam, beta.p, sd.eps, NA))
  }

  # Case 5 = Model 3 (extended Markov with heterogeneity in lambda)
  cat(" - 5") ; flush.console()
  # MCMC settings
  na <- 5000  ;  ni <- 15000  ;  nt <- 15   ;  nb <- 0  ;  nc <- 3 # 7.5 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out5 <- try(jags(bdata, inits, params3, "simModel3.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence, save results
  if(!inherits(out5, "try-error") &&
      !any(out5$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 5, i] <- out5$mean$N
    estisumN[, 5, i] <- out5$mean$sumN
    estigammaX[, 5, i] <- out5$mean$gammaX
    estiEtc[, 5, i] <- with(out5$mean, c(beta.lam, beta.gam, beta.p,
        sd.eps, sd.alpha.lam))
  }

  # Scenario 4: Extended Markov simulation model with heterogeneity in gamma
  #             (cases 6 and 7)
  # ------------------------------------------------------------------------
  set.seed(seeds[i])
  betas <- 0.5 # Value of all covariate coefficients
  sd.log.gamma.survey <- 0.5 # Value of overdispersion in gamma
  data <- simPOP(mean.lam = 3, beta.lam = betas, mean.gamma = 1.0,
      beta.gamma = betas, sd.log.gamma.survey = sd.log.gamma.survey,
      sd.rho = 0.2, mean.p = 0.6, beta.p = betas, show.plot=FALSE)
  trueN[, , 4, i] <- data$N
  truesumN[, 4, i] <- data$sumN
  truegammaX[, 4, i] <- data$gammaX

  # Bundle data
  bdata <- list(C = data$C, M = nrow(data$C),  T = ncol(data$C),
      Xlam = data$Xsite1, Xgam = data$Xsiteyear1, Xp = data$Xsiteyear2)

  # Gets inits for N
  Nst <- data$C + 10
  inits <- function(){list(N = Nst)}

  # Case 6 : Model 2 (extended Markov without heterogeneity)
  cat(" - 6") ; flush.console()
  # MCMC settings
  na <- 1000  ;  ni <- 10000  ;  nt <- 10   ;  nb <- 0  ;  nc <- 3 # 5 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out6 <- try(jags(bdata, inits, params2, "simModel2.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence, save results
  if(!inherits(out6, "try-error") &&
      !any(out6$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 6, i] <- out6$mean$N
    estisumN[, 6, i] <- out6$mean$sumN
    estigammaX[, 6, i] <- out6$mean$gammaX
    estiEtc[, 6, i] <- with(out6$mean, c(beta.lam, beta.gam, beta.p, sd.eps, NA))
  }

  # Case 7 : Model 4 (extended Markov with heterogeneity in gamma)
  cat(" - 7") ; flush.console()
  # MCMC settings
  na <- 5000  ;  ni <- 30000  ;  nt <- 15   ;  nb <- 0  ;  nc <- 3 # 9 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out7 <- try(jags(bdata, inits, params4, "simModel4.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence, save results
  if(!inherits(out7, "try-error") &&
      !any(out7$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 7, i] <- out7$mean$N
    estisumN[, 7, i] <- out7$mean$sumN
    estigammaX[, 7, i] <- out7$mean$gammaX
    estiEtc[, 7, i] <- with(out7$mean, c(beta.lam, beta.gam, beta.p,
        sd.eps, sd.alpha.gam))
  }

  # Scenario 5: Extended Markov simulation model with site-by-occasion-level
  #             heterogeneity in p (cases 8 and 9)
  # ----------------------------------------------------------------------------
  set.seed(seeds[i])
  betas <- 0.5
  sd.logit.p.survey <- 1 # Value of overdispersion in p
  data <- simPOP(mean.lam = 3, beta.lam = betas, mean.gamma = 1.0,
      # beta.gamma = betas, sd.rho = 0.2, mean.p = 0.6,
      beta.gamma = betas, sd.rho = 0.2, mean.p = 0.6, beta.p = betas,
      sd.logit.p.survey = sd.logit.p.survey, show.plot=FALSE)
  trueN[, , 5, i] <- data$N
  truesumN[, 5, i] <- data$sumN
  truegammaX[, 5, i] <- data$gammaX

  # Bundle data
  bdata <- list(C = data$C, M = nrow(data$C),  T = ncol(data$C),
      Xlam = data$Xsite1, Xgam = data$Xsiteyear1, Xp = data$Xsiteyear2)

  # Gets inits for N
  Nst <- data$C + 10
  inits <- function(){list(N = Nst)}

  # Case 8 : Model 2 (extended Markov without heterogeneity)
  cat(" - 8") ; flush.console()
  # MCMC settings
  na <- 15000  ;  ni <- 40000  ;  nt <- 10   ;  nb <- 0  ;  nc <- 3  # 13 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out8 <- try(jags(bdata, inits, params2, "simModel2.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence, save results
  if(!inherits(out8, "try-error") &&
      !any(out8$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 8, i] <- out8$mean$N
    estisumN[, 8, i] <- out8$mean$sumN
    estigammaX[, 8, i] <- out8$mean$gammaX
    estiEtc[, 8, i] <- with(out8$mean, c(beta.lam, beta.gam, beta.p, sd.eps, NA))
  }

  # Case 9 : Model 5 (extended Markov with heterogeneity in p)
  cat(" - 9\n") ; flush.console()
  # MCMC settings
  na <- 10000  ;  ni <- 60000  ;  nt <- 10   ;  nb <- 0  ;  nc <- 3 # 15 mins

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out9 <- try(jags(bdata, inits, params5, "simModel5.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence, save results
  if(!inherits(out9, "try-error") &&
      !any(out9$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , 9, i] <- out9$mean$N
    estisumN[, 9, i] <- out9$mean$sumN
    estigammaX[, 9, i] <- out9$mean$gammaX
    estiEtc[, 9, i] <- with(out9$mean, c(beta.lam, beta.gam, beta.p,
        sd.eps, sd.alpha.p))
  }
} )  # 3 took 2.7 hrs; 20 reps took 22 hrs

save(seeds, trueN, truesumN, truegammaX, estiN, estisumN, estigammaX, estiEtc,
    file="sims_cases1-9_test.RData")

# Check for NAs
# =============
rbind(seeds, badRhat=round(colSums(is.na(estiEtc[1,,])), 2))
  # number of cases with NAs by seed

rbind(case=1:9, badRhat=round(rowSums(is.na(estiEtc[1,,])), 2))
  # number of simulations with poor convergence by case
# Cases 2, 7 and 9 have a high proportion of NAs due to non-convergence.

# Prepare arrays for plotting
# ===========================

# Fix true* arrays to match cases instead of scenarios
trueNp <- trueN[, , c(1,2,2,3,3,4,4,5,5), ]
str(trueNp)
truesumNp <- truesumN[, c(1,2,2,3,3,4,4,5,5), ]
truegammaXp <- truegammaX[, c(1,2,2,3,3,4,4,5,5), ]


# Figure 1.8
# ==========

op <- par(mfrow=c(3,5), ask=dev.interactive(orNone=TRUE))
for(case in 1:9) {
  plot(trueNp[,1,case,], estiN[,1,case,], pch=16,
      ylab=paste("Case", case), xlab="Truth", main="N in Year 1")
  abline(0, 1, lwd=2, col='red')
  abline(lm(c(estiN[,1,case,]) ~ c(trueNp[,1,case,])),
      lwd=2, lty=2, col='blue')

  plot(trueNp[,10,case,], estiN[,10,case,], pch=16,
      ylab=paste("Case", case), xlab="Truth", main="N in Year 10")
  abline(0, 1, lwd=2, col='red')
  abline(lm(c(estiN[,10,case,]) ~ c(trueNp[,10,case,])),
      lwd=2, lty=2, col='blue')

  plot(truesumNp[,case,], estisumN[,case,], pch=16,
      ylab=paste("Case", case), xlab="Truth", main="Total N")
  abline(0, 1, lwd=2, col='red')
  abline(lm(c(estisumN[,case,]) ~ c(truesumNp[,case,])),
      lwd=2, lty=2, col='blue')

  plot(truegammaXp[,case,], estigammaX[,case,], pch=16,
      ylab=paste("Case", case), xlab="Truth", main="Growth rate")
  abline(0, 1, lwd=2, col='red')
  abline(lm(c(estigammaX[,case,]) ~ c(truegammaXp[,case,])),
      lwd=2, lty=2, col='blue')

  if(case < 3) {
    plot.new()
  } else {
    plot(density(estiEtc['sd.eps', case, ], na.rm=TRUE), col = 'black',
        main = "Rand. immigration",
        xlab = '', ylab = '', lwd = 2, frame = FALSE)
    abline(v = 0.2, col = "red", lwd = 2)
    abline(v = mean(estiEtc['sd.eps', case, ], na.rm=TRUE),
        col = "blue", lwd = 2, lty = 2)
  }
}
par(op)

# Figure 1.9
# ==========

op <- par(mfrow=c(3,4), ask=dev.interactive(orNone=TRUE))
for(case in 1:9) {
  plot(density(estiEtc['beta.lam', case, ], na.rm=TRUE), col = 'black',
      main = "lambda covariate", xlim=c(0,1),
      xlab = '', ylab = paste("Case", case), lwd = 2, frame = FALSE)
  abline(v = 0.5, col = "red", lwd = 2)
  abline(v = mean(estiEtc['beta.lam', case, ], na.rm=TRUE),
      col = "blue", lwd = 2, lty = 2)

  plot(density(estiEtc['beta.gam', case, ], na.rm=TRUE), col = 'black',
      main = "gamma covariate", xlim=c(0,1),
      xlab = '', ylab = paste("Case", case), lwd = 2, frame = FALSE)
  abline(v = 0.5, col = "red", lwd = 2)
  abline(v = mean(estiEtc['beta.gam', case, ], na.rm=TRUE),
      col = "blue", lwd = 2, lty = 2)

  plot(density(estiEtc['beta.p', case, ], na.rm=TRUE), col = 'black',
      main = "p covariate", xlim=c(0,1),
      xlab = '', ylab = paste("Case", case), lwd = 2, frame = FALSE)
  abline(v = 0.5, col = "red", lwd = 2)
  abline(v = mean(estiEtc['beta.p', case, ], na.rm=TRUE),
      col = "blue", lwd = 2, lty = 2)

  if(case == 5) {
    plot(density(estiEtc['sd.alpha.x', case, ], na.rm=TRUE), col = 'black',
        main = "SD lambda", xlim=c(0,2),
        xlab = '', ylab = paste("Case", case), lwd = 2, frame = FALSE)
    abline(v = 1, col = "red", lwd = 2)
    abline(v = mean(estiEtc['sd.alpha.x', case, ], na.rm=TRUE),
        col = "blue", lwd = 2, lty = 2)
  } else if(case == 7) {
    plot(density(estiEtc['sd.alpha.x', case, ], na.rm=TRUE), col = 'black',
        main = "SD gamma", xlim=c(0,2),
        xlab = '', ylab = paste("Case", case), lwd = 2, frame = FALSE)
    abline(v = 0.5, col = "red", lwd = 2)
    abline(v = mean(estiEtc['sd.alpha.x', case, ], na.rm=TRUE),
        col = "blue", lwd = 2, lty = 2)
  } else if(case == 9) {
    plot(density(estiEtc['sd.alpha.x', case, ], na.rm=TRUE), col = 'black',
        main = "SD p", xlim=c(0,2),
        xlab = '', ylab = paste("Case", case), lwd = 2, frame = FALSE)
    abline(v = 1, col = "red", lwd = 2)
    abline(v = mean(estiEtc['sd.alpha.x', case, ], na.rm=TRUE),
        col = "blue", lwd = 2, lty = 2)
  } else {
    plot.new()
  }
}
par(op)
