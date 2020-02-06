#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================

library(AHMbook)

# ~~~~ need to follow on from section 1.3

# 1.7 “Demographic” State-Space Models for Inference About Relative Abundance
# ===========================================================================

# 1.7.2 Demographic State-Space Models for Swiss Crested Tits
# -----------------------------------------------------------

# 1.7.2.1 Demographic SSM with Generalized Markovian Dynamics and With Overdispersion
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Bundle data (same as for Gaussian SSMs)
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]), T = ncol(C[sel,])))
# Specify model in BUGS language
cat(file = "model10d.txt","
model {
  # Priors
  # Model for expected initial abundance
  for(i in 1:M){
    alpha.lam[i] ~ dnorm(mu.alpha.lam, tau.alpha.lam)
  }
  mu.alpha.lam <- log(mean.lambda)
  mean.lambda ~ dunif(0, 50) # Mean of lambda
  tau.alpha.lam <- pow(sd.alpha.lam, -2)
  sd.alpha.lam ~ dnorm(0, 2) I(0.001,) # Site-level OD
  # Model for ’immigration-free’ population growth rate
  for(i in 1:M){
    for(t in 1:(T-1)){
      alpha.gam[i,t] ~ dnorm(mu.alpha.gam, tau.alpha.gam)
    }
  }
  mu.alpha.gam <- log(mean.gamma)
  mean.gamma ~ dunif(0.9, 1.1) # Mean of gamma
  tau.alpha.gam <- pow(sd.alpha.gam, -2)
  sd.alpha.gam ~ dnorm(0, 100) I(0.001,) # Site/ year-level OD
  # curve(dnorm(x, 0, sqrt(1 / 100)), 0, 0.5)
  # Model for detection probability
  for(i in 1:M){
    for(t in 1:T){
      alpha.p[i,t] ~ dnorm(mu.alpha.p, tau.alpha.p)
    }
  }
  mu.alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1) # Mean of p
  tau.alpha.p <- pow(sd.alpha.p, -2)
  sd.alpha.p ~ dnorm(0, 10) I(0.001,) # Site/ year-level OD

  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,) # Half-normal prior for sd
  # ’Likelihood’
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- loglam[i]
    loglam[i] <- alpha.lam[i]
    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1] + rho[t-1])
      log(gamma[i, t-1]) <- loggam[i, t-1]
      loggam[i, t-1] <- alpha.gam[i, t-1]
    }
    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p[i,t], N[i,t])
      logit(p[i,t]) <- lp[i,t]
      lp[i,t] <- alpha.p[i,t]
    }
  }
  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(N[,t])
  }
  for(t in 1:(T-1)){
    gammaX[t] <- popindex[t+1] / popindex[t] # Derived growth rate
  }
}
")
# Initial values
Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
inits <- function(){list(N = Nst + 1)}
# Parameters monitored
params <- c("mean.lambda", "mean.gamma", "mean.p", "mu.alpha.lam", "mu.alpha.gam",
  "mu.alpha.p", "sd.alpha.lam", "sd.alpha.gam", "sd.alpha.p", "sd.rho", "rho",
  "popindex", "gammaX", "alpha.lam", "alpha.gam", "alpha.p", "N")

# MCMC settings
# na <- 5000 ; ni <- 1e6 ; nt <- 500 ; nb <- 5e5 ; nc <- 3
na <- 5000 ; ni <- 1e5 ; nt <- 50 ; nb <- 5e4 ; nc <- 3  # ~~~~~ for testing
# Call JAGS (ART 491 min), check convergence and summarize posteriors
out10d <- jags(bdata, inits, params, "model10d.txt", n.adapt = na, n.chains = nc,
  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2,2)) ; traceplot(out10d) ; par(mfrow = c(1,1))
summary(out10d) ; View(out10d) ; print(out10d$summary[1:100,-c(4:6)], 2)

# 1.7.2.2 Demographic SSM With Generalized Markovian Dynamics and With Covariates
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# We assume you still have these covariates in your workspace, otherwise, you have to get them first.
# Bundle and summarize data
# str(bdata <- list(C = C[sel,], M = nrow(C[sel,]), T = ncol(C[sel,]).,
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]), T = ncol(C[sel,]),
    elev = elev.sc[sel], forest = forest.sc[sel], date = date.sc[sel,], 
    dur = dur.sc[sel,]))
# List of 7
# $ C : int [1:97, 1:18] 3 5 9 13 12 2 3 3 6 10 ...
# $ M : int 97
# $ T : int 18
# $ elev : num [1:97] -0.0614 -0.8418 -0.0614 -0.6857 0.0947 ...
# $ forest: num [1:97] 0.0092 0.9126 0.5512 0.8042 1.1294 ...
# $ date : num [1:97, 1:18] -0.6342 -0.5746 -0.6143 -1.0118 0.0614 ...
# $ dur : num [1:97, 1:18] -0.5418 -0.5418 0 0.9873 -0.0359 ...

# Specify model in BUGS language
cat(file = "model11.txt","
model {
  # Priors
  # Model for expected initial abundance
  alpha.lam <- log(mean.lambda)
  mean.lambda ~ dunif(0, 50) # Mean of lambda
  for(v in 1:3){ # Covariate coefficients
    beta.lam[v] ~ dnorm(0, 1)
  }
  # Model for ’immigration-free’ population growth rate
  alpha.gam <- log(mean.gamma)
  mean.gamma ~ dunif(0.9, 1.1) # Mean of gamma
  for(v in 1:3){ # Covariate coefficients
    beta.gam[v] ~ dnorm(0, 1)
  }
  # Model for detection probability
  alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1) # Mean of p
  for(v in 1:3){ # Covariate coefficients
    beta.p[v] ~ dnorm(0, 0.1)
  }
  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,) # Half-normal prior for sd
  # ’Likelihood’
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- loglam[i]
    loglam[i] <- alpha.lam + beta.lam[1] * elev[i] + beta.lam[2] * pow(elev[i],2) +
    beta.lam[3] * forest[i]
    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1] + rho[t-1])
      log(gamma[i, t-1]) <- loggam[i, t-1]
      loggam[i, t-1] <- alpha.gam + beta.gam[1] * elev[i] + 
        beta.gam[2] * pow(elev[i],2) + beta.gam[3] * forest[i]
    }
    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p[i,t], N[i,t])
      logit(p[i,t]) <- lp[i,t]
      lp[i,t] <- alpha.p + beta.p[1] * date[i,t] + beta.p[2] * pow(date[i,t],2) + 
        beta.p[3] * dur[i,t]
    }
  }
  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(N[,t])
  }
  for(t in 1:(T-1)){
    gammaX[t] <- popindex[t+1] / popindex[t] # Derived growth rate
  }
}
")

# Parameters monitored
params <- c("mean.lambda", "mean.gamma", "mean.p", "sd.rho", "beta.lam", "beta.gam",
  "beta.p", "rho", "popindex", "gammaX", "alpha.lam", "alpha.gam", "alpha.p", "N")
# MCMC settings
na <- 5000 ; ni <- 100000 ; nt <- 50 ; nb <- 50000 ; nc <- 3
# Call JAGS (ART 47 min), check convergence and summarize posteriors
out11 <- jags(bdata, inits, params, "model11.txt", n.adapt = na, n.chains = nc,
  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2,2)) ; traceplot(out11) ; par(mfrow = c(1,1))
summary(out11) ; View(out11) ; print(out11$summary[1:100,-c(4:6)], 2)

# 1.7.2.3 Comparison of the Inferences Under the Demographic State-Space Models

# no code
