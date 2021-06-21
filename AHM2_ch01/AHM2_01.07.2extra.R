#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from MS dated 2018-08-14

# Approximate run time for this script: 46 mins
# Run time with the full number of iterations:

# This file has code for models 8, 9, 10a-c and 12 described in
#   section 1.7.2, pp.38-39 and figure 1.13.


library(AHMbook)
library(jagsUI)

# ~~~~~ Need to run 1.3 before this ~~~~~~~
source("AHM2_01.03.R")
# ~~~~~ and this from 1.4 ~~~~~~~~~~~~~~~~~
M <- nrow(C)
T <- ncol(C)
# ~~~~~ and this from 1.5.3 ~~~~~~~~~~~~~~~~~
# Scale some covariates and mean-impute missing values in them
elev.sc <- standardize(dat$elev) # elevation of site
forest.sc <- standardize(dat$forest) # forest cover of site
date.sc <- standardize(date)
date.sc[is.na(date.sc)] <- 0 # mean impute
dur.sc <- standardize(dur)
dur.sc[is.na(dur.sc)] <- 0 # mean impute
# ~~~~~ and this from 1.5.4 ~~~~~~~~~~~~~~~~~
nzero <- apply(C, 1, function(x) sum(x == 0, na.rm = TRUE))
sel <- nzero <= 1 # Select sites with <= 1 zero count
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1.7 “Demographic” state-space models for inference about relative abundance
# ===========================================================================

# 1.7.2 Demographic state-space models for Swiss crested tits
# -----------------------------------------------------------

# Model 8: Markov model with constant parameters
# ''''''''''''''''''''''''''''''''''''''''''''''
# Bundle data
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]),  T = ncol(C[sel,])))
# List of 3
# $ C: int [1:97, 1:18] 3 5 9 13 12 2 3 3 6 10 ...
# $ M: int 97
# $ T: int 18

# Specify model in BUGS language
cat(file = "model8.txt","
model {

  # Priors
  lambda ~ dunif(0, 100)     # Expected initial abundance
  gamma ~ dunif(0.5, 1.5)    # Population growth rate
  p ~ dunif(0, 1)            # Detection probability

  # Likelihood
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda)

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma)
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p, N[i,t])
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(N[,t])
  }
  for(t in 1:(T-1)){
    gammaX[t] <- popindex[t+1] / popindex[t]    # Derived growth rate
  }
}
")

# Initial values
Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
inits <- function(){list(N = Nst + 1)}

# Parameters monitored
params <- c("lambda", "gamma", "p", "popindex", "gammaX", "N")

# MCMC settings
na <- 5000  ;  ni <- 200000  ;  nt <- 10   ;  nb <- 110000  ;  nc <- 3  # 20 mins
# na <- 5000  ;  ni <- 20000  ;  nt <- 10   ;  nb <- 10000  ;  nc <- 3  # 2.5 mins

# Call JAGS from R, check convergence and summarize marginal posteriors
out8 <- jags(bdata, inits, params, "model8.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 2))  #  ~~~ replace with 'layout' argument
traceplot(out8, layout=c(2,2))
print(out8, 3)


# Model 9: Generalized Markov model with constant parameters
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Bundle data
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]),  T = ncol(C[sel,])))

# Specify model in BUGS language
cat(file = "model9.txt","
model {

  # Priors
  lambda ~ dunif(0, 100)     # Expected initial abundance
  gamma ~ dunif(0.5, 1.5)    # 'Immigration-free' population growth rate
  p ~ dunif(0, 1)            # Detection probability
  for(t in 1:(T-1)){         # Model for random immigration
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,)  # Half-normal prior for sd
  # curve(dnorm(x, 0, sqrt(1 / 0.5)), 0, 5)  # how does it look like ?

  # Likelihood
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda)

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma + rho[t-1])
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p, N[i,t])
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(N[,t])
  }
  for(t in 1:(T-1)){
    gammaX[t] <- popindex[t+1] / popindex[t]    # Derived growth rate
  }
}
")

# Initial values
Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
inits <- function(){list(N = Nst + 1)}

# Parameters monitored
params <- c("lambda", "gamma", "p", "sd.rho", "rho", "popindex", "gammaX", "N")

# MCMC settings
na <- 5000  ;  ni <- 40000  ;  nt <- 10   ;  nb <- 20000  ;  nc <- 3  # 5 mins
# na <- 5000  ;  ni <- 20000  ;  nt <- 10   ;  nb <- 10000  ;  nc <- 3  # 3 mins

# Call JAGS from R, check convergence and summarize marginal posteriors
out9 <- jags(bdata, inits, params, "model9.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 2))  #  ~~~ replace with 'layout' argument
traceplot(out9, layout=c(2,2))
print(out9, 3)


# Model 10a: Generalized Markov model with site-level overdispersion
#            in lambda only
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Bundle data
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]),  T = ncol(C[sel,])))

# Use weakly informative priors, for instance, sufficiently wide half-normal
# distributions to constrain the standard deviations of overdispersion
# parameters away from large values.

cat(file = "model10a.txt","
model {

  # Priors
  # Model for expected initial abundance
  for(i in 1:M){
    alpha.lam[i] ~ dnorm(mu.alpha.lam, tau.alpha.lam)
  }
  mu.alpha.lam <- log(mean.lambda)
  mean.lambda ~ dunif(0, 50)     # Mean
  tau.alpha.lam <- pow(sd.alpha.lam, -2)
  # sd.alpha.lam ~ dunif(0, 3)      # Site-level overdispersion
  sd.alpha.lam ~ dnorm(0, 2) I(0.001,)# Site-level overdispersion
  # curve(dnorm(x, 0, sqrt(1 / 2)), 0, 3)

  gamma ~ dunif(0.8, 1.2)         # 'Immigration-free' pop.growth rate
  p ~ dunif(0, 1)                 # Detection probability

  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,)  # Half-normal prior for sd
  # curve(dnorm(x, 0, sqrt(1 / 0.5)), 0, 5)  # how does it look like ?

  # Likelihood
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- loglam[i]
    loglam[i] <- alpha.lam[i]

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma + rho[t-1])
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p, N[i,t])
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(N[,t])
  }
  for(t in 1:(T-1)){
    gammaX[t] <- popindex[t+1] / popindex[t]    # Derived growth rate
  }
}
")

# Initial values
Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
inits <- function(){list(N = Nst + 1)}

# Parameters monitored
params <- c("mean.lambda", "mu.alpha.lam", "sd.alpha.lam", "gamma", "p",
    "sd.rho", "rho", "popindex", "gammaX", "N")

# MCMC settings
na <- 5000  ;  ni <- 40000  ;  nt <- 10   ;  nb <- 20000  ;  nc <- 3  # 5.5 mins
# na <- 100  ;  ni <- 150  ;  nt <- 1  ;  nb <- 50  ;  nc <- 3     # Test

# Call JAGS from R, check convergence and summarize marginal posteriors
out10a <- jags(bdata, inits, params, "model10a.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 2))  #  ~~~ replace with 'layout' argument
traceplot(out10a, layout=c(2,2))
print(out10a, 3)


# Model 10b: Generalized Markov model with site-by-year-level oversdispersion
#  in gamma only
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Bundle data
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]),  T = ncol(C[sel,])))

# Specify model in BUGS language
cat(file = "model10b.txt","
model {

  # Priors
  lambda ~ dunif(0, 50)          # Expected initial abundance

  # Model for 'immigration-free' population growth rate
  for(i in 1:M){
    for(t in 1:(T-1)){
      alpha.gam[i,t] ~ dnorm(mu.alpha.gam, tau.alpha.gam)
    }
  }
  mu.alpha.gam <- log(mean.gamma)
  mean.gamma ~ dunif(0.8, 1.2)    # Mean
  tau.alpha.gam <- pow(sd.alpha.gam, -2)
  #sd.alpha.gam ~ dunif(0, 1)      # Site/year-level overdispersion
  sd.alpha.gam ~ dnorm(0, 100) I(0.001,) # Site/year-level overdispersion
  # curve(dnorm(x, 0, sqrt(1 / 100)), 0, 0.3)

  p ~ dunif(0, 1)                 # Detection probability

  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,)  # Half-normal prior for sd
  # curve(dnorm(x, 0, sqrt(1 / 0.5)), 0, 5)  # how does it look like ?

  # Likelihood
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda)

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1] + rho[t-1])
      log(gamma[i, t-1]) <- loggam[i, t-1]
      loggam[i, t-1] <- alpha.gam[i, t-1]
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p, N[i,t])
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(N[,t])
  }
  for(t in 1:(T-1)){
    gammaX[t] <- popindex[t+1] / popindex[t]    # Derived growth rate
  }
}
")

# Initial values
Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
inits <- function(){list(N = Nst + 1)}

# Parameters monitored
params <- c("lambda", "mean.gamma", "mu.alpha.gam", "sd.alpha.gam", "p",
    "sd.rho", "rho", "popindex", "gammaX", "N")

# MCMC settings
na <- 1000  ;  ni <- 100000  ;  nt <- 50   ;  nb <- 50000  ;  nc <- 3  # 20 mins
# na <- 100  ;  ni <- 150  ;  nt <- 1  ;  nb <- 50  ;  nc <- 3     # Test

# Call JAGS from R, check convergence and summarize marginal posteriors
out10b <- jags(bdata, inits, params, "model10b.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 2))  #  ~~~ replace with 'layout' argument
traceplot(out10b, layout=c(2,2))
print(out10b, 3)

# Model 10c: Generalized Markov model with site-by-year-level
#            oversdispersion in p
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Bundle data
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]),  T = ncol(C[sel,])))

# Specify model in BUGS language
cat(file = "model10c.txt","
model {

  # Priors
  lambda ~ dunif(0, 50)      # Initial expected abundance
  gamma ~ dunif(0.8, 1.2)    # 'Immigration-free' population growth rate

  # Model for detection probability
  for(i in 1:M){
    for(t in 1:T){
      alpha.p[i,t] ~ dnorm(mu.alpha.p, tau.alpha.p)
    }
  }
  mu.alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1)            # Mean
  tau.alpha.p <- pow(sd.alpha.p, -2)
  #sd.alpha.p ~ dunif(0, 1)        # Site/year-level overdispersion
  sd.alpha.p ~ dnorm(0, 10) I(0.001,) # Site/year-level overdispersion
  # curve(dnorm(x, 0, sqrt(1 / 10)), 0, 1)


  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,)  # Half-normal prior for sd
  # curve(dnorm(x, 0, sqrt(1 / 0.5)), 0, 5)  # how does it look like ?

  # Likelihood
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda)

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma + rho[t-1])
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
    gammaX[t] <- popindex[t+1] / popindex[t]    # Derived growth rate
  }
}
")

# Initial values
Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
inits <- function(){list(N = Nst + 1)}

# Parameters monitored
params <- c("lambda", "gamma", "mean.p", "mu.alpha.p", "sd.alpha.p",
    "sd.rho", "rho", "popindex", "gammaX", "N")

# MCMC settings
na <- 1000  ;  ni <- 10^6  ;  nt <- 500   ;  nb <- 5*10^5  ;  nc <- 3
# na <- 1000  ;  ni <- 10^4  ;  nt <- 5   ;  nb <- 5*10^3  ;  nc <- 3

# Call JAGS from R, check convergence and summarize marginal posteriors
out10c <- jags(bdata, inits, params, "model10c.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 2))  #  ~~~ replace with 'layout' argument
traceplot(out10c, layout=c(2,2))
print(out10a, 3)


# Model 12: Generalized Markov model with site-level overdispersion in lambda
#  and site-by-year-level overdispersion in both gamma and p and with all the
#  covariates from model 4
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Bundle and summarize data
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]),  T = ncol(C[sel,]),
    elev = elev.sc[sel], forest = forest.sc[sel], date = date.sc[sel,],
    dur = dur.sc[sel,]))
# List of 7
# $ C     : int [1:97, 1:18] 3 5 9 13 12 2 3 3 6 10 ...
# $ M     : int 97
# $ T     : int 18
# $ elev  : num [1:97] -0.0614 -0.8418 -0.0614 -0.6857 0.0947 ...
# $ forest: num [1:97] 0.0092 0.9126 0.5512 0.8042 1.1294 ...
# $ date  : num [1:97, 1:18] -0.6342 -0.5746 -0.6143 -1.0118 0.0614 ...
# $ dur   : num [1:97, 1:18] -0.5418 -0.5418 0 0.9873 -0.0359 ...

# Continue to use weakly informative priors !

# Specify model in BUGS language
cat(file = "model12.txt","
model {

  # Priors
  # Model for expected initial abundance
  for(i in 1:M){
    alpha.lam[i] ~ dnorm(mu.alpha.lam, tau.alpha.lam)
  }
  mu.alpha.lam <- log(mean.lambda)
  mean.lambda ~ dunif(0, 100)     # Mean
  for(v in 1:3){     # Covariate coefficients
    beta.lam[v] ~ dnorm(0, 1)
  }
  tau.alpha.lam <- pow(sd.alpha.lam, -2)
  #sd.alpha.lam ~ dunif(0, 1)      # Site-level overdispersion
  sd.alpha.lam ~ dnorm(0, 2) I(0.001,)   # Site-level overdispersion
  # curve(dnorm(x, 0, sqrt(1 / 2)), 0, 3)

  # Model for 'immigration-free' population growth rate
  for(i in 1:M){
    for(t in 1:(T-1)){
      alpha.gam[i,t] ~ dnorm(mu.alpha.gam, tau.alpha.gam)
    }
  }
  mu.alpha.gam <- log(mean.gamma)
  mean.gamma ~ dunif(0.9, 1.1)    # Mean
  for(v in 1:3){     # Covariate coefficients
    beta.gam[v] ~ dnorm(0, 1)
  }

  tau.alpha.gam <- pow(sd.alpha.gam, -2)
  #sd.alpha.gam ~ dunif(0, 1)      # Site/year-level overdispersion
  sd.alpha.gam ~ dnorm(0, 100) I(0.001,) # Site/year-level overdispersion
  # curve(dnorm(x, 0, sqrt(1 / 100)), 0, 0.3)


  # Model for detection probability
  for(i in 1:M){
    for(t in 1:T){
      alpha.p[i,t] ~ dnorm(mu.alpha.p, tau.alpha.p)
    }
  }
  mu.alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1)            # Mean
  for(v in 1:3){     # Covariate coefficients
    beta.p[v] ~ dnorm(0, 0.1)
  }
  tau.alpha.p <- pow(sd.alpha.p, -2)
  #sd.alpha.p ~ dunif(0, 1)        # Site/year-level overdispersion
  sd.alpha.p ~ dnorm(0, 10) I(0.001,) # Site/year-level overdispersion
  # curve(dnorm(x, 0, sqrt(1 / 10)), 0, 1)

  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,)  # Half-normal prior for sd
  # curve(dnorm(x, 0, sqrt(1 / 0.5)), 0, 5)  # how does it look like ?

  # 'Likelihood'
  # State process
  for(i in 1:M){
    # Initial conditions
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- loglam[i]
    loglam[i] <- alpha.lam[i] + beta.lam[1] * elev[i] + beta.lam[2] * pow(elev[i],2) + beta.lam[3] * forest[i]

    # Transition model
    for(t in 2:T){
      N[i,t] ~ dpois(N[i,t-1] * gamma[i, t-1] + rho[t-1])
      log(gamma[i, t-1]) <- loggam[i, t-1]
      loggam[i, t-1] <- alpha.gam[i, t-1] + beta.gam[1] * elev[i] + beta.gam[2] * pow(elev[i],2) + beta.gam[3] * forest[i]
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p[i,t], N[i,t])
      logit(p[i,t]) <- lp[i,t]
      lp[i,t] <- alpha.p[i,t] + beta.p[1] * date[i,t] + beta.p[2] * pow(date[i,t],2) + beta.p[3] * dur[i,t]
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(N[,t])
  }
  for(t in 1:(T-1)){
    gammaX[t] <- popindex[t+1] / popindex[t]    # Derived growth rate
  }
}
")


# Initial values
Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
inits <- function(){list(N = Nst + 1)}

# Parameters monitored
params <- c("mean.lambda", "mean.gamma", "mean.p", "mu.alpha.lam",
    "mu.alpha.gam", "mu.alpha.p", "sd.alpha.lam", "sd.alpha.gam",
    "sd.alpha.p", "sd.rho", "beta.lam", "beta.gam", "beta.p", "rho",
    "popindex", "gammaX", "alpha.lam", "alpha.gam", "alpha.p", "N")

# MCMC settings
# na <- 10000  ;  ni <- 1.5 * 10^6  ;  nt <- 500   ;  nb <- 10^6  ;  nc <- 3
na <- 1000  ;  ni <- 15000  ;  nt <- 5  ;  nb <- 10000  ;  nc <- 3     # Test, 7 mins

# Call JAGS from R, check convergence and summarize marginal posteriors
out12 <- jags(bdata, inits, params, "model12.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 2))  #  ~~~ replace with 'layout' argument
traceplot(out12, layout=c(2,2))
print(out12$summary[1:100,], 3)

# ~~~ save ouput for future use ~~~~
save(out8, out9, out10a, out10b, out10c, out12, file="AHM2_01.07.2extra_outs.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ extra code for figure 1.13 ~~~~
# Compare population inferences among models (= sumN)
# Load saved output
load("AHM2_01.07.2_out10d.RData")
load("AHM2_01.07.2_out11.RData")

# Grab necessary estimates
predN <- cbind(out8$mean$popindex, out9$mean$popindex, out10a$mean$popindex,
    out10b$mean$popindex, out10c$mean$popindex, out10d$mean$popindex,
    out11$mean$popindex)
LCI <- cbind(out8$q2.5$popindex, out9$q2.5$popindex, out10a$q2.5$popindex,
    out10b$q2.5$popindex, out10c$q2.5$popindex, out10d$q2.5$popindex,
    out11$q2.5$popindex)
UCI <- cbind(out8$q97.5$popindex, out9$q97.5$popindex, out10a$q97.5$popindex,
    out10b$q97.5$popindex, out10c$q97.5$popindex, out10d$q97.5$popindex,
    out11$q97.5$popindex)
# predTrend <- cbind(out8$mean$gammaX[-1], out9$mean$gammaX[-1], out10a$mean$gammaX[-1],
    # out10b$mean$gammaX[-1], out10c$mean$gammaX[-1], out10d$mean$gammaX, out11$mean$gammaX)#[-1 ,] # toss out first year NA !
predTrend <- cbind(out8$mean$gammaX, out9$mean$gammaX, out10a$mean$gammaX,
    out10b$mean$gammaX, out10c$mean$gammaX, out10d$mean$gammaX, out11$mean$gammaX)

op <- par(mfrow = c(2, 1), mar = c(1,4,2,2)+0.1)
# Estimated trajectories
xlim = c(1998.5, 2016.5)  ;  ylim = c(500, 6000)
off <- 0.1
plot(year-3*off, predN[,1], type = 'b', col = 1, xlim = xlim, ylim = ylim,
    main = '', xlab = '', ylab = 'Population size', frame = FALSE, las = 1,
    pch=16, axes = FALSE)
axis(1, labels = FALSE)
axis(2, las = 1)
points(year-2*off, predN[,2], type = 'b', pch=16, col = 2)
points(year-1*off, predN[,3], type = 'b', pch=16, col = 3)
points(year, predN[,4], type = 'b', pch=16, col = 4)
points(year+1*off, predN[,5], type = 'b', pch=16, col = 5)
points(year+2*off, predN[,6], type = 'b', pch=16, col = 6)
points(year+3*off, predN[,7], type = 'b', pch=16, col = 7)
for(i in 1:7){
  segments((year-4*off+i*off), LCI[,i], (year-4*off+i*off), UCI[,i], col = i)
}

# Compare gammaX among models (= annual trends)
par(mar = c(5,4,1,2)+0.1)
ylim = c(0.9, 1.4)
matplot(year[-length(year)], predTrend, type = 'b', pch = 16, lty = 1,
    xlab = 'Year', col = 1:7, xlim = xlim, ylim = ylim, main = '',
    ylab = 'Growth rate', las = 1, frame = FALSE, las=1)
abline(h = 1, col = 'black')
legend('topright', c('Model 8: Null, Markov',
    'Model 9: Null, generalized Markov', 'Model 10a: OD in lambda',
    'Model 10b: OD in gamma',
    'Model 10c: OD in p', 'Model 10d: OD in all three',
    'Model 11: with covariates'), col=1:7, pch=16, bty = 'n', cex=0.8)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
