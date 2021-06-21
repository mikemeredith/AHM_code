#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 1 hr
# Run time with the full number of iterations: 7.4 hrs

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

# 1.7.2.1 Demographic SSM with generalized markovian dynamics and with overdispersion
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
  sd.alpha.lam ~ dnorm(0, 2) I(0.001,)    # Site-level OD

  # Model for 'immigration-free' population growth rate
  for(i in 1:M){
    for(t in 1:(T-1)){
      alpha.gam[i,t] ~ dnorm(mu.alpha.gam, tau.alpha.gam)
    }
  }
  mu.alpha.gam <- log(mean.gamma)
  mean.gamma ~ dunif(0.9, 1.1) # Mean of gamma
  tau.alpha.gam <- pow(sd.alpha.gam, -2)
  sd.alpha.gam ~ dnorm(0, 100) I(0.001,)   # Site/ year-level OD
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
  sd.alpha.p ~ dnorm(0, 10) I(0.001,)      # Site/ year-level OD

  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,)          # Half-normal prior for sd

  # 'Likelihood'
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
params <- c("mean.lambda", "mean.gamma", "mean.p", "mu.alpha.lam",
    "mu.alpha.gam", "mu.alpha.p", "sd.alpha.lam", "sd.alpha.gam",
    "sd.alpha.p", "sd.rho", "rho", "popindex", "gammaX", "alpha.lam",
    "alpha.gam", "alpha.p", "N")

# MCMC settings
# na <- 5000 ; ni <- 1e6 ; nt <- 500 ; nb <- 5e5 ; nc <- 3
na <- 5000 ; ni <- 1e5 ; nt <- 50 ; nb <- 5e4 ; nc <- 3  # ~~~~~ for testing, 40 mins

# Call JAGS (ART 491 min), check convergence and summarize posteriors
out10d <- jags(bdata, inits, params, "model10d.txt", n.adapt = na, n.chains = nc,
  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  # ~~~ replaced with 'layout' argument
traceplot(out10d, layout=c(2,2))
summary(out10d) ; jags.View(out10d) ; print(out10d$summary[1:100,-c(4:6)], 2)

# ~~~ Save output for use in subsequent sections ~~~
save(out10d, file="AHM2_01.07.2_out10d.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1.7.2.2 Demographic SSM with generalized Markovian dynamics and with covariates
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Bundle and summarize data
str(bdata <- list(C = C[sel,], M = nrow(C[sel,]), T = ncol(C[sel,]),
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

# Specify model in BUGS language
cat(file = "model11.txt","
model {

  # Priors
  # Model for expected initial abundance
  alpha.lam <- log(mean.lambda)
  mean.lambda ~ dunif(0, 50)             # Mean of lambda
  for(v in 1:3){                         # Covariate coefficients
    beta.lam[v] ~ dnorm(0, 1)
  }

  # Model for 'immigration-free' population growth rate
  alpha.gam <- log(mean.gamma)
  mean.gamma ~ dunif(0.9, 1.1)           # Mean of gamma
  for(v in 1:3){                         # Covariate coefficients
    beta.gam[v] ~ dnorm(0, 1)
  }

  # Model for detection probability
  alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1)                   # Mean of p
  for(v in 1:3){                         # Covariate coefficients
    beta.p[v] ~ dnorm(0, 0.1)
  }

  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ dnorm(0, 0.5)I(0.001,)        # Half-normal prior for sd

  # 'Likelihood'
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
    gammaX[t] <- popindex[t+1] / popindex[t]    # Derived growth rate
  }
}
")

# Parameters monitored
params <- c("mean.lambda", "mean.gamma", "mean.p", "sd.rho", "beta.lam",
    "beta.gam", "beta.p", "rho", "popindex", "gammaX", "alpha.lam",
    "alpha.gam", "alpha.p", "N")

# MCMC settings
na <- 5000 ; ni <- 100000 ; nt <- 50 ; nb <- 50000 ; nc <- 3

# Call JAGS (ART 47 min), check convergence and summarize posteriors
out11 <- jags(bdata, inits, params, "model11.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out11, layout=c(2,2))
summary(out11) ; jags.View(out11) ; print(out11$summary[1:100,-c(4:6)], 2)

# ~~~ Save output for use in subsequent sections ~~~
save(out11, file="AHM2_01.07.2_out11.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1.7.2.3 Comparison of the inferences under the demographic state-space models
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# ~~~~~~~~ code to produce Fig 1.11 ~~~~~~~~~~~~~~~
op <- par(mfrow = c(2,2))
plot(density(out10d$sims.list$sd.alpha.lam), col = 'black', lwd = 2,
    main = "SD of OD in lambda", xlab = 'sd.alpha.lam', ylab = 'Density',
    frame = FALSE, axes = FALSE, xlim = c(0, 1))
axis(1)  ;  axis(2)

plot(density(out10d$sims.list$sd.alpha.gam), col = 'black', lwd = 2,
    main = "SD of OD in gamma", xlab = 'sd.alpha.gam', ylab = 'Density',
    frame = FALSE, axes = FALSE, xlim = c(0, 1))
axis(1)  ;  axis(2)

plot(density(out10d$sims.list$sd.alpha.p), col = 'black', lwd = 2,
    main = "SD of OD in p", xlab = 'sd.alpha.p', ylab = 'Density',
    frame = FALSE, axes = FALSE, xlim = c(0, 1))
axis(1)  ;  axis(2)

plot(density(out10d$sims.list$sd.rho), col = 'black', lwd = 2,
    main = "SD of Random immigration", xlab = 'sd.rho', ylab = 'Density',
    frame = FALSE, axes = FALSE, xlim = c(0, 3))
axis(1)  ;  axis(2)
par(op)

# ~~~~~~~~ figure 1.12 ~~~~~~~~~~
# New covariates for prediction on the original scale
elevo <- seq(250, 2750,, 100)
foresto <- seq(0, 100,, 100)
dateo <- seq(115, 205,, 100)
duro <- seq(60, 675,, 100)

# Scale them all identically as we did for the actual covariates in the analyses
elevp <- standardize2match(elevo, dat$elev)
forestp <- standardize2match(foresto, dat$forest)
datep <- standardize2match(dateo, date)
durp <- standardize2match(duro, dur)


# Form predictions on the inverse-link = natural scale
#Expected initial abundance lambda
predlam.elev <- exp(out11$mean$alpha.lam + out11$mean$beta.lam[1] *elevp +
    out11$mean$beta.lam[2] *elevp^2)
predlam.forest <- exp(out11$mean$alpha.lam + out11$mean$beta.lam[3] *forestp)

# Expected population growth rate gamma
predgam.elev <- exp(out11$mean$alpha.gam + out11$mean$beta.gam[1] *elevp +
    out11$mean$beta.gam[2] *elevp^2)
predgam.forest <- exp(out11$mean$alpha.gam + out11$mean$beta.gam[3] *forestp)

# Expected detection probability p
predp.date <- plogis(out11$mean$alpha.p + out11$mean$beta.p[1] *datep +
    out11$mean$beta.p[2] *datep^2)
predp.dur <- plogis(out11$mean$alpha.p + out11$mean$beta.p[3] *durp)

# Plot predictions
op <- par(mfrow = c(3, 2))
plot(elevo, predlam.elev, xlab = "Elevation (m)", ylab = "Expected lambda",
    main = "Initial abundance ~ Elevation", type = 'l', lwd = 3, col = 'blue',
    frame = FALSE, ylim = c(5, 15))
plot(foresto, predlam.forest, xlab = "Forest cover (%)", ylab = "Expected lambda",
    main = "Initial abundance ~ Forest cover", type = 'l', lwd = 3, col = 'blue',
    frame = FALSE, ylim = c(0, 32))
plot(elevo, predgam.elev, xlab = "Elevation (m)", ylab = "Expected gamma",
    main = "Growth rate ~ Elevation", type = 'l', lwd = 3, col = 'blue',
    frame = FALSE, ylim = c(0.7, 1.1))
plot(foresto, predgam.forest, xlab = "Forest cover (%)", ylab = "Expected gamma",
    main = "Growth rate ~ Forest cover", type = 'l', lwd = 3, col = 'blue',
    frame = FALSE, ylim = c(0.7, 1.1))
plot(dateo, predp.date, xlab = "Julian date (1 = Jan 1)", ylab = "Expected p",
    main = "Detection ~ Date", type = 'l', lwd = 3, col = 'blue',
    frame = FALSE, ylim = c(0.2, 0.8))
plot(duro, predp.dur, xlab = "Duration (min)", ylab = "Expected p",
    main = "Detection ~ Duration", type = 'l', lwd = 3, col = 'blue',
    frame = FALSE, ylim = c(0.2, 0.8))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Code for figure 1.13 is at the end of the file "AHM2_01.07.2extra.R"
