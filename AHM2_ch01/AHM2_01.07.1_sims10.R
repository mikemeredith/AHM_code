#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code based on MS "Report on 2nd set of sims: focused comparisons just for the
#  DM model", Marc Kery, 1 August 2018 and later, largely refactored by Mike.

# 1.7 “Demographic” state-space models for inference about relative abundance
# 1.7.1 Simulation assessment of a demographic state-space model
# --------------------------------------------------------------

# This file has code to run the simulations for case 10 described on
#    p.34 and to plot figure 1.10.


library(AHMbook)
library(jagsUI)

# Specify model 6 in BUGS language
# ================================

# MODEL 6 : Extended Markov model with site-level OD in lambda, and
#           site-by-interval-specific OD in both gamma and p
# ------------------------------------------------------------------
cat(file = "simModel6.txt","
model {
  # Priors
  for(i in 1:M){
    alpha.lam[i] ~ dnorm(mu.alpha.lam, tau.alpha.lam)
  }
  mu.alpha.lam ~ dnorm(0, 0.1) # Expected initial abundance with X = 0
  mean.lambda <- exp(mu.alpha.lam)
  beta.lam ~ dnorm(0, 0.1)     # Coefficient of covariate on lambda
  tau.alpha.lam <- pow(sd.alpha.lam, -2)
  sd.alpha.lam ~ dunif(0, 10)  # Overdispersion in lambda

  for(i in 1:M){
    for(t in 1:(T-1)){
      alpha.gam[i,t] ~ dnorm(mu.alpha.gam, tau.alpha.gam)
    }
  }
  mu.alpha.gam ~ dunif(-0.2, 0.2) # Population growth rate with X = 0
  mean.gamma <- exp(mu.alpha.gam)
  beta.gam ~ dnorm(0, 0.1)     # Coefficient of covariate on gamma
  tau.alpha.gam <- pow(sd.alpha.gam, -2)
  sd.alpha.gam ~ dunif(0, 10)  # Overdispersion in gamma

  for(i in 1:M){
    for(t in 1:T){
      alpha.p[i,t] ~ dnorm(mu.alpha.p, tau.alpha.p)
    }
  }
  mu.alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1)         # Detection probability
  beta.p ~ dnorm(0, 0.1)       # Coefficient of covariate on p
  tau.alpha.p <- pow(sd.alpha.p, -2)
  sd.alpha.p ~ dunif(0, 10)    # Overdispersion in p

  for(t in 1:(T-1)){           # half-normal prior for random immigration
    log(eps[t]) <- logeps[t]
    logeps[t] ~ dnorm(0, tau.eps)
  }
  tau.eps <- pow(sd.eps, -2)
  sd.eps ~ dnorm(0, 5)I(0,)
  # curve(dnorm(x, 0, sqrt(1 / 5)), 0, 3)

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
      loggam[i, t-1] <- alpha.gam[i, t-1] + beta.gam * Xgam[i, t-1]
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

# Parameters monitored
params6 <- c("mean.lambda", "mean.gamma", "mean.p", "mu.alpha.lam",
    "mu.alpha.gam", "mu.alpha.p", "sd.alpha.lam", "sd.alpha.gam",
    "sd.alpha.p", "sd.eps", "beta.lam", "beta.gam", "beta.p", "sumN",
    "gammaX", "eps", "alpha.lam", "alpha.gam", "alpha.p", "N")

# Prepare stuff we need
# =====================
# simrep <- 100
simrep <- 3 # ~~~ for testing
# simrep <- 30 # ~~~ for testing
# Each rep takes about 20 mins

# For truth, we need the N's, sumN's and gammaX's for each simrep
trueN <- array(NA, dim=c(100, 10, simrep))
truesumN <- array(NA, dim=c(10, simrep))
truegammaX <- array(NA, dim=c(9, simrep))

# We need the estimates of those
estiN <- array(NA, dim=c(100, 10, simrep))
estisumN <- array(NA, dim=c(10, simrep))
estigammaX <- array(NA, dim=c(9, simrep))

# and some other stuff
estiEtc <- array(NA, dim=c(7, simrep))
dimnames(estiEtc)[[1]] <- c("beta.lam", "beta.gam", "beta.p",
    "sd.alpha.lam", "sd.alpha.gam", "sd.alpha.p", "sd.eps")

seeds <- c(10, 13, 18) #, 19, 29) # ~~~ known to converge, which we need for the plot.
# Other seeds are possible, provided
length(seeds) >= simrep

# Run the simulations for scenario 6
# ==================================

# Scenario 6: Extended Markov simulation model with site-level heterogeneity
#             in lambda and site-by-occasion-level heterogeneity in both gamma
#             and p (case 10)

system.time(
for(i in 1:simrep){
  cat(paste("\n\n***** Case 10 : Simrep", i, "of", simrep, "*****\n\n"))
  flush.console()

  # Scenario 6: Extended Markov dynamics (rescue effect or random
  #    immigration) and heterogeneity in lambda, gamma and p
  set.seed(seeds[i])
  betas <- 0.5
  sd.log.lam <- 1 # Value of overdispersion in lambda
  sd.log.gamma.survey <- 0.5 # Value of overdispersion in gamma
  sd.logit.p.survey <- 1 # Value of overdispersion in p
  data <- simPOP(mean.lam = 3, beta.lam = betas, sd.log.lam = sd.log.lam,
      mean.gamma = 1.0, beta.gamma = betas,
      sd.log.gamma.survey = sd.log.gamma.survey,
      sd.rho = 0.2, mean.p = 0.6, beta.p = betas,
      sd.logit.p.survey = sd.logit.p.survey, show.plot=FALSE)
  trueN[, , i] <- data$N
  truesumN[, i] <- data$sumN
  truegammaX[, i] <- data$gammaX

  # Bundle data
  bdata <- list(C = data$C, M = nrow(data$C),  T = ncol(data$C),
      Xlam = data$Xsite1, Xgam = data$Xsiteyear1, Xp = data$Xsiteyear2)

  # Initial values
  # Gets inits for N
  Nst <- data$C + 10
  inits <- function(){list(N = Nst)}


  # MCMC settings
  na <- 20000  ;  ni <- 1e5  ;  nt <- 50   ;  nb <- 5e4  ;  nc <- 3
  # Some data sets need more iterations, but it's more efficient to run more
  #   simulations (increase 'simrep') and discard those that don't converge.

  # Call JAGS from R, check convergence and summarize marginal posteriors
  out <- try(jags(bdata, inits, params6, "simModel6.txt", n.adapt = na,
      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC = FALSE,
      parallel = TRUE, verbose = FALSE), silent=TRUE)
  # Check for errors and convergence (of first 20 nodes), save results
  if(!inherits(out, "try-error") &&
      !any(out$summary[1:20, 'Rhat'] > 1.1)) {
    estiN[, , i] <- out$mean$N
    estisumN[, i] <- out$mean$sumN
    estigammaX[, i] <- out$mean$gammaX
    estiEtc[, i] <- with(out$mean, c(beta.lam, beta.gam, beta.p,
        sd.alpha.lam, sd.alpha.gam, sd.alpha.p, sd.eps))
  }
} )  # 3 took 1 hr
save(trueN, truesumN, truegammaX, estiN, estisumN, estigammaX, estiEtc,
    # file="sims_case10.RData")
    file="sims_case10_test.RData")  # ~~~ for testing
###########################

estiEtc
estisumN
truesumN


# Figure 1.10
# ===========

op <- par(mfrow=3:4)

plot(trueN[,1,], estiN[,1,], pch=16,
    ylab="", xlab="True N", main="N in Year 1")
abline(0, 1, lwd=2, col='red')
abline(lm(c(estiN[,1,]) ~ c(trueN[,1,])),
    lwd=2, lty=2, col='blue')

plot(trueN[,5,], estiN[,5,], pch=16,
    ylab="", xlab="True N", main="N in Year 5")
abline(0, 1, lwd=2, col='red')
abline(lm(c(estiN[,5,]) ~ c(trueN[,5,])),
    lwd=2, lty=2, col='blue')

plot(trueN[,10,], estiN[,10,], pch=16,
    ylab="", xlab="True N", main="N in Year 10")
abline(0, 1, lwd=2, col='red')
abline(lm(c(estiN[,10,]) ~ c(trueN[,10,])),
    lwd=2, lty=2, col='blue')

plot(truesumN, estisumN, pch=16,
    ylab="", xlab="True sum(N)", main="Total N (all sites)")
abline(0, 1, lwd=2, col='red')
abline(lm(c(estisumN) ~ c(truesumN)),
    lwd=2, lty=2, col='blue')

plot(truegammaX, estigammaX, pch=16,
    ylab="", xlab="True growth rate", main="Growth rate")
abline(0, 1, lwd=2, col='red')
abline(lm(c(estigammaX) ~ c(truegammaX)),
    lwd=2, lty=2, col='blue')

tmp <- c("lambda", "gamma", "p")
title <- c(paste("covariate in", tmp), paste0("SD of OD(", tmp, ")"),
    "SD of Ramdom Imm.")
data_gen <- c(0.5, 0.5, 0.5, 1, 0.5, 1, 0.2)
for(i in 1:7) {
  plot(density(estiEtc[i,], na.rm=TRUE), col = 'black',
      main = title[i],
      xlab = '', ylab = '', lwd = 2, frame = FALSE)
  abline(v = data_gen[i], col = "red", lwd = 2)
  abline(v = mean(estiEtc[i, ], na.rm=TRUE),
      col = "blue", lwd = 2, lty = 2)
}
par(op)

