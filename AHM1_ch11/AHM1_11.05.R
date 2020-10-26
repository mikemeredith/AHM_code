#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 11. Hierarchical models for communities
# =========================================================================

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14" # location of the "WinBUGS14.exe" application
library(jagsUI)

# ~~~~~~ this section requires the data prepared in section 11.3 ~~~~~~~~~~
source("AHM1_11.03.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 11.5 Community models that ignore species identity
# ==================================================

# 11.5.1 Simple Poisson regression for the observed community size
# ------------------------------------------------------------------------
# Get covariates (from AHMbook::MHB2014) and standardise them

# Quadrat elevation and forest cover
orig.ele <- MHB2014$sites$elev
(mean.ele <- mean(orig.ele, na.rm = TRUE))
(sd.ele <- sd(orig.ele, na.rm = TRUE))
ele <- (orig.ele - mean.ele) / sd.ele
orig.forest <- MHB2014$sites$forest
(mean.forest <- mean(orig.forest, na.rm = TRUE))
(sd.forest <- sd(orig.forest, na.rm = TRUE))
forest <- (orig.forest - mean.forest) / sd.forest

# Average date and duration of survey
orig.mdate <- apply(MHB2014$date, 1, mean, na.rm = TRUE)
(mean.mdate <- mean(orig.mdate[-NAsites]))   # drop unsurved site
(sd.mdate <- sd(orig.mdate[-NAsites]))
mdate <- (orig.mdate - mean.mdate) / sd.mdate
mdate[NAsites] <- 0                 # impute mean for missing

orig.mdur <- apply(MHB2014$dur, 1, mean, na.rm = TRUE)
(mean.mdur <- mean(orig.mdur[-NAsites]))
(sd.mdur <- sd(orig.mdur[-NAsites]))
mdur <- (orig.mdur - mean.mdur) / sd.mdur
mdur[NAsites] <- 0                  # impute mean for missing

# Bundle data and summarize input data for BUGS
str( win.data <- list(C = C, nsite = length(C), ele = ele, forest = forest,
    mdate = mdate, mdur = mdur) )

# Specify model in BUGS language
sink("model1.txt")
cat("
model {

  # Priors
  gamma0 ~ dnorm(0, 0.001)           # Regression intercept
  for(v in 1:6){                     # Loop over regression coef's
    gamma[v] ~ dnorm(0, 0.001)
  }

  # Likelihood for Poisson GLM
  for (i in 1:nsite){
    C[i] ~ dpois(lambda[i])
    log(lambda[i]) <- gamma0 + gamma[1] * ele[i] + gamma[2] * pow(ele[i],2) +
        gamma[3] * forest[i] + gamma[4] * mdate[i] + gamma[5] * pow(mdate[i],2) +
        gamma[6] * mdur[i]
  }
}
",fill = TRUE)
sink()

# Initial values
inits <- function() list(gamma0 = rnorm(1), gamma = rnorm(6))

# Parameters monitored
params <- c("gamma0", "gamma")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART <1 min)
out1 <- bugs(win.data, inits, params, "model1.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)  # ~~~~ for autotesting

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
library(jagsUI)
out1J <- jags(win.data, inits, params, "model1.txt",
  # n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)  # ~~~~~ for testing
traceplot(out1J)    ;    print(out1J, dig = 3)


# Plot posterior distributions for potentially 'ecological' parameters
op <- par(mfrow = c(1,3))
hist(out1J$sims.list$gamma[,1], breaks = 50, col = "grey", main = "",
    xlab = "Slope of elevation (linear)")
hist(out1J$sims.list$gamma[,2], breaks = 50, col = "grey", main = "",
    xlab = "Slope of elevation (squared)")
hist(out1J$sims.list$gamma[,3], breaks = 50, col = "grey", main = "",
    xlab = "Slope of forest")
par(op)

# Get covariate values for prediction
orig.pred.ele <- seq(250, 2750,, 500)   # 500 vals spread between 250 and 2750
p.ele <- (orig.pred.ele - mean.ele) / sd.ele
orig.pred.forest <- seq(1, 100,, 500)
p.forest <- (orig.pred.forest - mean.forest) / sd.forest

# Compute predictions
nsamp <- out1J$mcmc.info$n.samples
pred.ele <- pred.forest <- array(NA, dim = c(500, nsamp))
for(i in 1:nsamp){
   pred.ele[,i] <- exp(out1J$sims.list$gamma0[i] + out1J$sims.list$gamma[i,1] *
      p.ele + out1J$sims.list$gamma[i,2]* p.ele^2)
   pred.forest[,i] <- exp(out1J$sims.list$gamma0[i] +
      out1J$sims.list$gamma[i,3] * p.forest)
}

# Plot posterior mean and a random sample of 100 from posterior of regression
selection <- sample(1:nsamp, 100)
op <- par(mfrow = c(1,3))
matplot(orig.pred.ele, pred.ele[,selection], ylab = "Predicted species count",
    xlab = "Elevation (m a.s.l.)", type = "l", lty = 1, lwd = 1, col = "grey",
    ylim = c(0, 50), frame = FALSE)
lines(orig.pred.ele, apply(pred.ele, 1, mean), lwd = 3, col = "blue")
matplot(orig.pred.forest, pred.forest[,selection], ylab = "Predicted species count",
    xlab = "Forest cover (%)", type = "l", lty = 1, lwd = 1, col = "grey",
    ylim = c(0, 50), frame = FALSE)
lines(orig.pred.forest, apply(pred.forest, 1, mean), lwd = 3, col = "blue")


# Get observed species richness per site and rep and plot
CC <- apply(y, c(1,2), sum, na.rm = TRUE)
CC[CC == 0] <- NA            # 0 means not surveyed
matplot(t(CC), type = 'l', lty = 1, lwd = 2, xlab = "First to third survey",
    ylab = "Number of species detected", frame = FALSE)  # Fig. 11–6 right
par(op)

# Get survey date and survey duration and standardise both
# Survey date (this is Julian date, with day 1 being April 1)
orig.DAT <- MHB2014$date
(mean.date <- mean(orig.DAT, na.rm = TRUE))
(sd.date <- sd(c(orig.DAT), na.rm = TRUE))
DAT <- (orig.DAT - mean.date) / sd.date      # scale
DAT[is.na(DAT)] <- 0                         # impute missings
# Survey duration (in minutes)
orig.DUR <- MHB2014$dur
(mean.dur <- mean(orig.DUR, na.rm = TRUE))
(sd.dur <- sd(c(orig.DUR), na.rm = TRUE))
DUR <- (orig.DUR - mean.dur) / sd.dur        # scale
DUR[is.na(DUR)] <- 0                         # mean impute missings

# Bundle data and summarize
str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele,
    forest = forest, DAT = DAT, DUR = DUR) )


# Specify model in BUGS language
sink("model2.txt")
cat("
model {

  # Priors
  gamma0 ~ dnorm(0, 0.001)
  for(v in 1:6){
    gamma[v] ~ dnorm(0, 0.001)
  }

  # Likelihood for Poisson GLM
  for (i in 1:M){             # Loop over sites
    for(j in 1:J){            # Loop over occasions
      CC[i,j] ~ dpois(lambda[i,j])
      log(lambda[i,j]) <- gamma0 + gamma[1] * ele[i] + gamma[2] * pow(ele[i],2) +
        gamma[3] * forest[i] + gamma[4] * DAT[i,j] + gamma[5] * pow(DAT[i,j],2) +
        gamma[6] * DUR[i,j]
    }
  }
}
",fill = TRUE)
sink()

# Initial values
inits <- function() list(gamma0 = rnorm(1), gamma = rnorm(6))

# Parameters monitored
params <- c("gamma0", "gamma")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART 1.7 min)
out2 <- bugs(win.data, inits, params, "model2.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)  # ~~~~ for autotesting

# Call JAGS from R (ART 0.6 min)
out2J <- jags(win.data, inits, params, "model2.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
traceplot(out2J)   ;   print(out2J, dig = 2)


# 11.5.2 Poisson random effects model for the observed community size
# ------------------------------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele,
    forest = forest, DAT = DAT, DUR = DUR) )

# Specify model in BUGS language
sink("model3.txt")
cat("
model {

  # Priors
  mugamma0 ~ dnorm(0, 0.001)     # Hyperparameters
  taugamma0 <- pow(sd.gamma0,-2)
  sd.gamma0 ~ dunif(0, 10)
  for(v in 1:6){                 # Parameters
    gamma[v] ~ dnorm(0, 0.001)
  }

  # Likelihood for Poisson GLMM
  for (i in 1:M){                # Loop over sites
    gamma0[i] ~ dnorm(mugamma0, taugamma0)     # site intercepts random now
    for(j in 1:J){               # Loop over repeated measurements
      CC[i,j] ~ dpois(lambda[i,j])
      log(lambda[i,j]) <- gamma0[i] + gamma[1]*ele[i] + gamma[2] * pow(ele[i],2) +
        gamma[3] * forest[i] + gamma[4] * DAT[i,j] + gamma[5] * pow(DAT[i,j],2) +
        gamma[6] * DUR[i,j]
    }
  }
}
",fill = TRUE)
sink()

# Initial values
inits <- function() list(gamma0 = rnorm(nrow(CC)), gamma = rnorm(6))

# Parameters monitored
params <- c("mugamma0", "sd.gamma0", "gamma0", "gamma")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART 2.9 min)
out3 <- bugs(win.data, inits, params, "model3.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)  # ~~~~ for autotesting

# Call JAGS from R (ART 0.7 min)
out3J <- jags(win.data, inits, params, "model3.txt",
  # n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)  # ~~~~~ for testing
traceplot(out3J, c('mugamma0', 'sd.gamma0', 'gamma'))   ;    print(out3J, dig = 2)


# 11.5.3 N-mixture model for the observed community size
# -----------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele,
    forest = forest, DAT = DAT, DUR = DUR) )

# Specify model in BUGS language
sink("model4.txt")
cat("
model {

  # Priors
  alpha0 ~ dnorm(0, 0.01)      # Base-line community detection probability
  beta0 ~ dnorm(0, 0.01)       # Base-line community size (number of species)
  for(v in 1:3){
    alpha[v] ~ dnorm(0, 0.01) # Covariate effects on detection
    beta[v] ~ dnorm(0, 0.01)  # Covariate effects on community size
  }

  # Likelihood
  # Ecological model for true community size
  for (i in 1:M){              # Loop over sites
    N[i] ~ dpois(lambda[i])   # Community size
    lambda[i] <- exp(beta0 + beta[1] * ele[i] + beta[2] * pow(ele[i],2) +
        beta[3] * forest[i])

    # Observation model for repeated measurements
    for (j in 1:J){          # Loop over occasions
      CC[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))
      lp[i,j] <- alpha0 + alpha[1] * DAT[i,j] + alpha[2] * pow(DAT[i,j],2) +
           alpha[3] * DUR[i,j]
    # logit(p) = ... causes undefined real result in WinBUGS (but not JAGS)
    }
  }
}
",fill = TRUE)
sink()

# Define function to generate random initial values
Nst <- apply(CC, 1, max, na.rm = TRUE) + 1
Nst[Nst == -Inf] <- max(Nst, na.rm = TRUE)  # Some nonzero val. for unsurv. sites
inits <- function() list(N = Nst, alpha0 = rnorm(1), alpha = rnorm(3),
    beta0 = rnorm(1), beta = rnorm(3))

# Parameters monitored
params <- c("alpha0", "alpha", "beta0", "beta","N")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Run JAGS from R (ART 1.5 min) in parallel, look at traceplots
#   and summarize posteriors
out4 <- jags(win.data, inits, params, "model4.txt",
   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
traceplot(out4, c('alpha0', 'alpha', 'beta0', 'beta'))  ;  print(out4, 3)

print(tmp <- cbind(out3$summary[c(1, 270:272),c(1:3,7)], out4$summary[5:8, c(1:3, 7)]), 3)

plot(orig.ele, out4$summary[9:275, 1], pch = 16, xlab = "Elevation (m)",
    ylab = "Species richness", ylim = c(0, 60), frame = FALSE)
segments(orig.ele, out4$summary[9:275, 3], orig.ele, out4$summary[9:275, 7])
points(orig.ele+20, C)      # elevation jittered

