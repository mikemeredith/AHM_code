#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5 : MODELING METACOMMUNITY DYNAMICS USING DYNAMIC COMMUNITY MODELS
# ==========================================================================
# Code from proofs dated 2020-08-19

# Approximate run time for this script: 15 mins

library(AHMbook)
library(jagsUI)

# 5.3 Fitting the simplest dynamic community models
# =================================================

# Simulate a small data set
set.seed(1)
str(dat <- simDCM(nspec = 50, nsites = 100, nsurveys = 3, nyears = 10,
    mean.psi1 = 0.2, sig.lpsi1 = 2,
    range.mean.phi = c(0.6, 0.6), sig.lphi = 1,
    range.mean.gamma = c(0.1, 0.1), sig.lgamma = 1,
    range.mean.p = c(0.3, 0.3), sig.lp = 2) ) # Not all output shown

# ** Number of species ever occurring: 50
# ** Number of species ever detected: 50
# ** Average number of years of occurrence: 9.96
# ** Average number of years with detection: 9.38

str(dat$y) # This is our data set for input to BUGS
# int [1:100, 1:3, 1:10, 1:50] 1 0 1 0 0 0 1 0 0 0 ...
# - attr(*, "dimnames")=List of 4
# ..$ : chr [1:100] "Site1" "Site2" "Site3" "Site4" ...
# ..$ : chr [1:3] "Survey1" "Survey2" "Survey3"
# ..$ : chr [1:10] "Year1" "Year2" "Year3" "Year4" ...
# ..$ : chr [1:50] "Spec1" "Spec2" "Spec3" "Spec4" ...


# 5.3.1 The fixed-effects DCM
# ---------------------------

# Bundle and summarize data set
y <- dat$y # copy 4D array
str(bdata <- list(y = y, nsite = dim(y)[1], nsurvey = dim(y)[2],
    nyear = dim(y)[3], nspec = dim(y)[4]))
# List of 5
# $ y    : int [1:100, 1:3, 1:10, 1:50] 1 0 1 0 0 0 1 0 0 0 ...
# ..- attr(*, "dimnames")=List of 4
# .. ..$ : chr [1:100] "Site1" "Site2" "Site3" "Site4" ...
# .. ..$ : chr [1:3] "Survey1" "Survey2" "Survey3"
# .. ..$ : chr [1:10] "Year1" "Year2" "Year3" "Year4" ...
# .. ..$ : chr [1:50] "Spec1" "Spec2" "Spec3" "Spec4" ...
# $ nsite  : int 100
# $ nsurvey: int 3
# $ nyear  : int 10
# $ nspec  : int 50

# Specify model in BUGS language
cat(file = "DCM1.txt", "
model {

  # Specify priors
  for(k in 1:nspec){ # Loop over species
    psi1[k] ~ dunif(0, 1)          # Initial occupancy
    phi[k] ~ dunif(0, 1)           # Persistence
    gamma[k] ~ dunif(0, 1)         # Colonization
    p[k] ~ dunif(0, 1)             # Detection
  }

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){ # Loop over sites
    for(k in 1:nspec){ # Loop over species
      # Initial conditions of system
      z[i,1, k] ~ dbern(psi1[k])   # Presence/absence at start of study
      # State transitions
      for (t in 2:nyear){          # Loop over years
        z[i,t,k] ~ dbern(z[i,t-1,k] * phi[k] + (1-z[i,t-1, k]) * gamma[k])
      }
    }
  }

  # Observation model
  for (i in 1:nsite){              # Loop over sites
    for(k in 1:nspec){             # Loop over species
      for (j in 1:nsurvey){        # Loop over surveys
        for (t in 1:nyear){        # Loop over years
          y[i,j,t,k] ~ dbern(z[i,t,k] * p[k])
        }
      }
    }
  }

  # Derived parameters: Number of occupied sites and population occupancy
  for(k in 1:nspec){ # Loop over species
    n.occ[1, k] <- sum(z[,1,k])    # Number of occupied sites
    psi[1, k] <- psi1[k]           # Population occupancy
    for (t in 2:nyear){            # Loop over years
      n.occ[t, k] <- sum(z[,t,k])
      psi[t, k] <- psi[t-1, k] * phi[k] + (1-psi[t-1, k]) * gamma[k]
    }
  }
}
")

# Initial values
zst <- apply(y, c(1,3,4), max) # Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored (could also add "z")
params <- c("psi1", "phi", "gamma", "p", "n.occ", "psi")

# MCMC settings
# na <- 1000 ; ni <- 6000 ; nt <- 6 ; nb <- 3000 ; nc <- 3
na <- 1000 ; ni <- 600 ; nt <- 1 ; nb <- 300 ; nc <- 3  # ~~~~ for testing, 7 mins

# Call JAGS (ART 196 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "DCM1.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
summary(out1) ; jags.View(out1) ; print(out1, 3) # not shown


# 5.3.2 The random-effects DCM
# ----------------------------

# Specify model in BUGS language
cat(file = "DCM2.txt", "
model {

  # Specify priors: Declare species-level effects as random
  for(k in 1:nspec){               # Loop over species
    logit(psi1[k]) <- lpsi1[k]     # Initial occupancy
    lpsi1[k] ~ dnorm(mu.lpsi1, tau.lpsi1)
    logit(phi[k]) <- lphi[k]       # Persistence
    lphi[k] ~ dnorm(mu.lphi, tau.lphi)
    logit(gamma[k]) <- lgamma[k]   # Colonization
    lgamma[k] ~ dnorm(mu.lgamma, tau.lgamma)
    logit(p[k]) <- lp[k]           # Detection
    lp[k] ~ dnorm(mu.lp, tau.lp)
  }

  # Specify hyperpriors: Priors for the hyperparameters
  mu.lpsi1 <- logit(mean.psi1)     # Initial occupancy
  mean.psi1 ~ dunif(0, 1)
  tau.lpsi1 <- pow(sd.lpsi1, -2)
  sd.lpsi1 ~ dunif(0, 10)
  mu.lphi <- logit(mean.phi)       # Persistence
  mean.phi ~ dunif(0, 1)
  tau.lphi <- pow(sd.lphi, -2)
  sd.lphi ~ dunif(0, 10)
  mu.lgamma <- logit(mean.gamma)   # Colonization
  mean.gamma ~ dunif(0, 1)
  tau.lgamma <- pow(sd.lgamma, -2)
  sd.lgamma ~ dunif(0, 10)
  mu.lp <- logit(mean.p)           # Detection
  mean.p ~ dunif(0, 1)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 10)

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){              # Loop over sites
    for(k in 1:nspec){             # Loop over species
      # Initial conditions of system
      z[i,1, k] ~ dbern(psi1[k])   # Presence/absence at start of study
      # State transitions
      for (t in 2:nyear){          # Loop over years
        z[i,t,k] ~ dbern(z[i,t-1,k] * phi[k] + (1-z[i,t-1, k]) * gamma[k])
      }
    }
  }

  # Observation model
  for (i in 1:nsite){              # Loop over sites
    for(k in 1:nspec){             # Loop over species
      for (j in 1:nsurvey){        # Loop over surveys
        for (t in 1:nyear){        # Loop over years
          y[i,j,t,k] ~ dbern(z[i,t,k] * p[k])
        }
      }
    }
  }

  # Derived parameters: Number of occupied sites and population occupancy
  for(k in 1:nspec){               # Loop over species
    n.occ[1, k] <- sum(z[,1,k])    # Number of occupied sites
    psi[1, k] <- psi1[k]           # Population occupancy
    for (t in 2:nyear){            # Loop over years
      n.occ[t, k] <- sum(z[,t,k])
      psi[t, k] <- psi[t-1, k] * phi[k] + (1-psi[t-1, k]) * gamma[k]
    }
  }
}
")

# Parameters monitored (could also add "z")
params <- c("mu.lpsi1", "sd.lpsi1", "mu.lphi", "sd.lphi", "mu.lgamma",
    "sd.lgamma", "mu.lp", "sd.lp", "psi1", "phi", "gamma", "p", "n.occ", "psi")

# Call JAGS (ART 188 min), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "DCM2.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out2)
summary(out2) ; jags.View(out2) ; print(out2, 3) # not shown


# ~~~~~~~~~ Additional plotting code for figures 5.1 and 5.2 ~~~~~~~~~~~

# Figure 5.1
# ''''''''''
# Compare estimates of hyperparameters with truth
op <- par(mfrow = c(4,2))
plot(density(plogis(out2$sims.list$mu.lpsi1)), lwd = 2, col = 'gray',
    main = 'Community mean of lpsi1', xlab = 'mu.lpsi1', ylab = 'Density',
    xlim = c(0, 0.5), frame = F)
abline(v = dat$mean.psi1, col = 'black', lwd = 2)
plot(density(out2$sims.list$sd.lpsi1), lwd = 2, col = 'gray',
    main = 'Community SD of lpsi1', xlab = 'sd.lpsi1', ylab = 'Density',
    xlim = c(0, 4), frame = F)
abline(v = dat$sig.lpsi1, col = 'black', lwd = 2)

plot(density(plogis(out2$sims.list$mu.lphi)), lwd = 2, col = 'gray',
    main = 'Community mean of lpsi1', xlab = 'mu.lphi', ylab = 'Density',
    xlim = c(0.4, 1), frame = F)
abline(v = dat$range.mean.phi[1], col = 'black', lwd = 2)
plot(density(out2$sims.list$sd.lphi), lwd = 2, col = 'gray',
    main = 'Community SD of lphi', xlab = 'sd.lphi', ylab = 'Density',
    xlim = c(0, 2), frame = F)
abline(v = dat$sig.lphi, col = 'black', lwd = 2)

plot(density(plogis(out2$sims.list$mu.lgamma)), lwd = 2, col = 'gray',
    main = 'Community mean of lgamma', xlab = 'mu.lgamma', ylab = 'Density',
    xlim = c(0, 0.3), frame = F)
abline(v = dat$range.mean.gamma[1], col = 'black', lwd = 2)
plot(density(out2$sims.list$sd.lgamma), lwd = 2, col = 'gray',
    main = 'Community SD of lgamma', xlab = 'sd.lgamma', ylab = 'Density',
    xlim = c(0, 2), frame = F)
abline(v = dat$sig.lgamma, col = 'black', lwd = 2)

plot(density(plogis(out2$sims.list$mu.lp)), lwd = 2, col = 'gray',
    main = 'Community mean of lp', xlab = 'mu.lp', ylab = 'Density',
    xlim = c(0, 1), frame = F)
abline(v = dat$range.mean.p[1], col = 'black', lwd = 2)
plot(density(out2$sims.list$sd.lp), lwd = 2, col = 'gray',
    main = 'Community SD of lp', xlab = 'sd.lp', ylab = 'Density',
    xlim = c(0, 4), frame = F)
abline(v = dat$sig.lp, col = 'black', lwd = 2)
par(op)

# Figure 5.2
# ''''''''''
psi1.true <- plogis(dat$beta0.lpsi)
phi.true <- plogis(dat$beta0.lphi)[,1]
gamma.true <- plogis(dat$beta0.lgamma)[,1]
p.true <- plogis(dat$beta0.lp)[,1]
nocc.true <- dat$n.occ

op <- par(mfrow = c(4, 2))
lim <- c(0,1)
# Initial occupancy (psi1)
plot(psi1.true, out1$mean$psi1, main = 'psi1 (fixed effects)', xlim = lim,
    ylim = lim, pch = 16, frame = FALSE, xlab = 'True value', ylab = 'Estimate')
segments(psi1.true, out1$q2.5$psi1, psi1.true, out1$q97.5$psi1)
abline(0,1, lwd = 2, col = 'red')
abline(lm(out1$mean$psi1 ~ psi1.true), col = 'blue', lwd = 2, lty = 2)

plot(psi1.true, out2$mean$psi1, main = 'psi1 (random effects)', xlim = lim,
    ylim = lim, pch = 16, frame = FALSE, xlab = 'True value', ylab = 'Estimate')
segments(psi1.true, out2$q2.5$psi1, psi1.true, out2$q97.5$psi1)
abline(0,1, lwd = 2, col = 'red')
abline(lm(out2$mean$psi1 ~ psi1.true), col = 'blue', lwd = 2, lty = 2)

# Persistence
plot(phi.true, out1$mean$phi, main = 'phi (fixed effects)', xlim = lim,
    ylim = lim, pch = 16, frame = FALSE, xlab = 'True value', ylab = 'Estimate')
segments(phi.true, out1$q2.5$phi, phi.true, out1$q97.5$phi)
abline(0,1, lwd = 2, col = 'red')
abline(lm(out1$mean$phi ~ phi.true), col = 'blue', lwd = 2, lty = 2)

plot(phi.true, out2$mean$phi, main = 'phi (random effects)', xlim = lim,
    ylim = lim, pch = 16, frame = FALSE, xlab = 'True value', ylab = 'Estimate')
segments(phi.true, out2$q2.5$phi, phi.true, out2$q97.5$phi)
abline(0,1, lwd = 2, col = 'red')
abline(lm(out2$mean$phi ~ phi.true), col = 'blue', lwd = 2, lty = 2)

# Colonization
plot(gamma.true, out1$mean$gamma, main = 'gamma (fixed effects)', xlim = lim,
    ylim = lim, pch = 16, frame = FALSE, xlab = 'True value', ylab = 'Estimate')
segments(gamma.true, out1$q2.5$gamma, gamma.true, out1$q97.5$gamma)
abline(0,1, lwd = 2, col = 'red')
abline(lm(out1$mean$gamma ~ gamma.true), col = 'blue', lwd = 2, lty = 2)

plot(gamma.true, out2$mean$gamma, main = 'gamma (random effects)', xlim = lim,
    ylim = lim, pch = 16, frame = FALSE, xlab = 'True value', ylab = 'Estimate')
segments(gamma.true, out2$q2.5$gamma, gamma.true, out2$q97.5$gamma)
abline(0,1, lwd = 2, col = 'red')
abline(lm(out2$mean$gamma ~ gamma.true), col = 'blue', lwd = 2, lty = 2)

# Number of occupied sites (n.occ)
lim <- c(0, 100)
plot(nocc.true, out1$mean$n.occ, main = 'n.occ (fixed effects)', xlim = lim,
    ylim = lim, pch = 16, frame = FALSE, xlab = 'True value', ylab = 'Estimate')
segments(nocc.true, out1$q2.5$n.occ, nocc.true, out1$q97.5$n.occ)
abline(0,1, lwd = 2, col = 'red')
abline(lm(c(out1$mean$n.occ) ~ c(nocc.true)), col = 'blue', lwd = 2, lty = 2)

plot(nocc.true, out2$mean$n.occ, main = 'n.occ (random effects)', xlim = lim,
    ylim = lim, pch = 16, frame = FALSE, xlab = 'True value', ylab = 'Estimate')
segments(nocc.true, out2$q2.5$n.occ, nocc.true, out2$q97.5$n.occ)
abline(0,1, lwd = 2, col = 'red')
abline(lm(c(out2$mean$n.occ) ~ c(nocc.true)), col = 'blue', lwd = 2, lty = 2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~ code to produce the table shown ~~~~~~~
RMSE <- array(NA, dim = c(4, 3), dimnames = list(c('psi1', 'phi', 'gamma', 'n.occ'),
    c('Fixed-eff. (DCM1)', 'Random-eff. (DCM2)', 'Ratio DCM2/DCM1')))
RMSE[1, 1] <- sqrt(mean(out1$mean$psi1-psi1.true)^2)
RMSE[1, 2] <- sqrt(mean(out2$mean$psi1-psi1.true)^2)
RMSE[2, 1] <- sqrt(mean(out1$mean$phi-phi.true)^2)
RMSE[2, 2] <- sqrt(mean(out2$mean$phi-phi.true)^2)
RMSE[3, 1] <- sqrt(mean(out1$mean$gamma-gamma.true)^2)
RMSE[3, 2] <- sqrt(mean(out2$mean$gamma-gamma.true)^2)
RMSE[4, 1] <- sqrt(mean(out1$mean$n.occ-nocc.true)^2)
RMSE[4, 2] <- sqrt(mean(out2$mean$n.occ-nocc.true)^2)
RMSE[,3] <- RMSE[,2]/RMSE[,1]
print(RMSE, 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Table of RMSE for both models and RMSE ratio
#       Fixed-eff. (DCM1) Random-eff. (DCM2) Ratio DCM2/DCM1
# psi1             0.0504              0.027           0.529
# phi              0.0079              0.020           2.488
# gamma            0.0324              0.002           0.062
# n.occ            2.6752              1.067           0.399
