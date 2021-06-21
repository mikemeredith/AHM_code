#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6 : MULTISTATE OCCUPANCY MODELS
# =======================================
# Code from proofs dated 2020-08-19

# Approximate code execution time for this script: 1 hr
# With full number of iterations: 12 hrs

library(jagsUI)
library(HDInterval)

# ~~~~~~~~~~~ need data preparation from 6.4.1 ~~~~~~~~~~~~~~~~~~
source("AHM2_06.4.1.R")
# ~~~~~ and this from 6.4.3 ~~~~~~~~~~~~~~~~~
date.scaled <- (date - 90) / 30.5
date.scaled[is.na(date.scaled)] <- 0
region <- as.numeric(dat$sites$region)
elev <- dat$sites$elev
forest <- dat$sites$forest
elev.scaled <- standardize(elev)
forest.scaled <- standardize(forest)
# ~~~~ and this from 6.4.4 ~~~~~~~~~~~~~~~~~~~
V <- yms
V[V > 1] <- 1
V[is.na(V)] <- 0
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 6.4 Case study: Swiss eagle owls
# ================================

# 6.4.5 Fitting a static multiseason model with trends and covariates
# -------------------------------------------------------------------
table(REGION <- as.numeric(dat$sites$region))
REGION[REGION > 4] <- 0
REGION[REGION > 0] <- 1
REGION <- REGION+1 ; table(REGION) # 170 Alps, 104 outside of Alps

year <- (1:10)-5.5                       # Create centered year covariate
str(bdata <- list(y = yms, nsites = dim(yms)[1], nsurveys = nsurveys,
    nyears = dim(yms)[3], REGION = REGION, year = year, elev = elev.scaled,
    forest = forest.scaled, date = date.scaled, V = V))
# List of 10
# $ y       : num [1:274, 1:20, 1:10] NA NA 1 NA NA NA NA 2 2 2 ...
# $ nsites  : int 274
# $ nsurveys: num [1:274, 1:10] 1 1 3 1 1 1 1 20 12 20 ...
# $ nyears  : int 10
# $ REGION  : num [1:274] 2 1 1 2 2 2 2 1 2 2 ...
# $ year    : num [1:10] -4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 4.5
# $ elev    : num [1:274] -0.84 0.115 -0.554 -0.267 1.07 ...
# $ forest  : num [1:274] -0.647 0.862 -0.459 -0.459 -1.213 ...
# $ date    : num [1:274, 1:20, 1:10] 0 0 -1.9 0 0 ...
# $ V       : num [1:274, 1:20, 1:10] 0 0 1 0 0 0 0 1 1 1 ...

# Specify model in BUGS language
# ~~~ missing code inserted from MS dated 2019-03-26 ~~~~~~~~~~~~~~~
cat(file = "staticMS1.txt", "
model{

  # Priors and linear models
  # (1) State process
  # Linear models for annual psi (prob site occupied) and r (prob occupied by pair, given occupied)
  for (i in 1:nsites){
    for (t in 1:nyears){
      logit(psi[i,t]) <- alpha.lpsi[REGION[i]] + trend.lpsi[REGION[i]] *
          year[t] + beta.lpsi[1] * elev[i] + beta.lpsi[2] * forest[i]
      logit(r[i,t]) <- alpha.lr[REGION[i]] + trend.lr[REGION[i]] * year[t] +
          beta.lr[1] * elev[i] + beta.lr[2] * forest[i]
    }
  }

  # Priors for parameters in the linear models of psi and r
  for (reg in 1:2){    # Loop over two levels of REGION
    # Intercepts
    alpha.lpsi[reg] <- logit(mean.psi[reg])
    mean.psi[reg] ~ dunif(0, 1)
    alpha.lr[reg] <- logit(mean.r[reg])
    mean.r[reg] ~ dunif(0, 1)
    # Occupancy trends
    trend.lpsi[reg] ~ dnorm(0, 0.01)
    trend.lr[reg] ~ dnorm(0, 0.01)
  }
  # Coefficients of 2 covariates
  for(k in 1:2){
    beta.lpsi[k] ~ dnorm(0, 0.1)
    beta.lr[k] ~ dnorm(0, 0.1)
  }

  # (2) Observation process
  # Linear models in observation process
  for (i in 1:nsites){
    for(t in 1:nyears){
      for(j in 1:nsurveys[i,t]){
        # Observation model for sites in occupied state 1 (= single bird)
        logit(p2[i,j,t]) <- alpha.lp2[REGION[i]] + beta.lp2[1, REGION[i]] *
            date[i,j,t] + beta.lp2[2, REGION[i]] * pow(date[i,j,t],2)

        # Observation model for sites in occupied state 2 (= pairs)
        # Specify linear models
        mlogit.p3[2,i,j,t] <- alpha.lp32[REGION[i]] + beta.lp32[1, REGION[i]] *
            date[i,j,t] + beta.lp32[2, REGION[i]] * pow(date[i,j,t],2)
        mlogit.p3[3,i,j,t] <- alpha.lp33[REGION[i]] + beta.lp33[1, REGION[i]] *
            date[i,j,t] + beta.lp33[2, REGION[i]] * pow(date[i,j,t],2)
      }
    }
  }

  # Priors for parameters in the linear models of p2, p32 and p33
  # Intercepts
  for (reg in 1:2){                      # Loop over levels of REGION
    alpha.lp2[reg] <- logit(mean.alpha.p2[reg])
    mean.alpha.p2[reg] ~ dunif(0, 1)
    alpha.lp32[reg] <- logit(mean.alpha.p32[reg])
    mean.alpha.p32[reg] ~ dunif(0, 1)
    alpha.lp33[reg] <- logit(mean.alpha.p33[reg])
    mean.alpha.p33[reg] ~ dunif(0, 1)
  }

  # Coefficients of survey date (linear and squared)
  for(k in 1:2){
    for(reg in 1:2){
      beta.lp2[k, reg] ~ dnorm(0, 0.1)
      beta.lp32[k, reg] ~ dnorm(0, 0.1)
      beta.lp33[k, reg] ~ dnorm(0, 0.1)
    }
  }

  # Implement Multinomial logit link for p3[2:3]
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

  # Definition of state vector (Omega) and observation matrix (Theta)
  # State vector (Omega)
  for (i in 1:nsites){
    for (t in 1:nyears){
      Omega[i,t,1] <- 1 - psi[i,t]           # Prob. of non-occupation
      Omega[i,t,2] <- psi[i,t] * (1-r[i,t])  # Prob. of occ. by single bird
      Omega[i,t,3] <- psi[i,t] * r[i,t]      # Prob. of occ. by pair
    }
  }

  # Observation matrix (Theta)
  # Order of indices: true state, observed state, site, occasion, year
  for(i in 1:nsites){
    for (t in 1:nyears){
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

  # State-space likelihood
  # State equation: model of true states (z)
  for (i in 1:nsites){
    for (t in 1:nyears){
      z[i,t] ~ dcat(Omega[i,t,])
    }
  }

  # Observation equation: model for observed multistate detections
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        y[i,j,t] ~ dcat(Theta[z[i, t], ,i,j,t])
      }
    }
  }

  # Derived quantities
  # Number of sites in each state per year
  for (t in 1:nyears){
    for (i in 1:nsites){
      state1[i,t] <- equals(z[i,t], 1) # Indicator for site in state 1
      state2[i,t] <- equals(z[i,t], 2) # ... state 2
      state3[i,t] <- equals(z[i,t], 3) # ... state 3
      stateOcc[i,t] <- max(state2[i,t], state3[i,t]) # ... occupied site
    }
    n.occ[t,1] <- sum(state1[,t]) # Number of unoccupied sites
    n.occ[t,2] <- sum(state2[,t]) # Number of sites with single birds
    n.occ[t,3] <- sum(state3[,t]) # Number of sites with pairs
    n.occ.total[t] <- sum(stateOcc[,t]) # All occupied
  }
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(bdata$nsites, bdata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("mean.psi", "alpha.lpsi", "trend.lpsi", "beta.lpsi", "mean.r",
    "alpha.lr", "trend.lr", "beta.lr", "alpha.lp2", "beta.lp2", "alpha.lp32",
    "alpha.lp33", "beta.lp32", "beta.lp33", "n.occ", "n.occ.total", "z")

# MCMC settings
# na <- 1000  ;  ni <- 20000   ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3
na <- 1000  ;  ni <- 2000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3  # ~~~ for testing, 16 mins

# Call JAGS from R, check convergence and summarize posteriors
out4 <- jags(bdata, inits, params, "staticMS1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out4)
summary(out4)  ;  jags.View(out4)

# The PS-variant of the same model
# ''''''''''''''''''''''''''''''''

# Specify model in BUGS language
# NOTE: This model now does accommodate preferential sampling
cat(file = "staticMS2.txt", "
model{
  ### (A) The 'biological model'
  # ----------------------------
  # Priors and linear models
  # (1) State process
  # Linear models for annual psi (prob site occupied) and r (prob occupied by pair, given occupied)
  for (i in 1:nsites){
    for (t in 1:nyears){
      logit(psi[i,t]) <- alpha.lpsi[REGION[i]] + trend.lpsi[REGION[i]] * year[t] +
          beta.lpsi[1] * elev[i] + beta.lpsi[2] * forest[i]
      logit(r[i,t]) <- alpha.lr[REGION[i]] + trend.lr[REGION[i]] * year[t] +
          beta.lr[1] * elev[i] + beta.lr[2] * forest[i]
    }
  }

  # Priors for parameters in the linear models of psi and r
  for (reg in 1:2){    # Loop over levels of REGION
    # Intercepts
    alpha.lpsi[reg] <- logit(mean.psi[reg])
    mean.psi[reg] ~ dunif(0, 1)
    alpha.lr[reg] <- logit(mean.r[reg])
    mean.r[reg] ~ dunif(0, 1)
    # Occupancy trends
    trend.lpsi[reg] ~ dnorm(0, 0.01)
    trend.lr[reg] ~ dnorm(0, 0.01)
  }
  # Coefficients of 2 covariates
  for(k in 1:2){
    beta.lpsi[k] ~ dnorm(0, 0.1)
    beta.lr[k] ~ dnorm(0, 0.1)
  }

  # (2) Observation process
  # Linear models in observation process
  for (i in 1:nsites){
    for(t in 1:nyears){
      for(j in 1:nsurveys[i,t]){
        # Observation model for sites in occupied state 1 (= single bird)
        logit(p2[i,j,t]) <- alpha.lp2[REGION[i]] + beta.lp2[1, REGION[i]] * date[i,j,t] +
            beta.lp2[2, REGION[i]] * pow(date[i,j,t],2)

        # Observation model for sites in occupied state 2 (= pairs)
        # Specify linear models
        mlogit.p3[2,i,j,t] <- alpha.lp32[REGION[i]] + beta.lp32[1, REGION[i]] * date[i,j,t] +
            beta.lp32[2, REGION[i]] * pow(date[i,j,t],2)
        mlogit.p3[3,i,j,t] <- alpha.lp33[REGION[i]] + beta.lp33[1, REGION[i]] * date[i,j,t] +
            beta.lp33[2, REGION[i]] * pow(date[i,j,t],2)
      }
    }
  }

  # Priors for parameters in the linear models of p2, p32 and p33
  # Intercepts
  for (reg in 1:2){    # Loop over levels of REGION
    alpha.lp2[reg] <- logit(mean.alpha.p2[reg])
    mean.alpha.p2[reg] ~ dunif(0, 1)
    alpha.lp32[reg] <- logit(mean.alpha.p32[reg])
    mean.alpha.p32[reg] ~ dunif(0, 1)
    alpha.lp33[reg] <- logit(mean.alpha.p33[reg])
    mean.alpha.p33[reg] ~ dunif(0, 1)
  }

  # Coefficients of survey date (linear and squared)
  for(k in 1:2){
    for(reg in 1:2){
      beta.lp2[k, reg] ~ dnorm(0, 0.1)
      beta.lp32[k, reg] ~ dnorm(0, 0.1)
      beta.lp33[k, reg] ~ dnorm(0, 0.1)
    }
  }

  # Implement Multinomial logit link for p3[2:3]
  for (i in 1:nsites){
    for (t in 1:nyears){
      for(j in 1:nsurveys[i,t]){
        p3[2,i,j,t] <- exp(mlogit.p3[2,i,j,t]) /
            (1 + exp(mlogit.p3[2,i,j,t]) + exp(mlogit.p3[3,i,j,t]))
        p3[3,i,j,t] <- exp(mlogit.p3[3,i,j,t]) /
            (1 + exp(mlogit.p3[2,i,j,t]) + exp(mlogit.p3[3,i,j,t]))
      }
    }
  }

  # Definition of state vector (Omega) and observation matrix (Theta)
  # State vector (Omega)
  for (i in 1:nsites){
    for (t in 1:nyears){
      Omega[i,t,1] <- 1 - psi[i,t]           # Prob. of non-occupation
      Omega[i,t,2] <- psi[i,t] * (1-r[i,t])  # Prob. of occ. by single bird
      Omega[i,t,3] <- psi[i,t] * r[i,t]      # Prob. of occ. by pair
    }
  }

  # Observation matrix (Theta)
  # Order of indices: true state, observed state, site, occasion, year
  for(i in 1:nsites){
    for (t in 1:nyears){
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

  # State-space likelihood
  # State equation: model of true states (z)
  for (i in 1:nsites){
    for (t in 1:nyears){
      z[i,t] ~ dcat(Omega[i,t,])
    }
  }

  # Observation equation: model for observed multistate detections
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        y[i,j,t] ~ dcat(Theta[z[i, t], ,i,j,t])
      }
    }
  }

  # Derived quantities
  # Number of sites in each state per year
  for (t in 1:nyears){
    for (i in 1:nsites){
      state1[i,t] <- equals(z[i,t], 1) # Indicator for site in state 1
      state2[i,t] <- equals(z[i,t], 2) # ... state 2
      state3[i,t] <- equals(z[i,t], 3) # ... state 3
      stateOcc[i,t] <- max(state2[i,t], state3[i,t]) # occupied site
    }
    n.occ[t,1] <- sum(state1[,t]) # Number of unoccupied sites
    n.occ[t,2] <- sum(state2[,t]) # Number of sites with single birds
    n.occ[t,3] <- sum(state3[,t]) # Number of sites with pairs
    n.occ.total[t] <- sum(stateOcc[,t]) # Number of occupied sites
  }

  ### (B) Model for site visitation
  # -------------------------------
  # Priors
  alpha.pv <- logit(mean.pv)
  mean.pv ~ dunif(0, 1)
  beta.pv ~ dnorm(0, 0.1)
  kappa ~ dnorm(0, 0.1)

  # Observation equation
  # First year: no direct modeling of PS possible, average over PS effect
  #     using estimated proportion occupied as weight
  for (i in 1:nsites){
    for (j in 1:nsurveys[i,1]){
      V[i,j,1] ~ dbern(Pvisit[i,j,1])
      logit(Pvisit[i,j,1]) <- alpha.pv + beta.pv * (1 - 5.5) +
          n.occ.total[1]/274 * kappa
    }
  }

  # Later years: use of occupancy at t-1 as predictor of Pvisit
  for (i in 1:nsites){
    for (t in 2:nyears){
      for (j in 1:nsurveys[i,t]){
        V[i,j,t] ~ dbern(Pvisit[i,j,t])
        logit(Pvisit[i,j,t]) <- alpha.pv + beta.pv * (t - 5.5) +
            kappa * stateOcc[i,t-1]
      }
    }
  }
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(bdata$nsites, bdata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("kappa", "mean.psi", "alpha.lpsi", "trend.lpsi", "beta.lpsi", "mean.r",
    "alpha.lr", "trend.lr", "beta.lr", "alpha.lp2", "beta.lp2", "alpha.lp32",
    "alpha.lp33", "beta.lp32", "beta.lp33", "n.occ", "n.occ.total", "alpha.pv",
    "mean.pv", "beta.pv", "z")

# MCMC settings
# na <- 1000  ;  ni <- 60000   ;   nt <- 20   ;   nb <- 40000   ;   nc <- 3
na <- 1000  ;  ni <- 6000   ;   nt <- 2   ;   nb <- 4000   ;   nc <- 3  # ~~~~ testing, 1 hr

# Call JAGS from R, check convergence and summarize posteriors
out5 <- jags(bdata, inits, params, "staticMS2.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out5)
summary(out5)  ;   jags.View(out5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~ extra code for figures 6.10 and 6.11 ~~~~~~~~~~~~~~~~~~~~
# Figure 6.10
# '''''''''''
# Compare population size estimates with and without PS
year <- 2007:2016
op <- par(mfrow = c(1, 2), mar = c(5,5,3,3), cex.lab = 1.5, cex.axis = 1.5)
plot(year, out4$mean$n.occ[,2], xlab = 'Year', ylab = 'Number of sites',
    type = "b", lwd = 3, frame = FALSE, pch = 's', col = 'black', cex = 1.5,
    lty = 1, ylim = c(0, 274), main = 'No correction for PS (staticMS1)')
segments(year, out4$q2.5$n.occ[,2], year, out4$q97.5$n.occ[,2])
points(year, out4$mean$n.occ[,3], type = "b", lwd = 3, pch = 'p',
    col = 'black', cex = 1.5, lty = 1)
segments(year, out4$q2.5$n.occ[,3], year, out4$q97.5$n.occ[,3])
points(year, out4$mean$n.occ.total, type = "b", lwd = 3, pch = 't',
    col = 'black', cex = 1.5, lty = 1)
segments(year, out4$q2.5$n.occ.total, year, out4$q97.5$n.occ.total)

plot(year, out5$mean$n.occ[,2], xlab = 'Year', ylab = 'Number of sites',
    type = "b", lwd = 3, frame = FALSE, pch = 's', col = 'black', cex = 1.5,
    lty = 1, ylim = c(0, 274), main = 'With correction for PS (staticMS2)')
segments(year, out5$q2.5$n.occ[,2], year, out5$q97.5$n.occ[,2])
points(year, out5$mean$n.occ[,3], type = "b", lwd = 3, pch = 'p',
    col = 'black', cex = 1.5, lty = 1)
segments(year, out5$q2.5$n.occ[,3], year, out5$q97.5$n.occ[,3])
points(year, out5$mean$n.occ.total, type = "b", lwd = 3, pch = 't',
    col = 'black', cex = 1.5, lty = 1)
segments(year, out5$q2.5$n.occ.total, year, out5$q97.5$n.occ.total)
points(year, obsnocc, type = "b", lwd = 3, pch = 'o', col = 'black',
    cex = 1.5, lty = 1)
par(op)

# Figure 6.11
# '''''''''''
# Compare trends in psi and r in- and outside of the Alps
# And compare between models with and without PS
npred <- 500
tmp <- out4$mean     # Grab posterior means
tmpPS <- out5$mean   # Grab posterior means
pred.psi <- pred.r <- pred.psiPS <- pred.rPS <- array(NA, dim = c(npred, 2))
yearO <- seq(1, 10, length.out = npred)
year <- yearO -5.5

# Predict at the means of the other covariate (i.e., 0) for both models
for(r in 1:2){
  pred.psi[,r] <- plogis(tmp$alpha.lpsi[r] + tmp$trend.lpsi[r] * year)
  pred.r[,r] <- plogis(tmp$alpha.lr[r] + tmp$trend.lr[r] * year)
  pred.psiPS[,r] <- plogis(tmpPS$alpha.lpsi[r] + tmpPS$trend.lpsi[r] * year)
  pred.rPS[,r] <- plogis(tmpPS$alpha.lr[r] + tmpPS$trend.lr[r] * year)
}

op <- par(mfrow = c(2, 2), mar = c(5,5,5,4), cex.lab = 1.5, cex.axis = 1.5)
matplot(seq(2007, 2016,,npred), pred.psi, type = 'l', lty = c(2,1),
    lwd = 3, ylim = c(0, 1), xlab = "Year", ylab = "psi",
    main = "Trend in psi (without PS)", frame = FALSE)
legend('bottomleft', c("Inside of the Alps", "Outside of the Alps"),
    lty = c(1,2), lwd = 3, col = c('red', 'black'), bty = 'n', cex = 1.2)
matplot(seq(2007, 2016,,npred), pred.r, type = 'l', lty = c(2,1), lwd = 3,
    ylim = c(0, 1), xlab = "Year", ylab = "r",
    main = "Trend in r (without PS)", frame = FALSE)
matplot(seq(2007, 2016,,npred), pred.psiPS, type = 'l', lty = c(2,1),
    lwd = 3, ylim = c(0, 1), xlab = "Year", ylab = "psi",
    main = "Trend in psi (with PS)", frame = FALSE)
matplot(seq(2007, 2016,,npred), pred.rPS, type = 'l', lty = c(2,1),
    lwd = 3, ylim = c(0, 1), xlab = "Year", ylab = "r",
    main = "Trend in r (with PS)", frame = FALSE)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Inference on trends in psi outside/inside of the Alps
out5$summary[c('trend.lpsi[1]', 'trend.lpsi[2]'), c(1:3,7)]
#                    mean         sd      2.5%     97.5%
# trend.lpsi[1] 0.3344544 0.03176273 0.2731411 0.3977807 # outside
# trend.lpsi[2] 0.2809060 0.02378117 0.2359825 0.3283898 # inside

# Inference on trends in r outside/;inside of the Alps
out5$summary[c('trend.lr[1]', 'trend.lr[2]'), c(1:3,7)]
#                    mean         sd       2.5%       97.5%
# trend.lr[1] -0.02274265 0.05743434 -0.1315113  0.09133061
# trend.lr[2] -0.20243689 0.04395900 -0.2921581 -0.11656813

# ~~~~~~~~~~ extra code for figure 6.12 ~~~~~~~~~~~~~~~~
# Compute regional number of sites occupied by single birds and by
#   pairs for every year under the static model
pm.nsingles <- array(NA, dim = c(6, nyears))
CRI.nsingles <- array(NA, dim = c(6, nyears, 2))
pm.npairs <- array(NA, dim = c(6, nyears))
CRI.npairs <- array(NA, dim = c(6, nyears, 2))
for(r in 1:6){
  tmp1 <- out5$sims.list$z[,which(region == r),]
  for(t in 1:10){
    tmp2 <- tmp1[,,t]
    singles <- rowSums(tmp2 == 2)
    pm.nsingles[r, t] <- mean(singles)
    CRI.nsingles[r, t, ] <- hdi(singles) # could also use quantile()
    pairs <- rowSums(tmp2 == 3)
    pm.npairs[r, t] <- mean(pairs)
    CRI.npairs[r, t, ] <- hdi(pairs)
  }
}

op <- par(mfrow = c(1, 2), mar = c(5, 5, 2, 3), cex.lab = 1.5, cex.axis = 1.5)
year <- 2007:2016
matplot(year, t(pm.nsingles), type = 'b', lwd = 3, lty = 1, ylim = c(0, 45),
    frame = FALSE, ylab = 'Number of sites', xlab = 'Year',
    main = 'Single birds', las = 1)
segments(year, CRI.nsingles[1,,1], year, CRI.nsingles[1,,2], col = 1)
segments(year, CRI.nsingles[2,,1], year, CRI.nsingles[2,,2], col = 2)
segments(year, CRI.nsingles[3,,1], year, CRI.nsingles[3,,2], col = 3)
segments(year, CRI.nsingles[4,,1], year, CRI.nsingles[4,,2], col = 4)
segments(year, CRI.nsingles[5,,1], year, CRI.nsingles[5,,2], col = 5)
segments(year, CRI.nsingles[6,,1], year, CRI.nsingles[6,,2], col = 6)

matplot(year, t(pm.npairs), type = 'b', lwd = 3, lty = 1, ylim = c(0, 45),
    frame = FALSE, ylab = 'Number of sites', xlab = 'Year',
    main = 'Pairs', las = 1)
segments(year, CRI.npairs[1,,1], year, CRI.npairs[1,,2], col = 1)
segments(year, CRI.npairs[2,,1], year, CRI.npairs[2,,2], col = 2)
segments(year, CRI.npairs[3,,1], year, CRI.npairs[3,,2], col = 3)
segments(year, CRI.npairs[4,,1], year, CRI.npairs[4,,2], col = 4)
segments(year, CRI.npairs[5,,1], year, CRI.npairs[5,,2], col = 5)
segments(year, CRI.npairs[6,,1], year, CRI.npairs[6,,2], col = 6)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
