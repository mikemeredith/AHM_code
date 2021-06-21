#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6 : MULTISTATE OCCUPANCY MODELS
# =======================================
# Code from proofs dated 2020-08-19

# Approximate time to execute code in this script: 1.8 hrs
# With full number of iterations: 19 hrs

library(jagsUI)

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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 6.4 Case study: Swiss eagle owls
# ================================

# 6.4.4 Accounting for preferential sampling in the big dynamic model
# -------------------------------------------------------------------

# Get visitation data set V
V <- yms ; table(V)
V[V > 1] <- 1 ; table(V)
V[is.na(V)] <- 0 ; table(V)

# Data bundle
str(bdata <- list(y = yms, nsites = dim(yms)[1], nsurveys = nsurveys,
    nyears = dim(yms)[3], region = region, elev = elev.scaled,
    forest = forest.scaled, date = date.scaled, V = V))
# List of 9
# $ y       : num [1:274, 1:20, 1:10] NA NA 1 NA NA NA NA 2 2 2 ...
# $ nsites  : int 274
# $ nsurveys: num [1:274, 1:10] 1 1 3 1 1 1 1 20 12 20 ...
# $ nyears  : int 10
# $ region  : num [1:274] 2 5 6 1 1 2 3 5 3 2 ...
# $ elev    : num [1:274] -0.84 0.115 -0.554 -0.267 1.07 ...
# $ forest  : num [1:274] -0.647 0.862 -0.459 -0.459 -1.213 ...
# $ date    : num [1:274, 1:20, 1:10] 0 0 -1.9 0 0 ...
# $ V       : num [1:274, 1:20, 1:10] 0 0 1 0 0 0 0 1 1 1 ...

# Specify model in BUGS language
cat(file = "dynMS3.txt", "
model {
  ### (A) Model for multi-state Eagle Owl detections
  # -----------------------------------------------
  # ~~~~ code inserted from MS dated 2019-03-26 ~~~~~~~~~~~~~~~~~~~~~~~~
  ### (1) Linear models and priors
  # ------------------------------
  ## (a) State process
  # Initial state vector Omega
  # Linear models
  for (i in 1:nsites){
    logit(psi[i]) <- alpha.lpsi[region[i]] + beta.lpsi[1] * elev[i] + beta.lpsi[2] * forest[i]
    logit(r[i]) <- alpha.lr[region[i]] + beta.lr[1] * elev[i] + beta.lr[2] * forest[i]
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
  # Hyperpriors for hyperparameters governing these random year effects
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
  beta.region.lphi1[1] <- 0      # Avoid overparameterization
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
      for(j in 1:nsurveys[i,t]){   # note nsurvey function of site and year
        # Observation model for sites in occupied state 1 (= single bird)
        logit(p2[i,j,t]) <- alpha.lp2 + beta.region.lp2[region[i]] + beta.lp2[1] * date[i,j,t] + beta.lp2[2] * pow(date[i,j,t],2)
        # Observation model for sites in occupied state 2 (= pairs)
        mlogit.p3[2,i,j,t] <- alpha.lp32 + beta.region.lp32[region[i]] + beta.lp32[1] * date[i,j,t] + beta.lp32[2] * pow(date[i,j,t],2)
        mlogit.p3[3,i,j,t] <- alpha.lp33 + beta.region.lp33[region[i]] + beta.lp33[1] * date[i,j,t] + beta.lp33[2] * pow(date[i,j,t],2)
      }
    }
  }

  # Priors for parameters in the linear models in Theta
  # Intercepts
  alpha.lp2 <- logit(mean.alpha.p2)
  mean.alpha.p2 ~ dunif(0, 1)
  alpha.lp32 ~ dnorm(0, 0.01) # Must be normal for multinomial logit
  alpha.lp33 ~ dnorm(0, 0.01)
  # Fixed effects of region in parameters in Theta
  beta.region.lp2[1] <- 0      # Avoid overparameterization
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
        p3[2,i,j,t] <- exp(mlogit.p3[2,i,j,t]) / (1 + exp(mlogit.p3[2,i,j,t]) + exp(mlogit.p3[3,i,j,t]))
        p3[3,i,j,t] <- exp(mlogit.p3[3,i,j,t]) / (1 + exp(mlogit.p3[2,i,j,t]) + exp(mlogit.p3[3,i,j,t]))
      }
    }
  }


  ### (2) Define Initial state vector (Omega), state transition matrix (PhiMat) and observation matrix (Theta)
  # Initial state vector (Omega): Year 1
  for (i in 1:nsites){
    Omega[i,1] <- 1 - psi[i]         # Prob. of non-occupation
    Omega[i,2] <- psi[i] * (1-r[i])  # Prob. of occ. by single bird
    Omega[i,3] <- psi[i] * r[i]      # Prob. of occ. by pair
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
  # ~~~~ end of inserted code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

  ### (B) Model for site visitation
  # -------------------------------
  # Priors
  alpha.pv <- logit(mean.pv)
  mean.pv ~ dunif(0, 1)
  beta.pv ~ dnorm(0, 0.1)
  kappa ~ dnorm(0, 0.1)

  # Observation equation
  # First year: no direct modeling of PS possible, average over PS effect
  # using estimated proportion occupied as weight
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

# ~~~~~ Full list of parameters monitored inserted from MS ~~~~~~~~~~~~
params <- c("kappa", "alpha.lpsi", "mean.psi", "beta.lpsi", "alpha.lr", "mean.r",
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
    "beta.region.lp2", "beta.region.lp32", "beta.region.lp33", "beta.lp2", "beta.lp32", "beta.lp33", "PhiMat.avg", "n.occ", "n.occ.total",
    "alpha.pv", "mean.pv", "beta.pv", "z")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MCMC settings
# na <- 1000 ; ni <- 60000 ; nt <- 20 ; nb <- 40000 ; nc <- 3
na <- 1000 ; ni <- 6000 ; nt <- 2 ; nb <- 4000 ; nc <- 3  # ~~~ for testing, 2 hrs

# Call JAGS (ART 24 h), check convergence and summarize posteriors
odms3 <- jags(bdata, inits, params, "dynMS3.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(odms3)
summary(odms3) ; jags.View(odms3) # not shown

round(odms3$mean$PhiMat.avg, 2) # Transition probabilities
#      [,1] [,2] [,3]
# [1,] 0.67 0.18 0.15
# [2,] 0.20 0.70 0.09
# [3,] 0.10 0.06 0.84

# ~~~~~~~~~~~~~~ extra code for figures 6.5 to 6.9 ~~~~~~~~~~~~~~~~~~~~~~~

# Figure 6.5
# ''''''''''
# Compare population size estimates
load("AHM2_06.4.3_odms2.RData")  # Need output for comparison
year <- 2007:2016
op <- par(mfrow = c(1, 2), mar = c(5,5,3,3), cex.lab = 1.5, cex.axis = 1.5)
plot(year, odms2$mean$n.occ[,2], xlab = 'Year', ylab = 'Number of sites',
    type = "b", lwd = 3, frame = FALSE, pch = 's', col = 'black', cex = 1.5,
    lty = 1, ylim = c(0, 274), main = 'No correction for PS (dynMS2)')
segments(year, odms2$q2.5$n.occ[,2], year, odms2$q97.5$n.occ[,2])
points(year, odms2$mean$n.occ[,3], type = "b", lwd = 3, pch = 'p',
    col = 'black', cex = 1.5, lty = 1)
segments(year, odms2$q2.5$n.occ[,3], year, odms2$q97.5$n.occ[,3])
points(year, odms2$mean$n.occ.total, type = "b", lwd = 3, pch = 't',
    col = 'black', cex = 1.5, lty = 1)
segments(year, odms2$q2.5$n.occ.total, year, odms2$q97.5$n.occ.total)

plot(year, odms3$mean$n.occ[,2], xlab = 'Year', ylab = '',
    type = "b", lwd = 3, frame = FALSE, pch = 's', col = 'black', cex = 1.5,
    lty = 1, ylim = c(0, 274), main = 'With correction for PS (dynMS3)')
segments(year, odms3$q2.5$n.occ[,2], year, odms3$q97.5$n.occ[,2])
points(year, odms3$mean$n.occ[,3], type = "b", lwd = 3, pch = 'p',
    col = 'black', cex = 1.5, lty = 1)
segments(year, odms3$q2.5$n.occ[,3], year, odms3$q97.5$n.occ[,3])
points(year, odms3$mean$n.occ.total, type = "b", lwd = 3, pch = 't',
    col = 'black', cex = 1.5, lty = 1)
segments(year, odms3$q2.5$n.occ.total, year, odms3$q97.5$n.occ.total)
points(year, obsnocc, type = "b", lwd = 3, pch = 'o', col = 'black',
    cex = 1.3, lty = 1)
par(op)

# Figure 6.6
# ''''''''''
# Plot estimates of z matrix (rounded to integers)
op <- par(mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
mapPalette <- colorRampPalette(c("white", "black"))
image(x = year, y = 1:25, z = round(t(odms3$mean$z[1:25,])),
    col = mapPalette(10), axes = TRUE, xlab = "Year", ylab = "Site", main = '')
par(op)

# Figure 6.7
# ''''''''''
# Compute regional number of sites occupied by single birds and by pairs
#    for every year (code courtesy of Mike)
library(HDInterval)
str(odms3$sims.list$z)
pm.nsingles <- array(NA, dim = c(6, nyears))
CRI.nsingles <- array(NA, dim = c(6, nyears, 2))
pm.npairs <- array(NA, dim = c(6, nyears))
CRI.npairs <- array(NA, dim = c(6, nyears, 2))
for(r in 1:6){
  tmp1 <- odms3$sims.list$z[,which(region == r),]
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

op <- par(mfrow = c(1, 2), mar = c(5, 5, 5, 4), cex.lab = 1.5, cex.axis = 1.5,
    cex.main = 1.5)
matplot(year, t(pm.nsingles), type = 'b', lwd = 3, lty = 1, ylim = c(0, 45),
    frame = FALSE, ylab = 'Number of occupied sites', xlab = 'Year',
    main = 'Single birds', las = 1)
segments(year, CRI.nsingles[1,,1], year, CRI.nsingles[1,,2], col = 1)
segments(year, CRI.nsingles[2,,1], year, CRI.nsingles[2,,2], col = 2)
segments(year, CRI.nsingles[3,,1], year, CRI.nsingles[3,,2], col = 3)
segments(year, CRI.nsingles[4,,1], year, CRI.nsingles[4,,2], col = 4)
segments(year, CRI.nsingles[5,,1], year, CRI.nsingles[5,,2], col = 5)
segments(year, CRI.nsingles[6,,1], year, CRI.nsingles[6,,2], col = 6)

matplot(year, t(pm.npairs), type = 'b', lwd = 3, lty = 1, ylim = c(0, 45),
    frame = FALSE, ylab = 'Number of occupied sites', xlab = 'Year',
    main = 'Pairs', las = 1)
segments(year, CRI.npairs[1,,1], year, CRI.npairs[1,,2], col = 1)
segments(year, CRI.npairs[2,,1], year, CRI.npairs[2,,2], col = 2)
segments(year, CRI.npairs[3,,1], year, CRI.npairs[3,,2], col = 3)
segments(year, CRI.npairs[4,,1], year, CRI.npairs[4,,2], col = 4)
segments(year, CRI.npairs[5,,1], year, CRI.npairs[5,,2], col = 5)
segments(year, CRI.npairs[6,,1], year, CRI.npairs[6,,2], col = 6)
par(op)

round(100*odms3$mean$PhiMat.avg, 0)    # Transition probabilities in %

# Figure 6.8
# ''''''''''
# Some predictions for covariate effects in Omega
# Predictions of psi and r for region, elevation and forest cover
tmp <- odms3$mean    # Grab posterior means
npred <- 500
pred.psi.elev <- pred.r.elev <- pred.psi.forest <- pred.r.forest <-
    array(NA, dim = c(npred, 6))
elevO <- seq(min(elev), max(elev), length.out = npred)
forestO <- seq(0, 1, length.out = npred)
elevP <- standardize2match(elevO, elev)  # require(AHMbook)
forestP <- standardize2match(forestO, forest)

# Predict at the means of the other covariate (i.e., 0)
for(r in 1:6){
  pred.psi.elev[,r] <- plogis(tmp$alpha.lpsi[r] + tmp$beta.lpsi[1] * elevP)
  pred.r.elev[,r] <- plogis(tmp$alpha.lr[r] + tmp$beta.lr[1] * elevP)
  pred.psi.forest[,r] <- plogis(tmp$alpha.lpsi[r] + tmp$beta.lpsi[2] * forestP)
  pred.r.forest[,r] <- plogis(tmp$alpha.lr[r] + tmp$beta.lr[2] * forestP)
}

op <- par(mfrow = c(2, 2), mar = c(5,5,5,4), cex.lab = 1.5, cex.axis = 1.5)
matplot(elevO, pred.psi.elev, type = 'l', lty = 1, lwd = 3, ylim = c(0, 1),
    xlab = "Elevation (m a.s.l)", ylab = "psi",
    main = "Elevation effect on psi in Omega", frame = FALSE)
matplot(elevO, pred.r.elev, type = 'l', lty = 1, lwd = 3, ylim = c(0, 1),
    xlab = "Elevation (m a.s.l)", ylab = "r",
    main = "Elevation effect on r in Omega", frame = FALSE)
matplot(100*forestO, pred.psi.forest, type = 'l', lty = 1, lwd = 3,
    ylim = c(0, 1), xlab = "Forest (%)", ylab = "psi",
    main = "Forest effect on psi in Omega", frame = FALSE)
matplot(100*forestO, pred.r.forest, type = 'l', lty = 1, lwd = 3,
    ylim = c(0, 1), xlab = "Forest (%)", ylab = "r",
    main = "Forest effect on r in Omega", frame = FALSE)
par(op)

# Figure 6.9
# ''''''''''
# Some predictions for covariate effects in Observation matrix Theta
# Predictions of p2 and p31, p32 and p33 for region and survey date
pred.p2 <- mlogit.p32 <- mlogit.p33 <- pred.p31 <- pred.p32 <- pred.p33 <-
    array(NA, dim = c(npred, 6))
dateO <- seq(1, 281, length.out = npred)
dateP <- (dateO - 90) / 30.5

# Predict at the means of the other covariate (i.e., 0)
for(r in 1:6){
  pred.p2[,r] <- plogis(tmp$alpha.lp2 + tmp$beta.region.lp2[r] +
      tmp$beta.lp2[1] * dateP + tmp$beta.lp2[2] * dateP^2)
  mlogit.p32[,r] <- plogis(tmp$alpha.lp32 + tmp$beta.region.lp32[r] +
      tmp$beta.lp32[1] * dateP + tmp$beta.lp32[2] * dateP^2)
  mlogit.p33[,r] <- plogis(tmp$alpha.lp33 + tmp$beta.region.lp33[r] +
      tmp$beta.lp33[1] * dateP + tmp$beta.lp33[2] * dateP^2)
  pred.p32[,r] <- exp(mlogit.p32[,r]) / (1 + exp(mlogit.p32[,r]) +
      exp(mlogit.p33[,r]))
  pred.p33[,r] <- exp(mlogit.p33[,r]) / (1 + exp(mlogit.p32[,r]) +
      exp(mlogit.p33[,r]))
  pred.p31[,r] <- 1 - pred.p32[,r] - pred.p33[,r]
}

# Prediction of the visitation model
pred.pv <- array(NA, dim = c(10, 2))
pred.pv[1,] <- plogis(tmp$alpha.pv + tmp$beta.pv * (1-5.5) +
    tmp$n.occ.total[1]/274 * tmp$kappa)
for(t in 2:10){
  pred.pv[t, 1] <- plogis(tmp$alpha.pv + tmp$beta.pv * (t-5.5))
  pred.pv[t, 2] <- plogis(tmp$alpha.pv + tmp$beta.pv * (1 - 5.5) + tmp$kappa)
}

op <- par(mfrow = c(2, 2), mar = c(5, 5, 4, 4), cex.lab = 1.5, cex.axis = 1.5,
    cex.main = 2)
matplot(dateO, pred.p2, type = 'l', lty = 1, lwd = 3, ylim = c(0, 1.1),
    xlab = "Survey date (Julian day)", ylab = "p2", main = 'A', frame = FALSE)
    # main = "Date effect on p2 in Theta"
matplot(dateO, pred.p32, type = 'l', lty = 1, lwd = 3, ylim = c(0, 1.1),
    xlab = "Survey date (Julian day)", ylab = "p32", main = 'B', frame = FALSE)
    # main = "Date effect on p32 in Theta"
matplot(dateO, pred.p33, type = 'l', lty = 1, lwd = 3, ylim = c(0, 1.1),
    xlab = "Survey date (Julian day)", ylab = "p33", main = 'C', frame = FALSE)
    # main = "Date effect on p33 in Theta",
matplot(2007:2016, pred.pv, type = 'l', lty = c(2,1), lwd = 3, ylim = c(0, 1.1),
    xlab = "Year", ylab = "Visitation prob.", main = 'D', frame = FALSE)
    # main = "Visitation model"
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
