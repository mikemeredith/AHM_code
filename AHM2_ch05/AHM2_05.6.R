#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5 : MODELING METACOMMUNITY DYNAMICS USING DYNAMIC COMMUNITY MODELS
# ==========================================================================
# Code from proofs dated 2020-08-19

# Approximate run time for this script: 17 mins
# With full number of iterations: 2.1 hrs

library(AHMbook)
library(jagsUI)

# 5.6 Using the DCM in meta-analyses and comparative studies
# ==========================================================

# 5.6.1 DCM in a meta-analysis to gauge the effects of an intervention in a BACI design
# -------------------------------------------------------------------------------------

# Vectors of first year, study duration and year of impact
initial.year <- c(5, 1, 1, 6, 8, 8, 9, 12, 13, 18, 10, 15)
dur <- c(19, 26, 24, 14, 16, 19, 16, 11, 14, 9, 17, 12)
impact <- c(11, 23, 19, 11, 11, 11, 13, 8, 10, 6, 11, 5)

# Generate detection/nondetection data for all 12 species
library(AHMbook)
data <- list()
set.seed(1)
for(s in 1:12){
  data[[s]] <- simDynocc(nsites = 100, nsurveys = 2, nyears = dur[s],
    year.of.impact = impact[s], mean.psi1 = 0.5, range.phi = c(0.5, 1),
    impact.phi = 25, range.gamma = c(0, 0.5), impact.gamma = 25,
    range.p = c(0.1, 0.9), show.plot = FALSE)
}

# Data for all 12 species
y <- array(NA, dim = c(100, 2, 26, 12))      # Detection/nondetection
imat <- array(NA, dim = c(26, 12))           # Binary impact covariate

# Assemble detection/nondetection data
for(s in 1:12){
  (data.set <- paste('data[[',s,']]$y', sep = ''))
  y[,,initial.year[s]:(initial.year[s]+dur[s]-1),s] <-
      eval(parse(text = data.set)) # Thanx to Fränzi Korner for figuring this out
}

# Quick sum check .... looks OK
sum(y, na.rm = TRUE) ; sum(sapply(data, function(x) sum(x$y)))

# Assemble impact data matrix and look at it
for(t in 1:12){
  imat[1:(initial.year[t]+impact[t]-1),t] <- 0
  imat[(initial.year[t]+impact[t]-1):26,t] <- 1
}
colnames(imat) <- paste('Study',1:12)
imat                                         # Look at the impact matrix

# Bundle and summarize data set
str(bdata <- list(y = y, nsite = dim(y)[1], nsurvey = dim(y)[2], nyear = dim(y)[3],
    nspec = dim(y)[4], imat = imat))
# List of 6
# $ y      : int [1:100, 1:2, 1:26, 1:12] NA NA NA NA NA NA NA NA ...
# $ nsite  : int 100
# $ nsurvey: int 2
# $ nyear  : int 26
# $ nspec  : int 12
# $ imat   : num [1:26, 1:12] 0 0 0 0 0 0 0 0 0 0 ...

# Specify model in BUGS language
cat(file = "DCM3.txt", "
model {

  # Specify linear models for parameters
  for(k in 1:nspec){
    psi1[k] ~ dbeta(1, 1)
  }

  # Persistence and colonization have effects of impact plus random year effects
  for(k in 1:nspec){
    for(t in 1:(nyear-1)){
      logit(phi[t, k]) <- lphi[t, k]
      lphi[t, k] ~ dnorm(mu.lphi[imat[t, k]+1, k], tau.lphi.year[k])
      logit(gamma[t, k]) <- lgamma[t, k]
      lgamma[t, k] ~ dnorm(mu.lgamma[imat[t, k]+1, k], tau.lgamma.year[k])
    }
    tau.lphi.year[k] <- pow(sd.lphi.year[k],-2)
    sd.lphi.year[k] ~ dunif(0.01, 10)
    tau.lgamma.year[k] <- pow(sd.lgamma.year[k],-2)
    sd.lgamma.year[k] ~ dunif(0.01, 10)
    for(h in 1:2){ # Loop over period before and after
      mu.lphi[h, k] <- logit(mean.phi[h,k])
      mean.phi[h,k] ~ dbeta(1, 1)
      mu.lgamma[h, k] <- logit(mean.gamma[h,k])
      mean.gamma[h,k] ~ dbeta(1, 1)
    }
  }

  # Detection has only random year effects
  for(k in 1:nspec){
    for(t in 1:nyear){
      logit(p[t,k]) <- lp[t,k]
      lp[t, k] ~ dnorm(mu.lp[k], tau.lp.year[k])
    }
    mu.lp[k] <- logit(mean.p[k])
    mean.p[k] ~ dbeta(1, 1)
    tau.lp.year[k] <- pow(sd.lp.year[k],-2)
    sd.lp.year[k] ~ dunif(0.01, 10)
  }

  # Ecological and observation submodels
  for (i in 1:nsite){ # Loop over sites
      for(k in 1:nspec){ # Loop over species
      # Initial conditions of system
      z[i,1, k] ~ dbern(psi1[k])            # Presence/absence at start of study
      # State transitions
      for (t in 2:nyear){ # Loop over years
        z[i,t,k] ~ dbern(z[i,t-1,k] * phi[t-1, k] + (1-z[i,t-1, k]) * gamma[t-1, k])
        # Observation model
        for (j in 1:nsurvey){
          y[i,j,t,k] ~ dbern(z[i,t,k] * p[t, k])
        }
      }
    }
  }

  # Derived quantities
  # Estimate effect size of impact on persistence and colonization at
  # species level (probability scale)
  for(k in 1:nspec){                       # Loop over species
    effect.phi.spec[k] <- mean.phi[2,k] - mean.phi[1, k]
    effect.gamma.spec[k] <- mean.gamma[2,k] - mean.gamma[1,k]
  }
}
")

# Initial values
zst <- apply(y, c(1,3,4), max) # Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("effect.phi.spec", "effect.gamma.spec", "mean.phi",
    "sd.lphi.year", "mean.gamma", "sd.lgamma.year", "mean.p", "sd.lp.year")

# MCMC settings
# na <- 1000 ; ni <- 20000 ; nt <- 10 ; nb <- 10000 ; nc <- 3
na <- 1000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~~ for testing, 5 mins

# Call JAGS (ART 61 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "DCM3.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out3)
jags.View(out3) ; print(out3$summary[1:24, c(1:3,5,7:10)], 2) # not shown


# Specify model in BUGS language
cat(file = "DCM4.txt", "
model {

  # Specify linear models for parameters
  for(k in 1:nspec){
    psi1[k] ~ dbeta(1, 1)
  }

  # Persistence and colonization have effect of impact + random year effects
  for(k in 1:nspec){
    for(t in 1:(nyear-1)){
      logit(phi[t, k]) <- lphi[t, k]
      lphi[t, k] ~ dnorm(mu.lphi[imat[t, k]+1, k], tau.lphi.year[k])
      logit(gamma[t, k]) <- lgamma[t, k]
      lgamma[t, k] ~ dnorm(mu.lgamma[imat[t, k]+1, k], tau.lgamma.year[k])
    }

    tau.lphi.year[k] <- pow(sd.lphi.year[k],-2)
    sd.lphi.year[k] ~ dunif(0.01, 10)    # Random year effects for lphi
    tau.lgamma.year[k] <- pow(sd.lgamma.year[k],-2)
    sd.lgamma.year[k] ~ dunif(0.01, 10)  # Random year effects for lgamma
    # Specify priors for mu.lphi[1,k], mu.lphi[2,k] and for mu.lgamma[1,k]
    # and mu.lgamma[2,k] (assume same variance across species before and after)
    mu.lphi[1, k] ~ dnorm(mean.lphi.before, tau.lphi.spec)
    mu.lphi[2, k] ~ dnorm(mean.lphi.after, tau.lphi.spec)
    mu.lgamma[1, k] ~ dnorm(mean.lgamma.before, tau.lgamma.spec)
    mu.lgamma[2, k] ~ dnorm(mean.lgamma.after, tau.lgamma.spec)
  }

  # Hyperpriors for these hyperparameters
  mean.lphi.before <- logit(mean.phi.before)
  mean.phi.before ~ dbeta(1, 1)          # Mean across species of phi.before
  mean.lphi.after <- logit(mean.phi.after)
  mean.phi.after ~ dbeta(1, 1)           # Mean across species of phi.after
  mean.lgamma.before <- logit(mean.gamma.before)
  mean.gamma.before ~ dbeta(1, 1)        # Mean across species of gamma.before
  mean.lgamma.after <- logit(mean.gamma.after)
  mean.gamma.after ~ dbeta(1, 1)         # Mean across species of gamma.after
  tau.lphi.spec <- pow(sd.lphi.spec,-2)
  sd.lphi.spec ~ dunif(0.01, 10)         # Common SD among species of mu.lphi
  tau.lgamma.spec <- pow(sd.lgamma.spec,-2)
  sd.lgamma.spec ~ dunif(0.01, 10)       # Common SD among species of mu.lgamma
  # Detection same as for model 3
  for(k in 1:nspec){
    for(t in 1:nyear){
      logit(p[t,k]) <- lp[t,k]
      lp[t, k] ~ dnorm(mu.lp[k], tau.lp.year[k])
    }
    mu.lp[k] <- logit(mean.p[k])
    mean.p[k] ~ dbeta(1, 1)
    tau.lp.year[k] <- pow(sd.lp.year[k],-2)
    sd.lp.year[k] ~ dunif(0.01, 10)
  }

  # Ecological and observation submodels
  for (i in 1:nsite){                    # Loop over sites
    for(k in 1:nspec){                   # Loop over species
      # Initial conditions of system
      z[i,1, k] ~ dbern(psi1[k])         # Presence/absence at start of study
      # State transitions
      for (t in 2:nyear){ # Loop over years
        z[i,t,k] ~ dbern(z[i,t-1,k] * phi[t-1, k] + (1-z[i,t-1, k]) * gamma[t-1, k])
        # Observation model
        for (j in 1:nsurvey){
          y[i,j,t,k] ~ dbern(z[i,t,k] * p[t, k])
        }
      }
    }
  }

  # Derived quantities
  # Estimate effect size of impact on persistence and colonization at
  # community level (probability scale)
  effect.phi.comm <- mean.phi.after-mean.phi.before
  effect.gamma.comm <- mean.gamma.after-mean.gamma.before
  # Estimate effect size of impact on persistence and colonization at
  # species level (probability scale)
  for(k in 1:nspec){                  # Loop over species
    effect.phi.spec[k] <- ilogit(mu.lphi[2,k]) - ilogit(mu.lphi[1, k])
    effect.gamma.spec[k] <- ilogit(mu.lgamma[2,k]) - ilogit(mu.lgamma[1,k])
  }
}
")

# Parameters monitored
params <- c("effect.phi.comm", "effect.gamma.comm", "effect.phi.spec",
    "effect.gamma.spec", "mu.lphi", "sd.lphi.year", "mu.lgamma",
    "sd.lgamma.year", "mean.p", "sd.lp.year", "mean.phi.before",
    "mean.phi.after", "mean.gamma.before", "mean.gamma.after",
    "sd.lphi.spec", "sd.lgamma.spec")

# Call JAGS (ART 69 min), check convergence and summarize posteriors
out4 <- jags(bdata, inits, params, "DCM4.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out4)
jags.View(out4) ; print(out4$summary[1:24, c(1:3,5,7:10)], 2) # not shown

# ~~~~~~~ extra code for figure 5.4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualize difference before-after under fixed- and random-eff. models
off <- 0.1              # Graphical offset
op <- par(mfrow = c(1,2), mar=c(5,6,4,2)+0.1)
plot(out3$mean$effect.phi.spec, 1:12-off, frame = FALSE, yaxt='n',
    xlim = c(-0.5, 0.25), ylim = c(0, 12),
    main = 'Persistence', xlab = 'Difference after-before',
    ylab = '', pch = 1, cex = 1.2)
segments(out3$q2.5$effect.phi.spec, 1:12-off, out3$q97.5$effect.phi.spec,
    1:12-off)
axis(2, at = c(0,1:12), labels = c('Hypermean', paste('Study', 1:12)), las = 1)
points(out4$mean$effect.phi.comm, 0, pch = 16)
segments(out4$q2.5$effect.phi.comm, 0, out4$q97.5$effect.phi.comm, 0)
abline(h = 0.5, lty = 2)
points(out4$mean$effect.phi.spec, 1:12+off, pch = 16, cex = 1.2)
segments(out4$q2.5$effect.phi.spec, 1:12+off, out4$q97.5$effect.phi.spec,
    1:12+off)
abline(v = 0, lwd = 1, col = 'grey')

plot(out3$mean$effect.gamma.spec, 1:12-off, frame = FALSE, yaxt='n',
    xlim = c(-0.3, 0.3), ylim = c(0, 12),
    main = 'Colonization', xlab = 'Difference after-before',
    ylab = '', pch = 1, cex = 1.2)
segments(out3$q2.5$effect.gamma.spec, 1:12-off, out3$q97.5$effect.gamma.spec,
    1:12-off)
axis(2, at = c(0,1:12), labels = c('Hypermean', paste('Study', 1:12)), las = 1)
points(out4$mean$effect.gamma.comm, 0, pch = 16, cex = 1.5)
segments(out4$q2.5$effect.gamma.comm, 0, out4$q97.5$effect.gamma.comm, 0)
abline(h = 0.5, lty = 2)
points(out4$mean$effect.gamma.spec, 1:12+off, pch = 16, cex = 1.2)
segments(out4$q2.5$effect.gamma.spec, 1:12+off, out4$q97.5$effect.gamma.spec,
    1:12+off)
abline(v = 0, lwd = 1, col = 'grey')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.6.2 DCM in a comparative study to explain differences among multiple species
# ------------------------------------------------------------------------------

# Generate mass data for 12 species (sorted for graphic convenience)
set.seed(62016)
nspec <- 12
mass <- sort(rlnorm(nspec, meanlog = log(1), sdlog = 1.8)) # in kg

# Choose logit-linear allometric relationships for expectations
epsi1 <- plogis(0 - 1 * log(mass))       # Params (0, -1)
ephi <- plogis(1 + 1 * log(mass))        # Params (1, 1)
egamma <- plogis(-1 -1 * log(mass))      # Params (-1, -1)
ep <- plogis(0 - 0.5 * log(mass))        # Params (1, -0.5)

# Add a little species-specific noise on top of each regression
sd <- 0.5                        # Variation among species (on logit scale)
psi1 <- plogis(qlogis(epsi1) + rnorm(nspec, 0, sd))
phi <- plogis(qlogis(ephi) + rnorm(nspec, 0, sd))
gamma <- plogis(qlogis(egamma) + rnorm(nspec, 0, sd))
p <- plogis(qlogis(ep) + rnorm(nspec, 0, sd))
cbind(psi1, phi, gamma, p)               # Inspect species-specific values

# ~~~~~~~ extra code for figure 5.5 ~~~~~~~~~~~~~~~~~~~~
cx <- 1.2  ;  xlim <- c(0, 10)  ;  ylim = c(0, 1)
op <- par(mfrow = c(2,2))
curve(plogis(0 - 1 * log(x)), 0, 10, xlab = "Mass (kg)", ylab = "Probability",
    xlim = xlim, ylim = ylim, frame = FALSE, main = 'psi1', lwd = 2)
points(mass, psi1, pch = 16, cex = cx)
curve(plogis(1 + 1 * log(x)), 0, 10, xlab = "Mass (kg)", ylab = "Probability",
    xlim = xlim, ylim = ylim, frame = FALSE, main = 'phi', lwd = 2)
points(mass, phi, pch = 16, cex = cx)
curve(plogis(-1 -1 * log(x)), 0, 10, xlab = "Mass (kg)", ylab = "Probability",
    xlim = xlim, ylim = ylim, frame = FALSE, main = 'gamma', lwd = 2)
points(mass, gamma, pch = 16, cex = cx)
curve(plogis(0 - 0.5 * log(x)), 0, 10, xlab = "Mass (kg)", ylab = "Probability",
    xlim = xlim, ylim = ylim, frame = FALSE, main = 'p', lwd = 2)
points(mass, p, pch = 16, cex = cx)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Simulate a data set for each species with constant parameters
data <- list()
set.seed(1)
for(s in 1:12){
  data[[s]] <- simDynocc(nsites = 100, nsurveys = 3, nyear = 10,
      mean.psi1 = psi1[s], range.phi = c(phi[s], phi[s]),
      range.gamma = c(gamma[s], gamma[s]), range.p = c(p[s], p[s]),
      show.plot = FALSE)
}

# Data for all 12 species
y <- array(NA, dim = c(100, 3, 10, 12))  # Detection/nondetection

# Assemble detection/nondetection data in 4D array
for(s in 1:12){
  (data.set <- paste('data[[',s,']]$y', sep = ''))
  y[,,,s] <- eval(parse(text=data.set))  # Thanks, Fränzi !
}

# Quick sum check .... looks OK
sum(y, na.rm = TRUE) ; sum(sapply(data, function(x) sum(x$y)))

# Bundle and summarize data set
pred.logmass <- log(seq(0.1, 10, length.out = 100))
str(bdata <- list(y = y, nsite = dim(y)[1], nsurvey = dim(y)[2],
    nyear = dim(y)[3], nspec = dim(y)[4], logmass = log(mass),
    pred.logmass = pred.logmass) )
# List of 7
# $ y           : int [1:100, 1:3, 1:10, 1:12] 1 1 1 1 1 1 1 0 0 1 ...
# $ nsite       : int 100
# $ nsurvey     : int 3
# $ nyear       : int 10
# $ nspec       : int 12
# $ logmass     : num [1:12] -2.889 -1.429 0.57 0.817 0.871 ...
# $ pred.logmass: num [1:100] -2.303 -1.609 -1.204 -0.916 -0.693 ...

# Specify model in BUGS language
cat(file = "DCM5.txt", "
model {

  # Specify linear models for parameters
  # For initial occupancy
  for(k in 1:nspec){
    logit(psi1[k]) <- alpha.lpsi1 + beta.lpsi1 * logmass[k]
  }
  alpha.lpsi1 <- logit(int.psi1)
  int.psi1 ~ dbeta(1, 1)
  beta.lpsi1 ~ dnorm(0, 0.1)

  # For persistence and colonization
  for(k in 1:nspec){
    for(t in 1:(nyear-1)){
      logit(phi[t, k]) <- lphi[t, k]
      lphi[t, k] ~ dnorm(mu.lphi[k], tau.lphi)
      logit(gamma[t, k]) <- lgamma[t, k]
      lgamma[t, k] ~ dnorm(mu.lgamma[k], tau.lgamma)
    }
    mu.lphi[k] <- alpha.lphi + beta.lphi * logmass[k]
    mu.lgamma[k] <- alpha.lgamma + beta.lgamma * logmass[k]
  }
  alpha.lphi <- logit(int.phi)
  int.phi ~ dbeta(1, 1)
  beta.lphi ~ dnorm(0, 0.1)
  alpha.lgamma <- logit(int.gamma)
  int.gamma ~ dbeta(1, 1)
  beta.lgamma ~ dnorm(0, 0.1)
  tau.lphi <- pow(sd.lphi,-2)
  sd.lphi ~ dunif(0, 3)
  tau.lgamma <- pow(sd.lgamma,-2)
  sd.lgamma ~ dunif(0, 3)

  # For detection
  for(k in 1:nspec){
    for(t in 1:nyear){
      logit(p[t,k]) <- lp[t,k]
      lp[t, k] ~ dnorm(mu.lp[k], tau.lp)
    }
    mu.lp[k] <- alpha.lp + beta.lp * logmass[k]
  }
  alpha.lp <- logit(int.p)
  int.p ~ dbeta(1, 1)
  beta.lp ~ dnorm(0, 0.1)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 3)

  # Ecological and observation submodels
  for (i in 1:nsite){
    for(k in 1:nspec){
      # Initial conditions of system
      z[i,1, k] ~ dbern(psi1[k])
      # State transitions
      for (t in 2:nyear){
        z[i,t,k] ~ dbern(z[i,t-1,k]*phi[t-1,k]+(1-z[i,t-1,k])*gamma[t-1,k])
        # Observation model
        for (j in 1:nsurvey){
          y[i,j,t,k] ~ dbern(z[i,t,k] * p[t, k])
        }
      }
    }
  }

  # Predictions of the allometric relationships for all 4 params
  for (i in 1:100){
    logit(pred[i, 1]) <- alpha.lpsi1 + beta.lpsi1 * pred.logmass[i]
    logit(pred[i, 2]) <- alpha.lphi + beta.lphi * pred.logmass[i]
    logit(pred[i, 3]) <- alpha.lgamma + beta.lgamma * pred.logmass[i]
    logit(pred[i, 4]) <- alpha.lp + beta.lp * pred.logmass[i]
  }
}
")

# Initial values
inits <- function(){ list(z = apply(y, c(1,3,4), max))}

# Parameters monitored
params <- c("int.psi1", "alpha.lpsi1", "beta.lpsi1", "int.phi", "alpha.lphi",
    "beta.lphi", "sd.lphi", "int.gamma", "alpha.lgamma", "beta.lgamma",
    "sd.lgamma", "int.p", "alpha.lp", "beta.lp", "sd.lp", "pred")

# MCMC settings
# na <- 1000 ; ni <- 50000 ; nt <- 25 ; nb <- 25000 ; nc <- 3
na <- 1000 ; ni <- 5000 ; nt <- 2 ; nb <- 2500 ; nc <- 3  # ~~~ for testing, 6 mins

# Call JAGS (ART 70 min), check convergence and summarize posteriors
out5 <- jags(bdata, inits, params, "DCM5.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out5)
jags.View(out5) ; print(out5, 2) # not shown


# ~~~~~~~~~~~ extra code for figure 5.6 ~~~~~~~~~~~~~~~~~~
op <- par(mfrow = c(2,2))
pred.mass <- seq(0.1, 10,,100)
plot(pred.mass, out5$mean$pred[,1], type = 'l', xlab = "Mass (kg)",
    ylab = "Probability", ylim = c(0,1), frame = FALSE, xlim = c(0, 10),
    main = 'psi1', las = 1, lty = 2)
polygon(x = c(pred.mass, rev(pred.mass)),
    y = c(out5$q2.5$pred[,1], rev(out5$q97.5$pred[,1])),
    col = "gray90", border = "gray90")
lines(pred.mass, out5$mean$pred[,1], col = "blue", lwd = 2, lty = 2)
curve(plogis(0 - 1 * log(x)), 0.1, 10, lwd = 2, col = 'red', add = TRUE)

plot(pred.mass, out5$mean$pred[,2], type = 'l', xlab = "Mass (kg)",
    ylab = "Probability", ylim = c(0,1), frame = FALSE, xlim = c(0, 10),
    main = 'phi', las = 1, lty = 2)
polygon(x = c(pred.mass, rev(pred.mass)),
    y = c(out5$q2.5$pred[,2], rev(out5$q97.5$pred[,2])),
    col = "gray90", border = "gray90")
lines(pred.mass, out5$mean$pred[,2], col = "blue", lwd = 2, lty = 2)
curve(plogis(1+ 1 * log(x)), 0.1, 10, lwd = 2, col = 'red', add = TRUE)

plot(pred.mass, out5$mean$pred[,3], type = 'l', xlab = "Mass (kg)",
    ylab = "Probability", ylim = c(0,1), frame = FALSE, xlim = c(0, 10),
    main = 'gamma', las = 1, lty = 2)
polygon(x = c(pred.mass, rev(pred.mass)),
    y = c(out5$q2.5$pred[,3], rev(out5$q97.5$pred[,3])),
    col = "gray90", border = "gray90")
lines(pred.mass, out5$mean$pred[,3], col = "blue", lwd = 2, lty = 2)
curve(plogis(-1 - 1 * log(x)), 0.1, 10, lwd = 2, col = 'red', add = TRUE)

plot(pred.mass, out5$mean$pred[,4], type = 'l', xlab = "Mass (kg)",
    ylab = "Probability", ylim = c(0,1), frame = FALSE, xlim = c(0, 10),
    main = 'p', las = 1, lty = 2)
polygon(x = c(pred.mass, rev(pred.mass)),
    y = c(out5$q2.5$pred[,4], rev(out5$q97.5$pred[,4])),
    col = "gray90", border = "gray90")
lines(pred.mass, out5$mean$pred[,4], col = "blue", lwd = 2, lty = 2)
curve(plogis(0 - 0.5 * log(x)), 0.1, 10, lwd = 2, col = 'red', add = TRUE)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
