#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# Expected run time for this script: 3.6 hrs

library(AHMbook)
library(jagsUI)

# 4.13 Accounting for temporary emigration and modeling
#   phenologies using occupancy: estimation of arrival and
#   departure in insects or migratory animals
# ==========================================

# Read in Marbled White data and do some data management (not all shown -- see website)
data(SwissMarbledWhite)       # from AHMbook
str(dat <- SwissMarbledWhite) # rename and look at data overview

# Data preparation
y <- as.matrix(dat[,14:24])    # Grab detection/nondetection data
DATE <- as.matrix(dat[,3:13])  # Grab survey dates
for(t in 1:11) {               # Mean-impute date (but don't transform)
  DATE[is.na(DATE[,t]),t] <- mean(DATE[,t], na.rm=TRUE)
}
year <- dat$year
nsites <- length(unique(dat$site))
nyears <- length(unique(dat$year))
nsurveys <- ncol(y)
nobs <- nrow(y)

# Bundle and summarize data
# Use year as factor and yr as regressor
str(bdata <- list(y = y, DATE = DATE, year = year-1997, yr = year-2004,
    site = dat$site, nobs = nobs, nsites = nsites, nyears = nyears,
    nsurveys = nsurveys))
# List of 9
# $ y       : int [1:1337, 1:11] 0 0 0 0 0 0 0 0 0 0 ...
# $ DATE    : num [1:1337, 1:11] 23 23 23 23 23 23 23 23 23 25 ...
# $ year    : num [1:1337] 1 1 1 1 1 1 1 1 1 1 ...
# $ yr      : num [1:1337] -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 ...
# $ site    : int [1:1337] 1 2 3 4 5 6 7 8 9 10 ...
# $ nobs    : int 1337
# $ nsites  : int 519
# $ nyears  : int 13
# $ nsurveys: int 11

cat(file = "PhenoOcc.txt", "
model {

  # Linear model for annual site occupancy (psi) with its priors
  for(i in 1:nobs) {
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- min(10, max(-10, lpsi[i]))
    lpsi[i] <- beta.lpsi[1] + beta.lpsi[2] * yr[i] + eps.lpsi[year[i]]
  }
  for(t in 1:nyears){
    eps.lpsi[t] ~ dnorm(0, tau.lpsi)
  }

  # Priors for occupancy
  beta.lpsi[1] <- logit(mean.psi)
  mean.psi ~ dunif(0, 1)
  beta.lpsi[2] ~ dnorm(0, 0.1)
  tau.lpsi <- pow(sigma.lpsi, -2)
  sigma.lpsi ~ dnorm(0, 0.1)I(0.01,)

  # Logit-linear model of detection on year
  for(t in 1:nyears) {
    for(j in 1:nsurveys) {
      logit(p[t,j]) <- min(10, max(-10, lp[t,j]))
      lp[t,j] <- beta.lp[1,j] + beta.lp[2,j]* (t-7) # year centered
    }
  }

  # Priors for detection
  for (j in 1:nsurveys) { # visit-specific regression coefs on year
    beta.lp[1,j] ~ dnorm(mu.lp1, tau.lp1)
    beta.lp[2,j] ~ dnorm(mu.lp2, tau.lp2)
  }
  mu.lp1 ~ dnorm(0, 0.1)
  mu.lp2 ~ dnorm(0, 1)
  tau.lp1 <- pow(sigma.lp1, -2)
  tau.lp2 <- pow(sigma.lp2, -2)
  sigma.lp1 ~ dnorm(0, 1)I(0.01,)
  sigma.lp2 ~ dnorm(0, 1)I(0.01,)
  # curve(dnorm(x, 0, sqrt(1)), 0, 10) # how does this look like ?

  # Linear regression of arrival date (arr) on year
  for(i in 1:nobs) {
    arr[i] ~ dnorm(mu.arr1[i], tau.arr)
    mu.arr1[i] <- min(200, max(30, mu.arr[i]))
    mu.arr[i] <- beta.arr[1] + beta.arr[2] * yr[i]
  }

  # Priors for arrival model
  beta.arr[1] ~ dnorm(90, 0.1)
  beta.arr[2] ~ dnorm(0, 1)
  tau.arr <- pow(sigma.arr, -2)
  sigma.arr ~ dnorm(0, 0.01)I(0.01,)
  # curve(dnorm(x, 0, sqrt(100)), 0, 20) # how does this look like ?

  # Linear regression of departure date (dep) on year
  for(i in 1:nobs) {
    dep[i] ~ dnorm(mu.dep1[i], tau.dep)
    mu.dep1[i] <- min(500, max(0, mu.dep[i]))
    mu.dep[i] <- beta.dep[1] + beta.dep[2] * yr[i]
  }

  # Priors for departure model
  beta.dep[1] ~ dnorm(120, 0.1)
  beta.dep[2] ~ dnorm(0, 1)
  tau.dep <- pow(sigma.dep, -2)
  sigma.dep ~ dnorm(0, 0.01)I(0.01,)

  # Model for the observed data
  for(i in 1:nobs) {
    for(j in 1:nsurveys) {
      y[i,j] ~ dbern(mu1[i,j])
      mu1[i,j] <- min(0.99, max(0.01, mu[i,j]))
      mu[i,j] <- z[i] * step(DATE[i,j] - arr[i]) * step(dep[i] -
          DATE[i,j]) * p[year[i],j]
    }
  }

  # Derived quantities
  # Average occupancy per year
  for(t in 1:nyears){
    for(s in 1:nsites){
      logit(tmp[s, t]) <- beta.lpsi[1] + beta.lpsi[2] * (t-7) + eps.lpsi[t]
    }
    psi.pred[t] <- mean(tmp[,t])
  }

  # Average detection per year and visit
  for(t in 1:nyears){
    for(j in 1:nsurveys){
      logit(p.pred[t,j]) <- beta.lp[1,j] + beta.lp[2,j] * (t-7)
    }
  }

  # Average arrival and departure time per year and length of flight period
  for(t in 1:nyears){
    arr.pred[t] <- beta.arr[1] + beta.arr[2] * (t-7)
    dep.pred[t] <- beta.dep[1] + beta.dep[2] * (t-7)
    fp.pred[t] <- dep.pred[t] - arr.pred[t]
  }
}
")

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)
inits <- function() {list(z = zst)}

# Parameters monitored
params <- c('mean.psi', 'beta.lpsi', 'sigma.lpsi', 'beta.lp', 'mu.lp1',
    'mu.lp2', 'sigma.lp1', 'sigma.lp2', 'beta.arr', 'sigma.arr', 'beta.dep',
    'sigma.dep', 'psi.pred', 'p.pred', 'arr.pred', 'dep.pred', 'fp.pred')

# MCMC settings
# na <- 5000 ; ni <- 100000 ; nt <- 80 ; nb <- 20000 ; nc <- 3
na <- 5000 ; ni <- 10000 ; nt <- 8 ; nb <- 2000 ; nc <- 3 # ~~~ for testing

# Call JAGS (ART 400 min), check convergence and summarize posteriors
out <- jags(bdata, inits, params, "PhenoOcc.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out)
print(out, dig = 2)

# ~~~~~~~~~ code for figure 4.37 ~~~~~~~~~~~~~~~~
# Visualize results as trajectories over 1998-2010
op <- par(mfrow = c(2,2))
# Occupancy
plot(1998:2010, out$summary[37:49,1], xlab = 'Year',
    ylab = 'Occupancy probability',
    ylim = c(0.1, 0.4), type = 'l', lwd = 3, frame = FALSE, las = 1)
polygon(c(1998:2010, rev(1998:2010)), c(out$summary[37:49,3],
    rev(out$summary[37:49, 7])), col = 'grey', border = 'grey')
lines(1998:2010, out$summary[37:49,1], lwd = 3)

# Visit-specific detection probability
matplot(1998:2010, out$mean$p.pred, xlab = 'Year', ylab = 'Detection probability',
    ylim = c(0.5, 1), type = 'l', lwd = 2, frame = FALSE, las = 1, lty = 1, col = 'black')

# Arrival and departure time
plot(1998:2010, out$summary[193:205,1], xlab = 'Year', ylab = 'Date',
    frame = FALSE, las = 1, ylim = c(80, 130), pch = 16)
segments(1998:2010, out$summary[193:205,1], 1998:2010, out$summary[206:218,1],
    lwd = 2, col = 'grey', lend = 'butt')
points(1998:2010, out$summary[193:205,1], pch = 16)
segments(1998:2010, out$summary[193:205,3], 1998:2010, out$summary[193:205,7])
points(1998:2010, out$summary[206:218,1], pch = 16)
segments(1998:2010, out$summary[206:218,3], 1998:2010, out$summary[206:218,7])

# Length of flight period
plot(1998:2010, out$summary[219:231,1], xlab = 'Year', ylab = 'Length (days)',
    frame = FALSE, las = 1, ylim = c(20, 40), type='n')
polygon(c(1998:2010, rev(1998:2010)), c(out$summary[219:231,3],
    rev(out$summary[219:231,7])), col = 'grey', border = 'grey')
lines(1998:2010, out$summary[219:231,1], lwd = 2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
