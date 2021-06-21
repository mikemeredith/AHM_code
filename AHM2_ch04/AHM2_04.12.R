#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

library(AHMbook)
library(jagsUI)

# 4.12 A demographic dynamic occupancy model
# ==========================================

simDemoDynocc(nsites = 100, nyears = 10, nvisit = 5, psi1 = 0.6,
    range.phi = c(0.2, 0.9), range.r = c(0, 0.4), range.p = c(0.1, 0.9),
    show.plot = TRUE)

# Some functionality of function simDemoDynocc
str(data <- simDemoDynocc(psi1 = 1))        # All sites initially occupied
str(data <- simDemoDynocc(nsites = 1000))            # Plenty more sites
str(data <- simDemoDynocc(nyears = 100))             # Plenty more years
str(data <- simDemoDynocc(nvisit = 20))              # Plenty more visits
str(data <- simDemoDynocc(range.phi = c(0.8, 0.8)))  # Constant survival
str(data <- simDemoDynocc(range.phi = c(0.2,0.3), range.r = c(0,0.2))) # Decline
str(data <- simDemoDynocc(range.phi = c(0.8,1), range.r = c(0.5,0.7))) # Increase
str(data <- simDemoDynocc(nvisit = 1))               # Single visit
str(data <- simDemoDynocc(range.p = c(1,1)))         # Perfect detection

# Generate a data set with year-specific parameters
set.seed(24)
str(data <- simDemoDynocc(psi1 = 0.6, nsites = 100, nyears = 20, nvisit = 5,
    range.phi = c(0.1, 0.9), range.r = c(0, 0.5), range.p = c(0.1, 0.9)))

# Bundle and summarize data
str(bdata <- list(y = data$y, nsites = data$nsites, nyears = data$nyears,
    nvisit = data$nvisit, first = data$f) )
# List of 5
# $ y     : int [1:100, 1:5, 1:20] 1 1 1 0 1 1 1 0 1 1 ...
# $ nsites: num 100
# $ nyears: num 20
# $ nvisit: num 5
# $ first : num [1:100] 1 1 1 5 1 1 1 8 1 1 ...

# Specify model in BUGS language
cat(file = "DemoDynocc1.txt", "
model {

  # Priors
  for(t in 1:(nyears-1)){
    phi[t] ~ dunif(0,1)
    r[t] ~ dunif(0,1)
    p[t] ~ dunif(0,1) # only nyears-1 p params in conditional model !
  }

  # Likelihood
  # Alive/dead process specified conditional on first detection
  for(i in 1:nsites){
    z[i, first[i]] ~ dbern(1) # Condition on year of first detection
    for(t in (first[i]+1):nyears) {
      z[i, t] ~ dbern(z[i,t-1] * (phi[t-1] + (1-phi[t-1]) * r[t-1]) +
          (1-z[i,t-1]) * r[t-1])
    }
  }

  # Observations conditional on Alive/dead process
  for(i in 1:nsites){
    for(t in (first[i]+1):nyears) {
      for(j in 1:nvisit){
        y[i,j,t] ~ dbern(z[i,t] * p[t-1])
      }
    }
  }
}
")

# Initial values
zst <- zinit(apply(bdata$y, c(1,3), max))
inits <- function() list(z = zst)

# Parameters monitored
params <- c("phi", "r", "p")

# MCMC settings
na <- 1000 ; ni <- 6000 ; nt <- 4 ; nb <- 2000 ; nc <- 3

# Call JAGS (ART 2 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "DemoDynocc1.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out1, layout=c(2,2))
print(out1, dig = 2)

# ~~~~~~~~~ extra code for figure 4.35 ~~~~~~~~~~~~~~~~~~
# Compare estimates and truth
op <- par(mfrow = c(1, 3))
lim <- c(0,1)
plot(data$phi, out1$mean$phi, xlab = "True phi", ylab = "Estimated phi",
    main = "Survival probability", pch = 16, xlim = lim, ylim = lim,
    frame = FALSE)
segments(data$phi, out1$q2.5$phi, data$phi, out1$q97.5$phi, lwd = 2)
abline(0, 1, col = 'red', lwd = 2)
abline(lm(out1$mean$phi ~ data$phi), col =  'blue', lwd = 2)

plot(data$r, out1$mean$r, xlab = "True r", ylab = "Estimated r",
    main = "Recruitment probability", pch = 16, xlim = c(0, 0.5),
    ylim = c(0, 0.65), frame = FALSE)
segments(data$r, out1$q2.5$r, data$r, out1$q97.5$r, lwd = 2)
abline(0, 1, col = 'red', lwd = 2)
abline(lm(out1$mean$r ~ data$r), col =  'blue', lwd = 2)

plot(data$p[-1], out1$mean$p, xlab = "True p", ylab = "Estimated p",
    main = "Detection probability", pch = 16, xlim = lim, ylim = lim, frame = FALSE)
segments(data$p[-1], out1$q2.5$p, data$p[-1], out1$q97.5$p, lwd = 2)
abline(0, 1, col = 'red', lwd = 2)
abline(lm(out1$mean$p ~ data$p[-1]), col = 'blue', lwd = 2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data loading from AHMbook and preparation
data(FrenchPeregrines)
fp <- FrenchPeregrines
y <- as.matrix(fp[,4:56])
f <- apply(y, 1, function(x) min(which(x!=0)))
f[f == 'Inf'] <- ncol(y)

# Bundle and summarize data
str(bdata <- list(y = y, nsites = nrow(y), nyears = ncol(y), first = f) )
# List of 4
# $ y     : int [1:284, 1:53] 1 NA NA NA NA NA NA NA 1 1 ...
# $ nsites: int 284
# $ nyears: int 53
# $ first : int [1:284] 1 39 35 12 32 24 22 23 1 1 ...

# Specify model in BUGS language
cat(file = "DemoDynocc2.txt", "
model {

  # Priors
  phi[1] ~ dunif(0, 1) # Survival in first interval
  r[1] ~ dunif(0, 1) # Recruitment in first interval
  lphi[1] <- logit(phi[1])
  lr[1] <- logit(r[1])
  for(t in 2:(nyears-1)){
    logit(phi[t]) <- lphi[t]
    logit(r[t]) <- lr[t]
  }

  # Random-walk smoother
  for (t in 2:(nyears-1)){ # Survival and recruitment in later intervals
    lphi[t] ~ dnorm(lphi[t-1], tau.lphi)
    lr[t] ~ dnorm(lr[t-1], tau.lr)
  }
  tau.lphi <- pow(sd.lphi,-2) # Hyperpriors for variances
  sd.lphi ~ dunif(0, 1)
  tau.lr <- pow(sd.lr,-2)
  sd.lr ~ dunif(0, 1)

  # Likelihood
  # Alive/dead process specified conditional on first detection
  for(i in 1:nsites){
    y[i, first[i]] ~ dbern(1) # Condition on year of first detection
    for(t in (first[i]+1):nyears) {
      y[i, t] ~ dbern(y[i,t-1] * (phi[t-1] + (1-phi[t-1]) * r[t-1]) +
          (1-y[i,t-1]) * r[t-1])
    }
  }
}
")

# Initial values
inits <- function(){ list(sd.lphi = runif(1), sd.lr = runif(1))}

# Parameters monitored
params <- c("phi", "r", "sd.lphi", "sd.lr")

# MCMC settings
na <- 1000 ; ni <- 12000 ; nt <- 6 ; nb <- 6000 ; nc <- 3

# Call JAGS (ART 3 min), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "DemoDynocc2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out2, layout=c(2,2))
print(out2, dig = 2) # not shown

# ~~~~~~~~~ extra code for figure 4.36 ~~~~~~~~~~~~~~~~~
# Visualize 'local survival' and recruitment
plot(1964:2015, out2$mean$phi, xlab = 'Year', ylab = 'Probability',
    type = 'b', pch = 1, ylim = c(0,1), frame = FALSE)
segments(1964:2015, out2$q2.5$phi, 1964:2015, out2$q97.5$phi)
points(1964:2015, out2$mean$r, type = 'b', pch = 16)
segments(1964:2015, out2$q2.5$r, 1964:2015, out2$q97.5$r)
legend('bottomright', c('Local survival', 'Recruitment'), pch = c(1, 16),
    bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
