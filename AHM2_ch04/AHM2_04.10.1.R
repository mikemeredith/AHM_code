#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 15 mins
# Run time with the full number of iterations: 90 mins

library(AHMbook)
library(jagsUI)

# 4.10 Analysis of citizen-science data using occupancy models
# ============================================================

# 4.10.1 Effects of trends in the magnitude of unmodeled detection heterogeneity
# ------------------------------------------------------------------------------

# Generate data set with increase in site-level heterogeneity
set.seed(1)
str(data <- simDynocc(nsites = 250, nyears = 20, nsurveys = 3, mean.psi1 = 0.6,
    range.p = c(0.5, 0.5), range.phi = c(0.8, 0.8), range.gamma = c(0.3, 0.3),
    trend.sd.site = c(0, 2))) # library(AHMbook)

# ~~~~~~~ extra code from MS dated 2019-01-04 ~~~~~~~~~~~~~~

# Traditional, non-heterogeneity occupancy model
# ''''''''''''''''''''''''''''''''''''''''''''''
# Bundle data
str(bdata <- list(y = data$y, nsite = dim(data$y)[1], nsurvey = dim(data$y)[2],
    nyear = dim(data$y)[3]))

# Specify model in BUGS language
cat(file = "dynocc.txt", "
model {

  # Specify priors
  psi1 ~ dunif(0, 1)
  for (t in 1:(nyear-1)){
    phi[t] ~ dunif(0, 1)
    gamma[t] ~ dunif(0, 1)
    p[t] ~ dunif(0, 1)
  }
  p[nyear] ~ dunif(0, 1)

  # Ecological submodel
  for (i in 1:nsite){
    z[i,1] ~ dbern(psi1)
    for (t in 2:nyear){
      z[i,t] ~ dbern(z[i,t-1]*phi[t-1] + (1-z[i,t-1])*gamma[t-1])
    }
  }

  # Observation model
  for (i in 1:nsite){
    for (j in 1:nsurvey){
      for (t in 1:nyear){
        y[i,j,t] ~ dbern(z[i,t] * p[t])
      }
    }
  }

  # Compute population and sample occupancy
  psi[1] <- psi1                         # Population occupancy
  psi.fs[1] <- sum(z[1:nsite,1]) / 250   # Sample occupancy
  for (t in 2:nyear){
    psi[t] <- psi[t-1]*phi[t-1] + (1-psi[t-1])*gamma[t-1]
    psi.fs[t] <- sum(z[1:nsite,t]) /  250
  }
}
")

# Initial values
inits <- function(){ list(z = apply(data$y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "psi.fs", "phi", "gamma", "p")

# MCMC settings
na <- 1000  ;  ni <- 6000  ;  nb <- 4000  ;  nt <- 2  ;  nc <- 3 # 5 mins

# Call JAGS from R (ART 25 min)
out0 <- jags(bdata, inits, params, "dynocc.txt",
    n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(4,4))  # ~~~ replace with 'layout' argument
traceplot(out0, layout=c(4,4))
print(out0, 2)
# ~~~~~~~~~~~~~~~~~~~~~~ end of extra code ~~~~~~~~~~~~~~~~~~~

# Heterogeneity occupancy model
# '''''''''''''''''''''''''''''
# Bundle data
str(bdata <- list(y = data$y, nsites = dim(data$y)[1],
    nsurveys = dim(data$y)[2], nyears = dim(data$y)[3]))
# List of 4
# $ y       : int [1:250, 1:3, 1:20] 0 0 1 1 0 0 1 0 0 0 ...
# $ nsites  : int 250
# $ nsurveys: int 3
# $ nyears  : int 20

# Specify model in BUGS language
cat(file = "dynoccH.txt", "
model {

  # Specify priors
  psi1 ~ dunif(0, 1)
  for (t in 1:(nyears-1)){
    phi[t] ~ dunif(0, 1)
    gamma[t] ~ dunif(0, 1)
  }
  for (t in 1:nyears){
    lp[t] <- logit(mean.p[t])
    mean.p[t] ~ dunif(0, 1)
  }

  # Random effects priors
  for (t in 1:nyears){
    for (i in 1:nsites){
      eps[i,t] ~ dnorm(0, tau.eps[t])
    }
    tau.eps[t] <- pow(sd.eps[t], -2)
    sd.eps[t] ~ dunif(0, 10) # Note different variance in every year
  }

  # Ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      z[i,t] ~ dbern(z[i,t-1]*phi[t-1] + (1-z[i,t-1])*gamma[t-1])
    }
  }

  # Observation model
  for (i in 1:nsites){
    for (j in 1:nsurveys){
      for (t in 1:nyears){
        logit(p[i,j,t]) <- lp[t] + eps[i,t] # time + site.survey effects
        y[i,j,t] ~ dbern(z[i,t] * p[i,j,t])
      }
    }
  }

  # Compute population and sample occupancy
  psi[1] <- psi1 # Population occupancy
  psi.fs[1] <- sum(z[1:nsites,1]) / 250 # Sample occupancy
  for (t in 2:nyears){
    psi[t] <- psi[t-1]*phi[t-1] + (1-psi[t-1])*gamma[t-1]
    psi.fs[t] <- sum(z[1:nsites,t]) / 250
  }
}
")

# Initial values
inits <- function(){ list(z = apply(data$y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "psi.fs", "phi", "gamma", "mean.p", "sd.eps")

# MCMC settings
# na <- 5000 ; ni <- 100000 ; nb <- 20000 ; nt <- 80 ; nc <- 3
na <- 5000 ; ni <- 10000 ; nb <- 2000 ; nt <- 8 ; nc <- 3  # ~~~~ testing, 15 mins

# Call JAGS (ART 163 min), check convergence and summarize posteriors
out <- jags(bdata, inits, params, "dynoccH.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(4,4))  # ~~~ replace with 'layout' argument
traceplot(out, layout=c(4,4))
print(out, 2) # not shown

save(out0, out,  file="AHM2_10.1_output.RData")

# ~~~~~~~~~ code for Figure 4.22 ~~~~~~~~~~~~~
op <- par(mfrow = c(1, 3))

# (1) Plot realised and apparent proportion of occupied sites
year <- 1:data$nyear
plot(year, data$psi.fs, type = "o", xlab = "Year", ylab = "Occupancy",
    col = "red", pch=16,
    ylim = c(0.35,0.8), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
lines(year, data$psi.app, type = "o", col = "black", lwd = 2)
legend('topleft', c('True (with linear reg. line)',
    'Apparent (with linear reg. line)'),
    col = c('red', 'black'), lwd = 2, bty = 'n', pch=16)
abline(lm(data$psi.fs ~ year), col = "red", lwd = 1)
abline(lm(data$psi.app ~ year), col = "black", lwd = 1)

# (2) Estimates from the two models
plot(year-0.1, out0$mean$psi, type = "o", xlab = "Year", ylab = "Occupancy",
    col = "blue", ylim = c(0.35,0.8), lwd = 2, lty = 1, frame.plot = FALSE,
    las = 1, main = '')
segments(1:data$nyear-0.1, out0$q2.5$psi, 1:data$nyear-0.1, out0$q97.5$psi,
    col = "blue", lwd=1)
points(1:data$nyear+0.1, out$mean$psi, col = "brown", type = "o",
    lwd = 2, pch = 16)
segments(1:data$nyear+0.1, out$q2.5$psi, 1:data$nyear+0.1, out$q97.5$psi,
    col = "brown", lwd=1)
legend('bottomleft', c('Heterogeneity model', 'Simple model'),
    col = c('brown', 'blue'), lwd = 2, bty = 'n', pch=16)

# (3) Estimated magnitude of detection heterogeneity
plot(year, out$mean$sd.eps, type = "p", xlab = "Year",
    ylab = "SD of Site heterogeneity", col = "blue",
    ylim = c(0,3.5), pch = 16, frame = FALSE, las=1)
segments(1:data$nyear, out$q2.5$sd.eps, 1:data$nyear, out$q97.5$sd.eps,
    col = "blue")
points(1:data$nyear, data$sd.site, col = "red", pch = 16)
legend('topleft', pch = 16, col = c('blue', 'red'),
    legend = c('Estimated', 'True'), bty = 'n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
