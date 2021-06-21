#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6 : MULTISTATE OCCUPANCY MODELS
# =======================================
# Code from proofs dated 2020-08-19

library(jagsUI)
library(unmarked)

# 6.3 Simulation and analysis of the simplest static and dynamic multistate models
# ================================================================================

# Pick values for 11 parameters
# Parameters for initial conditions (Omega)
psi <- 0.8                # Expected proportion of occupied sites
r <- 0.5                  # Exp. proportion of sites (among occupied) with a pair

# Parameters for state transition matrix (Phi)
phi <- c(0.2, 0.8, 0.9)   # Prob. to become occupied
rho <- c(0.5, 0.25, 0.89) # Prob. to get a pair if become occupied

# Parameters for observation matrix (Theta)
p2 <- 0.5                 # Detection probability of site with single bird
p32 <- 0.2                # Classification probability of site with pair as single bird
p33 <- 0.6                # Classification probability of site with pair as pair

# Assemble initial state vector and the two matrices
# Populate initial state probability vector (Omega)
Omega <- c(1-psi, psi*(1-r), psi*r)

# Population state transition probability matrix (Phi)
Phi <- matrix(
    c(1-phi[1], phi[1]*(1-rho[1]), phi[1]*rho[1],
      1-phi[2], phi[2]*(1-rho[2]), phi[2]*rho[2],
      1-phi[3], phi[3]*(1-rho[3]), phi[3]*rho[3]), ncol = 3, byrow = TRUE)

# Populate observation probability matrix (Theta)
Theta <- matrix(
  c(1, 0, 0,
    1-p2, p2, 0,
    1-p32- p33, p32, p33), ncol = 3, byrow = TRUE)

# Inspect the three arrays
Omega # Initial state vector
Phi   # State transition matrix
Theta # Observation matrix
# > Omega                     # Initial state vector
# [1] 0.2 0.4 0.4
# > Phi                       # State transition matrix
#      [,1]  [,2]  [,3]
# [1,]  0.8 0.100 0.100
# [2,]  0.2 0.600 0.200
# [3,]  0.1 0.099 0.801
# > Theta                     # Observation matrix
#      [,1] [,2] [,3]
# [1,]  1.0  0.0  0.0
# [2,]  0.5  0.5  0.0
# [3,]  0.2  0.2  0.6

# Pick sample sizes (note use of names instead of letters)
nsites <- 100                 # denoted "M" above
nsurveys <- 3                 # ... "J" ...
nyears <- 5                   # ... "T" ...

# Generate structures for latent states (z) and for observations (y)
z <- array(NA, dim = c(nsites, nyears))
y <- array(NA, dim = c(nsites, nsurveys, nyears))

# Draw initial states in year 1 using initial state vector Omega
get1 <- function(x) which(x==1)   # Get positition of the sole one (1)
set.seed(1)
z[,1] <- apply(rmultinom(nsites, 1, Omega), 2, get1)

# Draw states in following years (2-10) using state transition matrix Phi
for(i in 1:nsites){
  for(t in 2:nyears){
    z[i,t] <- get1(rmultinom(1, 1, Phi[z[i,t-1],]))
  }
}

# Draw observations (all years) using observation matrix Theta
for(i in 1:nsites){
  for(t in 1:nyears){
    y[i,,t] <- apply(rmultinom(3, 1, Theta[z[i,t],]), 2, get1)
  }
}

# Latent states for first 6 sites
head(z)
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    3    3    3    2    2
# [2,]    2    2    1    1    1
# [3,]    2    3    3    3    3
# [4,]    1    1    1    1    1
# [5,]    2    2    1    1    1
# [6,]    3    1    3    2    1

# Observed data for site 1 (dimension is survey x year)
y[1, ,]
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    3    1    3    2    2
# [2,]    3    2    3    1    1
# [3,]    1    3    3    1    2


# 6.3.1 A static single-season model
# ----------------------------------

# Grab data from year 1
str(y1 <- y[,,1])

# Tabulate observed data for each site
ttab <- array(0, dim = c(nsites, 3)) # sites by number of states
colnames(ttab) <- c('no.bird.detected', 'single.detected', 'pair.detected')
for(i in 1:nsites){
  tt <- table(y1[i,])
  ttab[i, as.numeric(names(tt))] <- tt
}

# Compare true state and summary of observed states
data.frame('True state' = z[,1], ttab)
#     True.state no.bird.detected single.detected pair.detected
# 1            3                1               0             2
# 2            2                1               2             0
# 3            2                1               2             0
# ....
# 98           2                0               3             0
# 99           3                0               0             3
# 100          1                3               0             0

# Bundle data
str(bdata <- list(y = y1, nsites = nrow(y1), nsurveys = ncol(y1)))
# List of 3
# $ y       : int [1:100, 1:3] 3 2 2 1 2 3 2 2 2 1 ...
# $ nsites  : int 100
# $ nsurveys: int 3

# Specify model in BUGS language
cat(file = "static1.txt", "
model {

  # Priors
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)
  p2 ~ dunif(0, 1)

  # Multinomial logit link for observation model for state 3 (= pair)
  lp32 ~ dnorm(0, 0.001)
  lp33 ~ dnorm(0, 0.001)
  p32 <- exp(lp32) / (1 + exp(lp32) + exp(lp33))
  p33 <- exp(lp33) / (1 + exp(lp32) + exp(lp33))
  p31 <- 1-p32-p33                     # Nondetection prob for pairs by difference

  # Alternative: induce Dirichlet prior for p3
  # for (s in 1:3) {
  #   beta[s] ~ dgamma(1, 1)
  #   p3[s] <- beta[s] / sum(beta[])
  # }

  # Define initial state vector (Omega)
  Omega[1] <- 1 - psi                  # Prob. of non-occupation
  Omega[2] <- psi * (1-r)              # Prob. of occupancy (w/ single bird)
  Omega[3] <- psi * r                  # Prob. of occupancy (with pair)

  # Define observation matrix (Theta)
  # Order of indices: true state, observed state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-p2
  Theta[2,2] <- p2
  Theta[2,3] <- 0
  Theta[3,1] <- p31                    # = 1-p32-p33 as per prior section
  Theta[3,2] <- p32
  Theta[3,3] <- p33

  # State-space likelihood
  # State equation: model of true states (z)
  for (i in 1:nsites){
    z[i] ~ dcat(Omega[])
  }

  # Observation equation
  for (i in 1:nsites){
    for (j in 1:nsurveys){
      y[i,j] ~ dcat(Theta[z[i],])
    }
  }

  # Derived quantities
  for (i in 1:nsites){
    occ1[i] <- equals(z[i], 1)
    occ2[i] <- equals(z[i], 2)
    occ3[i] <- equals(z[i], 3)
  }
  n.occ[1] <- sum(occ1[]) # Sites in state 1
  n.occ[2] <- sum(occ2[]) # Sites in state 2
  n.occ[3] <- sum(occ3[]) # Sites in state 3
}
")

# Initial values
zst <- rep(3, nrow(bdata$y)) # Initialize at highest possible state
inits <- function(){list(z = zst)}

# Parameters monitored (could add "z")
params <- c("psi", "r", "p2", "p31", "p32", "p33", "Omega", "Theta", "n.occ")

# MCMC settings
na <- 1000 ; ni <- 2000 ; nt <- 2 ; nb <- 1000 ; nc <- 3

# Call JAGS, check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "static1.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
print(out1, 3)
#             mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# psi        0.866 0.053 0.762 0.866 0.966    FALSE 1 1.001  1500
# r          0.429 0.060 0.310 0.430 0.549    FALSE 1 1.001  1500
# p2         0.489 0.056 0.381 0.490 0.592    FALSE 1 1.005   445
# p31        0.212 0.043 0.135 0.212 0.300    FALSE 1 1.000  1500
# p32        0.171 0.040 0.098 0.169 0.253    FALSE 1 1.002  1222
# p33        0.616 0.054 0.508 0.618 0.716    FALSE 1 1.000  1500
# Omega[1]   0.134 0.053 0.034 0.134 0.238    FALSE 1 1.001  1500
# Omega[2]   0.495 0.067 0.371 0.491 0.634    FALSE 1 1.000  1500
# Omega[3]   0.370 0.052 0.273 0.372 0.475    FALSE 1 1.001  1500
# Theta[1,1] 1.000 0.000 1.000 1.000 1.000    FALSE 1    NA     1
# .....

cbind('truth' = table(z[,1]), 'estimates' = out1$summary[19:21,c(1,3,7)])
#   truth     mean 2.5% 97.5%
# 1    17 12.56067    3    18
# 2    43 50.18467   43    60
# 3    40 37.25467   35    42

# 6.3.2 A static multiseason model
# --------------------------------

# Bundle all data now
str(bdata <- list(y = y, nsites = dim(y)[1], nsurveys = dim(y)[2],
    nyears = dim(y)[3]))
# List of 4
# $ y       : int [1:100, 1:3, 1:5] 3 2 2 1 2 3 2 2 2 1 ...
# $ nsites  : int 100
# $ nsurveys: int 3
# $ nyears  : int 5

# Specify model in BUGS language
cat(file = "static2.txt", "
model {

  # Priors for each year
  for (t in 1:nyears){
    psi[t] ~ dunif(0, 1)
    r[t] ~ dunif(0, 1)
    p2[t] ~ dunif(0, 1)
    # Multinomial logit link for p3[1:3], i.e., p31, p32, p33
    lp32[t] ~ dnorm(0, 0.001)
    lp33[t] ~ dnorm(0, 0.001)
    p32[t] <- exp(lp32[t]) / (1 + exp(lp32[t]) + exp(lp33[t]))
    p33[t] <- exp(lp33[t]) / (1 + exp(lp32[t]) + exp(lp33[t]))
    p31[t] <- 1-p32[t]-p33[t]             # Nondetection prob for pairs by difference
  }

  # Alternative: Dirichlet prior for p3
  # for (t in 1:nyears){
  #   for (s in 1:3) {
  #     beta[t,s] ~ dgamma(1, 1)
  #     p3[t,s] <- beta[,ts]/sum(beta[t,])
  #   }
  # }

  # Define state vector for each year
  # Called Omega here; has different meaning in dynamic model
  for (t in 1:nyears){
    Omega[t,1] <- 1 - psi[t]              # Prob. of non-occupation
    Omega[t,2] <- psi[t] * (1-r[t])       # Prob. of occupancy (w/ single bird)
    Omega[t,3] <- psi[t] * r[t]           # Prob. of occupancy (with pair)
  }

  # Define observation matrix for each year
  # Order of indices: true state, year, observed state
  for (t in 1:nyears){
    Theta[1,t,1] <- 1
    Theta[1,t,2] <- 0
    Theta[1,t,3] <- 0
    Theta[2,t,1] <- 1-p2[t]
    Theta[2,t,2] <- p2[t]
    Theta[2,t,3] <- 0
    Theta[3,t,1] <- p31[t]                # = 1-p32[t]-p33[t]
    Theta[3,t,2] <- p32[t]
    Theta[3,t,3] <- p33[t]
  }

  # State-space likelihood
  # Define separate params of state and observation equation for each year
  for (t in 1:nyears){
    for (i in 1:nsites){
      z[i,t] ~ dcat(Omega[t,])            # State equation
      for (j in 1:nsurveys){
        y[i,j,t] ~ dcat(Theta[z[i,t],t,]) # Observation equation
      }
    }
  }

  # Derived quantities
  for (t in 1:nyears){
    for (i in 1:nsites){
      occ1[i,t] <- equals(z[i,t], 1)
      occ2[i,t] <- equals(z[i,t], 2)
      occ3[i,t] <- equals(z[i,t], 3)
    }
    n.occ[t,1] <- sum(occ1[,t])           # Sites in state 1
    n.occ[t,2] <- sum(occ2[,t])           # Sites in state 2
    n.occ[t,3] <- sum(occ3[,t])           # Sites in state 3
  }
}
")

# Initial values
zst <- array(3, dim = c(bdata$nsite, bdata$nyears))
inits <- function(){list(z = zst)}

# Parameters monitored (could add "z")
params <- c("psi", "r", "p2", "p31", "p32", "p33", "Omega", "Theta", "n.occ")

# MCMC settings
na <- 1000 ; ni <- 2000 ; nt <- 2 ; nb <- 1000 ; nc <- 3

# Call JAGS, check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "static2.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out2)
print(out2, 3) # not shown

true.sumZ <- array(NA, dim = c(3, 5), dimnames = list(c('unocc',
    'single', 'pair'), paste('year', 1:5, sep = '')))
for(t in 1:5){ # Compute true sum of sites in each state
  true.sumZ[,t] <- table(z[,t])
}
true.sumZ                                # True
t(out2$mean$n.occ)                       # Estimated (posterior means)


# 6.3.3 The dynamic multiseason model
# -----------------------------------

# Same data bundle as in last section
str(bdata <- list(y = y, nsites = dim(y)[1], nsurveys = dim(y)[2],
    nyears = dim(y)[3]))

# Specify model in BUGS language
cat(file = "dynamic1.txt", "
model {

  # Priors
  # Priors for parameters in initial state vector Omega
  psi ~ dunif(0,1)
  r ~ dunif(0,1)

  # Priors for params. of state transition mat. PhiMat (called Phi in text)
  for(s in 1:3){
    phi[s] ~ dunif(0,1)
    rho[s] ~ dunif(0,1)
  }

  # Priors for parameters in the observation matrix Theta
  p2 ~ dunif(0, 1)
  for (s in 1:3) {                       # # Dirichlet prior for p3
    beta[s] ~ dgamma(1, 1)
    p3[s] <- beta[s] / sum(beta[])
  }

  # Define initial state vector Omega (this is only for year 1)
  Omega[1] <- 1 - psi                    # Prob. of non-occupation
  Omega[2] <- psi * (1-r)                # Prob. of single bird
  Omega[3] <- psi * r                    # Prob. of pair

  # Define state transition matrix PhiMat (years 2:nyears)
  # Order of indices: state at t, state at t+1
  PhiMat[1, 1] <- 1 - phi[1]
  PhiMat[1, 2] <- phi[1] * (1 - rho[1])
  PhiMat[1, 3] <- phi[1] * rho[1]
  PhiMat[2, 1] <- 1 - phi[2]
  PhiMat[2, 2] <- phi[2] * (1 - rho[2])
  PhiMat[2, 3] <- phi[2] * rho[2]
  PhiMat[3, 1] <- 1 - phi[3]
  PhiMat[3, 2] <- phi[3] * (1 - rho[3])
  PhiMat[3, 3] <- phi[3] * rho[3]

  # Define observation matrix Theta
  # Order of indices: true state, observed state
  Theta[1, 1] <- 1
  Theta[1, 2] <- 0
  Theta[1, 3] <- 0
  Theta[2, 1] <- 1-p2
  Theta[2, 2] <- p2
  Theta[2, 3] <- 0
  Theta[3, 1] <- p3[1]
  Theta[3, 2] <- p3[2]
  Theta[3, 3] <- p3[3]

  # State-space likelihood
  # Initial state governed by vector Omega
  for (i in 1:nsites){
    z[i,1] ~ dcat(Omega[])
  }

  # State transitions governed by transition matrix PhiMat
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      z[i,t+1] ~ dcat(PhiMat[z[i,t],])
    }
  }

  # Observation equation with observation matrix Theta
  for(t in 1:nyears){
    for (i in 1:nsites){
      for (j in 1:nsurveys){
        y[i,j,t] ~ dcat(Theta[z[i,t],])
      }
    }
  }

  # Derived quantities
  for (t in 1:nyears){
    for (i in 1:nsites){
      occ1[i,t] <- equals(z[i,t], 1)
      occ2[i,t] <- equals(z[i,t], 2)
      occ3[i,t] <- equals(z[i,t], 3)
    }
    n.occ[t,1] <- sum(occ1[,t])          # Sites in state 1
    n.occ[t,2] <- sum(occ2[,t])          # Sites in state 2
    n.occ[t,3] <- sum(occ3[,t])          # Sites in state 3
  }
}
")

# Initial values
inits <- function(){list(z = array(3, dim = c(bdata$nsites,
    bdata$nyears)))}

# Parameters monitored
params <- c("psi", "r", "phi", "rho", "p2", "p3", "Omega", "PhiMat",
    "Theta", "n.occ", "z")

# MCMC settings
na <- 1000 ; ni <- 2000 ; nt <- 2 ; nb <- 1000 ; nc <- 3

# Call JAGS, check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "dynamic1.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out3)
print(out3, 3) # not shown

# ~~~~~ extra code for figures 6.1 and 6.2 ~~~~~~~~~~~~~~~~~~~~~~
# Figure 6.1: Image plot for latent states of these first 25 sites
# z transposed for greater viewing pleasure
mapPalette <- colorRampPalette(c("white", "black"))
op <- par(mfrow=1:2)
# left panel, true state
image(x = 1:nyears, y = 1:25, z = t(z[1:25,]), xlab = 'Year',
    ylab = 'Site', col = mapPalette(10), las = 1)
# right panel, estimates from the dynamic model
image(x = 1:5, y = 1:25, z = round(t(out3$mean$z[1:25,])), col = mapPalette(10),
    xlab = "Year", ylab = "Site")
par(op)

# Figure 6.2
op <- par(mfrow = c(1,2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
cx <- 2
matplot(1:5, t(true.sumZ[2:3,]), xlab = 'Year', ylab = 'Number of sites',
    type = 'b', lty = 1, pch = c(16,17), cex = cx, col = 'black',
    xlim = c(0.8, 5.2), ylim = c(5, 60), frame = F)
lines(1:5-0.1, out2$mean$n.occ[,2], type = 'b', lty = 2, pch = 1, cex = cx)
segments(1:5-0.1, out2$q2.5$n.occ[,2], 1:5-0.1, out2$q97.5$n.occ[,2])
lines(1:5+0.1, out2$mean$n.occ[,3], type = 'b', lty = 2, pch = 2, cex = cx)
segments(1:5+0.1, out2$q2.5$n.occ[,3], 1:5+0.1, out2$q97.5$n.occ[,3])

matplot(1:5, t(true.sumZ[2:3,]), xlab = 'Year', ylab = '',
    type = 'b', lty = 1, pch = c(16,17), cex = cx, col = 'black',
    xlim = c(0.8, 5.2), ylim = c(5, 60), frame = F)
lines(1:5-0.1, out3$mean$n.occ[,2], type = 'b', lty = 2, pch = 1, cex = cx)
segments(1:5-0.1, out3$q2.5$n.occ[,2], 1:5-0.1, out3$q97.5$n.occ[,2])
lines(1:5+0.1, out3$mean$n.occ[,3], type = 'b', lty = 2, pch = 2, cex = cx)
segments(1:5+0.1, out3$q2.5$n.occ[,3], 1:5+0.1, out3$q97.5$n.occ[,3])
legend('bottomleft', c('True singles', 'True pairs', 'Estimated singles',
    'Estimated pairs'), pch = c(16, 17, 1, 2), lty = c(1,1,2,2),
    bty = 'n', cex = 1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 6.3.4 Doing it in unmarked
# --------------------------
summary(umf <- unmarkedFrameOccuMS(y = y1 - 1))

# unmarkedFrame Object

# 100 sites
# Maximum number of observations per site: 3
# Mean number of observations per site: 3
# Number of primary survey periods: 1
# Number of secondary survey periods: 3
# Sites with at least one detection: 80

# Tabulation of y observations:
#   0  1  2
# 139 92 69

# Check unmarked has recognized correct number of states
umf@numStates
# [1] 3

# Occupancy formulas for state 1 and state 2
psiformulas <- c("~1", "~1")

# Detection formulas for p1, p2, and delta
detformulas <- c("~1", "~1", "~1")

# Fit static model with conditional binomial parameterization
fm1 <- occuMS(detformulas, psiformulas, data = umf, parameterization = "condbinom")
summary(fm1)

# Call:
# occuMS(detformulas = detformulas, psiformulas = psiformulas, data = umf,
# parameterization = "condbinom")

# Occupancy (logit-scale):
# Estimate SE z P(>|z|)
# psi (Intercept) 1.903 0.450 4.23 2.32e-05
# R (Intercept) -0.301 0.244 -1.24 2.16e-01

# Detection (logit-scale):
# Estimate SE z P(>|z|)
# p[1] (Intercept) -0.0447 0.224 -0.20 8.42e-01
# p[2] (Intercept) 1.3216 0.263 5.02 5.21e-07
# delta (Intercept) 1.3081 0.297 4.40 1.07e-05

# Compare to JAGS estimates (requires you still have out1 !)
jags_est <- out1$summary[1:6,1]

# Raw unmarked
unmarked_est_raw <- unmarked_est <- plogis(coef(fm1))

# Translate to JAGS parameterization (see above)
unmarked_est[4] <- unmarked_est_raw[4] * (1-unmarked_est_raw[5])
unmarked_est[5] <- unmarked_est_raw[4] * unmarked_est_raw[5]

# Good match!
round(cbind(JAGS=jags_est[-4], unmarked=unmarked_est), 3)
#      JAGS unmarked
# psi 0.866    0.870
# r   0.429    0.425
# p2  0.489    0.489
# p32 0.171    0.168
# p33 0.616    0.621

# add some dummy covariates into our static model
# '''''''''''''''''''''''''''''''''''''''''''''''
# Simulate a random site and observational covariate
set.seed(25)
site_covs <- data.frame(elevation = rnorm(nsites))
obs_covs <- data.frame(date = rnorm(length(y1)))

# Create unmarked frame (remember to subtract 1 from y1)
summary(umf2 <- unmarkedFrameOccuMS(y = y1-1, siteCovs = site_covs,
    obsCovs = obs_covs))

# unmarkedFrame Object

# 100 sites
# Maximum number of observations per site: 3
# Mean number of observations per site: 3
# Number of primary survey periods: 1
# Number of secondary survey periods: 3
# Sites with at least one detection: 80

# Tabulation of y observations:
#   0  1  2
# 139 92 69

# Site-level covariates:
# elevation
# Min. :-2.3449
# 1st Qu.:-1.0332
# Median :-0.1420
# Mean :-0.1822
# 3rd Qu.: 0.3572
# Max. : 2.3678

# Observation-level covariates:
# date
# Min. :-2.744085
# 1st Qu.:-0.696238
# Median :-0.002968
# Mean : 0.015081
# 3rd Qu.: 0.666917
# Max. : 2.874340

# Set up state formulas
psiformulas <- c("~elevation", "~1")

# Set up detection formulas
detformulas <- c("~1", "~1", "~date")

# Fit model
fm2 <- occuMS(detformulas, psiformulas, data = umf2, parameterization = "condbinom")
summary(fm2)

# Call:
# occuMS(detformulas = detformulas, psiformulas = psiformulas, data = umf2,
# parameterization = "condbinom")

# Occupancy (logit-scale):
# Estimate SE z P(>|z|)
# psi (Intercept) 1.940 0.473 4.106 4.02e-05
# psi elevation 0.160 0.387 0.413 6.79e-01
# R (Intercept) -0.301 0.243 -1.238 2.16e-01

# Detection (logit-scale):
# Estimate SE z P(>|z|)
# p[1] (Intercept) -0.0441 0.224 -0.197 8.44e-01
# p[2] (Intercept) 1.3222 0.263 5.018 5.23e-07
# delta (Intercept) 1.3145 0.301 4.373 1.22e-05
# delta date 0.2669 0.265 1.008 3.13e-01


# simplest dynamic multi-state occupancy model
# ''''''''''''''''''''''''''''''''''''''''''''
# Construct unmarked frame
ywide <- matrix(y, nrow = dim(y)[1], ncol = dim(y)[2]*dim(y)[3])
summary(umf3 <- unmarkedFrameOccuMS(y = ywide - 1, numPrimary = dim(y)[3]) )

# unmarkedFrame Object

# 100 sites
# Maximum number of observations per site: 15
# Mean number of observations per site: 15
# Number of primary survey periods: 5
# Number of secondary survey periods: 3
# Sites with at least one detection: 93

# Tabulation of y observations:
#   0   1   2
# 826 318 356

psiformulas <- c("~1", "~1")

# Detection formulas for p1, p2, and delta
detformulas <- c("~1", "~1", "~1")

# Formulas for phi transition matrix
umf3@phiOrder$cond_binom # Look at guide
# [1] "phi[0]" "phi[1]" "phi[2]" "R[0]" "R[1]" "R[2]"

# Intercept-only formulas for all
phiformulas <- rep("~1", 6)

# Fit model
fm3 <- occuMS(detformulas, psiformulas, phiformulas, data = umf3,
    parameterization = "condbinom")
summary(fm3)

# Call:
# occuMS(detformulas = detformulas, psiformulas = psiformulas,
# phiformulas = phiformulas, data = umf3, parameterization = "condbinom")

# Initial Occupancy (logit-scale):
# Estimate SE z P(>|z|)
# psi (Intercept) 1.66 0.306 5.42 5.91e-08
# R (Intercept) -0.24 0.234 -1.02 3.06e-01

# Transition Probabilities (logit-scale):
# Estimate SE z P(>|z|)
# phi[0] (Intercept) -1.730 0.285 -6.060 1.36e-09
# phi[1] (Intercept) 1.010 0.227 4.452 8.50e-06
# phi[2] (Intercept) 2.498 0.347 7.198 6.13e-13
# R[0] (Intercept) -0.322 0.578 -0.557 5.78e-01
# R[1] (Intercept) -0.939 0.252 -3.731 1.91e-04
# R[2] (Intercept) 1.735 0.283 6.132 8.65e-10

# Detection (logit-scale):
# Estimate SE z P(>|z|)
# p[1] (Intercept) 0.0481 0.115 0.418 6.76e-01
# p[2] (Intercept) 1.2374 0.112 11.076 1.63e-28
# delta (Intercept) 1.3424 0.127 10.538 5.76e-26

round( cbind(JAGS = unlist(out3$mean[1:4]), unmarked = plogis(coef(fm3))[1:8]), 3)
#       JAGS unmarked
# psi  0.833    0.840
# r    0.443    0.440
# phi1 0.158    0.151
# phi2 0.731    0.733
# phi3 0.918    0.924
# rho1 0.428    0.420
# rho2 0.287    0.281
# rho3 0.846    0.850

# dynamic multi-state occupancy model with some dummy covariates
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Simulate random site, yearly site and observational covariates
set.seed(39)
site_covs <- data.frame(elevation = rnorm(nsites))
ysite_covs <- data.frame(snow = rnorm(nsites * dim(y)[3]))
obs_covs <- data.frame(date = rnorm(length(ywide)))

# Create unmarked frame
summary(
umf4 <- unmarkedFrameOccuMS(ywide - 1, siteCovs = site_covs,
    yearlySiteCovs = ysite_covs, obsCovs = obs_covs, numPrimary = dim(y)[3]) )

# unmarkedFrame Object

# 100 sites
# Maximum number of observations per site: 15
# Mean number of observations per site: 15
# Number of primary survey periods: 5
# Number of secondary survey periods: 3
# Sites with at least one detection: 93

# Tabulation of y observations:
# 0 1 2
# 826 318 356

# Site-level covariates:
# elevation
# Min. :-2.32145
# 1st Qu.:-0.62779
# Median :-0.04193
# Mean : 0.05228
# 3rd Qu.: 0.86542
# Max. : 2.28150

# Observation-level covariates:
# date
# Min. :-3.19225
# 1st Qu.:-0.64407
# Median : 0.05337
# Mean : 0.02240
# 3rd Qu.: 0.68832
# Max. : 3.12377

# Yearly-site-level covariates:
# snow
# Min. :-2.55590
# 1st Qu.:-0.75018
# Median :-0.05207
# Mean :-0.06135
# 3rd Qu.: 0.61092
# Max. : 2.34549

# Set up state formulas
# Put covariate on occupancy for state 1, state 2 is intercept only
psiformulas <- c("~elevation", "~1")

# Set up detection formulas
# p1 and p2 are still intercept-only, delta gets date covariate
detformulas <- c("~1", "~1", "~date")

# Put "snow" covariate on transition prob. from state 0 to 1 (phi[0])
umf4@phiOrder$cond_binom # Remind ourselves of the order
phiformulas <- c("~snow", rep("~1", 5))

# Fit model
fm4 <- occuMS(detformulas, psiformulas, phiformulas, data = umf4,
    parameterization = "condbinom")
summary(fm4)

# Call:
# occuMS(detformulas = detformulas, psiformulas = psiformulas,
# phiformulas = phiformulas, data = umf4, parameterization = "condbinom")

# Initial Occupancy (logit-scale):
# Estimate SE z P(>|z|)
# psi (Intercept) 1.697 0.319 5.33 1.00e-07
# psi elevation -0.270 0.318 -0.85 3.95e-01
# R (Intercept) -0.247 0.234 -1.06 2.90e-01

# Transition Probabilities (logit-scale):
# Estimate SE z P(>|z|)
# phi[0] (Intercept) -1.750 0.293 -5.979 2.24e-09
# phi[0] snow -0.186 0.316 -0.587 5.57e-01
# phi[1] (Intercept) 1.021 0.228 4.472 7.74e-06
# phi[2] (Intercept) 2.493 0.347 7.182 6.89e-13
# R[0] (Intercept) -0.341 0.585 -0.583 5.60e-01
# R[1] (Intercept) -0.927 0.250 -3.711 2.07e-04
# R[2] (Intercept) 1.721 0.280 6.138 8.38e-10

# Detection (logit-scale):
# Estimate SE z P(>|z|)
# p[1] (Intercept) 0.0448 0.115 0.391 6.96e-01
# p[2] (Intercept) 1.2416 0.112 11.104 1.20e-28
# delta (Intercept) 1.3652 0.130 10.538 5.79e-26
# delta date -0.2006 0.124 -1.622 1.05e-01
