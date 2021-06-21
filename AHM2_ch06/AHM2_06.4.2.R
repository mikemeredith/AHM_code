#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6 : MULTISTATE OCCUPANCY MODELS
# =======================================
# Code from proofs dated 2020-08-19

# Approximate time with full number of iterations: 14 mins

library(jagsUI)

# ~~~~~~~~~~~ need data preparation from 6.4.1 ~~~~~~~~~~~~~~~~~~
source("AHM2_06.4.1.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 6.4 Case study: Swiss eagle owls
# ================================

# 6.4.2 Fitting the simplest possible dynamic multiseason model
# -------------------------------------------------------------

# Summarize and bundle data
str(bdata <- list(y = yms, nsites = dim(yms)[1], nsurveys = nsurveys,
    nyears = dim(yms)[3]))
# List of 4
# $ y       : num [1:274, 1:20, 1:10] NA NA 1 NA NA NA NA 2 2 2 ...
# $ nsites  : int 274
# $ nsurveys: num [1:274, 1:10] 1 1 3 1 1 1 1 20 12 20 ...
# $ nyears  : int 10

# Specify model in BUGS language
cat(file = "dynMS1.txt", "
model {

  ### (1) Priors for parameters
  # State process priors
  # Priors for parameters in initial state vector (Omega)
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)

  # Priors for parameters in state transition matrix (PhiMat)
  for(s in 1:3){
    phi[s] ~ dunif(0, 1)
    rho[s] ~ dunif(0, 1)
  }

  # Priors for parameters in observation process (Theta)
  p2 ~ dunif(0, 1)                 # Detection prob. when in state 2
  for (s in 1:3) {                 # Detection prob. when in state 3
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- beta[s] / sum(beta[])
  }

  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector: Year 1
  Omega[1] <- 1 - psi              # Prob. of non-occupation
  Omega[2] <- psi * (1-r)          # Prob. of occ. by single bird
  Omega[3] <- psi * r              # Prob. of occ. by pair

  # Define state transition probability matrix (PhiMat): years 2:nyears
  # Define probabilities of state S(t+1) given S(t)
  # For now, constant over sites and years
  # Note conditional Bernoulli parameterization of multinomial
  # Order of indices: Departing state, arrival state
  PhiMat[1,1] <- 1 - phi[1]
  PhiMat[1,2] <- phi[1] * (1 - rho[1])
  PhiMat[1,3] <- phi[1] * rho[1]
  PhiMat[2,1] <- 1 - phi[2]
  PhiMat[2,2] <- phi[2] * (1 - rho[2])
  PhiMat[2,3] <- phi[2] * rho[2]
  PhiMat[3,1] <- 1 - phi[3]
  PhiMat[3,2] <- phi[3] * (1 - rho[3])
  PhiMat[3,3] <- phi[3] * rho[3]

  # Define observation probability matrix (Theta)
  # Order of indices: true state, observed state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-p2
  Theta[2,2] <- p2
  Theta[2,3] <- 0
  Theta[3,1] <- p3[1]
  Theta[3,2] <- p3[2]
  Theta[3,3] <- p3[3]

  ### (3) Likelihood
  # Initial state: year 1
  for (i in 1:nsites){
    z[i,1] ~ dcat(Omega[])
  }

  # State transitions from yearly interval 1:(nyears-1)
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      z[i,t+1] ~ dcat(PhiMat[z[i,t],])
    }
  }

  # Observation equation
  for (i in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys[i,t]){
        y[i,j,t] ~ dcat(Theta[z[i, t], ])
      }
    }
  }

  ### (4) Derived quantities
  # Number of sites in each state per year
  for (t in 1:nyears){
    for (i in 1:nsites){
      state1[i,t] <- equals(z[i,t], 1)   # Indicator for site in state 1
      state2[i,t] <- equals(z[i,t], 2)   # ... state 2
      state3[i,t] <- equals(z[i,t], 3)   # ... state 3
    }
    n.occ[t,1] <- sum(state1[,t])        # Number of unoccupied sites
    n.occ[t,2] <- sum(state2[,t])        # Number of sites with single birds
    n.occ[t,3] <- sum(state3[,t])        # Number of sites with pairs
    n.occ.total[t] <- n.occ[t,2] + n.occ[t, 3] # All occupied
  }
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(bdata$nsites, bdata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("psi", "r", "phi", "rho", "p2", "p3", "Omega", "PhiMat",
    "Theta", "n.occ", "n.occ.total") # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
odms1 <- jags(bdata, inits, params, "dynMS1.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(odms1)
print(odms1, 3)
#           mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# psi      0.699 0.061 0.577 0.698 0.818    FALSE 1 1.001  3000
# r        0.925 0.063 0.775 0.942 0.997    FALSE 1 1.001  1592
# phi[1]   0.420 0.060 0.310 0.418 0.543    FALSE 1 1.002  2475
# phi[2]   0.762 0.058 0.638 0.768 0.865    FALSE 1 1.001  1546
# phi[3]   0.992 0.006 0.977 0.993 1.000    FALSE 1 1.004   619
# rho[1]   0.087 0.061 0.005 0.077 0.231    FALSE 1 1.006   890
# rho[2]   0.189 0.049 0.105 0.185 0.297    FALSE 1 1.003  1810
# rho[3]   0.886 0.017 0.851 0.887 0.918    FALSE 1 1.002   778
# p2       0.290 0.022 0.246 0.291 0.333    FALSE 1 1.001  2098
# p3[1]    0.196 0.007 0.183 0.196 0.209    FALSE 1 1.001  1663
# p3[2]    0.411 0.008 0.394 0.411 0.426    FALSE 1 1.000  3000
# p3[3]    0.393 0.008 0.378 0.393 0.410    FALSE 1 1.001  1494
# Omega[1] 0.301 0.061 0.182 0.302 0.423    FALSE 1 1.001  3000
# Omega[2] 0.054 0.048 0.002 0.040 0.173    FALSE 1 1.001  1684
# Omega[3] 0.645 0.057 0.528 0.647 0.749    FALSE 1 1.000  3000
# [ ... output truncated ... ]

round(odms1$mean$PhiMat, 2) # Transition probabilities
#      [,1] [,2] [,3]
# [1,] 0.58 0.38 0.04
# [2,] 0.24 0.62 0.14
# [3,] 0.01 0.11 0.88

round(odms1$mean$Theta, 2) # Observation probabilities
#      [,1] [,2] [,3]
# [1,] 1.00 0.00 0.00
# [2,] 0.71 0.29 0.00
# [3,] 0.20 0.41 0.39
