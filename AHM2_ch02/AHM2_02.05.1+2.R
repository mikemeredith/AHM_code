#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# run time for this script 3 mins

library(jagsUI)
library(AHMbook)

# 2.5 Dynamic N-mixture model of Dail-Madsen
# ==========================================

# 2.5.1 A Dail-Madsen simulator
# -----------------------------

simDM0 <- function(nsites = 50, nsurveys = 3, nyears = 5, lambda = 4,
    phi = 0.8, gamma = 1.5, p = 0.7){

  ## Simulation for multiple-visit data (from pcountOpen help file)
  ## No covariates, constant time intervals between primary periods
  # nsites: Number of sites
  # nsurveys: Number of rep. (secondary) samples within period of closure
  # nyears: Number of primary samples (= period of closure):
  # years, seasons etc.
  # lambda: Initial expected abundance
  # phi, gamma: apparent survival and recruitment rates, respectively
  # p: detection probability

  y <- array(NA, dim = c(nsites, nyears, nsurveys))
  N <- matrix(NA, nsites, nyears)
  S <- R <- matrix(NA, nsites, nyears-1)
  N[,1] <- rpois(nsites, lambda) # Initial state
  for(t in 1:(nyears-1)) { # State dynamics
    S[,t] <- rbinom(nsites, N[,t], phi) # Survival process
    R[,t] <- rpois(nsites, gamma) # Recruitment process
    N[,t+1] <- S[,t] + R[,t]
  }
  for(j in 1:nsurveys){ # Observation process
    y[,,j] <- rbinom(nsites*nyears, N, p)
  }

  # Put observed data into two dimensions
  yy <- array(NA, dim = c(nsites, nsurveys*nyears))
  for(t in 1:nyears){
    yy[,(nsurveys * t-(nsurveys-1)):(nsurveys*t)] <- y[,t,]
  }
  return(list(nsites = nsites, nsurveys = nsurveys, nyears = nyears,
      lambda = lambda, phi = phi, gamma = gamma, p = p, N = N, S = S, R = R,
      y = y, yy = yy))
}

# Execute function
set.seed(2017, kind = "L'Ecuyer")
str(data <- simDM0(nsites = 50, nsurveys = 3, nyears = 5, lambda = 4,
    phi = 0.8, gamma = 1.5, p = 0.7))


# 2.5.2 Fitting the DM model in BUGS
# ----------------------------------

# Bundle data set
str(bdata <- list(C = data$y, nsites = dim(data$y)[1], nsurveys = dim(data$y)[3],
    nyears = dim(data$y)[2]))
# List of 4
# $ C       : int [1:50, 1:5, 1:3] 3 4 4 1 2 5 3 5 7 5 ...
# $ nsites  : int 50
# $ nsurveys: int 3
# $ nyears  : int 5

# Specify model in BUGS language
cat(file = "DM1.txt","
model {
  # Priors
  lambda ~ dunif(0, 100)   # Initial site-specific abundance
  phi ~ dunif(0, 1)        # Apparent survival (omega in paper/unmarked)
  gamma ~ dunif(0, 5)      # Per-capita recruitment rate
  p ~ dunif(0, 1)          # Detection probability

  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t])     # Survival process
      # R[i,t+1] ~ dpois(gamma)        # 'absolute' recruitment = 'constant'
      R[i,t+1] ~ dpois(N[i,t] * gamma) # per-capita recruitment = 'autoreg'
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    }
    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,t,j] ~ dbin(p, N[i,t])
      }
    }
  }
}
")

# Initial values
Nst <- apply(data$y, c(1,2), max) + 2
Nst[, 2:5] <- NA                   # cols 2:5 of N are deterministic, N <- S + R.
R1 <- apply(data$y, c(1,2), max)   # Observed max. counts + 1 as inits
R1[,1] <- NA
inits <- function(){list( lambda = runif(1, 6, 16), phi = runif(1),
    gamma = runif(1), p = runif(1), N = Nst, R = R1 + 1 )}

# Parameters monitored
params <- c("lambda", "phi", "gamma", "p")

# MCMC settings
na <- 1000 ; ni <- 25000 ; nt <- 4 ; nb <- 5000 ; nc <- 3

# Call JAGS (ART 2 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "DM1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# par(mfrow = c(2, 3))  #  ~~~ replace with 'layout' argument
traceplot(out1, layout=c(2,3))
print(out1, 3)

# Per-capita recruitment parameterisation
#         mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# lambda 4.168 0.315 3.575 4.159 4.813    FALSE 1 1.000 15000
# phi    0.863 0.034 0.789 0.865 0.923    FALSE 1 1.002  1982
# gamma  0.232 0.036 0.168 0.230 0.311    FALSE 1 1.002  1682
# p      0.696 0.018 0.659 0.696 0.729    FALSE 1 1.001  5292

# To choose the absolute recruitment parameterization, we edit the BUGS code above to fit the
# absolute (or "constant") recruitment parameterization simply by commenting out the recruitment line
# for "per capita" and then uncommenting the line previous to that.
# ~~~~~~~~~~~~~~~~~~ here's the new JAGS code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat(file = "DM1b.txt","
model {
  # Priors
  lambda ~ dunif(0, 100) # Initial site-specific abundance
  phi ~ dunif(0, 1)      # Apparent survival (omega in paper/unmarked)
  gamma ~ dunif(0, 5)    # Per-capita recruitment rate
  p ~ dunif(0, 1)        # Detection probability

  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t])       # Survival process
      R[i,t+1] ~ dpois(gamma)            # 'absolute' recruitment = 'constant'
      # R[i,t+1] ~ dpois(N[i,t] * gamma) # per-capita recruitment = 'autoreg'
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    } # end t
    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,t,j] ~ dbin(p, N[i,t])
      } # end j
    } # end t
  } # end i
}
")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Call JAGS (ART 1.3 min), check convergence and summarize posteriors
out1b <- jags(bdata, inits, params, "DM1b.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  #  ~~~ replace with 'layout' argument
traceplot(out1b, layout=c(2,3))
print(out1b, 2)

# Absolute parameterisation
#        mean   sd 2.5%  50% 97.5% overlap0 f Rhat n.eff
# lambda 4.09 0.31 3.51 4.08  4.72    FALSE 1    1  8586
# phi    0.86 0.03 0.79 0.86  0.92    FALSE 1    1   748
# gamma  1.14 0.17 0.85 1.13  1.51    FALSE 1    1  1061
# p      0.70 0.02 0.66 0.70  0.73    FALSE 1    1  5387

# ~~~ save the work so far ~~~
save.image("AHM2_02.05.2.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~