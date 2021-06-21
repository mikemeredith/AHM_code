#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 30 mins
# Run time with the full number of iterations: 6 hrs

# library(AHMbook)
library(jagsUI)

# ~~~~~ Need to run 1.3 before this ~~~~~~~
source("AHM2_01.03.R")
# ~~~~~ and this from 1.4 ~~~~~~~~~~~~~~~~~
M <- nrow(C)
T <- ncol(C)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1.5 Generalized Linear Mixed Models
# ===================================

# 1.5.4 The “TRIM Model” with temporal autocorrelation
# ----------------------------------------------------

# Bundle data
str(bdata <- list(C = C, M = M, T = T))

# Specify model in BUGS language
cat(file = "model5.txt","
model {

  # Priors
  for(i in 1:M){
    site[i] ~ dnorm(0, 0.001) # Prior for site effects
  }
  year[1] <- 0                # Constraint on year effects
  for(t in 2:T){
    year[t] ~ dnorm(0, 0.001) # Prior for year effects 2:T
  }
  tau <- pow(sd, -2)
  sd ~ dunif(0, 3)
  rho ~ dunif(-1,1)           # Autoregressive param. for temp. autocorrelation

  # 'Likelihood'
  # First year
  for (i in 1:M){
    eps[i,1] ~ dnorm(0, tau) # unstructured random variation (= OD)
    C[i,1] ~ dpois(lambda[i,1])
    log(lambda[i,1]) <- site[i] + year[1] + w[i,1]
    w[i,1] <- eps[i,1] / sqrt(1 - rho * rho)

    # Later years
    for(t in 2:T){
      eps[i,t] ~ dnorm(0, tau) # (same) unstructured random variation
      C[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- site[i] + year[t] + w[i,t]
      w[i,t] <- rho * w[i,t-1] + eps[i,t]
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(lambda[,t])
  }
}
")

# Initial values
inits <- function() list(site = rnorm(M), year = c(NA, rnorm(T-1)),
    rho = runif(1), eps = array(0.1, dim=c(M, T)))

# Parameters monitored
params <- c("popindex", "site", "year", "lam.sel", "rho", "sd")

# MCMC settings
# na <- 10000 ; ni <- 250000 ; nt <- 200 ; nb <- 50000 ; nc <- 3
na <- 10000 ; ni <- 25000 ; nt <- 20 ; nb <- 5000 ; nc <- 3  # ~~~ for testing, 30 mins

# Call JAGS (ART 251 min), check convergence and summarize posteriors
out5 <- jags(bdata, inits, params, "model5.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,2))  # ~~~ replaced with 'layout' argument
traceplot(out5, layout=c(3,2))
summary(out5) ; jags.View(out5) ; print(out5, 3)
#      mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# rho 0.945 0.026 0.892 0.947 0.987    FALSE 1 1.032    71
# sd  0.211 0.012 0.188 0.210 0.235    FALSE 1 1.007   277
# [ ... ]

# ~~~~ Save output for use in subsequent sections ~~~~~
save(out5, file="AHM2_01.05.4_out5.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
