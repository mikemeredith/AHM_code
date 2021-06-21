#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

library(jagsUI)

# ~~~~ need the Green Woodpecker data prepared in 2.2 ~~~~~~~~
source("AHM2_02.02.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.3 Year-stratified N-mixture model
# ===========================================

# Bundle and summarize data set
str(bdata <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3]))
# List of 4
# $ C       : int [1:267, 1:3, 1:14] 0 3 0 0 0 0 0 0 0 0 ...
# $ nsites  : int 267
# $ nsurveys: int 3
# $ nyears  : int 14

# Specify model in BUGS language
cat(file = "Nmix1.txt","
model {

  # Priors
  for (t in 1:nyears){           # Loop over years ('seasons')
    lambda[t] ~ dunif(0, 100)    # Expected abundance
    p[t] ~ dunif(0, 1)           # Detection probability
  }
  # Ecological model for true abundance
  for (i in 1:nsites){           # Loop over 26 sites
    for(t in 1:nyears){          # Loop over 14 years
      N[i,t] ~ dpois(lambda[t])
      # Observation model for replicated counts
      for (j in 1:nsurveys){     # Loop over 3 occasions
        C[i,j,t] ~ dbin(p[t], N[i,t])
      }
    }
  }
  # Derived quantity: Total abundance across all surveyed sites
  for (t in 1:nyears){
    totalN[t] <- sum(N[,t])      # includes sites with missing surveys
  }
}
")

# Initial values: avoid data/prior/inits conflict
Nst <- apply(C, c(1,3), max, na.rm = TRUE)+1
Nst[Nst == '-Inf'] <- 1
inits <- function() list(N = Nst, lambda = runif(dim(C)[3]))

# Parameters monitored
params <- c("lambda", "p", "totalN")

# MCMC settings
na <- 100 ; ni <- 3000 ; nt <- 2 ; nb <- 1000 ; nc <- 3

# Run JAGS (ART 1 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "Nmix1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out1)
print(out1, digits = 2) # not shown

# ~~~~ extra code for figure 2.2 ~~~~~~~~~~~~~~~~
# Summarize counts and abundance estimates (with 95% CRI)
op <- par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(2004:2017, annual.mean, xlab = 'Year', type = 'b', ylab = 'Mean count or E(N)',
  frame = FALSE, pch = 16, cex = 2, col = 'red', ylim = c(0, 2))
points(2004:2017, out1$mean$lambda, type = 'b', pch = 16, cex = 2)
segments(2004:2017, out1$q2.5$lambda, 2004:2017, out1$q97.5$lambda)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
