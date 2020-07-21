#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

# 10.3 Simulation and analysis of the simplest possible site-occupancy model
# --------------------------------------------------------------------------


# Choose sample sizes and prepare observed data array y
set.seed(24)                  # So we all get same data set
M <- 100                      # Number of sites
J <- 2                        # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
psi <- 0.8                    # Probability of occupancy or presence
p <- 0.5                      # Probability of detection

# Generate presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi)  # R has no Bernoulli

# Generate detection/nondetection data (i.e. presence/absence measurements)
for(j in 1:J){
   y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}

# Look at data
sum(z)                        # True number of occupied sites
sum(apply(y, 1, max))         # Observed number of occupied sites
head(cbind(z=z, y))           # Truth and measurements for first 6 sites

library(unmarked)
umf <- unmarkedFrameOccu(y = y)  # Create unmarked data frame
summary(umf)                     # Summarize data frame
(fm1 <- occu(~1 ~1, data = umf)) # Fit model

backTransform(fm1, "state")      # Get estimates on probability scale
backTransform(fm1, "det")


# Bundle data and summarize data bundle
str( win.data <- list(y = y, M = nrow(y), J = ncol(y)) )

# Specify model in BUGS language
sink("model.txt")
cat("
model {
  # Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  # Likelihood
  for (i in 1:M) {    # Loop over sites
    z[i] ~ dbern(psi)         # State model
    for (j in 1:J) { # Loop over replicate surveys
      y[i,j] ~ dbern(z[i]*p)  # Observation model (only JAGS !)
      # y[i,j] ~ dbern(mu[i])  # For WinBUGS define 'straw man'
    }
    # mu[i] <- z[i]*p          # Only WinBUGS
  }
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, 1, max)       # Avoid data/model/inits conflict
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("psi", "p")

# MCMC settings
ni <- 5000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3

# Call JAGS and summarize posteriors
library(jagsUI)
fm2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb)
print(fm2, dig = 3)

