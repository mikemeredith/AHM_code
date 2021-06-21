#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7 : MODELING FALSE POSITIVES
# ====================================
# Code from proofs dated 2020-08-19
library(unmarked)
library(jagsUI)

# 7.5 Bayesian analysis of models with false positives in JAGS
# =============================================================

# Random seed and simulation settings
set.seed(129, kind = "Mersenne")
nsites <- 200     # number of sites (i = 1, ..., M)
nsurveys <- 7     # number of visits (k = 1, ..., J)
psi <- 0.6        # expected psi
p <- 0.7          # detection probability (p_11)
fp <- 0.05        # false positive error probability (p_10)

# Simulate the latent states and the data
z <- matrix(NA, nrow = nsites, ncol = 1)        # empty matrix for occ states
z[1:nsites] <- rbinom(nsites, 1, psi)           # occupancy states
y <- matrix(NA, nrow = nsites, ncol = nsurveys) # empty matrix for det.
for(i in 1:nsites){
  pr_yequals1 <- p*z[i] + fp*(1-z[i])           # p11 + p10
  y[i,] <- rbinom(nsurveys, 1, pr_yequals1)     # realized observations
}

# Bundle data and summarize data bundle
str( bdata <- list(y = y, nsites = nrow(y), nsurveys = ncol(y)) )

# Specify model in BUGS language
cat(file = "occufp.txt","
model {

  # Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  fp ~ dunif(0, 1)

  # Likelihood and process model
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # State model
    for (j in 1:nsurveys) { # Loop over replicate surveys
      y[i,j] ~ dbern(z[i]*p + (1-z[i])*fp) # Observation model
    }
  }
}
")

# Initial values
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, p = 0.7, fp = 0.05)}

# Parameters monitored
params <- c("psi", "p", "fp")

# MCMC settings
na <- 1000 ; ni <- 5000 ; nt <- 1 ; nb <- 1000 ; nc <- 3

# Call JAGS (ART <1 min), assess convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "occufp.txt", n.adapt = na,
    n.chains = nc,n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 2))  # ~~~ replaced with 'layout' argument
traceplot(out1, layout=c(2,2))
print(out1, 3)
#      mean    sd  2.5%   50% 97.5% overlap0 f Rhat n.eff
# psi 0.614 0.036 0.543 0.614 0.684    FALSE 1    1 12000
# p   0.718 0.017 0.684 0.718 0.750    FALSE 1    1  5606
# fp  0.057 0.012 0.035 0.057 0.082    FALSE 1    1  9407

# With both Type 1 and Type 3 data, i.e., what they call
#   the "multiple detection states design"
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Random seed and simulation settings
set.seed(129)       # RNG seed
nsites <- 200       # number of sites
nsurv1 <- 3         # number of occasions with Type 1 data
nsurv2 <- 4         # number of occasions with Type 3 data
psi <- 0.6          # expected proportion of are occupied
p <- c(0.7, 0.5)    # detection prob of method 1 and method 2
fp <- 0.05          # false-positive error probability (p_10)
b <- 0.2            # probability y is recorded as certain

# Simulate the latent states and the data
z <- rbinom(nsites, 1, psi)
y <- matrix(NA, nrow = nsites, ncol = nsurv1 + nsurv2)
for(i in 1:nsites){
  p1 <- p[1]*z[i]                   # certainly detection (method 1)
  p2 <- p[2]*z[i] + fp*(1-z[i])     # uncertainly detection (method 2)
  y[i,1:3] <- rbinom(nsurv1, 1, p1) # simulate method 1 data
  y[i,4:7] <- rbinom(nsurv2, 1, p2) # simulate method 2 data
  # now introduce certain observations:
  pr.certain <- z[i] * y[i,4:7] * b
  y[i, 4:7] <- y[i, 4:7] + rbinom(4, 1, pr.certain)
}
head(y)                             # Look at data to understand them !

# Define a covariate for method
Method <- matrix(c(rep("1", 3), rep("2", 4)), nrow = nsites,
    ncol = nsurv1 + nsurv2, byrow = TRUE)
head(Method)

# Make data categorical so non-detection = 1, certain detection = 3
y[, nsurv2:(nsurv1 + nsurv2)] <- 1 + y[,nsurv2:(nsurv1 + nsurv2)]

# Bundle data and summarize data bundle (not shown)
str(bdata <- list(y = y, nsites = nrow(y), nsurv1 = nsurv1, nsurv2 = nsurv2))

# Specify model in BUGS language
cat(file = "occufp2.txt","
model {

  # Priors
  psi ~ dunif(0, 1)
  fp ~ dunif(0, 1)
  b ~ dunif(0, 1)
  alpha0 ~ dnorm(0,0.01)
  alpha1 ~ dnorm(0,0.01) # Method effect

  # Likelihood and process model
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # State model
    # Define observation matrix (obsmat)
    for(j in 1:(nsurv1+nsurv2)) {
      obsmat[i,j,1,1] <- 1-fp # z = 0 obs probs
      obsmat[i,j,2,1] <- fp
      obsmat[i,j,3,1] <- 0
      obsmat[i,j,1,2] <- 1-p[i,j] # z = 1 obs probs
      obsmat[i,j,2,2] <- (1-b)*p[i,j]
      obsmat[i,j,3,2] <- p[i,j]*b
    }
    # Observation model: part 1 (for first 3 cols in y)
    for(j in 1:nsurv1) { # Loop over replicate surveys
      logit(p[i,j]) <- alpha0
      y[i,j] ~ dbern(z[i]*p[i,j] ) # ordinary occupancy data
    }
    # Observation model: part 2 (for last 4 cols in y)
    for (j in (nsurv1+1):(nsurv1+nsurv2)) { # Loop over replicates
      logit(p[i,j]) <- alpha0 + alpha1
      y[i,j] ~ dcat(obsmat[i,j,1:3,z[i]+1] )
    }
  }
}
")

# Initial values
zst <- apply(y[, 1:nsurv1], 1, max)
inits <- function(){list(z = zst, alpha0 = rnorm(1, -1, 1),
    alpha1 = rnorm(1, 0, 1), fp = 0.05, b = 0.1)}

# Parameters monitored
params <- c("psi", "alpha0", "alpha1", "fp", "b")

# MCMC settings
na <- 1000 ; ni <- 5000 ; nt <- 1 ; nb <- 1000 ; nc <- 3

# Call JAGS (ART 1 min), assess convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "occufp2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out2, layout=c(2,3))
print(out2, 3)
#          mean    sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
# psi     0.633 0.034  0.566  0.634  0.699    FALSE 1 1.000  5589
# alpha0  1.121 0.125  0.879  1.120  1.371    FALSE 1 1.001  3616
# alpha1 -1.095 0.154 -1.397 -1.094 -0.794    FALSE 1 1.000  5442
# fp      0.040 0.012  0.019  0.039  0.065    FALSE 1 1.000  8842
# b       0.189 0.025  0.143  0.189  0.238    FALSE 1 1.001  4655
