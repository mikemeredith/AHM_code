#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# 6.3. Simulation and analysis of the simplest possible N-mixture model
# ------------------------------------------------------------------------

# Choose sample sizes and prepare observed data array C
set.seed(24)                # So we all get same data set
M <- 150                    # Number of sites
J <- 2                      # Number of abu. measurements per site (rep. counts)
C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
lambda <- 2.5               # Expected abundance
p <- 0.4                    # Probability of detection (per individual)

# Generate local abundance data (the truth)
N <- rpois(n = M, lambda = lambda)

# Conduct repeated measurements (generate replicated counts)
for(j in 1:J){
   C[,j] <- rbinom(n = M, size = N, prob = p)
}

# Look at data
# The truth ....
table(N)                    # True abundance distribution
sum(N)                      # True total population size at M sites
sum(N>0)                    # True number of occupied sites
mean(N)                     # True mean abundance (estimate of lambda)

# ... and the observations
table(apply(C, 1, max))     # Observed abundance distribution (max count)
sum(apply(C, 1, max))       # Observed total population size at M sites
sum(apply(C, 1, max)>0)     # Observed number of occupied sites
mean(apply(C, 1, max))      # Observed mean "relative abundance"

head(cbind(N=N, count1=C[,1], count2=C[,2])) # First 6 sites

cor(C)[1,2]

library(unmarked)                  # Load package
umf <- unmarkedFramePCount(y = C)  # Create um data frame
summary(umf)                       # Summarize
(fm1 <- pcount(~1 ~1, data = umf)) # Fit model: get estimates on link scale
backTransform(fm1, "state")        # Get estimates on natural scale
backTransform(fm1, "det")


# Bundle and summarize data set
win.data <- list(C = C, M = nrow(C), J = ncol(C))
str(win.data)                      # Look at data

# Specify model in BUGS language
sink("model1.txt")
cat("
model {
  # Priors
  lambda ~ dgamma(0.001, 0.001)
  p ~ dunif(0, 1)
  # Likelihood
  for (i in 1:M) {
    N[i] ~ dpois(lambda)      # State model
    for (j in 1:J) {
      C[i,j] ~ dbin(p, N[i]) # Observation model
    }
  }
}
",fill = TRUE)
sink()

# Initial values
Nst <- apply(C, 1, max)       # Avoid data/model/inits conflict
inits <- function(){list(N = Nst)}

# Parameters monitored
params <- c("lambda", "p")

# MCMC settings
ni <- 25000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

# Call JAGS (ART 1 min) and summarize posteriors
library(jagsUI)
fm2 <- jags(win.data, inits, params, "model1.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb)
print(fm2, dig = 3)

