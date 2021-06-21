#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 15 mins

library(jagsUI)

# 2.10 Spatially dynamic Dail-Madsen models
# =========================================

# Set parameter values and sample sizes
nsites <- 25   # number of sites (must be integer^2)
nyears <- 20   # number of years
nsurveys <- 3  # number of surveys per year
lambda0 <- 5   # initial abundance
phi <- 0.8     # survival rate
gamma <- 0.2   # reproduction rate
kappa <- 0.3   # emigration rate
p <- 0.75      # detection probability

# Define a grid of sites and adjacency structure
lat <- rep(1:sqrt(nsites), times = sqrt(nsites))
lon <- rep(1:sqrt(nsites), each = sqrt(nsites))
plot(lat, lon, pch = 20, las = 1)
dist <- matrix(, nsites, nsites)
for (i in 1:nsites) {
  for (j in 1:nsites) {
    dist[i,j] <- sqrt((lat[i] - lat[j])^2 + (lon[i] - lon[j])^2)
  }
}
adj <- ifelse(dist <= 1.05, 1, 0) # Binary matrix that identifies neighbours
diag(adj) <- 0
nadj <- colSums(adj) # Number of neighbours

# Simulate a data set
N <- matrix(NA, nsites, nyears)
R <- S <- E <- I <- Ilambda <- matrix(NA, nsites, nyears - 1)
y <- array(NA, c(nsites, nyears, nsurveys))
N[,1] <- rpois(nsites, lambda0)
for (t in 2:nyears) {
  R[,t-1] <- rpois(nsites, gamma*N[,t-1])
  S[,t-1] <- rbinom(nsites, N[,t-1], phi)
  E[,t-1] <- rbinom(nsites, S[,t-1], kappa)
  for (i in 1:nsites) {
    Ilambda[i,t-1] <- sum(E[,t-1] / nadj * adj[,i])
  }
  I[,t-1] <- rpois(nsites, Ilambda[,t-1])
  N[,t] <- S[,t-1] - E[,t-1] + R[,t-1] + I[,t-1]
}
for (j in 1:nsurveys) {
  y[,,j] <- rbinom(nsites*nyears, N, p)
}

# Bundle data
str( bdata <- list(nsites = nsites, nyears = nyears, nsurveys = nsurveys,
    y = y, adj = adj, nadj = nadj))
# List of 6
# $ nsites   : num 25
# $ nyears   : num 20
# $ nsurveys : num 3
# $ y        : int [1:25, 1:20, 1:3] 3 3 3 4 5 4 2 1 4 2 ...
# $ adj      : num [1:25, 1:25] 0 1 0 0 0 1 0 0 0 0 ...
# $ nadj     : num [1:25] 2 3 3 3 2 3 4 4 4 3 ...

# Specify model in BUGS language
cat(file="spatialDMmodel.txt", "
model {
  # Prior distributions
  lambda0 ~ dgamma(0.001, 0.001)
  phi ~ dunif(0, 1)
  gamma ~ dunif(0, 1)
  kappa ~ dunif(0, 1)
  p ~ dunif(0, 1)

  # Process model
  for(i in 1:nsites) {
    N[i,1] ~ dpois(lambda0)
    for(t in 2:nyears) {
      R[i,t-1] ~ dpois(gamma*N[i,t-1])
      S[i,t-1] ~ dbin(phi, N[i,t-1])
      E[i,t-1] ~ dbin(kappa, S[i,t-1])
      Ilambda[i,t-1] <- sum(E[1:nsites,t-1]/ nadj[1:nsites] * adj[1:nsites,i])
      I[i,t-1] ~ dpois(Ilambda[i,t-1])
      N[i,t] <- S[i,t-1] - E[i,t-1] + R[i,t-1] + I[i,t-1]
    }
  }

  # Observation model
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      for (j in 1:nsurveys) {
        y[i,t,j] ~ dbin(p, N[i,t])
      }
    }
  }
}
")

# Initial values. cheap.inits = true values
Ni <- N + 20 ; Ni[,-1] <- NA ; Ri <- R + 10 ; Si <- S + 10
Ei <- E + 5 ; Ii <- I + 5
inits <- function() list(lambda0 = runif(1, 1, 5), phi = runif(1),
    gamma = runif(1), kappa = runif(1), p = runif(1), N = Ni, R = Ri, S = Si,
    E = Ei, I = Ii)

# Parameters monitored
params <- c("lambda0", "phi", "gamma", "kappa", "p")

# MCMC settings
# na <- 2000 ; ni <- 20000 ; nt <- 10 ; nb <- 10000 ; nc <- 3
na <- 2000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ testing

# Call JAGS (ART 22 min), check convergence and summarize posteriors
out11 <- jags(bdata, inits, params, model = "spatialDMmodel.txt",
    n.adapt = na, n.chains = nc, n.burnin = nb, n.iter = ni, n.thin = nt,
    parallel = TRUE)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out11, layout=c(2,2))
print(out11)
#          mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# lambda0 4.738 0.463 3.878 4.716 5.671    FALSE 1 1.000  3000
# phi     0.825 0.035 0.752 0.826 0.888    FALSE 1 1.010   236
# gamma   0.170 0.035 0.109 0.169 0.240    FALSE 1 1.009   234
# kappa   0.272 0.034 0.211 0.272 0.341    FALSE 1 1.006   533
# p       0.748 0.011 0.726 0.748 0.768    FALSE 1 1.000  3000
