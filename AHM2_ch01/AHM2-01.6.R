#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================

# library(AHMbook)

# 1.6 Gaussian State-Space Models for Inference About Relative Abundance
# ======================================================================

# 1.6.1 Modeling Multiple Time-Series With Fixed Effects for Site-Level Parameters
# --------------------------------------------------------------------------------

# Check number of zero years in all sites
nzero <- apply(C, 1, function(x) sum(x == 0, na.rm = TRUE))
plot(sort(nzero)) # Make a graph of number of zero years
sum(nzero <= 1) # 97 sites with at most 1 zero year
table(nzero)
# Bundle data with restriction on sites
sel <- nzero <= 1 # Select sites with <= 1 zero count
newM <- sum(sel) # Define new number of sites
str(bdata <- list(C = C[sel,], M = newM, T = ncol(C[sel,])))
# List of 3
# $ C: int [1:97, 1:18] 3 5 9 13 12 2 3 3 6 10 ...
# $ M: int 97
# $ T: int 18

# Specify model in BUGS language
cat(file = "model6.txt","
model {
  # Priors
  for(i in 1:M){
    n[i, 1] ~ dnorm(0, 0.01)I(0,) # Prior for initial pop. sizes
    # curve(dnorm(x, 0, sqrt(1/ 0.01)), 0, 50) # how does it look like ?
    mean.gamma[i] ~ dunif(0, 10) # Prior for mean growth rates
    sigma.proc[i] ~ dnorm(0, 1)I(0,) # Prior for sd of state process
    sigma2.proc[i] <- pow(sigma.proc[i], 2)
    tau.proc[i] <- pow(sigma.proc[i], -2)
    sigma.obs[i] ~ dnorm(0, 0.01)I(0,) # Prior for sd of obs. process
    sigma2.obs[i] <- pow(sigma.obs[i], 2)
    tau.obs[i] <- pow(sigma.obs[i], -2)
  }
  # ’Likelihood’
  # State process
  for (i in 1:M){
    for (t in 1:(T-1)){
      gamma[i, t] ~ dnorm(mean.gamma[i], tau.proc[i])
      n[i, t+1] <- n[i, t] * gamma[i, t]
    }
  }
  # Observation process
  for (i in 1:M){
    for (t in 1:T){
      C[i, t] ~ dnorm(n[i, t], tau.obs[i])
    }
  }
  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(n[,t])
  }
}
")

# Initial values
inits <- function(){list(sigma.proc = runif(newM, 0, 5),
  mean.gamma = runif(newM, 0.1, 2), sigma.obs = runif(newM, 0, 10),
  n = cbind(runif(newM, 0, 50), array(NA, dim = c(newM, ncol(C)-1))))}
# Parameters monitored
params <- c("mean.gamma", "sigma2.proc", "sigma2.obs", "popindex", "n")
# MCMC settings
# na <- 10000 ; ni <- 6e6 ; nt <- 1000 ; nb <- 5e6 ; nc <- 2
na <- 10000 ; ni <- 6e5 ; nt <- 100 ; nb <- 5e5 ; nc <- 2  # ~~~~~~ for testing
# Call JAGS (ART 505 min), check convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "model6.txt", n.adapt = na, n.chains = nc,
  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(4,4)) ; traceplot(out6) # all params
summary(out6) ; View(out6) ; print(out6$summary[1:320,-c(4:6)], 3)
# Check how many and which parameters have failed to converge
which(out6$summary[,8] > 1.1) # 7 derived quants or latent variables

# Produce Fig. 1.6
par(mfrow = c(1, 2))
graphSSM(out6, bdata$C)

# 1.6.2 MODELING MULTIPLE TIME-SERIES WITH RANDOM EFFECTS FOR SITELEVEL PARAMETERS

# no code

# 1.6.3 BRIEF COMMENTS ON GAUSSIAN STATE-SPACE MODELS

# no code
