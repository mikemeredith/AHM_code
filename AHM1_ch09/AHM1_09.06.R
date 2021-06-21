#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9. Advanced Hierarchical Distance Sampling
# =========================================================================

library(AHMbook)
library(jagsUI)

# 9.6 Open HDS models: Implicit Dynamics
# ======================================


# Obtain a data set
set.seed(1236)
str(tmp <- simHDSopen("point", nreps=7, nyears=5, nsites=100, beta.trend=0.2) )
attach(tmp)

apply(tmp$M.true,2,sum)  # True population size per year

# Define distance class information
delta <- 0.5
nD <- B%/%delta                 # Number of distance classes
midpt <- seq(delta/2, B, delta) # Mid-point of distance intervals

# Create the 4-d array
y4d <- array(0, dim=c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- tmp$data[[yr]][[rep]]
    site <- data[,1]
    dclass <- data[,"d"]%/%delta + 1
    ndclass <- B%/%delta
    dclass <- factor(dclass, levels=  1:ndclass)
# ~~~~~ this cannot work ~~~~~~~~~~
    # y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass)
# ~~~~~ use this instead ~~~~~~~~~~~
    ttt <- table(site, dclass)
    siteID <- as.numeric(rownames(ttt))
    y4d[siteID,1:nD,rep,yr] <- ttt
  }
}


# Bundle and summarize the data set
nobs <- apply(y4d, c(1,3,4), sum)  # Total detections per site and occasion
str( data <- list(y4d=y4d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta,
    habitat=habitat, B=B, nobs = nobs, T=tmp$nyears) )

# Define model in BUGS
cat("
model {

  # Prior distributions
  beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
  mean.lam <- exp(beta0)
  beta1 ~ dnorm(0, 0.01)  # Coefficient on habitat
  phi ~ dunif(0,1)        # Probability of availability
  sigma ~ dunif(0,5)      # Detection function parameter
  beta.trend ~ dnorm(0, 0.01)

  # Construct the multinomial cell probabilities
  for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal
    f[b] <- (2*midpt[b]*delta)/(B*B)                # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
  }
  cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])
  for (s in 1:nsites) {
    for (k in 1:K) {
      pdet[s,k] <- sum(cellprobs[1:nD]) # Distance class probabilities
      pmarg[s,k] <- pdet[s,k]*phi       # Marginal probability
    }
  }

  for(t in 1:T){                        # Years
    for (s in 1:nsites) {               # Sites
      for (k in 1:K) {                  # Replicates
        # Model part 4: distance class frequencies
        y4d[s,1:nD,k,t] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k,t])
        # Model part 3: total number of detections:
        nobs[s,k,t] ~ dbin(pmarg[s,k], M[s,t])
        # Model part 2: Availability. Not used in this model but simulated.
        Navail[s,k,t] ~ dbin(phi, M[s,t])
      }  # end k loop
      # Model part 1: Abundance model
      M[s,t] ~ dpois(lambda[s,t])
      log(lambda[s,t]) <- beta0 + beta1*habitat[s] + beta.trend*(t-2.5)
    }  # end s loop
  } # end t loop

  # Derived quantities
  for(t in 1:T){
    Mtot[t] <- sum(M[,t])
      for(k in 1:K){
        Ntot[k,t] <- sum(Navail[,k,t])
    }
  }
} # End model
",file="tempemig4d.txt")

# Inits and parameters to save
Navail.st <- apply(y4d, c(1,3,4),sum)
Mst <- apply(Navail.st, c( 1,3), max) +2
inits <- function(){
  list(M=Mst, Navail = Navail.st, sigma = 1.0, phi=.9,beta0=log(2),beta1=.5)
}
params <- c("sigma", "phi", "beta0", "mean.lam", "beta.trend",
   "beta1", "Mtot", "Ntot")

# MCMC settings
# ni <- 12000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3
ni <- 1200   ;   nb <- 200   ;   nt <- 1   ;   nc <- 3  # ~~~~ use for testing

# Run JAGS (ART 9 min), look at trace plots and summarize
# ~~~~~~~ this crashes reliably ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   with error "Error in node cellprobs.cond[6] Invalid parent values"
#   when script is run as-is, with the RNG seed specified on line 16.
# Set a new seed to overcome this:
set.seed(1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
outRD <- jags(data, inits, params, "tempemig4d.txt",
  # n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = FALSE)
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)  # ~~~~ for faster testing
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(outRD)
summary(outRD)

