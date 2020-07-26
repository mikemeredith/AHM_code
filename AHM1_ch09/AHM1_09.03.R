#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9. Advanced Hierarchical Distance Sampling
# =========================================================================

# Approximate execution time for this code: 5 mins

library(AHMbook)
library(jagsUI)

# ~~~~~~~ changes to RNG defaults ~~~~~~~~~~~~~~~~~~~~~~~~
# Use the old default random number generator to get the printed numbers
RNGversion("3.2.0")

# 9.3 Time-removal and distance sampling combined
# ===============================================


# 9.3.1 The four-part hierarchical model
# ------------------------------------------------------------------------


# 9.3.2 Simulating some time-removal/DS data
# ------------------------------------------------------------------------
# Obtain a data set and harvest the results
set.seed(1235)                 # so we all create the same data set
temp <- simHDStr(type="point") # Simulate point count-removal data set
data <- temp$data              # harvest data
B <- temp$B                    # upper limit of counting (maximum count distance)
nsites <- temp$nsites          # Number of sites
habitat <- temp$habitat        # habitat covariate
K <- temp$K                    # Number of removal periods


# Create the observed encounter frequencies per site (include the zeros! )
data <- data[!is.na(data[,2]),]   # Sites where detections did occur
n <- rep(0,nsites)                # The full site vector
names(n) <- 1:nsites
n[names(table(data[,1]))] <- table(data[,1])  # Put in the counts
site <- data[,1]
nobs <- nrow(data)

# Create the distance class data
nD <- 10             # Number of distance classes
delta <- B/nD        # bin size or width
mdpts <- seq(delta/2,B,delta) # midpoint distance of bins up to max distance
dclass <- data[,"d"] # distance class for each observation
dclass <- dclass%/%delta  +1
tint <- data[,"aux"]

# Bundle data and summarize
str( win.data<-list(n=n, site=site, dclass=as.numeric(dclass),nsites=nsites,
    nobs=nobs, delta=delta, nD=nD,mdpts=mdpts,B=B, K=K, tint=tint,
    habitat=habitat) )


cat("
model {
  # Prior distributions for basic parameters
  # Intercepts
  beta.a0 ~ dnorm(0,0.01)    # intercept for availability
  alpha0 ~ dnorm(0, 0.01)    # intercept for sigma
  alpha1 ~ dnorm(0,0.01)     # slope on sigma covariate
  # Coefficients
  # beta.a1 ~ dnorm(0,0.01)  # slope for availability covariate
  beta0 ~ dnorm(0,0.01)      # intercept for lambda
  beta1 ~dnorm(0,0.01)       # slope for lambda covariate

  for(s in 1:nsites){
    # Add covariates to scale parameter DISTANCE (perceptibility)
    log(sigma[s]) <- alpha0 +  alpha1*habitat[s]
    # Add covariates for availability here TIME-REMOVAL (availability)
    p.a[s] <- exp(beta.a0) / (1+exp(beta.a0))
    # Optional covariates on availability
    # exp(beta.a0 + beta.a1*date[s])/(1+exp(beta.a0+beta.a1*date[s]))
    # Distance sampling detection probability model
    for(b in 1:nD){
      log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  # Half-normal
      f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) # Radial density function
      pi.pd[b,s] <- g[b,s]*f[b,s]  #  Product Pr(detect)*Pr(distribution)
      pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]  # Conditional probabilities
    }
    pdet[s] <- sum(pi.pd[,s])  # Probability of detection at all

    # Time-removal probabilities
    for (k in 1:K){
      pi.pa[k,s] <- p.a[s] * pow(1-p.a[s], (k-1))
      pi.pa.c[k,s] <- pi.pa[k,s]/phi[s] # Conditional probabilities of availability
    }
    phi[s] <- sum(pi.pa[,s]) # Probability of ever available
  }
  # Conditional observation model for categorical covariates
  for(i in 1:nobs){
    dclass[i] ~ dcat(pi.pd.c[,site[i]])
    tint[i] ~ dcat(pi.pa.c[,site[i]])
  }
  # Abundance model
  for(s in 1:nsites){
    # Binomial model for # of captured individuals
    # n[s] ~ dbin(pmarg[s], M[s]) # Formulation b, see text
    # pmarg[s] <- pdet[s]*phi[s]
    n[s] ~ dbin(pdet[s], N[s])    # Formulation a, see text
    N[s] ~ dbin(phi[s],M[s])      # Number of available individuals
    M[s] ~ dpois(lambda[s])       # Abundance per survey/site/point
    # Add site-level covariates to lambda
    log(lambda[s]) <- beta0 + beta1*habitat[s]
  }
  # Derived quantities
  Mtot <- sum(M[])  # Total population size
  Ntot <- sum(N[])  # Total available population size
  PDETmean <- mean(pdet[]) # Mean perceptibility across sites
  PHImean <- mean(phi[]) # Mean availability across sites
}
", fill=TRUE, file="tr-ds.txt")

# Create initial values (including for M and N) and list parameters to save
Mst <- Nst <- n + 1
inits <- function(){list(M=Mst, N=Nst, alpha0=1, beta0=runif(1,-1,1),
   beta.a1=runif(1,-1,1), beta1=runif(1,-1,1), alpha1=runif(1,-1,1),
   beta.a0=runif(1,-1,1))}
params <- c("beta.a0", "beta.a1", "alpha0", "alpha1", "beta0", "beta1",
    "PDETmean", "PHImean", "Mtot", "Ntot")

# MCMC settings
ni <- 50000   ;   nb <- 10000   ;   nt <- 4   ;   nc <- 3

# Run JAGS in parallel (ART 7.3 min), check convergence and summarize posteriors
out2a <- jags(data=win.data, inits=inits, parameters=params,
   model.file ="tr-ds.txt",n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
   parallel = TRUE)
traceplot(out2a)   ;   print(out2a, 3)

sum(temp$M)

# print(out2b,3)  # ~~~ should surely be...
print(out2a,3)

