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

# 9.4 Mark-Recapture/Double observer Distance Sampling
# ====================================================


# 9.4.1 Simulating MRDS data
# ------------------------------------------------------------------------
# Simulate a double-observer sampling data set
set.seed(1235)
temp <- simHDStr(type="point", method="double") # simulate double observer point count data set
data <- temp$data         # harvest data
B <- temp$B               # upper limit of counting (maximum count distance)
nsites <-temp$nsites      # number of sites
habitat <-temp$habitat    # habitat covariate


# Processing of the data: pad the count vector with 0s etc.
data <- data[!is.na(data[,2]),]
n <- rep(0,nsites)
names(n) <- 1:nsites
n[names(table(data[,1]))] <- table(data[,1])
site <- data[,1]
dclass <- data[,"d"]      # categorical distance class for each observation
aux <- data[,"aux"]       # the auxiliary variable is capture history

# Create the categorical distance variable, use 10 classes here.
nD <- 10
delta <- B/nD # bin width
mdpts <-seq(delta/2,B,delta) # midpoint of bins up to max distance
nobs <- nrow(data)
dclass <- dclass%/%delta  +1

# Bundle data and look at overview of data
str( win.data <-list(n=n,site=site, dclass=as.numeric(dclass), nsites=nsites,
    nobs=nobs, delta=delta, nD=nD, mdpts=mdpts, B=B, aux=aux, habitat=habitat) )


# 9.4.2 Analysis in BUGS
# ------------------------------------------------------------------------
# Define model in BUGS langauge
cat("
model {

  #Priors for fixed detection parameters
  # 2 observer detection probability parameters
  logitp1 ~ dnorm(0, 0.01)
  logitp2 ~ dnorm(0, 0.01)
  # Intercepts
  alpha0 ~ dnorm(0, 0.01)    # intercept for sigma
  alpha1 ~ dnorm(0, 0.01)    # slope on sigma covariate
  # Coefficients
  beta0 ~ dnorm(0,0.01)      # intercept for lambda
  beta1 ~ dnorm(0,0.01)      # slope for lambda covariate

  # Detection scale parameter model
  for(s in 1:nsites){
    # Covariates on scale parameter (perceptibility)
    log(sigma[s]) <- alpha0 + alpha1*habitat[s]
    # Double observer cell probabilities here if there are covariates
    logit(pobs[1,s]) <- logitp1 # + covariates
    logit(pobs[2,s]) <- logitp2 # + covariates

    # Distance sampling model and cell probabilities
    for(b in 1:nD){
      log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  # half-normal
      f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) # Scaled radial density function
      pi.pd[b,s] <- g[b,s]*f[b,s]          #  Product Pr(detect)*Pr(distribution)
      pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]   # Conditional cell probabilities
    }
    pdet[s] <- sum(pi.pd[,s])              # Marginal probability of detection

    # Double observer cell probabilities and conditional probabilities
    doprobs[1,s] <- pobs[1,s]*(1-pobs[2,s])
    doprobs.condl[1,s] <- doprobs[1,s]/sum(doprobs[,s])
    doprobs[2,s] <- (1-pobs[1,s])*pobs[2,s]
    doprobs.condl[2,s] <- doprobs[2,s]/sum(doprobs[,s])
    doprobs[3,s] <- pobs[1,s]*pobs[2,s]
    doprobs.condl[3,s] <- doprobs[3,s]/sum(doprobs[,s])
    pavail[s] <- sum(doprobs[,s])  # probability of availability AT ALL
  }

  # Observation model for two categorical covariates
  for(i in 1:nobs){
    dclass[i] ~ dcat(pi.pd.c[,site[i]])
    aux[i] ~ dcat(doprobs.condl[,site[i]])
  }

  # Abundance model
  for(s in 1:nsites){
    # Binomial model for # of captured individuals
    n[s] ~ dbin(pdet[s], N[s])
    N[s] ~ dbin(pavail[s], M[s])   # binomial availability model
    # Abundance model
    M[s] ~ dpois(lambda[s])        # predicted abundance per survey/site/point
    # Add site-level covariates to lambda
    log(lambda[s])<- beta0 + beta1*habitat[s]
  }
  # Derived parameters
  Mtot <- sum(M[])
  Ntot <- sum(N[])
  logit(p1) <- logitp1
  logit(p2) <- logitp2
  sigma0 <- exp(alpha0)            # Baseline sigma
}
", fill=TRUE,file="do_model.txt")

# Inits function
Nst <- n + 1         # inits for N
inits <- function(){list(M=Nst+1, N=Nst, alpha0=runif(1,1,2),
   beta0=runif(1,-1,1), beta1=runif(1,-1,1), alpha1=runif(1,-1,1),
   logitp1=0, logitp2=0)}

# Parameters to monitor
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntot", "Mtot", "logitp1",
   "logitp2", "p1", "p2", "sigma0")

# MCMC settings
ni <- 50000   ;   nb <- 10000   ;   nt <- 4   ;   nc <- 3

# Run JAGS in parallel (ART 6.8 min), check convergence and summarize the results
out3 <- jags(data=win.data, inits=inits, parameters.to.save=params,
   model.file="do_model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
   parallel = TRUE)
traceplot(out3)   ;   print(out3, 3)


# Put true values into a vector
truth <- temp$parms
psi <- 1-(1-truth["p.double1"])*(1-truth["p.double2"])# Compute availability
truth <- c(truth[c("p.double1", "p.double2")], exp(truth["alpha0"]),
   truth["beta0"], "Mtot" = sum(temp$M),
   "Ntot" = sum(temp$M)*as.numeric(psi),
   truth[c("alpha0","alpha1","beta1")])

# Get posterior means and 2.5% and 97.5% percentiles (95% CRI)
post <- out3$summary[c("p1", "p2", "sigma0", "beta0", "Mtot", "Ntot",
   "alpha0", "alpha1", "beta1"), c(1,3,7)]

# Table compares truth with posterior mean and 95% CRI from JAGS
cbind(truth, posterior = round(post, 3))


# 9.4.3 Remarks (no code)
# ------------------------------------------------------------------------

