#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8. Modeling abundance using hierarchical distance sampling (HDS)
# =========================================================================

# Approximate execution time for this code: 45 mins

library(AHMbook)
library(R2WinBUGS)
bd <- "C:/WinBUGS14" # Never forget this for WinBUGS

# 8.5 Bayesian HDS
# ================

# 8.5.1 Simulating some HDS data
# ------------------------------------------------------------------------
set.seed(1234)
tmp1 <- simHDS("point")   # Point transect
tmp2 <- simHDS()          # Line transect (this is the default)
str(tmp1)                 # Look at function output

# 8.5.2 Bayesian HDS using data augmentation
# ------------------------------------------------------------------------
# Recreate line transect data set
set.seed(1234)
tmp <- simHDS()                  # Line transect (default)
attach(tmp)

# Data augmentation: add a bunch of "pseudo-individuals"
nz <- 500                        # Augment by 500
nind <- nrow(data)
y <- c(data[,2], rep(0, nz))     # Augmented detection indicator y
site <- c(data[,1], rep(NA, nz)) # Augmented site indicator,
                                 # unknown (i.e., NA) for augmented inds.
d <- c(data[,5], rep(NA,nz))     # Augmented distance data (with NAs)

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, habitat=habitat, wind=wind, B=B,
    nind=nind, nz=nz, y=y, d=d, site=site) )
win.data$site                    # unknown site cov. for augmented inds.


# BUGS model for line transect HDS (NOT point transects!)
cat("
model{
  # Prior distributions
  beta0 ~ dunif(-10,10)   # Intercept of lambda-habitat regression
  beta1 ~ dunif(-10,10)   # Slope of log(lambda) on habitat
  alpha0 ~ dunif(-10,10)  # Intercept of log(sigma) (half-normal scale)
  alpha1 ~ dunif(-10,10)  # Slope of log(sigma) on wind

  # psi is a derived parameter under DA for stratified populations
  psi <- sum(lambda[]) / (nind+nz)

  # 'Likelihood' (sort of...)
  for(i in 1:(nind+nz)){                 # i is index for individuals
    z[i] ~ dbern(psi)                    # Data augmentation variables
    d[i] ~ dunif(0, B)                   # distance uniformly distributed
    p[i] <- exp(-d[i]*d[i]/(2*sigma[site[i]]*sigma[site[i]])) # Det. function
    mu[i] <- z[i]* p[i]                  # 'straw man' for WinBUGS
    y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
  }

  # Linear models for abundance and for detection
  for(s in 1:nsites){                    # s is index for sites
    # Model for abundance
    # next line not necessary, but allows to make predictions
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta1*habitat[s] # Linear model abundance
    site.probs[s] <- lambda[s] / sum(lambda[])

    # Linear model for detection
     log(sigma[s]) <- alpha0 + alpha1*wind[s]
  }
  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[])
  area<- nsites*1*2*B   # Unit length == 1, half-width = B
  D<- Ntotal/area
}
",fill=TRUE , file = "model1.txt")


# Inits
zst <- c(rep(1, sum(y)), rep(0, nz)) # ... and for DA variables
inits <- function(){list(beta0=0, beta1=0, alpha0=0, alpha1=0, z=zst)}

# Parameters to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Call BUGS (ART 33 min) ...
# bd <- "c:/Program Files/WinBUGS14/" # Never forget this for WinBUGS
out1 <- bugs(win.data, inits, params, "model1.txt",
   # n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, debug=TRUE, bugs.dir = bd)
   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, debug=FALSE, bugs.dir = bd) #~~~~~ for automated testing

# ... or try JAGS for a change (ART 6 min)
library(jagsUI)       # never forget to load jagsUI
out1 <- jags(win.data, inits, params, "model1.txt",
  # n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel=TRUE)  # ~~~~ speeds up testing

# Summarize posterior output
print(out1, 2)
sum(tmp$N.true)

# Prepare data
delta <- 0.1                    # width of distance bins for approx.
midpt <- seq(delta/2, B, delta) # make mid-points and chop up data
dclass <- d %/% delta + 1       # convert distances to cat. distances
nD <- length(midpt)             # Number of distance intervals

# Bundle and summarize data set
str( win.data <- list (y=y, dclass=dclass, site=site, midpt=midpt, delta=delta,
    B=B, nind=nind, nz=nz, nsites=nsites, nD=nD, habitat=habitat, wind=wind) )


# BUGS model specification for line-transect HDS (NOT point transects!)
cat("
model{
  # Prior distributions
  alpha0 ~ dunif(-10,10)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-10,10)
  beta1 ~ dunif(-10,10)

  psi <- sum(lambda[])/(nind+nz)     # psi is a derived parameter

  for(i in 1:(nind+nz)){             # Loop over individuals
    z[i] ~ dbern(psi)               # DA variables
    dclass[i] ~ dcat(pi[site[i],])  # Population distribution of dist class
    mu[i] <- z[i] * p[site[i],dclass[i]] # p depends on site AND dist class
    y[i] ~ dbern(mu[i])             # Basic Bernoulli response in DS model
    site[i] ~ dcat(site.probs[1:nsites]) # Site membership of inds
  }

  for(s in 1:nsites){                # Loop over sites
    # Construct cell probabilities for nG cells
    for(g in 1:nD){                    # midpt = mid point of each cell
      log(p[s,g]) <- -midpt[g]*midpt[g]/(2*sigma[s]*sigma[s])
      pi[s,g] <- delta/B              # probability of x per interval
      f[s,g] <- p[s,g]*pi[s,g]        # pdf of observed distances
    }

    # not necessary   N[s]~dpois(lambda[s]) except for prediction
    N[s] ~ dpois(lambda[s])        # predict abundance at each site
    log(lambda[s]) <- beta0 + beta1 * habitat[s] # linear model for N
    site.probs[s] <- lambda[s]/sum(lambda[])
    log(sigma[s]) <- alpha0 + alpha1*wind[s] # linear model for sigma
  }

  # Derived parameter
  Ntotal <- sum(z[])   # Also sum(N[]) which is size of a new population
  area <- nsites*1*2*B  # Unit length == 1, half-width = B
  D <- Ntotal/area
}
",fill=TRUE, file = "model2.txt")

# Inits
zst <- c(rep(1, sum(y)), rep(0, nz))
inits <- function(){list (alpha0=0, alpha1=0, beta0=0, beta1=0, z=zst) }

# Params to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal","D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Run JAGS with parallel processing (ART 1 min)
library(jagsUI)
out2 <- jags(win.data, inits, params, "model2.txt",
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
print(out2,2)


# 8.5.3 Bayesian HDS using the 3-part conditional multinomial model
# ------------------------------------------------------------------------
set.seed(1234)
tmp <- simHDS(type="line", discard0=FALSE)
attach(tmp)

# Get number of individuals detected per site
# ncap = 1 plus number of detected individuals per site
ncap <- table(data[,1])            # ncap = 1 if no individuals captured
sites0 <- data[is.na(data[,2]),][,1] # sites where nothing detected
ncap[as.character(sites0)] <- 0    # Fill in 0 for sites with no detections
ncap <- as.vector(ncap)

# Prepare other data
site <- data[!is.na(data[,2]),1]   # site ID of each observation
delta <- 0.1                       # distance bin width for rect. approx.
midpt <- seq(delta/2, B, delta)    # make mid-points and chop up data
dclass <- data[,5] %/% delta + 1   # convert distances to cat. distances
nD <- length(midpt)                # Number of distance intervals
dclass <- dclass[!is.na(data[,2])] # Observed categorical observations
nind <- length(dclass)             # Total number of individuals detected

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, nind=nind, B=B, nD=nD, midpt=midpt,
    delta=delta, ncap=ncap, habitat=habitat, wind=wind, dclass=dclass,
    site=site) )

# BUGS model specification for line-transect HDS (NOT point transects!)
cat("
model{
# Priors
  alpha0 ~ dunif(-10,10)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-10,10)
  beta1 ~ dunif(-10,10)

  for(i in 1:nind){
    dclass[i] ~ dcat(fc[site[i],]) # Part 1 of HM
  }

  for(s in 1:nsites){
    # Construct cell probabilities for nD multinomial cells
    for(g in 1:nD){                 # midpt = mid-point of each cell
      log(p[s,g]) <- -midpt[g] * midpt[g] / (2*sigma[s]*sigma[s])
      pi[s,g] <- delta / B          # probability per interval
      f[s,g] <- p[s,g] * pi[s,g]
      fc[s,g] <- f[s,g] / pcap[s]
    }
    pcap[s] <- sum(f[s,])           # Pr(capture): sum of rectangular areas

    ncap[s] ~ dbin(pcap[s], N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])         # Part 3 of HM
    log(lambda[s]) <- beta0 + beta1 * habitat[s] # linear model abundance
    log(sigma[s])<- alpha0 + alpha1*wind[s]      # linear model detection
  }
  # Derived parameters
  Ntotal <- sum(N[])
  area<- nsites*1*2*B  # Unit length == 1, half-width = B
  D<- Ntotal/area
}
",fill=TRUE, file = "model3.txt")

# Inits
Nst <- ncap + 1
inits <- function(){list(alpha0=0, alpha1=0, beta0=0, beta1=0, N=Nst)}

# Params to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntotal","D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# Run JAGS (ART 1 min) and summarize posteriors
library(jagsUI)
out3 <- jags(win.data, inits, params, "model3.txt", n.thin=nt,
   # n.chains=nc, n.burnin=nb, n.iter=ni)
   n.chains=nc, n.burnin=nb, n.iter=ni, parallel=TRUE)  # ~~~~ faster testing
print(out3, 2)


# 8.5.4 Point transect HDS using the conditional multinomial formulation
# ------------------------------------------------------------------------
# Simulate a data set using our simHDS function
set.seed(1234)
tmp <- simHDS(type="point", discard0=FALSE)
attach(tmp)

# Prepare data
# Number of individuals detected per site
ncap <- table(data[,1])            # ncap = 1 if no individuals captured
sites0 <- data[is.na(data[,2]),][,1] # sites where nothing was seen
ncap[as.character(sites0)] <- 0    # Fill in 0 for sites with no detections
ncap <- as.vector(ncap)            # Number of individuals detected per site

# Other data
site <- data[!is.na(data[,2]),1]   # Site ID of each observation
delta <- 0.1                       # Distance bin width for rect. approx.
midpt <- seq(delta/2, B, delta)    # Make mid-points and chop up data
dclass <- data[,5] %/% delta + 1   # Convert distance to distance category
nD <- length(midpt)                # Number of distance intervals
dclass <- dclass[!is.na(data[,2])] # Observed categorical observations
nind <- length(dclass)             # Total number of individuals detected

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, nind=nind, B=B, nD=nD, midpt=midpt,
    delta=delta, ncap=ncap, habitat=habitat, wind=wind, dclass=dclass,
    site=site) )

# BUGS model specification for point transect data
cat("
model{
  # Priors
  alpha0 ~ dunif(-10,10)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-10,10)
  beta1 ~ dunif(-10,10)

  for(i in 1:nind){
    dclass[i] ~ dcat(fc[site[i],]) # Part 1 of HM
  }
  for(s in 1:nsites){
    # Construct cell probabilities for nD distance bands
    for(g in 1:nD){                # midpt = mid-point of each band
      log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s] * sigma[s])
      pi[s,g] <- ((2 * midpt[g] ) / (B * B)) * delta # prob. per interval
      f[s,g] <- p[s,g] * pi[s,g]
      fc[s,g] <- f[s,g] / pcap[s]
    }
    pcap[s] <- sum(f[s,])           # Pr(capture): sum of rectangular areas

    ncap[s] ~ dbin(pcap[s], N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])         # Part 3 of HM
    log(lambda[s]) <- beta0 + beta1 * habitat[s] # linear model abundance
    log(sigma[s]) <- alpha0 + alpha1*wind[s]     # linear model detection
  }

  # Derived parameters
  Ntotal <- sum(N[])
  area <- nsites*3.141*B*B
  D <- Ntotal/area
}
",fill=TRUE, file="model4.txt")


# Inits
Nst <- ncap + 1
inits <- function(){list(alpha0=0, alpha1=0, beta0=0, beta1=0, N=Nst)}

# Params to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntotal","D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# Run BUGS (not STAN !) (ART 2.3 min) and summarize posteriors
out4 <- bugs(win.data, inits, params, "model4.txt",
   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
   # debug=TRUE, bugs.dir = bd)
   debug=FALSE, bugs.dir = bd) #~~~~~~~~~ for automated testing
print(out4, 2)

sum(tmp$N.true)                  # True realized population size
sum(!is.na(tmp$data[,"y"]))      # Observed population size


# 8.5.5 Analysis of the ISSJ data
# ------------------------------------------------------------------------
# Load the ISSJ data
library(unmarked)
data(issj)

# Prepare some data
nD <- 3                          # Number of intervals
delta <- 100                     # Interval width
B <- 300                         # Upper bound (max. distance)
midpt <- c(50, 150, 250)         # mid points

# Convert vector frequencies to individual distance class
H <- as.matrix(issj[,1:3])
nsites <- nrow(H)
ncap <- apply(H, 1, sum)         # Number of individuals detected per site
dclass <- rep(col(H), H)         # Distance class of each individual
nind <- length(dclass)           # Number of individuals detected
elevation <- as.vector(scale(issj[,c("elevation")])) # Prepare covariates
forest <- as.vector(scale(issj[,"forest"]))
chaparral <- as.vector(scale(issj[,"chaparral"]))


# Bundle and summarize data set
str( win.data <- list(nsites=nsites, nind=nind, B=B, nD=nD, midpt=midpt,
    delta=delta, ncap=ncap, chaparral=chaparral, elevation=elevation,
    dclass=dclass) )

# BUGS model specification
cat("
model{
  # Priors
  sigma ~ dunif(0,1000)
  beta0 ~ dunif(-10,10)
  beta1 ~ dunif(-10,10)
  beta2 ~ dunif(-10,10)
  beta3 ~ dunif(-10,10)
  sigma.site ~ dunif(0,10)
  tau <- 1/(sigma.site*sigma.site)

  # Specify hierarchical model
  for(i in 1:nind){
    dclass[i] ~ dcat(fc[]) # Part 1 of HM
  }

  # construct cell probabilities for nG cells
  for(g in 1:nD){                # midpt = mid-point of each cell
    log(p[g]) <- -midpt[g] * midpt[g] / (2 * sigma * sigma)
    pi[g] <- ((2 * midpt[g]) / (B * B)) * delta # prob. per interval
    f[g] <- p[g] * pi[g]
    fc[g] <- f[g] / pcap
  }
  pcap <- sum(f[])               # Pr(capture): sum of rectangular areas

  for(s in 1:nsites){
    ncap[s] ~ dbin(pcap, N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])      # Part 3 of HM
    log(lambda[s]) <- beta0 + beta1*elevation[s] + beta2*chaparral[s] + beta3*chaparral[s]*chaparral[s] + site.eff[s] # linear model for abundance
    site.eff[s] ~ dnorm(0, tau)   # Site log-normal 'residuals'
  }
  # Derived params
  Ntotal <- sum(N[])
  area <- nsites*3.141*300*300/10000   # Total area sampled, ha
  D <- Ntotal/area
}
",fill=TRUE, file="model5.txt")


# Inits
Nst <- ncap + 1
inits <- function(){list (sigma = runif(1, 30, 100), beta0 = 0, beta1 = 0,
    beta2 = 0, beta3 = 0, N = Nst, sigma.site = 0.2)}

# Params to save
params <- c("sigma", "beta0", "beta1", "beta2", "beta3", "sigma.site",
    "Ntotal","D")

# MCMC settings
# ni <- 52000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3
ni <- 7000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3  # ~~~~ for testing

# Run BUGS (ART 0.9 min) and summarize posteriors
out5 <- bugs(win.data, inits, params, "model5.txt",
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
  # debug=TRUE, bugs.dir = bd)
  debug=FALSE, bugs.dir = bd) #~~~~~ for automated testing

# Run JAGS (ART 0.5 min) and summarize posteriors
out5 <- jags(win.data, inits, params, "model5.txt",
  # n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel=TRUE)  # ~~~~ faster testing

print(out5, 3)

