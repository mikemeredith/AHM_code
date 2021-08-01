#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9. Advanced Hierarchical Distance Sampling
# =========================================================================

# Approximate execution time for this code: 10 mins
# Run time with the full number of iterations: 1.6 hrs

library(AHMbook)
library(unmarked)
library(jagsUI)

# ~~~~~~~ changes to RNG defaults ~~~~~~~~~~~~~~~~~~~~~~~~
# Use the old default random number generator to get the printed numbers
RNGversion("3.2.0")

# 9.7 Open HDS models: modeling population dynamics
# ==================================================

# 9.7.1 Simulating the ISSJ data over multiple years
# ------------------------------------------------------------------------
# We load the ISSJ data analyzed in chapter 8, package into an unmarked frame
library(unmarked)
# library(rjags)
data(issj)
covs <- issj[,c("elevation","forest","chaparral")]
area <- pi*300^2 / 100^2             # Area in ha
jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
   siteCovs=data.frame(covs, area), numPrimary=1,
   dist.breaks=c(0, 100, 200, 300),
   unitsIn="m", survey="point")
sc <- siteCovs(jayumf)
sc.s <- scale(sc)
sc.s[,"area"] <- pi*300^2 / 10000  # Don't standardize area
covs<- siteCovs(jayumf) <- sc.s
summary(jayumf)

# Fit the model using gdistsamp and look at the fit summary
(nb.C2E.C <- gdistsamp( ~chaparral + I(chaparral^2) + elevation , ~1, ~chaparral,
   data =jayumf, output="abund", mixture="NB", K = 150))

# Get coefficient estimates to be used in data simulation
beta <- coef(nb.C2E.C)
betaFall <- beta[c("lambda(Int)", "lambda(chaparral)",
   "lambda(elevation)", "lambda(I(chaparral^2))")]

# Predict expected abundance per point count on log-scale for simulation
Xmat <- cbind(rep(1,307),covs[,3],covs[,3]^2,covs[,1]) # Order: chap, chap^2, elev
loglam <- Xmat%*%(betaFall)
lamnew <- exp(loglam)

# Parameters of the detection function
dparm <- beta[c("p(Int)", "p(chaparral)")]
sigma <- exp(Xmat[,c(1, 2)]%*%dparm)
J <- nsites <- 307 # number of sampling points

# Number of years
nyrs <- 6

# Set dynamics parameters to achieve a target growth rate of 0.95
phi <- 0.6       # Survival probability
gamma <- 0.35    # Recruitment rate

# Distance category info
db <- c(0,50, 100, 150, 200, 250, 300)
midpt <- c(25, 75, 125, 175, 225, 275)
nD <- length(midpt)
delta <- 50       # Distance interval width
B <- 300

# Simulate an ISSJ data set and harvest the data objects
set.seed(2015)
dat <- issj.sim(B=300, db = db, lam=lamnew, sigma=sigma, phi=phi, gamma=gamma,
    npoints=nsites, nyrs=nyrs)

y <- dat$y
dclass <- dat$dclass
site <- dat$site

# Bundle and summarize the data set
str(data1<-list(nsites=nsites, chap=as.vector(covs[,"chaparral"])[dat$cell],
   chap2=as.vector(covs[,"chaparral"]^2)[dat$cell],
   elev=as.vector(covs[,"elevation"])[dat$cell], T=nyrs, nD=nD, midpt=midpt,
   B=B, delta=delta, y=y, dclass=dclass, site=site, nind=sum(y)) )

# 9.7.2 Fitting a slurry of open population models
# ------------------------------------------------------------------------

# 9.7.2.1 The independence model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
model{

  # Prior distributions
  # Regression parameters
  alpha0 ~ dunif(0,20)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-20,20)
  beta1 ~ dunif(-20,20)
  beta2 ~ dunif(-20,20)
  beta3 ~ dunif(-20,20)
  beta4 ~ dunif(-20,20) # Population trend parameter
  r ~ dunif(0,5)        # NegBin dispersion parameter
  rout <- log(r)

  # 'Likelihood'
  for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]
    # Compute detection probability
    for(k in 1:nD){
      pi[k,s] <- (2*midpt[k]*delta )/(B*B)
      log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
      f[k,s] <- p[k,s]*pi[k,s]
      fc[k,s] <- f[k,s]/pcap[s]
      fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability

    # Process model
    for (t in 1:T){
      log(lambda[s,t]) <- beta0 + beta1*chap[s] + beta2*chap2[s] +
          beta3*elev[s] + beta4*(t - t/2)  # Note trend parameter here
      y[s,t] ~ dbin(pcap[s], N[s,t])
      # ~~~ Original Negative Binomial formulation ~~~~~~~~~~~~
      N[s,t] ~ dnegbin(prob[s,t], r)
      prob[s,t] <- r/(r+lambda[s,t])
      # ~~~ Alternative Poisson-gamma formulation ~~~~~~~~~~~~~
      # N[s,t] ~ dpois(lamx[s,t])
      # lamx[s,t] ~ dgamma(r, r/lambda[s,t])
      # This is mathematically equivalent but implemented differently in JAGS:
      #   no need to disable Conjugate samplers.
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    } # End loop over years
  } # End loop over sites

  # Distance sampling observation model for observed (binned) distance data
  for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]])
  }
  # Derived parameters
  for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)   # 300 m point = 28.27 ha
  }
}
", file="Sollmann1.txt")

# Set up initial values, parameters vector and MCMC settings
Nst <- y+1  # this is for trend model
inits <- function(){list(N=Nst, beta0=runif(1), beta1=runif(1), beta2=runif(1),
   beta3=runif(1), beta4=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), r = 1)}

params <-c('beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'alpha0', 'alpha1',
    'Ntot', 'D', 'r')

# ni <- 22000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3
ni <- 2200   ;   nb <- 200   ;   nt <- 1   ;   nc <- 3  # ~~~~~ for testing

## Execute JAGS, look at convergence and summarize the results

# ~~~~~~ jagsUI now has a factories argument ~~~~~~~~~~~~~~~~~~~~~
## JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
# set.factory("bugs::Conjugate", FALSE, type="sampler")  # ~~~ no longer works

# Setting factories with rjags::set.factory before calling
#   jagsUI::jags does not work on current versions.
# Use the 'factories' argument to jagsUI::jags instead.

# USE THIS with the original Negative Binomial formulation:
open1 <- jags (data1, inits, params, "Sollmann1.txt",
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
  factories="bugs::Conjugate sampler FALSE", parallel=TRUE)

# USE THIS with the Poisson-gamma formulation:
# open1 <- jags (data1, inits, params, "Sollmann1.txt", n.thin=nt, n.chains=nc,
   # n.burnin=nb, n.iter=ni, parallel=TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(open1)
print(open1, 2)


# 9.7.2.2 The reduced-dynamics model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
model{

  # Prior distributions
  # Regression parameters
  alpha0 ~ dunif(0,20)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-20,20)
  beta1 ~ dunif(-20,20)
  beta2 ~ dunif(-20,20)
  beta3 ~ dunif(-20,20)
  theta ~ dunif(0,5)
  # NegBin dispersion parameter
  r ~ dunif(0,5)
  rout <- log(r)

  # 'Likelihood'
  for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]
    # Compute detection probability
    for(k in 1:nD){
      log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
      f[k,s] <- p[k,s]*pi[k,s]
      fc[k,s] <- f[k,s]/pcap[s]
      fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
      pi[k,s] <- (2*midpt[k]*delta )/(B*B)
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability

    # Process model
    # Abundance model for Yr1 as in Sillett et al 2012
    log(lambda[s,1]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s]
    y[s,1] ~ dbin(pcap[s], N[s,1])
    N[s,1] ~ dnegbin(prob[s,1], r)
    prob[s,1] <- r/(r+lambda[s,1])

    # Population dynamics model for subsequent years
    for (t in 2:T){
      N[s,t] ~ dpois(N[s, t-1] * theta)
      y[s,t] ~ dbin(pcap[s], N[s,t])
    }
  }
  # Distance sampling observation model for observed (binned) distance data
  for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]])
  }

  # Derived parameters
  for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)  # 300 m point = 28.27 ha
  }
}
", file="Sollmann2.txt")

# Set up initial values, parameters vector and MCMC settings
Nst <- y+1 # this is for trend model
inits <- function(){list(N=Nst, beta0=runif(1), beta1=runif(1), beta2=runif(1),
   beta3=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), theta=runif(1,0.6,0.99))}

params <- c('beta0', 'beta1', 'beta2', 'beta3', 'alpha0', 'alpha1', 'theta',
   'rout', 'Ntot', 'D', 'r')

# ni <- 22000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3
ni <- 2200   ;   nb <- 200   ;   nt <- 1   ;   nc <- 3  # ~~~~ for testing

# Execute JAGS, look at convergence and summarize the results
open2 <- jags (data1, inits, params, "Sollmann2.txt",
  # n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel=TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(open2)
print(open2, 2)

# 9.7.2.3 The Glorious Integrated HDS/Dail-Madsen Model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
model{
  # Prior distributions
  # Regression parameters
  alpha0 ~ dunif(0,20)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-20,20)
  beta1 ~ dunif(-20,20)
  beta2 ~ dunif(-20,20)
  beta3 ~ dunif(-20,20)

  # Priors for dynamics parameters: here they are constant across years
  # We could add covariate models for logit(phi) and log(gamma)
  phi ~ dunif(0,1)
  gamma ~ dunif(0,5)

  # NegBin dispersion parameter
  r ~ dunif(0,5)
  rout <- log(r)

  # 'Likelihood'
  for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]

    # Compute detection probability
    for(k in 1:nD){
      log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
      f[k,s] <- p[k,s]*pi[k,s]
      fc[k,s] <- f[k,s]/pcap[s]
      fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
      pi[k,s] <- (2*midpt[k]*delta )/(B*B)
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability

    # Process model
    # Abundance model for year 1
    log(lambda[s,1]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s]
    y[s,1] ~ dbin(pcap[s], N[s,1])
    N[s,1] ~ dnegbin(prob[s,1], r)
    prob[s,1] <- r/(r+lambda[s,1])

    # Population dynamics model for subsequent years
    for (t in 2:T){                      # Loop over years
      S[s,t] ~ dbinom(phi, N[s, t-1])   # Survivors
      R[s,t] ~ dpois(gamma * N[s, t-1]) # Recruits
      N[s,t] <- S[s,t] + R[s,t]         # N = Survivors + Recruits
     y[s,t]~ dbin(pcap[s],N[s,t])       # Measurement error
    }
  }

  # Distance sampling observation model for observed (binned) distance data
  for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]])
  }

  # Derived parameters
  for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)     # 300 m point = 28.27 ha
  }
}
", file="Sollmann3.txt")


# Set up some sensible starting values for S and R
yin <- y+1
yin[,2:6] <- NA
Sin <- Rin <- matrix(NA, nrow=nsites, ncol=nyrs)
y1 <- y + 1
for(s in 1:nsites){
  for (t in 2:6){
    Sin[s,t] <- rbinom(1,y1[s,t-1], phi )
    Rin[s,t] <- ifelse((y1[s,t]-Sin[s,t])>0, y1[s,t]-Sin[s,t], 0)
  }
}

# Set up initial values, parameters vector and MCMC settings
inits <-function(){list(N=yin, beta0=runif(1), beta1=runif(1), beta2=runif(1),
   beta3=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), phi=0.6, gamma=0.3,
   R=Rin, S=Sin) }

params <- c('beta0', 'beta1', 'beta2', 'beta3', 'alpha0', 'alpha1', 'phi',
   'gamma', 'Ntot', 'D', 'r')

# ni <- 152000   ;   nb <- 2000   ;   nt <- 10   ;   nc <- 3
ni <- 15200   ;   nb <- 200   ;   nt <- 1   ;   nc <- 3  # ~~~~~ for testing

# Run JAGS, look at convergence and summarize the results
library(jagsUI)
set.seed(1) # ~~~ prevents "node incompatible..." error
open3  <- jags (data1, inits, params, "Sollmann3.txt", n.thin=nt,
    n.chains=nc, n.burnin=nb, n.iter=ni, parallel=TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(open3)
print(open3, 2)


# Compare inferences in graph .... (Fig. 9-6)
plot(apply(dat$N,2,sum),ylim=c(600,1300),xlab="Year",ylab="Population size (307 sample units)")
lines(apply(open1$sims.list$Ntot,2,mean), lty=1, col="blue", lwd=2)
lines(apply(open2$sims.list$Ntot,2,mean), lty=1, col="red", lwd=2)
lines(apply(open3$sims.list$Ntot,2,mean), lty=1, col="green", lwd=2)
open1.95cri <- apply(open1$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
open2.95cri <- apply(open2$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
open3.95cri <- apply(open3$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
legend(1,750,legend=c("Independence","Reduced dynamics","Full dynamics"), lty=1, col=c("blue","red","green"))

matlines(1:6, t(open1.95cri), type="l", lty=2, lwd=2, col="blue")
matlines(1:6, t(open2.95cri), type="l", lty=2, lwd=2, col="red")
matlines(1:6, t(open3.95cri), type="l", lty=2, lwd=2, col="green")

# .... and table
parms <- c(betaFall, dparm)
round(post <- cbind(parms, Independent=open1$summary[c(1:4,6,7),1],
       Partial=open2$summary[1:6,1], Full=open3$summary[1:6,1]), 3)

# 9.7.2.4 Summary remarks on modelling populations over time (no code)
