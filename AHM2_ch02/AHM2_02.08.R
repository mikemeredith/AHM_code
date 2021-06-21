#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

library(jagsUI)

# ~~~~ need this code from 2.7.2 ~~~~~~~~~~~
alfl <- read.csv(system.file("csv", "alfl.csv", package = "unmarked"))
alfl.covs <- read.csv(system.file("csv", "alflCovs.csv", package = "unmarked"),
    row.names = 1)
alfl$captureHistory <- paste(alfl$interval1, alfl$interval2, alfl$interval3, sep = "")
alfl$captureHistory <- factor(alfl$captureHistory,
    levels = c("001", "010", "011", "100", "101", "110", "111"))
alfl$id <- factor(alfl$id, levels = rownames(alfl.covs))
alfl.v1 <- alfl[alfl$survey == 1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
alfl.v2 <- alfl[alfl$survey == 2,]
alfl.H2 <- table(alfl.v2$id, alfl.v2$captureHistory)
alfl.v3 <- alfl[alfl$survey == 3,]
alfl.H3 <- table(alfl.v3$id, alfl.v3$captureHistory)
# Arrange the data in 3-d array format and also the wide format
Y <- array(NA, c(50, 3, 7))
Y[1:50,1,1:7] <- alfl.H1
Y[1:50,2,1:7] <- alfl.H2
Y[1:50,3,1:7] <- alfl.H3
Ywide <- cbind(alfl.H1, alfl.H2, alfl.H3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.8 Multinomial mixtures with full dynamics
# ===========================================

# Prepare time and date (incl. median-centering)
date <- as.matrix(alfl.covs[,c("date.1", "date.2", "date.3")])
date <- date - median(date)
time <- as.matrix(alfl.covs[,c("time.1", "time.2", "time.3")])
time <- time - median(time)

# Prepare response array
y3d <- array(NA, dim = c(nrow(Y), 7, 3)) # Create 3d array
y3d[,,1] <- Ywide[,1:7]
y3d[,,2] <- Ywide[,8:14]
y3d[,,3] <- Ywide[,15:21]

# Sample sizes
nseasons <- 3                   # Number of primary occasions
nsites <- nrow(Ywide)           # Number of sites
nobs <- apply(y3d, c(1,3), sum) # Total detections/site, occasion

# Bundle data (sames as before)
str(bdata <- list(y3d = y3d, nsites = nsites, nseasons = nseasons,
    nobs = nobs, woody = alfl.covs[,"woody"], struct = alfl.covs[,"struct"],
    date = date, time = time))

# Dail-Madsen-type, open model for multinomial sampling protocol
# Specify model in BUGS language
cat(file = "mnDM.txt", "
model{

  # Prior distributions, regression parameters
  alpha0 ~ dnorm(0,0.01)
  alpha1 ~ dnorm(0,0.01)
  beta0 ~ dunif(-20,20)
  beta1 ~ dunif(-20,20)
  beta2 ~ dunif(-20,20)
  # Priors for dynamics parameters: here they are constant across years
  # Here, we could add covariate models for logit(phi) and log(gamma)
  phi ~ dunif(0,1)
  gamma ~ dunif(0,5)

  # 'Likelihood'
  for (i in 1:nsites){
    for (t in 1:nseasons){
      # Linear model for detection function scale
      log(p[i,t]) <- alpha0 + alpha1*woody[i]
      # Multinomial cell probabilities
      cp[i,t,1] <- (1-p[i,t])*(1-p[i,t])*p[i,t]
      cp[i,t,2] <- (1-p[i,t])*p[i,t]*(1-p[i,t])
      cp[i,t,3] <- (1-p[i,t])*p[i,t]*p[i,t]
      cp[i,t,4] <- p[i,t]*(1-p[i,t])*(1-p[i,t])
      cp[i,t,5] <- p[i,t]*(1-p[i,t])*p[i,t]
      cp[i,t,6] <- p[i,t]*p[i,t]*(1-p[i,t])
      cp[i,t,7] <- p[i,t]*p[i,t]*p[i,t]
      cp[i,t,8] <- 1-sum(cp[i,t,1:7])
      cellprobs.cond[i,t,1] <- cp[i,t,1]/ sum(cp[i,t,1:7])
      cellprobs.cond[i,t,2] <- cp[i,t,2]/ sum(cp[i,t,1:7])
      cellprobs.cond[i,t,3] <- cp[i,t,3]/ sum(cp[i,t,1:7])
      cellprobs.cond[i,t,4] <- cp[i,t,4]/ sum(cp[i,t,1:7])
      cellprobs.cond[i,t,5] <- cp[i,t,5]/ sum(cp[i,t,1:7])
      cellprobs.cond[i,t,6] <- cp[i,t,6]/ sum(cp[i,t,1:7])
      cellprobs.cond[i,t,7] <- cp[i,t,7]/ sum(cp[i,t,1:7])

      # Conditional 4-part version of the model
      pdet[i,t] <- sum(cp[i,t, 1:7])
      y3d[i,1:7,t] ~ dmulti(cellprobs.cond[i,t,1:7], nobs[i,t])
      # Conditional observation model
      nobs[i,t] ~ dbin(pdet[i,t], N[i,t]) # Individuals detected
    }
     # Poisson regression model for abundance
    log(lambda0[i]) <- beta0 + beta1*woody[i] + beta2*struct[i]
    N[i,1] ~ dpois(lambda0[i]) # Population size

    # Population dynamics model for subsequent years
    for (t in 2:nseasons){ # Loop over primary periods
      S[i,t] ~ dbinom(phi, N[i, t-1]) # Survivors
      R[i,t] ~ dpois(gamma * N[i, t-1]) # Recruits
      N[i,t] <- S[i,t] + R[i,t] # N = Survivors + Recruits
      y[i,t] ~ dbin(pdet[i,t],N[i,t]) # Observation model
    }
  }
  # Derived parameters
  for(t in 1:nseasons){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (0.785*nsites) # 50 m point = 0.785 ha
  }
}
")

# Set up some sensible starting values for S and R
nseasons <- 3
yin <- nobs+1
yin[,2:3] <- NA
Sin <- Rin <- matrix(NA, nrow = nsites, ncol = nseasons)
y1 <- nobs + 1
for(i in 1:nsites){
  for (t in 2:3){
    Sin[i,t] <- rbinom(1,y1[i,t-1], 0.7)
    Rin[i,t] <- ifelse((y1[i,t]-Sin[i,t])>0, y1[i,t]-Sin[i,t], 0)
  }
}
# Initial values
inits <- function(){list(N = yin, beta0 = runif(1), beta1 = runif(1),
    beta2 = runif(1), beta3 = runif(1), alpha0 = runif(1, -3, -2),
    alpha1 = runif(1), phi = 0.6, gamma = 0.3, R = Rin, S = Sin) }

# Parameters monitored
params <- c('beta0', 'beta1', 'beta2', 'beta3', 'alpha0', 'alpha1',
    'phi', 'gamma', 'Ntot', 'D')

# MCMC settings
na <- 5000 ; ni <- 100000 ; nb <- 50000 ; nt <- 50 ; nc <- 3

# Run JAGS (ART 6 min), look at convergence and summarize posteriors
out8 <- jags (bdata, inits, params, "mnDM.txt", n.adapt = na, n.iter=ni,
    n.burnin=nb, n.thin=nt, n.chains=nc, parallel=TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out8)
print(out8, 3)
#           mean    sd   2.5%    50%  97.5% overlap0     f  Rhat n.eff
# beta0   -1.122 0.415 -1.929 -1.112 -0.314    FALSE 0.998 1.001  2834
# beta1    2.083 0.692  0.744  2.074  3.473    FALSE 0.999 1.000  3000
# beta2    0.046 0.037 -0.026  0.047  0.118     TRUE 0.892 1.000  3000
# alpha0  -0.546 0.116 -0.778 -0.541 -0.333    FALSE 1.000 1.000  3000
# alpha1   0.334 0.235 -0.132  0.335  0.778     TRUE 0.922 1.000  3000
# phi      0.392 0.105  0.139  0.407  0.556    FALSE 1.000 1.006   430
# gamma    0.199 0.106  0.059  0.176  0.481    FALSE 1.000 1.003   918
# Ntot[1] 55.781 1.466 54.000 56.000 59.000    FALSE 1.000 1.000  3000
# Ntot[2] 33.927 1.010 33.000 34.000 36.000    FALSE 1.000 1.000  3000
# Ntot[3] 17.953 1.038 17.000 18.000 20.000    FALSE 1.000 1.000  3000
# D[1]     1.421 0.037  1.376  1.427  1.503    FALSE 1.000 1.000  3000
# D[2]     0.864 0.026  0.841  0.866  0.917    FALSE 1.000 1.000  3000
# D[3]     0.457 0.026  0.433  0.459  0.510    FALSE 1.000 1.000  3000
