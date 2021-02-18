# Applied hierarchical modeling in ecology - vol.2 - 2021
# Marc Kéry & J. Andy Royle
#
# Chapter 11 : SPATIALLY EXPLICIT DISTANCE SAMPLING ALONG TRANSECTS
# =================================================================
# Code from proofs dated 2020-07-30

print("Approximate execution time for this code: 20 mins")
# Run time with the full number of iterations: 70 mins

library(AHMbook)
library(jagsUI)

# 11.3 Distance sampling models on linear transects
# =================================================

library(AHMbook)
RNGversion("3.5.0")
set.seed(1234, kind = "Mersenne-Twister")
# Simulation settings
N <- 200
M <- 500                # Data augmentation
sigma <- 0.20
alpha0 <- -2
W <- 1/2                # Transect dimension: half-width
L <- 4                  # ... length

# Locations of individuals
u1 <- runif(N, 0, 4)
u2 <- runif(N, 0, 2*W)
plot(u1, u2, pch = 20, col = 'black', cex = 1)  # Start building up Fig. 11.2
abline(0.5, 0, lwd = 2, col = 'grey')
title("Transect population subject to detection")

# Represent transect by a sequence of points with delta = 0.02
line.pts <- seq(0.01, L - 0.01, 0.02)

# Call these points "traps"
traplocs <- cbind(line.pts, 0.5)

# Compute some things:
# d.to.trap = distance from point on line to individual detected
# pbar = prob. an individual is detected at all during the survey
# obs.pos = observation point (on the line)
d.to.trap <- obs.pos <- pbar <- rep(NA, N)
dmat <- e2dist(cbind(u1, u2), traplocs)
for(i in 1:nrow(dmat)){
  haz <- exp(alpha0)*exp( -(dmat[i,]^2)/(2*sigma*sigma))
  probs <- 1 - exp(-haz)
  captured <- rbinom(nrow(traplocs), 1, probs)      # Bernoulli trials
  pbar[i] <- 1 - exp(-sum(haz))
  if(sum(captured)==0)
    next
  obs.pos[i] <- which(captured == 1)[1]
  d.to.trap[i] <- dmat[i,][obs.pos[i]]
  lines(c(u1[i], traplocs[obs.pos[i],1]), c(u2[i],
      traplocs[obs.pos[i],2]) )                     # Finish Fig. 11.2
}

# Subset to encountered individuals only
data <- cbind(obs.pos, u1, u2)[!is.na(obs.pos),]
obs.pos <- data[,1]
nind <- nrow(data)

# Distribution of individual detection probability (Fig. 11.3)
hist(pbar, nclass = 20, col = 'grey')

# Next we create the observation matrix Ymat by filling it with
# the Bernoulli trials for each observation
ntraps <- nrow(traplocs)
Ymat <- matrix(0, nrow = M, ncol = ntraps)
Ymat[cbind(1:nind, obs.pos)] <- 1

# Data augmentation
# Augment trap.id for individuals not captured....augmented
# individuals have the likelihood evaluated for all points along the
# line and so we set obs.pos = ntraps for those individuals
nz <- M - nind
obs.pos <- c(obs.pos, rep(ntraps, nz) )
uaug <- data[ ,c("u1", "u2")]
uaug <- rbind(uaug, matrix(NA, nrow = nz, ncol = 2))

# Bundle and summarize the data for BUGS
str(data <- list (obs.pos = obs.pos, ntraps = ntraps, traplocs = traplocs,
    nind = nind, y = Ymat, nz = nz, u = uaug) )
# List of 7
# $ obs.pos : num [1:500] 13 121 110 112 159 7 32 120 133 111 ...
# $ ntraps  : int 200
# $ traplocs: num [1:200, 1:2] 0.01 0.03 0.05 0.07 0.09 0.11 ...
# $ nind    : int 141
# $ y       : num [1:500, 1:200] 0 0 0 0 0 0 0 0 0 0 ...
# $ nz      : num 359
# $ u       : num [1:500, 1:2] 0.455 2.489 2.437 2.494 3.444 ...

# Write the BUGS model
cat(file = "transectDS.txt", "
model{

  # Prior distributions
  sigma ~ dunif(0, 10)
  psi ~ dunif(0, 1)
  alpha0 ~ dnorm(0, 0.01)

  for(i in 1:(nind+nz)){
    # Models for DA variables, location and observed data
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(0, 4)
    u[i,2] ~ dunif(0, 1)
    # Compute distance = derived quantity
    for(j in 1:obs.pos[i]){
      d[i,j] <- pow(pow(u[i,1] - traplocs[j,1],2) + pow(u[i,2] - traplocs[j,2],2), 0.5)
      haz[i,j] <- exp(alpha0)*exp(-d[i,j]*d[i,j]/(2*sigma*sigma))# hazard
      p[i,j] <- 1 - exp(-haz[i,j])
      mu[i,j] <- p[i,j]*z[i]
      y[i,j] ~ dbern(mu[i,j]) # Observation model
    }
  }

  # Derived parameters
  N <- sum(z[]) # N is a derived parameter
  D <- N/4 # area = 4 ha
}
")

# MCMC settings
# na <- 500 ; ni <- 2000 ; nb <- 1000 ; nt <- 1 ; nc <- 4
na <- 500 ; ni <- 200 ; nb <- 100 ; nt <- 1 ; nc <- 3  # ~~~~ for testing, 20 mins

# Create inits and define parameters to monitor
ust <- uaug
ust[1:nind,] <- NA
ust[(nind+1):M,] <- cbind( runif(nz, 0, 4), runif(nz, 0, 1) )
inits <- function(){ list (sigma = runif(1, 0.2, 1), psi = runif(1),
    alpha0 = runif(1, -5, -2), z = c(rep(1, nind), rep(0, nz)), u = ust ) }
params <- c("sigma", "N", "psi", "D", "alpha0")

# Run JAGS (ART 71 min), check convergence and summarize posteriors
library(jagsUI)
(out1a <- jags (data, inits, params, "transectDS.txt", n.thin = nt,
    n.chains = nc, n.burnin = nb, n.iter = ni, n.adapt = na, parallel = TRUE))
#           mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# sigma    0.208  0.011   0.189   0.207   0.231    FALSE 1 1.014   207
# N      217.331 18.977 186.000 215.000 262.000    FALSE 1 1.009   363
# psi      0.435  0.044   0.359   0.432   0.535    FALSE 1 1.008   352
# D       54.333  4.744  46.500  53.750  65.500    FALSE 1 1.009   363
# alpha0  -2.175  0.212  -2.611  -2.159  -1.796    FALSE 1 1.008   466

# We next fit the perpendicular distance model
# ''''''''''''''''''''''''''''''''''''''''''''

# Create the data list for JAGS
str(data <- list ( nind = nind, y = c(rep(1,nind), rep(0,nz)), nz = nz,
    u = uaug))
# List of 4
# $ nind: int 141
# $ y   : num [1:500] 1 1 1 1 1 1 1 1 1 1 ...
# $ nz  : num 359
# $ u   : num [1:500, 1:2] 0.455 2.489 2.437 2.494 3.444 ...

# Write BUGS model
cat(file = "transectDSb.txt", "
model{

  # Prior distributions
  sigma ~ dunif(0,10)
  psi ~ dunif(0,1)

  # Models for DA variables, location and observed data
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(0, 4)
    u[i,2] ~ dunif(0, 1)
    # compute distance to transect = derived quantity
    d[i] <- abs(u[i,2]-0.5)
    p[i] <- exp(-d[i]*d[i]/(2*sigma*sigma))
    y[i] ~ dbern(p[i]*z[i])
  }

  # Derived parameters
  N <- sum(z[]) # N is a derived parameter
  D <- N/4 # area = 4 ha
}
")

# Create inits and define parameters to monitor
inits <- function(){ list (sigma = runif(1, 0.2, 1), psi = runif(1),
    z = c(rep(1, nind), rep(0, nz)), u = ust ) }
params <- c("sigma", "N", "psi", "D")

na <- 500 ; ni <- 2000 ; nb <- 1000 ; nt <- 1 ; nc <- 3

# Run JAGS (ART < 1 min), check convergence and summarize posteriors
(out1b <- jags (data, inits, params, "transectDSb.txt", n.thin = nt,
    n.chains = nc, n.burnin = nb, n.iter = ni, n.adapt = na, parallel = TRUE))
#          mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# sigma   0.280  0.033   0.232   0.275   0.357    FALSE 1 1.004   604
# N     219.891 19.280 184.000 219.000 260.000    FALSE 1 1.007   354
# psi     0.440  0.044   0.357   0.439   0.528    FALSE 1 1.007   375
# D      54.973  4.820  46.000  54.750  65.000    FALSE 1 1.007   354
