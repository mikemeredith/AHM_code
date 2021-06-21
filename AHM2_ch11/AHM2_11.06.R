# Applied hierarchical modeling in ecology - vol.2 - 2021
# Marc Kéry & J. Andy Royle
#
# Chapter 11 : SPATIALLY EXPLICIT DISTANCE SAMPLING ALONG TRANSECTS
# =================================================================
# Code from proofs dated 2020-08-19

print("Approximate execution time for this code: 60 mins")
# Run time with the full number of iterations: 2.8 hrs

library(AHMbook)
library(jagsUI)

# 11.6 Hierarchical distance sampling (HDS) with multiple transects
# =================================================================

library(AHMbook)
RNGversion("3.5.3")
set.seed(1234, kind = "Mersenne-Twister")

# Have to know the value of nPix and ntraps. Be careful!
nPix <- 400
ntran <- 2    # Number of transects
ntraps <- 200 # Number of 'traps' per transect
M <- 400      # Data augmentation
y2d <- pixel <- matrix(NA, nrow = M, ncol = ntran) # note: Now matrices!
umat <- array(NA, dim = c(M, 2, ntran))
Ymat <- array(0, dim = c(M, ntraps, ntran))
obs.pos <- matrix(NA, nrow = M, ncol = ntran)
nind <- rep(NA, ntran)
Habitat <- matrix(NA, nrow = nPix, ncol = ntran)

# Simulate multiple transects (produces Figure 11.6)
for(tran in 1:ntran){
  tmp <- simSpatialDSline(N = 200, beta = 1, sigma = 0.25, alpha0 = -3)
  # Harvest data
  data <- tmp$data
  u <- data[ , c("u1", "u2")]
  trap <- data[,1]      # Integer for line segment
  dist.obs <- data[,2]  # Distance from observer position to object
  traplocs <- tmp$traps # Constant for standardized transects
  Habitat[,tran] <- as.vector(tmp$Habitat)
  Habitat[,tran] <- Habitat[,tran] - mean(Habitat[,tran])
  Habgrid <- tmp$grid   # Always same for standardized transects grid
  nind[tran] <- nrow(u)

  # Do data augmentation, including for pixel ID
  nz <- M - nind[tran]
  y2d[, tran] <- c(rep(1,nind[tran]), rep(0,nz))
  uaug <- rbind(u, matrix(NA, nrow = nz, ncol = 2))
  umat[1:M, 1:2, tran] <- uaug

  # Fill-out the Ymat with pixel of each observation
  pixel[1:nind[tran], tran] <- tmp$pixel
  Ymat[cbind(1:nind[tran], trap, tran)] <- 1

  # Augment for individuals not captured
  obs.pos[, tran] <- c(trap, rep(ntraps, M - nind[tran]) )
}

# Bundle and summarize the data for BUGS
str(data <- list (obs.pos = obs.pos, ntraps = ntraps, ntran = ntran,
    traplocs = traplocs, nind = nind, y = Ymat, M = M, Habitat = Habitat,
    Habgrid = Habgrid, nPix = nPix, pixel = pixel))
# List of 11
# $ obs.pos : num [1:400, 1:2] 28 99 105 72 31 111 153 143 25 2 ...
# $ ntraps  : num 200
# $ ntran   : num 2
# $ traplocs: num [1:200, 1:2] 0.01 0.03 0.05 0.07 0.09 0.11 0.13 0.15 0.17 ...
# $ nind    : int [1:2] 105 120
# $ y       : num [1:400, 1:200, 1:2] 0 0 0 0 0 0 0 0 0 0 ...
# $ M       : num 400
# $ Habitat : num [1:400, 1:2] -0.375 -0.142 0.413 -0.546 -0.232 ...
# $ Habgrid : num [1:400, 1:2] 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 ...
# $ nPix    : num 400
# $ pixel   : int [1:400, 1:2] 206 143 301 178 169 225 270 307 329 323 ...

# Write BUGS model
cat(file = "spatialDSmulti.txt", "
model{

  # Prior distributions
  sigma ~ dunif(0,10)
  alpha0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  psi[1] ~ dunif(0,1)
  psi[2] ~ dunif(0,1)

  # Note that probabilities are computed for each transect now
  for(tran in 1:ntran){
    for(g in 1:nPix){                           # g is the pixel index
      lam[g,tran] <- exp(beta1*Habitat[g,tran])
      probs[g,tran] <- lam[g,tran]/sum(lam[,tran])
    }
  }

  # Likelihood and spatial model for each transect
  for(tran in 1:ntran){
    for(i in 1:M){
      z[i,tran] ~ dbern(psi[tran])
      pixel[i,tran] ~ dcat(probs[,tran])
      s[i,1:2,tran] <- Habgrid[pixel[i,tran],]  # location = derived quantity
      # compute distance = derived quantity
      for(j in 1:obs.pos[i,tran]){
        d[i,j,tran] <- pow(pow(s[i,1,tran] - traplocs[j,1],2) +
        pow(s[i,2,tran] - traplocs[j,2],2), 0.5)
        haz[i,j,tran] <- exp(alpha0)*exp(-d[i,j,
        tran]*d[i,j,tran]/(2*sigma*sigma))      # Half-normal hazard det. fctn.
        p[i,j,tran] <- 1 - exp(-haz[i,j,tran])
        mu[i,j,tran] <- p[i,j,tran]*z[i,tran]
        y[i,j,tran] ~ dbern(mu[i,j,tran])       # Observation model
      }
    }
    # Derived parameters
    N[tran] <- sum(z[,tran])                    # N is a derived parameter
    D[tran] <- N[tran]/4                        # area = 4 units
  }
}
")

# MCMC settings
# na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nt <- 1 ; nc <- 6
na <- 1000 ; ni <- 300 ; nb <- 100 ; nt <- 1 ; nc <- 3  # ~~~~ for testing

# Create inits and define parameters to monitor
inits <- function(){ list (sigma = runif(1,0.2,1), beta1 = rnorm(1, 1, 0.4),
    alpha0 = runif(1, -5, -2), z = y2d ) }
params <- c("sigma", "N", "psi", "beta1", "D", "alpha0")

# Run JAGS (ART 160 min), check convergence and summarize posteriors
out5 <- jags (data, inits, params, "spatialDSmulti.txt", n.adapt = na,
    n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE,
    factories = "base::Finite sampler FALSE")
# par(mfrow = c(3,2))  # ~~~ replaced with 'layout' argument
traceplot(out5, layout=c(3,2))
print(out5, 3)
#          mean    sd   2.5%    50%  97.5% overlap0 f Rhat n.eff
# sigma    0.25  0.01   0.23   0.25   0.27    FALSE 1 1.01   505
# N[1]   167.77 16.97 140.00 166.00 206.00    FALSE 1 1.01   314
# N[2]   193.45 19.43 162.00 191.00 239.00    FALSE 1 1.01   274
# psi[1]   0.42  0.05   0.33   0.42   0.53    FALSE 1 1.01   410
# psi[2]   0.48  0.05   0.39   0.48   0.60    FALSE 1 1.01   307
# beta1    1.00  0.09   0.83   1.00   1.17    FALSE 1 1.02   179
# D[1]    41.94  4.24  35.00  41.50  51.50    FALSE 1 1.01   314
# D[2]    48.36  4.86  40.50  47.75  59.75    FALSE 1 1.01   274
# alpha0  -2.73  0.19  -3.13  -2.72  -2.38    FALSE 1 1.00   970
