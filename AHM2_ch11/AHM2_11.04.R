# Applied hierarchical modeling in ecology - vol.2 - 2021
# Marc Kéry & J. Andy Royle
#
# Chapter 11 : SPATIALLY EXPLICIT DISTANCE SAMPLING ALONG TRANSECTS
# =================================================================
# Code from proofs dated 2020-08-19

print("Approximate execution time for this code: 25 mins")
# Run time with the full number of iterations: 75 mins

library(AHMbook)
library(jagsUI)

# 11.4 Mark-recapture/distance sampling models on linear transects
# ================================================================

# Simulation settings
set.seed( 1234, kind = "Mersenne-Twister")
N <- 200
M <- 500
sigma <- 0.20
alpha0 <- -3
W <- 1/2                                  # Transect half-width
L <- 4
K <- 2                                    # Replicates

# Locations of individuals
u1 <- runif(N, 0, L)
u2 <- runif(N, 0, 2*W)
plot(u1, u2, pch = 20, cex = 1)           # plot points (not shown)
abline(0.5, 0, lwd = 2, col = 'grey')
title("Transect population subject to detection")

# Represent transect by a sequence of points with delta = 0.2
line.pts <- seq(0.01, L - 0.01, .02)
traplocs <- cbind(line.pts, 0.5) # Call these points "traps"

# Initialize some objects. These are matrices now
d.to.trap <- obs.pos <- pbar <- matrix(NA, nrow = N, ncol = K)

# Simulate the detections of each individual as a sequence of Bernoulli
# trials and record only the first detection (where it was detected from on
# the line)
dmat <- e2dist(cbind(u1, u2), traplocs)
for(k in 1:K){
  for(i in 1:nrow(dmat)){
    haz <- exp(alpha0)*exp( -(dmat[i,]^2)/(2*sigma*sigma))
    probs <- 1 - exp(-haz)
    captured <- rbinom(nrow(traplocs), 1, probs)    # Seq. of Bern trials
    pbar[i,k] <- 1 - exp(-sum(haz))                 # Average prob. of capture
    if(sum(captured)==0)
      next
    obs.pos[i,k] <- (1:length(captured))[captured==1][1]
    d.to.trap[i,k] <- dmat[i,][obs.pos[i,k]]
    lines(c(u1[i], traplocs[obs.pos[i,k], 1]), c(u2[i],
    traplocs[obs.pos[i,k], 2]) )
  }
}

# Subset to detected individuals
ncap <- apply(!is.na(d.to.trap), 1, sum)
nind <- sum(ncap > 0)
data <- cbind(obs.pos, u1, u2)[ncap>0, ]
obs.pos <- data[, 1:2]              # K = 2 vectors here, one for each survey

ntraps <- nrow(traplocs)
Yarr <- array(0, dim = c(M, ntraps, K) )
for(k in 1:K){
  for(i in 1:nind){
    if(!is.na(obs.pos[i,k])) {
      Yarr[i, obs.pos[i,k], k] <- 1
    }
  }
}

# Data augmentation
nz <- M - nind
obs.pos[is.na(obs.pos)] <- ntraps
obs.pos <- rbind(obs.pos, matrix(ntraps, nrow = nz, ncol = K) )

# Augment location matrix u
uaug <- rbind( data[,c("u1", "u2")], matrix(NA, nrow = nz, ncol = 2))

# Bundle and summarize the data for BUGS
str(data_haz <- list (obs.pos = obs.pos, ntraps = ntraps, traplocs = traplocs,
    nind = nind, y = Yarr, nz = nz, u = uaug, nsurveys = K) )# output omitted

# Model 2a: MRDS with hazard detection model
cat(file="MRDS.txt", "
model {

  # Prior distributions
  sigma ~ dunif(0,10)
  psi ~ dunif(0,1)
  alpha0 ~ dnorm(0,0.01)

  # Models for DA variables and location (pixel)
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(0, 4)
    u[i,2] ~ dunif(0, 1)
    for(k in 1:nsurveys){
      for(j in 1:obs.pos[i,k]){
        d[i,j,k] <- pow( pow( u[i,1] - traplocs[j,1],2) +
        pow(u[i,2] - traplocs[j,2],2),0.5)
        haz[i,j,k] <- exp(alpha0)*exp(-d[i,j,k]*d[i,j,k]/ (2*sigma*sigma))
        p[i,j,k] <- 1 - exp(-haz[i,j,k])
        mu[i,j,k] <- p[i,j,k]*z[i]
        y[i,j,k] ~ dbern(mu[i,j,k])               # Observation model
      }
    }
  }
  # Derived parameters
  N <- sum(z[]) # N is a derived parameter
  D <- N/4 # area = 4 ha
}
")

# MCMC settings
# na <- 500 ; ni <- 1500 ; nb <- 500 ; nt <- 1 ; nc <-5
na <- 500 ; ni <- 150 ; nb <- 50 ; nt <- 1 ; nc <- 3  # ~~~ for testing, 25 mins

# Create inits and define parameters to monitor
ust <- uaug
ust[1:nind, ] <- NA
ust[(nind+1):M, ] <- cbind( runif(nz, 0, L), runif(nz, 0, W) )
inits <- function(){ list (sigma = runif(1, 0.2, 1), psi = runif(1),
    alpha0 = runif(1, -5, -2), z = c(rep(1, nind), rep(0, nz)), u = ust ) }
params <- c("sigma", "N", "psi", "D", "alpha0")

# Run JAGS (ART 95 min), check convergence and summarize posteriors
out2a <- jags (data_haz, inits, params, "MRDS.txt", n.thin = nt,
    n.chains = nc, n.burnin = nb, n.iter = ni, n.adapt = na, parallel = TRUE)
# par(mfrow = c(3,2))  # ~~~ replaced with 'layout' argument
traceplot(out2a , layout=c(3,2))
print(out2a, 2)
#          mean    sd   2.5%    50%  97.5% overlap0 f Rhat n.eff
# sigma    0.20  0.01   0.18   0.19   0.21    FALSE 1 1.00   762
# N      214.98 19.07 181.00 214.00 256.00    FALSE 1 1.01   388
# psi      0.43  0.04   0.35   0.43   0.52    FALSE 1 1.01   471
# D       53.74  4.77  45.25  53.50  64.00    FALSE 1 1.01   388
# alpha0  -3.06  0.15  -3.37  -3.06  -2.78    FALSE 1 1.01   236
