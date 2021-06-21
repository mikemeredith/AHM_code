# Applied hierarchical modeling in ecology - vol.2 - 2021
# Marc Kéry & J. Andy Royle
#
# Chapter 11 : SPATIALLY EXPLICIT DISTANCE SAMPLING ALONG TRANSECTS
# =================================================================
# Code from proofs dated 2020-08-19

print("Approximate execution time for this code: 80 mins")
# Run time with the full number of iterations: 11.2 hrs

library(AHMbook)
library(jagsUI)

# 11.5 Spatial models for transect sampling
# =========================================

# Simulate a data set and harvest the output
library(AHMbook)
RNGversion("3.5.0")
set.seed(1234 , kind = "Mersenne-Twister")
str(tmp <- simSpatialDSline(N = 200, beta = 1, sigma = 0.25,
    alpha0 = -3))               # produces a figure like 11.1

# 11.5.1 Analyzing the spatial transect model in JAGS
# ---------------------------------------------------

# Harvest the data
data <- tmp$data
nind <- nrow(data)
traplocs <- tmp$traps
ntraps <- nrow(traplocs)
Habitat <- as.vector(tmp$Habitat)
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
nPix <- nrow(Habgrid)

# Do data augmentation, including for pixel ID
M <- 500
nz <- M - nind
obs.pos <- data[,1]             # Location on line where detection was made
dist.obs <- data[,2]            # Distance from observer to object

# Create the encounter matrix: 0 up to 1st capture, then 1
Ymat <- matrix(0, nrow = M, ncol = ntraps)
Ymat[cbind(1:nind, obs.pos)] <- 1

# Grab observation position and we use discrete "pixel ID" here
# instead of continuous location.
# Augment for individuals not captured....
obs.pos <- c(obs.pos, rep(ntraps, M-nind) )
pixel <- c(tmp$pixel, rep(NA, M - nind))

# Bundle and summarize the data for BUGS
str(bdata <- list(obs.pos = obs.pos, ntraps = ntraps, traplocs = traplocs,
    nind = nind, y = Ymat, nz = nz, Habitat = Habitat, Habgrid = Habgrid,
    nPix = nPix, pixel = pixel)) # not shown

# Write model in BUGS language
cat(file = "spatialDS.txt", "
model{

  # Prior distributions
  sigma ~ dunif(0,10)
  psi ~ dunif(0,1)
  beta1 ~ dnorm(0,0.01)
  alpha0 ~ dnorm(0,0.01)
  for(g in 1:nPix){ # Discrete point process model
    probs.num[g] <- exp(beta1*Habitat[g])
    probs[g] <- probs.num[g]/sum(probs.num[])
  }
  # Models for DA variables and location (pixel)
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    pixel[i] ~ dcat(probs[])
    s[i,1:2] <- Habgrid[pixel[i],]            # location = look-up in table
    for(j in 1:obs.pos[i]){
      # compute distance == a derived quantity
      d[i,j] <- pow( pow( s[i,1] - traplocs[j,1],2) +
      pow(s[i,2] - traplocs[j,2],2), 0.5)
      haz[i,j] <- exp(alpha0)*exp(-d[i,j]*d[i,j]/(2*sigma*sigma))
      p[i,j] <- 1 - exp(-haz[i,j])
      mu[i,j] <- p[i,j]*z[i]
      y[i,j] ~ dbern(mu[i,j])                  # Observation model
    }
  }
  # Derived parameters
  N <- sum(z[]) # N is a derived parameter
  D <- N/4 # area = 4 units
}
")

# Load jagsUI and specify MCMC settings
# na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nt <- 2 ; nc <- 5
na <- 1000 ; ni <- 1200 ; nb <- 200 ; nt <- 1 ; nc <- 3  # ~~~ for testing, 40 mins

# Create inits
inits <- function(){ list (sigma = runif(1, 0.2, 1), psi = runif(1),
    alpha0 = runif(1, -5, -2), z = c(rep(1, nind), rep(0, nz)) ) }

# Parameters to monitor
params1 <- c("sigma", "N", "psi", "beta1", "D", "alpha0")
params2 <- c("sigma", "N", "psi", "beta1", "D", "alpha0", "pixel", "z")

# Run JAGS (ART 377 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params2, "spatialDS.txt", n.thin = nt,
    n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE,
    factories = "base::Finite sampler FALSE")
round(out3$summary[params1, ], 2) # Just main structural parameters
# par(mfrow = c(3,2))  # ~~~ replaced with 'layout' argument
traceplot(out3, params1, layout=c(3,2))

#          mean    sd   2.5%    25%    50%    75%  97.5% Rhat n.eff overlap0 f
# sigma    0.25  0.02   0.22   0.23   0.24   0.26   0.28 1.00  4295        0 1
# N      236.32 74.53 154.00 185.00 213.00 262.00 444.00 1.03   156        0 1
# psi      0.47  0.15   0.30   0.37   0.43   0.53   0.89 1.03   157        0 1
# beta1    1.29  0.12   1.06   1.21   1.29   1.37   1.52 1.00  2047        0 1
# D       59.08 18.63  38.50  46.25  53.25  65.50 111.00 1.03   156        0 1
# alpha0  -3.22  0.44  -4.21  -3.49  -3.15  -2.90  -2.53 1.02   224        0 1

# ~~~~ extra code for figure 11.4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pixel <- out3$sims.list$pixel*out3$sims.list$z
post <- table(pixel)/nrow(pixel)
post <- post[-1] # remove the 0 values, z = 0 individuals

library(raster)
op <- par(mfrow=c(2,1),mar=c(2,2,2,6))
prior.mean <- mean(out3$sims.list$beta)*as.vector(Habitat)
prior.mean <- mean(out3$sims.list$psi)*M*exp(prior.mean)/sum(exp(prior.mean))
plot(rast.prior<- rasterFromXYZ(cbind(Habgrid,prior.mean)),axes=FALSE,box=FALSE,
    col=rampYOR(225, bias=2), ylim=c(0,1))
title("Prior mean density (estimated)", line=0)
axis(1)
axis(2)
plot(rast.post <- rasterFromXYZ(cbind(Habgrid,as.vector(post))),axes=FALSE,
    box=FALSE, col=rampYOR(225, bias=2))
title("Posterior mean density", line=0)
axis(1)
axis(2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 11.5.2 Formulating a continuous-space model with a discrete raster
# ------------------------------------------------------------------

library(AHMbook)
set.seed(1234, kind = "Mersenne-Twister")
# Simulate with perpendicular distance
str(tmp <- simSpatialDSline(N = 200, beta = 1, sigma = 0.25, alpha0 = -3,
    perp = TRUE))

# Harvest the data
xlim <- c(0, tmp$L)
ylim <- c(0, tmp$W*2)
delta <- tmp$delta
u <- tmp$data[, c("u1", "u2")]
nind <- nrow(u)

# Convert Habraster to a matrix we can use in JAGS
library(raster)
tmpmat <- as.matrix(tmp$Habraster)
Habmat <- t(tmpmat[nrow(tmpmat):1, ]) # flip and transpose as per Fig. 11.5.
# Centre (no need to scale)
Habmat <- Habmat - mean(Habmat)
str(Habmat)

# Write the BUGS model file
# DS model with continuous-space prior for locations and using
# perpendicular distance and with habitat as a raster
cat(file = "spatialDSfast_perp.txt","
model{

  # Prior distributions
  sigma ~ dunif(0,10)
  beta1 ~ dnorm(0,0.01)
  psi ~ dunif(0,1)
  # lambda and probs are matrices the same dimensions as Habmat
  for(d1 in 1:Dim[1]){
    for(d2 in 1:Dim[2]){
      lam[d1, d2] <- exp(beta1*Habmat[d1, d2])
    }
  }
  probs <- lam/sum(lam)

  # Models for DA variables and location (pixel)
  for(i in 1:M){
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(xlim[1], xlim[2])
    u[i,2] ~ dunif(ylim[1], ylim[2])
    # Converts continuous location to indices
    xindex[i] <- round(u[i,1]/delta+0.5)
    yindex[i] <- round(u[i,2]/delta+0.5)
    # This is the zeros trick
    negLogDen[i] <- -log(probs[xindex[i], yindex[i]]) # neg. log-density
    zeros[i] ~ dpois(negLogDen[i])
    # Compute distance. Centerline is at 0.5
    d[i] <- pow( pow( u[i,2] - 0.5, 2), 0.5)
    p[i] <- 1.0*exp(-d[i]*d[i]/(2*sigma*sigma)) # p(0) must be 1
    mu[i] <- p[i]*z[i]
    y[i] ~ dbern(mu[i]) # Observation model
  }
  # Derived parameters
  N <- sum(z[]) # N is a derived parameter
  D <- N/4 # area = 4 units
}
")

# Do data augmentation
M <- 500
nz <- M - nind
y <- c(rep(1,nind), rep(0,nz))

# Augment location and pick starting values and inits
uaug <- rbind(u, matrix(NA, nrow = nz, ncol = 2))
uinit <- rbind( matrix(NA, nrow = nind, ncol = 2),
cbind( runif(nz, xlim[1], xlim[2]), runif(nz, ylim[1], ylim[2]) ) )

# Specify MCMC settings
ni <- 2000 ; nb <- 500 ; nt <- 2 ; nc <- 4

# Create data bundle, inits and list of parameters to monitor
data_perp <- list (y = y, M = M, Habmat = Habmat, Dim = dim(Habmat),
    xlim = xlim, ylim = ylim, u = uaug, delta = delta, zeros = rep(0, M) )
inits <- function(){ list (sigma = runif(1, 0.2, 1), beta1 = rnorm(1, 0, 0.2),
    z = c(rep(1,nind), rep(0,nz)), u = uinit ) }
params <- c("sigma", "N", "psi", "beta1", "D")

# Run JAGS (ART < 1 min), check convergence and summarize posteriors
out1perp <- jags (data_perp, inits, params, "spatialDSfast_perp.txt",
    n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
# par(mfrow = c(3,2))  # ~~~ replaced with 'layout' argument
traceplot(out1perp, layout=c(3,2))
print(out1perp, 3)
#          mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# sigma   0.244  0.025   0.203   0.241   0.300    FALSE 1 1.002  2438
# N     208.888 20.223 173.000 208.000 254.000    FALSE 1 1.003  1334
# psi     0.418  0.046   0.336   0.415   0.516    FALSE 1 1.003   871
# beta1   1.133  0.106   0.932   1.134   1.335    FALSE 1 1.015   184
# D      52.222  5.056  43.250  52.000  63.500    FALSE 1 1.003  1334

# Next, the hazard encounter model.
# '''''''''''''''''''''''''''''''''

# Generate another data set using the hazard detection model
set.seed(1234, kind = "Mersenne-Twister")
str(tmp <- simSpatialDSline(N = 200, beta = 1, sigma = 0.25, alpha0 = -3,
    perp = FALSE))

# Harvest data
xlim <- c(0, tmp$L)
ylim <- c(0, tmp$W*2)
delta <- tmp$delta
data <- tmp$data
u <- data[, c("u1", "u2")]
nind <- nrow(u)
traplocs <- tmp$traps
ntraps <- nrow(traplocs)
obs.pos <- data[,1]
dist.obs <- data[,2]

# Convert Habraster to a matrix we can use in JAGS
tmpmat <- as.matrix(tmp$Habraster)
Habmat <- t(tmpmat[nrow(tmpmat):1, ])
Habmat <- Habmat - mean(Habmat) # Centre (no need to scale)
str(Habmat)

# Data augmentation, including for pixel ID
M <- 500
nz <- M - nind

# Data augmentation
uaug <- rbind(u, matrix(NA, nrow = nz, ncol = 2))
Ymat <- matrix(0, nrow = M, ncol = ntraps)
Ymat[cbind(1:nind, obs.pos)] <- 1

# Augment observer position for individuals not captured
obs.pos <- c(obs.pos, rep(ntraps, M - nind) )

# Write out the BUGS model
cat(file = "spatialDSfast_haz.txt", "
model{
  # Prior distributions
  sigma ~ dunif(0,10)
  alpha0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  psi ~ dunif(0, 1)
  for(d1 in 1:Dim[1]){
    for(d2 in 1:Dim[2]){
      lam[d1, d2] <- exp(beta1*Habmat[d1, d2])
    }
  }
  probs <- lam/sum(lam)

  # Models for DA variables and location
  for(i in 1:M){
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(xlim[1], xlim[2])
    u[i,2] ~ dunif(ylim[1], ylim[2])
    xindex[i] <- round(u[i,1]/delta+0.5)
    yindex[i] <- round(u[i,2]/delta+0.5)
    negLogDen[i] <- -log(probs[xindex[i], yindex[i]])     # neg.log-density
    zeros[i] ~ dpois(negLogDen[i])                        # zeros trick
    # Evaluate likelihood
    for(j in 1:obs.pos[i]){
      d[i,j] <- pow( pow( u[i,1] - traplocs[j,1],2) +
      pow(u[i,2] - traplocs[j,2],2), 0.5)
      haz[i,j] <- exp(alpha0)*exp(-d[i,j]*d[i,j]/(2*sigma*sigma))
      p[i,j] <- 1 - exp(-haz[i,j])
      mu[i,j] <- p[i,j]*z[i]
      y[i,j] ~ dbern(mu[i,j])                             # Observation model
    }
  }
  # Derived parameters
  N <- sum(z[])
  D <- N/4
}
")

# MCMC settings
# na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nt <- 2 ; nc <- 4
na <- 1000 ; ni <- 1200 ; nb <- 200 ; nt <- 2 ; nc <- 3  # ~~~~ for testing

# Bundle up the data. Include vector of 0s for the zeros trick.
str(data_haz <- list (y = Ymat, M = M, Habmat = Habmat, Dim = dim(Habmat),
    xlim = xlim, ylim = ylim, u = uaug, delta = delta, zeros = rep(0, M),
    traplocs = traplocs, obs.pos = obs.pos) )

# Create inits and define parameters to monitor
inits <- function(){ list (sigma = runif(1, 0.2, 1), beta1 = rnorm(1, 1, 0.4),
    alpha0 = runif(1, -5, -2), z = c(rep(1, nind), rep(0, nz)) ) }
params <- c("sigma", "N", "psi", "beta1", "D", "alpha0")

# Run JAGS (ART 342 min) and summarize. Note: factories setting not required.
out1haz <- jags (data_haz, inits, params, "spatialDSfast_haz.txt",
    n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
# par(mfrow = c(3,2))  # ~~~ replaced with 'layout' argument
traceplot(out1haz, layout=c(3,2))
print(out1haz, 3)
#           mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# sigma    0.256  0.017   0.226   0.255   0.292    FALSE 1 1.000 20000
# N      191.585 48.747 138.000 178.000 330.000    FALSE 1 1.018   684
# psi      0.383  0.099   0.266   0.359   0.662    FALSE 1 1.016   714
# beta1    0.926  0.111   0.707   0.927   1.141    FALSE 1 1.001 10268
# D       47.896 12.187  34.500  44.500  82.500    FALSE 1 1.018   684
# alpha0  -3.044  0.391  -3.983  -2.985  -2.421    FALSE 1 1.004  1568
