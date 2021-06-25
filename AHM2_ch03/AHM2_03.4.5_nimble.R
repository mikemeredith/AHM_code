#   Applied hierarchical modeling in ecology - vol.2 - 2021
#   Marc Kéry & J. Andy Royle
#
# Chapter 3 : HIERARCHICAL MODELS OF SURVIVAL
# ===========================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 14 mins
# Run time with the full number of iterations: 4.1 hrs

##### USING NIMBLE INSTEAD OF WinBUGS #####

# Code put together by Mike Meredith

library(AHMbook)
library(nimble)
library(mcmcOutput)
library(spdep)

# ~~~~~ data wrangling code from 3.4.1 ~~~~~~~~~~~~~~
data(willowWarbler)
ch <- as.matrix(willowWarbler$birds[, 1:11])
sitevec <- willowWarbler$birds$cesID
nyear <- ncol(ch) # Number years: 11
nsite <- nrow(willowWarbler$CES) # Number of CE sites: 193
nblock <- nrow(willowWarbler$blocks) # Number of blocks: 495
marr <- ch2marray(ch)
MARR <- array(NA, dim = c(10, nyear, nsite))
R <- array(NA, dim = c(10, nsite))
for(k in 1:nsite){
  sel.part <- ch[sitevec == k, ]
  ma <- ch2marray(sel.part)
  MARR[,,k] <- ma
  R[,k] <- apply(ma, 1, sum)
}
# ~~~~ and this bit from 3.4.3 ~~~~~
# Scale linear and squared GDD separately and latitude (for entire grid)
scaled.gdd1 <- standardize(willowWarbler$cells$gdd)
scaled.gdd2 <- standardize(willowWarbler$cells$gdd^2)
scaled.lat <- standardize(willowWarbler$cells$lat)
# Pull out GDD and LATITUDE values for 193 sites (using CellID)
gdd1.site <- scaled.gdd1[willowWarbler$CES$CellID]
gdd2.site <- scaled.gdd2[willowWarbler$CES$CellID]
lat.site <- scaled.lat[willowWarbler$CES$CellID]
# ~~~~ and this bit from 3.4.4 ~~~~~
MARRWB <- aperm (MARR, c(3, 1, 2)) # MARR for WinBUGS
blockcoordgrid <- cbind(as.matrix(willowWarbler$blocks))
neigh <- dnearneigh(blockcoordgrid, d1 = 0, d2 = sqrt(2 * 25^2)+ 0.1)
winnb <- nb2WB(neigh) # Function to get CAR ingredients for BUGS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 3.4 Spatial hierarchical CJS models
# ===================================

# 3.4.5 Fitting a spatial hierarchical CJS model with spatial autocorrelation and with covariates
# -----------------------------------------------------------------------------------------------

# Bundle and summarize data set for NIMBLE
bdata <- list(MARRWB = MARRWB)
str(bconst <- list(R = R, n.site = nsite, n.occ = nyear,
    n.block = nblock, BlockID = willowWarbler$CES$BlockID, adj = winnb$adj,
    weights = winnb$weights, num = winnb$num, gdd1.site = gdd1.site,
    gdd2.site = gdd2.site, lat.site = lat.site))
# List of 11
# $ R         : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site    : int 193
# $ n.occ     : int 11
# $ n.block   : int 495
# $ BlockID   : int [1:193] 25 77 204 110 222 283 119 234 295 152 ...
# $ adj       : int [1:3458] 2 3 4 1 3 5 6 1 2 4 ...
# $ weights   : num [1:3458] 1 1 1 1 1 1 1 1 1 1 ...
# $ num       : int [1:495] 3 4 6 5 3 6 6 6 5 5 ...
# $ gdd1.site : num [1:193] 1.13 0.726 0.46 0.538 0.686 ...
# $ gdd2.site : num [1:193] 1.228 0.703 0.38 0.473 0.653 ...
# $ lat.site  : num [1:193] -1.397 -1.11 -0.403 -0.881 -0.308 ...

# Specify model in NIMBLE dialect of BUGS language
# ''''''''''''''''''''''''''''''''''''''''''''''''
# Necessary changes:
# 1. Specify 3rd dimension of pr[t,s, ] -> pr[t,s, 1:n.occ]
# 2. Remove I(-12,12) as not compatible (or necessary)
# 3. dcar_normal needs zero_mean=1 argument to mimic WinBUGS version
# Optional improvements:
# 1. Use logit(.) <-
# 2. sd instead of tau in dnorm
# 3. dbeta for probabilities instead of dunif

cjs9.code <- nimbleCode({
  # Priors and linear models
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      logit(phi[t, s]) <- lphi[t, s] # survival
      logit(p[t, s]) <- lp[t, s]     # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    # Linear model for site-level effects: add covariates
    alpha.lphi.site[s] <- mu.lphi + beta1 * gdd1.site[s] + beta2 * gdd2.site[s] + beta3 *
    lat.site[s] + eta[BlockID[s]]
    alpha.lp.site[s] ~ dnorm(mu.lp, sd=sd.lp.site)

    # backtransform site means
    logit(mean.p.site[s]) <- alpha.lp.site[s]
  }
  for (t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, sd=sd.lphi.time)
    beta.lp.time[t] ~ dnorm(0, sd=sd.lp.time)

    # backtransform time means (see Errata dated 2021-06-23)
    logit(mean.phi.time[t]) <- mu.lphi + beta.lphi.time[t]
    logit(mean.p.time[t]) <- mu.lp + beta.lp.time[t]
  }

  # Hyperpriors for hyperparams
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dbeta(1, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  sd.lp.site ~ dunif(0, 2)
  sd.lphi.time ~ dunif(0, 1)
  sd.lp.time ~ dunif(0, 1)

  # Coefficients for gdd1, gdd2 and lat
  beta1 ~ dunif(-3, 3)
  beta2 ~ dunif(-3, 3)
  beta3 ~ dunif(-3, 3)

  # CAR prior distribution for spatial random effects
  # NOTE: this is defined on the entire grid of 495 blocks
  # eta[1:n.block] ~ car.normal(adj[], weights[], num[], tau)
  eta[1:n.block] ~ dcar_normal(adj[], weights[], num[], tau, zero_mean=1)
  tau ~ dgamma(0.5, 0.0005) # cf. Lunn et al. (2013)
  veta <- 1/tau
  sdeta <- sqrt(veta)

  # Multinomial likelihood for the m-array data (WinBUGS style)
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      MARRWB[s, t,1:n.occ] ~ dmulti(pr[t,s, 1:n.occ], R[t,s])
    }
  }

  # Define the cell probabilities of the m-array
  # Main diagonal
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      q[t,s] <- 1-p[t,s] # Probability of non-recapture
      pr[t,s,t] <- phi[t,s]*p[t,s]
      # Above main diagonal
      for (j in (t+1):(n.occ-1)){
        pr[t,s,j] <- prod(phi[t:j,s])*prod(q[t:(j-1),s])*p[j,s]
      } #j
      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,s,j] <- 0
      } #j
    } #t
  } #s
  # Last column of m-array: probability of non-recapture
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      pr[t,s,n.occ] <- 1-sum(pr[t,s,1:(n.occ-1)])
    } #t
  }#s
})

# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1),
    eta = rep(0, nblock), beta1 = rnorm(1)/3, beta2 = rnorm(1)/3,
    beta3 = rnorm(1) /3, tau = runif(1))}

# Parameters monitored
params <- c("mean.phi", "mean.p", "mu.lphi", "mu.lp", "sd.lp.site", "sd.lphi.time",
    "sd.lp.time", "mean.p.site", "mean.phi.time", "mean.p.time", "veta", "sdeta",
    "beta1", "beta2", "beta3", "eta")

# MCMC settings
ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 2 # ~~~ for testing, 10 mins

system.time(
out <- nimbleMCMC(cjs9.code, data=bdata, constants=bconst,
    inits=inits, monitors=params,
    niter=ni, nburnin=nb, thin=nt, nchains=nc,
    samplesAsCodaMCMC=TRUE) )

mc9 <- mcmcOutput(out)
View(summary(mc9, n.eff=TRUE))
diagPlot(mc9)


# Run multiple NIMBLE instances in parallel with 'foreach'
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''
library(foreach)
library(parallel)
library(doParallel)

detectCores()  # number of cores on your machine
ncore <- 3     # <-- adjust this for your machine (3 used for testing)

cl <- makeCluster(ncore)
registerDoParallel(cl)

# ni <- 2e5 ; nt <- 100 ; nb <- 1e5 # 4 hrs
ni <- 1000  ;  nt <- 1  ;  nb <- 500  # ~~~ testing, 3.5 mins

seeds <- 1:ncore

system.time(
res9 <- foreach(x = seeds, .packages="nimble",
      .errorhandling='remove', .inorder=FALSE) %dopar% {
  set.seed(x)
  nimbleMCMC(cjs9.code, data=bdata, constants=bconst,
      inits=inits, monitors=params,
      niter=ni, nburnin=nb, thin=nt, nchains=1,
      samplesAsCodaMCMC=TRUE)
} )
stopCluster(cl)

# How many chains completed successfully?
length(res9)

# Convert to an mcmcOutput object and look at diagnostics
mclist <- coda::mcmc.list(res9)
( mco9 <- mcmcOutput(mclist) )
diagPlot(mco9[c("mean.phi", "mean.p", "mu.lphi", "mu.lp", "sd.lp.site", "sd.lphi.time",
    "sd.lp.time", "veta", "sdeta", "beta1", "beta2", "beta3")])
View(summary(mco9))

mean(mco9$beta1 > 0)
mean(mco9$beta2 < 0)
mean(mco9$beta3 > 0)
# [1] 0.9696667 # Prob gdd1 has positive effect
# [1] 0.9423333 # Prob gdd2 has negative effect
# [1] 0.9993333 # Prob lat has positive effect
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~ code for figure 3.15 ~~~~~~~~~~~
sims <- mco9            # Parallel output
n.sims <- nrow(mco9)

# Figure 3.15
op <- par(mfrow = c(1, 3))#, mar = c(5,6,5,4), cex.main = 2, cex.axis = 2, cex.lab = 2)
plot(density(sims$beta1), main = 'gdd (linear)', col = 'gray', frame = FALSE, xlab = 'beta1')
abline(v = 0)
plot(density(sims$beta2), main = 'gdd (squared)', col = 'gray', frame = FALSE, xlab = 'beta2')
abline(v = 0)
plot(density(sims$beta3), main = 'Latitude', col = 'gray', frame = FALSE, xlab = 'beta3')
abline(v = 0)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for figures 3.16 and 3.17 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ncell <- nrow(willowWarbler$cells)

# Predictions 1: effects of covariates only
phi.pred1<- array(NA, dim = c(ncell, n.sims))
for (s in 1:ncell){
  phi.pred1[s,] <- plogis(sims$mu.lphi + sims$beta1 * scaled.gdd1[s] +
      sims$beta2 * scaled.gdd2[s] + sims$beta3 * scaled.lat[s])
}
# Predictions 2: effects of spatial autocorrelation only
phi.pred2<- array(NA, dim = c(nblock, n.sims))
for (b in 1:nblock){
  phi.pred2[b,] <- plogis(sims$mu.lphi + sims$eta[,b])
}
# Predictions 3: effects of both covariates and neighbourhood relations
phi.pred3<- array(NA, dim = c(ncell, n.sims))
for (s in 1:ncell){
  phi.pred3[s,] <- plogis(sims$mu.lphi + sims$beta1 * scaled.gdd1[s] +
      sims$beta2 * scaled.gdd2[s] + sims$beta3 * scaled.lat[s] +
      sims$eta[, willowWarbler$cells$blockID[s]])
}

# Compute posterior means and sds
phi.pred1.pm <- apply(phi.pred1, 1, mean)    # Posterior mean pred 1
phi.pred2.pm <- apply(phi.pred2, 1, mean)    # Posterior mean pred 2
phi.pred3.pm <- apply(phi.pred3, 1, mean)    # Posterior mean pred 3
phi.pred3.psd <- apply(phi.pred3, 1, sd)     # Posterior sd pred 3

library(raster)

# Here's the actual Fig. 3-16
op <- par(mfrow = c(1, 2))
# Plot posterior mean of predicted apparent survival (predictions 1)
r1 <- with(willowWarbler, rasterFromXYZ(data.frame(x = cells$lon,
    y = cells$lat, z = phi.pred1.pm)))
plot(r1, col = rampYOR(100), axes = FALSE, box = FALSE, main ="Covariates only",
    zlim = c(0.1, 0.4))
# Plot posterior mean of predicted apparent survival (predictions 2)
r2 <- with(willowWarbler, rasterFromXYZ(data.frame(x = blocks$blockX,
    y = blocks$blockY, z = phi.pred2.pm)))
plot(r2, col = rampYOR(100), axes = FALSE, box = FALSE,
    main ="Residual spatial field only", zlim = c(0.25, 0.3))

# Here's the actual Fig. 3-17
# Plot posterior mean of predicted apparent survival
r1 <- with(willowWarbler, rasterFromXYZ(data.frame(x = cells$lon, y = cells$lat,
    z = phi.pred3.pm)))
plot(r1, col = rampYOR(100), axes = FALSE, box = FALSE)
# Plot unicertainty in this estimate of predicted apparent survival
r2 <- with(willowWarbler, rasterFromXYZ(data.frame(x = cells$lon, y = cells$lat,
    z=phi.pred3.psd)))
plot(r2, col = rampYOR(100), axes = FALSE, box = FALSE, zlim = c(0, 0.08))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
