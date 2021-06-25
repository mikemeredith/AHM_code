#   Applied hierarchical modeling in ecology - vol.2 - 2021
#   Marc Kéry & J. Andy Royle
#
# Chapter 3 : HIERARCHICAL MODELS OF SURVIVAL
# ===========================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 9 mins
# Run time with the full number of iterations: 1.8 hrs

##### USING NIMBLE INSTEAD OF WinBUGS #####

# Code put together by Mike Meredith

library(AHMbook)
library(nimble)
library(mcmcOutput)

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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 3.4 Spatial hierarchical CJS models
# ===================================

# 3.4.4 Fitting a spatial hierarchical CJS model with spatial autocorrelation only
# --------------------------------------------------------------------------------

# Blocks with a distance between d1 and d2 comprise the neighborhood
library(spdep)
blockcoordgrid <- cbind(as.matrix(willowWarbler$blocks))
neigh <- dnearneigh(blockcoordgrid, d1 = 0, d2 = sqrt(2 * 25^2)+ 0.1)
str(winnb <- nb2WB(neigh)) # Function to get CAR ingredients for BUGS
# List of 3
# $ adj    : int [1:3458] 2 3 4 1 3 5 6 1 2 4 ... # ID of neighbors
# $ weights: num [1:3458] 1 1 1 1 1 1 1 1 1 1 ... # Weights: here, equal
# $ num    : int [1:495] 3 4 6 5 3 6 6 6 5 5 ...  # Number of neighbors

# Frequency distribution of the number of neighbors
table(winnb$num) # Every block is connected to at least two neighbors

# Reformat the m-array for WinBUGS (not necessary for NIMBLE, but keep it anyway)
dim(MARR) # The nyear = 11 dimension (now #2) must come last for WinBUGS
dim(MARRWB <- aperm (MARR, c(3, 1, 2))) # MARR for WinBUGS

# Bundle and summarize data set for NIMBLE
# 'data' = response, everything else = 'constants'
bdata <- list(MARRWB = MARRWB)
str(bconst <- list(R = R, n.site = nsite, n.occ = nyear,
    n.block = nblock, BlockID = willowWarbler$CES$BlockID, adj = winnb$adj,
    weights = winnb$weights, num = winnb$num))
# List of 8
# $ R      : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site : int 193
# $ n.occ  : int 11
# $ n.block: int 495
# $ BlockID: int [1:193] 25 77 204 110 222 283 119 234 295 152 ...
# $ adj    : int [1:3458] 2 3 4 1 3 5 6 1 2 4 ...
# $ weights: num [1:3458] 1 1 1 1 1 1 1 1 1 1 ...
# $ num    : int [1:495] 3 4 6 5 3 6 6 6 5 5 ...

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

cjs8.code <- nimbleCode({
  # Priors and linear models
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      logit(phi[t, s]) <- lphi[t, s] # survival
      logit(p[t, s]) <- lp[t, s]     # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    # eta is spatial effect at the block level
    alpha.lphi.site[s] <- mu.lphi + eta[BlockID[s]]
    alpha.lp.site[s] ~ dnorm(mu.lp, sd = sd.lp.site)

    # backtransform site means
    logit(mean.p.site[s]) <- alpha.lp.site[s]
  }
  for (t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, sd = sd.lphi.time)
    beta.lp.time[t] ~ dnorm(0, sd = sd.lp.time)

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

  # CAR prior distribution for spatial random effects eta
  # NOTE: this is defined on the entire grid of 495 blocks
  # eta[1:n.block] ~ car.normal(adj[], weights[], num[], tau) #####
  eta[1:n.block] ~ dcar_normal(adj[], weights[], num[], tau, zero_mean=1)
  tau ~ dgamma(0.5, 0.0005) # cf. Lunn et al. (2013)
  # curve(1/dgamma(x, 0.5, 0.0005), 0, 10) # howsit look like ?
  veta <- 1/tau
  sdeta <- sqrt(veta)
  # Multinomial likelihood for the m-array data (WinBUGS style)
  # Note 'open index' in pr[t,s,] comes last
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
      }
    }
  }
  # Last column: probability of non-recapture
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      pr[t,s,n.occ] <- 1-sum(pr[t,s,1:(n.occ-1)])
    }
  }
})

# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1),
    eta = rep(0, nblock))}

# Parameters monitored
params <- c("mean.phi", "mean.p", "mu.lphi", "mu.lp",
    "sdeta", "veta","sd.lp.time",  "sd.lp.site", "sd.lphi.time",
    "mean.p.site", "mean.phi.time", "mean.p.time", "eta")


# MCMC settings
ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 2 # ~~~ for testing, 5 mins

system.time(
out <- nimbleMCMC(cjs8.code, data=bdata, constants=bconst,
    inits=inits, monitors=params,
    niter=ni, nburnin=nb, thin=nt, nchains=nc,
    samplesAsCodaMCMC=TRUE) )

mc8 <- mcmcOutput(out)
View(summary(mc8, n.eff=TRUE))
diagPlot(mc8)


# Run multiple NIMBLE instances in parallel with 'foreach'
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''
library(foreach)
library(parallel)
library(doParallel)

detectCores()  # number of cores on your machine
ncore <- 3     # <-- adjust this for your machine (3 used for testing)

cl <- makeCluster(ncore)
registerDoParallel(cl)

ni <- 1000  ;  nt <- 1  ;  nb <- 500  # ~~~ testing, took 3.5 mins
# ni <- 30000  ;  nt <- 30  ;  nb <- 15000  # took 34 mins
# ni <- 10e4 ; nt <- 50 ; nb <- 5e4  # 1.7 hrs

seeds <- 1:ncore

system.time(
res8 <- foreach(x = seeds, .packages="nimble",
  .errorhandling='remove') %dopar% {
  set.seed(x)
  outx <- nimbleMCMC(cjs8.code, data=bdata, constants=bconst, inits=inits, monitors=params,
    niter=ni, nburnin=nb, thin=nt, nchains=1,
    samplesAsCodaMCMC=TRUE)
  outx
} )
stopCluster(cl)

length(res8)

mclist <- coda::mcmc.list(res8)
# Convert to an mcmcOutput object and look at diagnostics
( mco8 <- mcmcOutput(mclist) )
diagPlot(mco8[params[1:6]])
View(summary(mco8))
save(mco8, file="AHM2_03.4.4_nimble_mco8.RData")

# ~~~~~~~ code for figure 3.14 ~~~~~~~~~~~~~~~~~~~~~
library(raster)
phi.block <- array(NA, dim = c(nblock, nrow(mco8)))
sims <- mco8

for (b in 1:nblock){
  phi.block[b,] <- plogis(sims$mu.lphi + sims$eta[,b])
}
post.mean <- apply(phi.block, 1, mean)    # Posterior mean
post.sd <- apply(phi.block, 1, sd)        # Posterior standard deviation
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

op <- par(mfrow = c(1, 2))
# Plot posterior mean of predicted apparent survival
r1 <- rasterFromXYZ(data.frame(x = willowWarbler$blocks$blockX,
    y = willowWarbler$blocks$blockY, z = post.mean))
plot(r1, col = mapPalette(100), axes = FALSE, box = FALSE)
# points(CES[, 1:2], pch = 16, col='black', cex = 1) # Can add CES locations
# Plot uncertainty in this estimate of predicted apparent survival
r2 <- rasterFromXYZ(data.frame(x = willowWarbler$blocks$blockX,
    y = willowWarbler$blocks$blockY, z = post.sd))
plot(r2, col = mapPalette(100), axes = FALSE, box = FALSE, zlim = c(0, 0.08))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
