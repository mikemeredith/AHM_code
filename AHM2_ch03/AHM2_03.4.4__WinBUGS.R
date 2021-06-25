#   Applied hierarchical modeling in ecology - vol.2 - 2021
#   Marc Kéry & J. Andy Royle
#
# Chapter 3 : HIERARCHICAL MODELS OF SURVIVAL
# ===========================================
# Code from proofs dated 2020-08-18

cat("Approximate run time for this script: 55 mins \n")
# Run time with the full number of iterations: 59 hrs

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14" # the location of the WinBUGS14.exe file on your machine

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

# Reformat the m-array for WinBUGS
dim(MARR) # The nyear = 11 dimension (now #2) must come last for WinBUGS
dim(MARRWB <- aperm (MARR, c(3, 1, 2))) # MARR for WinBUGS

# Bundle and summarize data set for WinBUGS
str(bdata <- list(MARRWB = MARRWB, R = R, n.site = nsite, n.occ = nyear,
    n.block = nblock, BlockID = willowWarbler$CES$BlockID, adj = winnb$adj,
    weights = winnb$weights, num = winnb$num))
# List of 9
# $ MARRWB : num [1:193, 1:10, 1:11] 1 3 0 2 0 0 0 0 0 0 ...
# $ R      : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site : int 193
# $ n.occ  : int 11
# $ n.block: int 495
# $ BlockID: int [1:193] 25 77 204 110 222 283 119 234 295 152 ...
# $ adj    : int [1:3458] 2 3 4 1 3 5 6 1 2 4 ...
# $ weights: num [1:3458] 1 1 1 1 1 1 1 1 1 1 ...
# $ num    : int [1:495] 3 4 6 5 3 6 6 6 5 5 ...

# Specify model in BUGS language
cat(file = "cjs8.txt","
model {

  # Priors and linear models
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      phi[t, s] <- 1 / (1 + exp(-lphi[t, s])) # survival
      p[t, s] <- 1 / (1 + exp(-lp[t, s])) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    # eta is spatial effect at the block level
    alpha.lphi.site[s] <- mu.lphi + eta[BlockID[s]]
    alpha.lp.site[s] ~ dnorm(mu.lp, tau.lp.site) I(-12, 12)

    # backtransform site means
    mean.p.site[s] <- 1 / (1 + exp(-alpha.lp.site[s]))
  }
  for (t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, tau.lphi.time) I(-12, 12)
    beta.lp.time[t] ~ dnorm(0, tau.lp.time) I(-12, 12)

    # backtransform time means (see Errata dated 2021-06-23)
    mean.phi.time[t] <- 1 / (1 + exp(-(mu.lphi + beta.lphi.time[t])))
    mean.p.time[t] <- 1 / (1 + exp(-(mu.lp + beta.lp.time[t])))
  }

  # Hyperpriors for hyperparams
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  tau.lp.site <- pow(sd.lp.site, -2)
  sd.lp.site ~ dunif(0, 2)
  tau.lphi.time <- pow(sd.lphi.time, -2)
  sd.lphi.time ~ dunif(0, 1)
  tau.lp.time <- pow(sd.lp.time, -2)
  sd.lp.time ~ dunif(0, 1)

  # CAR prior distribution for spatial random effects eta
  # NOTE: this is defined on the entire grid of 495 blocks
  eta[1:n.block] ~ car.normal(adj[], weights[], num[], tau)
  tau ~ dgamma(0.5, 0.0005) # cf. Lunn et al. (2013)
  # curve(1/dgamma(x, 0.5, 0.0005), 0, 10) # howsit look like ?
  veta <- 1/tau
  sdeta <- sqrt(veta)
  # Multinomial likelihood for the m-array data (WinBUGS style)
  # Note 'open index' in pr[t,s,] comes last
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      MARRWB[s, t,1:n.occ] ~ dmulti(pr[t,s, ], R[t,s])
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
}
")

# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1),
    eta = rep(0, nblock))}

# Parameters monitored
params <- c("mean.phi", "mean.p", "mu.lphi", "mu.lp", "sd.lp.site", "sd.lphi.time",
    "sd.lp.time", "mean.p.site", "mean.phi.time", "mean.p.time", "veta", "sdeta", "eta")

# ~~~~ alternative code for running WinBUGS ~~~~~
# ~~~~ in parallel is given below ~~~~~~~~~~~~~~~

# MCMC settings
# ni <- 100000 ; nt <- 50 ; nb <- 50000 ; nc <- 3 # 52 hours
ni <- 1000 ; nt <- 5 ; nb <- 500 ; nc <- 3 # ~~~ for testing, 40 mins

# Call WinBUGS from R (ART 52 h!) and summarize posteriors
# bugs.dir must be set to WinBUGS location, e.g., "c:/WinBUGS14/"
out8 <- bugs(bdata, inits, params, "cjs8.txt", n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
print(out8$summary[c(1:7, 221,222),c(1:3,5,7:9)], 3)
#                mean     sd      2.5%     50%  97.5% Rhat n.eff
# mean.phi      0.287 0.0187  0.251397  0.2874  0.324 1.00   730
# mean.p        0.381 0.0356  0.310497  0.3797  0.450 1.01   380
# mu.lphi      -0.910 0.0918 -1.091025 -0.9079 -0.736 1.00   730
# mu.lp        -0.490 0.1519 -0.797902 -0.4907 -0.201 1.01   390
# sd.lp.site    1.129 0.1302  0.887167  1.1220  1.407 1.00  3000
# sd.lphi.time  0.220 0.0834  0.101800  0.2063  0.423 1.00  3000
# sd.lp.time    0.148 0.1117  0.006906  0.1266  0.424 1.00  3000
# veta          0.060 0.0765  0.000452  0.0327  0.262 1.08    38
# sdeta         0.204 0.1361  0.021259  0.1809  0.512 1.08    38

# ~~~~ run multiple WinBUGS instances in parallel with 'foreach' ~~~~~~~~
library(foreach)
library(parallel)
library(doParallel)

detectCores()  # number of cores on your machine
ncore <- 3     # <-- adjust this for your machine (3 used for testing)

cl <- makeCluster(ncore)
registerDoParallel(cl)

ni <- 1000  ;  nt <- 1  ;  nb <- 500  ;  nc <- 1  # ~~~ testing, took 15 mins
# ni <- 30000  ;  nt <- 30  ;  nb <- 15000  ;  nc <- 1  # took 9.25 hrs
seeds <- 1:ncore

# The call to foreach will open ncore WinBUGS windows. These will close on successful completion, but if one throws an error you will need to close it manually. The results from trouble-free instances will be retained.

system.time(
res3 <- foreach(x = seeds, .combine=rbind, .packages="R2WinBUGS",
  .errorhandling='remove') %dopar% {
  set.seed(x)
  outx <- bugs(bdata, inits, params, "cjs8.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
  outx$sims.matrix
} )
stopCluster(cl)

# How many chains completed successfully?
( itersPerChain <- (ni - nb)/nt )
( nChains <- nrow(res3)/itersPerChain )

# Convert to an mcmcOutput object and look at diagnostics
library(mcmcOutput)
( mco8 <- mcmcOutput(res3, nChains = nChains) )
diagPlot(mco8, c(1:7, 11, 12))
View(summary(mco8))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~ code for figure 3.14 ~~~~~~~~~~~~~~~~~~~~~
library(raster)
phi.block <- array(NA, dim = c(nblock, out8$n.sims))
sims <- out8$sims.list        # Grab the simulations first
# ~~~ after the foreach run do this ~~~~~~~~~~~~
phi.block <- array(NA, dim = c(nblock, nrow(mco8)))
sims <- mco8
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
