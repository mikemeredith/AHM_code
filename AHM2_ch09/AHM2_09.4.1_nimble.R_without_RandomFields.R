#   Applied hierarchical modeling in ecology - Vol. 2
#   Marc Kéry & J. Andy Royle
#
# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================

# The code below is modified to run without the RandomFields package; if
#   RandomFields is available, it will be used by AHMbook, and the results should
#   be the same.
# When RandomFields is not available, the 'fields' package is used instead, and
#   the results will be different.
if(requireNamespace("RandomFields"))
  stop("Package 'RandomFields' IS available; this script is not needed.")

# Approximate execution time for this code: 12 mins
# Run time with the full number of iterations: 4.2 hrs

##### USING NIMBLE INSTEAD OF WinBUGS #####

# Code put together by Mike Meredith

library(AHMbook)
library(spdep)
library(nimble)
library(mcmcOutput)
library(fields)
library(raster)

# ~~~~~~~ uses data from 9.2 ~~~~~~~~
data(BerneseOberland)
bo <- BerneseOberland
RNGversion("3.5.3")
dat <- simNmixSpatial(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
    mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
    theta.RF = 10,seeds = c(10, 100), truncN = 6, show.plots=FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 9.4 Descriptive models of spatial autocorrelation
# =================================================

# 9.4.1 Fitting an intrinsic conditional autoregressive model
# -----------------------------------------------------------

# Compute the queen's neighborhood

coordgrid <- cbind(bo$x, bo$y)
neigh <- dnearneigh(coordgrid, d1 = 0, d2 = sqrt(2) * 1000 + 1)
winnb <- nb2WB(neigh)      # Function to get CAR ingredients for BUGS
str(winnb)
# List of 3
# $ adj    : int [1:19404] 2 51 52 1 3 51 ...   # Index of neighbours
# $ weights: num [1:19404] 1 1 1 1 1 1 1 11 ... # Weights
# $ num    : int [1:2500] 3 5 5 5 5 5 5 5 ...   # Size of neighbourhood

# Bundle data
bdata <- list(y = dat$yobs)
str(bconst <- with(dat, list(nsites = dim(y)[1], nrep = dim(y)[2],
    adj = winnb$adj, weights = winnb$weights, num = winnb$num,
    elev = elevationS, forest = forestS, wind = wind)))
# List of 8
# $ nsites : int 2500
# $ nrep   : int 3
# $ adj    : int [1:19404] 2 51 52 1 3 51 52 53 2 4 ...
# $ weights: num [1:19404] 1 1 1 1 1 1 1 1 1 1 ...
# $ num    : int [1:2500] 3 5 5 5 5 5 5 5 5 5 ...
# $ elev   : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
# $ forest : num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
# $ wind   : num [1:2500, 1:3] 0.534 1.369 -0.426 0.747 -0.414 ...

# Specify model in NIMBLE dialect of BUGS language
# ''''''''''''''''''''''''''''''''''''''''''''''''
# Necessary changes:
# 1. Replace I(0,) with T(..., 0, )
# 2. dcar_normal needs zero_mean=1 argument to mimic WinBUGS version
# Optional improvements:
# 1. Use logit(.) <-
# 2. dbeta for probabilities instead of dunif

CAR.Nmix.code <- nimbleCode({

  # Specify priors
  beta0 <- log(mean.lam)
  mean.lam ~ dunif(0, 20)
  alpha0 <- logit(mean.p)
  mean.p ~ dbeta(1, 1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }

  # CAR prior distribution for spatial random effects
  # eta[1:nsites] ~ car.normal(adj[], weights[], num[], tau)
  eta[1:nsites] ~ dcar_normal(adj[], weights[], num[], tau, zero_mean=1)
  v.eta ~ T(dnorm(0, 0.01), 0, ) # I(0,)
  tau <- 1/v.eta

  # Model for abundance
  for (i in 1:nsites){
    loglam[i] <- beta0 + beta[1] * elev[i] + beta[2] * pow(elev[i],2) + eta[i]
    loglam.lim[i] <- min(1000, max(-1000, loglam[i])) # 'Stabilize' log
    lam[i] <- exp(loglam.lim[i])
    N[i] ~ dpois(lam[i])
  }

  # Measurement error model
  for (i in 1:nsites){
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      logit(p[i,j]) <- lp.lim[i,j]
      lp.lim[i,j] <- min(1000, max(-1000, lp[i,j])) # 'Stabilize' logit
      lp[i,j] <- alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j]
    }
  }

  # Derived parameters: Total population size on grid
  Ntotal <- sum(N[])
} )

# Initial values
Nst <- apply(dat$yobs, 1, max) # Max observed abundance as inits for N
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2
# Provide initial values too for the missing values in 'y'
# (not essential, but avoids a slew of ugly warnings).
yst <- dat$yobs
yst[is.na(dat$yobs)] <- 1
yst[!is.na(dat$yobs)] <- NA
inits <- function(){ list(y = yst, N = Nst, mean.lam = 1, beta = rep(0, 2), mean.p = 0.5,
    alpha = rep(0, 2), eta = rep(0, nrow(coordgrid)))}

# Parameters monitored
params <- c("mean.lam", "beta0", "beta", "mean.p", "alpha0", "alpha", "v.eta",
    "Ntotal", "eta", "lam")

# MCMC settings
ni <- 1500 ; nt <- 1 ; nb <- 1000 ; nc <- 2  # ~~~ for testing, 9.5 mins

system.time(
out <- nimbleMCMC(CAR.Nmix.code, data=bdata, constants=bconst,
    inits=inits, monitors=params,
    niter=ni, nburnin=nb, thin=nt, nchains=nc,
    samplesAsCodaMCMC=TRUE) )
mc <- mcmcOutput(out)
View(summary(mc, n.eff=TRUE))
diagPlot(mc)

# Use 'foreach' to do a long run in parallel
# ''''''''''''''''''''''''''''''''''''''''''
library(foreach)
library(parallel)
library(doParallel)

detectCores()  # number of cores on your machine
ncore <- 3     # <-- adjust this for your machine (3 used for testing)

cl <- makeCluster(ncore)
registerDoParallel(cl)

# ni <- 2e5 ; nt <- 100 ; nb <- 1e5  # 4.2 hrs
ni <- 1500 ; nt <- 1 ; nb <- 1000  # ~~~ for testing, 4.5 mins

seeds <- 1:ncore

system.time(
res2 <- foreach(x = seeds, .packages="nimble",
      .errorhandling='remove', .inorder=FALSE) %dopar% {
  set.seed(x)
  nimbleMCMC(CAR.Nmix.code, data=bdata, constants=bconst,
      inits=inits, monitors=params,
      niter=ni, nburnin=nb, thin=nt, nchains=1,
      samplesAsCodaMCMC=TRUE)
} )
stopCluster(cl)

# How many chains completed successfully?
length(res2)

# Convert to an mcmcOutput object and look at diagnostics
mclist <- coda::mcmc.list(res2)
( mco2 <- mcmcOutput(mclist) )
diagPlot(mco2, c(1:7, 11, 12)) ### FIXME
View(summary(mco2))

# ~~~~ save for comparison with other models ~~~~~~~
out2 <- sumryList(mco2)
out2$summary <- summary(mco2)
save(out2, file="AHM2_09.4.1_out2.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

with(dat, cbind(beta0, beta[1], beta[2], alpha0, alpha[1], alpha[2],
    Ntotal = sum(N), summaxC = sum(apply(y,1,max))))
#      beta0 beta1 beta2 alpha0 alpha1 alpha2 Ntotal summaxC
# [1,]     2     2    -2      0     -1     -1   9504    5922

# ~~~ extra code for figure 9.6 ~~~~~

# Compute average detection probability for each cell
phat <- array(NA, dim = c(2500, 3))
for(j in 1:3){
  phat[,j] <- plogis(out2$mean$alpha0 + out2$mean$alpha[1] * dat$forestS +
      out2$mean$alpha[2] * dat$wind[,j])
}
pmean <- apply(phat, 1, mean)

op <- par(mfrow = c(3, 2), mar = c(3,3,3,4))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = dat$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Density lambda (true)", zlim = c(0, 10))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = out2$mean$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Density lambda (estimate)", zlim = c(0, 10))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = c(dat$field)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Spatial effect eta (true)")
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = out2$mean$eta))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Spatial effect eta (estimate)")
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = out2$sd$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Estimation uncertainty (lambda)", zlim = c(0, 10))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = pmean))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Average detection probability)", zlim = c(0, 1))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ extra code for table of results ~~~~
# needs out1 object from section AHM2 9.3
load("AHM2_09.3_out1.RData")
truth <- with(dat, c(alpha0 = alpha0, alpha1 = alpha[1], alpha2 = alpha[2],
    beta0 = beta0, beta1 = beta[1], beta2 = beta[2], Ntotal = sum(N)))
Bayes.est.Nmix0 <- rbind(out1$summary[c('alpha0','alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
Bayes.est.CAR <- rbind(out2$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
print(cbind(truth, Bayes.est.Nmix0, Bayes.est.CAR), digits=2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        truth    mean      sd    mean      sd
# alpha0     0    0.21   0.072    0.10   0.083
# alpha1    -1   -0.99   0.056   -0.88   0.074
# alpha2    -1   -1.18   0.053   -1.09   0.055
# beta0      2    1.81   0.046    1.19   0.214
# beta1      2    2.06   0.080    1.49   0.195
# beta2     -2   -1.92   0.079   -1.48   0.172
# Ntotal  9504 8035.78 311.441 7774.10 612.859
