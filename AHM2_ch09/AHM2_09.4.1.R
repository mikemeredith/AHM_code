#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================
# Code from proofs dated 2020-07-15

# Approximate execution time for this code: 20 mins
# Run time with the full number of iterations: 3 hrs

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14"
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
library(spdep)
coordgrid <- cbind(bo$x, bo$y)
neigh <- dnearneigh(coordgrid, d1 = 0, d2 = sqrt(2) * 1000 + 1)
winnb <- nb2WB(neigh) # Function to get CAR ingredients for BUGS
str(winnb)
# List of 3
# $ adj : int [1:19404] 2 51 52 1 3 51 ... # Index of neighbours
# $ weights: num [1:19404] 1 1 1 1 1 1 1 11 ... # Weights
# $ num : int [1:2500] 3 5 5 5 5 5 5 5 ... # Size of neighbourhood

# Bundle data
# str(bdata <- list(y = yobs, nsites = dim(y)[1], nrep = dim(y)[2],
    # adj = winnb$adj, weights = winnb$weights, num = winnb$num,
    # elev = as.numeric(elev), forest = as.numeric(forest), wind = wind))
str(bdata <- with(dat, list(y = yobs, nsites = dim(y)[1], nrep = dim(y)[2],
    adj = winnb$adj, weights = winnb$weights, num = winnb$num,
    elev = elevationS, forest = forestS, wind = wind)))
# List of 9
# $ y : int [1:2500, 1:3] NA NA 0 NA NA NA NA 3 NA NA ...
# $ nsites : int 2500
# $ nrep : int 3
# $ adj : int [1:19404] 2 51 52 1 3 51 52 53 2 4 ...
# $ weights: num [1:19404] 1 1 1 1 1 1 1 1 1 1 ...
# $ num : int [1:2500] 3 5 5 5 5 5 5 5 5 5 ...
# $ elev : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
# $ forest : num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
# $ wind : num [1:2500, 1:3] 0.534 1.369 -0.426 0.747 -0.414 ...

# Specify model in BUGS language
cat(file = "CAR.Nmix.txt", "
model {
  # Specify priors
  beta0 <- log(mean.lam)
  mean.lam ~ dunif(0, 20)
  alpha0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }
  # CAR prior distribution for spatial random effects
  eta[1:nsites] ~ car.normal(adj[], weights[], num[], tau)
  v.eta ~ dnorm(0, 0.01)I(0,)
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
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(1000, max(-1000, lp[i,j])) # 'Stabilize' logit
      lp[i,j] <- alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j]
    }
  }
  # Derived parameters: Total population size on grid
  Ntotal <- sum(N[])
}
")

# Initial values
Nst <- apply(dat$yobs, 1, max) # Max observed abundance as inits for N
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2
inits <- function(){ list(N = Nst, mean.lam = 1, beta = rep(0, 2), mean.p = 0.5,
    alpha = rep(0, 2), eta = rep(0, nrow(coordgrid)))}

# Parameters monitored
params <- c("mean.lam", "beta0", "beta", "mean.p", "alpha0", "alpha", "v.eta",
    "Ntotal", "eta", "lam")

# MCMC settings
# ni <- 200000 ; nt <- 100 ; nb <- 100000 ; nc <- 3
ni <- 15000 ; nt <- 1 ; nb <- 10000 ; nc <- 3  # ~~~ for testing, 17 mins
    # nb = 1000 not enough for adaptation to complete

# Call WinBUGS from R (ART longish) and summarize posteriors
library(R2WinBUGS)
# out2 <- bugs(bdata, inits, params, "CAR.Nmix.txt", n.chains = nc, n.thin = nt,
    # n.iter = ni, n.burnin = nb, bugs.directory = bugs.dir,
    # debug = TRUE, working.directory = getwd()) # Close program manually !
out2 <- bugs(bdata, inits, params, "CAR.Nmix.txt", n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, bugs.directory = bugs.dir)  # ~~~~ for testing
print(out2$summary[1:10,], 2)
# mean sd 2.5% 25% 50% 75% 97.5% Rhat n.eff
# mean.lam 3.5569 0.386 2.86 3.289 3.5310 3.799 4.34 1 3000
# beta0 1.2630 0.109 1.05 1.191 1.2620 1.335 1.47 1 3000
# beta[1] 1.9796 0.164 1.67 1.865 1.9770 2.091 2.31 1 2100
# beta[2] -1.9704 0.136 -2.24 -2.061 -1.9690 -1.878 -1.71 1 2100
# mean.p 0.5010 0.026 0.45 0.484 0.5013 0.518 0.55 1 2300
# alpha0 0.0041 0.102 -0.20 -0.065 0.0051 0.071 0.20 1 2300
# alpha[1] -0.9065 0.088 -1.08 -0.965 -0.9063 -0.848 -0.73 1 2300
# alpha[2] -0.9531 0.066 -1.08 -0.996 -0.9522 -0.909 -0.82 1 3000
# v.eta 2.3962 0.404 1.69 2.119 2.3670 2.652 3.26 1 2000
# Ntotal 7235.9227 649.527 6113.90 6771.000 7187.5000 7650.250 8659.17 1 3000

# ~~~~ save for comparison with other models ~~~~~~~
save(out2, file="AHM2_09.4.1_out2.RData")

with(dat, cbind(beta0, beta[1], beta[2], alpha0, alpha[1], alpha[2], Ntotal = sum(N),
    summaxC = sum(apply(y,1,max))))
# beta0 beta1 beta2 alpha0 alpha1 alpha2 Ntotal summaxC
# [1,] 2 2 -2 0 -1 -1 7192 4371

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
# needs out1 object from section
load("AHM2_09.3_out1.RData")
truth <- with(dat, c(alpha0 = alpha0, alpha1 = alpha[1], alpha2 = alpha[2], beta0 = beta0,
    beta1 = beta[1], beta2 = beta[2], Ntotal = sum(N)))
Bayes.est.Nmix0 <- rbind(out1$summary[c('alpha0','alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
Bayes.est.CAR <- rbind(out2$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
print(cbind(truth, Bayes.est.Nmix0, Bayes.est.CAR), 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# truth mean sd mean sd
# alpha0 0 0.31 0.073 0.0041 0.102
# alpha1 -1 -0.90 0.067 -0.9065 0.088
# alpha2 -1 -1.14 0.060 -0.9531 0.066
# beta0 2 1.29 0.052 1.2630 0.109
# beta1 2 2.45 0.102 1.9796 0.164
# beta2 -2 -1.80 0.086 -1.9704 0.136
# Ntotal 7192 5863.50 238.104 7235.9227 649.527

