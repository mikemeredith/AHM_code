#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================
# Code from proofs dated 2020-08-19

if(!requireNamespace("RandomFields"))
  stop("Package 'RandomFields' is not available.")

# Approximate execution time for this code: 6 mins
# Run time with the full number of iterations: 1 hr

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14" # the location of the WinBUGS14.exe file on your machine


# 9.4 Descriptive models of spatial autocorrelation
# =================================================

# 9.4.3 Fitting an intrinsic CAR model to simulated occupancy data
# ----------------------------------------------------------------

# Simulate a detection/nondetection data set with spatial dependence
library(AHMbook)
RNGversion("3.5.3")
str(dat <- simOccSpatial(nsurveys = 3, mean.psi = 0.9, beta = c(2, -2),
    mean.p = 0.4, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
    theta.RF = 10, seeds = c(10, 100)))

# Prepare data
# Compute the neighborhood (2nd order)
library(spdep)
coordgrid <- cbind(dat$xcoord, dat$ycoord)
neigh <- dnearneigh(coordgrid, d1 = 0, d2 = sqrt(2) * 1000 + 1)
winnb <- nb2WB(neigh) # Function to get CAR ingredients for BUGS
str(winnb)

# Bundle data
str(bdata <- list(y = dat$yobs, nsites = dim(dat$y)[1], nrep = dim(dat$y)[2],
    adj = winnb$adj, weights = winnb$weights, num = winnb$num,
    elev = dat$elevationS, forest = dat$forestS, wind = dat$wind))
# List of 9
# $ y      : int [1:2500, 1:3] NA NA 0 NA NA NA NA 0 NA NA ...
# $ nsites : int 2500
# $ nrep   : int 3
# $ adj    : int [1:19404] 2 51 52 1 3 51 52 53 2 4 ...
# $ weights: num [1:19404] 1 1 1 1 1 1 1 1 1 1 ...
# $ num    : int [1:2500] 3 5 5 5 5 5 5 5 5 5 ...
# $ elev   : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
# $ forest : num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
# $ wind   : num [1:2500, 1:3] 1.2703 -0.8881 -0.0708 -0.5315 -1.4702 ...


# Specify model in BUGS language
cat(file = "CAR.occ.txt", "
model {
  # ~~~ code inserted from MS dated 2019-07-02 ~~~~~~~~~~~~~
  # Specify priors
  beta0 <- logit(mean.psi)
  mean.psi ~ dunif(0, 1)
  alpha0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # CAR prior distribution for spatial random effects
  eta[1:nsites] ~ car.normal(adj[], weights[], num[], tau)

  # Parameter expansion for spatial random effects for better mixing
  for (i in 1:nsites){
    S[i] <- Xi * eta[i]
  }

  # Hierarchical half-Cauchy prior for random effect standard deviations
  # from Gelman (2006), also see Chelgren et al. (Ecology 2011a)
  Xi ~ dnorm(0, tau.xi)
  tau.xi <- pow(sig.xi, -2)
  sig.xi <- u # Cauchy scale parameter u estimated from data
  u ~ dunif(0, 3) # Prior for the Cauchy scale parameter
  sd.eta <- abs(Xi)/sqrt(tau)
  tau ~ dgamma(0.5, 0.5) # This is a chisquare rv with 1 df
  v.eta <- pow(sd.eta, 2)

  # Model for occupancy
  for (i in 1:nsites){
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- lpsi.lim[i]
    lpsi.lim[i] <- min(1000, max(-1000, lpsi[i]))
    lpsi[i] <- beta0 + beta[1] * elev[i] + beta[2] * pow(elev[i],2) + S[i]
  }
  # ~~~ code inserted from MS dated 2019-07-02 ~~~~~~~~~~~~~
  # Measurement error model
  for (i in 1:nsites){
    for (j in 1:nrep){
      y[i,j] ~ dbern(mup[i,j])
      mup[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(1000, max(-1000, lp[i,j]))  # 'Stabilize' logit
      lp[i,j] <- alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j]
    }
  }

  # Derived parameters: Total number of occupied sites
  nocc <- sum(z[])
}
")

# Initial values
zst <- apply(dat$yobs, 1, max)
zst[is.na(zst)] <- 1
inits <- function(){ list(z = zst, mean.psi = 0.5, beta = rep(0, 2),
    mean.p = 0.5, alpha = rep(0, 2), eta = rep(0, nrow(coordgrid)))}

# Parameters monitored
params <- c("mean.psi", "beta0", "beta", "mean.p", "alpha0", "alpha",
    "v.eta", "nocc", "psi", "eta", "S")

# MCMC settings
# ni <- 100000  ;  nt <- 50  ;  nb <- 50000  ;  nc <- 3 # 111 min
ni <- 6000    ;    nt <- 1   ;    nb <- 5000    ;    nc <- 3 # ~~~ for testing, 4.5 mins

# Call WinBUGS (ART 111 min) and summarize posteriors
library(R2WinBUGS)
out4 <- bugs(bdata, inits, params, "CAR.occ.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE,
    bugs.directory = bugs.dir)
options(scipen = 10)
print(out4$summary[1:10,-(4:6)], dig = 3)

truth <- c(mean.psi = dat$mean.psi, beta0 = dat$beta0, beta1 = dat$beta[1],
    beta2 = dat$beta[2], mean.p = dat$mean.p, alpha0 = dat$alpha0,
    alpha1 = dat$alpha[1], alpha2 = dat$alpha[2], Nocc = dat$trueNocc)
estimates <- rbind(out4$summary[c('mean.psi', 'beta0', 'beta[1]', 'beta[2]',
    'mean.p', 'alpha0', 'alpha[1]', 'alpha[2]', 'nocc'), c(1:3,7)])
print(cbind(truth, estimates), 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#             truth     mean      sd    2.5%    97.5%
# mean.psi    0.900    0.790  0.0753   0.646    0.941
# beta0       2.197    1.405  0.5581   0.601    2.768
# beta1       2.000    1.931  0.5753   1.191    3.364
# beta2      -2.000   -2.011  0.6301  -3.731   -1.217
# mean.p      0.400    0.399  0.0287   0.346    0.455
# alpha0     -0.405   -0.413  0.1199  -0.635   -0.180
# alpha1     -1.000   -0.894  0.1241  -1.142   -0.654
# alpha2     -1.000   -1.063  0.1214  -1.303   -0.832
# Nocc     1166.000 1107.278 66.2683 983.000 1235.025


# ~~~~ extra code for figure 9.12 ~~~~~~~~~~~~~~~~
library(raster)
op <- par(mfrow = c(2, 2))
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2],
    z = c(dat$field)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Spatial effect (true)")
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2],
    z = out4$mean$S))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Spatial effect (estimate)")
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2],
    z = dat$psi))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Occupancy prob. (true)", zlim = c(0, 1))
title(main = "", line = -2)
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2],
    z = out4$mean$psi))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Occupancy prob. (estimate)", zlim = c(0, 1))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
