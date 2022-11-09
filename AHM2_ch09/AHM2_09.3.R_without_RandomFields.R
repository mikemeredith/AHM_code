#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
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

# Approximate execution time for this code: 13 mins

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14" # the location of the WinBUGS14.exe file on your machine

# ~~~ regenerate the data ~~~~~~
RNGversion("3.5.3")
dat <- simNmixSpatial(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
    mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
    theta.RF = 10,seeds = c(10, 100), truncN = 6, show.plots=FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 9.3 Fitting a nonspatial N-mixture model to the simulated data
# ==============================================================

# ~~~ extra WinBUGS code for Nmix fitting ~~~~~~~~~~~~~~~~~~~
# Bundle data
bdata <- with(dat, list(y = yobs, nsites = dim(y)[1], nrep = dim(y)[2],
    elev = elevationS, forest = forestS, wind = wind))
str(bdata)
# List of 6
 # $ y     : int [1:2500, 1:3] NA NA 0 NA NA NA NA 3 NA NA ...
 # $ nsites: int 2500
 # $ nrep  : int 3
 # $ elev  : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
 # $ forest: num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
 # $ wind  : num [1:2500, 1:3] 0.534 1.369 -0.426 0.747 -0.414 ...

# Specify model in BUGS language
cat(file = "Nmix.txt", "
model {

  # Priors
  beta0 <- log(mean.lam)
  mean.lam ~ dunif(0, 20)
  alpha0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }

  # Model for abundance
  for (i in 1:nsites){
   loglam[i] <- beta0 + beta[1] * elev[i] + beta[2] * pow(elev[i],2)
   loglam.lim[i] <- min(1000, max(-1000, loglam[i]))  # Stabilize log
    lam[i] <- exp(loglam.lim[i])
    N[i] ~ dpois(lam[i])
  }

  # Measurement error model
  for (i in 1:nsites){
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(1000, max(-1000, lp[i,j]))  # Stabilize logit
      lp[i,j] <- alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j]
    }
  }

  # Derived parameters: total population size in grid
  Ntotal <- sum(N[])
}
"
)

# Initial values
Nst <- apply(dat$yobs, 1, max)# Max observed abundance as inits for N
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2
inits <- function(){ list(N = Nst, mean.lam = 1, beta = rep(0, 2),
    mean.p = 0.5, alpha = rep(0, 2))}

# Parameters monitored
params <- c("mean.lam", "beta0", "beta", "mean.p", "alpha0", "alpha", "Ntotal", "lam")

# MCMC settings
ni <- 20000    ;    nt <- 10    ;    nb <- 10000    ;    nc <- 3  # 13 mins

# Call WinBUGS from R (ART 18 min) and summarize posteriors
library(R2WinBUGS)
out1 <- bugs(bdata, inits, params, "Nmix.txt", n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
print(out1, dig = 2)
# save result for comparison with other models
save(out1, file="AHM2_09.3_out1_without_RandomFields.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Inference for Bugs model at "Nmix.txt", fit using WinBUGS,
# 3 chains, each with 20000 iterations (first 10000 discarded), n.thin = 10
# n.sims = 3000 iterations saved # ART 18 min
#             mean     sd    2.5%     25%     50%     75%   97.5% Rhat n.eff
# mean.lam    3.62   0.19    3.26    3.50    3.62    3.74    4.00    1  3000
# beta0       1.29   0.05    1.18    1.25    1.29    1.32    1.39    1  3000
# beta[1]     2.45   0.10    2.25    2.38    2.45    2.52    2.65    1  1200
# beta[2]    -1.80   0.09   -1.98   -1.86   -1.80   -1.74   -1.64    1  1100
# mean.p      0.58   0.02    0.54    0.56    0.58    0.59    0.61    1  3000
# alpha0      0.31   0.07    0.17    0.26    0.31    0.36    0.45    1  2800
# alpha[1]   -0.90   0.07   -1.04   -0.95   -0.90   -0.86   -0.78    1  1500
# alpha[2]   -1.14   0.06   -1.26   -1.18   -1.14   -1.10   -1.03    1  2300
# Ntotal   5863.50 238.10 5423.97 5701.00 5857.50 6018.00 6350.02    1  3000
# [....]

with(dat, cbind(beta0, beta, alpha0, alpha, Ntotal = sum(N),
    summaxC = sum(apply(y,1,max)))) # Remember the truth
#      beta0 beta1 beta2 alpha0 alpha1 alpha2 Ntotal summaxC
# [1,]     2     2    -2      0     -1     -1   7192    4371

# ~~~ extra code for figure 9.5 ~~~~~~~~~~~
# Compare maps of true and estimated density (lambda)
library(raster)
op <- par(mfrow = c(1, 2), mar = c(3, 3, 3,3))
# r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = lam))
r <- with(dat, rasterFromXYZ(data.frame(x = xcoord, y = ycoord, z = lam)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Expected abundance lambda (true)", zlim = c(0, 10))
# r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = out1$mean$lam))
r <- with(dat, rasterFromXYZ(data.frame(x = xcoord, y = ycoord, z = out1$mean$lam)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Expected abundance lambda (estimate)", zlim = c(0, 10))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
