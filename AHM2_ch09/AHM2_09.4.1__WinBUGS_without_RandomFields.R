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

# Approximate execution time for this code: 20 mins
# Run time with the full number of iterations: 3 hrs

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14" # the location of the WinBUGS14.exe file on your machine
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
winnb <- nb2WB(neigh)      # Function to get CAR ingredients for BUGS
str(winnb)
# List of 3
# $ adj    : int [1:19404] 2 51 52 1 3 51 ...   # Index of neighbours
# $ weights: num [1:19404] 1 1 1 1 1 1 1 11 ... # Weights
# $ num    : int [1:2500] 3 5 5 5 5 5 5 5 ...   # Size of neighbourhood

# Bundle data
# str(bdata <- list(y = yobs, nsites = dim(y)[1], nrep = dim(y)[2],
    # adj = winnb$adj, weights = winnb$weights, num = winnb$num,
    # elev = as.numeric(elev), forest = as.numeric(forest), wind = wind))
str(bdata <- with(dat, list(y = yobs, nsites = dim(y)[1], nrep = dim(y)[2],
    adj = winnb$adj, weights = winnb$weights, num = winnb$num,
    elev = elevationS, forest = forestS, wind = wind)))
# List of 9
# $ y      : int [1:2500, 1:3] NA NA 0 NA NA NA NA 3 NA NA ...
# $ nsites : int 2500
# $ nrep   : int 3
# $ adj    : int [1:19404] 2 51 52 1 3 51 52 53 2 4 ...
# $ weights: num [1:19404] 1 1 1 1 1 1 1 1 1 1 ...
# $ num    : int [1:2500] 3 5 5 5 5 5 5 5 5 5 ...
# $ elev   : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
# $ forest : num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
# $ wind   : num [1:2500, 1:3] 0.534 1.369 -0.426 0.747 -0.414 ...

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

# Call WinBUGS from R (ART longish) and summarize posteriors
library(R2WinBUGS)
# out2 <- bugs(bdata, inits, params, "CAR.Nmix.txt", n.chains = nc, n.thin = nt,
    # n.iter = ni, n.burnin = nb, bugs.directory = bugs.dir,
    # debug = TRUE, working.directory = getwd()) # Close program manually !
out2 <- bugs(bdata, inits, params, "CAR.Nmix.txt", n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, bugs.directory = bugs.dir)  # ~~~~ for testing
print(out2$summary[1:10,], 2)
#              mean      sd    2.5%      25%      50%     75%   97.5% Rhat n.eff
# mean.lam    4.604   0.407    3.86    4.327    4.584    4.86    5.50  1.0    58
# beta0       1.523   0.088    1.35    1.465    1.522    1.58    1.71  1.0    60
# beta[1]     1.979   0.159    1.66    1.869    1.976    2.08    2.29  1.4     9
# beta[2]    -1.870   0.126   -2.13   -1.956   -1.863   -1.78   -1.64  1.1    19
# mean.p      0.511   0.022    0.47    0.496    0.511    0.53    0.55  1.0   190
# alpha0      0.045   0.087   -0.12   -0.016    0.045    0.11    0.21  1.0   200
# alpha[1]   -0.921   0.072   -1.05   -0.972   -0.921   -0.87   -0.78  1.0    95
# alpha[2]   -1.065   0.058   -1.18   -1.104   -1.065   -1.02   -0.95  1.0   340
# v.eta       1.827   0.308    1.32    1.602    1.795    2.02    2.51  1.0    89
# Ntotal   8272.814 511.584 7335.00 7920.750 8252.000 8599.00 9338.00  1.0   120

# ~~~~ save for comparison with other models ~~~~~~~
save(out2, file="AHM2_09.4.1_out2_without_RandomFields.RData")
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
# needs out1 object from section
load("AHM2_09.3_out1_without_RandomFields.RData")
truth <- with(dat, c(alpha0 = alpha0, alpha1 = alpha[1], alpha2 = alpha[2],
    beta0 = beta0, beta1 = beta[1], beta2 = beta[2], Ntotal = sum(N)))
Bayes.est.Nmix0 <- rbind(out1$summary[c('alpha0','alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
Bayes.est.CAR <- rbind(out2$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
print(cbind(truth, Bayes.est.Nmix0, Bayes.est.CAR), 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        truth    mean      sd     mean      sd
# alpha0     0    0.21   0.072    0.045   0.087
# alpha1    -1   -0.99   0.056   -0.921   0.072
# alpha2    -1   -1.18   0.053   -1.065   0.058
# beta0      2    1.81   0.046    1.523   0.088
# beta1      2    2.06   0.080    1.979   0.159
# beta2     -2   -1.92   0.079   -1.870   0.126
# Ntotal  9504 8035.78 311.441 8272.814 511.584
