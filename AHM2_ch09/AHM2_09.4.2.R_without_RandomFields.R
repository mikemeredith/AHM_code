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

# Approximate execution time for this code: 7 mins
# Run time with the full number of iterations: 30 mins

library(AHMbook)
library(raster)
library(jagsUI)

# 9.4 Descriptive models of spatial autocorrelation
# =================================================

# 9.4.2 Fitting a two-dimensional generalized additive model
# ----------------------------------------------------------

# Re-create data set and sample 500 sites (using the same seeds)
library(AHMbook)
RNGversion("3.5.3")
str(dat <- simNmixSpatial(nsurveys = 3, mean.lambda = exp(2),
    beta = c(2, -2), mean.p = 0.5, alpha = c(-1, -1), sample.size = 500,
    variance.RF = 1, theta.RF = 10, seeds = c(10, 100), truncN = 6,
    show.plots = TRUE))

# Scale both sets of coordinates
head(coordgrid <- scale(cbind(dat$xcoord, dat$ycoord))) # All 2500 cells
head(sitelocs <- coordgrid[dat$surveyed.sites,])        # 500 surveyed cells

library(fields)
system.time( knots <- cover.design(R = coordgrid, nd = 125, nruns = 10,
    num.nn = 200, max.loop=20) )                        # takes about 4 min

# Define the Z matrix for the random effects/knot coefficients
knotlocs <- knots$design
omega <- (e2dist(knotlocs, knotlocs)/10)^3
svd.omega <- svd(omega)
sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
Zk <- (e2dist(coordgrid, knotlocs)/10)^3
Zmat <- t(solve(sqrt.omega, t(Zk)))

# Visualize basis vectors (as in Fig. 9.7)
head(Zmat) # Look at first couple of values
library(raster)
devAskNewPage(ask = dev.interactive(orNone=TRUE))  # ~~~ better for testing
op <- par(mfrow = c(3, 3), mar = c(2,2,4,11), "ask")
# for(i in 1:125){
for(i in 1:9){  # ~~~ for testing
  r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2],
      z = Zmat[,i]))
  plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
      main = paste("Basis vector", i), legend = TRUE,
      axis.args = list(cex.axis = 2), legend.width = 2)
  points(knots$design[i,1], knots$design[i,2], col = 'black',
      pch = 16, cex = 1.5)
}
par(op)

# Bundle data
n.knots <- 125 # Fill in the desired spline bases and #knots
y <- dat$yobs # 2500 sites, but not surveyed NA'd out
str(bdata <- list(y = y, nsites = dim(y)[1], sample.size = dat$sample.size,
    nrep = dim(y)[2], elev = dat$elevationS, forest = dat$forestS,
    wind = dat$wind, n.knots = n.knots, Zmat = Zmat))
# List of 9
# $ y          : int [1:2500, 1:3] NA NA 0 NA NA NA NA 3 NA NA ...
# $ nsites     : int 2500
# $ sample.size: num 500
# $ nrep       : int 3
# $ elev       : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
# $ forest     : num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
# $ wind       : num [1:2500, 1:3] 0.534 1.369 -0.426 0.747 -0.414 ...
# $ n.knots    : num 125
# $ Zmat       : num [1:2500, 1:125] 0.00548 0.00545 0.00544 ...

# Specify model in BUGS language
cat(file = "2dSplines.Nmix.txt", "
model {

  # Specify priors: intercepts and slopes of lambda and p
  beta0 <- log(mean.lam)
  mean.lam ~ dunif(0, 30)
  alpha0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }

  # Priors on random effects parameters representing the splines
  for (k in 1:n.knots){
    b[k] ~ dnorm(0, tau.b)
  }
  # Prior on random effects dispersion
  tau.b ~ dgamma(0.1, 0.1)
  sd.b <- pow(tau.b, -2)

  # Model for abundance
  for (i in 1:nsites){
    N[i] ~ dpois(lam[i])
    log(lam[i]) <- beta0 + beta[1] * elev[i] + beta[2] * pow(elev[i],2) + smooth[i]
    smooth[i] <- smooth2[i] - mean(smooth2[])
    smooth2[i] <- inprod(Zmat[i,], b[])
  }

  # Measurement error model
  for (i in 1:nsites){
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      logit(p[i,j]) <- alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j]
    }
  }

  # Derived parameters: Total population size on grid
  Ntotal <- sum(N[])
}
")

# Initial values
Nst <- apply(dat$yobs, 1, max)
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2
inits <- function(){ list(N = Nst, mean.lam = exp(rnorm(1)), beta = rnorm(2),
    mean.p = runif(1), alpha = rnorm(2), b = runif(bdata$n.knots, -0.5, 0.5))}

# Parameters monitored
params <- c("mean.lam", "beta0", "beta", "mean.p", "alpha0", "alpha", "Ntotal",
    "sd.b", "b", "lam", "smooth")

# MCMC settings
# na <- 10000 ; ni <- 40000 ; nt <- 20 ; nb <- 20000 ; nc <- 3
na <- 1000 ; ni <- 4000 ; nt <- 2 ; nb <- 2000 ; nc <- 3  # ~~~ for testing, 3 mins

# Call JAGS (ART 40 min), gauge convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "2dSplines.Nmix.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
summary(out3)
which(out3$summary[,8] > 1.1) ; sum(out3$summary[,8] > 1.1)
options(scipen = 10)
print(out3$summary[1:18, -c(4:6, 8:11)], dig = 3) # not shown

# ~~~ extra code to do the table ~~~~~~~
load("AHM2_09.3_out1_without_RandomFields.RData")
load("AHM2_09.4.1_out2_without_RandomFields.RData")
truth <- with(dat, c(alpha0 = alpha0, alpha1 = alpha[1], alpha2 = alpha[2], beta0 = beta0,
    beta1 = beta[1], beta2 = beta[2], Ntotal = sum(N)))
Bayes.est.Nmix0 <- rbind(out1$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
Bayes.est.CAR <- rbind(out2$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
Bayes.est.spline <- rbind(out3$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
    'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
print(cbind(truth, Bayes.est.Nmix0, Bayes.est.CAR, Bayes.est.spline), 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Truth, nonspatial model, spatial Nmix with CAR, spatial Nmix with GAM
#        truth    mean      sd      mean      sd    mean      sd
# alpha0     0    0.21   0.072    0.045   0.087    0.13   0.083
# alpha1    -1   -0.99   0.056   -0.921   0.072   -0.90   0.063
# alpha2    -1   -1.18   0.053   -1.065   0.058   -1.12   0.053
# beta0      2    1.81   0.046    1.523   0.088    1.39   0.095
# beta1      2    2.06   0.080    1.979   0.159    1.89   0.130
# beta2     -2   -1.92   0.079   -1.870   0.126   -1.73   0.106
# Ntotal  9504 8035.78 311.441 8272.814 511.584 7819.17 346.527

# ~~~~ extra code for figures 9.8, 9.9 and 9.10 ~~~~~~~~
# Fig. 9.8
op <- par(mfrow = c(2, 2), mar = c(3,3,3,4))
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2], z = c(dat$field)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Spatial effect (true)")
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2], z = out3$mean$smooth))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Spatial effect (estimate)")
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2], z = dat$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Density lambda (Truth)", zlim = c(0, 10))
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2], z = out3$mean$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Density lambda (Estimate)", zlim = c(0, 10))
par(op)

# Fig. 9.9, a comparison of cell-wise density
xylim <- c(0, 50)
op <- par(mfrow = c(1, 2))
plot(c(dat$lam), out2$mean$lam, xlab = 'True', ylab = 'Estimated',
    main = 'Spatial Nmix (CAR)', frame = FALSE, pch = 16, cex = 1.5,
    col = rgb(0,0,0,0.3), xlim = xylim, ylim = xylim)
abline(0, 1)
abline(lm(out2$mean$lam ~ c(dat$lam)), col = 'blue', lwd = 3)
plot(c(dat$lam), out3$mean$lam, xlab = 'True', ylab = 'Estimated',
    main = 'Spatial Nmix (2d splines)', frame = FALSE, pch = 16, cex = 1.5,
    col = rgb(0,0,0,0.3), xlim = xylim, ylim = xylim)
abline(0, 1)
abline(lm(out3$mean$lam ~ c(dat$lam)), col = 'blue', lwd = 3)
par(op)

# Fig 9.10 Visualize contributions of piece-wise regressions vectors
op <- par(mfrow = c(3, 4), mar = c(1,1,4,1))
for(i in 1:12){
  piecewise.reg <- Zmat[,i] * out3$mean$b[i]
  r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2],
      z = piecewise.reg))
  plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
      main = paste("Piecewise regression", i) )
  points(knots$design[i,1], knots$design[i,2], col = 'black', pch = 16, cex = 1)
}
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
