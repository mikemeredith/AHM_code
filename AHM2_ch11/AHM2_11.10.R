# Applied hierarchical modeling in ecology - vol.2 - 2021
# Marc Kéry & J. Andy Royle
#
# Chapter 11 : SPATIALLY EXPLICIT DISTANCE SAMPLING ALONG TRANSECTS
# =================================================================
# Code from proofs dated 2020-08-19

print("Approximate execution time for this code: 7 mins")

library(AHMbook)
library(jagsUI)

# 11.10 Fully Bayesian density surface modeling (DSM)
# ===================================================

# 11.10.1 Data simulation
# -----------------------

library(raster)
library(AHMbook)
data(wigglyLine)

RNGversion("3.5.0")
set.seed(123, kind = "Mersenne-Twister" )

plot(wigglyLine[,1], wigglyLine[,2], type = "l", pch = " ", lwd = 2,
    cex.axis = 1.5, cex = 2, cex.lab = 1.5, asp = 1, frame = FALSE) # not shown
points <- SpatialPoints(wigglyLine)
sLine <- Line(points)
regpoints <- spsample(sLine, 100, type = "regular")
points(regpoints, col = "black", pch = 20, lwd = 2)

# Simulation settings
set.seed(2027, kind = "Mersenne-Twister")
tmp <- simDSM(X = regpoints@coords, Ntotal = 400, sigma = 0.65, beta1 = 1.0,
    nsurveys = 2, xlim = c(-0.5, 3.5), ylim = c(-0.5, 4.5)) # Produces figure 11.15
str(tmp)
Habitat <- tmp$Habitat
Habgrid <- tmp$Habgrid
nind <- tmp$nind
pixel <- tmp$pixel
N <- tmp$N

# ~~~ extra code for figure 11.14 ~~~~~~~~~~~~~
library(raster)
op <- par(mar = c(3,3,3,6))
image(r <- rasterFromXYZ(cbind(Habgrid, Habitat)), col = topo.colors(10))
image_scale(Habitat, col = topo.colors(10))
lines(wigglyLine, col = "black", pch = 20, lwd = 3)
points(tmp$U, pch = 16)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Harvest data objects
Habitat <- (Habitat - mean(Habitat))/sqrt(var(Habitat))
nPix <- length(Habitat)

# Create a vector of pixel counts and pad it with zeros
yg <- matrix(NA, nrow = nPix, ncol = 2)
yg[,1] <- tabulate(pixel[,1], nbins = nPix)
yg[,2] <- tabulate(pixel[,2], nbins = nPix)

# Show the distribution of observed pixel counts (Fig. 11.16)
plot(rasterFromXYZ(cbind(Habgrid, yg[,1]+yg[,2])))

# Tabulate the pixel frequencies
table(yg)
# yg
#   0   1   2   3   4   5   6   7   8
# 687 207  72  18   8   1   1   3   3

# Compute MINIMUM distance between each pixel and transect
dist <- e2dist(Habgrid, tmp$X)
mind <- apply(dist, 1 ,min)

# Create a nPix x K matrix of pixel counts.
yg <- matrix(NA, nrow = nPix, ncol = 2)
yg[,1] <- tabulate(pixel[,1], nbins = nPix)
yg[,2] <- tabulate(pixel[,2], nbins = nPix)


# 11.10.2 The DSM model (no code)

# 11.10.3 Fitting the models
# --------------------------

# Fit Model 0
# Bundle and summarize the data for BUGS
str(bdata <- list(y = yg, nsites = nPix, dist = mind ))
# List of 3
# $ y     : int [1:500, 1:2] 0 0 0 0 1 1 0 0 0 0 ...
# $ nsites: int 500
# $ dist  : num [1:500] 1.256 1.108 0.969 0.849 0.754 ...

# Write BUGS model
cat(file = "M0.txt","
model {

  # Prior distributions
  sigma ~ dunif(0,10)
  beta0 ~ dnorm(0,0.01)
  alpha0 ~ dnorm(0,0.01)
  logit(p0) <- alpha0

  # Define the likelihood
  for (g in 1:nsites){
    # Detection probability model
    lam[g] <- exp(beta0)
    log(p[g]) <- -(1/(2*sigma*sigma))*dist[g]*dist[g]
    for(k in 1:2){
      y[g,k] ~ dbinom(p[g], N[g])
    }
    N[g] ~ dpois(lam[g])
  }

  # Derived parameters: abundance and GoF calculations
  Ntot <- sum(N[])
  for(k in 1:2){
    for(g in 1:nsites){
      Err[g,k] <- ( sqrt(y[g,k]) - sqrt(lam[g]*p[g]) )^2
      ynew[g,k] ~ dbinom(p[g], N[g])
      Errnew[g,k] <- ( sqrt(ynew[g,k]) - sqrt(lam[g]*p[g]) )^2
    }
  }
  Fit <- sum(Err[,])
  Fitnew <- sum(Errnew[,])
}
")

# MCMC settings
na <- 1000 ; ni <- 5000 ; nb <- 2000 ; nt <- 3 ; nc <- 3

# Inits
Nst1 <- apply(yg,1,max)+1
inits <- function(){ list ( sigma = runif(1, 0.2, 0.7), alpha0 = runif(1, -5, -2),
    beta0 = 0, N = Nst1 ) }

# Parameters to monitor
params1 <- c("sigma", "Ntot", "beta0", "Fit", "Fitnew")

# Run JAGS (ART < 1min), check convergence and summarize
out0 <- jags (bdata, inits, params1, "M0.txt", n.thin = nt,
    n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE )
# par(mfrow = c(2,3))  # ~~~ replaced with 'layout' argument
traceplot(out0, layout=c(2,3))
print(out0, 3)
#           mean     sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# sigma    0.645  0.031   0.589   0.644   0.711    FALSE 1.000 1.003   804
# Ntot   424.126 20.488 386.000 423.000 467.000    FALSE 1.000 1.006   320
# beta0   -0.166  0.068  -0.303  -0.165  -0.036    FALSE 0.994 1.006   321
# Fit    372.992 10.808 353.454 372.269 395.681    FALSE 1.000 1.003   633
# Fitnew 369.701 12.271 346.482 369.315 393.956    FALSE 1.000 1.003   709

beta0 <- mean(out0$sims.list$beta0)
lamhat0 <- exp(beta0)
(RMSE0 <- sum <- sqrt(mean( (N - lamhat0)^2 ) ))
# [1] 1.16 1973
( pval0 <- mean( out0$sims.list$Fit < out0$sims.list$Fitnew ) )
# [1] 0.2516667


# Fit Model 1 (habitat)
# Bundle and summarize data set
str(bdata <- list(y = yg, nsites = nPix, dist = mind, Habitat = Habitat))

# Specify model in BUGS language
cat(file = "M1.txt","
model {

  # Prior distributions
  sigma ~ dunif(0,10)
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  alpha0 ~ dnorm(0,0.01)
  logit(p0) <- alpha0

  # Define the likelihood
  for (g in 1:nsites){
    # Linear model for site-level effects: add covariates
    lam[g] <- exp(beta0 + beta1*Habitat[g] )
    log(p[g]) <- -(1/(2*sigma*sigma))*dist[g]*dist[g]
    for(k in 1:2){
      y[g,k] ~ dbinom(p[g], N[g])
    }
    N[g] ~ dpois(lam[g])
  }

  # Derived parameters: abundance and GoF calculations
  Ntot <- sum(N[])
  for(k in 1:2){
    for(g in 1:nsites){
      Err[g,k] <- ( sqrt(y[g,k]) - sqrt(lam[g]*p[g]) )^2
      ynew[g,k] ~ dbinom(p[g], N[g])
      Errnew[g,k] <- ( sqrt(ynew[g,k]) - sqrt(lam[g]*p[g]) )^2
    }
  }
  Fit <- sum(Err[,])
  Fitnew <- sum(Errnew[,])
}
")

# MCMC settings
na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nt <- 3 ; nc <- 3

# Inits
Nst1 <- apply(yg, 1, max) + 1
inits <- function(){ list (sigma = runif(1, 0.2, 0.7),
    alpha0 = runif(1, -5, -2), beta0 = 0, beta1 = 1, N = Nst1 ) }

# Parameters to monitor
params1 <- c("sigma", "Ntot", "beta0", "beta1", "Fit", "Fitnew" )

# Run JAGS (ART < 1min), check convergence and summarize
out1 <- jags (bdata, inits, params1, "M1.txt", n.thin = nt, n.chains = nc,
    n.burnin = nb, n.iter = ni, n.adapt = na, parallel = TRUE )
# par(mfrow = c(3,3))  # ~~~ no longer necessary
traceplot(out1)
print(out1, 3)
#           mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# sigma    0.667  0.033   0.607   0.665   0.737    FALSE 1 1.001  2270
# Ntot   402.762 18.278 370.000 402.000 441.000    FALSE 1 1.002   915
# beta0   -0.495  0.088  -0.663  -0.495  -0.321    FALSE 1 1.004   540
# beta1    0.736  0.059   0.620   0.737   0.850    FALSE 1 1.000  9999
# Fit    279.227  9.482 262.380 278.649 299.788    FALSE 1 1.004   543
# Fitnew 277.626 11.500 256.343 277.188 301.214    FALSE 1 1.003   725

log(tmp$Ntotal / nrow(Habgrid))
# [1] -0.2231436

# Compute the RMSE and p-value
beta0 <- mean(out1$sims.list$beta0)
beta1 <- mean(out1$sims.list$beta1)
lamhat1 <- exp( beta0 + beta1*Habitat )

(RMSE1 <- sqrt(mean( (tmp$N - lamhat1)^2 )) )
# [1] 1.002903

(pval1 <- mean( out1$sims.list$Fit < out1$sims.list$Fitnew ) )
# [1] 0.3846385

# Set some number of knots and plot
nk <- 86    # Number of knots .... and 86 is a nice number
set.seed(1)
(knotid <- sort(sample(1:nPix, nk)))
head(knotlocs <- Habgrid[knotid,])
knotlocs <- as.matrix(knotlocs)

library(fields)
# Define knot locations using cover.design (do once)
x86 <- cover.design(Habgrid, nd = nk, nruns = 10, num.nn = 200,
    max.loop = 40)                    # runtime < 1 min
# ~~~ use the new knotlocs for subsequent code ~~~~
knotlocs <- x86$design

# ~~~ extra code for figure 11.17 ~~~~~~~~~~~~~~~~~~~~~~
op <- par(mar = c(3,3,3,6))
image(r <- rasterFromXYZ(cbind(Habgrid, Habitat)), col = topo.colors(10))
image_scale(Habitat, col = topo.colors(10))
lines(wigglyLine, col = "black", pch = 20, lwd = 3)
points(tmp$Ucap, pch = 16) # ACs of individuals detected
# Add lines from ACs to nearest point on transect
dd <- e2dist(wigglyLine, tmp$Ucap)
closest <- wigglyLine[apply(dd, 2, which.min), ]
segments(x0=tmp$Ucap[,1], y0=tmp$Ucap[,2], x1=closest[,1], y1=closest[,2])
# Add node locations
points(knotlocs, pch = 3, col = "black", lwd=3)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define the Z matrix for the random effects
Z.k <- (e2dist(Habgrid, knotlocs)/10 )^3
omega.all <- (dist(knotlocs)/10 )^3
svd.omega.all <- svd(omega.all)
sqrt.omega.all <- t(svd.omega.all$v %*%
    (t(svd.omega.all$u)*sqrt(svd.omega.all$d)))
Z.mat <- t(solve(sqrt.omega.all, t(Z.k)))

# Browse through those basis functions (not shown)
for(i in 1:nk){
  plot(rasterFromXYZ(cbind(Habgrid,Z.mat[,i])))
  # browser()           # unhash if you want to have a better look at each
}

# Bundle and summarize data set
str(bdata <- list(y = yg, nsites = nPix, n.knots = nk, Z.mat = Z.mat,
    dist = mind ))
# List of 5
# $ y      : int [1:500, 1:2] 0 0 0 0 1 1 0 0 0 0 ...
# $ nsites : int 500
# $ n.knots: num 86
# $ Z.mat  : num [1:500, 1:86] 0.0911 0.0852 0.0812 0.0791 0.0785 ...
# $ dist   : num [1:500] 1.256 1.108 0.969 0.849 0.754 ...

# Specify model in BUGS language.
cat(file = "DSM1.txt","
model {

  # Prior distributions
  sigma ~ dunif(0,10)
  beta0 ~ dnorm(0,0.01)
  alpha0 ~ dnorm(0,0.01)
  logit(p0) <- alpha0

  # Prior on random effects
  tau.b ~ dgamma(0.01,0.01)
  sigma.b <- sqrt(1/tau.b)
  for (k in 1:n.knots){
    b[k] ~ dnorm(0, tau.b)
  }

  # Define the likelihood
  for (g in 1:nsites){
    # Detection probability model
    log(p[g]) <- -(1/(2*sigma*sigma))*dist[g]*dist[g]
    # Observation model
    for(k in 1:2){
      y[g,k] ~ dbinom(p[g], N[g])
    }
    # Abundance model including GAM
    N[g] ~ dpois(lam[g])
    lam[g] <- exp(beta0 + smooth[g] )
    smooth[g] <- inprod(Z.mat[g,], b[])
  }

  # Derived parameters: abundance and fit statistics
  Ntot <- sum(N[])
  for(k in 1:2){
    for(g in 1:nsites){
      Err[g,k] <- ( sqrt(y[g,k]) - sqrt(lam[g]*p[g]) )^2
      ynew[g,k] ~ dbinom(p[g], N[g])
      Errnew[g,k] <- ( sqrt(ynew[g,k]) - sqrt(lam[g]*p[g]) )^2
    }
  }
  Fit <- sum(Err[,])
  Fitnew <- sum(Errnew[,])
}
")

# Inits
Nst1 <- apply(yg,1,max)+1
Nst2 <- yg + 1
inits <- function(){ list (b = rep(0, nk), sigma = runif(1, 0.2, 0.7),
    alpha0 = runif(1, -5, -2), beta0 = 0, N = Nst1 ) }

# Parameters to monitor
params1 <- c("sigma", "Ntot", "EN", "beta0", "b", "tau.b", "sigma.b",
    "Fit", "Fitnew")

# MCMC settings
na <- 5000 ; ni <- 22000 ; nb <- 2000 ; nt <- 10 ; nc <- 3

# Call JAGS (ART 11 min), assess convergence and summarize posteriors
outdsm1 <- jags (bdata, inits, params1, "DSM1.txt", n.thin = nt, n.chains = nc,
    n.burnin = nb, n.iter = ni, n.adapt = na, parallel = TRUE )
# par(mfrow = c(3,3))  # ~~~ no longer necessary
traceplot(outdsm1)
print(outdsm1, 3)
#            mean     sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# sigma     0.707  0.050   0.614   0.705   0.814    FALSE 1.000 1.010   323
# Ntot    365.986 36.379 323.000 360.000 441.000    FALSE 1.000 1.077   148
# beta0    -0.529  1.481  -3.587  -0.505   2.370     TRUE 0.648 1.001  4201
# b[1]     -7.375 13.562 -35.078  -7.412  19.711     TRUE 0.720 1.000  6000
# b[2]     -7.709 13.803 -36.722  -6.715  18.249     TRUE 0.718 1.001  6000
# b[3]      6.471 14.902 -22.540   6.154  37.049     TRUE 0.673 1.000  6000
# b[4]     -5.363 15.152 -35.726  -5.192  24.440     TRUE 0.653 1.000  4208
# ... Output trunacated ...
# b[82]    -6.761 15.347 -38.737  -6.282  22.858     TRUE 0.675 1.001  2441
# b[83]     5.927 14.723 -22.065   5.322  36.974     TRUE 0.650 1.001  2117
# b[84]    -2.512 15.258 -33.394  -2.290  27.253     TRUE 0.563 1.001  6000
# b[85]     9.873 14.472 -18.391   9.369  40.036     TRUE 0.766 1.001  3436
# b[86]    -0.126 15.299 -29.654  -0.159  31.412     TRUE 0.505 1.000  6000
# tau.b     0.005  0.003   0.002   0.004   0.012    FALSE 1.000 1.011   314
# sigma.b  16.127  4.383   9.018  15.732  25.816    FALSE 1.000 1.005   397
# Fit     307.864 11.782 285.922 307.460 331.636    FALSE 1.000 1.009   250
# Fitnew  305.433 12.883 281.317 305.101 331.283    FALSE 1.000 1.008   266

b <- outdsm1$sims.list$b
beta0 <- mean(outdsm1$sims.list$beta0)
pred <- exp( beta0 + rowMeans(Z.mat%*%t(b) ))

# Quick summary plot (Fig. 11.18)
op <- par(mar = c(3,3,3,6))
image( r <- rasterFromXYZ(cbind(Habgrid,pred)), col = topo.colors(10) )
image_scale(pred, col = topo.colors(10) )
points(tmp$U, pch = 20, cex = 1.5, col = "white")
par(op)

# RMSE and Bayesian p-value
(RMSE <- sqrt(mean( (tmp$N - pred)^2 ) ) )
# [1] 1.053792

(pval <- mean( outdsm1$sims.list$Fit < outdsm1$sims.list$Fitnew ))
# [1] 0.32

# Fig 11.19
plot(log(pred), Habitat, frame = FALSE, xlab = "linear predictor",
    ylab = "Habitat", col = rgb(0, 0, 0, 0.3), cex = 2, pch = 16)
abline(0, 1, lwd = 3)


# 11.10.4 Concluding remarks on density surface models (no code)
