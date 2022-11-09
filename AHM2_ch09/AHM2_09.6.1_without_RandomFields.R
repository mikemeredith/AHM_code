#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================
# Code from proofs dated 2020-08-19

# The code below is modified to run without the RandomFields package; if
#   RandomFields is available, it will be used by AHMbook, and the results should
#   be the same.
# When RandomFields is not available, the 'fields' package is used instead, and
#   the results will be different.
if(requireNamespace("RandomFields"))
  stop("Package 'RandomFields' IS available; this script is not needed.")

# Approximate execution time for this code: 16 mins
# Run time with the full number of iterations: 1.5 hrs

library(AHMbook)
library(jagsUI)


# 9.6 Mechanistic, or dynamic, models of spatial autocorrelation
# ==============================================================

# 9.6.1 Autologistic dynamic occupancy models
# -------------------------------------------

# 9.6.1.1 Data simulation under the autologistic dynocc model
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# The code below is modified to run without the RandomFields package; if
#   RandomFields is available, it will be used by AHMbook, and the results should
#   be the same.
# When RandomFields is not available, the 'fields' package is used instead, and
#   the results will be different.

library(AHMbook)
# This call to simDynoccSpatial will fail if the 'RandomFields' package is not available
try(str(dat <- simDynoccSpatial(
    side = 10, nyears = 10, nsurveys = 3,               # sample sizes
    mean.psi1 = 0.4, beta.Xpsi1 = 0,                    # linear model psi1
    range.phi = c(0.8, 0.8), beta.Xphi = 0,             # ... phi
    range.gamma = c(0.1, 0.1), beta.Xgamma = 0,         # ... gamma
    range.p = c(0.4, 0.4), beta.Xp = 0,                 # ... p
    theta.XAC = 5000, beta.XAC = c(0, 0, 0, 0),         # SAC effects 1
    beta.Xautolog = c(0, 0),                            # SAC effects 2
    trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0), # Heterogeneity p
    seed.XAC = NA, seed = NULL, ask.plot = TRUE) ) )

# Take function output and fit a non-spatial dynocc model in unmarked
library(unmarked)
# Example 1: model with 1 covariate in each main parameter
try(str( dat <- simDynoccSpatial(beta.Xpsi = 1, beta.Xphi = 1, beta.Xgamma = 1,
    beta.Xp = 1) ) )
# summary(dat$umf)
# summary(fm0 <- colext(~1, ~1, ~1, ~1, dat$umf))              # constant model
# summary(fm1 <- colext(~Xpsi1, ~Xgamma, ~Xphi, ~Xp, dat$umf)) # data-generating model

# Example 2: Generate data with autocovariate effects
try(str( dat <- simDynoccSpatial(mean.psi1 = 0.1, beta.Xautolog = c(3, 3),
    range.p = c(0.1, 0.1), ask.plot = FALSE, seed.XAC = 1, seed = 1)) )
# summary(dat$umf)
# summary(fm1 <- colext(~1, ~Xautoobs, ~Xautoobs, ~1, dat$umf)) # 'naive' autologistic
# summary(fm2 <- colext(~1, ~Xauto, ~Xauto, ~1, dat$umf)) # 'true' autologistic
# confint(fm1, type = 'col') ; confint(fm1, type = 'ext') # CIs
# confint(fm2, type = 'col') ; confint(fm2, type = 'ext')

# 9.6.1.2 fitting the simplest autologistic dynocc model to simulated data
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# This fails without 'RandomFields':
# str(dat <- simDynoccSpatial(side = 30, nyears = 10, nsurveys = 3,
    # mean.psi1 = 0.1, range.phi = c(0.5, 0.5), range.gamma = c(0.2, 0.2),
    # range.p = c(0.4, 0.4), beta.Xautolog = c(1, 1), seed.XAC = 1,
    # seed = 24, ask.plot = TRUE) )

# Use the precooked simulated data
data(simDynoccSpatialData)
dat <- simDynoccSpatialData
library(unmarked)
summary(dat$umf)
fm <- colext(~1, ~Xauto, ~Xauto, ~1, dat$umf)
summary(fm)
confint(fm, type = 'col')[2,]
confint(fm, type = 'ext')[2,]

# Format detection/nondetection data in a 3D array
nsites <- dat$side^2
nsurveys <- dat$nsurveys
nyears <- dat$nyears
y <- array(NA, dim = c(nsites, nsurveys, nyears))
for(i in 1:nyears){
  y[,,i] <- dat$umf@y[,(3*i-2):(3*i)]
}

# Grab and look at the adjacency matrix: marks neighbors with a 1
amat <- dat$amatrix
table(apply(amat, 2, sum))     # Frequency of number of neighbors

# Compute the neighborhood info
library(spdep)
neigh <- dnearneigh(dat$grid, d1 = 0, d2 = sqrt(2) + 0.1)
str(winnb <- nb2WB(neigh))    # Function to get CAR ingredients for BUGS
numN <- winnb$num             # Number of neighbors for each cell

# Put the neighbor IDs into a matrix
neighID <- array(NA, dim = c(nsites, 8))
for(i in 1:nsites){
  neighID[i, 1:numN[i]] <- unlist(neigh[i])
}
head(neighID) # Site 1 is a neighbor of 2, 31 and 32 etc.
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# [1,]    2   31   32   NA   NA   NA   NA   NA
# [2,]    1    3   31   32   33   NA   NA   NA
# [3,]    2    4   32   33   34   NA   NA   NA
# [4,]    3    5   33   34   35   NA   NA   NA
# [5,]    4    6   34   35   36   NA   NA   NA
# [6,]    5    7   35   36   37   NA   NA   NA

str(bdata <- list(y = y, nsites = nsites, nsurveys = nsurveys, nyears = nyears,
    neighID = neighID, numN = numN))
# List of 6
# $ y       : int [1:900, 1:3, 1:10] 0 0 0 0 1 0 0 0 0 0 ...
# $ nsites  : num 900
# $ nsurveys: num 3
# $ nyears  : num 10
# $ neighID : int [1:900, 1:8] 2 1 2 3 4 5 6 7 8 9 ...
# $ numN    : int [1:900] 3 5 5 5 5 5 5 5 5 5 ...


# Specify model in BUGS language
cat(file = "autologistic1.txt","
model {

  # Priors
  psi1 ~ dunif(0, 1)                   # Initial occupancy
  phi.int ~ dunif(0, 1)                # Persistence
  alpha.lphi <- logit(phi.int)
  beta.lphi ~ dnorm(0, 0.01)
  gamma.int ~ dunif(0, 1)              # Colonization
  alpha.lgamma <- logit(gamma.int)
  beta.lgamma ~ dnorm(0, 0.01)
  p ~ dunif(0, 1)                      # Detection

  # Likelihood
  # Ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1] + (1-z[i,t-1]) * gamma[i,t-1])
      # Compute autocovariate and specify its effects on phi and gamma
      autocov[i,t-1] <- sum(z[neighID[i,1:numN[i]], t-1]) / numN[i]
      logit(phi[i,t-1]) <- alpha.lphi + beta.lphi * autocov[i,t-1]
      logit(gamma[i,t-1]) <- alpha.lgamma + beta.lgamma * autocov[i,t-1]
    }
  }

  # Observation model
  for (i in 1:nsites){
    for (j in 1:nsurveys){
      for (t in 1:nyears){
        y[i,j,t] ~ dbern(z[i,t] * p)
      }
    }
  }
}
")

# Initial values
zst <- array(1, dim = c(nsites, nyears))
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi1", "phi.int", "alpha.lphi", "beta.lphi", "gamma.int",
    "alpha.lgamma", "beta.lgamma", "p") # could also monitor "z"

# MCMC settings
# na <- 1000 ; ni <- 6000 ; nt <- 3 ; nb <- 3000 ; nc <- 3
na <- 1000 ; ni <- 600 ; nt <- 1 ; nb <- 300 ; nc <- 3  # ~~~ for testing, 8 mins

# Call JAGS (ART 30 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "autologistic1.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
print(out1, 3)

# Compare estimates with truth
truth <- rbind('psi1' = dat$mean.psi1, 'phi.int' = dat$range.phi[1],
    'alpha.lphi' = qlogis(dat$range.phi[1]), 'beta.lphi' = dat$beta.Xautolog[1],
    'gamma.int' = dat$range.gamma[1], 'alpha.lgamma' = qlogis(dat$range.gamma[1]),
    'beta.lgamma' = dat$beta.Xautolog[2], 'p' = dat$range.p[1])
print(cbind(truth, out1$summary[1:8, c(1:3,7)]), 3)
#              truth   mean      sd    2.5%  97.5%
# psi1          0.10  0.088 0.01068  0.0682  0.110
# phi.int       0.50  0.531 0.03679  0.4633  0.601
# alpha.lphi    0.00  0.125 0.14869 -0.1472  0.411
# beta.lphi     1.00  1.020 0.34297  0.3644  1.676
# gamma.int     0.20  0.198 0.01264  0.1740  0.223
# alpha.lgamma -1.39 -1.400 0.07974 -1.5579 -1.250
# beta.lgamma   1.00  1.060 0.20979  0.6588  1.491
# p             0.40  0.400 0.00712  0.3861  0.414


# 9.6.1.3 Fitting an autologistic dynocc model to the Eurasian lynx data
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Get the data
library(AHMbook)
data(EurasianLynx)
str(datfull <- EurasianLynx)
year <- 1994:2016

# Select certain data in Switzerland and format det/nondet data
selection <- datfull$type == 'certain' & datfull$Cntry == 'Switzerland'
str(dat <- datfull[selection,])

# Format detection/nondetection data in a 3D array
( nsites <- length(unique(dat$site.nr)) ) # 367
( nyears <- length(unique(dat$Year)) )    # 23
( nsurveys <- 3)
ylong <- as.matrix(dat[,3:5])
y <- array(NA, dim = c(nsites, nsurveys, nyears))
for(i in 1:nyears){
  y[,,i] <- ylong[(367*i-366):(367*i),]
}

# Get coordinate info for the Swiss lynx (kilometre units)
head( grid <- cbind(dat$xcoord[1:nsites], dat$ycoord[1:nsites]) )

# Get neighborhood info
library(spdep)
neigh <- dnearneigh(grid, d1 = 0, d2 = sqrt(2)*10 + 0.1)
str(winnb <- nb2WB(neigh)) # Function to get CAR ingredients for BUGS
numN <- winnb$num          # Number of neighbors for each cell
table(numN)
# numN
# 1 2  3  4  5  6  7   8
# 2 2 11 21 29 27 34 241

# Put the neighbor IDs into a matrix
neighID <- array(NA, dim = c(nsites, 8))
for(i in 1:nsites){
  neighID[i, 1:numN[i]] <- unlist(neigh[i])
}
head(neighID)
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# [1,]    2    4    5   NA   NA   NA   NA   NA
# [2,]    1    4    5    6   NA   NA   NA   NA
# [3,]   11   NA   NA   NA   NA   NA   NA   NA
# [4,]    1    2    5   12   13   14   NA   NA
# [5,]    1    2    4    6   13   14   15   NA
# [6,]    2    5    7   14   15   16   NA   NA

# Grab environmental covariate (forest) and standardize
oforest <- dat$forest[1:nsites] # Original forest
forest <- standardize(oforest)

# Compute observed lynx occurence in every year
zobs <- apply(y, c(1,3), max, na.rm = TRUE)

# Plot observed occurrence on top of forest map (Fig. 9.16 shows part)
library(raster)
mapPalette <- colorRampPalette(c("grey", "lightgreen", "darkgreen"))
for(t in 1:nyears){
  r <- rasterFromXYZ(data.frame(x = grid[,1], y = grid[,2], z = oforest))
  plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, zlim = c(0, 100),
  main = paste('Observed lynx occurrence per 10x10km2 in', year[t],
      '\noverlaid on map of forest cover (%)'), axis.args=list(cex.axis=0.8),
  legend.width = 0.8)
  points(grid[,1][zobs[,t] == 0], y = grid[,2][zobs[,t] == 0], pch = 22)
  points(grid[,1][zobs[,t] == 1], y = grid[,2][zobs[,t] == 1], pch = 15)
  # browser() # ~~~ take out for testing
}

# Bundle and summmarize data set
str(bdata <- list(y = y, nsites = nsites, nsurveys = nsurveys, nyears = nyears,
    neighID = neighID, numN = numN, forest = forest))
# List of 7
# $ y       : int [1:367, 1:3, 1:23] NA NA NA NA NA NA NA NA NA NA ...
# $ nsites  : int 367
# $ nsurveys: num 3
# $ nyears  : int 23
# $ neighID : int [1:367, 1:8] 2 1 11 1 1 2 6 7 8 9 ...
# $ numN    : int [1:367] 3 4 1 6 7 6 5 5 5 4 ...
# $ forest  : num [1:367] -0.97 -0.865 1.448 -0.769 0.26 ...

# Specify model in BUGS language
cat(file = "autologistic2.txt","
model {

  # Priors
  psi1.int ~ dunif(0, 1) # Initial occupancy
  alpha.lpsi1 <- logit(psi1.int)
  beta.lpsi1.forest ~ dnorm(0, 0.1)
  phi.int ~ dunif(0, 1) # Persistence
  alpha.lphi <- logit(phi.int)
  beta.lphi.forest ~ dnorm(0, 0.1)
  beta.lphi.auto ~ dnorm(0, 0.01)
  gamma.int ~ dunif(0, 1) # Colonization
  alpha.lgamma <- logit(gamma.int)
  beta.lgamma.forest ~ dnorm(0, 0.1)
  beta.lgamma.auto ~ dnorm(0, 0.01)
  p.int ~ dunif(0, 1) # Detection
  alpha.lp <- logit(p.int)
  beta.lp.forest ~ dnorm(0, 0.1)

  # Likelihood
  # Ecological submodel
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1[i])
    logit(psi1[i]) <- alpha.lpsi1 + beta.lpsi1.forest * forest[i]
    for (t in 2:nyears){
      z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1] + (1-z[i,t-1]) * gamma[i,t-1])
      # Compute autocovariate and specify its effects on phi and gamma
      autocov[i,t-1] <- sum(z[neighID[i,1:numN[i]], t-1]) / numN[i]
      logit(phi[i,t-1]) <- alpha.lphi + beta.lphi.forest * forest[i] +
      beta.lphi.auto * autocov[i,t-1]
      logit(gamma[i,t-1]) <- alpha.lgamma + beta.lgamma.forest * forest[i] +
      beta.lgamma.auto * autocov[i,t-1]
    }
  }

  # Observation model
  for (i in 1:nsites){
    logit(p[i]) <- alpha.lp + beta.lp.forest * forest[i]
    for (j in 1:nsurveys){
      for (t in 1:nyears){
        y[i,j,t] ~ dbern(z[i,t] * p[i])
      }
    }
  }
}
")

# Initial values
zst <- array(1, dim = c(nsites, nyears)) # Cheap inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi1.int", "alpha.lpsi1", "beta.lpsi1.forest", "phi.int",
    "alpha.lphi", "beta.lphi.forest", "beta.lphi.auto", "gamma.int",
    "alpha.lgamma", "beta.lgamma.forest", "beta.lgamma.auto", "p.int",
    "alpha.lp", "beta.lp.forest", "z")

# MCMC settings
# na <- 1000 ; ni <- 6000 ; nt <- 3 ; nb <- 3000 ; nc <- 3
na <- 100 ; ni <- 600 ; nt <- 1 ; nb <- 300 ; nc <- 3  # ~~~ for testing, 6 mins

# Call JAGS (ART 48 min), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "autologistic2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out2)
print(out2$summary[1:14, c(1:4,7)], 2)
#                      mean     sd   2.5%    25%  97.5%
# psi1.int            0.062 0.0202  0.027  0.047  0.106
# alpha.lpsi1        -2.782 0.3701 -3.580 -3.016 -2.132
# beta.lpsi1.forest   0.777 0.3230  0.173  0.555  1.447
# phi.int             0.316 0.0449  0.228  0.286  0.404
# alpha.lphi         -0.782 0.2112 -1.222 -0.916 -0.387
# beta.lphi.forest    0.002 0.1606 -0.301 -0.110  0.321
# beta.lphi.auto      4.495 0.4586  3.642  4.175  5.432
# gamma.int           0.023 0.0026  0.019  0.022  0.029
# alpha.lgamma       -3.737 0.1118 -3.958 -3.810 -3.512
# beta.lgamma.forest  0.378 0.0818  0.214  0.323  0.540
# beta.lgamma.auto    5.157 0.3448  4.502  4.928  5.860
# p.int               0.337 0.0121  0.313  0.329  0.360
# alpha.lp           -0.675 0.0544 -0.784 -0.711 -0.574
# beta.lp.forest      0.600 0.0519  0.500  0.565  0.704

# ~~~~~~~ extra code for figure 9.17 ~~~~~~~~~~
# Compute posterior means and CRIs of autologistic predictions
tmp <- out2$sims.list
nsims <- out2$mcmc.info$n.samples
pred.autocov <- seq(0, 1, length.out = 1000)
pred <- array(NA, dim = c(1000, 2, nsims))
for(i in 1:nsims){
  pred[,1,i] <- plogis(tmp$alpha.lphi[i] + tmp$beta.lphi.auto[i] * pred.autocov)
  pred[,2,i] <- plogis(tmp$alpha.lgamma[i] +
      tmp$beta.lgamma.auto[i] * pred.autocov)
}
pm <- apply(pred, c(1,2), mean)
CRI <- apply(pred, c(1,2), quantile, prob = c(0.025, 0.975))

op <- par(mfrow = c(1,2))
plot(pred.autocov, pm[,1], type = 'l', lty = 1, lwd = 3, col = 'blue',
    xlab = 'Proportion of neighboring quadrats occupied',
    ylab = ' Persistence probability', ylim = c(0, 1), frame = FALSE)
polygon(c(pred.autocov, rev(pred.autocov)), c(CRI[1,,1], rev(CRI[2,,1])),
    col = 'grey', border = NA)
lines(pred.autocov, pm[,1], type = 'l', lty = 1, lwd = 3, col = 'blue')
plot(pred.autocov, pm[,2], type = 'l', lty = 1, lwd = 3, col = 'blue',
    xlab = 'Proportion of neighboring quadrats occupied',
    ylab = 'Colonization probability', ylim = c(0, 1), frame = FALSE)
polygon(c(pred.autocov, rev(pred.autocov)), c(CRI[1,,2], rev(CRI[2,,2])),
    col = 'grey', border = NA)
lines(pred.autocov, pm[,2], type = 'l', lty = 1, lwd = 3, col = 'blue')
par(op)

# ~~~~~~~ extra code for figure 9.18 ~~~~~~~~~~
# Plot estimated SDM in first and last years (i.e., 1994 and 2016)
library(raster)
op <- par(mfrow = c(1,2), mar = c(1,1,3,5))
r <- rasterFromXYZ(data.frame(x = grid[,1], y = grid[,2], z = oforest))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, zlim = c(0, 100),
    main = 'Winter 1994/1995')
points(grid[,1], y = grid[,2], pch = 15, col = rgb(0,0,0, out2$mean$z[,1]))
r <- rasterFromXYZ(data.frame(x = grid[,1], y = grid[,2], z = oforest))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, zlim =
  c(0, 100), main = 'Winter 2016/2017')
points(grid[,1], y = grid[,2], pch = 15, col = rgb(0,0,0, out2$mean$z[,23]))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
