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

library(AHMbook)

# 9.2 Data simulation for spatial N-mixture and occupancy models
# ==============================================================

library(AHMbook)
? simExpCorrRF
str(dat <- simExpCorrRF(variance = 1, theta = 1, size = 50, seed = 1))

str(tmp <- simExpCorrRF(theta = 0.0001, size = 200))
str(tmp <- simExpCorrRF(theta = 1, size = 200))
str(tmp <- simExpCorrRF(theta = 5, size = 200))
str(tmp <- simExpCorrRF(theta = 10, size = 200))
try(str(tmp <- simExpCorrRF(theta = 100, size = 200))) # fails
try(str(tmp <- simExpCorrRF(theta = 10000, size = 200))) # fails

data(BerneseOberland)
head(bo <- BerneseOberland)
str(bo)

# ~~~~ extra code for figure 9.1 ~~~~~~~~~~
library(raster)
op <- par(mfrow = c(1, 2), mar = c(3,3,3,5), cex.main = 2, cex.axis = 1.5)
r1 <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = bo$elevation))
plot(r1, col = topo.colors(100), axes = FALSE, box = FALSE,
    main = "Elevation (m a.s.l.)", zlim = c(0, 3000))
r1 <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = bo$forest))
plot(r1, col = topo.colors(100), axes = FALSE, box = FALSE,
    main = "Forest cover (%)", zlim = c(0, 100))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create Gaussian random field
set.seed(10) # Fig. 9.2
s <- simExpCorrRF(theta = 10, size = 50)

# Choose sample sizes: 2500 sites and 3 surveys
nsites <- 2500 # Number of sites (corresponding to our 50 by 50 grid)
nreps <- 3     # Number of replicate observations

# Scale the real Bernese Oberland covariates
elev <- standardize(bo$elevation)
forest <- standardize(bo$forest)

# Ecological process
beta0 <- 2           # Abundance model: intercept
beta1 <- 2           # Linear effect of elevation: positive
beta2 <- -2          # Quadratic effect of elevation: negative
loglam0 <- beta0 + beta1 * elev + beta2 * elev^2
loglam <- beta0 + beta1 * elev + beta2 * elev^2 + c(s$field)
lam0 <- exp(loglam0) # without spatial autocorrelation
lam <- exp(loglam)   # with spatial autocorrelation

# ~~~~ extra code for figure 9.3 ~~~~
# Plot expected counts (lambda) as a function of covariates only,
#   i.e., excluding spatial field, and including spatial field (Fig. 20-4)
op <- par(mfrow = c(1,2), mar = c(5,8,5,2), cex.lab = 1.5)
plot(bo$elevation, lam0, cex = 1, pch = 16, xlab = "Elevation",
    ylab = "Expected counts (lambda)", main = "Excluding spatial field",
    frame = FALSE, col = rgb(0, 0, 0, 0.3))
plot(bo$elevation, lam, cex = 1, pch = 16, xlab = "Elevation",
    ylab = "Expected counts (lambda)", main = "Including spatial field",
    frame = FALSE, col = rgb(0, 0, 0, 0.3))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine actual abundances as Poisson random variables with parameter lam
N <- rpois(n = nsites, lambda = lam)
table(N)            # Distribution of abundances across sites
sum(N > 0) / nsites # Finite-sample occupancy
(totalN <- sum(N))  # Total population size in all 2500 sites

# Create wind speed observational covariate
wind <- matrix(rnorm(nsites*nreps), nrow = nsites, ncol = nreps)

# Observation process
alpha0 <- 0   # logit-linear intercept
alpha1 <- -1  # slope on forest
alpha2 <- -1  # slope on wind speed
p <- array(NA, dim = c(nsites, nreps))
for(j in 1:nreps){
  p[,j] <- plogis(alpha0 + alpha1 * forest + alpha2 * wind[,j])
}

# Count things
y <- array(dim = c(nsites, nreps)) # Array for counts
for (j in 1:nreps){
  y[,j] <- rbinom(n = nsites, size = N, prob = p[,j])
}
str(y)
# int [1:2500, 1:3] 1 0 0 0 0 0 0 3 0 0 ...
summary(N)
summary(c(y))
#  Min. 1st Qu. Median  Mean 3rd Qu.   Max.
# 0.000   0.000  1.000 2.877   3.000 66.000
#  Min. 1st Qu. Median  Mean 3rd Qu.   Max.
# 0.000   0.000  0.000 1.174   1.000 49.000

# Compare true and observed total abundance
(true <- totalN)                  # True
(obs <- sum(apply(y, 1, max)))    # Observed
cat("Underestimation of total abundance:", round(100*(1-obs/true)), "%\n")
# [1] 7192
# [1] 4371
# Underestimation of total abundance: 39 %

# Select a sample of sites for surveys
# set.seed(100)
set.seed(100, sample.kind = "Rounding")
sample.size <- 500
sample.sites <- sort(sample(1:nsites, size = sample.size))

# ~~~~ extra code for figure 9.4 ~~~~
op <- par(mfrow = c(1,3), mar = c(3,3,3,6))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = N))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Abundance (N, truncated at 6)", zlim = c(0, 6))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = apply(p, 1, mean)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Average detection probability")
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = apply(y, 1, max)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
    main = "Max count (truncated at 6)", zlim = c(0, 6))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yobs <- y                    # Make a copy
yobs[-sample.sites,] <- NA   # Turn counts of unsurveyed sites into NAs
head(sample.sites)           # Look at the simulated data set
head(yobs)

simNmixSpatial(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
    mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
    theta.RF = 10, seeds = c(10, 100), truncN = 6, show.plots = TRUE)

simOccSpatial(nsurveys = 3, mean.psi = 0.6, beta = c(2, -2),
    mean.p = 0.4, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
    theta.RF = 10, seeds = c(10, 100), show.plots = TRUE)
