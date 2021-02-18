# Applied hierarchical modeling in ecology - vol.2 - 2021
# Marc Kéry & J. Andy Royle
#
# Chapter 11 : SPATIALLY EXPLICIT DISTANCE SAMPLING ALONG TRANSECTS
# =================================================================
# Code from proofs dated 2020-08-19

library(AHMbook)
library(jagsUI)

# 11.7 Distance sampling models based on pixel frequencies
# ========================================================

# 11.7.1 Fitting the pixel frequency distance sampling model in JAGS
# ------------------------------------------------------------------

# Simulate a data set
library(AHMbook) ; RNGversion("3.5.0")
set.seed(1234, kind = "Mersenne-Twister")
tmp <- simSpatialDSline(N = 200, beta = 1, sigma = 0.25, alpha0 = -3,
    perp = TRUE )                 # perp = TRUE forces p(0) = 1

# Harvest data objects
Habitat <- as.vector(tmp$Habitat)
Habitat <- as.numeric(scale(Habitat))
Habgrid <- tmp$grid
nind <- nrow(tmp$data)
npixels <- length(tmp$Habitat)

# Create a vector of pixel counts and pad it with zeros
yg <- tabulate(tmp$pixel, nbins = npixels)

# Create a covariate: distance between line and pixel center
dist <- abs(Habgrid[,2] - 0.5 )

# Bundle and summarize the data for BUGS
str(bdata <- list(y = yg, Habitat = Habitat, npixels = npixels, dist = dist))
# List of 4
# $ y      : int [1:400] 0 0 0 0 0 0 1 0 0 0 ...
# $ Habitat: num [1:400] -0.401 -0.152 0.441 -0.584 -0.248 ...
# $ npixels: int 400
# $ dist   : num [1:400] 0.45 0.45 0.45 0.45 0.45 0.45 0.45 0.45 ...

# Specify model in BUGS language, here we set p0 = 1.
cat(file = "spatialDSpixel.txt", "
model{

  # Prior distributions
  sigma ~ dunif(0,10)
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)

  p0 <- 1                           # Fix p0 to 1 to honor the CDS assumption
  for(g in 1:npixels){ # Point process intensity
    lam[g] <- exp(beta0 + beta1*Habitat[g])
    p[g] <- p0*exp(- (1/(2*sigma*sigma))*dist[g]*dist[g])
    # y[g] ~ dpois(p[g]*lam[g])     # equivalent to next two lines
    N[g] ~ dpois(lam[g])
    y[g] ~ dbinom(p[g], N[g])
  }
  # Derived parameters
  Ntot <- sum(N[])                  # N is a derived parameter
  D <- Ntot/4                       # Density, with area = 4 units
}
")

# MCMC settings
na <- 1000 ; ni <- 20000 ; nb <- 10000 ; nt <- 10 ; nc <- 3

# Inits
Nst <- yg + 1
inits <- function(){ list (sigma = runif(1, 0.2, 1), beta0 = -1, beta1 = 1,
    N = Nst) }

# Parameters to monitor
params <- c("sigma", "Ntot", "beta0", "beta1", "D")

# Run JAGS (ART < 1 min), check convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "spatialDSpixel.txt", n.adapt = na,
    n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE )
op <- par(mfrow = c(2,3)) ; traceplot(out6)
par(op)
print(out6, 3)
#          mean     sd    2.5%     50%   97.5%  Rhat n.eff
# sigma   0.265  0.035   0.215   0.260   0.348 1.009   457
# Ntot  189.548 18.999 154.000 189.000 228.000 1.008   229
# beta0  -1.063  0.153  -1.370  -1.054  -0.778 1.013   167
# beta1   0.796  0.099   0.604   0.797   0.998 1.002  1401
# D      47.387  4.750  38.500  47.250  57.000 1.008   229


# 11.7.2 Analysis of the pixel frequency distance sampling model in unmarked
#        using pcount.spHDS
# --------------------------------------------------------------------------

# Construct an unmarkedFrame
library(unmarked)
summary(umf <- unmarkedFramePCount(y = matrix(yg, ncol = 1),
    siteCovs = data.frame(dist = dist, Habitat = Habitat)) )

# unmarkedFrame Object

# 400 sites
# Maximum number of observations per site: 1
# Mean number of observations per site: 1
# Sites with at least one detection: 86
#
# Tabulation of y observations:
# 0 1 2 3 4 7
# 314 64 13 7 1 1
#
# Site-level covariates:
# dist Habitat
# Min. :0.05 Min. :-2.39613
# 1st Qu.:0.15 1st Qu.:-0.55058
# Median :0.25 Median : 0.05871
# Mean :0.25 Mean : 0.00000
# 3rd Qu.:0.35 3rd Qu.: 0.56853
# Max. :0.45 Max. : 2.52070
#
# Fit the hierarchical distance sampling model
(fm1 <- pcount.spHDS(~ -1 + I(dist^2) ~ Habitat, umf))

# Abundance:
# Estimate SE z P(>|z|)
# (Intercept) -1.25 0.1534 -8.12 4.78e-16
# Habitat 1.07 0.0973 11.00 3.66e-28
#
# Detection:
# Estimate SE z P(>|z|)
# 2.49 0.16 15.6 9.01e-55
#
# AIC: 434.9047

# Generate predictions for each pixel
pred <- predict(fm1, type = "state")

# Make a quick plot to visualize the result (Fig. 11.7)
library(raster)
M1 <- matrix(pred[,1], nrow = 10, byrow = T)
M2 <- M1[nrow(M1):1, ]
op <- par( mar = c(3, 3, 3, 6) )
graphics::image(x = 1:40, y = 1:10, t(M2), col = topo.colors(12))
image_scale(t(M2), col = topo.colors(12))           # AHMbook package
par(op)

# Estimated total population size in the landscape
( Nhat <- sum(pred[,1]) )
# [1] 199.7485

log(200/sum(exp( 1*Habitat)))
# [1] -1.175485
