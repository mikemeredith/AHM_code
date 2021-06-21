# Applied hierarchical modeling in ecology - vol.2 - 2021
# Marc Kéry & J. Andy Royle
#
# Chapter 11 : SPATIALLY EXPLICIT DISTANCE SAMPLING ALONG TRANSECTS
# =================================================================
# Code from proofs dated 2020-08-19

print("Approximate execution time for this code: 30 mins")
# Run time with the full number of iterations: 75 mins

library(AHMbook)
library(jagsUI)

# 11.8 Temporary emigration distance sampling models
# ==================================================

# 11.8.1 The TEDS/TPP model of Mizel et al. (2018) (no code)

# 11.8.2 Simulating data under the TEDS/TPP model
# ------------------------------------------------

set.seed(12345, kind = "Mersenne-Twister" )
str(tmp <- simSpatialDSte(nsites = 33, lam0 = 0.3, beta = 2, phi = 0.6,
    nsurveys = 4))                        # Produces (in part) Fig. 11.8.
# output omitted

# Harvest the data objects
y <- tmp$y
d <- tmp$d
B <- tmp$B
nobs <- tmp$Counts
Habitat <- tmp$Habitat
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
M <- tmp$M
nsites <- tmp$nsites
nsurveys <- tmp$nsurveys
npixels <- tmp$npixels
nind <- sum(M)                            # Total superpopulation size

# 11.8.3 Fitting the TEDS/TPP model in JAGS
# -----------------------------------------

# Bundle the data
str(bdata <- list(y = y, d = d, nsites = nsites, Habitat = Habitat, nobs = nobs,
    npixels = npixels, nsurveys = nsurveys))
# List of 7
# $ y       : int [1:33, 1:100, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
# $ d       : num [1:100] 4.5 3.5 2.5 1.5 0.5 0.5 1.5 2.5 3.5 4.5 ...
# $ nsites  : num 33
# $ Habitat : num [1:100, 1:33] 0.6253 0.959 0.5104 -0.0354 0.4759 ...
# $ nobs    : int [1:33, 1:4] 4 0 1 1 0 5 1 7 1 0 ...
# $ npixels : num 100
# $ nsurveys: num 4

# Specify HDS model with TE in BUGS language (from Mizel et al. 2018)
cat(file = "TEDSmodel.txt", "
model {
  # Prior distributions
  sigma ~ dunif(0,25)
  phi ~ dunif(0,1)
  beta1 ~ dnorm(0,.01)
  beta0 ~ dnorm(0,.01)
  for (j in 1:npixels){
    # Distance-based detection model
    log(g[j]) <- -d[j] * d[j] / (2*sigma*sigma)
    # Spatial abundance model
    for (i in 1:nsites){
      log(lam[j, i]) <- beta0 + beta1 * Habitat[j, i]
      pix.probs[j, i] <- lam[j, i] / sum(lam[, i])
      # conditional cell probabilities:
      # Pr( detected|pixel ) Pr(pixel | N)
      # Basic rule: [sum_{pix} Pr(det | pix)*Pr(pix | N)] = Pr(det | N)
      cellprobs[j, i] <- g[j] * pix.probs[j, i]         # pi in algebra above
      # Following line normalizes the probabilities but dmulti does this
      # cellprobs.cond[j,i] <-
      # cellprobs[j,i]/sum(cellprobs[1:npixels,i] )
    }
  }
  # Now define the availability process and observation model
  for (i in 1:nsites) {
    # Last cell probability (not used in conditional model)
    cellprobs[npixels + 1, i] <- 1-sum(cellprobs[1:npixels, i])
    # Part 1 of the hierarchical model: Plot level abundance model
    M.lam[i] <- sum(lam[, i])
    M[i] ~ dpois(M.lam[i])
    pdet[i] <- sum(cellprobs[1:npixels, i])         # Pr(detected | N)
    pmarg[i] <- pdet[i] * phi                       # Marginal probability
    for (k in 1:nsurveys) {
      # 2-part model
      # nobs[i,k] ~ dbin(pdet[i,k],N[i,k])          # same as model just below
      # N[i,k] ~ dbin(phi, M[i])
      # combined together produce:
      # Part 2 of the hierarchical model
      nobs[i,k] ~ dbin(pmarg[i], M[i])
      # Part 3 of the hierarchical model
      y[i,1:npixels,k] ~ dmulti(cellprobs[1:npixels,i], nobs[i,k])
      # Note: cell probs don't add to 1 but dmulti standardizes them
    }
  }
  # Total population size across all plots
  Mtotal <- sum(M[])
}
")

# Inits function
inits <- function(){ list (sigma = 3, M = M, beta1 = 2, phi = 0.6, beta0 = -5.8) }

# Parameters to estimate
params <- c("sigma", "beta0", "beta1", "phi", "Mtotal")

# MCMC settings
na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nt <- 5 ; nc <- 3

# Run JAGS (ART 6 min), check convergence and summarize posteriors
out7 <- jags(bdata, inits, params, "TEDSmodel.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out7)
print(out7, digits = 2)
#          mean     sd   2.5%    50%   97.5% overlap0 f  Rhat n.eff
# sigma   2.782  0.436  2.198  2.702   3.838    FALSE 1 1.011  1562
# beta0  -5.587  0.284 -6.127 -5.592  -5.013    FALSE 1 1.001  4866
# beta1   1.931  0.094  1.748  1.931   2.113    FALSE 1 1.000  6000
# phi     0.613  0.118  0.365  0.619   0.829    FALSE 1 1.002  1074
# Mtotal 91.680 20.723 68.000 87.000 146.000    FALSE 1 1.003  3009

sum(tmp$M)
# [1] 83


# 11.8.4 Analysis of the tree sparrow data from Mizel et al. (2018)
# -----------------------------------------------------------------

# Load the data set and harvest some data
library(AHMbook)
library(raster)
data("treeSparrow")
str(ts <- treeSparrow)
maxVisit <- 5                 # Max number of surveys
nsites <- 150                 # Plots with a transect running through them
totsurvs <- 466

ndvi <- ts$pixels[,c("X", "Y", "NDVI")]
plot(rasterFromXYZ(ndvi))

summary(ts$pixels)                # not shown
summary(ts$obs)
summary(ts$surveyData)
site <- ts$pixels$Site            # site id
sitepix <- as.vector(table(site)) # number pixels in each plot
( max.pix <- max(sitepix) )       # max pixels across plots
d <- round(ts$pixels$NEAR_DIST,2) # dist to pixel centroid FROM TRANSECT
pixID <- ts$pixels$pixID          # pixel id vectorized

# Landscape covariates at pixel level
elev <- round(standardize(ts$pixels$elev),5)
elev2 <- elev*elev
elev2 <- round(elev2,5)
NDVI <- round(standardize(ts$pixels$NDVI),5)
NDVI2 <- NDVI*NDVI
( npixels <- length(NDVI) )       #total number of pixels in study area: 49250

# Take a look at the NDVI landscape (cool plot !)
op <- par(mfrow = c(1,2))
library(raster)
plot(rasterFromXYZ(ndvi), frame = FALSE, axes = FALSE)
ndvi.sub <- ndvi[site==83,]
xbar <- mean(ndvi.sub[,1])
ybar <- mean(ndvi.sub[,2])
library(plotrix)
draw.circle(xbar, ybar, 600, nv = 1000, border = NULL, col = NA,
    lty = 1, lwd = 1)
plot(rasterFromXYZ(ndvi.sub), frame = FALSE, axes = FALSE) # Fig. 11.10
par(op)

# Survey covariates, scale, square if necessary
sc <- as.data.frame( scale(ts$surveyData[, 5:7]) )
julian <- sc$juldate            # integer day
julian2 <- julian^2
time <- sc$reltime              # time of day
effort <- sc$effort             # length of transect through plot
site2 <- ts$surveyData$Site
nobs <- ts$surveyData$count

# Individual bird observations (325 total detections)
ind.site <- ts$obs$Site         # site id for ATSP observations
pixel.cap <- ts$obs$Pixel       # pixel id for ATSP observations
ind.visit <- ts$obs$Visit       # visit id for ATSP observations
ind <- ts$obs$ind               # id of bird detection within the survey
( Mind <- max(ind) )            # max number of birds observed in a survey
survey <- ts$obs$SurveyID

# Encounter frequency array
pixel <- array(0, dim = c(Mind, totsurvs))
for (i in 1:length(pixel.cap)){
  pixel[ind[i], survey[i]] <- pixel.cap[i]
}

# Pixel-level detection frequencies (y) from detection locations
y <- array(0, dim = c(totsurvs, max.pix), dimnames = list( c(1:totsurvs),
    c(1:max.pix)))
for(j in 1:totsurvs){
  y[j,] <- tabulate(pixel[,j], nbins = max.pix)
}

# Initial values for super-population size
M <- tapply(ts$surveyData$count, ts$surveyData$Site, max) + 1

# Bundle and summarize data (data courtesy of Mizel et al. 2018)
str(bdata <- list(y = y, d = d, nsites = nsites, NDVI = NDVI, nobs = nobs,
    npixels = npixels, site = site, pixID = pixID, julian = julian,
    NDVI2 = NDVI2, julian2 = julian2, time = time, effort = effort,
    site2 = site2, totsurvs = totsurvs, sitepix = sitepix, elev = elev,
    elev2 = elev2) )         # Not shown

# Specify TEDS/TPP model in BUGS language
cat(file="TEDS.txt","
model {
  # Prior distributions
  sigma ~ dunif(0, 500)
  alpha0 ~ dnorm(0, 0.01)
  alpha1 ~ dnorm(0, 0.01)
  alpha2 ~ dnorm(0, 0.01)
  alpha3 ~ dnorm(0, 0.01)
  alpha4 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)
  beta4 ~ dnorm(0, 0.01)

  # Define the pixel frequency probabilities
  for (i in 1:npixels){
    # Detection model
    log(g[i]) <- -d[i] * d[i] / (2*sigma*sigma)
    # Density model
    log(lam[pixID[i], site[i]]) <- beta0 + beta1 * NDVI[i] +
    beta2 * NDVI2[i] + beta3 * elev[i] + beta4 * elev2[i]
    pix.probs[pixID[i], site[i]] <- lam[pixID[i], site[i]] /
    sum(lam[1:sitepix[site[i]], site[i]])
    cellprobs[pixID[i], site[i]] <- g[i] * pix.probs[pixID[i], site[i]]
  }
  # Super-population size model
  for (j in 1:nsites) {
    M.lam[j] <- sum(lam[1:sitepix[j], j])
    M[j] ~ dpois(M.lam[j])
  }
  # Categorical/multinomial probabilities
  for (i in 1:totsurvs) {
    # Availability model
    logit(phi[i]) <- alpha0 + alpha1 * time[i] + alpha2 * julian[i] +
    alpha3 * julian2[i] + alpha4 * effort[i]
    # Observation model
    pdet[i] <- sum(cellprobs[1:sitepix[site2[i]], site2[i]])
    pmarg[i] <- pdet[i] * phi[i] # Marginal probability
    nobs[i] ~ dbin(pmarg[i], M[site2[i]])
    y[i, 1:sitepix[site2[i]]] ~ dmulti(cellprobs[1:sitepix[site2[i]],
    site2[i]], nobs[i])
  }
  Mtotal <- sum(M[])
}
")

# Inits function
inits <- function(){ list (M = M+1, beta1 = rnorm(1), beta0 = runif(1,-6,-4),
    beta2 = rnorm(1), beta3 = rnorm(1), beta4 = rnorm(1), alpha0 = rnorm(1),
    alpha1 = rnorm(1), alpha2 = rnorm(1), alpha3 = rnorm(1), alpha4 = rnorm(1),
    sigma = runif(1, 50, 200) )}

# Parameters to monitor
params <- c("sigma", "beta0", "beta1", "beta2", "beta3", "beta4",
    "alpha0", "alpha1", "alpha2", "alpha3", "alpha4", "Mtotal")

# MCMC settings
# na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nt <- 2 ; nc <- 4
na <- 1000 ; ni <- 300 ; nb <- 100 ; nt <- 1 ; nc <- 3  # ~~~~ for testing, 26 mins

# Run JAGS (ART 83 min), check convergence and summarize posteriors
out8 <- jags(bdata, inits, params, "TEDS.txt", n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out8)
print(out8, 3)
#          mean    sd   2.5%    50%  97.5% overlap0    f Rhat n.eff
# sigma  152.52 14.49 129.79 150.74 186.38    FALSE 1.00 1.00  3239
# beta0   -6.17  0.18  -6.52  -6.17  -5.80    FALSE 1.00 1.00   587
# beta1    1.66  0.11   1.44   1.66   1.89    FALSE 1.00 1.00   658
# beta2    0.33  0.04   0.25   0.34   0.39    FALSE 1.00 1.00  1442
# beta3   -0.27  0.11  -0.50  -0.28  -0.06    FALSE 0.99 1.00  2446
# beta4   -0.38  0.13  -0.64  -0.37  -0.14    FALSE 1.00 1.00  1115
# alpha0   0.08  0.30  -0.50   0.07   0.67     TRUE 0.59 1.00  1658
# alpha1   0.30  0.14   0.03   0.30   0.59    FALSE 0.99 1.00  4000
# alpha2   0.29  0.15   0.01   0.29   0.61    FALSE 0.98 1.00  4000
# alpha3   0.07  0.15  -0.23   0.07   0.39     TRUE 0.69 1.00  3212
# alpha4   0.86  0.22   0.49   0.84   1.33    FALSE 1.00 1.00  1414
# Mtotal 228.14 20.93 195.00 226.00 276.00    FALSE 1.00 1.01  1885


# 11.8.5 Summary of the TEDS/TPP model (no code)
