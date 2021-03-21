#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================
# Code put together by Mike Meredith
# Based on code from proofs dated 2020-08-19
# NIMBLE adaptation by Mike Meredith with input from Chris Paciorek

# Approximate execution time for this script: 13 mins
# Run time with the full number of iterations: 37 mins

##### USING NIMBLE INSTEAD OF WinBUGS  #####

library(AHMbook)
library(unmarked)         # Contains the Swiss landscape data set
library(nimble)
library(mcmcOutput)

# 9.5 Fitting a simple SVC model: spatial modeling of population trend
# ====================================================================

# ~~~~ data wrangling code inserted from MS dated 2019-07-02
# Read data set
library(AHMbook)          # Contains the Green woodpecker data set
library(unmarked)         # Contains the Swiss landscape data set
data(greenWoodpecker)
str(wp <- greenWoodpecker)
data(Switzerland)
str(ch <- Switzerland)

# Map of the modelled domain (Switzerland) and of the MHB site locations (not shown)
plot(ch$x, ch$y, pch = '.', col = 'grey', asp = 1, frame = FALSE)
points(wp$x, wp$y, pch = 16)  # The 267 MHB quadrats

# Create unique quadrat identifiers based on x and y coordinates
wp_quadid <- paste(wp$x, wp$y, sep = '.') # For the 267 survey quadrats
ch_quadid <- paste(ch$x, ch$y, sep = '.') # For all of Switzerland

# Assign each one of the 267 MHB survey quadrats to the large sample of 42k quads in the entire modeled domain (Switzerland)
belongs.to<- pmatch(wp_quadid, ch_quadid)
sum(is.na(belongs.to)) == 0 # TRUE, ie, all wp quads match 1 ch quad
all(diff(belongs.to) > 0) # TRUE, ie, quads are in same order in wp and ch, makes life easier!
ch$MHBquad <- NA
ch$MHBquad[belongs.to] <- 1:length(belongs.to) # each ch$MHBquad has MHB quadrat ID or NA.
summary(ch)

# Grab elevation and forest cover covariates from Swiss landscape file
oelev <- ch$elev[belongs.to] # 'o' means 'original' or unscaled
oforest <- ch$forest[belongs.to]

# Standardize elevation and forest cover
elev <- standardize(oelev)
forest <- standardize(oforest)

# Also grab counts, survey dates and survey durations
head(counts <- as.matrix(wp[,7:48]))   # Counts
head(dates <- as.matrix(wp[,49:90]))   # Survey dates
head(durs <- as.matrix(wp[,91:132]))   # Survey durations

# Put these three into 3d arrays
nsite <- 267
nrep <- 3
nyear <- 14
C <- array(counts, dim = c(nsite, nrep, nyear))
ODATE <- array(dates, dim = c(nsite, nrep, nyear))
ODUR <- array(durs, dim = c(nsite, nrep, nyear))

# Standardize and mean-impute date and duration
DATE <- standardize(ODATE)
DATE[is.na(DATE)] <- 0          # mean-impute
DUR <- standardize(ODUR)
DUR[is.na(DUR)] <- 0            # mean-impute

# Generate 5 x 5 km2 blocks and add new columns to the 'ch' data frame with the x and y coordinates and an ID number for the block within which each of the 267 MHB survey quadrats lies.

# Round all Swiss quadrats by 5 km
block.side <- 5000   # Length of side of block for spatial aggregation
ch$xblock <- block.side * round(ch$x/block.side)
ch$yblock <- block.side * round(ch$y/block.side)
plot(ch$xblock, ch$yblock, pch = 0, asp = 1,
    main = "Switzerland by 5 x 5 km2 blocks, with location of the 267 MHB quads",
    frame = FALSE)
points(wp$x, wp$y, col = 'red', pch = 16)# Map not shown


# Create block identifiers from the x and y coords
block_x.y <- paste(ch$xblock, ch$yblock, sep='.')
unam <- unique(block_x.y) # unique block identifiers
(n.block <- length(unam)) # 1842

# Add block numbers to the 'ch' data frame
ch$blockNr <- as.numeric(factor(block_x.y, levels=unam))
head(ch, 20)

# Now pull out the x and y coords of the blocks
blockX <- tapply(ch$xblock, ch$blockNr, mean)
blockY <- tapply(ch$yblock, ch$blockNr, mean)
plot(blockX, blockY, asp = 1)  # Plot again to check it is right

# Compute the neighborhood (2nd order, queen) for CAR modeling
# all blocks with a distance between d1 and d2 comprise the neighbourhood
# for a cell (choice of d2 uses Pythagoras)
library(spdep)
blockcoordgrid <- cbind(blockX, blockY)
neigh <- dnearneigh(blockcoordgrid, d1 = 0, d2 = sqrt(2) * block.side+ 0.1)
winnb <- nb2WB(neigh)  # Function to get CAR ingredients for BUGS

# Summarize the neighborhood information needed by WinBUGS
str(winnb)
# List of 3
# $ adj     : int [1:13930] 2 4 5 1 3 4 5 6 2 5 ...
# $ weights : num [1:13930] 1 1 1 1 1 1 1 1 1 1 ...
# $ num     : int [1:1842] 3 5 3 5 8 5 5 8 6 4 ...

# Frequency distribution of the number of neighbours
table(winnb$num)
   # 2    3    4    5    6    7    8
   # 1   19   56   77   70  110 1509

# Membership indicator: tells us in which block each MHB quadrat lies
MHBblockID <- ch$blockNr[!is.na(ch$MHBquad)]

# Quadrats are in the same order in 'ch' and 'wp' data frames.
sum(duplicated(MHBblockID))  # 1 blockID occurs twice, ie,
#  one block contains 2 MHB quadrats.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bundle and summarize data set
bdata <- list(C = C)
str(bconst <- list(elev = elev, forest = forest, DATE = DATE, DUR = DUR,
    nsite = nsite, nsurvey = 3, nyear = nyear, n.block = n.block,
    MHBblockID = MHBblockID, adj = winnb$adj, weights = winnb$weights,
    num = winnb$num))
# List of 12
# $ elev      : num [1:267] -1.1539 -1.1539 -0.2175 -0.3735 -0.0614 ...
# $ forest    : num [1:267] -1.1661 -0.4426 -0.1379 -0.9376 -0.0618 ...
# $ DATE      : num [1:267, 1:3, 1:14] -1.09 -1.32 -1.23 -1.27 -1.36 ...
# $ DUR       : num [1:267, 1:3, 1:14] 0.128 -0.978 -0.297 -0.893 ...
# $ nsite     : num 267
# $ nsurvey   : num 3
# $ nyear     : num 14
# $ n.block   : int 1842
# $ MHBblockID: num [1:267] 15 29 34 43 46 61 64 87 98 101 ...
# $ adj       : int [1:13930] 2 4 5 1 3 4 5 6 2 5 ...
# $ weights   : num [1:13930] 1 1 1 1 1 1 1 1 1 1 ...
# $ num       : int [1:1842] 3 5 3 5 8 5 5 8 6 4 ...

# Specify model in NIMBLE dialect of BUGS language
# ''''''''''''''''''''''''''''''''''''''''''''''''
# Necessary changes:
# 1. Replace I(0,) with T(..., 0, )
# 2. dcar_normal output is not zero-centered (unlike WinBUGS version), so should
#    not be used with an intercept: `mean.trend` removed.
#    (`trend` retained for consistency with the text, p.558)
# Optional improvements:
# 1. Use logit(.) <-
# 2. dbeta for probabilities instead of dunif
# 3. Use ^ notation instead of pow()

SVC.code <- nimbleCode({
  # Priors
  for(t in 1:nyear){
    alpha0[t] <- logit(p0[t])
    p0[t] ~ dbeta(1, 1) # Intercept of detection
  }
  for(v in 1:4){ # Priors for coefficients in lambda model
    beta1[v] ~ dnorm(0, 0.1)
  }
  # curve(dnorm(x, 0, sqrt(1/1)), -3, 3) # Look at prior (when executed in R)
  for(w in 1:3){ # Priors for coefficients in p model
    alpha1[w] ~ dnorm(0, 0.1)
  }
  beta0 ~ dnorm(0, 0.1)
  lam0 <- exp(beta0)

  # Linear model for the trend
  for(i in 1:nsite){
    trend[i] <- eta[MHBblockID[i]]
  }

  # CAR prior distribution for spatial random effects in the trend
  eta[1:n.block] ~ dcar_normal(adj[], weights[], num[], tau)
  tau <- 1/v.eta
  sd.eta <- sqrt(v.eta)
  v.eta ~ T(dnorm(0, 0.01),0, )

  # Likelihood
  for (i in 1:nsite){
    for(t in 1:nyear){
      # State process
      N[i,t] ~ dpois(lambda[i, t])
      log(lambda[i,t]) <- min(10, max(-10, loglam.lim[i,t]))
      loglam.lim[i,t] <- beta0 + beta1[1] * elev[i] +
          beta1[2] * elev[i]^2 + beta1[3] * elev[i]^3 +
          beta1[4] * forest[i] + trend[i] * (t-7.5)
      for(j in 1:nsurvey){
        # Observation process
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t])
        logit(p[i,j,t]) <- lp.lim[i,j,t]
        lp.lim[i,j,t] <- min(500, max(-500, lp[i,j,t]))
        lp[i,j,t] <- alpha0[t] + alpha1[1] * DATE[i,j,t] +
            alpha1[2] * DATE[i,j,t]^2 + alpha1[3] * DUR[i,j,t]
      }
    }
  }
  # Derived quantities
  for(t in 1:nyear){
    Ntotal[t] <- sum(N[,t])
    meanPopLevel[t] <- exp(beta0 + mean(trend[1:nsite]) * (t-7.5))
  }
} )


# Initial values
Nst <- apply(C, c(1,3), max, na.rm = TRUE) + 1
Nst[Nst == -Inf] <- 1
# Provide initial values too for the missing values in 'C'
# (not essential, but avoids a bunch of ugly warnings).
Cst <- C
Cst[is.na(C)] <- 1
Cst[!is.na(C)] <- NA

inits <- function(){list(N = Nst, C = Cst, eta = rep(0, n.block), v.eta = abs(rnorm(1)),
    # additional starting values added to fix issues with non-converging chains:
    p0 = runif(14, 0.1, 0.4), beta1 = runif(4, -0.5, 0.5), alpha1 = runif(3, -0.5, 0.5),
    beta0 = runif(1, 0.1, 0.5))}

# Parameters monitored
params <- c("beta0", "v.eta", "alpha1", "beta1", "p0", # <-- top-level nodes
    "lam0", "Ntotal", "meanPopLevel", "trend", "sd.eta", "eta")

# Run 2 chains in series to check for issues
# ''''''''''''''''''''''''''''''''''''''''''
# MCMC settings
ni <- 2000 ; nt <- 1 ; nb <- 1200 ; nc <- 2  # ~~~ for testing, 7 mins

system.time(
out <- nimbleMCMC(SVC.code, data=bdata, constants=bconst,
    inits=inits, monitors=params,
    niter=ni, nburnin=nb, thin=nt, nchains=nc,
    samplesAsCodaMCMC=TRUE) )
( mco <- mcmcOutput(out) )
View(summary(mco, n.eff=TRUE))
# check top-level nodes
diagPlot(mco, c("alpha", "beta", "v.eta", "p0"))

# Use 'foreach' to do a long run in parallel
# ''''''''''''''''''''''''''''''''''''''''''
library(foreach)
library(parallel)
library(doParallel)

detectCores()  # number of cores on your machine
# ncore <- 7   # <-- adjust this for your task
ncore <- 3     # ~~~ testing

cl <- makeCluster(ncore)
registerDoParallel(cl)

# ni <- 22000 ; nt <- 5 ; nb <- 2000  # 30  mins, enough with 3 chains
ni <- 2000 ; nt <- 1 ; nb <- 1200   # ~~~ for testing, 4.5 mins

seeds <- 1:ncore

system.time(
res <- foreach(x = seeds, .packages="nimble",
      .errorhandling='remove', .inorder=FALSE) %dopar% {
  set.seed(x)
  nimbleMCMC(SVC.code, data=bdata, constants=bconst,
      inits=inits, monitors=params,
      niter=ni, nburnin=nb, thin=nt, nchains=1,
      samplesAsCodaMCMC=TRUE)
} )
stopCluster(cl)

# Convert to an mcmcOutput object and look at diagnostics
mclist <- coda::mcmc.list(res)
( mco <- mcmcOutput(mclist) )
diagPlot(mco, c("alpha", "beta", "mean.trend", "v.eta", "p0"))
View(summary(mco))

# ~~~~~ extra code for figure 9.13 ~~~~~~~~~~~~~
out <- sumryList(mco)
# Look at map of trends
block.trend <- exp(out$mean$eta) # 1 is stable
summary(block.trend)

# Map trend at 5x5km2 scale
library(raster)
op <- par(mfrow = c(1, 1), mar = c(1,2,4,4), cex.main = 1.5)
mapPalette<- colorRampPalette(c("grey", "yellow", "orange", "red"))
r <- rasterFromXYZ(data.frame(x = blockX, y = blockY, z = block.trend))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, main = "",
    zlim = c(0.7, 1.5))
par(op)

# ~~~~~ extra code for figure 9.14 ~~~~~~~~~~~~~
# Distinguish areas with positive and with negative trend
tx <- block.trend
tx[tx >= 1] <- 1
tx[tx < 1] <- -1

op <- par(mfrow = c(1, 1), mar = c(1,2,4,1), cex.main = 1.5)
mapPalette<- colorRampPalette(c("blue", "yellow"))
r <- rasterFromXYZ(data.frame(x = blockX, y = blockY, z = tx))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, main = "",
    zlim = c(-1, 1))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
