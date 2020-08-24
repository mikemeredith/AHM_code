#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 50 mins
# Run time with the full number of iterations: 1.5 hrs

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14" # the location of the WinBUGS14.exe file on your machine


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
oelev<- ch$elev[belongs.to] # 'o' means 'original' or unscaled
oforest<- ch$forest[belongs.to]

# Standardize elevation and forest cover
elev<- standardize(oelev)
forest<- standardize(oforest)

# Also grab counts, survey dates and survey durations
head(counts <- as.matrix(wp[,7:48]))   # Counts
head(dates<- as.matrix(wp[,49:90]))   # Survey dates
head(durs <- as.matrix(wp[,91:132]))   # Survey durations

# Put these three into 3d arrays
nsite<- 267
nrep<- 3
nyear<- 14
C<- array(counts, dim = c(nsite, nrep, nyear))
ODATE<- array(dates, dim = c(nsite, nrep, nyear))
ODUR<- array(durs, dim = c(nsite, nrep, nyear))

# Standardize and mean-impute date and duration
DATE <- standardize(ODATE)
DATE[is.na(DATE)] <- 0          # mean-impute
DUR <- standardize(ODUR)
DUR[is.na(DUR)] <- 0            # mean-impute

# Generate 5 x 5 km2 blocks and add new columns to the 'ch' data frame with the x and y coordinates and an ID number for the block within which each of the 267 MHB survey quadrats lies.

# Round all Swiss quadrats by 5 km
block.side <- 5000   # Length of side of block for spatial aggregation
ch$xblock<- block.side *round(ch$x/block.side)
ch$yblock<- block.side *round(ch$y/block.side)
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
 # $ adj    : int [1:13930] 2 4 5 1 3 4 5 6 2 5 ...
 # $ weights: num [1:13930] 1 1 1 1 1 1 1 1 1 1 ...
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
str(bdata <- list(C = C, elev = elev, forest = forest, DATE = DATE, DUR = DUR,
    nsite = nsite, nsurvey = 3, nyear = nyear, n.block = n.block,
    MHBblockID = MHBblockID, adj = winnb$adj, weights = winnb$weights,
    num = winnb$num))
# List of 13
# $ C         : int [1:267, 1:3, 1:14] 0 3 0 0 0 0 0 0 0 0 ...
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

# Specify model in BUGS language
cat(file = "SVC.txt","
model {

  # Priors
  for(t in 1:nyear){
    alpha0[t] <- logit(p0[t])
    p0[t] ~ dunif(0, 1) # Intercept of detection
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

  # Linear model for the trend, with its priors
  for(i in 1:nsite){
    trend[i] <- mean.trend + eta[MHBblockID[i]]
  }
  mean.trend ~ dnorm(0, 1)

  # CAR prior distribution for spatial random effects in the trend
  eta[1:n.block] ~ car.normal(adj[], weights[], num[], tau)
  tau <- 1/v.eta
  sd.eta <- sqrt(v.eta)
  v.eta ~ dnorm(0, 0.01)I(0,)

  # Likelihood
  for (i in 1:nsite){
    for(t in 1:nyear){
      # State process
      N[i,t] ~ dpois(lambda[i, t])
      log(lambda[i,t]) <- min(10, max(-10, loglam.lim[i,t]))
      loglam.lim[i,t] <- beta0 + beta1[1] * elev[i] +
      beta1[2] * pow(elev[i],2) + beta1[3] * pow(elev[i],3) +
      beta1[4] * forest[i] + trend[i] * (t-7.5)
      for(j in 1:nsurvey){
        # Observation process
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t])
        logit(p[i,j,t]) <- lp.lim[i,j,t]
        lp.lim[i,j,t] <- min(500, max(-500, lp[i,j,t]))
        lp[i,j,t] <- alpha0[t] + alpha1[1] * DATE[i,j,t] +
        alpha1[2] * pow(DATE[i,j,t],2) + alpha1[3] * DUR[i,j,t]
      }
    }
  }
  # Derived quantities
  for(t in 1:nyear){
    Ntotal[t] <- sum(N[,t])
    meanPopLevel[t] <- exp(beta0 + mean.trend * (t-7.5))
  }
}
")

# Initial values
tmp <- apply(C, c(1,3), max, na.rm = TRUE) + 1
tmp[tmp == '-Inf'] <- 1
inits <- function(){list(N = tmp, beta0 = 1, v.eta = abs(rnorm(1)),
    eta = rep(0, n.block))}

# Parameters monitored
params <- c("lam0", "beta1", "mean.trend", "p0", "alpha1", "Ntotal",
    "meanPopLevel", "trend", "v.eta", "sd.eta", "eta")

# MCMC settings
# ni <- 12000 ; nt <- 6 ; nb <- 6000 ; nc <- 3
ni <- 6000 ; nt <- 1 ; nb <- 5000 ; nc <- 3  # ~~~ for testing, 42 mins

# Call WinBUGS (ART 100 min), check convergence and summarize posteriors
library(R2WinBUGS)
# out <- bugs(bdata, inits, params, "SVC.txt", n.chains = nc, n.thin = nt,
    # n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir,
    # DIC = FALSE, working.directory = getwd()) # Must close program by hand
# ~~~~ for testing, use the following call:
out <- bugs(bdata, inits, params, "SVC.txt", n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir,
    DIC = FALSE)
print(out$summary[1:20,], 3)

# ~~~~~ extra code for figure 9.13 ~~~~~~~~~~~~~
# Look at map of trends
block.trend <- exp(out$mean$mean.trend+out$mean$eta) # 1 is stable
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
