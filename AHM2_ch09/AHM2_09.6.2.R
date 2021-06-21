#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 65 mins
# Run time with the full number of iterations: 7.6 hrs

library(AHMbook)
library(jagsUI)

# 9.6 Mechanistic, or dynamic, models of spatial autocorrelation
# ==============================================================

# 9.6.2 Continuous representation of space in dynamic occupancy models
# --------------------------------------------------------------------

# Get the data
library(AHMbook)
data(EurasianLynx)
str(datfull <- EurasianLynx)
selection <- datfull$type == 'certain' & datfull$Cntry == 'Switzerland' &
    datfull$xcoord < 4185
str(dat <- datfull[selection,])      # Swiss, certain and west of x = 4185
year <- 1994:2016

# Format detection/nondetection data in a 3D array
( nsites <- length(unique(dat$site.nr)) )         # 160
( nyears <- length(unique(dat$Year)) )            # 23
( nsurveys <- 3)
ylong <- as.matrix(dat[,3:5])
head(dat)              # site changes fastest so must be 1st dimension
tmp <- array(ylong, c(nsites, nyears, nsurveys))
y <- aperm(tmp, c(1, 3, 2))

# Compute observed lynx occurence in every year
zobs <- apply(y, c(1,3), max, na.rm = TRUE)

# Get coordinate info for the Swiss lynx (kilometres)
head( grid <- cbind(dat$xcoord[1:nsites], dat$ycoord[1:nsites]) )

# Get forest for restricted set of 160 sites
oforest <- dat$forest[1:nsites]
forest <- standardize(oforest)

# Get distance matrix
d <- e2dist(grid/100, grid/100)              # In units of 100 km
hist(d, main = 'Distribution of inter-quadrat distance (in 100 km units)') # not shown
distSq <- d^2 # Need distance squared in half-normal kernel

# Bundle and summarize data set
str(bdata <- list(y = y, nsites = nsites, nsurveys = nsurveys, nyears = nyears,
    forest = forest, distSq = distSq))
# List of 6
# $ y       : int [1:160, 1:3, 1:23] NA NA NA NA NA NA NA NA NA NA ...
# $ nsites  : int 160
# $ nsurveys: num 3
# $ nyears  : int 23
# $ forest  : num [1:160] -0.746 -0.646 -0.556 0.416 -1.025 ...
# $ distSq  : num [1:160, 1:160] 0 0.01 0.01 0.02 0.05 ...

# Specify model in BUGS language
cat(file = "chandler.txt","
model {

  # Priors
  psi1.int ~ dunif(0, 1) # Initial occupancy
  alpha.lpsi1 <- logit(psi1.int)
  beta.lpsi1.forest ~ dnorm(0, 0.1)
  phi.int ~ dunif(0, 1) # Persistence
  alpha.lphi <- logit(phi.int)
  beta.lphi.forest ~ dnorm(0, 0.1)
  gamma.int ~ dunif(0, 1) # Colonization
  alpha.lgamma <- logit(gamma.int)
  beta.lgamma.forest ~ dnorm(0, 0.1)
  sigma ~ dgamma(1, 0.1)
  p.int ~ dunif(0, 1) # Detection
  alpha.lp <- logit(p.int)
  beta.lp.forest ~ dnorm(0, 0.1)

  # Likelihood
  # Ecological submodel
  for (i in 1:nsites){
  z[i,1] ~ dbern(psi1[i])
    logit(psi1[i]) <- alpha.lpsi1 + beta.lpsi1.forest * forest[i]
    logit(phi[i]) <- alpha.lphi + beta.lphi.forest * forest[i]
    logit(gamma0[i]) <- alpha.lgamma + beta.lgamma.forest * forest[i]
    for(t in 2:nyears) {
      for(n in 1:nsites) {
        # Pairwise colonization probability
        gammaDistPairs[i,n,t-1] <- 1 - gamma0[i] *
        exp( -distSq[i,n] / (2*sigma^2)) * z[n,t-1]
      }
      # Colonization probability
      gamma[i,t-1] <- 1 - prod(gammaDistPairs[i,1:nsites,t-1])
      ## Pr(z=1|z=1) = 1 - Pr(extinction)*Pr(not colonized) = 1 - ((1-phi)*(1-gamma))
      muz[i,t-1] <- gamma[i,t-1]*(1-z[i,t-1]) + (1-(1-phi[i])*(1-gamma[i,t-1]))* z[i,t-1]
      z[i,t] ~ dbern(muz[i,t-1])
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
zst <- array(1, dim = c(nsites, nyears))
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi1.int", "alpha.lpsi1", "beta.lpsi1.forest", "phi.int",
    "alpha.lphi", "beta.lphi.forest", "gamma.int", "alpha.lgamma",
    "beta.lgamma.forest", "sigma", "p.int", "alpha.lp", "beta.lp.forest",
    "z", "gamma")

# MCMC settings
# na <- 1000 ; ni <- 4000 ; nt <- 2 ; nb <- 2000 ; nc <- 3
na <- 100 ; ni <- 400 ; nt <- 1 ; nb <- 200 ; nc <- 3  # ~~~ for testing

# Call JAGS (ART 485 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "chandler.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out3)
print(out3$summary[1:13, c(1:3,7)], 3)
#                       mean     sd    2.5%  97.5%
# psi1.int            0.1189 0.0367  0.0543  0.199
# alpha.lpsi1        -2.0529 0.3706 -2.8576 -1.394
# beta.lpsi1.forest   1.0608 0.3533  0.4211  1.833
# phi.int             0.5118 0.0728  0.3634  0.646
# alpha.lphi          0.0479 0.2972 -0.5606  0.602
# beta.lphi.forest   -0.0978 0.2589 -0.6195  0.379
# gamma.int           0.1423 0.0374  0.0842  0.230
# alpha.lgamma       -1.8286 0.3015 -2.3868 -1.210
# beta.lgamma.forest  0.6396 0.0950  0.4636  0.828
# sigma               0.0988 0.0092  0.0819  0.118
# p.int               0.3202 0.0163  0.2894  0.352
# alpha.lp           -0.7541 0.0750 -0.8982 -0.609
# beta.lp.forest      0.6360 0.0649  0.5051  0.761


# ~~~ extra code for figure 9.19 ~~~~~~~~
# Compute posterior mean and CRI of distance effect on gamma
tmp <- out3$sims.list
nsims <- out3$mcmc.info$n.samples
pred <- array(NA, dim = c(1000, nsims))
preddist <- seq(0, 0.4, length.out = 1000)
preddistSq <- preddist^2
for(i in 1:nsims){
  pred[,i] <- tmp$gamma.int[i] * exp( -preddistSq / (2*tmp$sigma[i]^2))
}
pm <- apply(pred, 1, mean)
CRI <- apply(pred, 1, function(x) (quantile(x, prob = c(0.025, 0.975))))

# Plot
plot(100*preddist, pm, type = 'l', lty = 1, lwd = 3, col = 'blue',
    xlab = 'Intersite distance (km)', ylab = ' Colonizationprobability',
    ylim = c(0, 0.3), frame = FALSE)
polygon(c(100*preddist, rev(100*preddist)), c(CRI[1,], rev(CRI[2,])),
    col = 'grey', border = NA)
lines(100*preddist, pm, type = 'l', lty = 1, lwd = 3, col = 'blue')

# ~~~~ extra code for figure 9.20 ~~~~~~~~~~~
# Plot estimated colonization and SDM in first and last years (i.e., 1994 and 2016)
library(raster)
op <- par(mfrow = c(2,2), mar = c(1,1,4,6))
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
r <- rasterFromXYZ(data.frame(x = grid[,1], y = grid[,2], z = out3$mean$gamma[,1]))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, zlim = c(0, 1),
    main = 'Colonization in winter 1994/1995', legend = FALSE)
r <- rasterFromXYZ(data.frame(x = grid[,1], y = grid[,2], z = out3$mean$gamma[,22]))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, zlim = c(0, 1),
    main = 'Colonization in winter 2016/2017')

mapPalette <- colorRampPalette(c("grey", "lightgreen", "darkgreen"))
r <- rasterFromXYZ(data.frame(x = grid[,1], y = grid[,2], z = oforest))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, zlim = c(0, 100),
    main = 'Conditional occupancy in winter 1994/1995', legend = FALSE)
points(grid[,1], y = grid[,2], pch = 15, col = rgb(0,0,0, out3$mean$z[,1]), cex = 1)
r <- rasterFromXYZ(data.frame(x = grid[,1], y = grid[,2], z = oforest))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, zlim = c(0, 100),
    main = 'Conditional occupancy in winter 2016/2017')
points(grid[,1], y = grid[,2], pch = 15, col = rgb(0,0,0, out3$mean$z[,23]), cex = 1)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
