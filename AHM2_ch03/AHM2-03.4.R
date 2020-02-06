#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 3 : HIERARCHICAL MODELS OF SURVIVAL
# ===========================================

library(AHMbook)
library(jagsUI)
bugs.dir <- "C:/WinBUGS14"

# 3.4 SPATIAL HIERARCHICAL CJS MODELS
# ===================================

# 3.4.1 BRITISH WILLOW WARBLER DATA
# ---------------------------------

# Load the data set (in library AHMbook)
data(willowWarbler)
str(willowWarbler)
# List of 4
# $ birds :'data.frame': 10551 obs. of 12 variables:
# ..$ 1986 : num [1:10551] 1 1 1 1 1 1 1 1 1 1 ...
# ..$ 1987 : num [1:10551] 0 0 0 1 0 0 0 0 0 0 ...
# [ ... truncated ... ]
# ..$ 1995 : num [1:10551] 0 0 0 0 0 0 0 0 0 0 ...
# ..$ 1996 : num [1:10551] 0 0 0 0 0 0 0 0 0 0 ...
# ..$ cesID: num [1:10551] 1 1 1 1 1 1 1 1 1 1 ...
# $ cells :'data.frame': 9667 obs. of 4 variables:
# ..$ lon : num [1:9667] 167500 172500 167500 172500 177500 ...
# ..$ lat : num [1:9667] 12500 12500 17500 17500 17500 22500 22500 ...
# ..$ gdd : num [1:9667] 2136 2120 1945 1943 1956 ...
# ..$ blockID: num [1:9667] 1 1 3 3 3 2 2 3 3 3 ...
# $ CES :'data.frame': 193 obs. of 4 variables:
# ..$ cesx : num [1:193] 307500 432500 457500 312500 477500 ...
# ..$ cesy : num [1:193] 92500 167500 352500 227500 377500 ...
# ..$ BlockID: num [1:193] 25 77 204 110 222 283 119 234 295 152 ...
# ..$ CellID : num [1:193] 319 1366 4302 2266 4627 ...
# $ blocks:'data.frame': 495 obs. of 2 variables:
# ..$ blockX: num [1:495] 175 150 175 200 150 175 200 225 250 275 ...
# ..$ blockY: num [1:495] 0 25 25 25 50 50 50 50 50 50 ...
attach (willowWarbler)
ch <- as.matrix(birds[, 1:11])
sitevec <- birds$cesID

# Number of captures per bird, year and site
table(apply(ch, 1, sum)) # Table of capture frequency per bird
apply(ch, 2, sum) # Number of birds per year
plot(table(table(sitevec))) # Frequency distribution of number of birds captured per site (not shown)
summary(as.numeric(table(sitevec)))
# 1 2 3 4 5 6
# 8997 1214 279 48 12 1
# [1] 692 795 1022 1423 1372 1231 1283 1331 1575 1516 280
# Min. 1st Qu. Median Mean 3rd Qu. Max.
# 2.00 14.00 30.00 54.67 76.00 338.00

# Map of GDD covariate values and 193 CES locations (Fig. 3.11)
library(raster)
mapPalette <- colorRampPalette(c("gray", "yellow", "orange", "red"))
r1 <- with(cells, rasterFromXYZ(data.frame(x = lon, y = lat, z = gdd)))
plot(r1, col = mapPalette(100), axes = F, box = F,
    main ="Map of GDD covariate with 193 CES locations")
with(CES, points(cesx, cesy, pch=16, col='blue', cex = 0.8))

# Get sample sizes
(nyear <- ncol(ch)) # Number years: 11
(nsite <- nrow(CES)) # Number of CE sites: 193
(nblock <- nrow(blocks)) # Number of blocks: 495

(marr <- ch2marray(ch))
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# [1,] 55 13 1 1 0 0 0 0 0 0 622
# [2,] 0 123 17 2 1 0 0 0 0 0 652
# [3,] 0 0 143 20 4 0 0 0 0 0 855
# [4,] 0 0 0 179 16 6 0 0 0 0 1222
# [5,] 0 0 0 0 183 14 1 0 0 0 1174
# [6,] 0 0 0 0 0 161 16 5 0 0 1049
# [7,] 0 0 0 0 0 0 172 20 4 0 1087
# [8,] 0 0 0 0 0 0 0 239 41 9 1042
# [9,] 0 0 0 0 0 0 0 0 252 24 1299
# [10,] 0 0 0 0 0 0 0 0 0 247 1269

# Calculate the number of birds released each year 1–10
(r <- apply(marr, 1, sum))
# [1] 692 795 1022 1423 1372 1231 1283 1331 1575 1516

# Create 3d (or multi-site) m-array and r array: MARR and R
MARR <- array(NA, dim = c(10, nyear, nsite))
R <- array(NA, dim = c(10, nsite))
for(k in 1:nsite){
  sel.part <- ch[sitevec == k, ]
  ma <- ch2marray(sel.part)
  MARR[,,k] <- ma
  R[,k] <- apply(ma, 1, sum)
}
MARR ; R # Look at them and make sure you understand them

# 3.4.2 A HIERARCHICAL CJS MODEL WITH RANDOM SITE AND RANDOM YEAR EFFECTS
# -----------------------------------------------------------------------

# Bundle and summarize data set
str(bdata <- list(MARR = MARR, R = R, n.site = nsite, n.occ = nyear))
# List of 4
# $ MARR : num [1:10, 1:11, 1:193] 1 0 0 0 0 0 0 0 0 0 ...
# $ R : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site : int 193
# $ n.occ : int 11
# Specify model in BUGS language
cat(file = "cjs6.txt","
model {
  # Priors and linear models
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      phi[t, s] <- ilogit(lphi[t, s]) # survival
      p[t, s] <- ilogit(lp[t, s]) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    # Define random site effects
    alpha.lphi.site[s] ~ dnorm(mu.lphi, tau.lphi.site)
    alpha.lp.site[s] ~ dnorm(mu.lp, tau.lp.site)
    mean.phi.site[s] <- ilogit(alpha.lphi.site[s])
    mean.p.site[s] <- ilogit(alpha.lp.site[s])
  }
  # Define random year effects
  for (t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, tau.lphi.time)
    beta.lp.time[t] ~ dnorm(0, tau.lp.time)
    mean.phi.time[t] <- ilogit(mu.lphi + beta.lphi.time[t])
    mean.p.time[t] <- ilogit(mu.lp + beta.lp.time[t])
  }
  # Hyperpriors for hyperparams
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  tau.lphi.site <- pow(sd.lphi.site, -2)
  sd.lphi.site ~ dunif(0, 3)
  tau.lp.site <- pow(sd.lp.site, -2)
  sd.lp.site ~ dunif(0, 3)
  tau.lphi.time <- pow(sd.lphi.time, -2)
  sd.lphi.time ~ dunif(0, 3)
  tau.lp.time <- pow(sd.lp.time, -2)
  sd.lp.time ~ dunif(0, 3)
  # Multinomial likelihood for the m-array data (JAGS style)
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      MARR[t,1:n.occ,s] ~ dmulti(pr[t, , s], R[t,s])
    }
  }
  # Define the cell probabilities of the m-array
  # Main diagonal
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      q[t,s] <- 1-p[t,s] # Probability of non-recapture
      pr[t,t,s] <- phi[t,s]*p[t,s]
      # Above main diagonal
      for (j in (t+1):(n.occ-1)){
        pr[t,j,s] <- prod(phi[t:j,s])*prod(q[t:(j-1),s])*p[j,s]
      } #j
      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,j,s] <- 0
      } #j
    } #t
  } #s
  # Last column: probability of non-recapture
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      pr[t,n.occ,s] <- 1-sum(pr[t,1:(n.occ-1),s])
    } #t
  }#s
}
")

# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1))}
# Parameters monitored
params <- c("mean.phi", "mean.p", "sd.lphi.site", "sd.lp.site", "sd.lphi.time",
    "sd.lp.time", "mean.phi.site", "mean.p.site", "mean.phi.time", "mean.p.time")
# MCMC settings# MCMC settings
na <- 1000 ; ni <- 30000 ; nt <- 10 ; nb <- 20000 ; nc <- 3
# Call JAGS (ART 192 min), check convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "cjs6.txt", n.adapt = na, n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2,3)) ; traceplot(out6)
print(out6, 3)

# mean sd 2.5% 50% 97.5% overlap0 f Rhat n.eff
# mean.phi 0.281 0.020 0.243 0.281 0.321 FALSE 1 1.031 104
# mean.p 0.380 0.038 0.310 0.379 0.458 FALSE 1 1.002 779
# sd.lphi.site 0.158 0.078 0.031 0.153 0.329 FALSE 1 1.007 279
# sd.lp.site 1.138 0.135 0.888 1.135 1.413 FALSE 1 1.001 3000
# sd.lphi.time 0.223 0.086 0.100 0.208 0.431 FALSE 1 1.006 1007
# sd.lp.time 0.149 0.108 0.005 0.130 0.389 FALSE 1 1.014 226
# mean.phi.site[1] 0.277 0.040 0.187 0.280 0.354 FALSE 1 1.016 175
# mean.phi.site[2] 0.274 0.034 0.205 0.275 0.335 FALSE 1 1.006 356
# [ ... output truncated ... ]
# Visualizations (Fig. 3.12)
par(mfrow = c(2, 1), mar = c(3,5,2,1))
plot(1:193, out6$mean$mean.phi.site, xlab = 'Site', ylab = 'Apparent survival', frame = F,
    pch = 16, ylim = c(0.14, 0.45))
segments(1:193, out6$q2.5$mean.phi.site, 1:193, out6$q97.5$mean.phi.site)
abline(h = out6$mean$mean.phi, lty = 1, lwd = 2)
plot(1:193, out6$mean$mean.p.site, xlab = 'Site', ylab = 'Recapture', frame = F, pch = 16,
    ylim = c(0, 1))
segments(1:193, out6$q2.5$mean.p.site, 1:193, out6$q97.5$mean.p.site)
abline(h = out6$mean$mean.p, lty = 1, lwd = 2)

# 3.4.3 ADDING SPATIAL COVARIATES INTO THE HIERARCHICAL CJS MODEL
# ---------------------------------------------------------------

# Scale linear and squared GDD separately and latitude (for entire grid)
scaled.gdd1 <- standardize(cells$gdd)  #################
scaled.gdd2 <- standardize(cells$gdd^2)
scaled.lat <- standardize(cells$lat)
# Pull out GDD and LATITUDE values for 193 sites (using CellID)
gdd1.site <- scaled.gdd1[CES$CellID]
gdd2.site <- scaled.gdd2[CES$CellID]
lat.site <- scaled.lat[CES$CellID]
# Bundle and summarize data set
str(bdata <- list(MARR = MARR, R = R, n.site = nsite, n.occ = nyear, gdd1.site = gdd1.site,
    gdd2.site = gdd2.site, lat.site = lat.site))
    # , n.cell = n.cell, gdd1.grid = scaled.gdd, gdd2.grid = scaled.gdd2, lat.grid = lat # if want to make predictions inside BUGS
# List of 7
# $ MARR : num [1:10, 1:11, 1:193] 1 0 0 0 0 0 0 0 0 0 ...
# $ R : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site : int 193
# $ n.occ : int 11
# $ gdd1.site: num [1:193] 1.13 0.726 0.46 0.538 0.686 ...
# $ gdd2.site: num [1:193] 1.228 0.703 0.38 0.473 0.653 ...
# $ lat.site : num [1:193] -1.397 -1.11 -0.403 -0.881 -0.308 ...

# Specify model in BUGS language
cat(file = "cjs7.txt","
model {
  # Priors and linear models
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      phi[t, s] <- ilogit(lphi[t, s]) # survival
      p[t, s] <- ilogit(lp[t, s]) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    # Linear model for site-level effects: add covariates
    alpha.lphi.site[s] ~ dnorm(mu.lphi.site[s], tau.lphi.site)
    mu.lphi.site[s] <- alpha.mu.lphi + beta1 * gdd1.site[s] + beta2 * gdd2.site[s] + beta3 *
    lat.site[s]
    alpha.lp.site[s] ~ dnorm(mu.lp, tau.lp.site)
    # backtransform site means
    mean.phi.site[s] <- ilogit(alpha.lphi.site[s])
    mean.p.site[s] <- ilogit(alpha.lp.site[s])
  }
  # Define year random effects
  for (t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, tau.lphi.time)
    beta.lp.time[t] ~ dnorm(0, tau.lp.time)
    # backtransform time means
    mean.phi.time[t] <- ilogit(alpha.mu.lphi + beta.lphi.time[t])
    mean.p.time[t] <- ilogit(mu.lp + beta.lp.time[t])
  }
  # Hyperpriors for hyperparams
  alpha.mu.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  tau.lphi.site <- pow(sd.lphi.site, -2)
  sd.lphi.site ~ dunif(0, 3)
  tau.lp.site <- pow(sd.lp.site, -2)
  sd.lp.site ~ dunif(0, 3)
  tau.lphi.time <- pow(sd.lphi.time, -2)
  sd.lphi.time ~ dunif(0, 3)
  tau.lp.time <- pow(sd.lp.time, -2)
  sd.lp.time ~ dunif(0, 3)
  # Coefficients for gdd1, gdd2 and lat
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  beta3 ~ dnorm(0, 0.1)
  # Multinomial likelihood for the m-array data (JAGS style)
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      MARR[t,1:n.occ,s] ~ dmulti(pr[t, , s], R[t,s])
    }
  }
  # Define the cell probabilities of the m-array
  # Main diagonal
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      q[t,s] <- 1-p[t,s] # Probability of non-recapture
      pr[t,t,s] <- phi[t,s]*p[t,s]
      # Above main diagonal
      for (j in (t+1):(n.occ-1)){
        pr[t,j,s] <- prod(phi[t:j,s])*prod(q[t:(j-1),s])*p[j,s]
      } #j
      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,j,s] <- 0
      } #j
    } #t
  } #s
  # Last column of m-array: probability of non-recapture
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      pr[t,n.occ,s] <- 1-sum(pr[t,1:(n.occ-1),s])
    } #t
  }#s
  # Derived quantities: predictions for entire grid
  #for (s in 1:n.cell){
  # phi.grid[s] <- ilogit(alpha.mu.lphi + beta1 * gdd1.grid[s] + beta2 * #gdd2.grid[s] + beta3 * lat.grid[s])
  #}
}
")
# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1), beta1 = rnorm(1),
  beta2 = rnorm(1), beta3 = rnorm(1))}
# Parameters monitored
params <- c("mean.phi", "mean.p", "alpha.mu.lphi", "mu.lp",
    "sd.lphi.site", "sd.lp.site", "sd.lphi.time", "sd.lp.time", "mean.phi.site",
    "mean.p.site", "mean.phi.time", "mean.p.time", "beta1", "beta2", "beta3")
    # , "phi.grid" # for predictions within JAGS
    
# MCMC settings
# na <- 5000 ; ni <- 60000 ; nt <- 30 ; nb <- 30000 ; nc <- 3
na <- 5000 ; ni <- 6000 ; nt <- 3 ; nb <- 3000 ; nc <- 3 # ~~~~~~~~~ for testing
# Call JAGS (ART 330 min), check convergence and summarize posteriors
out7 <- jags(bdata, inits, params, "cjs7.txt", n.adapt = na, n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2,3)) ; traceplot(out7)
print(out7, 3)
# mean sd 2.5% 50% 97.5% overlap0 f Rhat n.eff
# mean.phi 0.260 0.023 0.216 0.259 0.309 FALSE 1.000 1.029 94
# mean.p 0.401 0.040 0.327 0.399 0.484 FALSE 1.000 1.007 401
# alpha.mu.lphi -1.048 0.120 -1.291 -1.052 -0.806 FALSE 1.000 1.028 97
# mu.lp -0.404 0.168 -0.721 -0.408 -0.065 FALSE 0.989 1.006 412
# sd.lphi.site 0.171 0.094 0.005 0.169 0.369 FALSE 1.000 1.029 105
# sd.lp.site 1.065 0.141 0.793 1.064 1.350 FALSE 1.000 1.002 695
# sd.lphi.time 0.223 0.087 0.100 0.210 0.436 FALSE 1.000 1.003 754
# sd.lp.time 0.158 0.117 0.007 0.135 0.439 FALSE 1.000 1.002 1794
# mean.phi.site[1] 0.227 0.039 0.142 0.229 0.303 FALSE 1.000 1.007 979
# mean.phi.site[2] 0.243 0.032 0.178 0.243 0.310 FALSE 1.000 1.008 359
# [ ... output truncated ... ]
# mean.p.time[9] 0.409 0.046 0.326 0.408 0.503 FALSE 1.000 1.003 788
# mean.p.time[10] 0.403 0.050 0.310 0.400 0.508 FALSE 1.000 1.005 448
# beta1 1.043 0.600 -0.113 1.030 2.188 TRUE 0.960 1.021 125
# beta2 -0.735 0.506 -1.683 -0.730 0.248 TRUE 0.931 1.016 153
# beta3 0.320 0.102 0.123 0.323 0.519 FALSE 0.999 1.003 664

ncell <- nrow(cells) 
phi.grid <- array(NA, dim = c(ncell, out7$mcmc.info$n.samples))
# Derived quantities: predictions for entire grid
sims <- out7$sims.list # Grab the posterior simulations
for (s in 1:ncell){
  phi.grid[s,] <- plogis(sims$alpha.mu.lphi + sims$beta1 * scaled.gdd1[s] + sims$beta2 *
      scaled.gdd2[s] + sims$beta3 * scaled.lat[s])
}
post.mean <- apply(phi.grid, 1, mean) # Posterior mean
post.sd <- apply(phi.grid, 1, sd) # Posterior standard deviation

# 3.4.4 FITTING A SPATIAL HIERARCHICAL CJS MODEL WITH SPATIAL AUTOCORRELATION ONLY
# --------------------------------------------------------------------------------

# Blocks with a distance between d1 and d2 comprise the neighborhood
library(spdep)
blockcoordgrid <- cbind(as.matrix(blocks))
neigh <- dnearneigh(blockcoordgrid, d1 = 0, d2 = sqrt(2 * 25^2)+ 0.1)
str(winnb <- nb2WB(neigh)) # Function to get CAR ingredients for BUGS
# List of 3
# $ adj : int [1:3458] 2 3 4 1 3 5 6 1 2 4 ... # ID of neighbors
# $ weights: num [1:3458] 1 1 1 1 1 1 1 1 1 1 ... # Weights: here, equal
# $ num : int [1:495] 3 4 6 5 3 6 6 6 5 5 ... # Number of neighbors

# Frequency distribution of the number of neighbors
table(winnb$num) # Every block is connected to at least two neighbors

# Reformat the m-array for WinBUGS
dim(MARR) # The nyear = 11 dimension (now #2) must come last for WinBUGS
dim(MARRWB <- aperm (MARR, c(3, 1, 2))) # MARR for WinBUGS
# Bundle and summarize data set for WinBUGS
str(bdata <- list(MARRWB = MARRWB, R = R, n.site = nsite, n.occ = nyear, n.block = nblock,
    BlockID = CES$BlockID, adj = winnb$adj, weights = winnb$weights, num = winnb$num))
# List of 9
# $ MARRWB : num [1:193, 1:10, 1:11] 1 3 0 2 0 0 0 0 0 0 ...
# $ R : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site : int 193
# $ n.occ : int 11
# $ n.block: int 495
# $ BlockID: int [1:193] 25 77 204 110 222 283 119 234 295 152 ...
# $ adj : int [1:3458] 2 3 4 1 3 5 6 1 2 4 ...
# $ weights: num [1:3458] 1 1 1 1 1 1 1 1 1 1 ...
# $ num : int [1:495] 3 4 6 5 3 6 6 6 5 5 ...

# Specify model in BUGS language
cat(file = "cjs8.txt","
model {
  # Priors and linear models
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      phi[t, s] <- 1 / (1 + exp(-lphi[t, s])) # survival
      p[t, s] <- 1 / (1 + exp(-lp[t, s])) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    # rho is spatial effect at the block level
    alpha.lphi.site[s] <- mu.lphi + eta[BlockID[s]]
    alpha.lp.site[s] ~ dnorm(mu.lp, tau.lp.site) I(-12, 12)
    # backtransform site means
    mean.p.site[s] <- 1 / (1 + exp(-alpha.lp.site[s]))
  }
  for (t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, tau.lphi.time) I(-12, 12)
    beta.lp.time[t] ~ dnorm(0, tau.lp.time) I(-12, 12)
    # backtransform time means
    mean.phi.time[t] <- 1 / (1 + exp(-mu.lphi + beta.lphi.time[t]))
    mean.p.time[t] <- 1 / (1 + exp(-mu.lphi + beta.lp.time[t]))
  }
  # Hyperpriors for hyperparams
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  tau.lp.site <- pow(sd.lp.site, -2)
  sd.lp.site ~ dunif(0, 2)
  tau.lphi.time <- pow(sd.lphi.time, -2)
  sd.lphi.time ~ dunif(0, 1)
  tau.lp.time <- pow(sd.lp.time, -2)
  sd.lp.time ~ dunif(0, 1)
  # CAR prior distribution for spatial random effects eta
  # NOTE: this is defined on the entire grid of 495 blocks
  eta[1:n.block] ~ car.normal(adj[], weights[], num[], tau)
  tau ~ dgamma(0.5, 0.0005) # cf. Lunn et al. (2013)
  # curve(1/dgamma(x, 0.5, 0.0005), 0, 10) # howsit look like ?
  veta <- 1/tau
  sdeta <- sqrt(veta)
  # Multinomial likelihood for the m-array data (WinBUGS style)
  # Note ’open index’ in pr[t,s,] comes last
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      MARRWB[s, t,1:n.occ] ~ dmulti(pr[t,s, ], R[t,s])
    }
  }
  # Define the cell probabilities of the m-array
  # Main diagonal
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      q[t,s] <- 1-p[t,s] # Probability of non-recapture
      pr[t,s,t] <- phi[t,s]*p[t,s]
      # Above main diagonal
      for (j in (t+1):(n.occ-1)){
        pr[t,s,j] <- prod(phi[t:j,s])*prod(q[t:(j-1),s])*p[j,s]
      } #j
      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,s,j] <- 0
      } #j
    } #t
  } #s
  # Last column: probability of non-recapture
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      pr[t,s,n.occ] <- 1-sum(pr[t,s,1:(n.occ-1)])
    } #t
  }#s
}
")

# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1), eta = rep(0, nblock))}
# Parameters monitored
params <- c("mean.phi", "mean.p", "mu.lphi", "mu.lp",
    "sd.lp.site", "sd.lphi.time", "sd.lp.time", "mean.p.site", "mean.phi.time", "mean.p.time",
    "veta", "sdeta", "eta")
# MCMC settings
# ni <- 100000 ; nt <- 50 ; nb <- 50000 ; nc <- 3 # 52 hours
ni <- 1000 ; nt <- 5 ; nb <- 500 ; nc <- 3 # ~~~~~~~ for testing
# Call WinBUGS from R (ART 52 h!) and summarize posteriors
# bugs.dir must be set to WinBUGS location, e.g., "c:/WinBUGS14/"
library(R2WinBUGS)
out8 <- bugs(bdata, inits, params, "cjs8.txt", n.chains = nc, n.thin = nt, n.iter = ni,
    n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
print(out8$summary[c(1:7, 221,222),c(1:3,5,7:9)], 3)
# mean sd 2.5% 50% 97.5% Rhat n.eff
# mean.phi 0.287 0.0187 0.251397 0.2874 0.324 1.00 730
# mean.p 0.381 0.0356 0.310497 0.3797 0.450 1.01 380
# mu.lphi -0.910 0.0918 -1.091025 -0.9079 -0.736 1.00 730
# mu.lp -0.490 0.1519 -0.797902 -0.4907 -0.201 1.01 390
# sd.lp.site 1.129 0.1302 0.887167 1.1220 1.407 1.00 3000
# sd.lphi.time 0.220 0.0834 0.101800 0.2063 0.423 1.00 3000
# sd.lp.time 0.148 0.1117 0.006906 0.1266 0.424 1.00 3000
# veta 0.060 0.0765 0.000452 0.0327 0.262 1.08 38
# sdeta 0.204 0.1361 0.021259 0.1809 0.512 1.08 38

# 3.4.5 FITTING A SPATIAL HIERARCHICAL CJS MODEL WITH SPATIAL AUTOCORRELATION AND WITH COVARIATES
# -----------------------------------------------------------------------------------------------
# Bundle and summarize data set for WinBUGS
str(bdata <- list(MARRWB = MARRWB, R = R, n.site = nsite, n.occ = nyear, n.block = nblock,
    BlockID = CES$BlockID, adj = winnb$adj, weights = winnb$weights, num = winnb$num,
    gdd1.site = gdd1.site, gdd2.site = gdd2.site, lat.site = lat.site))
# List of 12
# $ MARRWB : num [1:193, 1:10, 1:11] 1 3 0 2 0 0 0 0 0 0 ...
# $ R : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site : int 193
# $ n.occ : int 11
# $ n.block : int 495
# $ BlockID : int [1:193] 25 77 204 110 222 283 119 234 295 152 ...
# $ adj : int [1:3458] 2 3 4 1 3 5 6 1 2 4 ...
# $ weights : num [1:3458] 1 1 1 1 1 1 1 1 1 1 ...
# $ num : int [1:495] 3 4 6 5 3 6 6 6 5 5 ...
# $ gdd1.site : num [1:193] 1.13 0.726 0.46 0.538 0.686 ...
# $ gdd2.site : num [1:193] 1.228 0.703 0.38 0.473 0.653 ...
# $ lat.site : num [1:193] -1.397 -1.11 -0.403 -0.881 -0.308 ...

# Specify model in BUGS language
cat(file = "cjs9.txt","
model {
  # Priors and linear models
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      phi[t, s] <- 1 / (1 + exp(-lphi[t, s])) # survival
      p[t, s] <- 1 / (1 + exp(-lp[t, s])) # recapture
      lphi[t, s] <- alpha.lphi.site[s] + beta.lphi.time[t]
      lp[t, s] <- alpha.lp.site[s] + beta.lp.time[t]
    }
    # Linear model for site-level effects: add covariates
    alpha.lphi.site[s] <- mu.lphi + beta1 * gdd1.site[s] + beta2 * gdd2.site[s] + beta3 *
    lat.site[s] + eta[BlockID[s]]
    alpha.lp.site[s] ~ dnorm(mu.lp, tau.lp.site) I(-12, 12)
    # backtransform site means
    mean.p.site[s] <- 1 / (1 + exp(-alpha.lp.site[s]))
  }
  for (t in 1:(n.occ-1)){
    beta.lphi.time[t] ~ dnorm(0, tau.lphi.time) I(-12, 12)
    beta.lp.time[t] ~ dnorm(0, tau.lp.time) I(-12, 12)
    # backtransform time means
    mean.phi.time[t] <- 1 / (1 + exp(-mu.lphi + beta.lphi.time[t]))
    mean.p.time[t] <- 1 / (1 + exp(-mu.lphi + beta.lp.time[t]))
  }
  # Hyperpriors for hyperparams
  mu.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0, 1)
  mu.lp <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  tau.lp.site <- pow(sd.lp.site, -2)
  sd.lp.site ~ dunif(0, 2)
  tau.lphi.time <- pow(sd.lphi.time, -2)
  sd.lphi.time ~ dunif(0, 1)
  tau.lp.time <- pow(sd.lp.time, -2)
  sd.lp.time ~ dunif(0, 1)
  # Coefficients for gdd1, gdd2 and lat
  beta1 ~ dunif(-3, 3)
  beta2 ~ dunif(-3, 3)
  beta3 ~ dunif(-3, 3)
  # CAR prior distribution for spatial random effects
  # NOTE: this is defined on the entire grid of 495 blocks
  eta[1:n.block] ~ car.normal(adj[], weights[], num[], tau)
  tau ~ dgamma(0.5, 0.0005) # cf. Lunn et al. (2013)
  veta <- 1/tau
  sdeta <- sqrt(veta)
  # Multinomial likelihood for the m-array data (WinBUGS style)
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      MARRWB[s, t,1:n.occ] ~ dmulti(pr[t,s, ], R[t,s])
    }
  }
  # Define the cell probabilities of the m-array
  # Main diagonal
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      q[t,s] <- 1-p[t,s] # Probability of non-recapture
      pr[t,s,t] <- phi[t,s]*p[t,s]
      # Above main diagonal
      for (j in (t+1):(n.occ-1)){
        pr[t,s,j] <- prod(phi[t:j,s])*prod(q[t:(j-1),s])*p[j,s]
      } #j
      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,s,j] <- 0
      } #j
    } #t
  } #s
  # Last column of m-array: probability of non-recapture
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      pr[t,s,n.occ] <- 1-sum(pr[t,s,1:(n.occ-1)])
    } #t
  }#s
}
")

# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1), eta = rep(0, nblock),
    beta1 = rnorm(1)/3, beta2 = rnorm(1)/3, beta3 = rnorm(1) /3, tau = runif(1))}
# Parameters monitored
params <- c("mean.phi", "mean.p", "mu.lphi", "mu.lp",
    "sd.lp.site", "sd.lphi.time", "sd.lp.time", "mean.p.site", "mean.phi.time", "mean.p.time",
    "veta", "sdeta", "beta1", "beta2", "beta3", "eta")

# MCMC settings
# ni <- 200000 ; nt <- 100 ; nb <- 100000 ; nc <- 3 # c. 2 wk
ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3 # ~~~~~~~~~~~~ for testing
# You may have to launch WinBUGS a couple of times until you don’t get an “undefined real result”
# crash (or alternatively, you could “stabilize” the logit, see Trick 15 in Appendix 1
#  in Kéry and Schaub, 2012).
# Call WinBUGS from R (ART = loooong) and summarize posteriors (for a subset of parameters only)
library(R2WinBUGS)
out9 <- bugs(bdata, inits, params, "cjs9.txt", n.chains = nc, n.thin = nt, n.iter = ni,
    n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
print(out9$summary[c(1:7, 221:225),c(1:3,5,7:9)], 3)
# mean sd 2.5% 50% 97.5% Rhat n.eff
# mean.phi 0.2675 0.0206 0.229485 0.26670 0.310 1.00 3000
# mean.p 0.3824 0.0338 0.316297 0.38200 0.450 1.00 1100
# mu.lphi -1.0099 0.1051 -1.211075 -1.01100 -0.799 1.00 3000
# mu.lp -0.4820 0.1440 -0.770615 -0.48100 -0.202 1.00 1000
# sd.lp.site 1.1052 0.1283 0.864600 1.09900 1.379 1.00 3000
# sd.lphi.time 0.2186 0.0845 0.097788 0.20405 0.417 1.00 3000
# sd.lp.time 0.1554 0.1140 0.009250 0.13480 0.430 1.00 1500
# veta 0.0137 0.0380 0.000201 0.00212 0.118 1.06 48
# sdeta 0.0788 0.0867 0.014189 0.04600 0.344 1.06 48
# beta1 1.0107 0.5355 0.045398 1.01600 2.008 1.00 3000
# beta2 -0.6984 0.4467 -1.527000 -0.70535 0.195 1.00 3000
# beta3 0.2984 0.0976 0.115275 0.29945 0.490 1.01 320

mean(out9$sims.list$beta1>0)
mean(out9$sims.list$beta2<0)
mean(out9$sims.list$beta3>0)
# [1] 0.9696667 # Prob gdd1 has positive effect
# [1] 0.9423333 # Prob gdd2 has negative effect
# [1] 0.9993333 # Prob lat has positive effect






