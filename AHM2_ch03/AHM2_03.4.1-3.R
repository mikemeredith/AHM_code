#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 3 : HIERARCHICAL MODELS OF SURVIVAL
# ===========================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 1 hr
# Run time with the full number of iterations: 10.5 hrs

library(AHMbook)
library(jagsUI)

# 3.4 Spatial hierarchical CJS models
# ===================================

# 3.4.1 British willow warbler data
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
# ..$ cesx   : num [1:193] 307500 432500 457500 312500 477500 ...
# ..$ cesy   : num [1:193] 92500 167500 352500 227500 377500 ...
# ..$ BlockID: num [1:193] 25 77 204 110 222 283 119 234 295 152 ...
# ..$ CellID : num [1:193] 319 1366 4302 2266 4627 ...
# $ blocks:'data.frame': 495 obs. of 2 variables:
# ..$ blockX: num [1:495] 175 150 175 200 150 175 200 225 250 275 ...
# ..$ blockY: num [1:495] 0 25 25 25 50 50 50 50 50 50 ...

attach (willowWarbler)
ch <- as.matrix(birds[, 1:11])
sitevec <- birds$cesID

# Number of captures per bird, year and site
table(apply(ch, 1, sum))    # Table of capture frequency per bird
apply(ch, 2, sum)           # Number of birds per year
plot(table(table(sitevec))) # Frequency distribution of number of
                            # birds captured per site (not shown)
summary(as.numeric(table(sitevec)))
#    1    2   3  4  5 6
# 8997 1214 279 48 12 1

#     1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996
# [1]  692  795 1022 1423 1372 1231 1283 1331 1575 1516  280

#  Min. 1st Qu. Median  Mean 3rd Qu.   Max.
# 2.00    14.00  30.00 54.67   76.00 338.00

# Map of GDD covariate values and 193 CES locations (Fig. 3.11)
library(raster)
mapPalette <- colorRampPalette(c("gray", "yellow", "orange", "red"))
r1 <- with(cells, rasterFromXYZ(data.frame(x = lon, y = lat, z = gdd)))
plot(r1, col = mapPalette(100), axes = FALSE, box = FALSE,
    main ="Map of GDD covariate with 193 CES locations")
with(CES, points(cesx, cesy, pch=16, col='blue', cex = 0.8))

# Get sample sizes
(nyear <- ncol(ch))        # Number years: 11
(nsite <- nrow(CES))       # Number of CE sites: 193
(nblock <- nrow(blocks))   # Number of blocks: 495

(marr <- ch2marray(ch))
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# [1,]    55   13    1    1    0    0    0    0    0     0   622
# [2,]     0  123   17    2    1    0    0    0    0     0   652
# [3,]     0    0  143   20    4    0    0    0    0     0   855
# [4,]     0    0    0  179   16    6    0    0    0     0  1222
# [5,]     0    0    0    0  183   14    1    0    0     0  1174
# [6,]     0    0    0    0    0  161   16    5    0     0  1049
# [7,]     0    0    0    0    0    0  172   20    4     0  1087
# [8,]     0    0    0    0    0    0    0  239   41     9  1042
# [9,]     0    0    0    0    0    0    0    0  252    24  1299
# [10,]    0    0    0    0    0    0    0    0    0   247  1269

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

# 3.4.2 A hierarchical CJS model with random site and random year effects
# -----------------------------------------------------------------------

# Bundle and summarize data set
str(bdata <- list(MARR = MARR, R = R, n.site = nsite, n.occ = nyear))
# List of 4
# $ MARR   : num [1:10, 1:11, 1:193] 1 0 0 0 0 0 0 0 0 0 ...
# $ R      : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site : int 193
# $ n.occ  : int 11

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
      }
    }
  }
  # Last column: probability of non-recapture
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      pr[t,n.occ,s] <- 1-sum(pr[t,1:(n.occ-1),s])
    }
  }
}
")

# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1))}

# Parameters monitored
params <- c("mean.phi", "mean.p", "sd.lphi.site", "sd.lp.site",
    "sd.lphi.time", "sd.lp.time", "mean.phi.site", "mean.p.site",
    "mean.phi.time", "mean.p.time")

# MCMC settings
# na <- 1000 ; ni <- 30000 ; nt <- 10 ; nb <- 20000 ; nc <- 3
na <- 1000 ; ni <- 3000 ; nt <- 1 ; nb <- 2000 ; nc <- 3  # ~~~ for testing

# Call JAGS (ART 192 min), check convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "cjs6.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,3))  #  ~~~ replace with 'layout' argument
traceplot(out6, layout=c(2,3))
print(out6, 3)
#                   mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# mean.phi         0.281 0.020 0.243 0.281 0.321    FALSE 1 1.031   104
# mean.p           0.380 0.038 0.310 0.379 0.458    FALSE 1 1.002   779
# sd.lphi.site     0.158 0.078 0.031 0.153 0.329    FALSE 1 1.007   279
# sd.lp.site       1.138 0.135 0.888 1.135 1.413    FALSE 1 1.001  3000
# sd.lphi.time     0.223 0.086 0.100 0.208 0.431    FALSE 1 1.006  1007
# sd.lp.time       0.149 0.108 0.005 0.130 0.389    FALSE 1 1.014   226
# mean.phi.site[1] 0.277 0.040 0.187 0.280 0.354    FALSE 1 1.016   175
# mean.phi.site[2] 0.274 0.034 0.205 0.275 0.335    FALSE 1 1.006   356
# [ ... output truncated ... ]

# Visualizations (Fig. 3.12)
op <- par(mfrow = c(2, 1), mar = c(3,5,2,1))
plot(1:193, out6$mean$mean.phi.site, xlab = 'Site', ylab = 'Apparent survival',
    frame = FALSE, pch = 16, ylim = c(0.14, 0.45))
segments(1:193, out6$q2.5$mean.phi.site, 1:193, out6$q97.5$mean.phi.site)
abline(h = out6$mean$mean.phi, lty = 1, lwd = 2)
plot(1:193, out6$mean$mean.p.site, xlab = 'Site', ylab = 'Recapture',
    frame = FALSE, pch = 16, ylim = c(0, 1))
segments(1:193, out6$q2.5$mean.p.site, 1:193, out6$q97.5$mean.p.site)
abline(h = out6$mean$mean.p, lty = 1, lwd = 2)
par(op)

# 3.4.3 Adding spatial covariates into the hierarchical CJS model
# ---------------------------------------------------------------

# Scale linear and squared GDD separately and latitude (for entire grid)
scaled.gdd1 <- standardize(cells$gdd)
scaled.gdd2 <- standardize(cells$gdd^2)
scaled.lat <- standardize(cells$lat)

# Pull out GDD and LATITUDE values for 193 sites (using CellID)
gdd1.site <- scaled.gdd1[CES$CellID]
gdd2.site <- scaled.gdd2[CES$CellID]
lat.site <- scaled.lat[CES$CellID]

# Bundle and summarize data set
str(bdata <- list(MARR = MARR, R = R, n.site = nsite, n.occ = nyear,
    gdd1.site = gdd1.site, gdd2.site = gdd2.site, lat.site = lat.site))
    # , n.cell = n.cell, gdd1.grid = scaled.gdd, gdd2.grid = scaled.gdd2, lat.grid = lat # if want to make predictions inside BUGS
# List of 7
# $ MARR     : num [1:10, 1:11, 1:193] 1 0 0 0 0 0 0 0 0 0 ...
# $ R        : num [1:10, 1:193] 13 5 0 0 0 0 0 0 0 0 ...
# $ n.site   : int 193
# $ n.occ    : int 11
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
      }
    }
  }

  # Last column of m-array: probability of non-recapture
  for (s in 1:n.site){
    for (t in 1:(n.occ-1)){
      pr[t,n.occ,s] <- 1-sum(pr[t,1:(n.occ-1),s])
    }
  }

  # Derived quantities: predictions for entire grid
  #for (s in 1:n.cell){
  # phi.grid[s] <- ilogit(alpha.mu.lphi + beta1 * gdd1.grid[s] + beta2 * #gdd2.grid[s] + beta3 * lat.grid[s])
  #}
}
")

# Initial values
inits <- function(){list(mean.phi = runif(1), mean.p = runif(1),
    beta1 = rnorm(1), beta2 = rnorm(1), beta3 = rnorm(1))}

# Parameters monitored
params <- c("mean.phi", "mean.p", "alpha.mu.lphi", "mu.lp", "sd.lphi.site",
    "sd.lp.site", "sd.lphi.time", "sd.lp.time", "mean.phi.site", "mean.p.site",
    "mean.phi.time", "mean.p.time", "beta1", "beta2", "beta3")
    # , "phi.grid" # for predictions within JAGS

# MCMC settings
# na <- 5000 ; ni <- 60000 ; nt <- 30 ; nb <- 30000 ; nc <- 3
na <- 5000 ; ni <- 6000 ; nt <- 3 ; nb <- 3000 ; nc <- 3 # ~~~~~~~~~ for testing

# Call JAGS (ART 330 min), check convergence and summarize posteriors
out7 <- jags(bdata, inits, params, "cjs7.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,3))  #  ~~~ replace with 'layout' argument
traceplot(out7, layout=c(2,3))
print(out7, 3)
#                    mean    sd   2.5%    50%  97.5% overlap0     f  Rhat n.eff
# mean.phi          0.260 0.023  0.216  0.259  0.309    FALSE 1.000 1.029    94
# mean.p            0.401 0.040  0.327  0.399  0.484    FALSE 1.000 1.007   401
# alpha.mu.lphi    -1.048 0.120 -1.291 -1.052 -0.806    FALSE 1.000 1.028    97
# mu.lp            -0.404 0.168 -0.721 -0.408 -0.065    FALSE 0.989 1.006   412
# sd.lphi.site      0.171 0.094  0.005  0.169  0.369    FALSE 1.000 1.029   105
# sd.lp.site        1.065 0.141  0.793  1.064  1.350    FALSE 1.000 1.002   695
# sd.lphi.time      0.223 0.087  0.100  0.210  0.436    FALSE 1.000 1.003   754
# sd.lp.time        0.158 0.117  0.007  0.135  0.439    FALSE 1.000 1.002  1794
# mean.phi.site[1]  0.227 0.039  0.142  0.229  0.303    FALSE 1.000 1.007   979
# mean.phi.site[2]  0.243 0.032  0.178  0.243  0.310    FALSE 1.000 1.008   359
# [ ... output truncated ... ]
# mean.p.time[9]    0.409 0.046  0.326  0.408  0.503    FALSE 1.000 1.003   788
# mean.p.time[10]   0.403 0.050  0.310  0.400  0.508    FALSE 1.000 1.005   448
# beta1             1.043 0.600 -0.113  1.030  2.188     TRUE 0.960 1.021   125
# beta2            -0.735 0.506 -1.683 -0.730  0.248     TRUE 0.931 1.016   153
# beta3             0.320 0.102  0.123  0.323  0.519    FALSE 0.999 1.003   664

ncell <- nrow(cells)
phi.grid <- array(NA, dim = c(ncell, out7$mcmc.info$n.samples))

# Derived quantities: predictions for entire grid
sims <- out7$sims.list                # Grab the posterior simulations
for (s in 1:ncell){
  phi.grid[s,] <- plogis(sims$alpha.mu.lphi + sims$beta1 * scaled.gdd1[s] + sims$beta2 *
      scaled.gdd2[s] + sims$beta3 * scaled.lat[s])
}
post.mean <- apply(phi.grid, 1, mean) # Posterior mean
post.sd <- apply(phi.grid, 1, sd)     # Posterior standard deviation


# ~~~~~~~~~ code for figure 3.13 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(raster)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

op <- par(mfrow = c(1, 2))#, mar = c(2, 6, 2, 8), cex.lab = 2, cex.axis = 2)
# Plot posterior mean of predicted apparent survival
r1 <- rasterFromXYZ(data.frame(x = cells$lon, y = cells$lat, z = post.mean))
plot(r1, col = mapPalette(100), axes = FALSE, box = FALSE)
# Plot uncertainty in this estimate of predicted apparent survival
r2 <- rasterFromXYZ(data.frame(x = cells$lon, y = cells$lat, z = post.sd))
plot(r2, col = mapPalette(100), axes = FALSE, box = FALSE, zlim = c(0, 0.1))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
