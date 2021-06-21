#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10 : INTEGRATED MODELS FOR MULTIPLE TYPES OF DATA
# =========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 68 mins
# Run time with the full number of iterations: 19 hrs

library(AHMbook)
library(jagsUI)

# 10.6 Example 4: Combination of counts plus opportunistic presence-only data
# ===========================================================================

# 10.6.1 Data simulation under the model of Dorazio (2014)
# --------------------------------------------------------

# Call function with default values for arguments
str(dat <- simDataDK(sqrt.npix = 100, alpha = c(-1,-1), beta = c(6,0.5),
    drop.out.prop.pb = 0.7, quadrat.size = 4, gamma = c(0,-1.5),
    nquadrats = 250, nsurveys = 3, show.plot = TRUE), 1)

# Finer pixel approximation of the continuous landscape
# yields 1 Million pixels (as in Dorazio (2014)) ... too big for BUGS !
str(dat <- simDataDK(sqrt.npix = 1000), 1) # Yields 1 Million pixels

# No spatial bias in point pattern data set 1: no effect of covariate W
str(dat <- simDataDK(alpha = c(-1, 0)), 1)

# Homogeneous point pattern: no effect of covariate X
str(dat <- simDataDK(beta = c(6, 0)), 1) # don't make intercept too big

# No effect of covariate W on detection in point count sampling
str(dat <- simDataDK(gamma = c(0, 0)), 1)

# No thinning of the original point pattern
str(dat <- simDataDK(alpha = c(20, 0), drop.out.prop.pb = 0), 1)

# Much larger quadrats for point count sampling
str(dat <- simDataDK(quadrat.size = 20, nquadrats = 12), 1)

RNGkind(sample.kind = "Rounding")
# Call function with default values for arguments
set.seed(111)
str(dat <- simDataDK(sqrt.npix = 100, alpha = c(-2,-1), beta = c(6,0.5),
    drop.out.prop.pb = 0, quadrat.size = 4, gamma = c(0,-1.5),
    nquadrats = 250, nsurveys = 3, show.plot = TRUE), 1)

# Inspect the presence-only data generated
head(dat$loc.det) # Actual coordinates of 294 points detected
#               x          y
# [1,]  0.5388595  0.2530584
# [2,]  0.2553583  0.6808913
# [3,] -0.8700825  0.3739989
# [4,]  0.9571221  0.2905681
# [5,]  0.1677716 -0.1535967
# [6,]  0.8483838  0.4965985

head(dat$pixel.id.det) # Pixel ID of 294 points detected
# [1] 3777 1563 3107 3598 5759 2593

# Are there any pixels with more than one point ?
table(table(dat$pixel.id.det))
#   1 2
# 286 4

# Overview of the point count data set (with 250 quads)
head(dat$countData)
#      quadID          x           w N
# [1,]      2  0.4697062  0.42953716 6 3 2 4
# [2,]     10 -0.1762401  0.91133985 1 0 0 0
# [3,]     11 -0.2730764  0.72495689 4 2 1 1
# [4,]     14 -0.4757779  0.02358506 4 2 3 1
# [5,]     15 -0.5031330 -0.21286719 1 0 0 0
# [6,]     17 -0.4987504 -0.62829456 2 1 1 1


# 10.6.2 Fitting a poisson point pattern (PPP) model to the opportunistic
#        data alone: the point process variant of a relative-abundance model
# --------------------------------------------------------------------------

# Get point indicators: turn presence-only data into 0/1 data
y <- numeric(dat$npix)
y[dat$pixel.id.det] <- 1

# Get pixel area (for use as an offset)
logarea <- log(dat$s.area / dat$npix)

# Bundle data
str(bdata <- list(y = y, logarea = logarea, npix = dat$npix, xcov = dat$xcov))
# List of 4
# $ y      : num [1:10000] 0 0 0 0 0 0 0 0 0 0 ...
# $ logarea: num -7.82
# $ npix   : num 10000
# $ xcov   : num [1:10000] 0.374 0.372 0.369 0.365 0.36 ...

# Specify model in BUGS language
cat(file = "ipp.txt", "
model {

  # Priors
  beta0 ~ dnorm(0, 0.01)              # Regression params for intensity
  beta1 ~ dnorm(0, 0.1)
  # Howsit look like (execute in R)?
  # curve(dnorm(x, mean=0, sd=sqrt(1/0.01)), -20, 20)

  # Likelihood: binary regression with cloglog link and area offset
  # as an approximation to the Poisson point process (PPP) model
  for(i in 1:npix){
    y[i] ~ dbern(theta[i])
    cloglog(theta[i]) <- logarea + log(lambda[i]) # cloglog with offset logarea
    log(lambda[i]) <- beta0 + beta1 * xcov[i]
  }
  # Derived quantity
  N <- sum(lambda[]* exp(logarea))     # Population size
}
")

# Initial values
inits <- function() {list(beta0 = runif(1), beta1 = runif(1)) }

# Parameters monitored
params <- c("beta0", "beta1", "N")

# MCMC settings
# na <- 1000 ; ni <- 5000 ; nt <- 2 ; nb <- 3000 ; nc <- 3
na <- 1000 ; ni <- 500 ; nt <- 1 ; nb <- 300 ; nc <- 3  # ~~~~ for testing, 3 mins

# Call JAGS (ART 8 min), assess convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "ipp.txt", n.adapt = na, n.thin = nt,
    n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
# par(mfrow = c(2,2))  # ~~~ replaced with 'layout' argument
traceplot(out1, layout=c(2,2))
print(out1, 3)
#          mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# beta0   4.144  0.069   4.006   4.146   4.276    FALSE 1 1.009   256
# beta1   0.587  0.074   0.439   0.585   0.731    FALSE 1 1.014   152
# N     295.373 17.408 261.973 295.407 327.695    FALSE 1 1.000  1500


# 10.6.3 Accounting for spatial sampling and detection bias in the PPP model
#        by accommodation of thinning
# --------------------------------------------------------------------------

# Bundle data
str(bdata <- list(y = y, logarea = logarea, npix = dat$npix, xcov = dat$xcov,
    wcov = dat$wcov))
# List of 5
# $ y      : num [1:10000] 0 0 0 0 0 0 0 0 0 0 ...
# $ logarea: num -7.82
# $ npix   : num 10000
# $ xcov   : num [1:10000] 0.374 0.372 0.369 0.365 0.36 ...
# $ wcov   : num [1:10000] 0.0631 0.123 0.1829 0.2427 0.3022 ...

# Specify model in BUGS language
cat(file = "thinned.ipp.txt", "
model {

  # Priors
  beta0 ~ dnorm(0, 0.01)      # Regression params for intensity
  beta1 ~ dnorm(0, 0.1)
  alpha0 ~ dnorm(0, 0.01)     # Regression parameters for thinning
  alpha1 ~ dnorm(0, 0.1)

  # Likelihood: binary regression with cloglog link, area offset
  # and thinning parameter as an approximation to a thinned PPP model
  for(i in 1:npix){
    # # Variant 1: # Model occupancy through underlying abundance process
    # y[i] ~ dbern(1-exp(-exp(theta[i]) * p[i]))
    # theta[i] <- logarea + log(lambda[i])
    # log(lambda[i]) <- beta0 + beta1 * xcov[i] # model for intensity
    # logit(p[i]) <- alpha0 + alpha1 * wcov[i] # model for thinning
    # Variant 2: # Same, simply via the cloglog link !
    y[i] ~ dbern(psi[i])
    cloglog(psi[i]) <- logarea + log(lambda[i]) + log(p[i])
    log(lambda[i]) <- beta0 + beta1 * xcov[i] # model for intensity
    logit(p[i]) <- alpha0 + alpha1 * wcov[i] # model for thinning
  }
  # Derived quantity
  mean.p <- ilogit(alpha0)          # Intercept thinning probability
  N <- sum(lambda[]* exp(logarea))  # Population size
}
")

# Initial values
inits <- function() {
  list(alpha0 = runif(1), alpha1 = runif(1), beta0 = runif(1),
        beta1 = runif(1))
}

# Parameters monitored
params <- c("mean.p", "alpha0", "alpha1", "beta0", "beta1", "N")

# MCMC settings
# na <- 5000 ; ni <- 200000 ; nt <- 100 ; nb <- 100000 ; nc <- 3
na <- 5000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ for testing, 30 mins

# Call JAGS (ART 892 min), assess convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "thinned.ipp.txt", n.adapt = na, n.thin = nt,
    n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
# par(mfrow = c(2,2))  # ~~~ replaced with 'layout' argument
traceplot(out2, layout=c(2,2))
print(out2, 3)
#              mean         sd    2.5%      50%      97.5% overlap0    f  Rhat n.eff
# mean.p      0.154      0.147   0.000    0.110      0.482    FALSE 1.00 1.023    93
# alpha0     -2.640      2.101  -8.233   -2.092     -0.071    FALSE 0.98 1.104    38
# alpha1     -0.973      0.247  -1.588   -0.903     -0.660    FALSE 1.00 1.018   129
# beta0       6.756      1.907   4.828    6.124     12.101    FALSE 1.00 1.119    36
# beta1       0.580      0.076   0.437    0.577      0.732    FALSE 1.00 1.001  2424
# N      120506.472 927113.374 581.071 2133.640 911248.463    FALSE 1.00 1.310    81


# Compare estimates of thinned PPM with truth in data simulation
# ~~~ code to produce the table ~~~~~~~~~~
truth <- c(dat$alpha, dat$beta, dat$N.ipp)
esti <- out2$summary[2:6, c(1,3,7)]
print(cbind(truth, esti), 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         truth       mean    2.5%       97.5%
# alpha0   -2.0     -2.640  -8.233     -0.0714
# alpha1   -1.0     -0.973  -1.588     -0.6596
# beta0     6.0      6.756   4.828     12.1007
# beta1     0.5      0.580   0.437      0.7324
# N      1831.0 120506.472 581.071 911248.4629


# 10.6.4 Fitting the integrated model: a bugs implementation of the
#        model of Dorazio (2014)
# ------------------------------------------------------------------

# Bundle and summarize data set
str(bdata <- list(
    # Data for thinned point process submodel:
    y = y, logarea = logarea, npix = dat$npix, xcov = dat$xcov, wcov = dat$wcov,
    # Data for replicated counts (Nmix) submodel:
    C = dat$countData[,5:7], covarX = dat$countData[,'x'],
        covarW = dat$countData[,'w'],
    nsite = dat$nquadrats, nsurveys = dat$nsurveys,
    area = dat$s.area * (dat$quadrat.size^2 / dat$npix) ) )
# List of 11
# $ y       : num [1:10000] 0 0 0 0 0 0 0 0 0 0 ...
# $ logarea : num -7.82
# $ npix    : num 10000
# $ xcov    : num [1:10000] 0.374 0.372 0.369 0.365 0.36 ...
# $ wcov    : num [1:10000] 0.0631 0.123 0.1829 0.2427 0.3022 ...
# $ C       : num [1:250, 1:3] 3 0 2 2 0 1 1 0 2 1 ...
# $ covarX  : num [1:250] 0.47 -0.176 -0.273 -0.476 -0.503 ...
# $ covarW  : num [1:250] 0.4295 0.9113 0.725 0.0236 -0.2129 ...
# $ nsite   : num 250
# $ nsurveys: num 3
# $ area    : num 0.0064

# Specify model in BUGS language
cat(file = "JointPPPNmix.txt", "
model {

  # Priors for shared intensity parameters: this is the SDM
  beta0 ~ dnorm(0, 0.01) # intensity intercept
  beta1 ~ dnorm(0, 0.1) # slope of intensity on X

  # Submodel 1 for presence-only data (thinned PPP)
  # -----------------------------------------------
  # Priors
  alpha0 ~ dnorm(0,0.01)                # thinning intercept
  alpha1 ~ dnorm(0, 0.1)                # slope of thinning on W

  # Likelihood
  # logistic regression with cloglog link models intensity of underlying Poisson PP
  for(i in 1:npix){
    y[i] ~ dbern(psi[i])
    # cloglog link with offset logarea
    cloglog(psi[i]) <- logarea + log(lambda[i]) + log(p1[i])
    # linear model for intensity
    log(lambda[i]) <- beta0 + beta1 * xcov[i] # Note same beta as below
    # linear model for thinning probability: we call detection differently,
    # b/c it is different from the p in the Nmix part of the model below
    logit(p1[i]) <- alpha0 + alpha1 * wcov[i]
  }

  # Submodel 2 for replicated counts
  # --------------------------------
  # Priors
  gamma0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)           # Detection intercepts
  gamma1 ~ dnorm(0, 0.1)         # Slope of detection on quadrat average of W
  # Likelihood
  # Ecological model for true abundance
  for (i in 1:nsite){
    N[i] ~ dpois(area * lam[i])   # Here's the area scaling
    log(lam[i]) <- beta0 + beta1 * covarX[i] # note same beta's as above
    # Observation model for replicated counts
    for (j in 1:nsurveys){
    C[i,j] ~ dbin(p2[i], N[i])
  }
  logit(p2[i]) <- gamma0 + gamma1 * covarW[i]
  }

  # Derived quantity: can estimate abundance in two ways
  Nppm <- sum(lambda[]* exp(logarea))     # from the IPP
  Nnmix <- (sum(N) / 250) * 625           # from the Nmix
}
")

# Initial values
Nst <- apply(bdata$C, 1, max)+1 # usual Nmix inits
inits <- function() list(N = Nst, beta1 = rnorm(1), gamma1 = rnorm(1))

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "gamma0", "gamma1",
    "mean.p", "Nppm", "Nnmix")

# MCMC settings
# na <- 1000 ; ni <- 50000 ; nt <- 40 ; nb <- 10000 ; nc <- 3
na <- 1000 ; ni <- 5000 ; nt <- 4 ; nb <- 1000 ; nc <- 3  # ~~~~ for testing, 27 mins

# Call JAGS (ART 282 min), assess convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "JointPPPNmix.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  # ~~~ replaced with 'layout' argument
traceplot(out3, layout=c(2,2))
print(out3, 2)
#           mean    sd    2.5%     50%   97.5% overlap0    f Rhat n.eff
# alpha0   -1.86  0.11   -2.07   -1.86   -1.64    FALSE 1.00    1  1113
# alpha1   -0.92  0.09   -1.10   -0.92   -0.75    FALSE 1.00    1  3000
# beta0     5.93  0.06    5.82    5.93    6.04    FALSE 1.00    1  3000
# beta1     0.55  0.04    0.46    0.55    0.63    FALSE 1.00    1  3000
# gamma0    0.14  0.09   -0.03    0.14    0.33     TRUE 0.94    1  1965
# gamma1   -1.55  0.09   -1.73   -1.55   -1.37    FALSE 1.00    1  3000
# mean.p    0.54  0.02    0.49    0.54    0.58    FALSE 1.00    1  1969
# Nppm   1729.81 93.60 1549.71 1726.98 1918.92    FALSE 1.00    1  3000
# Nnmix  1664.60 63.50 1545.00 1662.50 1795.00    FALSE 1.00    1  1711

# Compare estimates of joint model with truth in data simulation
# ~~~~ extra code for the table ~~~~~~~~~~~
truth <- c(dat$alpha, dat$beta, dat$gamma, dat$N.ipp)
esti <- out3$summary[c(1:6,8), c(1,3,7)]
print(cbind(truth, esti), 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         truth     mean      2.5%    97.5%
# alpha0   -2.0   -1.859   -2.0743   -1.636
# alpha1   -1.0   -0.924   -1.1045   -0.745
# beta0     6.0    5.933    5.8199    6.045
# beta1     0.5    0.546    0.4630    0.630
# gamma0    0.0    0.143   -0.0349    0.330
# gamma1   -1.5   -1.548   -1.7293   -1.373
# Nppm   1831.0 1729.814 1549.7149 1918.916


# ~~~ extra code for figure 10.6 ~~~~~~~~~~~~~~
# Form predictions for intensity lambda: this is the SDM !
library(raster)
nsims <- length(out3$sims.list$beta0)
sims <- out3$sims.list
pred.lam <- array(NA, dim = c(nsims, dat$npix))
for(i in 1:nsims){
  pred.lam[i, ] <- exp(sims$beta0[i] + sims$beta1[i] * values(dat$s)[,'x'])
}
pred.lam.pm <- apply(pred.lam, 2, mean)
pred.lam.psd <- apply(pred.lam, 2, sd)

# Define a raster object to be filled with data
pred <- raster(ncol=dat$sqrt.npix, nrow=dat$sqrt.npix, xmn=-1, xmx=1,
    ymn=-1, ymx=1)
pred.loc <- xyFromCell(pred, 1:ncell(pred)) # Coordinates of every cell

# Fill raster with values of posterior mean and sd
pred.pm <- raster(pred)
values(pred.pm) <- pred.lam.pm
names(pred.pm) <- 'posterior_mean'
pred <- addLayer(pred, pred.pm)
pred.sd <- raster(pred)
values(pred.sd) <- pred.lam.psd
names(pred.sd) <- 'posterior_standardDeviation'
pred <- addLayer(pred, pred.sd)
plot(pred, asp = 1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 10.6.5 Assumptions and implications of the thinned PPP-based SDM model
# (no code)

# 10.6.6 Relationship between the thinned PPP-SDM and the single-visit
#   occupancy model of Lele et al. (2012) (no code)

# 10.6.7 Concluding remarks on the analysis of citizen science data using
#   integrated models (no code)
