#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 10 mins

library(AHMbook)
library(jagsUI)
library(unmarked)

# 2.7 Modeling dynamics with multinomial N-mixture models
# =======================================================

# 2.7.1 Modeling TE with a three-level multinomial-mixture model in unmarked using function gmultmix
# ------------------------------------------------------------------------------------------

# Define function for data simulation
simMultMix <- function(nsites = 100, nyears = 4, nsurveys = 3, lambda = 3,
    theta = 0.5, p = 0.3){
  # Simulate data using the multinomial-Poisson model with a
  # repeated constant-interval removal design (written by R.B. Chandler)
  #
  # lambda, theta and p: expected abundance, availability and detection prob.
  y <- array(NA, c(nsites, nyears, nsurveys))
  M <- rpois(nsites, lambda)      # Local population size
  N <- matrix(NA, nsites, nyears) # Individuals available for detection

  for(i in 1:nsites) {
    N[i,] <- rbinom(nyears, M[i], theta)
    y[i,,1] <- rbinom(nyears, N[i,], p) # Observe some
    Nleft1 <- N[i,] - y[i,,1] # Remove them
    y[i,,2] <- rbinom(nyears, Nleft1, p) # ...
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(nyears, Nleft2, p)
  }
  y2d <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
  return(list(nsites = nsites, nyears = nyears, nsurveys = nsurveys,
      lambda = lambda, theta = theta, p = p, M = M, N = N, y = y, y2d = y2d))
}

# Execute function
set.seed(24)
str(data <- simMultMix(nsites = 100, nyears = 4, nsurveys = 3, lambda = 3,
    theta = 0.5, p = 0.3) )

# Package data in unmarked format GMM, fit model and inspect estimates
library(unmarked)
umf <- unmarkedFrameGMM(y = data$y2d, numPrimary = data$nyears, type = "removal")
(fm <- gmultmix(~1, ~1, ~1, data = umf, K = 30))
backTransform(fm, type = "lambda")     # Individuals per plot
backTransform(fm, type = "phi")        # Prob(available) (phi)
(p <- backTransform(fm, type = "det")) # Prob. of detection
p <- coef(p)

# Multinomial cell probabilities under removal design
c(p, (1-p) * p, (1-p)^2 * p)

# Or more generally:
head(getP(fm))

# Empirical Bayes estimates of super-population size (see Fig. 2.11)
re <- ranef(fm)
plot(re, layout=c(5, 1), xlim = c(-1, 20), subset = site%in%1:5, lwd = 5)


# 2.7.2 Chandler’s alder flycatcher data
# --------------------------------------

alfl <- read.csv(system.file("csv", "alfl.csv", package = "unmarked"))
alfl.covs <- read.csv(system.file("csv", "alflCovs.csv",
    package = "unmarked"), row.names = 1)
head(alfl.covs)

alfl$captureHistory <- paste(alfl$interval1, alfl$interval2, alfl$interval3, sep = "")
alfl$captureHistory <- factor(alfl$captureHistory,
    levels = c("001", "010", "011", "100", "101", "110", "111"))
alfl$id <- factor(alfl$id, levels = rownames(alfl.covs))

head(alfl, 5)
#          id survey interval1 interval2 interval3 captureHistory
# 1 crick1_05      1         1         1         1            111
# 2 crick1_05      3         1         0         1            101
# 3   his1_05      1         0         1         1             11
# 4   his1_05      1         1         1         1            111
# 5   his1_05      2         0         1         1             11

# Note: we analyzed survey == 1 in AHM1 ch. 7, here use all surveys
alfl.v1 <- alfl[alfl$survey == 1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
alfl.v2 <- alfl[alfl$survey == 2,]
alfl.H2 <- table(alfl.v2$id, alfl.v2$captureHistory)
alfl.v3 <- alfl[alfl$survey == 3,]
alfl.H3 <- table(alfl.v3$id, alfl.v3$captureHistory)

# Arrange the data in 3-d array format and also the wide format
Y <- array(NA, c(50, 3, 7))
Y[1:50,1,1:7] <- alfl.H1
Y[1:50,2,1:7] <- alfl.H2
Y[1:50,3,1:7] <- alfl.H3
Ywide <- cbind(alfl.H1, alfl.H2, alfl.H3)

# unmarkedFrame for a static model just using 1 primary sample
intervalMat <- matrix(c('1', '2', '3'), 50, 3, byrow = TRUE)
class(alfl.H1) <- "matrix"
o2y <- matrix(1, 3, 7)

summary(umf.cr1 <- unmarkedFrameMPois(y = alfl.H1,
    siteCovs = alfl.covs[,c("woody", "struct", "time.1")],
    obsCovs = list(interval = intervalMat), obsToY = o2y,
    piFun = "crPiFun") )

# Standardize some obsCovs
time <- as.matrix(alfl.covs[,c("time.1", "time.2", "time.3")] )
time <- time - median(time)
date <- as.matrix(alfl.covs[,c("date.1", "date.2", "date.3")] )
date <- date - median(date)

occ <- matrix(NA, nrow = 50, ncol = 9)
occ <- col(occ)

summary(alfl.umf <- unmarkedFrameGMM(y = Ywide,
    siteCovs = alfl.covs[,c("woody", "struct")], obsCovs = list(occ = occ),
    numPrimary = 3, yearlySiteCovs = list(time = time, date = date),
    obsToY = o2y, piFun = "crPiFun") )

m0 <- gmultmix( ~1, ~1, ~1, data = alfl.umf, mixture = "P")
m1 <- gmultmix( ~1, ~time, ~1, data = alfl.umf, mixture = "P")
m1b <- gmultmix( ~1, ~time, ~time, data = alfl.umf, mixture = "P")
m1c <- gmultmix( ~1, ~1, ~time, data = alfl.umf, mixture = "P")
m2 <- gmultmix( ~1, ~date, ~1, data = alfl.umf, mixture = "P")
m2b <- gmultmix( ~1, ~date, ~date, data = alfl.umf, mixture = "P")
m2c <- gmultmix( ~1, ~1, ~date, data = alfl.umf, mixture = "P")
m3 <- gmultmix( ~1, ~time + date, ~1, data = alfl.umf, mixture = "P")
m3b <- gmultmix( ~1, ~time + date, ~time, data = alfl.umf, mixture = "P")
m3c <- gmultmix( ~1, ~date, ~time, data = alfl.umf, mixture = "P")
m3d <- gmultmix( ~1, ~time, ~date, data = alfl.umf, mixture = "P")
m3e <- gmultmix( ~1, ~time, ~time + date, data = alfl.umf, mixture = "P")
m3f <- gmultmix( ~1, ~1, ~time + date, data = alfl.umf, mixture = "P")
m3g <- gmultmix( ~1, ~time + date, ~date, data = alfl.umf, mixture = "P")
m3h <- gmultmix( ~1, ~time + date, ~time + date, data = alfl.umf, mixture = "P")
m3i <- gmultmix( ~1, ~date, ~time + date, data = alfl.umf, mixture = "P")

(fl <- modSel(fitList(m0, m1, m1b, m1c, m2, m2b, m2c, m3, m3b, m3c, m3d,
    m3e, m3f, m3g, m3h, m3i)) )
#     nPars    AIC delta  AICwt cumltvWt
# m3g     6 565.79  0.00 0.6600     0.66
# m3h     7 567.46  1.66 0.2900     0.95
# m2b     5 572.27  6.47 0.0260     0.97
# m3      5 573.93  8.14 0.0110     0.98
# m3i     6 574.21  8.41 0.0098     0.99
# m3b     6 575.28  9.48 0.0058     1.00
# [... output truncated ... ]

m3g

# Abundance:
# Estimate SE z P(>|z|)
# 0.406 0.154 2.64 0.00818

# Availability:
# Estimate SE z P(>|z|)
# (Intercept) -0.2607 0.2635 -0.989 3.22e-01
# time -0.3985 0.1423 -2.800 5.11e-03
# date -0.0615 0.0144 -4.277 1.89e-05

# Detection:
# Estimate SE z P(>|z|)
# (Intercept) 0.7173 0.1618 4.43 9.33e-06
# date -0.0369 0.0117 -3.14 1.69e-03

m3g1 <- gmultmix(~woody, ~time + date, ~date, data = alfl.umf, mixture = "P")
m3g2 <- gmultmix(~struct, ~time + date, ~date, data = alfl.umf, mixture = "P")
m3g3 <- gmultmix(~woody + struct, ~time + date, ~date, data = alfl.umf, mixture = "P")
(fl <- modSel(fitList(m0, m1, m1b, m1c, m2, m2b, m2c, m3, m3b, m3c, m3d, m3e,
    m3f, m3g, m3h, m3i, m3g1, m3g2, m3g3)) )
#      nPars    AIC delta   AICwt cumltvWt
# m3g3     8 554.80  0.00 0.51000     0.51
# m3g1     7 554.97  0.17 0.47000     0.99
# m3g2     7 562.60  7.81 0.01000     1.00
# m3g      6 565.79 11.00 0.00210     1.00
# m3h      7 567.46 12.66 0.00092     1.00
# . . . [output truncated] . . .

m3g3

# Abundance:
# Estimate SE z P(>|z|)
# (Intercept) -0.7886 0.4054 -1.95 0.05177
# woody 1.9647 0.6275 3.13 0.00174
# struct 0.0554 0.0364 1.52 0.12849

# Availability:
# Estimate SE z P(>|z|)
# (Intercept) -0.3644 0.2886 -1.26 2.07e-01
# time -0.4020 0.1442 -2.79 5.30e-03
# date -0.0575 0.0144 -4.00 6.46e-05

# Detection:
# Estimate SE z P(>|z|)
# (Intercept) 0.7177 0.1617 4.44 9.13e-06
# date -0.0369 0.0117 -3.14 1.69e-03

# Carry out a parametric bootstrap goodness-of-fit test (ART 7 min)
(pb <- parboot(m3g3, fitstats, nsim = 1000, report = 1, ncores=3))

# Parametric Bootstrap Statistics:

# t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
# SSE 87.1 -14.16 15.4 0.823
# Chisq 1046.5 6.79 92.5 0.434
# freemanTukey 115.8 -0.73 11.9 0.513

# t_B quantiles:
# 0% 2.5% 25% 50% 75% 97.5% 100%
# SSE 58 73 91 100 111 135 158
# Chisq 819 880 977 1030 1095 1244 1477
# freemanTukey 83 93 109 116 125 139 155

# t0 = Original statistic compuated from data
# t_B = Vector of bootstrap samples

# For comparison, we will now stack the data and fit a similar set of models
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Ystacked <- rbind(Ywide[,1:7], Ywide[,8:14], Ywide[,15:21])
sc <- alfl.covs[,c("woody", "struct")]
sc <- rbind(sc, sc, sc)                         # Repeat the site covs matrix
sc <- cbind(sc, date = c(date), time = c(time)) # YearlySiteCovs become siteCovs
summary(umf2 <- unmarkedFrameGMM(y = Ystacked, siteCovs = sc, numPrimary = 1,
    obsToY = o2y, piFun = "crPiFun") )

# Models with time on lambda or p
m0 <- gmultmix( ~1, ~1, ~1, data = umf2, mixture = "P")
m1 <- gmultmix( ~1, ~1, ~time, data = umf2, mixture = "P")
m1b <- gmultmix( ~time, ~1, ~1, data = umf2, mixture = "P")
m1c <- gmultmix( ~time, ~1, ~time, data = umf2, mixture = "P")

# Models with date
m2 <- gmultmix( ~1, ~1, ~date, data = umf2, mixture = "P")
m2b <- gmultmix( ~date, ~1, ~1, data = umf2, mixture = "P")
m2c <- gmultmix( ~date, ~1, ~date, data = umf2, mixture = "P")

# Models with both time and date .... all run in no time
m3 <- gmultmix( ~1, ~1, ~time + date, data = umf2, mixture = "P")
# ~~~~ extra models not shown in book ~~~~~~~~~~~~~~~~~~~~~
m3b <- gmultmix( ~time, ~1, ~date, data=umf2, mixture = "P")
m3c <- gmultmix( ~time, ~1, ~time + date, data=umf2, mixture = "P")
m3d <- gmultmix( ~date,  ~1, ~time, data=umf2, mixture = "P")
m3e <- gmultmix( ~time + date, ~1, ~1, data=umf2, mixture = "P")
m3f <- gmultmix( ~time + date, ~1, ~time, data=umf2, mixture = "P")
m3g <- gmultmix( ~date,  ~1, ~time + date, data=umf2, mixture = "P")
m3h <- gmultmix( ~time + date, ~1, ~date, data=umf2, mixture = "P")
m3i <- gmultmix( ~time + date, ~1, ~time + date, data=umf2, mixture = "P")

( fl <- modSel(fitList(m0, m1, m1b, m1c, m2, m2b, m2c, m3,
    m3b, m3c, m3d, m3e, m3f, m3g, m3h, m3i)) )
#     nPars    AIC delta  AICwt cumltvWt
# m2c     4 597.88  0.00 0.4100     0.41
# m3h     5 598.61  0.73 0.2800     0.69
# m3g     5 599.70  1.83 0.1600     0.85
# m3i     6 600.34  2.47 0.1200     0.97
# m2b     3 606.31  8.43 0.0060     0.98
# m3b     4 606.63  8.75 0.0051     0.98
# m2      3 606.65  8.77 0.0051     0.99
# ... truncated ...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m2c

# Abundance:
# Estimate SE z P(>|z|)
# (Intercept) -0.5061 0.11643 -4.35 1.38e-05
# date -0.0258 0.00785 -3.28 1.04e-03

# Detection:
# Estimate SE z P(>|z|)
# (Intercept) 0.7122 0.1627 4.38 1.21e-05
# date -0.0374 0.0118 -3.17 1.51e-03

# Models with habitat covariates. ART < 1 min
m3g1 <- gmultmix( ~date + woody, ~1, ~date, data = umf2, mixture = "P")
m3g2 <- gmultmix( ~date + struct, ~1, ~date, data = umf2, mixture = "P")
m3g3 <- gmultmix( ~date + woody + struct, ~1, ~date, data = umf2, mixture = "P")

# ~~~~ model selection not shown ~~~~~~~~~~~
(fl <- modSel(fitList(m0, m1, m1b, m1c, m2, m2b, m2c, m3, m3g1, m3g2, m3g3)) )
#      nPars    AIC   delta   AICwt cumltvWt
# m3g3     6 576.56  0.0000 5.0e-01      0.5
# m3g1     5 576.56  0.0047 5.0e-01      1.0
# m3g2     5 592.40 15.8383 1.8e-04      1.0
# m2c      4 597.88 21.3200 1.2e-05      1.0
# m3h      5 598.61 22.0527 8.1e-06      1.0
# ... output truncated ...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Top model
m3g3

# Abundance:
# Estimate SE z P(>|z|)
# (Intercept) -1.6709 0.31285 -5.34 9.26e-08
# date -0.0259 0.00807 -3.21 1.33e-03
# woody 2.1163 0.50174 4.22 2.47e-05
# struct 0.0403 0.02793 1.44 1.49e-01

# Detection:
# Estimate SE z P(>|z|)
# (Intercept) 0.7156 0.1621 4.41 1.02e-05
# date -0.0371 0.0117 -3.16 1.60e-03

# Goodness-of-fit analysis (ART 7 min)
(pb <- parboot(m3g3, fitstats, nsim = 1000, ncores = 3))

# Parametric Bootstrap Statistics:

# t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
# SSE 86.2 -10.30 12.2 0.803
# Chisq 1083.3 41.84 75.9 0.285
# freemanTukey 114.0 -1.42 10.3 0.564

# t_B quantiles:
# 0% 2.5% 25% 50% 75% 97.5% 100%
# SSE 61 74 89 96 104 122 135
# Chisq 819 909 986 1037 1090 1195 1323
# freemanTukey 81 94 108 116 122 135 149

# t0 = Original statistic compuated from data
# t_B = Vector of bootstrap samples


# 2.7.3 Temporary emigration models in BUGS
# -----------------------------------------

# Remember, Y is the 3-d flycatcher data created in Section 2.7.2
y3d <- array(NA, dim = c(nrow(Y), 7, 3) ) # Create 3d array and fill it
y3d[,,1] <- Ywide[,1:7]
y3d[,,2] <- Ywide[,8:14]
y3d[,,3] <- Ywide[,15:21]
nseasons <- 3                   # Number of primary occasions
nsites <- nrow(Ywide)           # Number of sites
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion

# Bundle data
str(bdata <- list(y3d = y3d, nsites = nsites, nseasons = nseasons,
    nobs = nobs, woody = alfl.covs[,"woody"], struct = alfl.covs[,"struct"],
    date = date, time = time, pi = pi))
# List of 9
# $ y3d     : int [1:50, 1:7, 1:3] 0 0 0 0 0 0 0 0 1 0 ...
# $ nsites  : int 50
# $ nseasons: num 3
# $ nobs    : int [1:50, 1:3] 1 2 0 1 0 3 0 0 1 0 ...
# $ woody   : num [1:50] 0.3 0.05 0.35 0.3 0.1 0.4 0.2 0 0.2 0.15 ...
# $ struct  : num [1:50] 5.45 4.75 14.7 5.05 4.15 9.75 9.6 15.7 9.2 ...
# $ date    : num [1:50, 1:3] -21 -7 -7 -7 -19 -19 -19 -26 -26 -26 ...
# $ time    : num [1:50, 1:3] 1.12 1.87 0.69 0.21 2.01 ...
# $ pi      : num 3.14

# Specify model in BUGS language
cat(file="CR_TE.txt", "
model {

  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  # Availability parameters
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)
  gamma2 ~ dnorm(0, 0.01)
  # Detection parameters
  alpha0 ~ dnorm(0, 0.01)
  alpha1 ~ dnorm(0, 0.01)

  # Likelihood
  for (i in 1:nsites) {
    # Linear model for lambda
    log(lambda[i]) <- beta0 + beta1*woody[i] + beta2*struct[i]
    for(k in 1:nseasons){
      # Linear models for phi and p
      logit(phi[i,k]) <- gamma0 + gamma1*date[i,k] + gamma2*time[i,k]
      logit(p[i,k]) <- alpha0 + alpha1*date[i,k]
      # Define multinomial cell probabilities
      cp[i,k,1] <- (1-p[i,k])*(1-p[i,k])*p[i,k]
      cp[i,k,2] <- (1-p[i,k])*p[i,k]*(1-p[i,k])
      cp[i,k,3] <- (1-p[i,k])*p[i,k]*p[i,k]
      cp[i,k,4] <- p[i,k]*(1-p[i,k])*(1-p[i,k])
      cp[i,k,5] <- p[i,k]*(1-p[i,k])*p[i,k]
      cp[i,k,6] <- p[i,k]*p[i,k]*(1-p[i,k])
      cp[i,k,7] <- p[i,k]*p[i,k]*p[i,k]
      cp[i,k,8] <- 1-sum(cp[i,k,1:7])
      cellprobs.cond[i,k,1] <- cp[i,k,1]/ sum(cp[i,k,1:7])
      cellprobs.cond[i,k,2] <- cp[i,k,2]/ sum(cp[i,k,1:7])
      cellprobs.cond[i,k,3] <- cp[i,k,3]/ sum(cp[i,k,1:7])
      cellprobs.cond[i,k,4] <- cp[i,k,4]/ sum(cp[i,k,1:7])
      cellprobs.cond[i,k,5] <- cp[i,k,5]/ sum(cp[i,k,1:7])
      cellprobs.cond[i,k,6] <- cp[i,k,6]/ sum(cp[i,k,1:7])
      cellprobs.cond[i,k,7] <- cp[i,k,7]/ sum(cp[i,k,1:7])
      # Conditional 4-part version of the model
      pdet[i,k] <- sum(cp[i,k, 1:7])
      pmarg[i,k] <- pdet[i,k]*phi[i,k]
      #Part 4: multinomial
      y3d[i,1:7,k] ~ dmulti(cellprobs.cond[i,k,1:7], nobs[i,k])
      #Part 3: Number of detected individuals
      nobs[i,k] ~ dbin(pmarg[i,k], M[i])
      #Part 2: Number of available individuals
      Navail[i,k] ~ dbin(phi[i,k], M[i])
    }
    M[i] ~ dpois(lambda[i]) #Part 1: Abundance model
  }

  # Derived quantities
  for(k in 1:nseasons){
    Davail[k] <- mean(phi[,k])*exp(beta0)/ area
  }
  Mtotal <- sum(M[])
  area <- (pi*50*50)/ 10000 # 50 m point counts, D in units per ha
  Dtotal <- exp(beta0)/ area
}
")

# Initial values
Navail.st <- apply(y3d, c(1, 3), sum)
Mst <- apply(Navail.st, 1, max, na.rm = TRUE) + 2
inits <- function() list(M = Mst, sigma = 100)

# Parameters monitored
params <- c("beta0" , "beta1", "beta2", "alpha0", "alpha1", "gamma0",
    "gamma1", "gamma2", "Mtotal", "Davail", "Dtotal")

# MCMC settings
na <- 1000 ; ni <- 10000 ; nb <- 5000 ; nt <- 5 ; nc <- 3

# Call JAGS (ART 0.2 min), gauge convergence, summarize posteriors
out7 <- jags(bdata, inits, params, "CR_TE.txt", n.adapt = na, n.iter = ni,
    n.burnin = nb, n.thin = nt, n.chains = nc, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out7)
print(out7, 3)
#             mean     sd   2.5%    50%   97.5% overlap0     f  Rhat n.eff
# beta0     -0.797  0.429 -1.637 -0.797   0.047     TRUE 0.966 1.003   849
# beta1      1.964  0.628  0.735  1.955   3.212    FALSE 0.998 1.001  1914
# beta2      0.054  0.038 -0.026  0.055   0.126     TRUE 0.924 1.004   466
# alpha0     0.716  0.161  0.397  0.715   1.033    FALSE 1.000 1.000  3000
# alpha1    -0.038  0.012 -0.061 -0.038  -0.016    FALSE 1.000 1.000  3000
# gamma0    -0.351  0.300 -0.985 -0.337   0.204     TRUE 0.877 1.014   164
# gamma1    -0.059  0.015 -0.091 -0.059  -0.031    FALSE 1.000 1.005   397
# gamma2    -0.415  0.149 -0.706 -0.412  -0.125    FALSE 0.997 1.001  1535
# Mtotal    80.120 10.959 66.000 78.000 106.000    FALSE 1.000 1.035    99
# Davail[1]  0.393  0.176  0.157  0.362   0.848    FALSE 1.000 1.004   753
# Davail[2]  0.287  0.128  0.117  0.264   0.618    FALSE 1.000 1.004   879
# Davail[3]  0.162  0.076  0.060  0.148   0.353    FALSE 1.000 1.002  1410
# Dtotal     0.629  0.280  0.248  0.574   1.334    FALSE 1.000 1.003  1918
