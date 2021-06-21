#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8 : MODELING INTERACTIONS AMONG SPECIES
# ===============================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 55 mins
# Run time with the full number of iterations: 4 hrs

library(AHMbook)
library(jagsUI)

# 8.3 Joint occupancy models for few species with directional interactions
# ========================================================================

# 8.3.1 The Hubbard Brook data set
# --------------------------------

library(AHMbook)
data(HubbardBrook)
str(HubbardBrook)

# Get four-letter codes for the 13 species
speclist <- dimnames(HubbardBrook$counts)[[4]]

# Restrict all data to last 10 years: 2009-2018
year <- 2009:2018
str(counts <- HubbardBrook$counts[,,11:20,]) # Counts
str(dates <- HubbardBrook$dates[,,11:20])    # Survey dates
str(hours <- HubbardBrook$times[,,11:20])    # Survey hour (in hours)

# Date of full canopy leaf expansion as a phenological measure
budburst <- c(138.4, 133.8, 143.6, 133.9, 137.3, 146.6, 142.0, 148.3, 143.3, 141.0)
plot(year, budburst, xlab = 'Year', ylab = 'Day of budburst',
    ylim = c(130, 150), type = 'b', pch = 16, cex = 2) # not shown

# Get subset of six of the more common species
old.speclist <- speclist               # Make copy of original
sel.species <- c(3:4, 6:7, 12:13)      # Select those 6 species
speclist[sel.species]                  # Check them out
speclist <- speclist[sel.species]      # Update species list

# Restrict count data to the selected six species
str(counts <- counts[,,,sel.species])

# Scale and mean-impute both survey dates and survey hour (times)
DATES <- standardize(dates)
DATES[is.na(DATES)] <- 0
HOURS <- standardize(hours)
HOURS[is.na(HOURS)] <- 0

# Grab and scale elevation and aspect covariates
hist(elev <- HubbardBrook$sitecov$Elev, col = 'gold')
hist(aspect <- HubbardBrook$sitecov$Aspect, col = 'gold')
elev <- (elev - 500) / 100                 # center on 500 and scale by 100
north <- abs(HubbardBrook$sitecov$Aspect-180)/180 * pi

# Get sample sizes
nsites <- dim(counts)[1]
nreps <- dim(counts)[2]
nyears <- dim(counts)[3]
nspec <- dim(counts)[4]

# Look at relative abundance and occurrence patterns of those species
maxC <- apply(counts, c(1,3,4), max, na.rm = TRUE)
maxC[maxC == '-Inf'] <- NA
obsz <- maxC
obsz[obsz > 1] <- 1           # Observed p/a array for site, year and species

# Number of surveyed plots per year (for last 10 years)
tmp <- maxC[,,1]
tmp[!is.na(tmp)] <- 1
(nplots.per.year <- apply(tmp, 2, sum, na.rm = TRUE))

# Observed proportion of occupied sites per year
nocc.obs <- apply(obsz, c(2,3), sum, na.rm = TRUE)
psi.obs <- nocc.obs / nplots.per.year

# Plot observed proportion of occupied sites (points) (Fig. 8.11)
op <- par(cex.lab = 1.5, cex.axis = 1.5)
matplot(year, psi.obs, type = 'b', lty = 1,
    main = 'Observed occupancy per year and species', pch = 1:6, col = 1:6,
    xlab = 'Year', ylab = 'Observed occupancy', frame = FALSE, las = 1,
    ylim = c(0, 1.15), cex = 1.5, lwd = 3)
legend('topleft', legend = speclist, pch = 1:6, lwd = 2, col = 1:6,
    bty = 'n', horiz = TRUE, cex = 1.5)
par(op)

# 8.3.2 Directional dependence in a static occupancy model
# --------------------------------------------------------

# Quantize counts to become detection/nondetection data
y <- counts ; y[y > 1] <- 1

# Pull out REVI and BTNW data from 2016
str(yR <- y[, , '2016', 'REVI'])
str(yB <- y[, , '2016', 'BTNW'])

# Bundle and summarize data set
str(bdata <- list(yR = yR, yB = yB, nsites = nsites, nreps = nreps,
    elev = elev, north = north, DATES = DATES[,,'2016'],
    HOURS = HOURS[,,'2016']))
# List of 8
# $ yR   : num [1:373, 1:3] 0 1 1 1 1 1 1 1 0 0 ...
# $ yB   : num [1:373, 1:3] 0 0 0 1 0 0 0 1 0 0 ...
# $ nsite: num 373
# $ nrep : num 3
# $ elev : num [1:373] -0.26 -0.08 0.46 0.72 0.88 0.78 0.74 0.74 0.85 ...
# $ north: num [1:373] 0.1107 0.3218 0.0907 0.0588 1.9757 ...
# $ DATES: num [1:373, 1:3] -1.45 -1.45 -1.45 -1.45 -1.45 ...
# $ HOURS: num [1:373, 1:3] -1.969 -1.778 -1.464 -1.218 -0.958 ...

# Specify model in BUGS language
cat(file = "staticWaddle.txt", "
model {

  # Model for Red-eyed vireo (REVI): the 'dominant' species
  # Priors for intercepts and coefficients
  mean.psiR <- ilogit(beta0R)
  beta0R ~ dnorm(0, 0.1)
  alpha0R ~ dnorm(0, 0.1)
  mean.pR <- ilogit(alpha0R)
  beta1R ~ dnorm(0, 0.1)
  beta2R ~ dnorm(0, 0.1)
  alpha1R ~ dnorm(0, 0.1)
  alpha2R ~ dnorm(0, 0.1)
  alpha3R ~ dnorm(0, 0.1)

  # Likelihood for REVI ('dominant')
  # Ecological model
  for (i in 1:nsites){
    zR[i] ~ dbern(psiR[i])
    logit(psiR[i]) <- beta0R + beta1R * elev[i] + beta2R * north[i]
    # Observation model
    for (j in 1:nreps){
      yR[i,j] ~ dbern(pR[i,j] * zR[i])
      logit(pR[i,j]) <- alpha0R + alpha1R * DATES[i,j] + alpha2R *
      pow(DATES[i,j],2) + alpha3R * HOURS[i,j]
    }
  }

  # Model for Black-throated Green Warbler (BTNW): the 'subordinate' sp.
  # Priors for intercepts and coefficients
  for(k in 1:2){ # Note all are stratified by presence/absence of REVI !
    mean.psiB[k] <- ilogit(beta0B[k])
    beta0B[k] ~ dnorm(0, 0.1)
    beta1B[k] ~ dnorm(0, 0.1)
    beta2B[k] ~ dnorm(0, 0.1)
    alpha0B[k] ~ dnorm(0, 0.1)
    mean.pB[k] <- ilogit(alpha0B[k])
    alpha1B[k] ~ dnorm(0, 0.1)
    alpha2B[k] ~ dnorm(0, 0.1)
    alpha3B[k] ~ dnorm(0, 0.1)
  }

  # Likelihood for BTNW ('subordinate')
  # Ecological model
  for (i in 1:nsites){
    zB[i] ~ dbern(psiB[i])
    logit(psiB[i]) <- beta0B[zR[i]+1] + beta1B[zR[i]+1] * elev[i] +
    beta2B[zR[i]+1] * north[i] # note nested indexing as a 'switch' !
    # Observation model
    for (j in 1:nreps){
      yB[i,j] ~ dbern(pB[i,j] * zB[i])
      logit(pB[i,j]) <- alpha0B[zR[i]+1] + alpha1B[zR[i]+1] * DATES[i,j] +
      alpha2B[zR[i]+1] * pow(DATES[i,j],2) + alpha3B[zR[i]+1] * HOURS[i,j]
    }
  }

  # Derived quantities
  # Differences in BTNW parameters between REVI presence and absence
  diffbeta0B <- beta0B[2] - beta0B[1] # 'presence of R minus absence ...'
  diffbeta1B <- beta1B[2] - beta1B[1]
  diffbeta2B <- beta2B[2] - beta2B[1]
  diffalpha0B <- alpha0B[2] - alpha0B[1]
  diffalpha1B <- alpha1B[2] - alpha1B[1]
  diffalpha2B <- alpha2B[2] - alpha2B[1]
  diffalpha3B <- alpha3B[2] - alpha3B[1]
  # Tally up number of sites in each of 4 states
  for(i in 1:nsites){
    tmpRB[i] <- zR[i] * zB[i]
    tmpRb[i] <- zR[i] * (1-zB[i])
    tmprB[i] <- (1-zR[i]) * zB[i]
    tmprb[i] <- (1-zR[i]) * (1-zB[i])
  }
  nRB <- sum(tmpRB[]) # Number of sites with both species present
  nRb <- sum(tmpRb[]) # Number of sites with only the vireo present
  nrB <- sum(tmprB[]) # Number of sites with only the warbler present
  nrb <- sum(tmprb[]) # Number of sites with both species absent
}
")

# Initial values
zst <- rep(1, nsites)
inits <- function() list(zR = zst, zB = zst)

# Parameters monitored
params <- c('mean.psiR', 'beta1R', 'beta2R', 'mean.pR', 'alpha1R', 'alpha2R',
    'alpha3R', 'mean.psiB', 'beta1B', 'beta2B', 'mean.pB', 'alpha1B', 'alpha2B',
    'alpha3B', 'diffbeta0B', 'diffbeta1B', 'diffbeta2B', 'diffalpha0B',
    'diffalpha1B', 'diffalpha2B', 'diffalpha3B', 'nRB', 'nRb', 'nrB', 'nrb')

# MCMC settings
# na <- 5000 ; nc <- 3 ; ni <- 30000 ; nb <- 10000 ; nt <- 20
na <- 500 ; nc <- 3 ; ni <- 3000 ; nb <- 1000 ; nt <- 2  # ~~~ for testing, 2.2 mins

# Call JAGS (ART 39 min), assess convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "staticWaddle.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
which(out1$summary[,8] > 1.1)
print(out1, 2)
#               mean    sd   2.5%    50%  97.5% overlap0    f Rhat n.eff
# mean.psiR     0.95  0.04   0.86   0.96   1.00    FALSE 1.00 1.00  1769
# beta1R       -1.64  0.35  -2.44  -1.61  -1.04    FALSE 1.00 1.00  3000
# beta2R        0.20  0.21  -0.26   0.20   0.60     TRUE 0.83 1.00  3000
# .....
# diffbeta0B    5.74  3.53  -2.15   5.76  12.77     TRUE 0.94 1.01   684
# diffbeta1B   -2.33  1.52  -5.82  -2.20   0.40     TRUE 0.96 1.01   313
# diffbeta2B    1.74  2.44  -2.94   1.53   7.12     TRUE 0.78 1.03    86
# diffalpha0B  -0.15  0.34  -0.79  -0.15   0.55     TRUE 0.68 1.00  3000
# diffalpha1B   0.67  0.24   0.22   0.66   1.15    FALSE 1.00 1.00  3000
# diffalpha2B   0.15  0.29  -0.44   0.16   0.69     TRUE 0.71 1.00   977
# diffalpha3B   0.33  0.20  -0.05   0.33   0.71     TRUE 0.96 1.00   533
# nRB         247.96 15.25 215.00 250.00 275.00    FALSE 1.00 1.05    52
# nRb          17.02 14.34   0.00  14.00  52.00     TRUE 1.00 1.05    51
# nrB          90.19  8.71  75.00  90.00 109.00    FALSE 1.00 1.00   603
# nrb          17.82  8.93   0.00  19.00  34.00     TRUE 1.00 1.00  1993

# ~~~ extra code for figure 8.12 ~~~
tmp <- out1$mean # Grab the point estimates (posterior means)
# Predictions of elevation profile of occupancy
epred.orig <- 200:1000     # Approximate elevation range in data
epred <- (epred.orig - 500) / 100
op <- par(mfrow = c(1,2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot(epred.orig, plogis(qlogis(tmp$mean.psiB[1]) + tmp$beta1B[1] * epred),
    type = 'l', lwd = 3, lty = 2, col = 'red', frame = FALSE,
    xlab = 'Elevation (metres)', ylab = 'Predicted occupancy prob. of BTNW',
    ylim = c(0,1))
lines(epred.orig, plogis(qlogis(tmp$mean.psiB[2]) + tmp$beta1B[2] * epred),
    type = 'l', lwd = 3, col = 'blue')

# Predictions of seasonal profile of detection
dpred.orig <- 120:200   # Approximate date range in data
dpred <- (dpred.orig - mean(dates, na.rm = TRUE)) / sd(dates, na.rm = TRUE)
plot(dpred.orig, plogis(qlogis(tmp$mean.pB[1]) + tmp$alpha1B[1] * dpred +
    tmp$alpha2B[1] * dpred^2), type = 'l', lwd = 3, lty = 2, col = 'red',
    frame = FALSE, xlab = 'Survey date (Julian day)',
    ylab = 'Predicted detection prob. of BTNW', ylim = c(0,1))
lines(dpred.orig, plogis(qlogis(tmp$mean.pB[2]) + tmp$alpha1B[2] * dpred +
    tmp$alpha2B[2] * dpred^2), type = 'l', lwd = 3, col = 'blue')
legend('bottomleft', c('Red-eyed Vireo present', 'Red-eyed Vireo absent'),
    lwd = 3, col = c('blue', 'red'), lty = c(1,2), cex = 1.5, bty = 'n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 8.3.3 Directional dependence in a dynamic occupancy model
# ---------------------------------------------------------

# Grab REVI and BTNW data for all 10 years
str(yR <- y[,,, 'REVI'])
str(yB <- y[,,, 'BTNW'])

# Bundle and summarize data set
bbcov <- as.numeric(scale(budburst))
str(bdata <- list(yR = yR, yB = yB, nsites = nsites, nyears = nyears,
    nreps = nreps, elev = elev, north = north, DATES = DATES, HOURS = HOURS,
    bbcov = bbcov) )
# List of 10
# $ yR    : num [1:373, 1:3, 1:10] 1 1 0 1 1 0 1 0 0 1 ...
# $ yB    : num [1:373, 1:3, 1:10] 1 0 1 0 0 0 0 1 1 0 ...
# $ nsites: int 373
# $ nyears: int 10
# $ nreps : int 3
# $ elev  : num [1:373] -0.26 -0.08 0.46 0.72 0.88 0.78 0.74 0.74 ...
# $ north : num [1:373] 0.1107 0.3218 0.0907 0.0588 1.9757 ...
# $ DATES : num [1:373, 1:3, 1:10] -1.36 -1.36 -1.36 -1.36 -1.36 ...
# $ HOURS : num [1:373, 1:3, 1:10] -1.58 -1.378 -1.039 -0.798 -0.572 ...
# $ bbcov : num [1:10] -0.489 -1.42 0.562 -1.399 -0.712 ...

# Specify model in BUGS language
cat(file = "dynamicWaddle.txt", "
model {

  # Dynocc model for Red-eyed Vireo (REVI; 'dominant')
  # Priors for initial occupancy
  mean.psiR ~ dunif(0, 1)
  alpha.lpsiR <- logit(mean.psiR)
  # Linear models and priors for phi and gamma
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      alpha.lphiR[i,t] <- alpha0.lphiR + beta.lphiR[1] * elev[i] +
      beta.lphiR[2] * north[i] + beta.lphiR[3] * bbcov[t]
      alpha.lgammaR[i, t] <- alpha0.lgammaR + beta.lgammaR[1] * elev[i] +
      beta.lgammaR[2] * north[i] + beta.lgammaR[3] * bbcov[t]
    }
  }
  for(v in 1:3){
    beta.lphiR[v] ~ dnorm(0, 0.1)
    beta.lgammaR[v] ~ dnorm(0, 0.1)
  }
  mean.phiR ~ dunif(0,1)
  alpha0.lphiR <- logit(mean.phiR)
  mean.gammaR ~ dunif(0,1)
  alpha0.lgammaR <- logit(mean.gammaR)
  # Priors for p
  for(t in 1:nyears){
    alpha.lpR[t] <- mean.pR[t]
    mean.pR[t] ~ dunif(0, 1)
  }
  alpha1R ~ dnorm(0, 0.1)
  alpha2R ~ dnorm(0, 0.1)
  alpha3R ~ dnorm(0, 0.1)

  # Likelihood for 'REVI'
  # Ecological submodel
  for (i in 1:nsites){
    # Initial conditions of system
    zR[i,1] ~ dbern(psi1R[i])
    logit(psi1R[i]) <- alpha.lpsiR
    # State transitions
    for (t in 2:nyears){
      zR[i,t] ~ dbern(zR[i,t-1]*phiR[i,t-1] + (1-zR[i,t-1])*gammaR[i,t-1])
      logit(phiR[i,t-1]) <- alpha.lphiR[i, t-1]
      logit(gammaR[i,t-1]) <- alpha.lgammaR[i, t-1]
    }
  }

  # Observation model
  for (i in 1:nsites){
    for (j in 1:nreps){
      for (t in 1:nyears){
        yR[i,j,t] ~ dbern(zR[i,t] * pR[i,j,t])
        logit(pR[i,j,t]) <- alpha.lpR[t] + alpha1R * DATES[i,j,t] +
        alpha2R * pow(DATES[i,j,t],2) + alpha3R * HOURS[i,j,t]
      }
    }
  }

  # Dynocc model for Black-throated Green Warbler (BTNW; 'subordinate')
  # Priors for initial occupancy
  mean.psiB ~ dunif(0, 1)
  alpha.lpsiB <- logit(mean.psiB)
  # Linear models and priors for phi and gamma
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      alpha.lphiB[i,t] <- alpha0.lphiB[zR[i,t]+1] +
          beta.lphiB1[zR[i,t]+1] * elev[i] +
          beta.lphiB2[zR[i,t]+1] * north[i] + beta.lphiB3[zR[i,t]+1] * bbcov[t]
      alpha.lgammaB[i,t] <- alpha0.lgammaB[zR[i,t]+1] +
          beta.lgammaB1[zR[i,t]+1] * elev[i] +
          beta.lgammaB2[zR[i,t]+1] * north[i] +
          beta.lgammaB3[zR[i,t]+1] * bbcov[t]
    }
  }
  for(k in 1:2){ # Note here stratification by vireo presence/absence
    mean.phiB[k] ~ dunif(0,1)
    alpha0.lphiB[k] <- logit(mean.phiB[k])
    mean.gammaB[k] ~ dunif(0,1)
    alpha0.lgammaB[k] <- logit(mean.gammaB[k])
    beta.lphiB1[k] ~ dnorm(0, 0.1)
    beta.lphiB2[k] ~ dnorm(0, 0.1)
    beta.lphiB3[k] ~ dnorm(0, 0.1)
    beta.lgammaB1[k] ~ dnorm(0, 0.1)
    beta.lgammaB2[k] ~ dnorm(0, 0.1)
    beta.lgammaB3[k] ~ dnorm(0, 0.1)
  }

  # Priors for model of p
  for(t in 1:nyears){
    alpha.lpB[t] <- mean.pB[t]
    mean.pB[t] ~ dunif(0, 1)
  }
  alpha1B ~ dnorm(0, 0.1)
  alpha2B ~ dnorm(0, 0.1)
  alpha3B ~ dnorm(0, 0.1)

  # Likelihood for 'BTNW'
  # Ecological submodel
  for (i in 1:nsites){
    # Initial conditions of system
    zB[i,1] ~ dbern(psi1B[i])
    logit(psi1B[i]) <- alpha.lpsiB
    # State transitions
    for (t in 2:nyears){
      zB[i,t] ~ dbern(zB[i,t-1]*phiB[i,t-1] + (1-zB[i,t-1])*gammaB[i,t-1])
      logit(phiB[i,t-1]) <- alpha.lphiB[i,t-1]
      logit(gammaB[i,t-1]) <- alpha.lgammaB[i,t-1]
    }
  }

  # Observation model
  for (i in 1:nsites){
    for (j in 1:nreps){
      for (t in 1:nyears){
        yB[i,j,t] ~ dbern(zB[i,t] * pB[i,j,t])
        logit(pB[i,j,t]) <- alpha.lpB[t] + alpha1B * DATES[i,j,t] +
        alpha2B * pow(DATES[i,j,t],2) + alpha3B * HOURS[i,j,t]
      }
    }
  }

  # Derived quantities
  # Differences in BTNW parameters between REVI presence and absence
  diffalpha0.lphiB <- alpha0.lphiB[2] - alpha0.lphiB[1]
  diffalpha0.lgammaB <- alpha0.lgammaB[2] - alpha0.lgammaB[1]
  diffbeta.lphiB1 <- beta.lphiB1[2] - beta.lphiB1[1]
  diffbeta.lphiB2 <- beta.lphiB2[2] - beta.lphiB2[1]
  diffbeta.lphiB3 <- beta.lphiB3[2] - beta.lphiB3[1]
  diffbeta.lgammaB1 <- beta.lgammaB1[2] - beta.lgammaB1[1]
  diffbeta.lgammaB2 <- beta.lgammaB2[2] - beta.lgammaB2[1]
  diffbeta.lgammaB3 <- beta.lgammaB3[2] - beta.lgammaB3[1]

  # Tally up number of sites in each of 4 states
  for (t in 1:nyears){
    for(i in 1:nsites){
      tmpRB[i,t] <- zR[i,t] * zB[i,t]
      tmpRb[i,t] <- zR[i,t] * (1-zB[i,t])
      tmprB[i,t] <- (1-zR[i,t]) * zB[i,t]
      tmprb[i,t] <- (1-zR[i,t]) * (1-zB[i,t])
    }
    nRB[t] <- sum(tmpRB[,t]) # Number of sites with both species
    nRb[t] <- sum(tmpRb[,t]) # Number of sites with only the vireo
    nrB[t] <- sum(tmprB[,t]) # Number of sites with only the warbler
    nrb[t] <- sum(tmprb[,t]) # Number of sites with neither species
  }
}
")

# Initial values
zst <- array(1, dim = c(nsites, nyears))
inits <- function(){ list(zR = zst, zB = zst)}

# Parameters monitored
params <- c('mean.psiR', 'mean.phiR', 'mean.gammaR', 'mean.pR',
    'alpha.lpsiR', 'alpha0.lphiR', 'alpha0.lgammaR', 'beta.lphiR',
    'beta.lgammaR', 'alpha.lpR', 'alpha1R', 'alpha2R', 'alpha3R',
    'mean.psiB', 'mean.phiB', 'mean.gammaB', 'mean.pB', 'alpha.lpsiB',
    'alpha0.lphiB', 'alpha0.lgammaB', 'beta.lphiB1', 'beta.lphiB2',
    'beta.lphiB3', 'beta.lgammaB1', 'beta.lgammaB2', 'beta.lgammaB3',
    'alpha.lpB', 'alpha1B', 'alpha2B', 'alpha3B', 'diffalpha0.lphiB',
    'diffalpha0.lgammaB', 'diffbeta.lphiB1', 'diffbeta.lphiB2',
    'diffbeta.lphiB3', 'diffbeta.lgammaB1', 'diffbeta.lgammaB2',
    'diffbeta.lgammaB3', 'nRB', 'nRb', 'nrB', 'nrb')

# MCMC settings
# na <- 5000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 500 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~ for testing, 110 mins

# Call JAGS (ART 323 min), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "dynamicWaddle.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out2)
which(out2$summary[,8] > 1.1)
print(out2$summary[, -c(4:6, 10:11)], 2)
#                        mean      sd      2.5%    97.5% Rhat n.eff
# mean.psiR           0.58788 0.03053  0.527482  0.64761 1.00  3000
# mean.phiR           0.97604 0.00893  0.954922  0.98890 1.03    84
# mean.gammaR         0.62962 0.05668  0.517656  0.73867 1.00   831
# ....
# diffalpha0.lphiB    0.01327 0.67091 -1.454793  1.20533 1.03    72
# diffalpha0.lgammaB -0.84677 1.37805 -3.826917  1.41830 1.05    56
# diffbeta.lphiB1    -0.03145 0.22120 -0.480559  0.39371 1.00   514
# diffbeta.lphiB2     0.56479 0.24648  0.121004  1.08160 1.02   131
# diffbeta.lphiB3    -0.14117 0.30585 -0.706631  0.48144 1.00  3000
# diffbeta.lgammaB1   0.49208 0.39742 -0.236849  1.32898 1.03    79
# diffbeta.lgammaB2  -0.43980 0.37352 -1.189054  0.29883 1.01   240
# diffbeta.lgammaB3  -1.19232 0.50848 -2.292762 -0.30256 1.00   684
# ... output truncated ...

# ~~~~~~ extra code for figure 8.13 ~~~~~~~~~~~
tmp <- out2$mean # Grab the point estimates (posterior means)

# Predictions of aspect profile of persistence
apred <- seq(0, pi, length.out = 1000)  # Aspect range in data
op <- par(mfrow = c(1,2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot(apred, plogis(tmp$alpha0.lphiB[1] + tmp$beta.lphiB2[1] * apred),
    type = 'l', lwd = 3, lty = 2, col = 'red', frame = FALSE,
    xlab = 'Aspect (north)', ylab = 'Predicted persistence prob. of BTNW',
    ylim = c(0.6,1))
lines(apred, plogis(tmp$alpha0.lphiB[2] + tmp$beta.lphiB2[2] * apred),
    type = 'l', lwd = 3, col = 'blue')
legend('bottomleft', c('Red-eyed Vireo present', 'Red-eyed Vireo absent'),
    lwd = 3, col = c('blue', 'red'), lty = c(1,2), cex = 1.5, bty = 'n')
# Predictions of budburst date profile of colonization
bpred.orig <- 130:150   # Approximate budburst range in data
bpred <- (bpred.orig - mean(budburst[-1])) / sd(budburst[-1])
plot(bpred.orig, plogis(tmp$alpha0.lgammaB[1] + tmp$beta.lgammaB3[1] * bpred),
    type = 'l', lwd = 3, lty = 2, col = 'red', frame = FALSE,
    xlab = 'Date of full budburst',
    ylab = 'Predicted colonization prob. of BTNW',
    ylim = c(0.4,1))
lines(bpred.orig, plogis(tmp$alpha0.lgammaB[2] + tmp$beta.lgammaB3[2] * bpred),
    type = 'l', lwd = 3, col = 'blue')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ extra code for figure 8.14 ~~~~~~~~~~
# Plot estimated number of sites in 4 occupancy states (Fig. 19–17)
op <- par(mar = c(6, 6, 4, 3), cex.lab = 1.5, cex.axis = 1.5)
matplot(year, cbind(tmp$nRB, tmp$nRb, tmp$nrB, tmp$nrb), type = 'b',
    lty = 1, main = 'Number of sites per occupancy state', pch = 1:4,
    col = 1:4, xlab = 'Year', ylab = 'Number of sites', frame = FALSE,
    las = 1, ylim = c(0, 320), cex = 1.5, lwd = 3)
legend('topleft', legend = c('both', 'REVI only', 'BTNW only', 'none'),
    pch = 1:4, lwd = 2, col = 1:4, bty = 'n', horiz = TRUE, cex = 1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 8.3.4 Structural equation modeling of species interaction networks
# (no code)
