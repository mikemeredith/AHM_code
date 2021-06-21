#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8 : MODELING INTERACTIONS AMONG SPECIES
# ===============================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 9 mins
# Run time with the full number of iterations: 9 hrs

library(AHMbook)
library(abind)
library(jagsUI)
library(corrplot)

# ~~~ Need to prepare data as in 8.3.1 ~~~~~~~~~~~~~
data(HubbardBrook)
speclist <- dimnames(HubbardBrook$counts)[[4]]
year <- 2009:2018
counts <- HubbardBrook$counts[,,11:20,] # Counts
dates <- HubbardBrook$dates[,,11:20] # Survey dates
hours <- HubbardBrook$times[,,11:20] # Survey hour (in hours)
sel.species <- c(3:4, 6:7, 12:13) # Select those 6 species
speclist[sel.species] # Check them out
speclist <- speclist[sel.species] # Update species list
counts <- counts[,,,sel.species]
DATES <- standardize(dates)
DATES[is.na(DATES)] <- 0
HOURS <- standardize(hours)
HOURS[is.na(HOURS)] <- 0
elev <- HubbardBrook$sitecov$Elev
aspect <- HubbardBrook$sitecov$Aspect
elev <- (elev - 500) / 100 # center on 500 and scale by 100
north <- abs(HubbardBrook$sitecov$Aspect-180)/180 * pi
nsites <- dim(counts)[1]
nreps <- dim(counts)[2]
nyears <- dim(counts)[3]
nspec <- dim(counts)[4]
maxC <- apply(counts, c(1,3,4), max, na.rm = TRUE)
maxC[maxC == -Inf] <- NA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 8.5 Joint models for abundance
# ==============================

# 8.5.1 Static binomial N-mixture model with few species and
#       symmetric interactions
# -----------------------------------------------------------

Rmat <- diag(nspec)      # Identity matrix
df <- nspec + 1

# Bundle and summarize data set
str(bdata <- list(C = counts[,,'2016',], nsites = nsites, nspec = nspec,
    nreps = nreps, elev = elev, north = north, DATES = DATES[,,'2016'],
    HOURS = HOURS[,,'2016'], R = Rmat, df = df) )
# List of 10
# $ C     : num [1:373, 1:3, 1:6] 0 0 0 0 0 0 0 0 0 0 ...
# $ nsites: num 373
# $ nspec : num 6
# $ nreps : num 3
# $ elev  : num [1:373] -0.26 -0.08 0.46 0.72 0.88 0.78 0.74 0.74 ...
# $ north : num [1:373] 0.1107 0.3218 0.0907 0.0588 1.9757 ...
# $ DATES : num [1:373, 1:3] -1.45 -1.45 -1.45 -1.45 -1.45 ...
# $ HOURS : num [1:373, 1:3] -1.969 -1.778 -1.464 -1.218 -0.958 ...
# $ R     : num [1:6, 1:6] 1 0 0 0 0 0 0 1 0 0 ...
# $ df    : num 7

# Specify model in BUGS language
cat(file = "Nmix1.txt", "
model {
  # Priors
  # Intercepts and coefficients all fixed effects
  for(k in 1:nspec){
    mean.lambda[k] <- exp(beta0[k])
    beta0[k] ~ dnorm(0, 0.1)
    alpha0[k] <- logit(mean.p[k])
    mean.p[k] ~ dunif(0,1)
    beta1[k] ~ dnorm(0, 0.1)
    beta2[k] ~ dnorm(0, 0.1)
    alpha1[k] ~ dnorm(0, 0.1)
    alpha2[k] ~ dnorm(0, 0.1)
    alpha3[k] ~ dnorm(0, 0.1)
  }

  # Specify MVN prior for random site effects in lambda for each species
  for (i in 1:nsites){
    eta.lam[i,1:nspec] ~ dmnorm(mu.eta[], Omega[,])
  }
  for (k in 1:nspec){
    mu.eta[k] <- 0
  }

  # Vague inverse Wishart prior for variance-covariance matrix
  Omega[1:nspec,1:nspec] ~ dwish(R[,], df)
  Sigma2[1:nspec,1:nspec] <- inverse(Omega[,])

  # Scale var/covar matrix to become the correlation matrix
  for (i in 1:nspec){
    for (k in 1:nspec){
      rho[i,k] <- Sigma2[i,k] / (sqrt(Sigma2[i,i]) * sqrt(Sigma2[k,k]))
    }
  }

  # Likelihood
  # Ecological model for true abundance
  for (i in 1:nsites){
    for(k in 1:nspec){
      N[i,k] ~ dpois(lambda[i,k])
      log(lambda[i,k]) <- beta0[k] + beta1[k] * elev[i] +
          beta2[k] * north[i] + eta.lam[i,k]

      # Observation model for replicated counts
      for (j in 1:nreps){
        C[i,j,k] ~ dbin(p[i,j,k], N[i,k])
        logit(p[i,j,k]) <- alpha0[k] + alpha1[k] * DATES[i,j] +
            alpha2[k] * pow(DATES[i,j],2) + alpha3[k] * HOURS[i,j]
      }
    }
  }

  # Derived quantities
  for(k in 1:nspec){
    Ntot[k] <- sum(N[,k])         # Total abundance per species
  }
}
")

# Initial values
Nst <- maxC[,'2016',] + 1
Nst[is.na(Nst)] <- 1
inits <- function() list(N = Nst, mean.p = rep(0.2, 6), Omega = diag(6))

# Parameters monitored
params <- c('mean.lambda', 'beta1', 'beta2', 'mean.p', 'alpha1', 'alpha2',
    'alpha3', 'Sigma2', 'rho', 'Ntot')

# MCMC settings
# na <- 5000 ; nc <- 3 ; ni <- 300000 ; nb <- 100000 ; nt <- 200
na <- 500 ; nc <- 3 ; ni <- 3000 ; nb <- 1000 ; nt <- 2  # ~~~ for testing, 3 mins

# Call JAGS (ART 400 min), assess convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "Nmix1.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
which(out1$summary[,8] > 1.1)
print(out1, 2)

# Get residual correlations
tmp2 <- round(out1$mean$rho, 2)
dimnames(tmp2) <- list(speclist, speclist)
tmp2             # This is with adjustment for two environmental covariates

# ~~~ extra code for figure 8.17 ~~~~~~~~~~~
op <- par(mar = c(10, 4, 4, 3))
xlabnames <- c(paste(speclist[1], speclist[2:6], sep='-'),
    paste(speclist[2], speclist[3:6], sep='-'),
    paste(speclist[3], speclist[4:6], sep='-'),
    paste(speclist[4], speclist[5:6], sep='-'),
    paste(speclist[5], speclist[6], sep='-'))
# model 1
mns1 <- c(out1$mean$rho[1, 2:6], out1$mean$rho[2, 3:6], out1$mean$rho[3, 4:6],
    out1$mean$rho[4, 5:6], out1$mean$rho[5, 6])
lcl1 <- c(out1$q2.5$rho[1, 2:6], out1$q2.5$rho[2, 3:6], out1$q2.5$rho[3, 4:6],
    out1$q2.5$rho[4, 5:6], out1$q2.5$rho[5, 6])
ucl1 <- c(out1$q97.5$rho[1, 2:6], out1$q97.5$rho[2, 3:6],
    out1$q97.5$rho[3, 4:6], out1$q97.5$rho[4, 5:6], out1$q97.5$rho[5, 6])
plot(1:15, mns1, xlab = '', ylab = 'Residual correlation', ylim = c(-1,1),
    pch = 16, cex = 2, frame = FALSE, axes = FALSE)
axis(1, at = 1:15, labels = xlabnames, las = 2)
axis(2)
segments(1:15, lcl1, 1:15, ucl1)
abline(h = 0)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 8.5.2 Static binomial N-mixture model with two species and
#       directional interactions
# ----------------------------------------------------------

# Select REVI and BTNW data from 2016
str(CR <- counts[, , '2016', 'REVI'])
str(CB <- counts[, , '2016', 'BTNW'])

# Bundle and summarize data set
str(bdata <- list(CR = CR, CB = CB, nsites = nsites, nreps = nreps, elev = elev,
    north = north, DATES = DATES[,,'2016'], HOURS = HOURS[,,'2016']) )
# List of 8
# $ CR    : num [1:373, 1:3] 0 1 1 2 1 1 1 1 0 0 ...
# $ CB    : num [1:373, 1:3] 0 0 0 1 0 0 0 1 0 0 ...
# $ nsites: num 373
# $ nreps : num 3
# $ elev  : num [1:373] -0.26 -0.08 0.46 0.72 0.88 0.78 0.74 0.74 ...
# $ north : num [1:373] 0.1107 0.3218 0.0907 0.0588 1.9757 ...
# $ DATES : num [1:373, 1:3] -1.45 -1.45 -1.45 -1.45 -1.45 ...
# $ HOURS : num [1:373, 1:3] -1.969 -1.778 -1.464 -1.218 -0.958 ...

# Specify model in BUGS language
cat(file = "Nmix2.txt", "
model {

  # Model for Red-eyed Vireo (REVI): the 'dominant' species
  # Priors for intercepts and coefficients
  # mean.lambdar ~ dunif(0, 5)
  mean.lambdaR ~ dunif(0, 5)
  beta0R <- log(mean.lambdaR)
  alpha0R <- logit(mean.pR)
  mean.pR ~ dunif(0,1)
  beta1R ~ dnorm(0, 0.1)
  beta2R ~ dnorm(0, 0.1)
  alpha1R ~ dnorm(0, 0.1)
  alpha2R ~ dnorm(0, 0.1)
  alpha3R ~ dnorm(0, 0.1)

  # Likelihood for REVI ('dominant')
  # Ecological model
  for (i in 1:nsites){
    NR[i] ~ dpois(lambdaR[i])
    log(lambdaR[i]) <- beta0R + beta1R * elev[i] + beta2R * north[i]
    # Observation model
    for (j in 1:nreps){
      CR[i,j] ~ dbin(pR[i,j], NR[i])
      logit(pR[i,j]) <- alpha0R + alpha1R * DATES[i,j] +
      alpha2R * pow(DATES[i,j],2) + alpha3R * HOURS[i,j]
    }
  }

  # Model for Black-throated Green Warbler (BTNW): the 'subordinate' sp.
  # Priors for intercepts and coefficients
  mean.lambdaB ~ dunif(0, 5)
  beta0B <- log(mean.lambdaB)
  beta1B ~ dnorm(0, 0.1)
  beta2B ~ dnorm(0, 0.1)
  alpha0B <- logit(mean.pB)
  mean.pB ~ dunif(0,1)
  alpha1B ~ dnorm(0, 0.1)
  alpha2B ~ dnorm(0, 0.1)
  alpha3B ~ dnorm(0, 0.1)
  gamma0 ~ dnorm(0, 0.1) # These are the 'interaction coefficients' !
  gamma1 ~ dnorm(0, 0.1)
  gamma2 ~ dnorm(0, 0.1)

  # Likelihood for BTNW ('subordinate')
  # Ecological model
  for (i in 1:nsites){
    NB[i] ~ dpois(lambdaB[i])
    lambdaB[i] <- exp(beta0B + beta1B * elev[i] + beta2B * north[i] +
    gamma0 * NR[i] + gamma1 * elev[i] * NR[i] + gamma2 * north[i] * NR[i])
    # Observation model
    for (j in 1:nreps){
      CB[i,j] ~ dbin(pB[i,j], NB[i])
      logit(pB[i,j]) <- alpha0B + alpha1B * DATES[i,j] +
      alpha2B * pow(DATES[i,j],2) + alpha3B * HOURS[i,j]
    }
  }
}
")

# Initial values
Nst <- rep(10, nsites)
inits <- function() list(NR = Nst, NB = Nst)

# Parameters monitored
params <- c('mean.lambdaR', 'beta1R', 'beta2R', 'mean.pR', 'alpha1R',
    'alpha2R', 'alpha3R', 'mean.lambdaB', 'beta1B', 'beta2B', 'mean.pB',
    'alpha1B', 'alpha2B', 'alpha3B', 'gamma0', 'gamma1', 'gamma2')

# MCMC settings
# na <- 5000 ; nc <- 3 ; ni <- 100000 ; nb <- 50000 ; nt <- 50
na <- 500 ; nc <- 3 ; ni <- 10000 ; nb <- 5000 ; nt <- 5  # ~~~ for testing, 3 mins

# Call JAGS (ART 40 min), assess convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "Nmix2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out2)
which(out2$summary[,8] > 1.1)
print(out2, 2)
#               mean   sd  2.5%   50% 97.5% overlap0    f Rhat n.eff
# mean.lambdaR  2.45 0.62  1.53  2.36  3.94    FALSE 1.00 1.00  1314
# ... output truncated ...
# gamma0        0.09 0.10 -0.12  0.10  0.27     TRUE 0.85 1.02   110
# gamma1       -0.01 0.03 -0.07 -0.01  0.04     TRUE 0.70 1.00   528
# gamma2        0.03 0.04 -0.05  0.03  0.12     TRUE 0.76 1.00   702

# 8.5.3 An open binomial N-mixture for two species with directional interactions
# ------------------------------------------------------------------------------

# Select all REVI and BTNW data
str(CR <- counts[,,, 'REVI'])
str(CB <- counts[,,, 'BTNW'])
table(CR) ; table(CB) # Descriptifs of these counts

# Bundle and summarize data set
str(bdata <- list(CR = CR, CB = CB, nsites = nsites, nyears = nyears,
    nreps = nreps, elev = elev, north = north, DATES = DATES, HOURS = HOURS) )
# List of 9
# $ CR    : num [1:373, 1:3, 1:10] 2 1 0 1 1 0 1 0 0 1 ...
# $ CB    : num [1:373, 1:3, 1:10] 1 0 1 0 0 0 0 1 1 0 ...
# $ nsites: num 373
# $ nyears: num 10
# $ nreps : num 3
# $ elev  : num [1:373] -0.26 -0.08 0.46 0.72 0.88 0.78 0.74 0.74 ...
# $ north : num [1:373] 0.1107 0.3218 0.0907 0.0588 1.9757 ...
# $ DATES : num [1:373, 1:3, 1:10] -1.36 -1.36 -1.36 -1.36 -1.36 ...
# $ HOURS : num [1:373, 1:3, 1:20] -1.655 -1.45 -1.109 -0.753 -1.204 ...

# Specify model in BUGS language
cat(file = "Nmix3.txt", "
model {

  # Model for Red-eyed Vireo (REVI): the 'dominant' species
  # Priors for intercepts and coefficients
  for(t in 1:nyears){
    mean.lambdaR[t] ~ dunif(0, 10)
    beta0R[t] <- log(mean.lambdaR[t])
    alpha0R[t] <- logit(mean.pR[t])
    mean.pR[t] ~ dunif(0,1)
  }
  beta1R ~ dnorm(0, 0.1)
  beta2R ~ dnorm(0, 0.1)
  alpha1R ~ dnorm(0, 0.1)
  alpha2R ~ dnorm(0, 0.1)
  alpha3R ~ dnorm(0, 0.1)

  # Likelihood for REVI ('dominant')
  # Ecological model
  for(t in 1:nyears){
    for (i in 1:nsites){
      NR[i,t] ~ dpois(lambdaR[i,t])
      lambdaR[i,t] <- exp(beta0R[t] + beta1R * elev[i] + beta2R * north[i])
      # Observation model
      for (j in 1:nreps){
        CR[i,j,t] ~ dbin(pR[i,j,t], NR[i,t])
        logit(pR[i,j,t]) <- alpha0R[t] + alpha1R * DATES[i,j,t] +
        alpha2R * pow(DATES[i,j,t],2) + alpha3R * HOURS[i,j,t]
      }
    }
  }

  # Model for Black-throated Green Warbler (BTNW): the 'subordinate' sp.
  # Priors for intercepts and coefficients
  for(t in 1:nyears){
    mean.lambdaB[t] ~ dunif(0, 10)
    beta0B[t] <- log(mean.lambdaB[t])
    alpha0B[t] <- logit(mean.pB[t])
    mean.pB[t] ~ dunif(0,1)
  }
  beta1B ~ dnorm(0, 0.1)
  beta2B ~ dnorm(0, 0.1)
  alpha1B ~ dnorm(0, 0.1)
  alpha2B ~ dnorm(0, 0.1)
  alpha3B ~ dnorm(0, 0.1)
  gamma ~ dnorm(0, 0.1)           # This is the 'interaction coefficient'

  # Likelihood for BTNW ('subordinate')
  # For year 1 = 2009: no lagged effect of NR
  for (i in 1:nsites){
    NB[i,1] ~ dpois(lambdaB[i,1])
    lambdaB[i,1] <- exp(beta0B[1] + beta1B * elev[i] + beta2B * north[i])
    for (j in 1:nreps){
      CB[i,j,1] ~ dbin(pB[i,j,1], NB[i,1])
      logit(pB[i,j,1]) <- alpha0B[1] + alpha1B * DATES[i,j,1] +
          alpha2B * pow(DATES[i,j,1],2) + alpha3B * HOURS[i,j,1]
    }
    # For years 2010-2018: including a lagged effect of NR
    for(t in 2:nyears){
      NB[i,t] ~ dpois(lambdaB[i,t])
      lambdaB[i,t] <- exp(beta0B[t] + beta1B * elev[i] +
          beta2B * north[i] + gamma * NR[i,t-1])
      for (j in 1:nreps){
        CB[i,j,t] ~ dbin(pB[i,j,t], NB[i,t])
        logit(pB[i,j,t]) <- alpha0B[t] + alpha1B * DATES[i,j,t] +
            alpha2B * pow(DATES[i,j,t],2) + alpha3B * HOURS[i,j,t]
      }
    }
  }
}
")

# Cheap initial values: make sure N is initialized high enough
Nst <- array(10, dim = c(nsites, nyears))
inits <- function() list(NR = Nst, NB = Nst)

# Parameters monitored
params <- c('mean.lambdaR', 'beta1R', 'beta2R', 'mean.pR', 'alpha1R', 'alpha2R',
    'alpha3R', 'mean.lambdaB', 'beta1B', 'beta2B', 'mean.pB', 'alpha1B',
    'alpha2B', 'alpha3B', 'gamma')

# MCMC settings
# na <- 5000 ; nc <- 3 ; ni <- 100000 ; nb <- 50000 ; nt <- 50
na <- 500 ; nc <- 3 ; ni <- 1000 ; nb <- 500 ; nt <- 1  # ~~~ for testing

# Call JAGS (ART 509 min), assess convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "Nmix3.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out3)
which(out3$summary[,8] > 1.1)
print(out3$summary[,-(4:6)], 3)
#                    mean     sd     2.5%  97.5% Rhat n.eff overlap0     f
# mean.lambdaR[1] 2.29138 0.4938  1.55938 3.4651    1  3000        0 1.000
# ... output heavily truncated ...
# gamma           0.02872 0.0153 -0.00229 0.0588    1  2586        1 0.964
