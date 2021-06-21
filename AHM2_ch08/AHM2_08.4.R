#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8 : MODELING INTERACTIONS AMONG SPECIES
# ===============================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 75 mins
# Run time with the full number of iterations: 4 days

library(AHMbook)
library(abind)
library(jagsUI)

# 8.4 Joint occupancy models for “many” species:
#       joint species distribution models (JSDMS)
# ===============================================
# 8.4.1 (no code)

# 8.4.2 The latent-variable occupancy model of Tobler et al. (2019)
# -----------------------------------------------------------------

# Load data set from Tobler et al. (2019)
library(AHMbook)
data(SwissAtlasHa)
str(dat <- SwissAtlasHa)

# Grab counts and threshold them to become detection/nondetection data
str(counts <- dat$counts)   # Array of counts (2318 x 3 x 78)
y <- counts
y[y > 1] <- 1               # Array of detection/nondetection data
str(y)

# Subset sites to a random sample of 1200 (about half)
set.seed(1, sample.kind="Rounding")
sel.sites <- sort(sample(1:2318, 1200, replace = FALSE))
str(cc <- counts[sel.sites,,])
str(yy <- y[sel.sites,,])

# Calculate sum of max counts across sites and obs. number of occ. sites
tmp <- apply(cc, c(1,3), max, na.rm = TRUE)
summax <- apply(tmp, 2, sum)        # p-ignorant estimate of Ntotal
tmp <- apply(yy, c(1,3), max, na.rm = TRUE)
nobs <- apply(tmp, 2, sum)          # p-ignorant estimate of sum(z)
sort(nobs)                          # Look at species ordered by nobs

# Restrict data set according to number of observed occurrences
sel.species <- which(nobs > 39)     # yields about 30 species left
cc <- cc[,, sel.species]
yy <- yy[,,sel.species]
str(cc) ; str(yy) # Counts and det/nondet. for 1200 sites and approx. 30 spec.

# Determine sample sizes
(nsites <- dim(yy)[1])
(nspec <- dim(yy)[3])
table(nreps <- dat$sitecovs[sel.sites,'nsurveys']) # 2 or 3 surveys per site

# Grab, restrict, scale and mean-impute covariates
library(abind)
str(xocc <- as.matrix(dat$sitecovs[sel.sites,3:6])) # Occ. covariates
str(xocc <- cbind(xocc, xocc^2))      # Add squares of covariates
xocc <- scale(xocc)                   # Scale column-wise
str(xdet <- dat$dates[sel.sites,])    # Detection covariates
xdettmp <- standardize(xdet)          # Scale matrix-wide
xdettmp[is.na(xdettmp)] <- 0          # Mean-impute dates of 3rd survey
str(xdet <- abind(xdettmp, xdettmp^2, along = 3) )

# Provide for different numbers of LVs: could try 2, 5, 10, 15
NLV <- c(2, 5, 10, 15)                # Here we will just take 2

# Bundle data (incl. choice of number of LVs)
str(bdata <- list(y = aperm(yy, c(1,3,2)), Xocc = xocc, Xdet = xdet,
    ncov.occ = ncol(xocc), ncov.det = 2, nlv = NLV[1], nsites = nsites,
    nspec = nspec, nreps = nreps))
# List of 9
# $ y       : num [1:1200, 1:30, 1:3] 0 0 0 1 0 0 0 0 0 0 ...
# $ Xocc    : num [1:1200, 1:8] -1.033 -0.724 -0.413 -1.055 -0.726 ...
# $ Xdet    : num [1:1200, 1:3, 1:2] -1.34 -1.38 -1.42 -1.68 -1.34 ...
# $ ncov.occ: int 8
# $ ncov.det: num 2
# $ nlv     : num 2
# $ nsites  : int 1200
# $ nspec   : int 30
# $ nreps   : num [1:1200] 3 3 3 3 3 3 3 3 3 3 ...

# Specify model in BUGS language
cat(file = "JSDMocc.txt", "
model{
  # Community priors for occupancy
  mu.beta0 <- logit(mean.psi0) # Intercept
  mean.psi0 ~ dunif(0, 1)
  tau.beta0 <- pow(sd.beta0, -2)
  sd.beta0 ~ dunif(0, 2)
  for(v in 1:ncov.occ) { # Coefficients
    mu.beta[v] ~ dnorm(0, 0.2)
    tau.beta[v] <- pow(sd.beta[v], -2)
    sd.beta[v] ~ dunif(0, 2)
  }
  # Community priors for detection
  mu.alpha0 <- logit(mean.p) # Intercept
  mean.p ~ dunif(0, 1)
  tau.alpha0 <- pow(sd.alpha0, -2)
  sd.alpha0 ~ dunif(0, 1)
  for(v in 1:ncov.det) { # Coefficients
    mu.alpha[v] ~ dnorm(0, 0.1)
    tau.alpha[v] <- pow(sd.alpha[v], -2)
    sd.alpha[v] ~ dunif(0, 1)
  }

  # Define species random effects for all coefficients
  for (k in 1:nspec) {
    # Random species effects in the occupancy model
    beta0[k] ~ dnorm(mu.beta0, tau.beta0) # Intercepts
    for(v in 1:ncov.occ) { # Coefficients
      beta[k, v] ~ dnorm(mu.beta[v], tau.beta[v])
    }
    # Random effects for detection
    alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0) # Intercepts
    for(v in 1:ncov.det) { # Coefficients
      alpha[k, v] ~ dnorm(mu.alpha[v], tau.alpha[v])
    }
  }

  # Priors for latent variables: standard Normal rv
  for(i in 1:nsites) {
    for(l in 1:nlv){
      LV[i,l] ~ dnorm(0, 1)
    }
  }

  # Latent variable coefficients with constraints
  # Diagonal elements positive, upper diagonal equal to 0
  for(l in 1:(nlv-1)){
    for(l2 in (l+1):nlv){
      lv.coef[l,l2] <- 0
    }
  }
  ## Sign constraints on diagonal elements
  for(l in 1:nlv) {
    lv.coef[l,l] ~ dunif(0, 1)
  }
  # Lower diagonal free
  for(l in 2:nlv){
    for(l2 in 1:(l-1)){
      lv.coef[l,l2] ~ dunif(-1, 1)
    }
  }
  # Other elements free
  for(l in (nlv+1):nspec) {
    for(l2 in 1:nlv){
      lv.coef[l,l2] ~ dunif(-1, 1)
    }
  }

  # Define the multi-species occupancy model
  for (i in 1:nsites) {                          # Loop over sites
    for (k in 1:nspec) {                         # Loop over species
      # Probit link GLM for occupancy via auxiliary variable approach
      eta[i,k] <- beta0[k] + inprod(beta[k, ], Xocc[i, ]) +
      inprod(lv.coef[k,], LV[i, ])
      # Draw Gaussian auxiliary variable, with variance constrained to 1
      mu.psi[i,k] ~ dnorm(eta[i,k], 1/(1-sum(lv.coef[k,1:nlv]^2)))
      z[i,k] <- step(mu.psi[i,k])
      # Bernoulli GLM for detection
      for (j in 1:nreps[i]) {                    # Loop over 2 or 3 surveys
        logit(p[i,k,j]) <- alpha0[k] + alpha[k, 1] *
        Xdet[i,j,1] + alpha[k, 2] * Xdet[i,j,2]
        y[i,k,j] ~ dbern(p[i,k,j]*z[i,k])
      }
    }
  }
}
")

# 'Blind' initial values
inits <- function() { # nsites = nsites; nspec = nspec;
  lv.coef <- matrix(0, nspec, bdata$nlv)
  lv.coef[1:bdata$nlv, 1:bdata$nlv] <- 0
  for(l in 1:bdata$nlv-1){ lv.coef[l, (l+1):bdata$nlv] <- NA}
  LV <- matrix(rnorm(bdata$nlv * nsites), nsites, bdata$nlv)
  lv.coef <- matrix(runif(bdata$nlv * nspec, -sqrt(1/(bdata$nlv+1)),
  sqrt(1/(bdata$nlv+1))), nspec, bdata$nlv) * lv.coef
  mu.psi <- array(1, dim = c(nsites, nspec))         # yields psi = 0.84
  list(LV = LV, lv.coef = lv.coef , mu.psi = mu.psi)
}

# Parameters to be monitored
params <- c('mean.psi0', 'mu.beta0', 'sd.beta0', 'mu.beta', 'sd.beta',
    'mean.p', 'mu.alpha0', 'sd.alpha0', 'mu.alpha', 'sd.alpha', 'beta0',
    'beta', 'alpha0', 'alpha', 'LV', 'lv.coef') # could add 'z'

# MCMC settings
# na <- 10000 ; ni <- 250000 ; nt <- 100 ; nb <- 150000 ; nc <- 3  # 3.5 days
na <- 1000 ; ni <- 2500 ; nt <- 1 ; nb <- 1500 ; nc <- 3 # ~~~ for testing, 90 mins

# Call JAGS (ART is long !), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "JSDMocc.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
summary(out1)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1) # All estimands
traceplot(out1, 'lv.coef') # Only lv coefficients
which(out1$summary[,8] > 1.1)
# out1X <- update(out1, n.iter = 50000) # Can update additional 50k, 15 hrs
out1X <- update(out1, n.iter = 500) # Can update additional 50k  # ~~~ for testing, 15 mins

# Compute the posterior mean of the correlation matrix
R <- getLVcorrMat(lv.coef = out1$sims.list$lv.coef) # From AHMbook
colnames(R) <- rownames(R) <- names(sel.species)
R                                                   # unwieldy

# Plot the correlation matrix (Fig. 8.16)
library(corrplot)
corrplot(R, type = "lower", diag = FALSE, mar = c(1,0.5,5,1), tl.col = 'black',
    tl.pos = 'ld', tl.srt = 45, xpd = TRUE, main = 'Residual correlations')
