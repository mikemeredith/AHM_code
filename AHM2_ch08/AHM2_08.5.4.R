#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8 : MODELING INTERACTIONS AMONG SPECIES
# ===============================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 20 mins
# Run time with the full number of iterations: 28 hrs

library(AHMbook)
library(abind)
library(jagsUI)
library(corrplot)

# ~~~ Need to prepare data as in 8.4.2 ~~~~~~~~~~~~~
data(SwissAtlasHa)
dat <- SwissAtlasHa
counts <- dat$counts # Array of counts (2318 x 3 x 78)
y <- counts
y[y > 1] <- 1 # Array of detection/nondetection data
set.seed(1, sample.kind="Rounding")
sel.sites <- sort(sample(1:2318, 1200, replace = FALSE))
cc <- counts[sel.sites,,]
yy <- y[sel.sites,,]
tmp <- apply(yy, c(1,3), max, na.rm = TRUE)
nobs <- apply(tmp, 2, sum) # p-ignorant estimate of sum(z)
sel.species <- which(nobs > 39) # yields about 30 species left
cc <- cc[,, sel.species]
yy <- yy[,,sel.species]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 8.5 Joint models for abundance
# ==============================

# 8.5.4 The n-mixture model variant of the latent-variable model of
#       Tobler et al. (2019)
# -----------------------------------------------------------------

# Determine sample sizes
(nsites <- dim(cc)[1]) ; (nspec <- dim(cc)[3]) ; (nreps <- dim(cc)[2])
table(nreps <- dat$sitecovs[sel.sites,'nsurveys']) # 2/3 surveys per site

# Prepare occupancy and detection covariates and choose latent vars
library(abind)
str(xocc <- as.matrix(dat$sitecovs[sel.sites,3:6]))# Occ. covariates
str(xocc <- cbind(xocc, xocc^2))                   # Add squares of covariates
xocc <- scale(xocc)                                # Scale column-wise
str(xdet <- dat$dates[sel.sites,])                 # Detection covariates
xdettmp <- standardize(xdet)                       # Scale matrix-wide
xdettmp[is.na(xdettmp)] <- 0                       # Mean-impute dates of 3rd survey
str(xdet <- abind(xdettmp, xdettmp^2, along = 3) )

# Bundle data
nlv <- 2 # Choose number of latent variables
str(bdata <- list(C = aperm(cc, c(1,3,2)), Xocc = xocc, Xdet = xdet,
    ncov.occ = ncol(xocc), ncov.det = 2, nlv = nlv, nsites = nsites,
    nspec = nspec, nreps = nreps))
# List of 9
# $ C       : num [1:1200, 1:30, 1:3] 0 0 0 1 0 0 0 0 0 0 ...
# $ Xocc    : num [1:1200, 1:8] -1.033 -0.724 -0.413 -1.055 -0.726 ...
# $ Xdet    : num [1:1200, 1:3, 1:2] -1.34 -1.38 -1.42 -1.68 -1.34 ...
# $ ncov.occ: int 8
# $ ncov.det: num 2
# $ nlv     : num 2
# $ nsites  : int 1200
# $ nspec   : int 30
# $ nreps   : num [1:1200] 3 3 3 3 3 3 3 3 3 3 ...

# Specify model in BUGS language
cat(file = "JSDMnmix.txt", "
model{
  # Community priors for abundance
  mu.beta0 <- log(mean.lambda) # Intercepts
  mean.lambda ~ dunif(0, 1)
  tau.beta0 <- pow(sd.beta0, -2)
  sd.beta0 ~ dunif(0, 3)
  for(v in 1:ncov.occ) { # Coefficients
    mu.beta[v] ~ dnorm(0, 0.1)
    tau.beta[v] <- pow(sd.beta[v], -2)
    sd.beta[v] ~ dunif(0, 3)
  }

  # Community priors for detection
  mu.alpha0 <- logit(mean.p) # Intercept
  mean.p ~ dunif(0, 1)
  tau.alpha0 <- pow(sd.alpha0, -2)
  sd.alpha0 ~ dunif(0, 2)
  for(v in 1:ncov.det) { # Coefficients
    mu.alpha[v] ~ dnorm(0, 1)
    tau.alpha[v] <- pow(sd.alpha[v], -2)
    sd.alpha[v] ~ dunif(0, 1)
  }

  # Define species random effects for all coefficients
  for (k in 1:nspec) {
    # Random species effects in the abundance model
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

  # Priors for latent variables: standard Normal rv #
  for(i in 1:nsites) {
    for(l in 1:nlv){
      LV[i,l] ~ dnorm(0, 1)
    }
  }

  # Latent variable coefficients with constraints
  # (slightly different from occupancy version of the model)
  # Upper diagonal equal to 0
  for(l in 1:(nlv-1)){
    for(l2 in (l+1):nlv){
      lv.coef[l,l2] <- 0
    }
  }

  # Sign constraints on diagonal elements (positive)
  for(l in 1:nlv) {
    lv.coef[l,l] ~ dnorm(0, 0.1)I(0, )
  }

  # Lower diagonal free
  for(l in 2:nlv){
    for(l2 in 1:(l-1)){
      lv.coef[l,l2] ~ dnorm(0, 0.1)
    }
  }

  # Other elements also free
  for(l in (nlv+1):nspec) {
    for(l2 in 1:nlv){
      lv.coef[l,l2] ~ dnorm(0, 0.1)
    }
  }

  # Define the multi-species binomial N-mixture model
  for (i in 1:nsites) { # Loop over sites
    for (k in 1:nspec) { # Loop over species
      # State model
      N[i,k] ~ dpois(exp(loglam[i,k]))
      loglam[i,k] <- beta0[k] + inprod(beta[k, ], Xocc[i, ]) +
      inprod(lv.coef[k,], LV[i, ])
      # Observation model
      for (j in 1:nreps[i]) { # Loop over 2 or 3 surveys
        logit(p[i,k,j]) <- alpha0[k] + alpha[k, 1] * Xdet[i,j,1] +
        alpha[k, 2] * Xdet[i,j,2]
        C[i,k,j] ~ dbinom(p[i,k,j], N[i,k])
      }
    }
  }
}
")

# (Original) Initial values
inits <- function() { # nsites = nsites; nspec = nspec;
  tmp <- apply(bdata$C, c(1,2), max)
  tmp[is.na(tmp)] <- 0
  Nst <- tmp + 5
  lv.coef <- matrix(0, nspec, bdata$nlv)
  lv.coef[1:bdata$nlv, 1:bdata$nlv] <- 0
  for(l in 1:bdata$nlv-1){ lv.coef[l, (l+1):bdata$nlv] <- NA}
  LV <- matrix(rnorm(bdata$nlv * nsites), nsites, bdata$nlv)
  lv.coef <- matrix(runif(bdata$nlv * nspec, -sqrt(1/(bdata$nlv+1)),
  sqrt(1/(bdata$nlv+1))), nspec, bdata$nlv) * lv.coef
  beta0 <- rep(0, nspec)                          # yields lambda = 1
  list(N = Nst, LV = LV, lv.coef = lv.coef , beta0 = beta0)
}

# Parameters to be monitored
params <- c('mean.lambda', 'mu.beta0', 'sd.beta0', 'mu.beta', 'sd.beta',
    'mean.p', 'mu.alpha0', 'sd.alpha0', 'mu.alpha', 'sd.alpha', 'beta0',
    'beta', 'alpha0', 'alpha', 'LV', 'lv.coef') # could add 'N'

# MCMC settings
# na <- 10000 ; ni <- 150000 ; nt <- 50 ; nb <- 100000 ; nc <- 2
na <- 100 ; ni <- 1500 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ for testing, 21 mins

# Call JAGS (ART 36 h), check convergence and summarize posteriors
out2X <- jags(bdata, inits, params, "JSDMnmix.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(4,4))  # ~~~ replaced with 'layout' argument
traceplot(out2X, layout=c(4,4)) # All params
traceplot(out2X,'lv.coef', layout=c(4,4)) # Only LV coefficients
traceplot(out2X, 'deviance') # Only the deviance
which(out2X$summary[,8] > 1.1)


# Compute the posterior mean of the correlation matrix
R <- getLVcorrMat(lv.coef = out2X$sims.list$lv.coef, type= "Nmix")
colnames(R) <- rownames(R) <- names(sel.species)
R                           # unwieldy

# Plot the correlation matrix (Fig. 8.18) and dendrogram (Fig. 8.19)
library(corrplot)
corrplot(R, type = "lower", diag = FALSE, mar = c(1,0.5,5,1),
    tl.col = 'black', tl.pos = 'ld', tl.srt = 45, xpd = TRUE,
    main = 'Residual correlations')
dist <- as.dist(1 - R) # so R = 1 -> dist = 0
plot(hclust(dist), xlab = "", sub = "", ylab = "Correlation coefficient",
    yaxt = 'n') # Fig. 8.19
axis(2, at = c(0, 0.5, 1, 1.5, 2), labels = c(1, 0.5, 0, -0.5, -1), las = 1)
