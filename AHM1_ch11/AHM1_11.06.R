#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 11. Hierarchical models for communities
# =========================================================================

# Approximate execution time for this code: 1.7 hrs
# Run time with the full number of iterations: 2.5 hrs

library(jagsUI)

# ~~~~~~ this section requires the data prepared in section 11.3 ~~~~~~~~~~
source("AHM1_11.03.R")
# ~~~~~~ and this code from section 11.5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
orig.ele <- MHB2014$sites$elev
(mean.ele <- mean(orig.ele, na.rm = TRUE))
(sd.ele <- sd(orig.ele, na.rm = TRUE))
ele <- (orig.ele - mean.ele) / sd.ele
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 11.6 Community models that retain species identity
# ==================================================


# 11.6.1 Simplest community occupancy model: n-fold single species
#        occupancy model with species treated as fixed effects
# ------------------------------------------------------------------------
# Collapse 3D detection/nondetection data to 2D detection frequencies
ysum <- apply(y, c(1,3), sum, na.rm = TRUE) # Collapse to detection frequency
ysum[NAsites,] <- NA                     # Have to NA out sites with NA data

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey,
    nspec = dim(ysum)[2]) )

# Specify model in BUGS language
sink("model5.txt")
cat("
model {

  # Priors
  for(k in 1:nspec){          # Loop over species
    psi[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }

  # Ecological model for latent occurrence z (process model)
  for(k in 1:nspec){          # Loop over species
    for (i in 1:M) {         # Loop over sites
      z[i,k] ~ dbern(psi[k])
    }
  }

  # Observation model for observed data y
  for(k in 1:nspec){          # Loop over species
    for (i in 1:M) {
      mup[i,k] <- z[i,k] * p[k]
      ysum[i,k] ~ dbin(mup[i,k], J[i])
    }
  }

  # Derived quantities
  for(k in 1:nspec){          # Loop over species
    Nocc.fs[k] <- sum(z[,k]) # Add up number of occupied sites among the 267
  }
  for (i in 1:M) {            # Loop over sites
    Nsite[i] <- sum(z[i,])   # Add up number of occurring species at each site
  }
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as inits for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, psi = rep(0.4, nspec), p = rep(0.4, nspec))

# Parameters monitored
params <- c("psi", "p", "Nsite", "Nocc.fs")

# MCMC settings
ni <- 2500   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3

# Call JAGS from R (ART 2.1 min)
out5 <- jags(win.data, inits, params, "model5.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(4,4))  #  ~~~ replace with 'layout' argument
traceplot(out5, layout=c(4,4))   ;   print(out5, dig = 3)

# Compare observed and estimated site species richness
op <- par(cex = 1.3)
twice <- which(MHB2014$sites$nsurvey == 2)
plot(C[twice], out5$summary[291:557,1][twice], xlab = "Observed number of species",
    ylab = "Estimated number of species", frame = FALSE, xlim = c(0, 60),
    ylim = c(0, 70), col = "red", pch = 16)
segments(C[twice], out5$summary[291:557,3][twice], C[twice],
    out5$summary[291:557,7][twice], col = "red")
points(C[-twice], out5$summary[291:557,1][-twice], col = "blue", pch = 16)
segments(C[-twice], out5$summary[291:557,3][-twice], C[-twice],
    out5$summary[291:557,7][-twice], col = "blue")
abline(0,1)
par(op)

# Observed and estimated number of occupied sites for each species
# in a table
cbind(obs.occu = obs.occ, out5$summary[558:702, c(1,3,7)])

# and in a plot
plot(obs.occ, out5$summary[558:702, 1], xlab = "Observed number of occupied sites",
    ylab = "Estimated version of quantity", ylim = c(0, 267), frame = FALSE, pch = 16)
abline(0,1)
segments(obs.occ, out5$summary[558:702,3], obs.occ, out5$summary[558:702,7],
    col = "grey", lwd = 2)


# Estimated occupancy and detection probability for each species
plot(out5$summary[1:145,1], out5$summary[146:290,1], xlab = "Occupancy estimate",
    ylab = "Detection estimate", xlim = c(0,1), ylim = c(0,1), frame = FALSE, pch = 16)
segments(out5$summary[1:145,3], out5$summary[146:290,1], out5$summary[1:145,7],
    out5$summary[146:290,1], col = "grey", lwd = 2)
segments(out5$summary[1:145,1], out5$summary[146:290,3], out5$summary[1:145,1],
    out5$summary[146:290,7], col = "grey", lwd = 2)


# 11.6.2 Community occupancy model with bivariate species-specific random effects
# -----------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey,
    nspec = dim(ysum)[2], R = matrix(c(5,0,0,1), ncol = 2), df = 3) )


# Specify model in BUGS language
sink("model6.txt")
cat("
model {

  # Priors
  for(k in 1:nspec){  # Group lpsi and lp together in array eta
    lpsi[k] <- eta[k,1]
    lp[k] <- eta[k,2]
    eta[k, 1:2] ~ dmnorm(mu.eta[], Omega[,])
  }
  # Hyperpriors
  # Priors for mu.lpsi=mu.eta[1] and mu.lp=mu.eta[2]
  # probs = community means of occupancy and detection probability
  for(v in 1:2){
    mu.eta[v] <- log(probs[v] / (1-probs[v]))
    probs[v] ~ dunif(0,1)
  }
  # Prior for variance-covariance matrix
  Omega[1:2, 1:2] ~ dwish(R[,], df)
  Sigma[1:2, 1:2] <- inverse(Omega[,])

  # Ecological model for latent occurrence z (process model)
  for(k in 1:nspec){
    logit(psi[k]) <- lpsi[k]   # Must take outside of i loop
    for (i in 1:M) {
      z[i,k] ~ dbern(psi[k])
    }
  }

  # Observation model for observed data y
  for(k in 1:nspec){
    logit(p[k]) <- lp[k]       # Needs to be outside of i loop
    for (i in 1:M) {
      mu.p[i,k] <- z[i,k] * p[k]
      ysum[i,k] ~ dbin(mu.p[i,k], J[i])
    }
  }

  # Derived quantities
  rho <- Sigma[1,2] / sqrt(Sigma[1,1] * Sigma[2,2])  # Correlation coefficient
  for(k in 1:nspec){
    Nocc.fs[k] <- sum(z[,k])   # Number of occupied sites among the 267
  }
  for (i in 1:M) {
    Nsite[i] <- sum(z[i,])     # Number of occurring species
  }
}
",fill = TRUE)
sink()


# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as starting values for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, Omega = matrix(c(1,0,0,1), ncol = 2),
    eta = matrix(0, nrow = nspec, ncol = 2))

# Parameters monitored
params <- c("mu.eta", "probs", "psi", "p", "Nsite", "Nocc.fs", "Sigma", "rho")

# MCMC settings
# ni <- 20000   ;   nt <- 15   ;   nb <- 5000   ;   nc <- 3
ni <- 2000   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3  # ~~~~ for testing

# Call JAGS from R (ART 12 min), check traceplots and summarize posteriors
out6 <- jags(win.data, inits, params, "model6.txt", n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out6, c('mu.eta', 'probs', 'Sigma', 'rho'))
print(out6, 3)

# Graphically compare some estimates between fixed- and random-effects model
op <- par(mfrow = c(2,2))      # not shown
# Species-specific occupancy (probability scale)
plot(out5$summary[1:145,1], out6$summary[5:149,1],
    main = "Species-specific occupancy probability")   ;   abline(0,1)
# Species-specific detection (probability scale)
plot(out5$summary[146:290,1], out6$summary[150:294,1],
    main = "Species-specific detection probability")   ;   abline(0,1)
# Site-specific species richness
plot(out5$summary[291:557,1], out6$summary[295:561,1],
    main = "Site-specific species richness (conditional on list of 145 detected)")
abline(0,1)
# Species-specific number of presences
plot(out5$summary[558:702,1], out6$summary[562:706,1],
    main = "Species-specific number of presences (in 267 sites)")
abline(0,1)
par(op)

# Estimated occupancy and detection probability for each species
plot(out6$summary[5:149,1], out6$summary[150:294,1],
    xlab = "Occupancy estimate", ylab = "Detection estimate",
    xlim = c(0,1), ylim = c(0,1), frame = FALSE, pch = 16)
segments(out6$summary[5:149,3], out6$summary[150:294,1], out6$summary[5:149,7],
    out6$summary[150:294,1], col = "grey", lwd = 2)
segments(out6$summary[5:149,1], out6$summary[150:294,3], out6$summary[5:149,1],
    out6$summary[150:294,7], col = "grey", lwd = 2)



# 11.6.3 Modeling species-specific effects in community occupancy models
# ------------------------------------------------------------------------
# Look at distribution of body mass among 145 observed species (see errata)
mass <- MHB2014$species$body.mass[-toss.out] # Get  species mass of observed species
hist(log10(mass), breaks = 40, col = "grey")      # Look at log10
gmass <- as.numeric(log10(mass) %/% 1.3 + 1)      # size groups 1, 2 and 3
gmass[gmass == 4] <- 3                            # Mute swan is group 3, too

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, g = gmass, M = nrow(ysum),
    J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2]) )

# Specify model in BUGS language
sink("model7.txt")
cat("
model {

  # Priors
  for(k in 1:nspec){      # loop over species
    lpsi[k] ~ dnorm(mu.lpsi[g[k]], tau.lpsi[g[k]]) # note g-dependence now
    lp[k] ~ dnorm(mu.lp[g[k]], tau.lp[g[k]])
  }

  # Hyperpriors
  for(g in 1:3){          # loop over groups (g)
    mu.lpsi[g] <- logit(mu.psi[g])      # everything is indexed g now
    mu.lp[g] <- logit(mu.p[g])
    mu.psi[g] ~ dunif(0,1)
    mu.p[g] ~ dunif(0,1)
    tau.lpsi[g] <- pow(sd.lpsi[g], -2)
    sd.lpsi[g] ~ dunif(0,5)
    tau.lp[g] <- pow(sd.lp[g], -2)
    sd.lp[g] ~ dunif(0,5)
  }

  # Ecological model for latent occurrence z (process model)
  for(k in 1:nspec){      # no change at all down here in model
    logit(psi[k]) <- lpsi[k]
    for (i in 1:M) {
      z[i,k] ~ dbern(psi[k])
    }
  }

  # Observation model for observed data ysum
  for(k in 1:nspec){      # Loop over species
    logit(p[k]) <- lp[k]
    for (i in 1:M) {
      mu.px[i,k] <- z[i,k] * p[k]  # call mu.px to avoid conflict with above
      ysum[i,k] ~ dbin(mu.px[i,k], J[i])
    }
  }

  # Derived quantities
  for(k in 1:nspec){          # Loop over species
    Nocc.fs[k] <- sum(z[,k]) # Number of occupied sites among the 267
  }
  for (i in 1:M) {            # Loop over sites
    Nsite[i] <- sum(z[i,])   # Number of occurring species at each site
  }
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst)

# Parameters monitored
params <- c("mu.psi", "mu.lpsi", "sd.lpsi", "mu.p", "mu.lp", "sd.lp")

# MCMC settings
# ni <- 6000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
ni <- 600   ;   nt <- 1   ;   nb <- 200   ;   nc <- 3  # ~~~~~ for testing

# Call JAGS from R (ART 6 min), look at convergence and summarize posteriors
out7 <- jags(win.data, inits, params, "model7.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out7)
print(out7, dig = 3)


# Bundle and summarize data set
logmass <- as.numeric(log10(mass))         # Take log10 of body mass
str( win.data <- list(ysum = ysum, logmass = logmass, M = nrow(ysum),
    nsite = nrow(ysum), J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2]) )


# Specify model in BUGS language
sink("model8.txt")
cat("
model {

  # Priors
  for(k in 1:nspec){              # loop over species
    lpsi[k] ~ dnorm(mu.lpsi[k], tau.lpsi[k]) # now all indexed by k, not g
    tau.lpsi[k] <- 1/var.lpsi[k]
    lp[k] ~ dnorm(mu.lp[k], tau.lp[k])
    tau.lp[k] <- 1/var.lp[k]
    mu.lpsi[k] <- delta0.lpsi + delta1.lpsi * logmass[k]
    mu.lp[k] <- delta0.lp + delta1.lp * logmass[k]
    log(var.lpsi[k]) <- phi0.lpsi + phi1.lpsi * logmass[k]
    log(var.lp[k]) <- phi0.lp + phi1.lp * logmass[k]
  }
  # Priors for regression params for means
  delta0.lpsi ~ dnorm(0, 0.01)
  delta1.lpsi ~ dnorm(0, 0.01)
  delta0.lp ~ dnorm(0, 0.01)
  delta1.lp ~ dnorm(0, 0.01)
  # Priors for regression params for variances
  phi0.lpsi ~ dnorm(0, 0.01)
  phi1.lpsi ~ dnorm(0, 0.01)
  phi0.lp ~ dnorm(0, 0.01)
  phi1.lp ~ dnorm(0, 0.01)

  # Ecological model for latent occurrence z (process model)
  for(k in 1:nspec){
    logit(psi[k]) <- lpsi[k]
    for (i in 1:M) {
      z[i,k] ~ dbern(psi[k])
    }
  }

  # Observation model for observed data ysum
  for(k in 1:nspec){              # Loop over species
    logit(p[k]) <- lp[k]
    for (i in 1:M) {
      mu.p[i,k] <- z[i,k] * p[k]
      ysum[i,k] ~ dbin(mu.p[i,k], J[i])
    }
  }

  # Derived quantities
  for(k in 1:nspec){          # Loop over species
    Nocc.fs[k] <- sum(z[,k]) # Number of occupied sites among the 267
  }
  for (i in 1:M) {            # Loop over sites ## see errata
    Nsite[i] <- sum(z[i,])   # Number of occurring species at each site
  }
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, delta0.lpsi = rnorm(1), delta1.lpsi = rnorm(1),
    delta0.lp = rnorm(1), delta1.lp = rnorm(1), phi0.lpsi = rnorm(1),
    phi1.lpsi = rnorm(1), phi0.lp = rnorm(1), phi1.lp = rnorm(1))

# Parameters monitored
params <- c("delta0.lpsi", "delta1.lpsi", "delta0.lp", "delta1.lp", "phi0.lpsi",
    "phi1.lpsi", "phi0.lp", "phi1.lp", "psi", "p", "Nocc.fs", "Nsite")

# MCMC settings
# ni <- 12000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
ni <- 1200   ;   nt <- 1   ;   nb <- 200   ;   nc <- 3  # ~~~~ for testing

# Call JAGS from R (ART 12 min), look at convergence and summarize posteriors
out8 <- jags(win.data, inits, params, "model8.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out8, c('delta0.lpsi', 'delta1.lpsi', 'delta0.lp', 'delta1.lp',
    'phi0.lpsi', 'phi1.lpsi', 'phi0.lp', 'phi1.lp'))
print(out8, dig = 3)

# ~~~~~ suggest saving for use later ~~~~~~~~~~~~~
save(out5, out6, out7, out8, file="AHM1_11.06_JAGSoutput.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get covariate values for prediction
predm <- seq(10, 10000,,500)        # Predict for mass of 10g to 10 kg
pred.logm <- log10(predm)

# Compute predictions (all in one array)
tmp <- out8$sims.list               # Grab simulation list
nsamp <- out8$mcmc.info$n.samples   # Number of MCMC samples
pred <- array(NA, dim = c(500, nsamp, 4)) # Array for predictions
# [,,1] mu.psi, [,,2] mu.p, [,,3] var.lpsi, [,,4] var.lp
for(i in 1:nsamp){                  # Fill array
  pred[,i,1] <- plogis(tmp$delta0.lpsi[i] + tmp$delta1.lpsi[i] * pred.logm)
  pred[,i,2] <- plogis(tmp$delta0.lp[i] + tmp$delta1.lp[i] * pred.logm)
  pred[,i,3] <- exp(tmp$phi0.lpsi[i] + tmp$phi1.lpsi[i] * pred.logm)
  pred[,i,4] <- exp(tmp$phi0.lp[i] + tmp$phi1.lp[i] * pred.logm)
}

# Plot posterior mean and a random sample of 100 from posterior of regression
selection <- sample(1:nsamp, 100)   # Choose random sample of MCMC output
op <- par(mfrow = c(2,2), mar = c(5,5,2,2))
matplot(predm, pred[,selection,1], ylab = "Occupancy mean",
    xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey",
    ylim = c(0, 0.4), frame = FALSE)
lines(predm, apply(pred[,,1], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,2], ylab = "Detection mean",
    xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey",
    ylim = c(0, 0.8), frame = FALSE)
lines(predm, apply(pred[,,2], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,3], ylab = "Occupancy variance",
    xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey",
    ylim = c(0, 8), frame = FALSE)
lines(predm, apply(pred[,,3], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,4], ylab = "Detection variance",
    xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey",
    ylim = c(0, 8), frame = FALSE)
lines(predm, apply(pred[,,4], 1, mean), lwd = 3, col = "blue")
par(op)


# 11.6.4 Modeling species richness in a two-step analysis
# ------------------------------------------------------------------------
# Extract estimates of N from model 5
N.pm <- out5$summary[291:557, 1]       # Posterior means of Nsite
N.psd <- out5$summary[291:557, 2]      # ... posterior sd's of Nsite
N.cri <- out5$summary[291:557, c(3,7)] # ... CRL's of Nsite

# Plot estimates as a function of elevation
elev <- MHB2014$sites$elev
plot(elev, N.pm, xlab = "Altitude (m a.s.l.)",
    ylab = "Estimated avian species richness", ylim = c(0, 70), frame = FALSE)
segments(elev, N.cri[,1], elev, N.cri[,2], col = "grey")
lines(smooth.spline(N.pm ~ elev, w = 1 / N.psd), col = "grey", lwd = 3)

# Bundle and summarize data set
pred.ele <- (seq(200, 2750,5) - mean.ele) / sd.ele # elevation standardised
str(win.data <- list(ele = ele, N = N.pm, psd = N.psd, n = length(N.pm),
    pred.ele = pred.ele, npred = length(pred.ele)))

# Define model in BUGS language
sink("meta.analysis.txt")
cat("
model{

  # Priors
  for(v in 1:4){         # Priors for intercept and polynomial coefficients
    beta[v] ~ dnorm(0, 0.0001)
  }
  tau.site <- pow(sd.site, -2)
  sd.site ~ dunif(0,10)

  # Likelihood
  for (i in 1:n){
    N[i] ~ dnorm(muN[i], tau.psd[i]) # Measurement error model for estimated N
    tau.psd[i] <- pow(psd[i], -2)    # 'Known' part of residual: meas. error
    muN[i] <- beta[1] + beta[2] * ele[i] + beta[3] * pow(ele[i],2) +
    beta[4] * pow(ele[i],3) + eps.site[i] # add another source of uncertainty
    eps.site[i] ~ dnorm(0, tau.site) # this is the usual 'residual'
  }
  # Get predictions for plot
  for(i in 1:npred){
    Npred[i] <- beta[1] + beta[2] * pred.ele[i] + beta[3] * pow(pred.ele[i],2) + beta[4] * pow(pred.ele[i],3)
  }
} # end model
",fill=TRUE)
sink()

# Initial values, params monitored, and MCMC settings
inits <- function() list(beta = rnorm(4))
params <- c("beta", "sd.site", "Npred")
# ni <- 12000   ;   nt <- 10   ;   nb <- 2000   ;   nc <- 3
ni <- 1200   ;   nt <- 1   ;   nb <- 200   ;   nc <- 3  # ~~~~~ for testing

# Call JAGS and summarize posterior
out <- jags(win.data, inits, params, "meta.analysis.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # parallel=FALSE)
  parallel=TRUE)  # ~~~~~ for testing
print(out, 3)

lines(seq(200, 2750,5), out$mean$Npred, col = "blue", lwd = 3)
matlines(seq(200,2750,5), out$summary[6:516,c(3, 7)], col = "blue",
    lwd = 2, lty= "dashed")


op <- par(mfrow = c(3, 3), mar = c(5,4,3,2))
oldask <- devAskNewPage(ask=dev.interactive(orNone=TRUE)) #~~~~ to replace browser call
for(i in 1:267){
   for(j in 1:267){
      plot(jitter(out5$sims.list$Nsite[,i]), jitter(out5$sims.list$Nsite[,j]),
      main = paste("Joint posterior sites", i, "and", j))
#   browser() #~~~~~ incompatible with automated testing
   }
}
devAskNewPage(oldask) #~~~~ clean up
par(op)
