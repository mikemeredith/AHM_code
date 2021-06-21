#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 11. Hierarchical models for communities
# =========================================================================

# Approximate execution time for this code: 50 mins
# Run time with the full number of iterations: 7.7 hrs

library(jagsUI)

# ~~~~~~ this section requires the data prepared in section 11.3 ~~~~~~~~~~
source("AHM1_11.03.R")
# ~~~~~~ and this code from section 11.6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collapse 3D detection/nondetection data to 2D detection frequencies
ysum <- apply(y, c(1,3), sum, na.rm = TRUE) # Collapse to detection frequency
ysum[NAsites,] <- NA                     # Have to NA out sites with NA data
# ~~~~~~ and this code from section 11.5 ~~~~~~~~~~
# Quadrat elevation and forest cover
orig.ele <- MHB2014$sites$elev
(mean.ele <- mean(orig.ele, na.rm = TRUE))
(sd.ele <- sd(orig.ele, na.rm = TRUE))
ele <- (orig.ele - mean.ele) / sd.ele
orig.forest <- MHB2014$sites$forest
(mean.forest <- mean(orig.forest, na.rm = TRUE))
(sd.forest <- sd(orig.forest, na.rm = TRUE))
forest <- (orig.forest - mean.forest) / sd.forest

# Get survey date and survey duration and standardise both
# Survey date (this is Julian date, with day 1 being April 1)
orig.DAT <- MHB2014$date
(mean.date <- mean(orig.DAT, na.rm = TRUE))
(sd.date <- sd(c(orig.DAT), na.rm = TRUE))
DAT <- (orig.DAT - mean.date) / sd.date      # scale
DAT[is.na(DAT)] <- 0                         # impute missings
# Survey duration (in minutes)
orig.DUR <- MHB2014$dur
(mean.dur <- mean(orig.DUR, na.rm = TRUE))
(sd.dur <- sd(c(orig.DUR), na.rm = TRUE))
DUR <- (orig.DUR - mean.dur) / sd.dur        # scale
DUR[is.na(DUR)] <- 0                         # mean impute missings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 11.7 The Dorazio-Royle (DR) community occupancy model with data augmentation
# ============================================================================


# 11.7.1 The simplest DR community model with data augmentation
# ------------------------------------------------------------------------
# Augment data set (DA part)
nz <- 150                # Number of potential species in superpopulation
M <- nspec + nz          # Size of augmented data set ('superpopulation')
yaug <- cbind(ysum, array(0, dim=c(nsite, nz))) # Add all zero histories

# Bundle and summarize data set
str( win.data <- list(yaug = yaug, nsite = nrow(ysum),
    nrep = MHB2014$sites$nsurvey, M = M, nspec = nspec, nz = nz) )

# Specify model in BUGS language
sink("model9.txt")
cat("
model {

  # Priors to describe heterogeneity among species in community
  for(k in 1:M){                  # Loop over all species in augmented list
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm(mu.lp, tau.lp)
  }

  # Hyperpriors to describe full community
  omega ~ dunif(0,1)              # Data augmentation or 'occupancy' parameter
  mu.lpsi ~ dnorm(0,0.001)        # Community mean of occupancy (logit)
  mu.lp ~ dnorm(0,0.001)          # Community mean of detection (logit)
  tau.lpsi <- pow(sd.lpsi, -2)
  sd.lpsi ~ dunif(0,5)            # Species heterogeneity in logit(psi)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0,5)              # Species heterogeneity in logit(p)

  # Superpopulation process:this is the 'paramater expansion' part of PX-DA
  for(k in 1:M){
    w[k] ~ dbern(omega)           # Metacommunity membership indicator
  }                               # (or data augmentation variable)

  # Ecological model for latent occurrence z (process model)
  for(k in 1:M){
    mu.psi[k] <- w[k] * psi[k]    # species not part of community zeroed out for z
    logit(psi[k]) <- lpsi[k]
    for (i in 1:nsite) {
      z[i,k] ~ dbern(mu.psi[k])
    }
  }

  # Observation model for observed detection frequencies
  for(k in 1:M){
    logit(p[k]) <- lp[k]
    for (i in 1:nsite) {
      mu.p[i,k] <- z[i,k] * p[k]  # non-occurring species are zeroed out for p
      yaug[i,k] ~ dbin(mu.p[i,k], nrep[i])
    }
  }

  # Derived quantities
  for(k in 1:M){
     Nocc.fs[k] <- sum(z[,k])     # Number of occupied sites among the 267
  }
  for (i in 1:nsite) {
    Nsite[i] <- sum(z[i,])        # Number of occurring species at each site
  }
  n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species in metacommunity
  Ntotal <- sum(w[])              # Total metacommunity size (= nspec + n0)
}
",fill = TRUE)
sink()

# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at 'occurring'
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto for z
inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz),
    lp = rnorm(n = nspec+nz))

# Parameters monitored
params <- c("mu.lpsi", "sd.lpsi", "mu.lp", "sd.lp", "psi", "p", "Nsite",
    "Ntotal", "omega", "n0")

# MCMC settings
# ni <- 22000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
ni <- 2200   ;   nt <- 2   ;   nb <- 200   ;   nc <- 3  # ~~~ for testing

# Call JAGS from R (ART 62 min), check convergence and summarize posteriors
out9 <- jags(win.data, inits, params, "model9.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out9, c('mu.lpsi', 'sd.lpsi', 'mu.lp', 'sd.lp'), layout=c(2,2))
print(out9, dig = 3)
# ~~~~~ suggest saving for use later ~~~~~~~~~~~~~
save(out9, file="AHM1_11.07_out9.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot posterior distribution of site-specific species richness (Nsite)
op <- par(mfrow = c(3,3), mar = c(5,4,3,2))
oldask <- devAskNewPage(ask=dev.interactive(orNone=TRUE)) #~~~~ to replace browser call
for(i in 1:267){
   plot(table(out9$sims.list$Nsite[,i]), main = paste("Quadrat", i),
   xlab = "Local species richness", ylab = "", frame = FALSE,
   xlim = c((min(C[i], out9$sims.list$Nsite[,i], na.rm = TRUE)-2),
   max(out9$sims.list$Nsite[,i]) ))
   abline(v = C[i], col = "grey", lwd = 4)
   # browser()  # ~~~~~ incompatible with automated testing
}
devAskNewPage(ask=oldask) #~~~~ clean up
par(op)

# Plot it only for a selection of sites
op <- par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
   plot(table(out9$sims.list$Nsite[,i]), main = paste("Quadrat", i),
       xlab = "Local species richness", ylab = "", frame = FALSE,
       xlim = c((min(C[i], out9$sims.list$Nsite[,i], na.rm = TRUE)-2),
       max(out9$sims.list$Nsite[,i]) ))
   abline(v = C[i], col = "grey", lwd = 4)
}
par(op)

# Plot posterior distribution of total species richness (Ntotal)
plot(table(out9$sims.list$Ntotal), main = "", ylab = "",
    xlab = "Avian metacommunity size in Swiss MHB survey (267 1km2 quadrats)",
    frame = FALSE, xlim = c(144, 245))
abline(v = nspec, col = "grey", lwd = 4)



# 11.7.2 Dorazio-Royle community model with covariates
# ------------------------------------------------------------------------
# Augment data set: choose one of two different priors on Ntotal
nz <- 250                 # Use for vague prior on Ntotal: M = 395
nz <- 215 - nspec         # Use for informative prior on Ntotal: M = 215
yaug <- array(0, dim=c(nsite, nrep, nspec+nz)) # array with only zeroes
yaug[,,1:nspec] <- y      # copy into it the observed data

# Create same NA pattern in augmented species as in the observed species
missings <- is.na(yaug[,,1]) # e.g., third survey in high-elevation quads
for(k in (nspec+1):(nspec+nz)){
   yaug[,,k][missings] <- NA
}

# Bundle and summarize data
str(win.data <- list(y = yaug, nsite = dim(y)[1], nrep = dim(y)[2],
    nspec = dim(y)[3], nz = nz, M = nspec + nz, ele = ele, forest = forest,
    DAT = DAT, DUR = DUR) )

# Specify model in BUGS language
sink("model10.txt")
cat("
model {

  # Priors
  omega ~ dunif(0,1)
  # Priors for species-specific effects in occupancy and detection
  for(k in 1:M){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)    # Hyperparams describe community
    betalpsi1[k] ~ dnorm(mu.betalpsi1, tau.betalpsi1)
    betalpsi2[k] ~ dnorm(mu.betalpsi2, tau.betalpsi2)
    betalpsi3[k] ~ dnorm(mu.betalpsi3, tau.betalpsi3)
    lp[k] ~ dnorm(mu.lp, tau.lp)
    betalp1[k] ~ dnorm(mu.betalp1, tau.betalp1)
    betalp2[k] ~ dnorm(mu.betalp2, tau.betalp2)
    betalp3[k] ~ dnorm(mu.betalp3, tau.betalp3)
  }

  # Hyperpriors
  # For the model of occupancy
  mu.lpsi ~ dnorm(0,0.01)
  tau.lpsi <- pow(sd.lpsi, -2)
  sd.lpsi ~ dunif(0,8)   # as always, bounds of uniform chosen by trial and error
  mu.betalpsi1 ~ dnorm(0,0.1)
  tau.betalpsi1 <- pow(sd.betalpsi1, -2)
  sd.betalpsi1 ~ dunif(0, 4)
  mu.betalpsi2 ~ dnorm(0,0.1)
  tau.betalpsi2 <- pow(sd.betalpsi2, -2)
  sd.betalpsi2 ~ dunif(0,2)
  mu.betalpsi3 ~ dnorm(0,0.1)
  tau.betalpsi3 <- pow(sd.betalpsi3, -2)
  sd.betalpsi3 ~ dunif(0,2)

  # For the model of detection
  mu.lp ~ dnorm(0,0.1)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 2)
  mu.betalp1 ~ dnorm(0,0.1)
  tau.betalp1 <- pow(sd.betalp1, -2)
  sd.betalp1 ~ dunif(0,1)
  mu.betalp2 ~ dnorm(0,0.1)
  tau.betalp2 <- pow(sd.betalp2, -2)
  sd.betalp2 ~ dunif(0,1)
  mu.betalp3 ~ dnorm(0,0.1)
  tau.betalp3 <- pow(sd.betalp3, -2)
  sd.betalp3 ~ dunif(0,1)

  # Superpopulation process: Ntotal species sampled out of M available
  for(k in 1:M){
    w[k] ~ dbern(omega)
  }

  # Ecological model for true occurrence (process model)
  for(k in 1:M){
    for (i in 1:nsite) {
      logit(psi[i,k]) <- lpsi[k] + betalpsi1[k] * ele[i] +
        betalpsi2[k] * pow(ele[i],2) + betalpsi3[k] * forest[i]
      mu.psi[i,k] <- w[k] * psi[i,k]
      z[i,k] ~ dbern(mu.psi[i,k])
    }
  }

  # Observation model for replicated detection/nondetection observations
  for(k in 1:M){
    for (i in 1:nsite){
      for(j in 1:nrep){
        logit(p[i,j,k]) <- lp[k] + betalp1[k] * DAT[i,j] +
            betalp2[k] * pow(DAT[i,j],2) + betalp3[k] * DUR[i,j]
        mu.p[i,j,k] <- z[i,k] * p[i,j,k]
        y[i,j,k] ~ dbern(mu.p[i,j,k])
      }
    }
  }

  # Derived quantities
  #for(k in 1:M){
  #   Nocc.fs[k] <- sum(z[,k])       # Number of occupied sites among the 267
  #}
  for (i in 1:nsite){
     Nsite[i] <- sum(z[i,])          # Number of occurring species at each site
  }
  n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species
  Ntotal <- sum(w[])                 # Total metacommunity size

  # Vectors to save (S for 'save'; discard posterior samples for
  # all minus 1 of the potential species to save disk space)
  # we do this for nz = 250 (i.e., M = 395)
  lpsiS[1:(nspec+1)] <- lpsi[1:(nspec+1)]
  betalpsi1S[1:(nspec+1)] <- betalpsi1[1:(nspec+1)]
  betalpsi2S[1:(nspec+1)] <- betalpsi2[1:(nspec+1)]
  betalpsi3S[1:(nspec+1)] <- betalpsi3[1:(nspec+1)]
  lpS[1:(nspec+1)] <- lp[1:(nspec+1)]
  betalp1S[1:(nspec+1)] <- betalp1[1:(nspec+1)]
  betalp2S[1:(nspec+1)] <- betalp2[1:(nspec+1)]
  betalp3S[1:(nspec+1)] <- betalp3[1:(nspec+1)]
}
",fill = TRUE)
sink()


# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at occurring
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto
inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz),
    betalpsi1 = rnorm(n = nspec+nz), betalpsi2 = rnorm(n = nspec+nz),
    betalpsi3 = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz),
    betalp1 = rnorm(n = nspec+nz), betalp2 = rnorm(n = nspec+nz),
    betalp3 = rnorm(n = nspec+nz))

# Set 1
params1 <- c("omega", "mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1",
    "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3", "mu.lp",
    "sd.lp", "mu.betalp1", "sd.betalp1", "mu.betalp2", "sd.betalp2",
    "mu.betalp3", "sd.betalp3", "Ntotal", "Nsite")

# MCMC settings
# ni <- 15000   ;   nt <- 10   ;   nb <- 5000   ;   nc <- 3
ni <- 1500   ;   nt <- 1   ;   nb <- 500   ;   nc <- 3  # ~~~~~ use for testing

# Run JAGS, check convergence and summarize posteriors
out101 <- jags(win.data, inits, params1, "model10.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 2))  #  ~~~ replace with 'layout' argument
traceplot(out101, c("omega", "mu.lpsi", "sd.lpsi", "mu.betalpsi1",
    "sd.betalpsi1", "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3",
    "sd.betalpsi3", "mu.lp", "sd.lp", "mu.betalp1", "sd.betalp1",
    "mu.betalp2", "sd.betalp2", "mu.betalp3", "sd.betalp3", "Ntotal"), layout=c(2,2) )

# ~~~~~ save for use in sections 11.8 and 11.9 ~~~~~~~~~~~~~
out10 <- out101
save(out10, file="AHM1_11.07_out10.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set 2
params2 <- c("mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1",
    "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3",
    "lpsi", "betalpsi1", "betalpsi2", "betalpsi3", "lp", "betalp1",
    "betalp2", "betalp3", "z", "w")
# ni <- 12000   ;   nt <- 20   ;   nb <- 2000   ;   nc <- 3  # 3 hrs
ni <- 1200   ;   nt <- 2   ;   nb <- 200   ;   nc <- 3  # ~~~~~ use for testing
out102 <- jags.basic(win.data, inits, params2, "model10.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# library(coda) ~~~ not necessary, causes conflict with traceplot function in jagsUI
all10 <- as.matrix(out102) # Put output from 3 chains into a matrix
# ~~~~~ save for use in sections 11.8 and 11.9 ~~~~~~~~~~~~~
save(all10, file="AHM1_11.07_all10.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~ won't work with 59,349 parameters ~~~~~~~~~~~~~~~~~~~~~~~~
# summary(out102)            # Pointless, can't display 59,000*2 rows
# gelman.diag(out102)        # Error: Cannot allocate vector of size 78.7 Gb
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Comparison of main hyperparameters when M = 215 and with M = 395
# (not all code to produce this output is shown)
# ~~~~~ does not work without the rest of the code ~~~~~~~~~~~~~~~~~~~~~~~~~
# print(cbind(out10.215$summary[1:17,c(1:3, 7)], out10.395$summary[1:17, c(1:3, 7)]), 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out10 <- out101


op <- par(mfrow = c(1,2))       # Fig. 11-16
psi.sample <- plogis(rnorm(10^6, mean = out10$mean$mu.lpsi, sd = out10$mean$sd.lpsi))
p.sample <- plogis(rnorm(10^6, mean = out10$mean$mu.lp, sd = out10$mean$sd.lp))
hist(psi.sample, freq = FALSE, breaks = 50, col = "grey",
    xlab = "Species occupancy probability", ylab = "Density", main = "")
hist(p.sample, freq = FALSE, breaks = 50, col = "grey",
    xlab = "Species detection probability", ylab = "Density", main = "")
par(op)
summary(psi.sample)   ;   summary(p.sample)

op <- par(mfrow = c(2,4))  # Among-species variability in parameters (not shown)
hist(out10$sims.list$sd.lpsi, breaks = 100, col = "grey", xlim = c(0,6),
    main = "Occupancy: intercept")
abline(v = mean(out10$sims.list$sd.lpsi), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalpsi1, breaks = 100, col = "grey", xlim = c(0,3),
    main = "Occupancy: linear effect of elevation")
abline(v = mean(out10$sims.list$sd.betalpsi1), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalpsi2, breaks = 100, col = "grey", xlim = c(0,3),
    main = "Occupancy: quadratic effect of elevation")
abline(v = mean(out10$sims.list$sd.betalpsi2), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalpsi3, breaks = 100, col = "grey", xlim = c(0,3),
    main = "Occupancy: linear effect of forest cover")
abline(v = mean(out10$sims.list$sd.betalpsi3), col = "blue", lwd = 3)
hist(out10$sims.list$sd.lp, breaks = 100, col = "grey", xlim = c(0,2),
    main = "Detection: intercept")
abline(v = mean(out10$sims.list$sd.lp), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalp1, breaks = 100, col = "grey", xlim = c(0,1),
    main = "Detection: linear effect of survey date")
abline(v = mean(out10$sims.list$sd.betalp1), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalp2, breaks = 100, col = "grey", xlim = c(0,1),
    main = "Detection: quadratic linear effect of survey date")
abline(v = mean(out10$sims.list$sd.betalp2), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalp3, breaks = 100, col = "grey", xlim = c(0,1),
    main = "Detection: linear effect of survey duration")
abline(v = mean(out10$sims.list$sd.betalp3), col = "blue", lwd = 3)
par(op)

# Visualize covariate mean relationships for the average species
o.ele <- seq(200, 2500,,500)               # Get covariate values for prediction
o.for <- seq(0, 100,,500)
o.dat <- seq(15, 120,,500)
o.dur <- seq(100, 420,,500)
ele.pred <- (o.ele - mean.ele) / sd.ele
for.pred <- (o.for - mean.forest) / sd.forest
dat.pred <- (o.dat - mean.date) / sd.date
dur.pred <- (o.dur - mean.dur) / sd.dur

# Predict occupancy for elevation and forest and detection for date and duration
# Put all fourpredictions into a single
str( tmp <- out10$sims.list )              # grab MCMC samples
nsamp <- length(tmp[[1]])    # number of mcmc samples
predC <- array(NA, dim = c(500, nsamp, 4)) # "C" for 'community mean'
for(i in 1:nsamp){
   predC[,i,1] <- plogis(tmp$mu.lpsi[i] + tmp$mu.betalpsi1[i] * ele.pred +
     tmp$mu.betalpsi2[i] * ele.pred^2 )
   predC[,i,2] <- plogis(tmp$mu.lpsi[i] + tmp$mu.betalpsi3[i] * for.pred)
   predC[,i,3] <- plogis(tmp$mu.lp[i] + tmp$mu.betalp1[i] * dat.pred +
     tmp$mu.betalp2[i] * dat.pred^2 )
   predC[,i,4] <- plogis(tmp$mu.lp[i] + tmp$mu.betalp3[i] * dur.pred)
}

# Get posterior means and 95% CRIs and plot (Fig. 11–17)
pmC <- apply(predC, c(1,3), mean)
criC <- apply(predC, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))

op <- par(mfrow = c(2, 2))
plot(o.ele, pmC[,1], col = "blue", lwd = 3, type = 'l', lty = 1, frame = FALSE,
    ylim = c(0, 0.05), xlab = "Elevation (m a.s.l)",
    ylab = "Community mean occupancy")
matlines(o.ele, t(criC[,,1]), col = "grey", lty = 1)
plot(o.for, pmC[,2], col = "blue", lwd = 3, type = 'l', lty = 1, frame = FALSE,
    ylim = c(0, 0.05), xlab = "Forest cover",
    ylab = "Community mean occupancy")
matlines(o.for, t(criC[,,2]), col = "grey", lty = 1)
plot(o.dat, pmC[,3], col = "blue", lwd = 3, type = 'l', lty = 1, frame = FALSE,
    ylim = c(0.2, 0.8), xlab = "Survey date",
    ylab = "Community mean detection")
matlines(o.dat, t(criC[,,3]), col = "grey", lty = 1)
plot(o.dur, pmC[,4], col = "blue", lwd = 3, type = 'l', lty = 1, frame = FALSE,
    ylim = c(0.2, 0.8), xlab = "Survey duration",
    ylab = "Community mean detection")
matlines(o.dur, t(criC[,,4]), col = "grey", lty = 1)
par(op)

# Plot posterior distribution of site-specific species richness (Nsite)
op <- par(mfrow = c(3,3), mar = c(5,4,3,2))
oldask <- devAskNewPage(ask=dev.interactive(orNone=TRUE)) #~~~~ to replace browser call
for(i in 1:267){
   plot(table(out10$sims.list$Nsite[,i]), main = paste("Quadrat", i),
   xlab = "Local species richness", ylab = "", frame = FALSE,
   xlim = c((min(C[i], out10$sims.list$Nsite[,i], na.rm = TRUE)-2),
   max(out10$sims.list$Nsite[,i]) ))
   abline(v = C[i], col = "grey", lwd = 4)
   # browser() #~~~~ incompatible with automated testing
}
devAskNewPage(ask=oldask) #~~~~ clean up
par(op)

# Plot it only for a selection of sites (Fig. 11-18)
op <- par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
   plot(table(out10$sims.list$Nsite[,i]), main = paste("Quadrat", i),
       xlab = "Local species richness", ylab = "", frame = FALSE,
       xlim = c((min(C[i], out10$sims.list$Nsite[,i], na.rm = TRUE)-2),
       max(out10$sims.list$Nsite[,i]) ))
   abline(v = C[i], col = "grey", lwd = 4)
}
par(op)

# Plot Nsite estimates under models 9 & 10 vs. elevation (Fig. 11-19)
elev <- orig.ele
offset <- 30    # Set off elevation for better visibility
plot(elev, out9$mean$Nsite, xlab = "Elevation (metres)",
    ylab = "Community size estimate (Nsite)", frame = FALSE, ylim = c(0,60),
    pch = 16) # black: model 9
lines(smooth.spline(out9$mean$Nsite ~ elev), lwd = 3)
points(elev+offset, out10$mean$Nsite, pch = 16, col = "blue") # red: model 10
lines(smooth.spline(out10$mean$Nsite ~ elev), lwd = 3, col = "blue")


str(all10)                    # look at the MCMC output
# ~~~ note that with 'jags.basic' node names are in alphabetical order, not in
# ~~~ the order given in 'params2'
pm <- apply(all10, 2, mean)    # Get posterior means and 95% CRIs
cri <- apply(all10, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs

# Effects of date (linear and quadratic) and of duration on detection - figure 11.20
op <- par(mfrow = c(1,3), cex.lab = 1.3, cex.axis = 1.3) # Can put all three in one
# op <- par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
# Date linear (Fig. 11 – 20 left)
plot(pm[1:145], 1:145, xlim = c(-1.5, 1.5), xlab = "Parameter estimate",
    ylab = "Species number", main = "Effect of date (linear) on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 1:145], 1:145, cri[2, 1:145], 1:145, col = "grey", lwd = 1)
sig1 <- (cri[1, 1:145] * cri[2, 1:145]) > 0
segments(cri[1, 1:145][sig1 == 1], (1:145)[sig1 == 1], cri[2, 1:145][sig1 == 1],
    (1:145)[sig1 == 1], col = "blue", lwd = 2)
# ~~~~ The order of elements in jags output has changed; use names instead of numerical indices ~~~~
# abline(v = out101$summary[11,1], lwd = 3, col = "red")
# abline(v = out101$summary[11,c(3,7)], lwd = 2, col = "red", lty = 2)
abline(v = out101$summary['mu.betalp1','mean'], lwd = 3, col = "red")
abline(v = out101$summary['mu.betalp1',c('2.5%', '97.5%')], lwd = 2,
    col = "red", lty = 2)

# Date quadratic (not shown)
plot(pm[216:360], 1:145, xlim = c(-1.5, 1.5), xlab = "Parameter estimate",
    ylab = "Species number", main = "Effect of date (quadratic) on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 216:360], 1:145, cri[2, 216:360], 1:145, col = "grey", lwd = 1)
sig2 <- (cri[1, 216:360] * cri[2, 216:360]) > 0
segments(cri[1, 216:360][sig2 == 1], (1:145)[sig2 == 1], cri[2, 216:360][sig2 == 1],
    (1:145)[sig2 == 1], col = "blue", lwd = 2)
# abline(v = out101$summary[13,1], lwd = 3, col = "red")
# abline(v = out101$summary[13, c(3,7)], lwd = 3, col = "red", lty = 2)
abline(v = out101$summary['mu.betalp2','mean'], lwd = 3, col = "red")
abline(v = out101$summary['mu.betalp2',c('2.5%', '97.5%')], lwd = 2,
    col = "red", lty = 2)

# Survey duration (Fig. 11-20 right)
plot(pm[431:575], 1:145, xlim = c(-0.5, 1), xlab = "Parameter estimate",
    ylab = "Species number", main = "Effect of survey duration on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 431:575], 1:145, cri[2, 431:575], 1:145, col = "grey", lwd = 1)
sig3 <- (cri[1, 431:575] * cri[2, 431:575]) > 0
segments(cri[1, 431:575][sig3 == 1], (1:145)[sig3 == 1], cri[2, 431:575][sig3 == 1],
    (1:145)[sig3 == 1], col = "blue", lwd = 2)
# abline(v = out101$summary[15,1], lwd = 3, col = "red")
# abline(v = out101$summary[15, c(3,7)], lwd = 3, col = "red", lty = 2)
abline(v = out101$summary['mu.betalp3','mean'], lwd = 3, col = "red")
abline(v = out101$summary['mu.betalp3',c('2.5%', '97.5%')], lwd = 2,
    col = "red", lty = 2)
par(op)

# Effects of elevation (linear and quadratic) and of forest on occupancy - figures 11.21 - 11.23
op <- par(mfrow = c(1,3), cex.lab = 1.3, cex.axis = 1.3) # can do all in one
# Effect of elevation (linear) on occupancy probability (Fig. 11-21)
plot(pm[646:790], 1:145, xlim = c(-8, 8), xlab = "Parameter estimate",
    ylab = "Species number", main = "Effect of elevation (linear) on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 646:790], 1:145, cri[2, 646:790], 1:145, col = "grey", lwd = 1)
sig4 <- (cri[1, 646:790] * cri[2, 646:790]) > 0
segments(cri[1, 646:790][sig4 == 1], (1:145)[sig4 == 1], cri[2, 646:790][sig4 == 1],
    (1:145)[sig4 == 1], col = "blue", lwd = 2)
# abline(v = out101$summary[3,1], lwd = 3, col = "red")
# abline(v = out101$summary[3,c(3,7)], lwd = 3, col = "red", lty = 2)
abline(v = out101$summary['mu.betalpsi1','mean'], lwd = 3, col = "red")
abline(v = out101$summary['mu.betalpsi1',c('2.5%', '97.5%')], lwd = 2, col = "red", lty = 2)

# Effect of elevation (quadratic) on occupancy probability (Fig. 11-22)
plot(pm[861:1005], 1:145, xlim = c(-4, 2), xlab = "Parameter estimate",
    ylab = "Species number", main = "Effect of elevation (quadratic) on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 861:1005], 1:145, cri[2, 861:1005], 1:145, col = "grey", lwd=1)
sig5 <- (cri[1, 861:1005] * cri[2, 861:1005]) > 0
segments(cri[1, 861:1005][sig5 == 1], (1:145)[sig5 == 1], cri[2, 861:1005][sig5 == 1],
    (1:145)[sig5 == 1], col = "blue", lwd = 2)
# abline(v = out101$summary[5,1], lwd = 3, col = "red")
# abline(v = out101$summary[5,c(3,7)], lwd = 3, col = "red", lty = 2)
abline(v = out101$summary['mu.betalpsi2','mean'], lwd = 3, col = "red")
abline(v = out101$summary['mu.betalpsi2',c('2.5%', '97.5%')], lwd = 2,
    col = "red", lty = 2)


# Effect of forest (linear) on occupancy probability (Fig. 11-23)
plot(pm[1076:1220], 1:145, xlim = c(-3, 4), xlab = "Parameter estimate",
    ylab = "Species number", main = "Effect of forest cover on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 1076:1220], 1:145, cri[2, 1076:1220],1:145, col = "grey", lwd=1)
sig6 <- (cri[1, 1076:1220] * cri[2, 1076:1220]) > 0
segments(cri[1, 1076:1220][sig6 == 1], (1:145)[sig6 == 1], cri[2, 1076:1220][sig6 == 1],
    (1:145)[sig6 == 1], col = "blue", lwd = 2)
# abline(v = out101$summary[7,1], lwd = 3, col = "red")
# abline(v = out101$summary[7,c(3,7)], lwd = 3, col = "red", lty = 2)
abline(v = out101$summary['mu.betalpsi3','mean'], lwd = 3, col = "red")
abline(v = out101$summary['mu.betalpsi3',c('2.5%', '97.5%')], lwd = 2,
    col = "red", lty = 2)
par(op)
negsig6 <- (cri[1, 1076:1220] < 0 & cri[2, 1076:1220] < 0) == 1 # sig negative
possig6 <- (cri[1, 1076:1220] > 0 & cri[2, 1076:1220] > 0) == 1 # sig positive

# Predict detection for date and duration and occupancy for elevation and forest
# for each of the 145 observed species
predS <- array(NA, dim = c(500, nspec, 4))   # covariate value x species x response, "S" for 'species'
p.coef <- cbind(lp=pm[1292:1436], betalp1 = pm[1:145], betalp2 = pm[216:360], betalp3 = pm[431:575])
psi.coef <- cbind(lpsi=pm[1507:1651], betalpsi1 = pm[646:790], betalpsi2 = pm[861:1005], betalpsi3 = pm[1076:1220])

for(i in 1:nspec){          # Loop over 145 observed species
   predS[,i,1] <- plogis(p.coef[i,1] + p.coef[i,2] * dat.pred +
     p.coef[i,3] * dat.pred^2 )     # p ~ date
   predS[,i,2] <- plogis(p.coef[i,1] + p.coef[i,4] * dur.pred) # p ~ duration
   predS[,i,3] <- plogis(psi.coef[i,1] + psi.coef[i,2] * ele.pred +
     psi.coef[i,3] * ele.pred^2 )     # psi ~ elevation
   predS[,i,4] <- plogis(psi.coef[i,1] + psi.coef[i,4] * for.pred) # psi ~ forest
}

# Plots for detection probability and survey date and duration (Fig. 11-24)
op <- par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.dat, predS[,1,1], lwd = 3, type = 'l', lty = 1, frame = FALSE,
   ylim = c(0, 1), xlab = "Survey date (1 = 1 April)",
   ylab = "Detection probability")
for(i in 2:145){
   lines(o.dat, predS[,i,1], col = i, lwd = 3)
}

plot(o.dur, predS[,1,2], lwd = 3, type = 'l', lty = 1, frame = FALSE,
   ylim = c(0, 1), xlab = "Survey duration (min)",
   ylab = "Detection probability")
for(i in 2:145){
   lines(o.dur, predS[,i,2], col = i, lwd = 3)
}
par(op)

# Plots for occupancy probability and elevation and forest cover (Fig. 11-25)
op <- par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.ele, predS[,1,3], lwd = 3, type = 'l', lty = 1, frame = FALSE,
   ylim = c(0, 1), xlab = "Elevation (m a.s.l.)",
   ylab = "Occupancy probability")
for(i in 2:145){
   lines(o.ele, predS[,i,3], col = i, lwd = 3)
}

plot(o.for, predS[,1,4], lwd = 3, type = 'l', lty = 1, frame = FALSE,
   ylim = c(0, 1), xlab = "Forest cover (%)", ylab = "Occupancy probability")
for(i in 2:145){
   lines(o.for, predS[,i,4], col = i, lwd = 3)
}
par(op)
