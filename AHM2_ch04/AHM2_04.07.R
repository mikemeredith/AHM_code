#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from MS dated 2019-08-28, no code in the book.

library(AHMbook)
library(unmarked)

# 4.7 Study design, bias, and precision of estimators
# ===================================================

# 4.7.1 Can we fit dynocc models to single-visit data ?
# ------------------------------------------------------

# Case 1 : Fully time-dependent model
# '''''''''''''''''''''''''''''''''''

# Do simulation with 100 reps (takes about 2 hours)
# simrep <- 100
simrep <- 10  # ~~~ for testing

# Define arrays to hold true values and the results
true.vals1 <- array(dim = c(39, simrep))
estimates1 <- array(dim = c(39, simrep))

# Start simulation
system.time(                  # time whole thing: 101 mins
for(k in 1:simrep){           # Loop over k simreps
  # Counter
  cat("** simrep", k, "***\n")
  # Generate a data set using simDynocc
  pick.psi1 <- runif(1, 0.01, 0.99)
  data <- simDynocc(nsites = 267, nsurveys = 1, nyears = 10,
    mean.psi1 = pick.psi1, range.phi = c(0.01, 0.99),
    range.gamma = c(0.01, 0.99), range.p = c(0.01, 0.99),
    beta.Xp = 0, trend.sd.site = c(1, 1), trend.sd.survey = c(1, 1),
    show.plot = FALSE)
  # Fit model
  yy <- matrix(data$y, data$nsites, data$nsurveys * data$nyears)
  yr <- matrix(c('01','02','03','04','05','06','07','08','09','10'),
    nrow(yy), data$nyears, byrow=TRUE)
  umf <- unmarkedMultFrame(y=yy, yearlySiteCovs=list(year=yr),
        numPrimary=data$nyears)
  summary(fm <- colext(~1, ~year-1, ~year-1, ~year-1, data = umf,
    control=list(trace=TRUE, REPORT=20, maxit = 250), se = F) )
  # Save results (true parameter values, MLEs and projected psi)
  true.vals1[,k] <- c(data$mean.psi1, data$mean.gamma, 1 - data$mean.phi,
      data$mean.p, data$mean.psi)
  estimates1[,k] <- c(coef(fm), projected(fm)[2,])
}
)  # 10 took 4.5 mins
rownames(true.vals1) <- rownames(estimates1) <- c(names(coef(fm)),
    paste0('proj', 1:10))

# Case 2 : Intercept only model
# '''''''''''''''''''''''''''''

# Do simulation with 1000 reps (takes about 6 mins)
# simrep <- 100
simrep <- 10  # ~~~ for testing

# Define arrays to hold true values and the results
true.vals2 <- array(dim = c(4, simrep))
estimates2 <- array(dim = c(4, simrep))

# Start simulation
system.time(                  # time whole thing: about 6 minutes
for(k in 1:simrep){           # Loop k over simreps
  # Counter
  cat("** simrep", k, "***\n")
  # Generate a data set: pick p and use p in simDynocc()
  pick.psi1 <- runif(1, 0.01, 0.99)
  pick.phi <- sort(runif(1, 0.01, 0.99))
  pick.gamma <- sort(runif(1, 0.01, 0.99))
  pick.p <- sort(runif(1, 0.01, 0.99))

  data <- simDynocc(nsites = 267, nsurveys = 1, nyears = 10,
    mean.psi1 = pick.psi1, range.phi = c(pick.phi, pick.phi),
    range.gamma = c(pick.gamma, pick.gamma), range.p = c(pick.p, pick.p),
    beta.Xp = 0, trend.sd.site = c(1, 1), trend.sd.survey = c(1, 1),
    show.plot = FALSE)

  # Fit model
  yy <- matrix(data$y, data$nsites, data$nsurveys * data$nyears)
  yr <- matrix(c('01','02','03','04','05','06','07','08','09','10'),
    nrow(yy), data$nyears, byrow=TRUE)
  umf <- unmarkedMultFrame(y=yy, yearlySiteCovs=list(year=yr),
        numPrimary=data$nyears)
  summary(fm <- colext(~1, ~1, ~1, ~1, data = umf,
    control=list(trace=TRUE, REPORT=20, maxit = 250), se = F) )

  # Save results (true parameter values, MLEs and projected psi)
  true.vals2[,k] <- c(data$mean.psi1, data$mean.gamma[1],
      1 - data$mean.phi[1], data$mean.p[1])
  estimates2[,k] <-coef(fm)
}
)  # 10 took 14 secs
rownames(true.vals2) <- rownames(estimates2) <- names(coef(fm))

# Case 3 : Additional random noise in p
# '''''''''''''''''''''''''''''''''''''

# Do simulation with 100 reps (takes about 40 minutes)
# simrep <- 100
simrep <- 10  # ~~~ for testing

# Define arrays to hold true values and the results
true.vals3 <- array(dim = c(40, simrep))
estimates3 <- array(dim = c(40, simrep))

# Start simulation
system.time(                  # time whole thing
for(k in 1:simrep){           # Loop k over simreps
  # Counter
  cat("** simrep", k, "***\n")
  pick.psi1 <- runif(1, 0.01, 0.99)

  # Generate a data set using simDynocc()
  data <- simDynocc(nsites = 267, nsurveys = 1, nyears = 10,
    mean.psi1 = pick.psi1, range.phi = c(0.01, 0.99),
    range.gamma = c(0.01, 0.99), range.p = c(0.01, 0.99),
    beta.Xp = 1, trend.sd.site = c(1, 1), trend.sd.survey = c(1, 1),
    show.plot = FALSE)

  # Fit model
  yy <- matrix(data$y, data$nsites, data$nsurveys * data$nyears)
  Xp <- matrix(data$Xp, data$nsites, data$nsurveys * data$nyears)
  yr <- matrix(c('01','02','03','04','05','06','07','08','09','10'),
    nrow(yy), data$nyears, byrow=TRUE)
  umf <- unmarkedMultFrame(y=yy, obsCovs = list(Xp = Xp),
   yearlySiteCovs=list(year=yr), numPrimary=data$nyears)
  summary(fm <- colext(~1, ~year-1, ~year-1, ~ year + Xp - 1, data = umf,
    control=list(trace=TRUE, REPORT=20, maxit = 250), se = F) )

  # Save results (true parameter values, MLEs and projected psi)
  true.vals3[,k] <- c(data$mean.psi1, data$mean.gamma, 1 - data$mean.phi,
      data$mean.p, data$beta.Xp, data$mean.psi)
  estimates3[,k] <- c(coef(fm), projected(fm)[2,])
}
)  # 10 took 3 mins
rownames(true.vals3) <- rownames(estimates3) <- c(names(coef(fm)),
    paste0('proj', 1:10))

# Figure 4.8
# ''''''''''
op <- par(mfrow = c(3, 3), mar = c(5,5,3,2), cex.lab = 1.2)

# Full time-dependent model
# Occupancy
plot(true.vals1[30:39,], estimates1[30:39,], xlab = "True occupancy prob.",
    ylab = "Occupancy estimate", xlim = c(0, 1), ylim = c(0,1),
    main = "Full time-dep model", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(estimates1[30:39,] ~ true.vals1[30:39,], df = 5),
    col = "blue", lwd = 2)

# Colonization
plot(true.vals1[2:10,], plogis(estimates1[2:10,]),
    xlab = "True colonization prob.", ylab = "Colonization estimate",
    xlim = c(0, 1), ylim = c(0,1), main = "Full time-dep model", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(plogis(estimates1[2:10,]) ~ true.vals1[2:10,], df = 5),
    col = "blue", lwd = 2)

# Detection
plot(true.vals1[20:29,], plogis(estimates1[20:29,]),
    xlab = "True detection prob.", ylab = "Detection estimate",
    xlim = c(0, 1), ylim = c(0,1), main = "Full time-dep model", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(plogis(estimates1[20:29,]) ~ true.vals1[20:29,], df = 5),
    col = "blue", lwd = 2)


# Intercepts-only model
# Occupancy
plot(true.vals2[1,], plogis(estimates2[1,]),
    xlab = "True initial occ. prob.", ylab = "Initial occupancy estimate",
    xlim = c(0, 1), ylim = c(0,1), main = "Intercept-only model", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(plogis(estimates2[1,]) ~ true.vals2[1,], df = 5),
    col = "blue", lwd = 2)

# Colonization
plot(true.vals2[2,], plogis(estimates2[2,]),
    xlab = "True colonization prob.", ylab = "Colonization estimate",
    xlim = c(0, 1), ylim = c(0,1), main = "Intercept-only model", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(plogis(estimates2[2,]) ~ true.vals2[2,], df = 5),
    col = "blue", lwd = 2)

# Detection
plot(true.vals2[4,], plogis(estimates2[4,]),
    xlab = "True detection prob.", ylab = "Detection estimate",
    xlim = c(0, 1), ylim = c(0,1), main = "Intercept-only model", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(plogis(estimates2[4,]) ~ true.vals2[4,], df = 5),
    col = "blue", lwd = 2)


# Full time-dependent model with one site-covariate in p
# Occupancy
plot(true.vals3[31:40,], estimates3[31:40,], xlab = "True occupancy prob.",
    ylab = "Occupancy estimate", xlim = c(0, 1), ylim = c(0,1),
    main = "Full time-dep model with covariate in p", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(estimates3[31:40,] ~ true.vals3[31:40,], df = 5),
    col = "blue", lwd = 2)

# Colonization
plot(true.vals3[2:10,], plogis(estimates3[2:10,]),
    xlab = "True colonization prob.", ylab = "Colonization estimate",
    xlim = c(0, 1), ylim = c(0,1),
    main = "Full time-dep model with covariate in p", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(plogis(estimates3[2:10,]) ~ true.vals3[2:10,], df = 5),
    col = "blue", lwd = 2)

# Detection
plot(true.vals3[20:29,], plogis(estimates3[20:29,]),
    xlab = "True detection prob.", ylab = "Detection estimate",
    xlim = c(0, 1), ylim = c(0,1),
    main = "Full time-dep model with covariate in p", frame = FALSE)
abline(0, 1, col='red', lwd = 2)
lines(smooth.spline(plogis(estimates3[20:29,]) ~ true.vals3[20:29,], df = 5),
    col = "blue", lwd = 2)
par(op)


# 4.7.2 Bias and precision as a function of nsite, nsurvey and p
# ---------------------------------------------------------------

# Do simulation with 100 reps (takes about 2 hours)
# simrep <- 100
simrep <- 10  # ~~~ for testing

# Choose number and levels of simulation factors
Msim <- c(20, 120, 250)     # number of sites
Jsim <- c(2, 5, 10)         # number of repeat surveys

# Define arrays to hold true values and the results (previous fit req'd)
true.vals <- array(dim = c(39, 3, 3, simrep))
estimates <- array(dim = c(39, 3, 3, simrep))

# Start simulation
system.time(                  # time whole thing
for(k in 1:simrep){           # Loop k over simreps
  for(i in 1:3) {             # Loop i over site factor (Msim)
    for(j in 1:3) {           # Loop j over survey factor (Jsim)
      # Counter
      cat("** nsites", Msim[i], ", nsurveys", Jsim[j], ", simrep", k, "***\n")
      # Generate a data set: pick p and use p in simDynocc()
      p.draw <- runif(1, 0.01, 0.99)
      data <- simDynocc(nsite = Msim[i], nsurvey = Jsim[j], nyear = 10,
        mean.psi1 = 0.6, range.phi = c(0.1, 0.9),
        range.gamma = c(0.1, 0.9), range.p = c(p.draw, p.draw),
        show.plot = FALSE)
      # Fit model
      yy <- matrix(data$y, data$nsite, data$nsurvey * data$nyear)
      yr <- matrix(c('01','02','03','04','05','06','07','08','09','10'),
        nrow(yy), data$nyear, byrow=TRUE)
      umf <- unmarkedMultFrame(y=yy, yearlySiteCovs=list(year=yr),
        numPrimary=data$nyear)
      summary(fm <- colext(~1, ~year-1, ~year-1, ~year-1, data = umf,
        control=list(trace=TRUE, REPORT=20, maxit = 250), se = F) )

      # Save results (true parameter values, MLEs and projected psi)
      true.vals[,i,j,k] <- c(data$mean.psi1, data$mean.gamma,
          1 - data$mean.phi, data$mean.p, data$mean.psi)
      estimates[,i,j,k] <- c(coef(fm), projected(fm)[2,])
    }
  }
}
)
dimnames(true.vals) <- dimnames(estimates) <- list(c(names(coef(fm)),
    paste0('proj', 1:10)), Msim, Jsim, NULL)

# Compute relative error in percent (for occupancy probability)
error.psi.hat <- 100 * (estimates[30:39,,,] - true.vals[30:39,,,]) /
    true.vals[30:39,,,]

# Figure 4.9
op <- par(mfrow = c(3,3), mar = c(5,5,3,1), cex.main = 1.2,
    cex.lab = 1.2, cex.axis = 1.2)
for(i in 1:3){
  for(j in 1:3){
    lab <- paste(Msim[i],"sites,", Jsim[j],"surveys")
    plot(true.vals[20:29,i,j,], error.psi.hat[,i,j,],
    xlab = "Detection prob.", ylab = "% Error in Occ. prob.",
    main = lab, xlim = c(0, 1), ylim = c(-200,200))
    abline(h = 0, col = "red", lwd = 2)
    lines(smooth.spline(error.psi.hat[,i,j,] ~ true.vals[20:29,i,j,],
    df = 5), col = "blue", lwd = 2)
   }
}
par(op)

# 4.7.3 A power analysis for occupancy trend estimation
# -----------------------------------------------------
# Design points of factorial simulation design and number of sims
Msim <- c(20, 50, 100, 250)
Jsim <- c(2, 5, 10, 25)
# simrep <- 1000          # takes about 70 min
simrep <- 100  # ~~~ for testing

# Data structure for results
results <- array(NA, dim = c(length(Msim), length(Jsim), 2, simrep))
dimnames(results) <- list(Msim, Jsim, c("Det p", "LRT p"), NULL)

system.time(       # Being Swiss .... and obsessed with time...
# Loop over settings for M
for(i in 1:length(Msim)){
  # Loop over settings for J
  for(j in 1:length(Jsim)){
    cat(paste("\n*** nsite:", Msim[i], "; nsurvey:", Jsim[j], "***\n"))
    # Loop over simulation reps
    for(k in 1:simrep){
      cat(paste("survey:", k, "\n"))
      # Simulate data set as a realization of above process
      p.draw <- runif(1, 0.01, 0.99)   # Draw a value of detection prob.
      data <- simDynocc(nsite = Msim[i], nsurvey = Jsim[j], nyear = 10,
        mean.psi1 = 0.62, range.phi = c(0.9, 0.9),
        range.gamma = c(0.1, 0.1), range.p = c(p.draw, p.draw),
        show.plot = FALSE)

      # Estimate trend using stacked data and static occupancy model
      yrv <- rep(1:data$nyear, each = data$nsite)
      ystack <- array(NA, dim = c(data$nsite * data$nyear, data$nsurvey))
      for(t in 1:data$nyear){
        ystack[((t-1)*data$nsite+1):(data$nsite*t),] <- data$y[,,t]
      }
      umf <- unmarkedFrameOccu(y=ystack, siteCovs=data.frame(year=yrv))
      fm0 <- occu(~1 ~ 1, data=umf, se = F)     # Null model
      fmT <- occu(~1 ~ year, data=umf, se = F)  # Trend model
      lrt <- LRT(fm0, fmT)  # likelihood ratio test for trend

      # Save value of detection probability and LRT result
      results[i,j,1,k] <- p.draw
      results[i,j,2,k] <- lrt[['Pr(>Chisq)']]
    }
  }
} )  # 75 mins (100 took 7.6 mins)

power <- array(NA, c(4,4))
rownames(power) <- c('nsite = 20', 'nsite = 50', 'nsite = 100', 'nsite = 250')
colnames(power) <- c('nsurvey = 2', 'nsurvey = 5', ' nsurvey = 10',
    ' nsurvey = 25')

for(i in 1:length(Msim)){
  for(j in 1:length(Jsim)){
    power[i,j] <- mean(results[i,j,2,] < 0.05)
  }
}
power
            # nsurvey = 2 nsurvey = 5  nsurvey = 10  nsurvey = 25
# nsite = 20        0.173       0.219         0.233         0.247
# nsite = 50        0.218       0.298         0.361         0.378
# nsite = 100       0.317       0.427         0.469         0.523
# nsite = 250       0.563       0.687         0.752         0.795


# 4.7.4 Effects of unmodelled detection heterogeneity
# ---------------------------------------------------

# Figure 4.10
base.p <- c(0.2, 0.5, 0.8)
e1 <- rnorm(10^7, 0, 0.2)
e2 <- rnorm(10^7, 0, 1)
e3 <- rnorm(10^7, 0, 2)

op <- par(mfrow = c(1, 3))
set1 <- plogis(qlogis(base.p[1]) + e1)
set2 <- plogis(qlogis(base.p[2]) + e1)
set3 <- plogis(qlogis(base.p[3]) + e1)
dens1 <- density(set1) ; dens2 <- density(set2) ; dens3 <- density(set3)
plot(dens1$x, dens1$y / max(dens1$y), col = 'black', xlim = c(0,1),
    ylim = c(0,1), main = 'Heterogeneity SD = 0.2',
    xlab = 'Detection probability', ylab = '',
    type = 'l', yaxt = 'n', frame=FALSE, lwd=2)
lines(dens2$x, dens2$y / max(dens2$y), col = 'black', lwd=2)
lines(dens3$x, dens3$y / max(dens3$y), col = 'black', lwd=2)

set1 <- plogis(qlogis(base.p[1]) + e2)
set2 <- plogis(qlogis(base.p[2]) + e2)
set3 <- plogis(qlogis(base.p[3]) + e2)
dens1 <- density(set1) ; dens2 <- density(set2) ; dens3 <- density(set3)
plot(dens1$x, dens1$y / max(dens1$y), col = 'black', xlim = c(0,1),
    ylim = c(0,1), main = 'Heterogeneity SD = 1',
    xlab = 'Detection probability', ylab = '',
    type = 'l', frame = FALSE, lwd = 2, yaxt = 'n')
lines(dens2$x, dens2$y / max(dens2$y), col = 'black', lwd = 2)
lines(dens3$x, dens3$y / max(dens3$y), col = 'black', lwd = 2)

set1 <- plogis(qlogis(base.p[1]) + e3)
set2 <- plogis(qlogis(base.p[2]) + e3)
set3 <- plogis(qlogis(base.p[3]) + e3)
dens1 <- density(set1) ; dens2 <- density(set2) ; dens3 <- density(set3)
plot(dens1$x, dens1$y / max(dens1$y), col = 'black', xlim = c(0,1),
    ylim = c(0,1), main = 'Heterogeneity SD = 2',
    xlab = 'Detection probability', ylab = '', type = 'l',
    frame = FALSE, lwd = 2, yaxt = 'n')
lines(dens2$x, dens2$y / max(dens2$y), col = 'black', lwd = 2)
lines(dens3$x, dens3$y / max(dens3$y), col = 'black', lwd = 2)
par(op)


# Set up simulation (sd.site = 0.2, i.e., little heterogeneity)
# ------------------------------------------------------------
nsim <- 100                         # Number of sim reps
psi0.2 <- eps0.2 <- gamma0.2 <- p0.2 <- psi.hat0.2.ml <- eps.hat0.2.ml <-
    gamma.hat0.2.ml <- p.hat0.2.ml <- array(NA, dim = c(8, nsim))

system.time(
for(i in 1:nsim){
  cat("\n\n*** Simrep number", i, "***\n\n")

  # Generate data set and save time-specific parameters
  pick.psi1 <- runif(1, 0.01, 0.99)
  data <- simDynocc(nsite = 250, nsurvey = 3, nyear = 8, mean.psi1 = pick.psi1, range.p = c(0.01, 0.99), range.phi = c(0.01, 0.99), range.gamma = c(0.01, 0.99), trend.sd.site = c(0.2, 0.2), show.plot = F)
  yy <- matrix(data$y,
  data$nsite, data$nsurvey* data$nyear) # Format data sideways
  year <- matrix(factor(1:data$nyear), nrow(yy), data$nyear, byrow=TRUE)
  psi0.2[,i] <- data$mean.psi
  eps0.2[1:(data$nyear-1),i] <- 1-data$mean.phi
  gamma0.2[1:(data$nyear-1),i] <- data$mean.gamma
  p0.2[,i] <- data$mean.p

  # Fit standard model by ML
  simUMF <- unmarkedMultFrame(y = yy, yearlySiteCovs = list(year = year), numPrimary=data$nyear)
  summary(fm <- colext(~1, ~ year-1, ~ year-1, ~ year-1, data = simUMF,
      control=list(trace=TRUE, REPORT=10), se = FALSE))

  # Get time-specific MLEs and save them
  nd1 <- data.frame(year= factor(1:(data$nyear-1))) # Year cov for eps, gamma
  nd2 <- data.frame(year=factor(1:data$nyear)) # Year cov for psi, p
  psi.hat0.2.ml[,i] <- projected(fm)[2,]
  eps.hat0.2.ml[1:(data$nyear-1),i] <- predict(fm, type='ext', newdata=nd1)[,1]
  gamma.hat0.2.ml[1:(data$nyear-1),i] <-
      predict(fm, type='col', newdata=nd1)[,1]
  p.hat0.2.ml[,i] <- predict(fm, type='det', newdata=nd2)[,1]
} )  # 17 mins, 2 did not converge


# Set up simulation (sd.site = 1)
# ---------------------------------------
nsim <- 100                         # Number of sim reps
psi1 <- eps1 <- gamma1 <- p1 <- psi.hat1.ml <- eps.hat1.ml <- gamma.hat1.ml <-
    p.hat1.ml <- array(NA, dim = c(8, nsim))
    # ml for MLEs and b for Bayesian MCMC posterior means
system.time(
for(i in 1:nsim){
  cat("\n\n*** Simrep number", i, "***\n\n")

  # Generate data set and save time-specific parameters
  pick.psi1 <- runif(1, 0.01, 0.99)
  data <- simDynocc(nsite = 250, nsurvey = 3, nyear = 8, mean.psi1 = pick.psi1,
      range.p = c(0.01, 0.99), range.phi = c(0.01, 0.99), range.gamma = c(0.01, 0.99),
      trend.sd.site = c(1, 1), show.plot = F)
  yy <- matrix(data$y,
  data$nsite, data$nsurvey* data$nyear) # Format data sideways
  year <- matrix(factor(1:data$nyear), nrow(yy), data$nyear, byrow=TRUE)
  psi1[,i] <- data$mean.psi
  eps1[1:(data$nyear-1),i] <- 1-data$mean.phi
  gamma1[1:(data$nyear-1),i] <- data$mean.gamma
  p1[,i] <- data$mean.p

  # Fit standard model by ML
  simUMF <- unmarkedMultFrame(y = yy, yearlySiteCovs = list(year = year),
      numPrimary=data$nyear)
  summary(fm <- colext(~1, ~ year-1, ~ year-1, ~ year-1, data = simUMF,
      control=list(trace=TRUE, REPORT=10), se = FALSE))

  # Get time-specific MLEs and save them
  nd1 <- data.frame(year= factor(1:(data$nyear-1))) # Year cov for eps, gamma
  nd2 <- data.frame(year=factor(1:data$nyear)) # Year cov for psi, p
  psi.hat1.ml[,i] <- projected(fm)[2,]
  eps.hat1.ml[1:(data$nyear-1),i] <- predict(fm, type='ext', newdata=nd1)[,1]
  gamma.hat1.ml[1:(data$nyear-1),i] <- predict(fm, type='col', newdata=nd1)[,1]
  p.hat1.ml[,i] <- predict(fm, type='det', newdata=nd2)[,1]
} )  # 16 mins, 2 convergence failures


# Set up simulation (sd.site = 2)
# ---------------------------------------
nsim <- 100                         # Number of sim reps
psi2 <- eps2 <- gamma2 <- p2 <- psi.hat2.ml <- eps.hat2.ml <- gamma.hat2.ml <-
    p.hat2.ml <- array(NA, dim = c(8, nsim))
system.time(
for(i in 1:nsim){
  cat("\n\n*** Simrep number", i, "***\n\n")

  # Generate data set and save time-specific parameters
  pick.psi1 <- runif(1, 0.01, 0.99)
  data <- simDynocc(nsite = 250, nsurvey = 3, nyear = 8, mean.psi1 = pick.psi1,
      range.p = c(0.01, 0.99), range.phi = c(0.01, 0.99),
      range.gamma = c(0.01, 0.99),
      trend.sd.site = c(2, 2), show.plot = F)
  yy <- matrix(data$y,
  data$nsite, data$nsurvey* data$nyear) # Format data sideways
  year <- matrix(factor(1:data$nyear), nrow(yy), data$nyear, byrow=TRUE)
  psi2[,i] <- data$mean.psi
  eps2[1:(data$nyear-1),i] <- 1-data$mean.phi
  gamma2[1:(data$nyear-1),i] <- data$mean.gamma
  p2[,i] <- data$mean.p

  # Fit standard model by ML
  simUMF <- unmarkedMultFrame(y = yy, yearlySiteCovs = list(year = year), numPrimary=data$nyear)
  summary(fm <- colext(~1, ~ year-1, ~ year-1, ~ year-1, data = simUMF,
      control=list(trace=TRUE, REPORT=10), se = FALSE))

  # Get time-specific MLEs and save them
  nd1 <- data.frame(year= factor(1:(data$nyear-1))) # Year cov for eps, gamma
  nd2 <- data.frame(year=factor(1:data$nyear)) # Year cov for psi, p
  psi.hat2.ml[,i] <- projected(fm)[2,]
  eps.hat2.ml[1:(data$nyear-1),i] <- predict(fm, type='ext', newdata=nd1)[,1]
  gamma.hat2.ml[1:(data$nyear-1),i] <- predict(fm, type='col', newdata=nd1)[,1]
  p.hat2.ml[,i] <- predict(fm, type='det', newdata=nd2)[,1]
} )  # 13 mins

# Figure 4.11
# -----------
# Visualisation of simulation results
# Plots of estimates versus truth (all 4 params, 3 levels of heterogeneity)
op <- par(mfrow = c(4,3), mar = c(5,5,3,3), cex.lab = 1.4, cex.axis = 1.4)
plot(psi0.2, psi.hat0.2.ml, xlab = "Truth", ylab = "Occupancy",
    main = "Heterogeneity sd = 0.2", xlim = c(0,1), ylim = c(0,1),
    frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(psi0.2, psi.hat0.2.ml), col = "blue", lwd = 2)
plot(psi1, psi.hat1.ml, xlab = "Truth", ylab = "Occupancy",
    main = "Heterogeneity sd = 1", xlim = c(0,1), ylim = c(0,1),
    frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(psi1, psi.hat1.ml), col = "blue", lwd = 2)
plot(psi2, psi.hat2.ml, xlab = "Truth", ylab = "Occupancy",
    main = "Heterogeneity sd = 2", xlim = c(0,1), ylim = c(0,1),
    frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(psi2, psi.hat2.ml), col = "blue", lwd = 2)

plot(gamma0.2, gamma.hat0.2.ml, xlab = "Truth", ylab = "Colonization",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(gamma0.2[1:4,], gamma.hat0.2.ml[1:4,]),
    col = "blue", lwd = 2)
plot(gamma1, gamma.hat1.ml, xlab = "Truth", ylab = "Colonization",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(gamma1[1:4,], gamma.hat1.ml[1:4,]), col = "blue", lwd = 2)
plot(gamma2, gamma.hat2.ml, xlab = "Truth", ylab = "Colonization",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(gamma2[1:4,], gamma.hat2.ml[1:4,]), col = "blue", lwd = 2)

plot(eps0.2, eps.hat0.2.ml, xlab = "Truth", ylab = "Extinction",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(eps0.2[1:4,], eps.hat0.2.ml[1:4,]), col = "blue", lwd = 2)
plot(eps1, eps.hat1.ml, xlab = "Truth", ylab = "Extinction",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(eps1[1:4,], eps.hat1.ml[1:4,]), col = "blue", lwd = 2)
plot(eps2, eps.hat2.ml, xlab = "Truth", ylab = "Extinction",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(eps2[1:4,], eps.hat2.ml[1:4,]), col = "blue", lwd = 2)

plot(p0.2, p.hat0.2.ml, xlab = "Truth", ylab = "Detection",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(p0.2, p.hat0.2.ml), col = "blue", lwd = 2)
plot(p1, p.hat1.ml, xlab = "Truth", ylab = "Detection",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(p1, p.hat1.ml), col = "blue", lwd = 2)
plot(p2, p.hat2.ml, xlab = "Truth", ylab = "Detection",
    main = "", xlim = c(0,1), ylim = c(0,1), frame = FALSE)
abline(0, 1, col = "red", lwd = 1)
lines(smooth.spline(p2, p.hat2.ml), col = "blue", lwd = 2)
par(op)
