#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

# Approximate execution time for this code: 10 mins
# Run time with the full number of iterations: 85 mins

library(AHMbook) # Requires version ‘0.1.4.9088’ or later
library(jagsUI)

# ~~~~~~~ changes to RNG defaults ~~~~~~~~~~~~~~~~~~~~~~~~
# Use the old default random number generator to get the printed numbers
RNGversion("3.2.0")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 10.14 Modeling wiggly covariate relationships: penalized splines in hierarchical models
# =======================================================================================


# Execute the function and inspect file produced
data <- wigglyOcc(seed = 1)
str(data)

# Convert matrix data into vectors and prepare stuff
y <- c(data$y)                      # Detection/nondetection data (response)
Xsite <- data$Xsite                 # Fine as is
Xsurvey <- c(data$Xsurvey)          # Survey covariate
site <- rep(1:data$M, data$J)       # Site index

# tmp1 <-spline.prep(Xsite, 20) # This would choose 20 knots for Xsite
tmp1 <-spline.prep(Xsite, NA)   # Choose variable default number of knots
tmp2 <-spline.prep(Xsurvey, NA)
Xocc <- tmp1$X         # Fixed-effects part of covariate Xsite in occ
Zocc <- tmp1$Z         # Random-effects part of covariate Xsite in occ
Xdet <- tmp2$X         # Fixed-effects part of covariate Xsite in det
Zdet <- tmp2$Z         # Random-effects part of covariate Xsite in det
nk.occ <- length(tmp1$knots)    # Number of knots in occupancy spline
nk.det <- length(tmp2$knots)    # Number of knots in detection spline


# Bundle and summarize data set
win.data <- list(y1 = y, site = site, M = data$M, Xocc = Xocc, Zocc = Zocc,
    nk.occ = nk.occ, Xdet = Xdet, Zdet = Zdet, nk.det = nk.det,
    nobs = data$M*data$J, y2 = y, onesSite = rep(1, 240),
    onesSurvey = rep(1, 720), Xsite = Xsite, Xsite2 = Xsite^2,
    Xsurvey = Xsurvey, Xsurvey2 = Xsurvey^2)
str(win.data)    # onesSite and onesSurvey are for occ and det intercepts

# Specify two models in one in BUGS language
cat(file = "hypermodel.txt",
"model {

  # *** Spline model for the data***
  # --------------------------------
  # Priors
  for(k in 1:3){                 # Regression coefficients
    alpha1[k] ~ dnorm(0, 0.1)   # Detection model
    beta1[k] ~ dnorm(0, 0.1)    # Occupancy model
  }
  for(k in 1:nk.occ){ # Random effects at specified knots (occupancy)
    b.occ[k] ~ dnorm(0, tau.b.occ)
  }
  for(k in 1:nk.det){ # Random effects at specified knots (detection)
    b.det[k] ~ dnorm(0, tau.b.det)
  }
  tau.b.occ ~ dgamma(0.01, 0.01)
  tau.b.det ~ dgamma(0.01, 0.01)

  # Likelihood
  # Model for latent occupancy state
  for (i in 1:M) {
    z1[i] ~ dbern(psi1[i])
    logit(psi1[i]) <- fix.terms.occ[i] + smooth.terms.occ[i]
    fix.terms.occ[i] <- beta1[1]*Xocc[i,1] + beta1[2]*Xocc[i,2] + beta1[3]*Xocc[i,2]
    smooth.terms.occ[i] <- inprod(b.occ[], Zocc[i,])
  }
  # Model for observations
  for(i in 1:nobs){
    y1[i] ~ dbern(mu.y1[i])
    mu.y1[i] <- z1[site[i]] * p1[i]
    logit(p1[i]) <- fix.terms.det[i] + smooth.terms.det[i]
    fix.terms.det[i] <- alpha1[1]*Xdet[i,1] + alpha1[2]*Xdet[i,2] +
        alpha1[3]*Xdet[i,2]
    smooth.terms.det[i] <- inprod(b.det[], Zdet[i,])
  }

  # Derived quantities
  sum.z1 <- sum(z1[])            # Number of occupied sites in sample
  sd.b.occ <- sqrt(1/tau.b.occ)  # SD of spline random effects variance Occ.
  sd.b.det <- sqrt(1/tau.b.det)  # SD of spline random effects variance Det.


  # *** Polynomial model for same data ***
  # --------------------------------------
  # Priors
  for(k in 1:3){                 # Regression coefficients
    alpha2[k] ~ dnorm(0, 0.1)  # Detection model
    beta2[k] ~ dnorm(0, 0.1)   # Occupancy model
  }

  # Likelihood
  # Model for latent occupancy state
  for (i in 1:M) {
    z2[i] ~ dbern(psi2[i])
    logit(psi2[i]) <- beta2[1]*onesSite[i] + beta2[2]*Xsite[i] + beta2[3]*Xsite2[i]
  }
  # Model for observations
  for(i in 1:nobs){
    y2[i] ~ dbern(mu.y2[i])
    mu.y2[i] <- z2[site[i]] * p2[i]
    logit(p2[i]) <- alpha2[1]*onesSurvey[i] + alpha2[2]*Xsurvey[i] +
        alpha2[3] * Xsurvey2[i]
  }

  # Derived quantities
  sum.z2 <- sum(z2[])            # Number of occupied sites in sample
}
")

# Initial values
zst <- apply(data$y, 1, max)
inits <- function(){list(z1=zst, alpha1=rnorm(3), beta1=rnorm(3),
    b.occ = rnorm(nk.occ), b.det = rnorm(nk.det), tau.b.occ = runif(1),
    tau.b.det = runif(1), z2=zst, alpha2=rnorm(3), beta2=rnorm(3))}

# Parameters monitored
params <- c("alpha1", "beta1", "psi1", "p1", "fix.terms.occ",
    "smooth.terms.occ", "b.occ", "fix.terms.det", "smooth.terms.det", "b.det",
    "sum.z1", "sd.b.occ", "sd.b.det", "alpha2", "beta2", "psi2", "p2", "sum.z2")

# MCMC settings
# ni <- 100000   ;   nb <- 10000   ;   nt <- 90   ;   nc <- 3
ni <- 10000   ;   nb <- 1000   ;   nt <- 10   ;   nc <- 3  # ~~~~ for testing

# Call JAGS from R (ART 95 min) and summarize posteriors
system.time(fhm <- jags(win.data, inits, params, "hypermodel.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE))
traceplot(fhm)   ;   print(fhm, 3)
print(fhm$summary[c(1:6, 2957:2965, 3926),], 3)  # Compare some key estimands

# Plot prediction of psi and p
op <- par(mfrow = c(1,2), mar = c(5,4,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(Xsite, data$psi, main = "Occupancy probability",type = "l",
    ylim = c(-0.1, 1.1), col = "red", xlab = "Site covariate (Xsite)",
    ylab = "", lwd = 3)
#abline(v = tmp1$knots, col = "grey")
points(Xsite, jitter(data$z, amount = 0.02))
lines(Xsite, fhm$mean$psi1, col = "blue", lty = 1, lwd = 3)
lines(Xsite, fhm$mean$psi2, col = "brown", lty = 1, lwd = 3)

plot(Xsurvey[order(data$x.index)], data$p.ordered,
    main = "Detection probability ",type = "l", ylim = c(-0.1, 1.1),
    col = "red", xlab = "Survey covariate (Xsurvey)", ylab = "", lwd = 3)
points(Xsurvey, jitter(y, amount = 0.02))
#abline(v = tmp2$knots, col = "grey")
lines(Xsurvey[order(data$x.index)], fhm$mean$p1[order(data$x.index)],
    col = "blue", lwd = 3)
lines(Xsurvey[order(data$x.index)], fhm$mean$p2[order(data$x.index)],
    col = "brown", lwd = 3)
par(op)
