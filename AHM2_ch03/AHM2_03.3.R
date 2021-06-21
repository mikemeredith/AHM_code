#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 3 : HIERARCHICAL MODELS OF SURVIVAL
# ===========================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 25 mins

library(AHMbook)
library(jagsUI)

# 3.3 Hierarchical CJS models to combine information across, and compare, groups
# ==============================================================================

# 3.3.1 Community CJS models for inference at both the species and the community scale
# ------------------------------------------------------------------------------------

# Pick "hyperdata" for 50 species: hyperparameters which describe the community
nspec <- 50       # choose number of species
mu.lphi <- 0      # mean survival of all species on logit scale
sigma.lphi <- 0.5 # standard deviation of survival on logit scale
mu.lp <- -1       # mean recapture of all species on logit scale
sigma.lp <- 0.5   # standard deviation of recapture on logit scale

# Draw values of survival and recapture (on the logit scale) to describe each species
set.seed(123)
lphi <- rnorm(n = nspec, mean = mu.lphi, sd = sigma.lphi)
lp <- rnorm(n = nspec, mean = mu.lp, sd = sigma.lp)

# Transform to probability scale and print species-level values
phi <- plogis(lphi)
p <- plogis(lp)
round(sort(phi),2)   # not shown
round(sort(p),2)     # not shown

str(ch <- array(NA, dim = c(100, 6, nspec)))

# Simulate data set for 50 species
for(s in 1:nspec){
  # Plug each pair (phi, p) into simCJS to create a data set
  data <- simCJS(phi = phi[s], p = p[s], show.plot = FALSE)
  # Save the data in slice s of the 3D results array called ch
  ch[,,s] <- data$ch
}

# Look at simulated data (not shown)
ch                                  # unwieldy
ch[1:5,,]                           # First 5 individuals only
head(tmp <- apply(ch, c(1,3), sum)) # Number of captures (first 6 ind.)

op <- par(mfrow = c(3,3))
devAskNewPage(dev.interactive(orNone=TRUE)) # ~~~~ preferred to 'browser()'
for(k in 1:nspec){ # Capture-frequency for each species
  plot(table(tmp[,k]), main = paste('Capture frequencies (species', k, ')'), xlim = c(0,
  data$n.occ), ylab = 'Number of ind.', frame = FALSE)
  # browser() # ~~~~ take out for autotesting
}
par(op)
plot(table(tmp), main = 'Overall capture frequencies\n (all species)',
    xlim = c(0, data$n.occ), ylab = 'Number of ind.', frame = FALSE,
    type = 'h', lend = 'butt', lwd = 10)


# Bundle and summarize data set
str(bdata <- list(y = ch, f = data$f, n.ind = data$n.ind,
    n.occ = data$n.occ, nspec = dim(ch)[3]))
# List of 5
# $ y    : num [1:100, 1:6, 1:50] 1 1 1 1 1 1 1 1 1 1 ...
# $ f    : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
# $ n.ind: num 100
# $ n.occ: num 6
# $ nspec: int 50

# Specify model in BUGS language
cat(file = "cjs3.txt","
model {

  # Priors and hyperpriors
  # Define phi and p as random effects from a prior distribution
  # This is the submodel for how species vary in the community
  for(s in 1:nspec){
    phi[s] <- ilogit(lphi[s])
    lphi[s] ~ dnorm(mu.lphi, tau.lphi)
    p[s] <- ilogit(lp[s])
    lp[s] ~ dnorm(mu.lp, tau.lp)
  }

  # Give (hyper-)priors for the hyperparameters that
  # characterise the community
  mu.lphi <- logit(mean.phi)                 # Hyperpriors for survival hyperparams
  mean.phi ~ dunif(0,1)                      # mean hyperparam. (community average)
  tau.lphi <- pow(sigma.lphi, -2)
  sigma.lphi ~ dunif(0, 3)                   # sd hyperparam. (community heterogeneity)
  mu.lp <- logit(mean.p)                     # Hyperpriors for recapture hyperparams
  mean.p ~ dunif(0,1)                        # mean hyperparam.
  tau.lp <- pow(sigma.lp, -2)
  sigma.lp ~ dunif(0, 3)                     # sd hyperparam.

  # 'Likelihood'
  for(s in 1:nspec){                         # Loop over species
    for(i in 1:n.ind){                       # Loop over individuals
      # Define latent state at first capture
      z[i,f[i], s] <- 1
      for(t in (f[i]+1):n.occ){              # Loop over occasions
        # State process: the latent alive/dead state
        z[i,t,s] ~ dbern(z[i,t-1,s] * phi[s])# phi indexed by species now
        # Obs. process: relates true state to observed state, y [ ch
        y[i,t,s] ~ dbern(z[i,t,s] * p[s])    # p also indexed by species
      }
    }
  }
}
")

# Initial values
zst <- ch
for(s in 1:50){
  zst[,,s] <- zinit(ch[,,s])
}
inits <- function(){list(z = zst, mean.phi = runif(1), sigma.lphi = runif(1),
    mean.p = runif(1), sigma.lp = runif(1))}

# Parameters monitored
params<- c("mean.phi", "mu.lphi", "sigma.lphi", "mean.p", "mu.lp", "sigma.lp",
    "phi", "p") # could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~ testing

# Call JAGS (ART 12 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "cjs3.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out3)
print(out3, 3)
#              mean    sd   2.5%    50%  97.5% overlap0     f  Rhat n.eff
# mean.phi    0.486 0.023  0.439  0.486  0.530    FALSE 1.000 1.001  1401
# mu.lphi    -0.056 0.092 -0.246 -0.056  0.120     TRUE 0.721 1.001  1402
# sigma.lphi  0.441 0.083  0.286  0.437  0.614    FALSE 1.000 1.000  3000
# mean.p      0.287 0.020  0.246  0.287  0.325    FALSE 1.000 1.002  1015
# mu.lp      -0.914 0.100 -1.121 -0.909 -0.732    FALSE 1.000 1.002   972
# sigma.lp    0.423 0.088  0.255  0.421  0.601    FALSE 1.000 1.003  2642
# phi[1]      0.465 0.070  0.336  0.463  0.611    FALSE 1.000 1.000  2588
# phi[2]      0.470 0.064  0.350  0.468  0.604    FALSE 1.000 1.000  3000
# [ ... Output truncated ... ]
# p[49]       0.282 0.059  0.177  0.280  0.410    FALSE 1.000 1.000  3000
# p[50]       0.216 0.057  0.119  0.211  0.337    FALSE 1.000 1.000  3000

# ~~~~ plotting code for figs 3.3 and 3.4 ~~~~~~~~~
# Fig 3.3 Compare truth and estimates (marginal posteriors) for hyperparams
op <- par(mfrow = c(2, 2), las = 1)
plot(density(out3$sims.list$mu.lphi), col = "gray",
    main = "Community survival mean", xlab = "mu.lphi",
    frame = FALSE, axes = FALSE)
axis(1)  ;  axis(2)
abline(v = mu.lphi, lwd = 2)
plot(density(out3$sims.list$sigma.lphi), col = "gray",
    main = "Community survival heterogeneity", xlab = "sigma.lphi",
    frame = FALSE, axes = FALSE)
axis(1)  ;  axis(2)
abline(v = sigma.lphi, lwd = 2)
plot(density(out3$sims.list$mu.lp), col = "gray",
    main = "Community recapture mean", xlab = "mu.lp",
    frame = FALSE, axes = FALSE)
axis(1)  ;  axis(2)
abline(v = mu.lp, lwd = 2)
plot(density(out3$sims.list$sigma.lp), col = "gray",
    main = "Community recapture heterogeneity", xlab = "sigma.lp",
    frame = FALSE, axes = FALSE)
axis(1)  ;  axis(2)
abline(v = sigma.lp, lwd = 2)
par(op)

# Fig 3.4 Compare truth and estimates for species-level params
op <- par(mfrow = c(3, 4), las = 1)
oldask <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
for(s in 1:nspec){ # Print 6 species per page
  plot(density(out3$sims.list$phi[,s]), col = "gray",
      main = paste("Apparent survival species", s), xlab = "phi", xlim = c(0,1),
      ylab = '', axes = FALSE)
  axis(1)
  abline(v = phi[s], lwd = 2)
  plot(density(out3$sims.list$p[,s]), col = "gray",
      main = paste("Recapture species", s), xlab = "p", xlim = c(0,1),
      ylab = '', axes = FALSE)
  axis(1)
  abline(v = p[s], lwd = 2)
  # browser()  # ~~~ take out browser for testing
}
devAskNewPage(oldask)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 3.3.2 Community CJS models for improved estimates of species-level parameters
# -----------------------------------------------------------------------------

# Create "hyperdata" for the 50 species: species-level phi and p
set.seed(123)
nspec <- 50              # choose number of species
range.phi <- c(0.3, 0.7) # Draw survival from a uniform with these bounds
range.p <- c(0.1, 0.5)   # Draw recapture from another uniform

# Draw values of survival and recapture on the probability scale
phi <- runif(n = nspec, range.phi[1], range.phi[2]) # survival
p <- runif(n = nspec, range.p[1], range.p[2])       # recapture
round(sort(phi),2)       # Look at simulated values for phi (not shown)
round(sort(p),2)         # ... and for p (also not shown)

ch <- array(NA, dim = c(20, 6, nspec))

# Simulate data set for 50 species
for(s in 1:nspec){
  # Plug each pair (phi, p) into simCJS to create a data set
  data <- simCJS(n.marked = 4, phi = phi[s], p = p[s], show.plot = FALSE)
  # Save the data in slice i of the 3D results array called ch
  ch[,,s] <- data$ch
}

# Bundle and summarize data set
str(bdata <- list(y = ch, f = data$f, n.ind = data$n.ind,
    n.occ = data$n.occ, nspec = dim(ch)[3]))
# List of 5
# $ y    : num [1:20, 1:6, 1:50] 1 1 1 1 0 0 0 0 0 0 ...
# $ f    : int [1:20] 1 1 1 1 2 2 2 2 3 3 ...
# $ n.ind: num 20
# $ n.occ: num 6
# $ nspec: int 50

# Specify model in BUGS language
cat(file = "cjs4f.txt","
model {

# Priors
  # These are independent for each phi and p,
  # and there are no shared parameters that are estimated
  for(s in 1:nspec){
    phi[s] ~ dunif(0, 1)
    p[s] ~ dunif(0, 1)
  }

  # 'Likelihood'
  for(s in 1:nspec){                          # Loop over species
    for(i in 1:n.ind){                        # Loop over individuals
      # Define latent state at first capture
      z[i,f[i], s] <- 1
      for(t in (f[i]+1):n.occ){               # Loop over occasions
        # State process: the latent alive/dead state
        z[i,t,s] ~ dbern(z[i,t-1,s] * phi[s]) # phi indexed by species now
        # Obs. process: relates true state to observed state, y [ ch
        y[i,t,s] ~ dbern(z[i,t,s] * p[s])     # p also indexed by species
      }
    }
  }
}
")

# Initial values
zst <- ch
for(s in 1:50){
  zst[,,s] <- zinit(ch[,,s])
}
inits <- function(){list(z = zst, phi = runif(50), p = runif(50))}

# Parameters monitored
params <- c("phi", "p") # could add "z"

# MCMC settings
# na <- 1000 ; ni <- 20000 ; nt <- 10 ; nb <- 10000 ; nc <- 3
na <- 1000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3   # ~~~ testing

# Call JAGS (ART 3 min), check convergence and summarize posteriors
out4f <- jags(bdata, inits, params, "cjs4f.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out4f)
print(out4f, 3) # not shown

# Next, we fit the random-effects CJS model.
# Initial values
inits <- function(){list(z = zst, mean.phi = runif(1), sigma.lphi = runif(1),
    mean.p = runif(1), sigma.lp = runif(1))}

# Parameters monitored
params <- c("mean.phi", "mu.lphi", "sigma.lphi", "mean.p", "mu.lp",
    "sigma.lp", "phi", "p")

# MCMC settings
# na <- 1000 ; ni <- 20000 ; nt <- 10 ; nb <- 10000 ; nc <- 3
na <- 1000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ testing

# Call JAGS (ART 3 min), check convergence and summarize posteriors
out4r <- jags(bdata, inits, params, "cjs3.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out4r)
print(out4r, 3) # not shown

(avge.perc.error.fixed <- mean(100*abs(phi-out4f$mean$phi) / phi))
(avge.perc.error.random <- mean(100*abs(phi-out4r$mean$phi) / phi))
(RMSE.fixed <- sqrt(mean((phi-out4f$mean$phi)^2)))
(RMSE.random <- sqrt(mean((phi-out4r$mean$phi)^2)))
# [1] 23.796
# [1] 16.05593
# [1] 0.1572338
# [1] 0.09872242

# ~~~ extra code for figure 3.5 and 3.6 ~~~~~~~~~~~
# Compare truth and estimates (marginal posteriors) for all 50 species
op <- par(mfrow = c(2, 2), las = 1)

# Fixed-effects estimates
plot(phi, out4f$mean$phi, main = "Fixed-effects estimates", ylim = c(0,1),
    xlab = "True apparent survival (phi)", ylab = "", pch = 16, frame = FALSE)
segments(phi, out4f$q2.5$phi, phi, out4f$q97.5$phi)
abline(0,1, col = "red")
abline(lm(out4f$mean$phi ~ phi), col = 'blue', lty = 2)
points(phi, out4f$mean$phi, pch = 16)

plot(p, out4f$mean$p, main = "Fixed-effects estimates", ylim = c(0, 1),
    xlab = "True recapture (p)", ylab = "", pch = 16, frame = FALSE)
segments(p, out4f$q2.5$p, p, out4f$q97.5$p)
abline(0, 1, col = "red")
abline(lm(out4f$mean$p ~ p), col = 'blue', lty = 2)
points(p, out4f$mean$p, pch = 16)

# Random-effects estimates
plot(phi, out4r$mean$phi, main = "Random-effects estimates", ylim = c(0, 1), xlab = "True apparent survival (phi)", ylab = "", pch = 16, frame = FALSE)
segments(phi, out4r$q2.5$phi, phi, out4r$q97.5$phi)
abline(0, 1, col = "red")
abline(lm(out4r$mean$phi ~ phi), col = 'blue', lty = 2)
abline(h = out4r$mean$mean.phi, col = 'black', lty = 1)
points(phi, out4r$mean$phi, pch = 16)

plot(p, out4r$mean$p, main = "Random-effects estimates", ylim = c(0, 1), xlab = "True recapture (p)", ylab = "", pch = 16, frame = FALSE)
segments(p, out4r$q2.5$p, p, out4r$q97.5$p)
abline(0, 1, col = "red")
abline(lm(out4r$mean$p ~ p), col = 'blue', lty = 2)
abline(h = out4r$mean$mean.p, col = 'black',lty = 1)
points(p, out4r$mean$p, pch = 16)
par(op)

# Fig 3.6
op <- par(mfrow = c(3, 4), las = 1)
for(s in 1:6){
  plot(density(out4f$sims.list$phi[,s]), xlim = c(0,1), xlab = "phi", ylab = "",
      main = paste("Species", s, "(fixed effects)"), col = "blue", frame = FALSE)
  abline(v = phi[s], col = "red")
  plot(density(out4r$sims.list$phi[,s]), xlim = c(0,1), xlab = "phi", ylab = "",
      main = paste("Species", s, "(random effects)"), col = "gold", frame = FALSE)
  abline(v = phi[s], col = "red")
}
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 3.3.3 Community CJS models for exploration of patterns among species
# --------------------------------------------------------------------

# Create the species-level data
set.seed(24)
nspec <- 50                               # Number of species (need even number here)
beija <- c(rep('yes', nspec/2), rep('no',nspec/2)) # beija-flor indicator
mass <- round(runif(nspec, 10, 100))               # average mass of species 10-100 g

# Look at covariate data: beija-flor indicator and body mass
head(data.frame(beija, mass)) # First 6 species
#   beija mass
# 1   yes   36
# 2   yes   30
# 3   yes   73
# 4   yes   57
# 5   yes   70
# 6   yes   93

# express effect of mass for a centered version of original mass
massC <- mass - mean(mass)

# (1) Get design matrix (DM) and inspect how it looks like
(DM <- model.matrix(~ beija * massC-1-massC))
#    beijano beijayes beijano:mass beijayes:mass
# 1        0        1         0.00        -12.96
# 2        0        1         0.00        -18.96
# [ ... output truncated ... ]
# 49       1        0       -15.96          0.00
# 50       1        0        -8.96          0.00

# (2) Pick values for parameter vector (PV)
# Orders of DM and PV must match
mu.alpha.nonbeija <- 0 # logit-linear intercept beijaflor=no
mu.alpha.beija <- 0.7  # logit-linear intercept beijaflor=yes
beta.nonbeija <- 0.03  # logit-linear slope for beijaflor=no
beta.beija <- 0.01     # logit-linear slope for beijaflor=yes
(PV <- c(mu.alpha.nonbeija, mu.alpha.beija, beta.nonbeija, beta.beija))

# (3) Expected survival at link scale [ linear predictor (LinP)
(LinP <- DM %*% PV) # %*% does matrix multiplication
#       [,1]
# 1   0.5704
# 2   0.5104
# [ ... output truncated ... ]
# 49 -0.4788
# 50 -0.2688

# (4) Add random species effects and pick value of hyperparameter sigma_phi
sigma.phi <- 0.4
alpha0 <- rnorm(50, 0, sigma.phi)
lphi <- LinP + alpha0

# (5) Inverse-logit transformation yields survival for every species
exp.phi <- plogis(LinP) # Expected survival (w/o species random effects)
phi <- plogis(lphi)     # Realized survival (with species random effects)

# All the simulated data (not shown)
data.frame(beija, mass, 'expected survival' = round(exp.phi,3),
  'realized survival' = round(phi, 3))

# ~~~~ extra code for figure 3.7 ~~~~~~~~~~
# Plots of survival data alone and as function of tax. group and body mass
op <- par(mfrow = c(1, 2), las = 1)
hist(phi, col = 'grey', xlim = c(0,1), xlab = 'Realized apparent survival (phi)', main = '')
oo1 <- order(massC[beija == 'yes'])  # get order of mass for plotting
oo2 <- order(massC[beija == 'no'])   # get order of mass for plotting
plot(sort(massC[beija == 'yes']), exp.phi[beija == 'yes'][oo1], type = 'l',
    col = 'black', lty = 'solid', xlab = 'Body mass (in g, centered)',
    ylab = 'Apparent survival (phi)', ylim = c(0,1), frame = FALSE, xlim = c(-40, 50))
lines(sort(massC[beija == 'no']), exp.phi[beija == 'no'][oo2], col = 'black', lty = 2)
points(massC[beija == 'yes'], phi[beija == 'yes'], pch = 16, col = 'black')
points(massC[beija == 'no'], phi[beija == 'no'], pch = 1, col = 'black')
legend('bottomright', c('Beija-flores', 'Non-beija-flores'),
    pch = c(16,1), lty = c(1, 2), bty = 'n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Simulate recapture probabilities for 50 species
mu.lp <- -1     # mean recapture of all species on logit scale
sigma.lp <- 0.5 # standard deviation of recapture on logit scale
lp <- rnorm(n = nspec, mean = mu.lp, sd = sigma.lp) # logit(recapture)
p <- plogis(lp)
round(p, 3)     # Recapture prob. of the 50 species (not shown)

str(ch <- array(NA, dim = c(100, 6, 50)))
for(s in 1:nspec){
  data <- simCJS(phi = phi[s], p = p[s], show.plot = FALSE)
  ch[,,s] <- data$ch
}


# Bundle and summarize data set
new.beija <- rep(c(1, 2), each = 25) # 1 is non-beija, 2 is beija-flor
str(bdata <- list(y = ch, f = data$f, n.ind = data$n.ind, n.occ = data$n.occ,
  nspec = dim(ch)[3], beija = new.beija, massC = massC))
# List of 7
# $ y    : num [1:100, 1:6, 1:50] 1 1 1 1 1 1 1 1 1 1 ...
# $ f    : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
# $ n.ind: num 100
# $ n.occ: num 6
# $ nspec: int 50
# $ beija: num [1:50] 1 1 1 1 1 1 1 1 1 1 ...
# $ massC: num [1:50] -12.96 -18.96 24.04 8.04 21.04 ...

# Specify model in BUGS language
cat(file = "cjs5.txt","
model {

  # Priors and hyperpriors
  # Submodel for how species vary in the community
  # logit(phi) ~ beija * mass
  for(s in 1:nspec){ # Define priors for phi and p parameters
    phi[s] <- ilogit(lphi[s])
    lphi[s] ~ dnorm(mu.lphi[s], tau.lphi)
    mu.lphi[s] <- alpha.lphi[beija[s]] + beta.lphi[beija[s]] * massC[s]
    p[s] <- ilogit(lp[s])
    lp[s] ~ dnorm(mu.lp, tau.lp)
  }

  # Submodel for community hyperparameters
  for(k in 1:2){
    alpha.lphi[k] <- logit(mean.phi[k])
    mean.phi[k] ~ dunif(0,1)
    beta.lphi[k] ~ dnorm(0, 0.1)
  }
  tau.lphi <- pow(sigma.lphi, -2)
  sigma.lphi ~ dunif(0, 1)   # Community heterogeneity in survival
  mu.lp <- logit(mean.p)     # Hyperpriors for recapture hyperparams
  mean.p ~ dunif(0,1)        # Community average recapture
  tau.lp <- pow(sigma.lp, -2)
  sigma.lp ~ dunif(0, 1)     # Community heterogeneity in recapture

  # 'Likelihood'
  for(s in 1:nspec){                          # Loop over species
    for(i in 1:n.ind){                        # Loop over individuals
      # Define latent state at first capture
      z[i,f[i], s] <- 1
      for(t in (f[i]+1):n.occ){               # Loop over occasions
        # State process: the latent alive/dead state
        z[i,t,s] ~ dbern(z[i,t-1,s] * phi[s]) # phi indexed by species now
        # Obs. process: relates true state to observed state, y [ ch
        y[i,t,s] ~ dbern(z[i,t,s] * p[s])     # p also indexed by species
      }
    }
  }
}
")

# Initial values
zst <- ch
for(s in 1:50){
  zst[,,s] <- zinit(ch[,,s])
}
inits <- function(){list(z = zst, mean.phi = runif(2), sigma.lphi = runif(1),
    mean.p = runif(1), sigma.lp = runif(1))}

# Parameters monitored
params <- c('mean.phi', 'alpha.lphi', 'beta.lphi', 'sigma.lphi', 'mean.p', 'mu.lp',
    'sigma.lp', 'mu.lphi', 'phi', 'p')

# MCMC settings
# na <- 5000 ; ni <- 20000 ; nt <- 10 ; nb <- 10000 ; nc <- 3
na <- 5000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ testing

# Call JAGS (ART 14 min), check convergence and summarize posteriors
out5 <- jags(bdata, inits, params, "cjs5.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out5)
print(out5, 2)
#                mean   sd  2.5%   50% 97.5% overlap0    f Rhat n.eff
# mean.phi[1]    0.66 0.03  0.60  0.66  0.72    FALSE 1.00 1.00  1029
# mean.phi[2]    0.54 0.04  0.46  0.54  0.61    FALSE 1.00 1.01   285
# alpha.lphi[1]  0.67 0.14  0.40  0.67  0.94    FALSE 1.00 1.00   990
# alpha.lphi[2]  0.14 0.15 -0.15  0.14  0.43     TRUE 0.83 1.01   286
# beta.lphi[1]   0.01 0.01  0.00  0.01  0.02    FALSE 0.98 1.00  3000
# beta.lphi[2]   0.04 0.01  0.03  0.04  0.05    FALSE 1.00 1.00   398
# sigma.lphi     0.52 0.11  0.31  0.51  0.74    FALSE 1.00 1.01   181
# mean.p         0.25 0.02  0.21  0.25  0.28    FALSE 1.00 1.01   282
# mu.lp         -1.13 0.10 -1.32 -1.12 -0.95    FALSE 1.00 1.01   290
# sigma.lp       0.45 0.09  0.29  0.45  0.65    FALSE 1.00 1.01   162
# [ ... output truncated ... ]

# ~~~ code for figures 3.8 and 3.9 ~~~~~~~~~~
# Fig 3.8
oo1 <- order(massC[beija == 'yes'])
oo2 <- order(massC[beija == 'no'])
plot(sort(massC[beija == 'yes']), exp.phi[beija == 'yes'][oo1], type = 'l',
    col = 'black', lty = 'solid', xlab = 'Body mass (in g, centered)',
    ylab = 'Apparent survival (phi)', ylim = c(0, 1), frame = FALSE,
    xlim = c(-40, 50))
lines(sort(massC[beija == 'no']), exp.phi[beija == 'no'][oo2],
    col = 'black', lty = 2)

# Add realized values of phi
points(massC[beija=='yes'], phi[beija=='yes'], col='black', pch=16)
points(massC[beija=='no'], phi[beija=='no'], col='black', pch=1)

# Add estimates from the hierarchical CJS model in BUGS
lines(sort(massC[beija == "yes"]), plogis(out5$mean$mu.lphi[beija == "yes"])[oo1],
    col = "blue", lty = 1)
lines(sort(massC[beija == "no"]), plogis(out5$mean$mu.lphi[beija == "no"])[oo2],
    col = "blue", lty = 2)
legend('bottomright', legend = c('True expected (beija-flor)',
    'Estimated expected (beija-flor)', 'True expected (passerine)',
    'Estimated expected (passerine)'), lty = c(1, 1, 2, 2),
    col = c('black', 'blue', 'black', 'blue'), bty = 'n')

# Fig 3.9
off <- 0.5
plot(1:nspec, phi, xlab = "Species number", ylab = "phi", ylim = c(0,1),
    pch = 1,frame = F)

# Bayesian posterior means (estimated as random effects)
points(1:nspec+off, out5$mean$phi, pch = 16)
segments(1:nspec+off, out5$q2.5$phi, 1:nspec+off, out5$q97.5$phi)
legend('bottomleft', pch = c(1, 16), c('True', 'Estimated'), bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
