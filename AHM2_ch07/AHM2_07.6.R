#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7 : MODELING FALSE POSITIVES
# ====================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 6 mins

library(jagsUI)
library(AHMbook)

# 7.6 Modeling  false positives from acoustic monitoring data
# ===========================================================

# 7.6.1 A variation of the Chambert et al. (2017) model
# -----------------------------------------------------

# Simulation settings
set.seed(2019, kind = "Mersenne")
nsites <- 100     # Number of sites
nsurveys <- 5     # Number of replicates/occasions
psi <- 0.7        # Occupancy
p11 <- 0.5        # Detection probability at an occupied site
p10 <- 0.05       # False detection probability
lam <- 3          # Rate of true positives from ARU
ome <- 0.50       # Rate of false positives from ARU

# Simulate true occupancy states
z <- rbinom(nsites, 1, psi)

# Define detection probability
p <- z * p11 + (1-z) * p10

# Simulate occupancy data and ARU count frequencies
yARU <- y <- K <- Q <- matrix(NA, nsites, nsurveys)
for(i in 1:nsites){
  y[i,] <- rbinom(nsurveys, 1, p[i]) # Detection/nondetection data
  K[i,] <- rpois(nsurveys, lam*z[i]) # True positive detection frequency
  Q[i,] <- rpois(nsurveys, ome)      # False-positive detection frequency
  yARU[i,] <- K[i,] + Q[i,]          # Number of ARU detections
}

# Bundle and summarize data
str( bdata <- list(y = y, yARU = yARU, nsites = nsites, nsurveys = nsurveys ))
# List of 4
# $ y       : int [1:100, 1:5] 0 0 0 0 1 1 0 0 1 0 ...
# $ yARU    : int [1:100, 1:5] 1 2 5 5 2 2 0 4 7 4 ...
# $ nsites  : num 100
# $ nsurveys: num 5

# Specify Model A in BUGS language
cat(file = "modelA.txt","
model {

  # Priors
  psi ~ dunif(0, 1) # psi = Pr(Occupancy)
  p10 ~ dunif(0, 1) # p10 = Pr(y = 1 | z = 0)
  p11 ~ dunif(0, 1) # p11 = Pr(y = 1 | z = 1)
  lam ~ dunif(0, 1000)
  ome ~ dunif(0, 1000)

  # Likelihood:process and observation models
  for (i in 1:nsites) {
    z[i] ~ dbern(psi) # Occupancy status of site i
    p[i] <- z[i] * p11 + (1-z[i]) * p10 # false-positive detection model
    for(j in 1:nsurveys) {
      y[i,j] ~ dbern(p[i]) # Binary occupancy data
      yARU[i,j] ~ dpois(lam * z[i] + ome) # ARU detection frequency data
    }
  }
}
")

# Initial values
inits <- function(){list(z = apply(y, 1, max), psi = runif(1),
    p10 = runif(1, 0, 0.05), p11 = runif(1, 0.5, 0.8), lam = runif(1, 1, 2),
    ome = runif(1, 0, 0.4) )}

# Parameters monitored
params <- c("psi", "p10", "p11", "lam", "ome")

# MCMC settings
na <- 1000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3

# Call JAGS (tiny ART), gauge convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "modelA.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out1, layout=c(2,3))
print(out1, 3)
#      mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# psi 0.706 0.044 0.616 0.706 0.789    FALSE 1 1.001  1662
# p10 0.069 0.021 0.034 0.067 0.115    FALSE 1 1.001  2027
# p11 0.476 0.025 0.428 0.477 0.526    FALSE 1 1.000  3000
# lam 3.077 0.112 2.863 3.075 3.296    FALSE 1 1.005   540
# ome 0.353 0.050 0.264 0.350 0.457    FALSE 1 1.002  1238


# 7.6.2 Making use of acoustic validation data
# --------------------------------------------

k <- n <- matrix(NA, nrow = nsites, ncol = nsurveys)
for(j in 1:nsurveys){
  bleen <- valid_data(yARU[,j], tp = K[,j], n.valid = 50)
  n[,j] <- bleen$n
  k[,j] <- bleen$k
}

# Bundle and summarize data
str(bdata <- list(y = y, yARU = yARU, nsites = nsites, nsurveys = nsurveys,
    n = n, k = k))

# Specify Model B in BUGS language
cat(file = "modelB.txt","
model {

  # Priors
  psi ~ dunif(0, 1) # psi = Pr(Occupancy)
  p10 ~ dunif(0, 1) # p10 = Pr(y = 1 | z = 0)
  p11 ~ dunif(0, 1) # p11 = Pr(y = 1 | z = 1)
  lam ~ dunif(0, 1000)
  ome ~ dunif(0, 1000)

  # Likelihood and process model
  for (i in 1:nsites){
    z[i] ~ dbern(psi) # Occupancy status of site i
    p[i] <- z[i] * p11 + (1-z[i]) * p10 # False-positive detection model
    true.positives[i] <- lam*z[i] / (lam * z[i] + ome) # Pr(true positive)
    for(j in 1:nsurveys){
      y[i,j] ~ dbern(p[i]) # Binary occupancy data
      yARU[i,j] ~ dpois(lam*z[i] + ome) # ARU detection freq. data
      K[i,j] ~ dbin(true.positives[i], yARU[i,j]) # K = true positives
      Q[i,j] <- yARU[i,j] - K[i,j] # False positives, derived
      k[i,j] ~ dhyper(K[i,j], Q[i,j], n[i,j], 1)
    }
  }
}
")

# Initial values, must satisfy some constraints
zst <- apply(y, 1, max); zst[apply(k,1,sum)>0] <- 1
Kst <- k
inits <- function(){list(z = zst, psi = runif(1), p10 = runif(1, 0, 0.05),
    p11 = runif(1, 0.5, 0.8), lam = runif(1, 1, 2), ome = runif(1, 0, 0.4),
    K = Kst)}

# Parameters monitored
params <- c("psi", "p10", "p11", "lam", "ome")

# MCMC settings. Mixing is worse so we run longer.
na <- 5000 ; ni <- 33000 ; nt <- 2 ; nb <- 3000 ; nc <- 3

# Call JAGS (ART 2 min), gauge convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "modelB.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out2, layout=c(2,3))
print(out2, 3)
#      mean    sd  2.5%   50% 97.5% overlap0 f Rhat n.eff
# psi 0.706 0.045 0.614 0.707 0.790    FALSE 1    1 38245
# p10 0.068 0.021 0.033 0.066 0.114    FALSE 1    1 11329
# p11 0.476 0.026 0.425 0.476 0.528    FALSE 1    1 45000
# lam 3.020 0.102 2.819 3.019 3.221    FALSE 1    1  8432
# ome 0.392 0.043 0.313 0.391 0.480    FALSE 1    1  5412


# 7.6.3. The coupled classification model: A new model for integrating
#        species classification from bioacoustics data with
#        occupancy models
# ---------------------------------------------------------------------

# Simulation settings
set.seed(2019)
nsites <- 100          # Number of sites
nsurveys <- 5          # Number of replicates
psi <- 0.7             # Occupancy
p11 <- 0.5             # Detection probability at an occupied site
p10 <- 0.05            # False-positive rate for occupancy surveys
lam <- 1               # Rate at which target species vocalizations are obtained
ome <- 3               # Rate at which non-target species are obtained
mu.mix.target <- 1     # mean of feature for target species
mu.mix.nontarget <- -1 # ... for non-target species
sd.mix <- 0.5          # sd of feature for both

# Simulate true occupancy states
z <- rbinom(nsites, 1, psi)

# Define detection probability, with false positives
p <- z * p11 + (1-z) * p10

# Simulate occupancy data and ARU count frequencies
# Main simulation loop
covdata <- NULL
yARU <- y <- K <- Q <- matrix(NA, nsites, nsurveys)
for(i in 1:nsites){
  for(j in 1:nsurveys){
    y[i,j] <- rbinom(1, 1 ,p[i]) # Detection/nondetection data
    K[i,j] <- rpois(1, lam*z[i]) # True detection frequency
    Q[i,j] <- rpois(1, ome)      # False-positive detection frequency
    yARU[i,j] <- K[i,j] + Q[i,j] # Number of ARU detections
    if(yARU[i,j] > 0){ # now simulate some scores and store them for output
      K.score <- rnorm(K[i,j], mu.mix.target, sd.mix)
      Q.score <- rnorm(Q[i,j], mu.mix.nontarget, sd.mix)
      covdata <- rbind(covdata, cbind(rep(i, yARU[i,j]), rep(j, yARU[i,j]),
      c(K.score, Q.score), c(rep(1, K[i,j]), rep(2, Q[i,j]) ) ))
    }
  }
}
colnames(covdata) <- c("siteID", "occID", "score", "truth")

# Harvest the data
siteid <- covdata[,1]
occid <- covdata[,2]
score <- covdata[,3]
truth <- covdata[,4]

# Bundle and summarize data
str(bdata <- list(y = y, score = score, yARU = yARU, siteid = siteid,
    occid = occid, nsamples = length(score), nsites = nsites,
    nsurveys = nsurveys))
# List of 8
# $ y       : int [1:100, 1:5] 0 0 0 1 1 0 0 1 0 0 ...
# $ score   : num [1:1822] -1.12 -1.17 -1.38 -1.43 -1.39 ...
# $ yARU    : int [1:100, 1:5] 4 2 5 7 8 2 1 1 7 6 ...
# $ siteid  : num [1:1822] 1 1 1 1 1 1 1 1 1 1 ...
# $ occid   : num [1:1822] 1 1 1 1 2 2 2 2 2 2 ...
# $ nsamples: int 1822
# $ nsites  : num 100
# $ nsurveys: num 5

# Specify Model C in BUGS language
# Occupancy with false-positives, bioacoustics data with built-in
# species classification using Gaussian mixtures
cat(file="modelC.txt","
model {

  # Priors
  psi ~ dunif(0, 1) # psi = Pr(Occupancy)
  p10 ~ dunif(0, 1) # p10 = Pr(y = 1 | z = 0)
  p11 ~ dunif(0, 1) # p11 = Pr(y = 1 | z = 1)
  lam ~ dunif(0, 1000) # lambda: rate of target-species calls detected
  ome ~ dunif(0, 1000) # omega: rate of non-target detections

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(0, 0.01)
  mu[2] ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 10)
  tau <- 1 / (sigma * sigma)

  # Likelihood part 1: detection data and ARU counts
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # Latent occupancy states
    p[i] <- z[i]*p11 + (1-z[i])*p10 # Detection probability
    site.prob[i] <- lam*z[i]/(lam*z[i]+ome) # Pr(sample is target species)
    for(j in 1:nsurveys) { # Loop over occasions
      y[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
      yARU[i,j] ~ dpois(lam*z[i] + ome) # Total samples processed
    }
  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    # Sample specific covariate
    score[k] ~ dnorm(mu[g[k]], tau) # parameters are group specific
    probs[k,1] <- site.prob[siteid[k]]
    probs[k,2] <- 1 - site.prob[siteid[k]] # the prior class probabilities
    g[k] ~ dcat(probs[k,])
    N1[k] <- ifelse(g[k]==1, 1, 0)
  }

  # Derived quantities
  Npos <- sum(N1[])
}
")

# Initial values
zst <- rep(1, nrow(y))
gst <- sample(1:2, length(score), replace = TRUE)
gst[score>0] <- 1
gst[score<=0] <- 2
inits <- function(){ list(mu = c(1, -1), sigma = 0.2, z = zst,
    psi = runif(1), p10 = runif(1, 0, 0.05), p11 = runif(1, 0.5, 0.8),
    lam = runif(1, 1, 2), ome = runif(1, 0, 0.4), g = gst) }

# Parameters monitored
params <- c("psi", "p10", "p11", "lam", "ome", "mu", "sigma", "Npos")

# MCMC settings
na <- 1000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3

# Call JAGS (ART 1 min), gauge convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "modelC.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out3, layout=c(2,3))
print(out3, 3)
#          mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# psi     0.704 0.045   0.612   0.705   0.788    FALSE 1 1.003   570
# p10     0.055 0.019   0.023   0.053   0.099    FALSE 1 1.002  1673
# p11     0.470 0.026   0.420   0.470   0.521    FALSE 1 1.005   476
# lam     0.927 0.053   0.821   0.926   1.034    FALSE 1 1.004   432
# ome     2.991 0.078   2.840   2.991   3.145    FALSE 1 1.001  2817
# mu[1]   1.006 0.029   0.948   1.006   1.064    FALSE 1 1.000  2636
# mu[2]  -0.987 0.013  -1.012  -0.988  -0.961    FALSE 1 1.002   904
# sigma   0.480 0.008   0.464   0.480   0.498    FALSE 1 1.003   619
# Npos  327.099 5.284 317.000 327.000 338.000    FALSE 1 1.001  1286

table(truth)
#   1    2
# 334 1488


# 7.6.4. Predicting the group membership (species) of a given sample
# ------------------------------------------------------------------

# Parameters monitored
params <- c("psi", "p10", "p11", "lam", "ome", "mu", "sigma", "g")

# Refit model 3 (ART 1 min)
out3g <- jags(bdata, inits, params, "modelC.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Process the output to compute Pr(g=1) for each sample
g <- out3g$sims.list$g # Extract posterior samples for p
gbar <- out3g$mean$g
prob <- apply(g, 2, function(x) mean(x == 1))
plot(score, prob, ylab = "Posterior class probability", xlab = "Feature score",
    frame = FALSE, pch = 20, cex=2.5) # Produces Fig. 7.4 without the red bits;

# ~~~~ inserted code for the red bits ~~~~~~~~
look <- (score > -0.20 & prob < 0.046) # trial and error is necessary
points(score[look], prob[look], pch=20, cex=2.5, col="red")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 7.6.5  The uncoupled classification model
# -------------------------------------------

# Model C v2, uncoupled
# ~~~~ inserted from previous model C ~~~~~~~~~~
cat(file="modelCv2.txt", "
model{

  # Priors
  psi ~ dunif(0,1)	# psi = Pr(Occupancy)
  p10 ~ dunif(0,1)	# p10 = Pr(y = 1 | z = 0)
  p11 ~ dunif(0,1)	# p11 = Pr(y = 1 | z = 1)
  lam ~ dunif(0,1000)	# lambda: rate of target-species calls detected
  ome ~ dunif(0,1000)	# omega: rate of non-target detections
  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(0, 0.01)
  mu[2] ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 10)
  tau <- 1/(sigma*sigma)

  # Likelihood
  for (i in 1:nsites){
    z[i] ~ dbern(psi)			# Latent occupancy states
    p[i] <- z[i]*p11 + (1-z[i])*p10 	# Detection probability
    site.prob[i]<- lam*z[i]/(lam*z[i] + ome)   # Pr(a given call is target species)
    for(j in 1:nsurveys){
      y[i,j] ~ dbern(p[i])	      # Observed occupancy data (if available)
      yARU[i,j] ~ dpois(lam*z[i] + ome)  # Total samples processed
    } # j
  }
  # ~~~~ end of inserted  code ~~~~~~~~~~
  pi0 ~ dunif(0,1) # New class probability parameter
  for(s in 1:nsamples) {
    # Sample specific covariate
    score[s] ~ dnorm(mu[g[s]], tau) # Mixture distribution for score
    probs[s,1] <- pi0
    probs[s,2] <- 1-pi0
    g[s] ~ dcat(probs[s,])
    N1[s] <- ifelse(g[s]==1,1,0)
  }
  Npos <- sum(N1[])
}
")

# ~~~~ code inserted from MS ~~~~~~~~
# Initial values
zst <- apply(y, 1, max)
zst <- rep(1, nrow(y))
gst <- sample(1:2, length(score), replace=TRUE)
gst[score>0] <- 1
gst[score<=0] <- 2
inits <- function(){list(  mu=c(1, -1), sigma = 0.2, z = zst,
    psi=runif(1), p10=runif(1, 0, 0.05), p11=runif(1, 0.5, 0.8),
    lam=runif(1,1,2),ome=runif(1, 0, 0.4), g = gst ) }

# Parameters monitored
params <- c("psi", "p10","p11","lam","ome" , "sigma", "mu","Npos","pi0")

# MCMC settings
ni <- 2000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3; na<- 1000

# Call JAGS and summarize posteriors
library(jagsUI)
fit3b <- jags(bdata, inits, params, "modelCv2.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE)
print(fit3b, dig = 3)
# ~~~~ end of inserted code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# psi     0.688 0.072   0.528   0.694   0.814    FALSE 1 1.009   295
# p10     0.060 0.033   0.010   0.055   0.138    FALSE 1 1.010   817
# p11     0.476 0.036   0.410   0.474   0.551    FALSE 1 1.005   422
# lam     0.962 0.216   0.533   0.967   1.371    FALSE 1 1.014   150
# ome     2.979 0.183   2.617   2.978   3.329    FALSE 1 1.017   129
# sigma   0.479 0.008   0.462   0.478   0.496    FALSE 1 1.001  1892
# mu[1]   1.009 0.031   0.950   1.009   1.071    FALSE 1 1.002   981
# mu[2]  -0.988 0.013  -1.013  -0.988  -0.962    FALSE 1 1.001  1837
# Npos  326.789 5.703 316.000 327.000 338.000    FALSE 1 1.003   837
# pi0     0.180 0.010   0.161   0.179   0.199    FALSE 1 1.001  2023


# 7.6.6  A simulation study to evaluate the effect of occupancy on classification
# -------------------------------------------------------------------------------

# anova(out <- lm(aucDIFF ~ beta0 + p0 + beta1 + sigma + beta0*p0 +
#     beta0*beta1 + beta0*sigma + p0*beta1 + p0*sigma + beta1*sigma,
#     data = data))

# Analysis of Variance Table

# Response: aucDIFF
#              Df   Sum Sq   Mean Sq    F value     Pr(>F)
# beta0         4   118.76    29.689   576.9596  < 2.2e-16 ***
# p0            1     1.07     1.072    20.8306  1.077e-05 ***
# beta1         2     1.66     0.832    16.1599  4.752e-07 ***
# sigma         2   613.46   306.728  5960.8267  < 2.2e-16 ***
# beta0:p0      4     0.38     0.096     1.8579     0.1211
# beta0:beta1   8     7.04     0.879    17.0909  < 2.2e-16 ***
# beta0:sigma   8     2.09     0.261     5.0747  1.458e-05 ***
# p0:beta1      2     0.02     0.009     0.1845     0.8317
# p0:sigma      2     0.15     0.077     1.4960     0.2275
# beta1:sigma   4     0.35     0.088     1.7130     0.1504
# Residuals   142     7.31     0.051

# 7.6.7 Validation and digital archive records
# --------------------------------------------

# Choose some validation sample
nvalid <- 100
valid.id <- sample(1:length(score), nvalid)

# Group membership (target species = 1, false-positive = 2).
g <- rep(NA, length(score))
g[valid.id] <- truth[valid.id] # Set some to 1 or 2

# Bundle and summarize data. Note g is data now.
str(bdata <- list(y = y, score = score, yARU = yARU, siteid = siteid,
    occid = occid, nsamples = length(score), nsites = nsites,
    nsurveys = nsurveys, g = g) ) # not shown

# ~~~~ code inserted from MS ~~~~~~~~~~~~~~~~~~~
# Model D, built-in species classification
cat(file = "modelD.txt","
model {

  # Priors
  psi ~ dunif(0, 1)     # psi = Pr(Occupancy)
  p10 ~ dunif(0, 1)     # p10 = Pr(y = 1 | z = 0)
  p11 ~ dunif(0, 1)     # p11 = Pr(y = 1 | z = 1)
  lam ~ dunif(0, 1000)  # lambda: rate of target-species calls
  ome ~ dunif(0, 1000)  # omega: rate of non-target calls

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(0, 0.01)
  mu[2] ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 10)
  tau <- 1 / (sigma * sigma)

  # Likelihoodpart 1
  # For basic occupancy and ARU data sets
  for (i in 1:nsites){    # Loop over sites
    z[i] ~ dbern(psi)                     # Latent occupancy states
    p[i] <- z[i] * p11 + (1 - z[i]) * p10 # Detection probability
    site.prob[i] <- lam * z[i] / (lam * z[i] + ome) # Pr(a given call is target species)
    for(j in 1:nsurveys){   # Loop over occasions
      y[i,j] ~ dbern(p[i])      # Observed occupancy data (if available)
      yARU[i,j] ~ dpois(lam * z[i] + ome)  # Total samples processed
    }
  }

  # Likelihood part 2
  # Mixture model for group membership and score distributions
  for(s in 1:nsamples){
    # Sample specific covariate model
    score[s] ~  dnorm(mu[g[s]], tau)
    probs[s,1] <- site.prob[siteid[s]]
    probs[s,2] <- 1-site.prob[siteid[s]]
    g[s] ~ dcat(probs[s,])
    N1[s] <- ifelse(g[s]==1,1,0)
  }
  # Derived quantities
  Npos <- sum(N1[])
} # end 'model'
")

# Initial values
zst <- apply(y, 1, max)
zst <- rep(1, nrow(y))
gst <- sample(1:2, length(score), replace=TRUE)
gst[score>0] <- 1
gst[score<=0] <- 2
gst[valid.id] <- NA  # validation samples must not be given inits !
inits <- function(){list(mu = c(1, -1), sigma = 0.2, z = zst,
  psi = runif(1), p10 = runif(1, 0, 0.05), p11 = runif(1, 0.5, 0.8),
  lam = runif(1, 1, 2), ome = runif(1, 0, 0.4))} # removed g from inits, no reason

# Parameters monitored
params <- c("psi", "p10", "p11", "lam", "ome", "mu","sigma", "Npos")

# MCMC settings, more iters because poor mixing
na <- 1000  ;  ni <- 12000  ;  nt <- 2  ;  nb <- 2000  ;  nc <- 3
# ~~~ end of inserted code ~~~~~~~~~~~

# Call JAGS (ART 4 min), gauge convergence and summarize posteriors
(out4 <- jags(bdata, inits, params, "modelD.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE) )
#          mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# psi     0.702 0.046   0.610   0.704   0.788    FALSE 1 1.001  3164
# p10     0.056 0.020   0.023   0.054   0.100    FALSE 1 1.000  9058
# p11     0.469 0.027   0.417   0.469   0.521    FALSE 1 1.000 15000
# lam     0.930 0.054   0.825   0.929   1.040    FALSE 1 1.001  3041
# ome     2.991 0.077   2.841   2.990   3.144    FALSE 1 1.000 15000
# sigma   0.480 0.009   0.463   0.480   0.497    FALSE 1 1.000 12148
# mu[1]   1.006 0.030   0.947   1.006   1.064    FALSE 1 1.000 15000
# mu[2]  -0.988 0.013  -1.013  -0.988  -0.962    FALSE 1 1.000 15000
# Npos  327.404 5.198 317.000 327.000 338.000    FALSE 1 1.000 15000

round( (unlist(out3$sd)/unlist(out4$sd) ), 3)
#   psi   p10   p11   lam   ome sigma   mu1   mu2  Npos
# 0.984 0.996 0.964 0.974 1.011 0.993 0.973 1.010 1.017



