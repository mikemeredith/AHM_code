#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5 : MODELING METACOMMUNITY DYNAMICS USING DYNAMIC COMMUNITY MODELS
# ==========================================================================
# Code from proofs dated 2020-08-19

# Approximate run time for this script: 20 mins
# With the full number of iterations: 52 hrs

library(AHMbook)
library(jagsUI)

# 5.5 Estimating species richness via the random-effects DCM with data augmentation
# =================================================================================

set.seed(1)
dat <- simDCM(nspec = 100, nsites = 50, nsurveys = 2, nyears = 6,
    mean.psi1 = 0.1, sig.lpsi1 = 3,
    range.mean.phi = c(0.6, 0.6), sig.lphi = 1,
    range.mean.gamma = c(0.1, 0.1), sig.lgamma = 1,
    range.mean.p = c(0.1, 0.1), sig.lp = 3)

# ** Number of species ever occurring: 100
# ** Number of species ever detected: 81
# ** Average number of years of occurrence: 5.59
# ** Average number of years with detection: 3.53

(missed <- which(dat$nyears.det == 0))       # species detected in 0 years
# Spec2 Spec9 Spec13 Spec20 Spec24 Spec35 Spec38 Spec40 Spec43 Spec48
#     2     9     13     20     24     35     38     40     43     48
# Spec50 Spec51 Spec72 Spec74 Spec81 Spec83 Spec88 Spec89 Spec94
#     50     51     72     74     81     83     88     89     94

# Toss out species never detected
y <- dat$y                   # copy 4D array
y <- y[,,,-missed]           # Drop data from 19 missed species
str(y)                       # 50 sites x 2 reps x 6 years x 81 species

# Augment the data set with 100 potential species
# Create zero 4D array for M species
nz <- 100                              # Number of 'zero species'
M <- dim(y)[4] + nz                    # Size of augmented data set
yaug <- array(0, dim = c(50, 2, 6, M)) # Prefill with zeroes
dim(yaug)                              # Check if it went well: and it did !

# Fill in the observed data into this larger array
yaug[,,,1:dim(y)[4]] <- y
str(yaug)
sum(y) ; sum(yaug)                     # Quick sum check: should give same sum

# Bundle and summarize data set (M is nspec for augmented data set)
str(bdata <- list(yaug = yaug, nsite = dim(yaug)[1], nsurvey = dim(yaug)[2],
    nyear = dim(yaug)[3], M = dim(yaug)[4]))
# List of 5
# $ yaug   : num [1:50, 1:2, 1:6, 1:181] 0 0 0 0 0 0 0 0 0 0 ...
# $ nsite  : int 50
# $ nsurvey: int 2
# $ nyear  : int 6
# $ M      : int 181

# ~~~~~~~~~~~ model code inserted from MS dated 2019-01-01 ~~~~~~~~~~

# First variant, "natural" parameterization
# '''''''''''''''''''''''''''''''''''''''''
# Specify model in BUGS language
cat(file = "DCM2A.txt", "
model {

  # Specify hyperpriors: Priors for the hyperparameters
  mu.lpsi1 <- logit(mean.psi1)    # Initial occupancy
  mean.psi1 ~ dunif(0, 1)
  tau.lpsi1 <- pow(sd.lpsi1, -2)
  sd.lpsi1 ~ dunif(0, 10)
  mu.lphi <- logit(mean.phi)      # Persistence
  mean.phi ~ dunif(0, 1)
  tau.lphi <- pow(sd.lphi, -2)
  sd.lphi ~ dunif(0, 10)
  mu.lgamma <- logit(mean.gamma)  # Colonization
  mean.gamma ~ dunif(0, 1)
  tau.lgamma <- pow(sd.lgamma, -2)
  sd.lgamma ~ dunif(0, 10)
  mu.lp <- logit(mean.p)          # Detection
  mean.p ~ dunif(0, 1)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 10)

  # Specify priors: Declare species-level effects as random
  for(k in 1:M){  # Loop over species
    logit(psi1[k]) <- lpsi1[k]    # Initial occupancy
    lpsi1[k] ~ dnorm(mu.lpsi1, tau.lpsi1)
    logit(phi[k]) <- lphi[k]      # Persistence
    lphi[k] ~ dnorm(mu.lphi, tau.lphi)
    logit(gamma[k]) <- lgamma[k]  # Colonization
    lgamma[k] ~ dnorm(mu.lgamma, tau.lgamma)
    logit(p[k]) <- lp[k]          # Detection
    lp[k] ~ dnorm(mu.lp, tau.lp)
  }

  # Data augmentation model
  omega ~ dunif(0,1)  # Royle et al.(JCGS, 2007) prior for data augmentation parameter
  # omega ~ dbeta(0.001, 1)   # Scale prior of Link (Ecology, 2013)
  for(k in 1:M){      # Loop over all M species, including all-zero species
    w[k] ~ dbern(omega)
  }

  # Ecological submodel: Define state conditional on parameters
  # Note every expression for z has premultiplication with w
  for (i in 1:nsite){                    # Loop over sites
    for(k in 1:M){                       # Loop over species
      # Initial conditions of system
      z[i,1, k] ~ dbern(w[k] * psi1[k])  # P/A at start of study
      # State transitions
      for (t in 2:nyear){                # Loop over years
        z[i,t,k] ~ dbern(w[k] * (z[i,t-1,k] * phi[k] + (1- z[i,t-1, k]) * gamma[k]) )
      }
    }
  }

  # Observation submodel, for augmented data set
  for (i in 1:nsite){                    # Loop over sites
    for(k in 1:M){                       # Loop over all M species
      for (j in 1:nsurvey){              # Loop over surveys
        for (t in 1:nyear){              # Loop over years
          yaug[i,j,t,k] ~ dbern(z[i,t,k] * p[k])
        }
      }
    }
  }

  # Derived parameters: Number of occupied sites and population occupancy
  for(k in 1:M){
    n.occ[1, k] <- sum(z[,1,k])          # Number of occupied sites
    psi[1, k] <- psi1[k]                 # Population occupancy
    for (t in 2:nyear){# Loop over years
      n.occ[t, k] <- sum(z[,t,k])
      psi[t, k] <- psi[t-1, k] * phi[k] + (1-psi[t-1, k]) * gamma[k]
    }
  }

  # Species richness: total and for each site individually
  Ntotal <- sum(w[])                     # Estimate of community size
  for(i in 1:nsite){
    for(t in 1:nyear){
      Nspec[i,t] <- sum(z[i,t,])         # Site-specific species richness
    }
  }
}
")

# Initial values (needed for both w and z)
zst <- apply(yaug, c(1,3,4), max) # Observed occurrence as inits for z
wst <- apply(zst, 3, max)         # Need inits for w also
inits <- function(){ list(w = wst, z = zst)}

# Parameters monitored
params <- c("omega", "Ntotal", "mu.lpsi1", "sd.lpsi1", "mu.lphi", "sd.lphi",
    "mu.lgamma", "sd.lgamma", "mu.lp", "sd.lp", "psi1", "phi", "gamma", "p",
    "n.occ", "Nspec")

# MCMC settings
# na <- 1000 ; ni <- 50000 ; nt <- 25 ; nb <- 25000 ; nc <- 20 # 4.8 hrs
# na <- 1000 ; ni <- 350000 ; nt <- 25 ; nb <- 25000 ; nc <- 3 # ~~~ keep to 3 cores, do more iterations, 26 hrs
na <- 1000 ; ni <- 500 ; nt <- 25 ; nb <- 250 ; nc <- 3 # ~~~ for testing, 9 mins

# Call JAGS from R, check convergence and summarize posteriors
out2A <- jags(bdata, inits, params, "DCM2A.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out2A)
summary(out2A, 3)
jags.View(out2A)

# Second variant, Tobler et al. (2015) and Yamaura et al. (2016) parameterization
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Specify model in BUGS language
cat(file = "DCM2B.txt", "
model {

  # Specify hyperpriors: Priors for the hyperparameters
  mu.lpsi1 <- logit(mean.psi1)    # Initial occupancy
  mean.psi1 ~ dunif(0, 1)
  tau.lpsi1 <- pow(sd.lpsi1, -2)
  sd.lpsi1 ~ dunif(0, 10)
  mu.lphi <- logit(mean.phi)      # Persistence
  mean.phi ~ dunif(0, 1)
  tau.lphi <- pow(sd.lphi, -2)
  sd.lphi ~ dunif(0, 10)
  mu.lgamma <- logit(mean.gamma)  # Colonization
  mean.gamma ~ dunif(0, 1)
  tau.lgamma <- pow(sd.lgamma, -2)
  sd.lgamma ~ dunif(0, 10)
  mu.lp <- logit(mean.p)          # Detection
  mean.p ~ dunif(0, 1)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 10)

  # Specify priors: Declare species-level effects as random
  for(k in 1:M){# Loop over species
    logit(psi1[k]) <- lpsi1[k]    # Initial occupancy
    lpsi1[k] ~ dnorm(mu.lpsi1, tau.lpsi1)
    logit(phi[k]) <- lphi[k]      # Persistence
    lphi[k] ~ dnorm(mu.lphi, tau.lphi)
    logit(gamma[k]) <- lgamma[k]  # Colonization
    lgamma[k] ~ dnorm(mu.lgamma, tau.lgamma)
    logit(p[k]) <- lp[k]          # Detection
    lp[k] ~ dnorm(mu.lp, tau.lp)
  }

  # Data augmentation model
  omega ~ dunif(0,1)  # Royle et al. (JCGS, 2007) prior for data augmentation parameter
  # omega ~ dbeta(0.001, 1)   # Scale prior of Link (Ecology, 2013)
  for(k in 1:M){      # Loop over all M species, including all-zero species
    w[k] ~ dbern(omega)
  }

  # Ecological submodel: Define state conditional on parameters
  # No premultiplication with w
  for (i in 1:nsite){             # Loop over sites
    for(k in 1:M){                # Loop over species
      # Initial conditions of system
      z[i,1, k] ~ dbern(psi1[k])
      # State transitions
      for (t in 2:nyear){         # Loop over years
        z[i,t,k] ~ dbern(z[i,t-1,k] * phi[k] + (1- z[i,t-1, k]) * gamma[k])
      }
    }
  }

  # Observation model, for augmented data set !
  # Note premultiplication with w comes only here
  for (i in 1:nsite){             # Loop over sites
    for(k in 1:M){                # Loop over all M species
      for (j in 1:nsurvey){       # Loop over surveys
        for (t in 1:nyear){       # Loop over years
          yaug[i,j,t,k] ~ dbern(w[k] * z[i,t,k] * p[k])
        }
      }
    }
  }

  # Derived parameters: Number of occupied sites and population occupancy
  # Note multiplication of z with w here
  for(k in 1:M){
    n.occ[1, k] <- sum(w[k] * z[,1,k])  # Number of occupied sites
    psi[1, k] <- psi1[k]          # Population occupancy
    for (t in 2:nyear){# Loop over years
      n.occ[t, k] <- sum(w[k] * z[,t,k])
      psi[t, k] <- psi[t-1, k] * phi[k] + (1-psi[t-1, k]) * gamma[k]
    }
  }

  # Species richness: total and for each site individually
  Ntotal <- sum(w[])              # Estimate of community size
  for(i in 1:nsite){
    for(t in 1:nyear){
      for(k in 1:M){
        tmp[i,t,k] <- w[k] * z[i,t,k] # zero out species with w = 0
      }
      Nspec[i,t] <- sum(tmp[i,t,])# Site-specific species richness
    }
  }
}
")

# Call JAGS from R, check convergence and summarize posteriors
out2B <- jags(bdata, inits, params, "DCM2B.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out2B)
# print(out2B, 3)
summary(out2B)
jags.View(out2B)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Comparison of the two DA parameterisations
print(cbind(out2A$summary[1:10,c(1:2,8:9)],
    out2B$summary[1:10,c(1:2,8:9)]), 3)
#             mean     sd Rhat n.eff   mean      sd Rhat n.eff
# omega      0.537 0.0519 1.04    56  0.545  0.0658 1.00   733
# Ntotal    97.079 6.7752 1.10    25 98.894 10.1222 1.01   676
# mu.lpsi1  -2.918 0.6525 1.00   679 -2.936  0.6752 1.03    85
# sd.lpsi1   4.648 0.9177 1.00   577  4.760  0.9474 1.00   651
# mu.lphi    0.243 0.1827 1.00  3000  0.247  0.1866 1.02   103
# sd.lphi    1.103 0.1670 1.00   817  1.107  0.1597 1.00  1357
# mu.lgamma -1.914 0.1765 1.00  1432 -1.930  0.1802 1.02   115
# sd.lgamma  1.185 0.1587 1.00  3000  1.192  0.1584 1.01   497
# mu.lp     -2.022 0.4869 1.04    50 -2.097  0.5988 1.00  2334
# sd.lp      3.065 0.3859 1.05    46  3.099  0.4561 1.00  1452

# Tally up true Nspec in the simulated data set by aggregating z
Nspec <- apply(dat$z, 1:2, sum)

# Summarize average proportion of species detected per site and year
Nspec.obs <- apply(apply(dat$y, c(1,3,4), max), 1:2, sum)
summary(P <- c(Nspec.obs) / c(out2B$mean$Nspec))
#   Min. 1st Qu. Median   Mean 3rd Qu.   Max.
# 0.1974  0.3097 0.3576 0.3569  0.4041 0.5527


# ~~~~~~~ extra code for figure 5.3 ~~~~~~~~~~~~
# Visualize estimates of community size (Ntotal) and site/year specific estimates
op <- par(mfrow = c(1,2))
hist(out2B$sims.list$Ntotal, col = 'grey', breaks = 30,
    xlab = 'Total number of species in community', ylab = 'Density',
    main = 'Overall species richness', freq = FALSE, xlim = c(50, 181))
abline(v = dat$nspec.det, lwd = 2)
matplot(1:6, t(out2B$mean$Nspec), xlab = 'Year', ylab = 'Species richness',
    main = 'Per site and year', type = 'l',
    lty = 1, lwd = 2, frame = FALSE, col = 'gray50')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
