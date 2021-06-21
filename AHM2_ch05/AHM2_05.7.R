#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5 : MODELING METACOMMUNITY DYNAMICS USING DYNAMIC COMMUNITY MODELS
# ==========================================================================
# Code from proofs dated 2020-08-19

# Approximate run time for this script: 5 mins
# With full number of iterations: 5.4 hrs

library(AHMbook)
library(jagsUI)

# 5.7 Case study: Norwegian birds
# ===============================

# Get the data set
library(AHMbook)
data(Finnmark)
str(dat <- Finnmark) # Be sure to inspect these data
# List of 4
# $ counts : int [1:37, 1:5, 1:4, 1:17] 1 0 4 3 3 0 2 1 1 0 ...
# ..- attr(*, "dimnames")=List of 4
# .. ..$ sites  : chr [1:37] "AE1K" "AE2K" "AE3K" "AE4K" ...
# .. ..$ visits : chr [1:5] "1" "2" "3" "4" ...
# .. ..$ years  : chr [1:4] "2005" "2006" "2007" "2008"
# .. ..$ species: chr [1:17] "Common Redpoll" "Bluethroat" "Long-tailed Skua"...
# $ timeOfDay: int [1:37, 1:5, 1:4] 1280 1340 1265 1245 1220 165 1230 1260 1280 1305 ...
# ..- attr(*, "dimnames")=List of 3
# .. ..$ site: chr [1:37] "AE1K" "AE2K" "AE3K" "AE4K" ...
# .. ..$ rep : chr [1:5] "1" "2" "3" "4" ...
# .. ..$ year: chr [1:4] "2005" "2006" "2007" "2008"
# $ sites :'data.frame': 37 obs. of 8 variables:
# ..$ region   : Factor w/ 3 levels "Ifjord","Komag",..: 1 1 1 1 1 1 2 2 2 2 ...
# ..$ catchment: Factor w/ 11 levels "IF1","IF2","IF3",..: 1 2 3 4 5 4 6 6 6 6 ...
# ..$ plot     : Factor w/ 37 levels "AE1K","AE2K",..: 1 2 3 4 5 6 7 8 9 10 ...
# ..$ plotnr   : int [1:37] 1 2 3 4 5 6 7 8 9 10 ...
# ..$ area     : num [1:37] 17.9 52.2 27.4 58.9 52.5 ...
# ..$ edge     : num [1:37] 1090 498 571 402 809 ...
# ..$ height   : num [1:37] 77.5 162.5 142.5 220 160 ...
# ..$ density  : num [1:37] 2.5 3.25 3 3 3.25 2.75 5.5 5.5 2.5 0.75 ...
# $ species :'data.frame': 17 obs. of 3 variables:
# ..$ species   : Factor w/ 17 levels "Bluethroat","Brambling",..: 3 1 7 13 5 4 8 ...
# ..$ latin     : Factor w/ 17 levels "Anthus cervinus",..: 6 10 15 3 17 14 2 4 11 9 ...
# ..$ assemblage: Factor w/ 3 levels "OT","WCB","WGB": 2 3 1 1 2 1 1 1 1 3 ...

# Turn counts into detection/nondetection data
str(y <- dat$counts)                     # Copy counts
y[y > 1] <- 1                            # Quantize

# Look at annual number of sites where each species observed
tmp <- apply(y, c(1,3,4), max, na.rm = TRUE)
tmp[tmp == '-Inf'] <- NA
t(apply(tmp, c(2,3), sum, na.rm = TRUE))
#                         years
#                 species 2005 2006 2007 2008
#          Common Redpoll   34   36   34   28
#              Bluethroat   11   15   16   21
#        Long-tailed Skua    2    4    9    1
#    Rough-legged Buzzard    1    2    2    1
#               Fieldfare   13   12   18   15
#  Eurasian Golden Plover    7    6   11   10
#            Meadow Pipit   34   35   34   33
#          Lapland Buntin   22   11   18   14
#           White Wagtail    7   13   13   11
#        Willow Ptarmigan    5    2    7    4
#          Willow Warbler   18   21   24   26
#                 Redwing   26   23   24   28
#            Reed Bunting    3    3    2    2
#       Northern Wheatear    7    5    3    5
#        Temminck's Stint    4    4    5    6
#               Brambling    1    2    6    1
#      Red-throated Pipit    3    3    0    5

# Get sample sizes
nsites <- dim(y)[1]
nsurveys <- dim(y)[2]
nyears <- dim(y)[3]
nspec <- dim(y)[4]

# Get covariates, plot and standardize them
area.o <- dat$sites$area
ragged.o <- dat$sites$edge / area.o
op <- par(mfrow = c(1,3))               # Describe covariates: plot not shown
hist(area.o, breaks = 30, col = 'grey', main = 'Thicket area')
hist(ragged.o, breaks = 30, col = 'grey', main = 'Thicket raggedness')
plot(area.o, ragged.o, pch = 16, frame = FALSE)
area <- standardize(area.o)
ragged <- standardize(ragged.o)
par(op)

# Fill the 3D array Time of day covariate
str(tod <- dat$timeOfDay)                # Grab and inspect survey time data
tod[is.na(tod)] <- mean(tod, na.rm = TRUE)

# Create zero 4D array for M species
nz <- 5                                  # Number of potential 'zero species'
M <- nspec + nz                          # M is the new 'nspec'
yaug <- array(0, dim = c(nsites, nsurveys, nyears, M))
dim(yaug)                                # Check if it went well

# Fill in the observed data into this larger array
yaug[,,,1:dim(y)[4]] <- y
str(yaug)
sum(y, na.rm = TRUE) ; sum(yaug, na.rm = TRUE) # Quick sum check

# Create matrix of site covs and inspect result
head(covs <- cbind(area = area, area2 = area^2, ragged = ragged,
    ragged2 = ragged^2) )

# Bundle and summarize data set
str(bdata <- list(yaug = yaug, nsites = nsites, nsurveys = nsurveys,
    nyears = nyears, M = M, covs = covs, tod = tod, pi = pi, Tday = 1440) )
# List of 9
# $ yaug    : num [1:37, 1:5, 1:4, 1:22] 1 0 1 1 1 0 1 1 1 0 ...
# $ nsites  : int 37
# $ nsurveys: int 5
# $ nyears  : int 4
# $ M       : num 22
# $ covs    : num [1:37, 1:4] -0.5912 1.6329 0.0228 2.0663 1.6547 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:4] "area" "area2" "ragged" "ragged2"
# $ tod     : num [1:37, 1:5, 1:4] 1280 1340 1265 1245 1220 ...
# $ pi      : num 3.14
# $ Tday    : num 1440

# Specify model in BUGS language
cat(file = "DCM6.txt", "
model {

  # *** Priors and hyperpriors for model on psi1 ***
  # Priors
  for(k in 1:M){ # Loop over 22 species in augmented list
    alpha.lpsi1[k] ~ dnorm(mu.alpha.lpsi1, tau.alpha.lpsi1)
      for(g in 1:4){ # Loop over 4 coefficients
      beta.lpsi1[g, k] ~ dnorm(mu.beta.lpsi1[g], tau.beta.lpsi1[g])
    }
  }

  # Hyperpriors
  mu.alpha.lpsi1 <- logit(mean.alpha.psi1)
  mean.alpha.psi1 ~ dunif(0, 1)
  tau.alpha.lpsi1 <- pow(sd.alpha.lpsi1, -2)
  sd.alpha.lpsi1 ~ dunif(0, 10)
  for(g in 1:4){ # Loop over 4 coefficients
    mu.beta.lpsi1[g] ~ dnorm(0, 0.1)
    tau.beta.lpsi1[g] <- pow(sd.beta.lpsi1[g], -2)
    sd.beta.lpsi1[g] ~ dnorm(0, 0.1)I(0,) # Half-Normal prior
    # curve(dnorm(x, 0, sqrt(1 / 0.1)), 0, 20) # howsit look like ?
  }

  # *** Priors and hyperpriors for model on phi ***
  # Priors
  for(k in 1:M){ # Loop over all 22 species
    for(t in 1:(nyears-1)){ # Loop over 3 intervals
      alpha.lphi[t,k] ~ dnorm(mu.alpha.lphi[t], tau.alpha.lphi)
      # phi intercept for 3 intervals, different mean, same variance
    }
    for(g in 1:4){ # Loop over 4 coefficients
      beta.lphi[g, k] ~ dnorm(mu.beta.lphi[g], tau.beta.lphi[g])
    }
  }

  # Hyperpriors
  for(t in 1:(nyears-1)){ # Loop over 3 intervals
    mu.alpha.lphi[t] <- logit(mean.alpha.phi[t])
    mean.alpha.phi[t] ~ dunif(0, 1)
  }
  tau.alpha.lphi <- pow(sd.alpha.lphi, -2)
  sd.alpha.lphi ~ dnorm(0, 0.1)I(0,)
  for(g in 1:4){ # Loop over 4 coefficients
    mu.beta.lphi[g] ~ dnorm(0, 0.01)
    tau.beta.lphi[g] <- pow(sd.beta.lphi[g], -2)
    sd.beta.lphi[g] ~ dnorm(0, 0.1)I(0,)
  }

  # *** Priors and hyperpriors for model on gamma ***
  # Priors
  for(k in 1:M){ # Loop over all 22 species
    for(t in 1:(nyears-1)){ # Loop over 3 intervals
      alpha.lgamma[t,k] ~ dnorm(mu.alpha.lgamma[t], tau.alpha.lgamma)
      # gamma intercept for 3 intervals, different mean, same variance
    }
    for(g in 1:4){ # Loop over 4 coefficients
      beta.lgamma[g, k] ~ dnorm(mu.beta.lgamma[g], tau.beta.lgamma[g])
    }
  }

  # Hyperpriors
  for(t in 1:(nyears-1)){ # Loop over 3 intervals
    mu.alpha.lgamma[t] <- logit(mean.alpha.gamma[t])
    mean.alpha.gamma[t] ~ dunif(0, 1)
  }
  tau.alpha.lgamma <- pow(sd.alpha.lgamma, -2)
  sd.alpha.lgamma ~ dnorm(0, 0.1)I(0,)
  for(g in 1:4){ # Loop over 4 coefficients
    mu.beta.lgamma[g] ~ dnorm(0, 0.1)
    tau.beta.lgamma[g] <- pow(sd.beta.lgamma[g], -2)
    sd.beta.lgamma[g] ~ dnorm(0, 0.1)I(0,)
  }

  # *** Priors and hyperpriors for model on p ***
  # Priors
  for(k in 1:M){ # Loop over all 22 species
    for(t in 1:nyears){ # Loop over 4 years
      alpha.lp[t,k] ~ dnorm(mu.alpha.lp[t], tau.alpha.lp)
      # p intercept for 4 years, different mean, same variance
    }
    for(g in 1:6){ # Loop over (now) 6 coefficients
      beta.lp[g, k] ~ dnorm(mu.beta.lp[g], tau.beta.lp[g]) # coefficients
    }
  }

  # Hyperpriors
  for(t in 1:nyears){ # Loop over 4 years
    mu.alpha.lp[t] <- logit(mean.alpha.p[t])
    mean.alpha.p[t] ~ dunif(0, 1)
  }
  tau.alpha.lp <- pow(sd.alpha.lp, -2)
  sd.alpha.lp ~ dnorm(0, 0.1)I(0,)
  for(g in 1:6){ # Loop over 6 coefficients
    mu.beta.lp[g] ~ dnorm(0, 0.1)
    tau.beta.lp[g] <- pow(sd.beta.lp[g], -2)
    sd.beta.lp[g] ~ dnorm(0, 0.1)I(0,)
  }

  # Likelihood of the model
  # Data augmentation submodel
  omega ~ dunif(0, 1) # Prior for data augmentation parameter
  # omega ~ dbeta(0.001, 1) # Scale prior (Link, 2013) as an alternative
  for(k in 1:M){ # Loop over all 22 species, including 5 all-zero species
    w[k] ~ dbern(omega)
  }

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){ # Loop over 37 sites
    for(k in 1:M){ # Loop over 22 species
      # Initial conditions of system (incl. covariate effects)
      z[i,1, k] ~ dbern(psi1[i, k])
      logit(psi1[i,k]) <- alpha.lpsi1[k] +
      beta.lpsi1[1,k] * covs[i,1] + beta.lpsi1[2,k] * covs[i,2] +
      beta.lpsi1[3,k] * covs[i,3] + beta.lpsi1[4,k] * covs[i,4]

      # State transitions (incl. covariate effects)
      for (t in 2:nyears){ # Loop over years
        z[i,t,k] ~ dbern(z[i,t-1,k]*phi[i,t-1,k] + (1-z[i,t-1, k])*gamma[i,t-1,k])
        logit(phi[i,t-1,k]) <- alpha.lphi[t-1,k] +
        beta.lphi[1,k] * covs[i,1] + beta.lphi[2,k] * covs[i,2] +
        beta.lphi[3,k] * covs[i,3] + beta.lphi[4,k] * covs[i,4]
        logit(gamma[i,t-1,k]) <- alpha.lgamma[t-1,k] +
        beta.lgamma[1,k] * covs[i,1] + beta.lgamma[2,k] * covs[i,2] +
        beta.lgamma[3,k] * covs[i,3] + beta.lgamma[4,k] * covs[i,4]
      }
    }
  }

  # Observation model (incl. covariate effects)
  # Multiplication with w is now here (see Section 5.5)
  # Also note the cosinor function for cyclic timeOfDay (= tod) covariate
  for (i in 1:nsites){
    for(k in 1:M){
      for (j in 1:nsurveys){
        for (t in 1:nyears){
          yaug[i,j,t,k] ~ dbern(w[k] * z[i,t,k] * p[i,j,t,k])
          logit(p[i,j,t,k]) <- alpha.lp[t,k] +
          beta.lp[1,k] * covs[i,1] + beta.lp[2,k] * covs[i,2] +
          beta.lp[3,k] * covs[i,3] + beta.lp[4,k] * covs[i,4] +
          beta.lp[5,k] * cos(2*pi*tod[i,j,t]/Tday) +
          beta.lp[6,k] * sin(2*pi*tod[i,j,t]/Tday)
        }
      }
    }
  }

  # Derived parameters (note multiplication with w)
  # Number of occupied sites
  for(k in 1:M){
    for (t in 1:nyears){
      n.occ[t, k] <- sum(w[k] * z[,t,k])
    }
  }
  # Species richness: total and per site/year
  Ntotal <- sum(w[]) # Total species richness (community size)
  for(i in 1:nsites){
    for(t in 1:nyears){
      for(k in 1:M){
        tmp[i,t,k] <- w[k] * z[i,t,k]
      }
      Nspec[i,t] <- sum(tmp[i,t,]) # Species richness per site and year
    }
  }
}
")

# Initial values (simply initialize all at 1)
zst <- array(1, dim = c(nsites, nyears, M))
wst <- rep(1, M)
inits <- function(){ list(w = wst, z = zst)}

# Parameters monitored
params <- c("omega", "mu.alpha.lpsi1", "sd.alpha.lpsi1", "mu.beta.lpsi1",
    "sd.beta.lpsi1", "mu.alpha.lphi", "sd.alpha.lphi", "mu.beta.lphi",
    "sd.beta.lphi", "mu.alpha.lgamma", "sd.alpha.lgamma", "mu.beta.lgamma",
    "sd.beta.lgamma", "mu.alpha.lp", "sd.alpha.lp", "mu.beta.lp",
    "sd.beta.lp", "alpha.lpsi1", "beta.lpsi1", "alpha.lphi", "beta.lphi",
    "alpha.lgamma", "beta.lgamma", "alpha.lp", "beta.lp", "Ntotal", "Nspec",
    "n.occ")                             # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 150000 ; nt <- 50 ; nb <- 100000 ; nc <- 3
na <- 1000 ; ni <- 1500 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ for testing, 5 mins

# Call JAGS (ART 386 min), check convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "DCM6.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out6)
jags.View(out6) ; print(out6, 2) # not shown

# ~~~~~~~~~ extra code for figure 5.9 ~~~~~~~~~~~~~~~~
# Compare observed and estimated number of species per site and year
tmp <- apply(y, c(1,3,4), max, na.rm = TRUE)
tmp[tmp == '-Inf'] <- NA
Nobs <- apply(tmp, 1:2, sum)
plot(Nobs, out6$mean$Nspec, xlab = 'Observed N', ylab = 'Estimated true N',
    xlim = c(0, 12), ylim = c(0, 16), pch = 16, cex = 1.5, col = 'gray30', frame = FALSE)
abline(0, 1)
(mean.P <- 100*mean(Nobs/out6$mean$Nspec, na.rm = TRUE)) # % detected
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Generate five covariates for prediction
area.predo <- seq(min(area.o), max(area.o), , 100)
ragged.predo <- seq(min(ragged.o), max(ragged.o), , 100)
area.pred <- standardize2match(area.predo, area.o)
ragged.pred <- standardize2match(ragged.predo, ragged.o)
tod.pred <- seq(min(tod), max(tod), , 100)

# Create matrix of prediction site covs with quadratic terms
head(pred.covs <- cbind(area = area.pred, area2 = area.pred^2,
    ragged = ragged.pred, ragged2 = ragged.pred^2))

# Generate arrays to hold the predictions
comm.pred <- array(NA, dim = c(length(area.pred), 3))       # Community
spec.pred <- array(NA, dim = c(length(area.pred), 3, 17))   # Species
str(tmp <- out6$mean)                    # Grab the posterior means
lin <- c(1,3)                            # Index variable for linears
quad <- c(2,4)                           # Index variable for quadratics

# Compute predictions for initial occupancy
for(g in 1:2){                           # For the community mean first ...
  comm.pred[,g] <- plogis(tmp$mu.alpha.lpsi +
  tmp$mu.beta.lpsi1[lin[g]] * pred.covs[,lin[g]] +
  tmp$mu.beta.lpsi1[quad[g]] * pred.covs[,quad[g]])
  for(s in 1:17){                        # ... and then for each species
    spec.pred[,g,s] <- plogis(tmp$alpha.lpsi1[s] +
    tmp$beta.lpsi1[lin[g],s] * pred.covs[,lin[g]] + tmp$beta.lpsi1[quad[g],s] *
    pred.covs[,quad[g]])
  }
}

# Spagetti-plot predictions for initial occupancy
# ( Code shown here only for top row of 4x2 matrix in Fig. 5.10)
cx <- 1.3
op <- par(mfrow = c(1, 2), mar = c(5,5,2,1), cex.lab = cx, cex.axis = cx,
    cex.main = cx)
matplot(area.predo, spec.pred[,1,], xlab = 'Willow thicket area',
    ylab = "Occupancy (psi1)", ylim = c(0, 1), frame = FALSE, las = 1, lty = 1,
    type = 'l', lwd = 3, col = rgb(0,0,0,0.4))
lines(area.predo, comm.pred[,1], lwd = 5)
matplot(ragged.predo, spec.pred[,2,], xlab = ' Willow thicket raggedness',
    ylab = "Occupancy (psi1)", ylim = c(0, 1), frame = FALSE, las = 1, lty = 1,
    type = 'l', lwd = 3, col = rgb(0,0,0,0.4))
lines(ragged.predo, comm.pred[,2], lwd = 5)
par(op)

# ~~~ Extra code for all the plots in figure 5.10 ~~~~~~~~~~~
# Generate arrays to hold the predictions
comm.pred <- array(NA, dim = c(length(area.pred), 2))     # Community
spec.pred <- array(NA, dim = c(length(area.pred), 2, 17)) # Species
str(tmp <- out6$mean)    # Grab the posterior means
lin <- c(1,3)            # To pull out linears
quad <- c(2,4)           # To pull out quads

op <- par(mfrow = c(4,2))

# Predictions for initial occupancy
for(g in 1:2){       # For the community mean first ...
  comm.pred[,g] <- plogis(tmp$mu.alpha.lpsi + tmp$mu.beta.lpsi1[lin[g]] *
      pred.covs[,lin[g]] + tmp$mu.beta.lpsi1[quad[g]] * pred.covs[,quad[g]])
  for(s in 1:17){    # ... and then for each species
    spec.pred[,g,s] <- plogis(tmp$alpha.lpsi1[s] + tmp$beta.lpsi1[lin[g],s] *
        pred.covs[,lin[g]] + tmp$beta.lpsi1[quad[g], s] * pred.covs[,quad[g]])
  }
}
plot(area.predo, comm.pred[,1], xlab = '', ylab = 'Initial occupancy',
    type = 'l', col = 'black', ylim = c(0, 1), frame = FALSE, las = 1)
matlines(area.predo, spec.pred[,1,], lty = 1, type = 'l', col = rgb(0,0,0,0.4))
lines(area.predo, comm.pred[,1], lwd = 2)
plot(ragged.predo, comm.pred[,2], xlab = '', ylab = '', type = 'l',
    col = 'black', ylim = c(0, 1), frame = FALSE, las = 1)
matlines(ragged.predo, spec.pred[,2,], lty = 1, type = 'l',
    col = rgb(0,0,0,0.4))
lines(ragged.predo, comm.pred[,2], lwd = 2)

# Predictions for persistence
for(g in 1:2){     # For the community mean first ...
  comm.pred[,g] <- plogis(mean(tmp$mu.alpha.lphi) + tmp$mu.beta.lphi[lin[g]] *
      pred.covs[,lin[g]] + tmp$mu.beta.lphi[quad[g]] * pred.covs[,quad[g]])
  for(s in 1:17){  # ... and then for each species
    spec.pred[,g,s] <- plogis(mean(tmp$alpha.lphi[,s]) +
        tmp$beta.lphi[lin[g],s] * pred.covs[,lin[g]] +
        tmp$beta.lphi[quad[g], s] * pred.covs[,quad[g]])
  }
}
plot(area.predo, comm.pred[,1], xlab = '', ylab = 'Persistence',
    type = 'l', col = 'black', ylim = c(0.4, 1), frame = FALSE, las = 1)
matlines(area.predo, spec.pred[,1,], lty = 1, type = 'l', col = rgb(0,0,0,0.4))
lines(area.predo, comm.pred[,1], lwd = 2)
plot(ragged.predo, comm.pred[,2], xlab = '', ylab = '', main = '', type = 'l',
    col = 'black', ylim = c(0.4, 1), frame = FALSE, las = 1)
matlines(ragged.predo, spec.pred[,2,], lty = 1, type = 'l',
    col = rgb(0,0,0,0.4))
lines(ragged.predo, comm.pred[,2], lwd = 2)

# Predictions for colonization
for(g in 1:2){     # For the community mean first ...
  comm.pred[,g] <- plogis(mean(tmp$mu.alpha.lgamma) +
      tmp$mu.beta.lgamma[lin[g]] * pred.covs[,lin[g]] +
      tmp$mu.beta.lgamma[quad[g]] * pred.covs[,quad[g]])
  for(s in 1:17){  # ... and then for each species
    spec.pred[,g,s] <- plogis(mean(tmp$alpha.lgamma[,s]) +
        tmp$beta.lgamma[lin[g],s] * pred.covs[,lin[g]] +
        tmp$beta.lgamma[quad[g], s] * pred.covs[,quad[g]])
  }
}
# Plot predictions for colonization
plot(area.predo, comm.pred[,1], xlab = '', ylab = 'Colonization',
    type = 'l', col = 'black', ylim = c(0, 1), frame = FALSE, las = 1)
matlines(area.predo, spec.pred[,1,], lty = 1, type = 'l', col = rgb(0,0,0,0.4))
lines(area.predo, comm.pred[,1], lwd = 2)
plot(ragged.predo, comm.pred[,2], xlab = '', ylab = '', type = 'l',
    col = 'black', ylim = c(0, 1), frame = FALSE, las = 1)
matlines(ragged.predo, spec.pred[,2,], lty = 1, type = 'l',
    col = rgb(0,0,0,0.4))
lines(ragged.predo, comm.pred[,2], lwd = 2)

# Predictions for detection
for(g in 1:2){       # For the community mean first ...
  comm.pred[,g] <- plogis(mean(tmp$mu.alpha.lp) + tmp$mu.beta.lp[lin[g]] *
      pred.covs[,lin[g]] + tmp$mu.beta.lp[quad[g]] * pred.covs[,quad[g]] +
      tmp$mu.beta.lp[5] * cos(2*pi*720/1440) +
      tmp$mu.beta.lp[6] * sin(2*pi*720/1440) )
  for(s in 1:17){    # ... and then for each species
    spec.pred[,g,s] <- plogis(mean(tmp$alpha.lp[,s]) +
        tmp$beta.lp[lin[g],s] * pred.covs[,lin[g]] +
        tmp$beta.lp[quad[g], s] * pred.covs[,quad[g]] +
        tmp$beta.lp[5, s] * cos(2*pi*720/1440) +
        tmp$beta.lp[6, s] * sin(2*pi*720/1440))
  }
}
# Plot predictions for detection
plot(area.predo, comm.pred[,1], xlab = 'Willow patch area',
    ylab = 'Detection', type = 'l', col = 'black', ylim = c(0, 1),
    frame = FALSE, las = 1)
matlines(area.predo, spec.pred[,1,], lty = 1, type = 'l', col = rgb(0,0,0,0.4))
lines(area.predo, comm.pred[,1], lwd = 2)
plot(ragged.predo, comm.pred[,2], xlab = 'Willow patch raggedness',
    ylab = '', type = 'l', col = 'black', ylim = c(0, 1),
    frame = FALSE, las = 1)
matlines(ragged.predo, spec.pred[,2,], lty = 1, type = 'l',
    col = rgb(0,0,0,0.4))
lines(ragged.predo, comm.pred[,2], lwd = 2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~ code for figure 5.11 ~~~~~~~~~~~~~~~~
# Generate arrays to hold the predictions
comm.pred <- array(NA, dim = c(length(tod.pred)))     # Community
spec.pred <- array(NA, dim = c(length(tod.pred), 17)) # Species

# Compute predictions for p (averaging over 4 years
comm.pred <- plogis(mean(tmp$mu.alpha.lp) + tmp$mu.beta.lp[5] *
    cos(2*pi*tod.pred/1440) + tmp$mu.beta.lp[6] * sin(2*pi*tod.pred/1440))
for(s in 1:17){      # Make predictions for every species
  spec.pred[,s] <- plogis(mean(tmp$alpha.lp[,s]) + tmp$beta.lp[5,s] *
      cos(2*pi*tod.pred/1440) + tmp$beta.lp[6,s] * sin(2*pi*tod.pred/1440))
}

# Plot prediction for detection ~ time of day for all 17 species
matplot(tod.pred, spec.pred, xlab = 'Time of Day',
    ylab = "Detection probability (p)", type = 'l', lty = 1,
    col = 'blue', ylim = c(0, 1), frame = FALSE, las = 1, xaxt='n')
lines(tod.pred, comm.pred, lwd = 2)
axis(1, at=c(0,6,12,18,24)*60,
    labels=c("00:00", "06:00", "12:00", "18:00", "24:00"))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
