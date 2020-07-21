#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

library(AHMbook)
library(jagsUI)

# 10.6 A model with lots of covariates: use of R function model.matrix with BUGS
# ==============================================================================


set.seed(148)
data <- simOcc(time.effects = c(0,0), sd.lp = 0, b = 0)
str(data)                         # Look at data object

# Create habitat factor
hab <- c(rep("A", 90), rep("B", 90), rep("C", 87)) # must have M = 267 sites

# Load library, format data and summarize unmarked data frame
library(unmarked)
umf <- unmarkedFrameOccu(
   y = data$y,
   siteCovs = data.frame(elev = data$elev, forest = data$forest, hab = hab),
   obsCovs = list(wind = data$wind))
summary(umf)

summary(fm <- occu(~elev+wind ~elev*forest*hab, data=umf))

# Bundle and summarize data set
HAB <- as.numeric(as.factor(hab))  # Get numeric habitat factor
str( win.data <- list(y = data$y, M = nrow(data$y), J = ncol(data$y),
    elev = data$elev, forest = data$forest, wind = data$wind, HAB = HAB) )

# Specify model in BUGS language
sink("modelA.txt")
cat("
model {

# Priors
  mean.p ~ dunif(0, 1)          # Detection intercept on prob. scale
  alpha0 <- logit(mean.p)       #   same on logit scale
  mean.psi ~ dunif(0, 1)        # Occupancy intercept on prob. scale
  beta0 <- logit(mean.psi)      #   same on logit scale
  for(k in 1:2){                # 2 terms in detection model
    alpha[k] ~ dnorm(0, 0.1)   # Covariates on logit(detection)
  }
  for(k in 1:11){               # 11 terms in occupancy model
    beta[k] ~ dnorm(0, 0.1)    # Covariates on logit(occupancy)
  }

  # Likelihood
  for (i in 1:M) {              # Loop over sites
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- beta0 +                 # occupancy (psi) intercept
      beta[1] * elev[i] +                    # effect of elev
      beta[2] * forest[i] +                  # effect of forest
      beta[3] * equals(HAB[i],2) +           # effect of habitat 2 (= B)
      beta[4] * equals(HAB[i],3) +           # effect of habitat 3 (= C)
      beta[5] * elev[i] * forest[i] +                     # elev:forest
      beta[6] * elev[i] * equals(HAB[i],2) +              # elev:habB
      beta[7] * elev[i] * equals(HAB[i],3) +              # elev:habC
      beta[8] * forest[i] * equals(HAB[i],2) +            # forest:habB
      beta[9] * forest[i] * equals(HAB[i],3) +            # forest:habC
      beta[10] * elev[i] * forest[i] * equals(HAB[i],2) + # elev:forest:habB
      beta[11] * elev[i] * forest[i] * equals(HAB[i],3)   # elev:forest:habC
    for (j in 1:J) {           # Loop over replicates
      y[i,j] ~ dbern(z[i] * p[i,j])        # WinBUGS would need 'straw man' !
      logit(p[i,j]) <- alpha0 +            # detection (p) intercept
           alpha[1] * elev[i] +              # effect of elevation on p
           alpha[2] * wind[i,j]              # effect of wind on p
    }
  }
}
",fill = TRUE)
sink()

# Inits
inits <- function(){list(z = apply(data$y, 1, max), mean.psi = runif(1),
    mean.p = runif(1), alpha = rnorm(2), beta = rnorm(11))}

# Parameters monitored
params <- c("alpha0", "alpha", "beta0", "beta")

# MCMC settings
ni <- 50000   ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3

# Run JAGS (ART 4 min), look at convergence and summarize posteriors
outA <- jags(win.data, inits, params, "modelA.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
traceplot(outA)    ;    print(outA, 3)

# Compare MLEs and SEs with posterior means and sd's
tmp <- summary(fm)
cbind(rbind(tmp$state[1:2], tmp$det[1:2]),
    Post.mean = outA$summary[c(4:15, 1:3), 1],
    Post.sd = outA$summary[c(4:15, 1:3), 2])


# Create design matrix for occupancy covariates and look at it
occDM <- model.matrix(~ data$elev * data$forest * hab)[,-1] # Drop first col.
head(occDM)             # Look at design matrix
str(occDM)

# Bundle and summarize data set
str( win.data <- list(y = data$y, M = nrow(data$y), J = ncol(data$y),
    elev = data$elev, wind = data$wind, occDM = occDM) )

# Specify model in BUGS language
sink("modelB.txt")
cat("
model {

  # Priors
  mean.p ~ dunif(0, 1)          # Detection intercept on prob. scale
  alpha0 <- logit(mean.p)       #   same on logit scale
  mean.psi ~ dunif(0, 1)        # Occupancy intercept on prob. scale
  beta0 <- logit(mean.psi)      #   same on logit scale
  for(k in 1:2){                # 2 terms in detection model
    alpha[k] ~ dnorm(0, 0.1)   # Covariates on logit(detection)
  }
  for(k in 1:11){               # 11 terms in occupancy model
    beta[k] ~ dnorm(0, 0.1)    # Covariates on logit(occupancy)
  }

  # Likelihood
  for (i in 1:M) {
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- beta0 + inprod(beta[], occDM[i,])  # slick !
    for (j in 1:J) {
      y[i,j] ~ dbern(z[i] * p[i,j]) # In WinBUGS need 'straw man'
      logit(p[i,j]) <- alpha0 +     # detection (p) intercept
           alpha[1] * elev[i] +       # effect of elevation on p
           alpha[2] * wind[i,j]       # effect of wind on p
    }
  }
}
",fill = TRUE)
sink()

# Call JAGS from R (ART 3.3 min) and summarize posteriors
outB <- jags(win.data, inits, params, "modelB.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
traceplot(outB)    ;    print(outB, 3)

