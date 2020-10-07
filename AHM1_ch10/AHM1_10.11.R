#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

# Approximate execution time for this code: 5 mins
# Run time with the full number of iterations: 55 mins

library(AHMbook)
library(jagsUI)

# 10.11 Space-for-time substitution
# =================================


# 10.11.1 A magical covariate
# ------------------------------------------------------------------------
set.seed(1)
data <- sim3Occ(nunit = 500, nsubunit = 5, nrep = 1, mean.psi = 0.8,
    beta.Xpsi = 1, sd.logit.psi = 0.4, mean.theta = 0.6,
    theta.time.range = c(0, 0), beta.Xtheta = 1, sd.logit.theta = 0.6,
    mean.p = 0.4, p.time.range = c(0,0), beta.Xp = -1, sd.logit.p = 0.8)


# Bundle and summarize data set
y <- data$y
str( win.data <- list(y = y, nunit = dim(y)[1], nsubunit = dim(y)[2],
    nrep = dim(y)[3], covA = data$covA, covB = data$covB, covC = data$covC) )

# Define model in BUGS langauge
sink("model.txt")
cat("
model {

  # Priors
  int.psi ~ dunif(0,1)   # Occupancy probability
  int.theta ~ dunif(0,1) # Availability probability
  int.p ~ dunif(0,1)     # Detection probability
  beta.lpsi ~ dnorm(0, 0.01)
  beta.ltheta ~ dnorm(0, 0.01)
  beta.lp ~ dnorm(0, 0.01)

  # Likelihood
  for (i in 1:nunit){
    # Occupancy model for quad i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi) + beta.lpsi * covA[i]
    for (j in 1:nsubunit){
      # Availability in subquad j
      a[i,j] ~ dbern(mu.a[i,j])
      mu.a[i,j] <- z[i] * theta[i,j]
      logit(theta[i,j]) <- logit(int.theta) + beta.ltheta * covB[i,j]
      for (k in 1:nrep){
        # PCR detection error process in replicate k
        y[i,j,k] ~ dbern(mu.y[i,j,k])
        mu.y[i,j,k] <- a[i,j] * p[i,j]
        logit(p[i,j]) <- logit(int.p) + beta.lp * covC[i,j,1]
      }
    }
    tmp[i] <- step(sum(a[i,])-0.1)
  }

  # Derived quantities
  sum.z <- sum(z[])      # Total number of occupied quadrats
  sum.a <- sum(tmp[])    # Total number of quads with presence in samples
  p.theta <- int.p * int.theta # What a 2-level model estimates as 'p'
}
",fill=TRUE)
sink()

# Initial values
inits <- function() list(z = array(1, dim = data$nunit),
    a = array(1, dim =c(data$nunit, data$nsubunit)) )    # Set all to 1 to avoid conflict

# Parameters monitored
params <- c("int.psi", "int.theta", "int.p", "beta.lpsi", "beta.ltheta",
    "beta.lp","p.theta", "sum.z", "sum.a")

# MCMC settings
# ni <- 25000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
ni <- 2500   ;   nt <- 1   ;   nb <- 200   ;   nc <- 3  # ~~~ for testing

# Call JAGS (ART 15 min) and summarize posterior
out <- jags(win.data, inits, params, "model.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
traceplot(out)   ;   print(out, 3)

# Compare truth and estimate in table
tmp <- cbind(rbind(mean.psi = data$mean.psi, mean.theta = data$mean.theta,
    mean.p = data$mean.p, beta.lpsi = data$beta.Xpsi,
    beta.ltheta = data$beta.Xtheta, beta.lp = data$beta.Xp,
    product.theta.p = data$mean.theta* data$mean.p, sum.z = data$sum.z),
    out$summary[c(1:8),c(1, 3, 7)])
colnames(tmp) <- c("Truth", "Post.mean", "LCRL", "UCRL")
print(tmp, 3)


# 10.11.2 No magical covariate known: theta and p are confounded
# ------------------------------------------------------------------------
# Load unmarked and define a function to fit the model with unmarked
library(unmarked)
occUM.fn <- function(data = data, inits = c(1, 1, -1)){
   umf <- unmarkedFrameOccu(y = data$y[,,1], siteCovs =
      data.frame(covA = data$covA))
   tmp1 <- summary(fm <- occu(~1 ~covA, data=umf, starts=inits))
   tmp <- matrix(unlist(c(tmp1$state[,1], tmp1$det[,1], tmp1$state[,2],
      tmp1$det[,2])), ncol = 2, byrow = FALSE)
   dimnames(tmp) <- list(c("Occ_Int", "Occ_A", "Det_Int"), c("MLE", "SE"))
   return(MLE = tmp)
}

# Choose number of simulations and create structures to hold results
# simreps <- 10000                # takes about 30 min
simreps <- 1000                # ~~~~ enough for testing
obs.stats <- array(NA, dim = c(simreps, 3))
dimnames(obs.stats) <- list(NULL, c("sum.z", "obs.sum.z", "sum.z.x"))
MLEs <- array(NA, dim = c(3,2,simreps))
dimnames(MLEs) <- list(c("Occ_Int", "Occ_A", "Det_Int"), c("MLE", "SE"), NULL)

# Set timer and launch simulation
system.time(
for(i in 1:simreps){
  cat("\n\n*** Simrep Number:", i, "***\n\n")
  # Generate data set
  data <- sim3Occ(nunit = 500, nsubunit = 5, nrep = 1, mean.psi = 0.8,
      beta.Xpsi = 1, sd.logit.psi = 0.4, mean.theta = 0.6,
      theta.time.range = c(-1, 1), beta.Xtheta = 0, sd.logit.theta = 0.6,
      mean.p = 0.4, p.time.range = c(0,0), beta.Xp = 0, sd.logit.p = 0.8)
  # Save stats
  obs.stats[i,] <- unlist(data[23:25])
  # Get MLEs of occupancy model and save them
  UMmle <- try(occUM.fn(data = data, inits = c(1,1,-1)))
  # if (class(UMmle) == "try-error") {v<-1} else {  # ~~~ better to use 'inherits'
  if (!inherits(UMmle, "try-error")) {
    MLEs[,,i] <- UMmle
  }
  rm(data, UMmle)
}
)

# ~~~~~ the 'data' object was 'rm'ed in the loop ~~~~~~
# Create a new one with the same arguments
   data <- sim3Occ(nunit = 500, nsubunit = 5, nrep = 1, mean.psi = 0.8,
      beta.Xpsi = 1, sd.logit.psi = 0.4, mean.theta = 0.6,
      theta.time.range = c(-1, 1), beta.Xtheta = 0, sd.logit.theta = 0.6,
      mean.p = 0.4, p.time.range = c(0,0), beta.Xp = 0, sd.logit.p = 0.8)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Visualize results
op <- par(mfrow = c(1,3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5,
    cex.main = 1.5)
# Estimate of occupancy (psi)
hist(plogis(MLEs[1,1,]), breaks = 40, col = "grey",
    main = "Quadrat occupancy (psi)")
abline(v = mean(plogis(MLEs[1,1,]), na.rm = T), col = "blue", lwd = 3)
abline(v = 0.8, col = "red", lwd = 3)
# Estimate of occupancy covariate (A)
hist(MLEs[2,1,], breaks = 40, col = "grey",
    main = "Quadrat occupancy covariate (covA)")
abline(v = mean(MLEs[2,1,], na.rm = T), col = "blue", lwd = 3)
abline(v = data$beta.Xpsi, col = "red", lwd = 3)
# Estimate of "detection": product of theta and p
hist(plogis(MLEs[3,1,]), breaks = 40, col = "grey",
    main = "'Detection probability' = \n theta * p")
abline(v = mean(plogis(MLEs[3,1,]), na.rm = T), col = "blue", lwd = 3)
abline(v = data$mean.theta * data$mean.p, col = "red", lwd = 3, lty = 2)
par(op)
