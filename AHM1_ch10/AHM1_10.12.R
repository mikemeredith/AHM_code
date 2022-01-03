#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

library(AHMbook)
library(R2WinBUGS)
bd <- "C:/WinBUGS14" # location of the "WinBUGS14.exe" application

# 10.12 Models for data along transects: Poisson, exponential,
#       Weibull and removal observation models
# ============================================================


# 10.12.1 Occupancy models with "survival model" observation process:
#         exponential time-to-detection (TTD) model with simulated data
# ------------------------------------------------------------------------
simOccttd(M = 250, mean.psi = 0.4, mean.lambda = 0.3, beta1 = 1,
    alpha1 = -1, Tmax = 10)

set.seed(1)
data <- simOccttd()
str(data)

# Plot response (not shown)
hist(data$ttd, breaks = 50, col = "grey",
    main = "Observed distribution of time to detection",
    xlim = c(0, data$Tmax), xlab = "Measured time to detection")
abline(v = data$Tmax, col = "grey", lwd = 3)

# Bundle data
str( win.data <- list(ttd = data$ttd, d = data$d, covA = data$covA,
   covB = data$covB, nobs = data$M, Tmax = data$Tmax) )

# Define exponential observation model
cat(file = "model1.txt", "
model {

  # Priors
  int.psi ~ dunif(0, 1)               # Intercept occupancy on prob. scale
  beta1 ~ dnorm(0, 0.001)             # Slope coefficient in logit(occupancy)
  int.lambda ~ dgamma(0.0001, 0.0001) # Poisson rate parameter
  alpha1 ~ dnorm(0, 0.001)            # Slope coefficient in log(rate)

  # Likelihood
  for (i in 1:nobs){
    # Model for occurrence
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi) + beta1 * covB[i]

    # Observation model
    # Exponential model for time to detection ignoring censoring
    ttd[i] ~ dexp(lambda[i])
    log(lambda[i]) <- log(int.lambda) + alpha1 * covA[i]
    # Model for censoring due to species absence and ttd>=Tmax
    d[i] ~ dbern(theta[i])
    theta[i] <- z[i] * step(ttd[i] - Tmax) + (1 - z[i])
  }
  # Derived quantities
  n.occ <- sum(z[])                   # Number of occupied sites among M
}
")

# Inits function for some params
# Initialize with z = 1 throughout and
#   all missings due to censoring, rather than non-occurrence
zst <- rep(1, length(win.data$ttd))
ttdst <-rep(win.data$Tmax+1, data$M)
ttdst[win.data$d == 0] <- NA
inits <- function(){list(z =zst, ttd = ttdst, int.psi = runif(1),
    int.lambda = runif(1))}

# Parameters to estimate
params <- c("int.psi", "beta1", "int.lambda", "alpha1", "n.occ")

# MCMC settings
ni <- 12000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART 1.3 min) and summarize posteriors
out1 <- bugs(win.data, inits, params, "model1.txt",
  n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt,
  # debug = TRUE, bugs.directory = bd)
  debug = FALSE, bugs.directory = bd) # ~~~~~ for autotesting
print(out1, dig = 3)


# 10.12.2 Time-to-detection analysis with real data:
#         Weibull occupancy model for the peregrine spring survey
# ------------------------------------------------------------------------
# Code modified slightly to use the ttdPeregrine data set in the AHMbook package
?ttdPeregrine  # check the description of the data
data(ttdPeregrine)
data <- ttdPeregrine
# Manage data and standardize time of day
nobs <- length(data$SiteNumber)     # Number of observations
d <- as.numeric(is.na(data$ttd))    # Censoring indicator
mean.tod <- mean(data$MinOfDay)
sd.tod <- sd(data$MinOfDay)
tod <- (data$MinOfDay -mean.tod) / sd.tod

# Bundle and summarize data set
str( win.data <- list(M = max(data$SiteNumber), site = data$SiteNumber,
   tod = tod, male = as.numeric(data$sex)-1, ttd = data$ttd, d = d, nobs = nobs,
   Tmax = data$Tmax) )


# Define model
cat(file = "model2.txt", "
model {

  # Priors
  psi ~ dunif(0, 1)              # Occupancy intercept
  lambda.int[1] ~ dgamma(0.001, 0.001) # Poisson rate parameter for females
  lambda.int[2] ~ dgamma(0.001, 0.001) # Poisson rate parameter for males
  alpha1 ~ dnorm(0, 0.001)       # Coefficient of time of day (linear)
  alpha2 ~ dnorm(0, 0.001)       # Coefficient of time of day (squared)
  shape ~ dgamma(0.001,0.001)    # Weibull shape
  sexratio ~ dunif(0,1)          # Sex ratio (proportion males)

  # Likelihood
  for (i in 1:M){                # Model for occurrence at site level
    z[i] ~ dbern(psi)
  }

  for (i in 1:nobs){             # Observation model at observation level
    # Weibull model for time to detection ignoring censoring
    ttd[i] ~ dweib(shape, lambda[i])
    log(lambda[i]) <- (1-male[i])*log(lambda.int[1]) + male[i]*log(lambda.int[2])+
        alpha1 * tod[i] + alpha2 * pow(tod[i],2)
    # Model for censoring due to species absence and ttd>=Tmax
    d[i] ~ dbern(theta[i])
    theta[i] <- z[site[i]] * step(ttd[i] - Tmax[i]) + (1 - z[site[i]])
    # Model for sex of unobserved individuals
    male[i] ~ dbern(sexratio)   # Will impute sex for unobserved individuals
  }
  # Derived quantities
  n.occ <- sum(z[])              # Number of occupied sites among M
}
")

# Inits function
zst <- rep(1, win.data$M)
ttdst <-rep(win.data$Tmax+1)
ttdst[win.data$d == 0] <- NA
inits <- function(){list(z =zst, ttd = ttdst, psi = runif(1),
    lambda.int = runif(2), alpha1 = rnorm(1), alpha2 = rnorm(1), shape = runif(1))}

# Parameters to estimate
params <- c("psi", "lambda.int", "alpha1", "alpha2", "n.occ", "z", "sexratio", "shape")

# MCMC settings
ni <- 15000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART 0.6 min) and summarize posteriors
out2 <- bugs(win.data, inits, params, "model2.txt",
  n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt,
  # debug = TRUE, bugs.directory = bd)
  debug = FALSE, bugs.directory = bd)  # ~~~~~ for autotesting
print(out2, dig = 3)


# Predict detection over time of day: prediction cov. runs from 7h to 19h
minutes <- 60:780
pred.tod <- (minutes - mean.tod) / sd.tod   # Standardize as real data

# Predict p over time of day, averaging over sex, and for duration of 10 min
sex.mean <- apply(out2$sims.list$lambda.int, 1, mean)
p.pred1 <- 1 - exp(-exp(log(mean(sex.mean)) + out2$mean$alpha1 * pred.tod +
    out2$mean$alpha2 * pred.tod^2) * 10)

# Predict p for durations of 1-60 min, averaging over time of day and sex
duration <- 1:60
p.pred2 <- 1 - exp(-exp(log(mean(sex.mean))) * duration)

# Visualize analysis
op <- par(mfrow = c(2,2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(data$ttd, breaks = 40, col = "grey",
    xlab = "Time to first detection (min)", main = "")
plot(table(out2$sims.list$n.occ)/length(out2$sims.list$n.occ),
    xlab = "Number of occupied sites", ylab = "Density", frame = FALSE)
plot(minutes, p.pred1, xlab = "Minutes after 6.00 hours (i.e., 7.00 - 19.00h)",
    ylab = "Detection prob.", ylim = c(0.6, 1), type = "l", col = "blue",
    lwd = 3, frame = FALSE)
plot(duration, p.pred2, xlab = "Survey duration (min)",
    ylab = "Detection prob.", ylim = c(0, 1), type = "l", col = "blue",
    lwd = 3, frame = FALSE)
par(op)

# 10.12.3 Occupancy models with removal design observation process (no code)
# ------------------------------------------------------------------------

