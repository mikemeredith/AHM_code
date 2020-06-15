#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-01-09

library(AHMbook)
library(jagsUI)

# 1.8 MODELING POPULATION DYNAMICS AT TWO TEMPORAL SCALES
# =======================================================

# 1.8.2 MODELING THE TWO-SCALE DYNAMICS OF BRITISH BUTTERFLIES
# ------------------------------------------------------------

# Load Marbled White data and do some data summaries (fig not shown)
library(AHMbook)
data(UKmarbledWhite)
str(white <- UKmarbledWhite)
par(mfrow = c(2,3))
tmp <- tapply(white$C, list(white$year, white$site), function(x) length(x))
plot(table(c(tmp)), main = "Counts per site and year")
plot(table(white$site), main = "Counts per site")
plot(table(white$north), main = "Counts per latitude class")
plot(table(white$date), main = "Counts per date")
plot(table(white$year), main = "Counts per year")
plot(table(white$C), main = "Frequency distribution of counts")

# Variables for use later
year <- white$year - 1990
# Get latitude for each site (n [ 80) and scale it
latitude <- tapply(white$north, white$site, mean)
mean.lat <- mean(latitude)
sd.lat <- sd(latitude)
lat <- (latitude - mean.lat ) / sd.lat
# Get standardized date
jdate <- white$date - 150 # 1 = Day 151
date <- (jdate/ 10) - mean((jdate/ 10)) # Express in units of 10 days

# Data bundle
str(bdata <- list(C = white$C, pop = white$site, year = year, date = date, lat = lat,
    lat2 = lat^2, npop = max(white$site), nyear = max(year), nobs = length(year),
    pi = pi) )
# List of 10
# $ C : int [1:9651] 15 59 53 4 4 16 25 63 74 28 ...
# $ pop : num [1:9651] 1 1 1 1 1 1 1 1 1 1 ...
# $ year : num [1:9651] 1 1 1 1 1 2 2 2 2 2 ...
# $ date : num [1:9651] -0.546 0.654 1.254 2.954 3.454 ...
# $ lat : num [1:80(1d)] -0.415 -0.646 1.004 0.961 1.112 ...
# $ lat2 : num [1:80(1d)] 0.172 0.418 1.009 0.924 1.236 ...
# $ npop : num 80
# $ nyear: num 25
# $ nobs : int 9651
# $ pi : num 3.14

# Specify model in BUGS language
cat(file = "modelPH.txt","
model {
  ### Priors and linear models for parameters
  # Initial abundance
  for(i in 1:npop){
    expn1[i] <- exp(loglam[i])
    loglam[i] ~ dnorm(mu.llam[i], tau.llam)
    mu.llam[i] <- alpha.llam + beta.llam[1]*lat[i] + beta.llam[2]*lat2[i]
  }
  alpha.llam <- log(mean.lambda)
  mean.lambda ~ dunif(1, 300)
  beta.llam[1] ~ dnorm(0, 0.1)
  beta.llam[2] ~ dnorm(0, 0.1)
  tau.llam <- pow(sd.llam, -2)
  sd.llam ~ dunif(0.001, 5)
  # Interannual population growth rate
  for(t in 1:(nyear-1)){
    for(i in 1:npop){
      gamma[i, t] <- exp(loggam[i, t])
      loggam[i,t] <- alpha0.lgam + alpha.lgam.site[i] + alpha.lgam.time[t] +
          beta.lgam[1] * (t-13) + beta.lgam[2]*lat[i] + beta.lgam[3]*lat2[i] +
          beta.lgam[4]*lat[i] * (t-13) + beta.lgam[5]*lat2[i]*(t-13)
    }
  }
  alpha0.lgam <- log(mean.gamma)
  mean.gamma ~ dunif(0.001, 3)
  for(i in 1:npop){
    alpha.lgam.site[i] ~ dnorm(0, tau.lgam.site)
  }
  for(t in 1:(nyear-1)){
    alpha.lgam.time[t] ~ dnorm(0, tau.lgam.time)
  }
  for(v in 1:5){
    beta.lgam[v] ~ dnorm(0, 1)
  }
  tau.lgam.site <- pow(sd.lgam.site, -2)
  sd.lgam.site ~ dunif(0.001, 1)
  tau.lgam.time <- pow(sd.lgam.time, -2)
  sd.lgam.time ~ dunif(0.001, 1)
  # Mean (or peak) activity period
  for(t in 1:nyear){
    for(i in 1:npop){ # Peak activity period
      mu[i, t] <- alpha0.mu + alpha.mu.site[i] + alpha.mu.time[t] + beta.mu[1] *
      (t-13) + beta.mu[2]*lat[i] + beta.mu[3]*lat2[i] + beta.mu[4]*lat[i]*(t-13) +
      beta.mu[5]*lat2[i]*(t-13)
    }
  }
  alpha0.mu <- mean.mu
  mean.mu ~ dunif(-3, 3)
  for(i in 1:npop){
    alpha.mu.site[i] ~ dnorm(0, tau.mu.site)
  }
  for(t in 1:nyear){
    alpha.mu.time[t] ~ dnorm(0, tau.mu.time)
  }
  for(v in 1:5){
    beta.mu[v] ~ dnorm(0, 1)
  }
  tau.mu.site <- pow(sd.mu.site, -2)
  sd.mu.site ~ dunif(0.001, 1)
  tau.mu.time <- pow(sd.mu.time, -2)
  sd.mu.time ~ dunif(0.001, 1)
  # Half-length of activity period
  for(t in 1:nyear){
    for(i in 1:npop){ # Length of activity period
      sigma[i, t] <- exp(lsig[i,t]) # Width
      lsig[i,t] <- alpha0.lsig + alpha.lsig.site[i] + alpha.lsig.time[t] +
          beta.lsig[1] * (t-13) + beta.lsig[2]*lat[i] + beta.lsig[3]*lat2[i] +
          beta.lsig[4]*lat[i] * (t-13) + beta.lsig[5]*lat2[i]*(t-13)
    }
  }
  alpha0.lsig <- log(mean.sigma)
  mean.sigma ~ dunif(0.001, 3)
  for(i in 1:npop){
    alpha.lsig.site[i] ~ dnorm(0, tau.lsig.site)
  }
  for(t in 1:nyear){
    alpha.lsig.time[t] ~ dnorm(0, tau.lsig.time)
  }
  for(v in 1:5){
    beta.lsig[v] ~ dnorm(0, 1)
  }
  tau.lsig.site <- pow(sd.lsig.site, -2)
  sd.lsig.site ~ dunif(0.001, 1)
  tau.lsig.time <- pow(sd.lsig.time, -2)
  sd.lsig.time ~ dunif(0.001, 1)
  # ’Likelihood’
  # Model for between-year dynamics
  for(i in 1:npop){
    # Initial year
    n[i,1] ~ dpois(expn1[i])
    n1[i] <- n[i,1]
    # Autoregressive transitions from t to t+1
    for(t in 2:nyear){
      n[i,t] ~ dpois(gamma[i, t-1]*n[i,(t-1)])
    }
  }
  # Phenomenological within-season population model
  for(i in 1:nobs){
    C[i] ~ dpois(lambda[i])
    lambda[i] <- n[pop[i],year[i]]*(1 / (sigma[pop[i],year[i]]*sqrt(2*pi)) )*exp( -
        pow((date[i] - mu[pop[i], year[i]]),2) / (2*pow(sigma[pop[i],year[i]], 2)) )
  }
}
")
# Initial values
nst <- tapply(bdata$C, list(bdata$pop, bdata$year), max, na.rm = TRUE)
nst[is.na(nst)] <- round(mean(nst, na.rm = TRUE))
inits <- function() list(n = 2*nst)
# Parameters monitored
# Choose if want both hyperparams and latent variables
params <- c("mu.llam", "mean.lambda", "alpha.llam", "beta.llam", "sd.llam",
    "alpha0.lgam", "mean.gamma", "alpha.lgam", "beta.lgam", "sd.lgam.site",
    "sd.lgam.time", "mean.mu", "alpha0.mu", "beta.mu", "sd.mu.site", "sd.mu.time",
    "mean.sigma", "mu.lsig", "alpha0.lsig", "beta.lsig", "sd.lsig.site", "sd.lsig.time",
    "expn1", "n1", "n", "alpha.lgam.site", "alpha.lgam.time", "gamma", "alpha.mu.site",
    "alpha.mu.time", "mu", "alpha.lsig.site", "alpha.lsig.time", "sigma")
# MCMC settings
# na <- 10000 ; ni <- 150000 ; nt <- 50 ; nb <- 100000 ; nc <- 3
na <- 1000 ; ni <- 15000 ; nt <- 5 ; nb <- 10000 ; nc <- 3 # ~~~~ for testing
# Call JAGS (ART 33 hours)
out13 <- jags(bdata, inits, params, "modelPH.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
summary(out13) ; View(out13)
# Convergence check
par(mfrow = c(3,3), mar = c(3,3,3,2)) ; traceplot(out13) # All params
# Posterior summary only for parameters with Rhat > 1.1
print(out13$summary[which(out13$summary[,8] > 1.1), -c(4:6)], 3)

# Select some key parameters to summarize their posteriors
dim(out13$summary)
cbind(1:2000, out13$summary[1:2000, -c(4:6)])
sel.params <- c(81:112)
print(out13$summary[sel.params,-c(4:6)], 3)
