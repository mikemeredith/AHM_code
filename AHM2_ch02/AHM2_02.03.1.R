#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 20 mins

library(jagsUI)

# ~~~~ need the Green Woodpecker data prepared in 2.2 ~~~~~~~~
source("AHM2_02.02.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.3 Year-stratified N-mixture model
# ===========================================

# 2.3.1 Adding covariates and estimating a common trend over time
# ---------------------------------------------------------------

# Bundle data
str(bdata <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3], elev = elev, forest = forest, DATE = DATE,
    INT = INT))
# List of 8
# $ C       : int [1:267, 1:3, 1:14] 0 3 0 0 0 0 0 0 0 0 ...
# $ nsites  : int 267
# $ nsurveys: int 3
# $ nyears  : int 14
# $ elev    : num [1:267] -1.206 -1.16 -0.184 -0.439 -0.126 ...
# $ forest  : num [1:267] -1.1529 -0.467 0.0023 -0.9002 -0.106 ...
# $ DATE    : num [1:267, 1:3, 1:14] -1.09 -1.32 -1.23 -1.27 -1.36 ...
# $ INT     : num [1:267, 1:3, 1:14] -0.532 -0.959 0.168 -0.452 -0.865 ...

# Specify model in BUGS language
cat(file = "Nmix2.txt","
model {
  # Priors
  for (t in 1:nyears){
    alpha0[t] <- logit(mean.p[t])
    mean.p[t] ~ dunif(0, 1)
  }
  for (k in 1:3){
    alpha[k] ~ dnorm(0, 0.01)
  }
  beta0 ~ dnorm(0, 0.01)
  for (k in 1:4){
    beta[k] ~ dnorm(0, 0.01)
  }

  # Likelihood
  # Ecological model for true abundance
  for (i in 1:nsites){
    for(t in 1:nyears){
      N[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- beta0 + beta[1] * elev[i]+ beta[2] *
      pow(elev[i],2) + beta[3] * forest[i] + beta[4] * (t-7.5)
      # Observation model for replicated counts
      for (j in 1:nsurveys){
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t])
        logit(p[i,j,t]) <- alpha0[t] + alpha[1] * DATE[i,j,t] + alpha[2] *
        pow(DATE[i,j,t],2) + alpha[3] * INT[i,j,t]
      }
    }
  }

  # Derived quantity: Total abundance across all surveyed sites
  for (t in 1:nyears){
    totalN[t] <- sum(N[,t])
  }
}
")

# Initial values
Nst <- apply(C, c(1,3), max, na.rm = TRUE)+1 # Inits for latent N
Nst[Nst == '-Inf'] <- 1
inits <- function() list(N = Nst, alpha = c(runif(2), NA),
    beta0 = runif(1), beta = runif(4))

# Parameters monitored
params <- c("alpha0", "alpha", "beta0", "beta", "totalN")

# MCMC settings
na <- 1000 ; ni <- 15000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

# Run JAGS (ART 22 mins), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "Nmix2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out2)
print(out2, digits = 2) # shown partially only
#           mean   sd  2.5%   50% 97.5% overlap0    f Rhat n.eff
# [...output truncated...]
# alpha[1] -0.22 0.03 -0.27 -0.22 -0.17    FALSE 1.00 1.00   394
# alpha[2]  0.06 0.02  0.01  0.06  0.10    FALSE 0.99 1.00  3000
# alpha[3]  0.35 0.03  0.28  0.35  0.41    FALSE 1.00 1.02   182
# beta0     0.17 0.04  0.10  0.17  0.24    FALSE 1.00 1.00  1127
# beta[1]  -0.59 0.03 -0.64 -0.59 -0.54    FALSE 1.00 1.00  3000
# beta[2]  -0.17 0.03 -0.23 -0.17 -0.11    FALSE 1.00 1.00  2578
# beta[3]   0.13 0.02  0.10  0.13  0.17    FALSE 1.00 1.00  2422
# beta[4]   0.05 0.01  0.03  0.05  0.06    FALSE 1.00 1.01  2914
# [... output truncated ...]


# ~~~~~~~~~~~~ extra code for figure 2.3 ~~~~~~~~~~~~~~~~~~~~~
# Quick plots of response curves of lambda and p to all covariates
op <- par(mfrow = c(2, 3), mar = c(5, 5, 3, 3))
ylim <- c(0, 2)
curve(exp(out2$mean$beta0 + out2$mean$beta[1] * (x - mean(elevo)) / sd(elevo) +
    out2$mean$beta[2] * ((x - mean(elevo)) / sd(elevo))^2), 200, 3000,
    xlab = 'Elevation (m)', ylab = 'E(N)', main = 'lambda ~ elev',
    frame = FALSE, ylim = ylim)
curve(exp(out2$mean$beta0 + out2$mean$beta[3] * (x - mean(foresto)) / sd(foresto)),
    0, 100, xlab = 'Forest cover (%)', ylab = 'E(N)', main = 'lambda ~ forest',
    frame = FALSE, ylim = ylim)
curve(exp(out2$mean$beta0 + out2$mean$beta[4] * (x - 7.5)), 1, 15,
    xlab = 'Time (years)', ylab = 'E(N)',
    main = 'lambda ~ year (Time trend in abundance)',
    frame = FALSE, ylim = ylim)
ylim = c(0, 1)
plot(2004:2017, plogis(out2$mean$alpha0), xlab = 'Year', ylab = 'p',
    main = 'p ~ year (Detection intercepts)', cex = 2, pch = 16,
    frame = FALSE, ylim = ylim)
# curve(plogis(mean(out2$mean$alpha0) + out2$mean$alpha[1] * (x - mean.date) / sd.date + out2$mean$alpha[2] * ((x - mean.date) / sd.date)^2), 100, 200, xlab = 'Julian date', ylab = 'p', main = 'p ~ date', frame = FALSE, ylim = ylim)
curve(plogis(mean(out2$mean$alpha0) + out2$mean$alpha[1] * standardize2match(x, dateso) +
    out2$mean$alpha[2] * standardize2match(x, dateso)^2), 100, 200,
    xlab = 'Julian date', ylab = 'p', main = 'p ~ date', frame = FALSE,
    ylim = ylim)
curve(plogis(mean(out2$mean$alpha0) + out2$mean$alpha[3] * standardize2match(x, dateso)),
    0, 200, xlab = 'Survey intensity (min / km)', ylab = 'p',
    main = 'p ~ intensity', frame = FALSE, ylim = ylim)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
