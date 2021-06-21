#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

# Run time with the full number of iterations: 8 mins

library(AHMbook)
library(jagsUI)

# ~~~~~ Need to run 1.3 before this ~~~~~~~
source("AHM2_01.03.R")
# ~~~~~ and this from 1.4 ~~~~~~~~~~~~~~~~~
M <- nrow(C)
T <- ncol(C)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1.5 Generalized Linear Mixed Models
# ===================================

# 1.5.3 Adding covariates into a GLMM for relative abundance and modeling a trend
# -------------------------------------------------------------------------------

# Scale some covariates and mean-impute missing values in them
elev.sc <- standardize(dat$elev)       # elevation of site
forest.sc <- standardize(dat$forest)   # forest cover of site
date.sc <- standardize(date)
date.sc[is.na(date.sc)] <- 0           # mean impute
dur.sc <- standardize(dur)
dur.sc[is.na(dur.sc)] <- 0             # mean impute

# Bundle and summarize data
str(bdata <- list(C = C, yr = year - mean(year), elev = elev.sc,
    forest = forest.sc, date = date.sc, dur = dur.sc,
    twosurveys = as.numeric(dat$nsurveys == 2), M = M, T = T) )
# List of 9
# $ C         : int [1:267, 1:18] 1 0 NA 0 3 NA NA 5 0 0 ...
# $ yr        : num [1:18] -8.5 -7.5 -6.5 -5.5 -4.5 -3.5 -2.5 -1.5 ...
# $ elev      : num [1:267, 1] -1.1539 -1.1539 -0.2175 -0.3735 -0.0614 ...
# $ forest    : num [1:267, 1] -1.1471 -0.4967 -0.0992 -0.9303 0.0092 ...
# $ date      : num [1:267, 1:18] -0.157 -1.012 0 -0.257 -0.634 ...
# $ dur       : num [1:267, 1:18] -0.8792 -1.1884 0 0.0485 -0.5418 ...
# $ twosurveys: num [1:267] 0 0 0 0 0 0 0 0 0 0 ...
# $ M         : int 267
# $ T         : int 18

# Specify model in BUGS language
cat(file = "model4.txt","
model {

  # 'Priors' and linear models
  mu ~ dnorm(0, 0.1)               # Grand mean (intercept)
  for(i in 1:M){
    site[i] ~ dnorm(0, tau.site)   # Random site effects
  }

  # Linear model for effect of elevation on expectation of trends
  for(i in 1:M){ # NOTE: here we model the trends
    gamma[i] ~ dnorm(mu.gamma[i], tau.gamma) # Random site-level trends
    mu.gamma[i] <- alpha.mu.gamma + beta1.mu.gamma * elev[i] +
        beta2.mu.gamma * pow(elev[i],2)
  }
  alpha.mu.gamma ~ dnorm(0, 0.1)   # intercept of mean trend on elev
  beta1.mu.gamma ~ dnorm(0, 0.1)   # lin effect of elev on trend
  beta2.mu.gamma ~ dnorm(0, 0.1)   # quad effect of elev on trend
  tau.gamma <- pow(sd.gamma, -2)
  sd.gamma ~ dunif(0, 0.2)         # Variability of trends

  # Other priors
  tau.site <- pow(sd.site, -2)
  sd.site ~ dunif(0, 3)
  for(i in 1:7){
    theta[i] ~ dnorm(0, 0.1)       # Covariate effects
  }
  for(t in 1:T){
    year[t] ~ dnorm(0, tau.year)   # Random year effects
  }
  tau.year <- pow(sd.year, -2)
  sd.year ~ dunif(0, 2)
  tau <- pow(sd, -2)
  sd ~ dunif(0, 1)

  # 'Likelihood'
  for (i in 1:M){
    for(t in 1:T){
      C[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- mu + gamma[i] * yr[t] +
          theta[1] * elev[i] + theta[2] * pow(elev[i],2) +
          theta[3] * forest[i] + theta[4] * date[i,t] +
          theta[5] * pow(date[i,t],2) + theta[6] * dur[i,t] +
          theta[7] * twosurveys[i] + site[i] + year[t] + eps[i,t]
      eps[i,t] ~ dnorm(0, tau)
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(lambda[,t])
  }
}
")

# Initial values
inits <- function() list(mu = rnorm(1), gamma = rnorm(M), theta = rnorm(7),
    site = rnorm(M), year = rnorm(T), eps = array(1, dim=c(M, T)))

# Parameters monitored
params <- c("mu", "alpha.mu.gamma", "beta1.mu.gamma", "beta2.mu.gamma",
    "sd.beta", "theta", "sd.site", "sd.year", "sd", "popindex")
# could also monitor some random effects: "gamma", "site", "year",

# MCMC settings
na <- 5000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3

# Call JAGS (ART 7 min), check convergence and summarize posteriors
out4 <- jags(bdata, inits, params, "model4.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,2))  # ~~~ replaced with 'layout' argument
traceplot(out4, layout=c(3,2))
print(out4, 2)
#                  mean    sd   2.5%    50%  97.5% overlap0     f  Rhat n.eff
# mu              1.018 0.137  0.754  1.016  1.295    FALSE 1.000 1.008   240
# alpha.mu.gamma  0.013 0.007 -0.002  0.013  0.028     TRUE 0.960 1.001  1652
# beta1.mu.gamma  0.019 0.005  0.009  0.019  0.029    FALSE 1.000 1.015   145
# beta2.mu.gamma -0.001 0.007 -0.015 -0.001  0.013     TRUE 0.558 1.003  2368
# theta[1]        0.929 0.117  0.701  0.930  1.158    FALSE 1.000 1.000  3000
# theta[2]       -0.919 0.136 -1.192 -0.913 -0.649    FALSE 1.000 1.029    77
# theta[3]        0.667 0.101  0.467  0.667  0.866    FALSE 1.000 1.002   832
# theta[4]       -0.035 0.027 -0.089 -0.035  0.020     TRUE 0.903 1.003   744
# theta[5]        0.010 0.019 -0.027  0.011  0.049     TRUE 0.710 1.002   914
# theta[6]        0.125 0.020  0.084  0.124  0.165    FALSE 1.000 1.000  3000
# theta[7]       -3.039 0.393 -3.803 -3.041 -2.273    FALSE 1.000 1.000  3000
# sd.site         1.199 0.070  1.069  1.195  1.341    FALSE 1.000 1.011   185
# sd.year         0.103 0.023  0.067  0.100  0.158    FALSE 1.000 1.001  3000
# sd              0.193 0.017  0.161  0.193  0.226    FALSE 1.000 1.039   180
# ........

# ~~~~~~ extra code for these results ~~~~~~~~~~~~~~~~~~~~
load("AHM2_01.05.1_out2.RData")
(R2site <- 100* (out2$mean$sd.site - out4$mean$sd.site) / out2$mean$sd.site)
(R2year <- 100* (out2$mean$sd.year - out4$mean$sd.year) / out2$mean$sd.year)
(R2resi <- 100* (out2$mean$sd - out4$mean$sd) / out2$mean$sd)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# [1] 49.35703
# [1] 10.07
# [1] 31.38612

# ~~~ Save output for use in subsequent sections ~~~
save(out4, file="AHM2_01.05.3_out4.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
