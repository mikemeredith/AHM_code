#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

# Run time with the full number of iterations: 1 hr

library(jagsUI)

# ~~~~~ Need to run 1.3 before this ~~~~~~~
source("AHM2_01.03.R")
# ~~~~~ and this from 1.4 ~~~~~~~~~~~~~~~~~
M <- nrow(C)
T <- ncol(C)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1.5 Generalized Linear Mixed Models
# ===================================

# 1.5.2 A GLMM with trends in relative abundance
# ----------------------------------------------

# Bundle data (and center Year)
str(bdata <- list(C = C, yr = year - mean(year), M = M, T = T) )
# List of 4
# $ C : int [1:267, 1:18] 1 0 NA 0 3 NA NA 5 0 0 ...
# $ yr: num [1:18] -8.5 -7.5 -6.5 -5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5 ...
# $ M : int 267
# $ T : int 18

# Specify model in BUGS language
cat(file = "model3.txt","
model {

  # 'Priors'
  mu ~ dnorm(0, 0.1) # Grand mean (intercept)
  for(i in 1:M){
    gamma[i] ~ dnorm(mu.gamma, tau.gamma) # Random site-level trends
    site[i] ~ dnorm(0, tau.site) # Random site effects
  }
  mu.gamma ~ dnorm(0, 0.1) # Mean trend
  tau.gamma <- pow(sd.gamma, -2)
  sd.gamma ~ dunif(0, 0.2) # Variability of trends
  tau.site <- pow(sd.site, -2)
  sd.site ~ dunif(0, 3)

  for(t in 1:T){
    year[t] ~ dnorm(0, tau.year) # Random year effects
  }
  tau.year <- pow(sd.year, -2)
  sd.year ~ dunif(0, 2)
  tau <- pow(sd, -2)
  sd ~ dunif(0, 1)

  # 'Likelihood'
  for (i in 1:M){
    for(t in 1:T){
      C[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- mu + gamma[i] * yr[t] + site[i] + year[t] + eps[i,t]
      eps[i,t] ~ dnorm(0, tau) # Overdispersion
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(lambda[,t]) # Population index
    for (i in 1:M){ # Site-specific trends
      pred.lam1[i,t] <- exp(mu + site[i] + gamma[i] * yr[t])
      pred.lam2[i,t] <- exp(mu + site[i] + gamma[i] * yr[t]) / exp(mu + site[i])
    }
  }
}
")

# Initial values
inits <- function() list(mu = rnorm(1), gamma = rnorm(M), site = rnorm(M),
    year = rnorm(T), eps = array(1, dim=c(M, T)))

# Parameters monitored
params <- c("mu", "mu.gamma", "sd.gamma", "sd.site", "sd.year", "sd", "gamma",
    "site", "year", "popindex", "pred.lam1", "pred.lam2")

# MCMC settings
# na <- 5000 ; ni <- 60000 ; nt <- 40 ; nb <- 20000 ; nc <- 3
na <- 5000 ; ni <- 6000 ; nt <- 4 ; nb <- 2000 ; nc <- 3  # ~~~~ for testing

# Call JAGS (ART 88 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "model3.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# par(mfrow = c(3,2))  # ~~~ replaced with 'layout' argument
traceplot(out3, layout=c(3,2))
summary(out3) ; jags.View(out3) ; print(out3$summary[1:800,-c(4:6)], 3)
#                      mean       sd       2.5%     97.5% Rhat n.eff overlap0     f
# mu              -0.456599  0.15677  -0.748280  -0.15097 1.04    55        0 0.998
# mu.gamma         0.013879  0.00612   0.001937   0.02570 1.00  3000        0 0.986
# sd.gamma         0.045594  0.00388   0.038483   0.05369 1.00  2744        0 1.000
# sd.site          2.380585  0.13583   2.134211   2.66447 1.00   443        0 1.000
# sd.year          0.099117  0.02263   0.063627   0.14971 1.00  3000        0 1.000
# sd               0.197952  0.01660   0.164070   0.22954 1.00   511        0 1.000
# gamma[1]        -0.020510  0.04022  -0.100570   0.05683 1.00  3000        1 0.690
# [ ..... ]
# year[18]        -0.081874  0.05552  -0.190255   0.02882 1.00  3000        1 0.936
# popindex[1]    676.204662 26.18911 626.416755 727.11841 1.00  3000        0 1.000
# [ ..... ]
# popindex[18]   888.755919 28.49865 832.911551 944.49179 1.00  1297        0 1.000
# pred.lam1[1,1]   0.465903  0.23323   0.147836   1.01941 1.00  1844        0 1.000

# ~~~~~ code for figure 1.5 ~~~~~~~~~~~~~~~~
op <- par(mfrow = c(1, 2), mar = c(5,5,3,3))
matplot(year, t(out3$mean$pred.lam1), type = 'l', lty = 1, lwd = 2, xlab = 'Year',
    ylab = 'Expected abundance', frame = FALSE, ylim = c(0, 32), las = 1)
matplot(year, t(out3$mean$pred.lam2), type = 'l', lty = 1, lwd = 2, xlab = 'Year',
    ylab = 'Site-specific trends (standardized)', frame = FALSE,
    ylim = c(0, 2.6), las = 1)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

table(out3$summary[7:273,'overlap0'])
# 0 1
# 45 222

summary(apply(out3$sims.list$gamma > 0, 1, sum))
#  Min. 1st Qu. Median  Mean 3rd Qu.  Max.
# 111.0   157.0  166.0 166.2   175.0 216.0

# ~~~ for the comparison, need to load model out2 ~~~~
load("AHM2_01.05.1_out2.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(R2site <- 100* (out2$mean$sd.site - out3$mean$sd.site) / out2$mean$sd.site)
(R2year <- 100* (out2$mean$sd.year - out3$mean$sd.year) / out2$mean$sd.year)
(R2resi <- 100* (out2$mean$sd - out3$mean$sd) / out2$mean$sd)
# [1] -0.5243542
# [1] 13.56628
# [1] 29.6783

# ~~~ Save output for use in subsequent sections ~~~
save(out3, file="AHM2_01.05.2_out3.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
