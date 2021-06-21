#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

# Approximate execution time for this code: 6 mins

library(jagsUI)

# ~~~~~ Need to run 1.3 before this ~~~~~~~
source("AHM2_01.03.R")
# ~~~~~ and this from 1.4 ~~~~~~~~~~~~~~~~~
M <- nrow(C)
T <- ncol(C)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1.5 Generalized Linear Mixed Models
# ===================================

# 1.5.1 A GLMM to partition the total variability in relative abundance
# ---------------------------------------------------------------------

# Bundle and summarize data (same as before)
str(bdata <- list(C = C, M = M, T = T))

# Specify model in BUGS language
cat(file = "model2.txt","
model {
  # 'Priors'
  mu ~ dnorm(0, 0.001) # Intercept
  for(i in 1:M){
    site[i] ~ dnorm(0, tau.site) # Random site effects
  }
  tau.site <- pow(sd.site, -2)
  sd.site ~ dunif(0, 5)
  for(t in 1:T){
    year[t] ~ dnorm(0, tau.year) # Random year effects
  }
  tau.year <- pow(sd.year, -2)
  sd.year ~ dunif(0, 3)
  tau <- pow(sd, -2)
  sd ~ dunif(0, 3)

  # 'Likelihood'
  for (i in 1:M){
    for(t in 1:T){
      C[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- mu + site[i] + year[t] + eps[i,t]
      eps[i,t] ~ dnorm(0, tau) # 'Overdispersion'
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(lambda[,t])
  }
}
")

# Initial values
inits <- function() list(mu = rnorm(1), site = rnorm(M), year = rnorm(T),
    eps = array(1, dim=c(M, T)))

# Parameters monitored
params <- c("mu", "sd.site", "sd.year", "sd", "site", "year", "popindex")

# MCMC settings
na <- 1000 ; ni <- 15000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

# Call JAGS (ART 4 min), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "model2.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# par(mfrow = c(3, 2))  #  ~~~ replace with 'layout' argument
traceplot(out2, layout=c(3,2))
summary(out2) ; jags.View(out2) ; print(out2, 3)
#           mean    sd   2.5%    50%  97.5% overlap0     f  Rhat n.eff
# mu      -0.434 0.157 -0.749 -0.431 -0.130    FALSE 0.997 1.023    92
# sd.site  2.378 0.141  2.115  2.371  2.659    FALSE 1.000 1.084    29
# sd.year  0.116 0.025  0.075  0.112  0.177    FALSE 1.000 1.001  2168
# sd       0.282 0.015  0.252  0.282  0.310    FALSE 1.000 1.007   444
# site[1] -0.603 0.416 -1.466 -0.594  0.149     TRUE 0.929 1.001  3000
# [ ...... ]

# Save output for use in subsequent sections
save(out2, file="AHM2_01.05.1_out2.RData")

# ~~~~~~~~~ plot to compare results of out1 and out2 (figures not shown) ~~~~~~~~~~~~~~~~~
load("AHM2_01.04_out1.RData")

# Plot population index under models 1 and 2
plot(year-0.1, out1$mean$popindex, cex = 2, pch = 16, xlab = 'Year',
    ylab = 'Population index', col = 'red', type = 'b', frame = FALSE,
    ylim = c(500, 1000), main = 'red: fixed-effects, blue: random-effects')
segments(year-0.1, out1$q2.5$popindex, year-0.1, out1$q97.5$popindex, col = 'red')
points(year+0.1, out2$mean$popindex, cex = 2, pch = 16, col = 'blue', type = 'b')
segments(year+0.1, out2$q2.5$popindex, year+0.1, out2$q97.5$popindex, col = 'blue')

# Graphical partitioning of the variance
library(denstrip)
op <- par(mar = c(4, 8, 4, 2))
plot(out2$sims.list$sd, xlim = c(0, 3), ylim = c(0.5, 3.5), xlab="", ylab="",
    type="n", axes = FALSE,
    main = "Variance partitioning in space, time and space-time")
axis(1)
axis(2, at = 1:3, labels = c('Space-Time (sd)', 'Time (sd.year)', 'Space (sd.site)'), las = 1)
for(k in 1:3){
   denstrip(unlist(out2$sims.list[k+1]), at = k, ticks = out2$summary[k+1, c(3,5,7)])
}
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~