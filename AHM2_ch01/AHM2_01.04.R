#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

# Approximate time to execute this code: 4 mins

# ~~~~~ Need to run 1.3 before this ~~~~~~~
source("AHM2_01.03.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(coda)

# 1.4 Generalized Linear Models: A site-by-year model
# ===================================================

# Bundle and summarize data
M <- nrow(C)
T <- ncol(C)
str(bdata <- list(C = C, M = M, T = T) )
# List of 3
# $ C: int [1:267, 1:18] 1 0 NA 0 3 NA NA 5 0 0 ...
# $ M: int 267
# $ T: int 18

# Specify model in BUGS language
cat(file = "model1.txt","
model {
  # Priors
  for(i in 1:M){
    site[i] ~ dnorm(0, 0.001) # Priors for site effects
  }
  year[1] <- 0 # Constraint on year effects
  for(t in 2:T){
    year[t] ~ dnorm(0, 0.001) # Priors for year effects 2:T
  }
  # Likelihood
  for (i in 1:M){
    for(t in 1:T){
      C[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- site[i] + year[t]
    }
  }
  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(lambda[,t])
  }
}
")

# Initial values
inits <- function() list(site = rnorm(nrow(C)), year = c(NA, rnorm(ncol(C)-1)))

# Parameters monitored
params <- c("site", "year", "popindex")

# MCMC settings
na <- 1000 ; ni <- 15000; nt <- 10 ; nb <- 5000 ; nc <- 3

# Call JAGS (ART 5 min), check convergence and summarize posteriors
library(jagsUI)
out1 <- jags(bdata, inits, params, "model1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
jags.View(out1) ; print(out1, 2) # Two formats for posterior summaries

# ~~~~~ save output for use in subsequent sections
save(out1, file="AHM2_01.04_out1.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Reformat data for analysis in R using glm
cvec <- c(C)
sitevec <- rep(1:267, 18)
yrvec <- rep(1:18, each = 267)
cbind(sitevec, yrvec, cvec)   # Look at data in this format
summary(fm <- glm(cvec ~ as.factor(sitevec) + as.factor(yrvec) - 1,
    family = 'poisson'))

cor(exp(fm$coef[1:nsite]), exp(out1$summary[1:nsite,1])) # > 0.999
cor(fm$coef[268:284], out1$summary[269:285,1])           # > 0.999

# ~~~~~ code to plot figure 1.4 ~~~~~~~~~~~~~~
# Plot population size index
plot(year, out1$mean$popindex, xlab = 'Year', ylab = 'Population index',
    ylim = c(500, 1000), pch = 16, frame = FALSE, type= 'b', )
segments(year, out1$q2.5$popindex, year, out1$q97.5$popindex)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


params <- c("site", "year", "popindex", "C")
na <- 100 ; ni <- 1500 ; nt <- 2 ; nb <- 500 ; nc <- 2
out1X <- jags(bdata, inits, params, "model1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

(n.missing <- rowSums(is.na(C)) )
(na.sites <- which(n.missing > 0) )
head(C[160:169, 1:13]) # Some of the data with many NAs ...
# ... and the corresponding estimates
head(round(out1X$mean$C[160:169, 1:13], 1))
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
# [1,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0     0     0     0
# [2,]  1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0   0.0     0     0     0
# [3,]  5.0  5.1  4.6  5.5  6.1  6.3  6.6  6.1  6.3   6.6     1     0    10
# [4,]  1.6  4.0  1.0  1.0  2.0  1.0  2.0  1.0  2.0   1.0     2     3     1
# [5,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0     1     1     0
# [6,]  8.0  6.0  3.0  7.0  9.0  8.0 11.0  8.0 19.0  11.0    13    12    17

head(round(out1X$sd$C[160:169, 1:13], 1))
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
# [1,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0     0     0     0
# [2,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0     0     0     0
# [3,]  2.4  2.4  2.2  2.5  2.6  2.7  2.9  2.7  2.6   2.7     0     0     0
# [4,]  1.3  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0     0     0     0
# [5,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0     0     0     0
# [6,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0     0     0     0

