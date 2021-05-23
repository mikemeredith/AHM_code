#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# Approximate execution time for this code: 22 mins

library(AHMbook)
library(unmarked)
library(jagsUI)

# 6.8 Goodness of fit
# ===================
# ~~~~~~~~~~~ specify number of cores for Nmix.gof.test ~~~~~~~~~~
# Default is to use all-but-one available cores, but that leads to
# a crash if other applications are using multiple cores
ncores <- 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Case 1: Test GoF of correct model
library(AICcmodavg)
op <- par(mfrow = c(3,3))
for(i in 1:9){
  data <- simNmix(show.plot = FALSE)          # Create data set
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100, ncores = ncores) # 100 bootstrap reps
}

# Case 2: Simulate data with zero inflation and analyse without
val.range <- seq(0.1, 1,,9)            # Much to no zero-inflation
for(i in 1:9){
  data <- simNmix(mean.theta = val.range[i], show.plot = FALSE)
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100, ncores = ncores)
}

# Case 3: Extra-Poisson dispersion in lambda
val.range <- seq(1, 0,,9)            # Some to no extra-Poisson dispersion
for(i in 1:9){
   data <- simNmix(sigma.lam = val.range[i], show.plot = FALSE)
   fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
   pb.gof <- Nmix.gof.test(fm, nsim = 100, ncores = ncores)
}

# Case 4: Site covariate in lambda
val.range <- seq(3, 0,,9)            # Strong to no effect of covariate
for(i in 1:9){
   data <- simNmix(beta3.lam = val.range[i], show.plot = FALSE)
   fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
   pb.gof <- Nmix.gof.test(fm, nsim = 100, ncores = ncores)
}

# Case 5: Extra-binomial dispersion in p (survey random effect)
val.range <- seq(1, 0,,9)            # Strong to no effect extra-dispersion
for(i in 1:9){
   data <- simNmix(sigma.p.survey = val.range[i], show.plot = FALSE)
   fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
   pb.gof <- Nmix.gof.test(fm, nsim = 100, ncores = ncores)
}

# Case 6: Site covariate in p
val.range <- seq(3, 0,,9)            # Strong to no covariate effect
for(i in 1:9){
   data <- simNmix(beta3.p = val.range[i], show.plot = FALSE)
   fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
   pb.gof <- Nmix.gof.test(fm, nsim = 100, ncores = ncores)
}

# Case 7: Observational covariate in p
val.range <- seq(3, 0,,9)            # Strong to no covariate effect
for(i in 1:9){
   data <- simNmix(beta.p.survey = val.range[i], show.plot = FALSE)
   fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
   pb.gof <- Nmix.gof.test(fm, nsim = 100, ncores = ncores)
}
par(op)

# Bundle and summarize data set
str( win.data <- list(C = data$C, M = nrow(data$C), J = ncol(data$C), e = 0.001))

# Specify model in BUGS language
sink("model.txt")
cat("
model {
  # Priors
  lambda ~ dgamma(0.001, 0.001)
  p ~ dunif(0, 1)
  # Likelihood
  for (i in 1:M) {
    N[i] ~ dpois(lambda)      # State model
    for (j in 1:J) {
      C[i,j] ~ dbin(p, N[i]) # Observation model
    }
  }
  # Posterior predictive distributions of chi2 discrepancy
  for (i in 1:M) {
    for (j in 1:J) {
      C.sim[i,j] ~ dbin(p, N[i]) # Create new data set under model
      e.count[i,j] <- N[i] * p   # Expected datum
      # Chi-square discrepancy for the actual data
      chi2.actual[i,j] <- pow((C[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
      # Chi-square discrepancy for the simulated ('perfect') data
      chi2.sim[i,j] <- pow((C.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
      # Add small value e to denominator to avoid division by zero
    }
  }
  # Add up individual chi2 values for overall fit statistic
  fit.actual <- sum(chi2.actual[,])  # Fit statistic for actual data set
  fit.sim <- sum(chi2.sim[,])        # Fit statistic for a fitting model
  c.hat <- fit.actual / fit.sim      # c-hat estimate
  bpv <- step(fit.sim-fit.actual)    # Bayesian p-value
}
",fill = TRUE)
sink()

# Do other preps and run model with JAGS
inits <- function(){list(N = apply(data$C, 1, max)+1)}
params <- c("lambda", "p", "fit.actual", "fit.sim", "c.hat", "bpv")
ni <- 2500   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3
fm <- jags(win.data, inits, params, "model.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb)
print(fm, dig = 3)

ppc.plot(fm)               # Produces Fig. 6–7

