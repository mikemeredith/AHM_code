#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle

# Chapter 5. Fitting models using the Bayesian modeling software BUGS and JAGS
# =========================================================================

library(AHMbook)
library(jagsUI)

# ~~~~~ this section requires the following code from section 5.3 ~~~~~~~~~~
# Generate data with data.fn from chapter 4
set.seed(24)
data <- data.fn(show.plot=FALSE)
attach(data)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.12 Moment-matching in a binomial GLM to accommodate underdispersion
# =====================================================================


# Bundle data
win.data <- list(y = data$C, M = nrow(data$C), J = ncol(data$C), elev = elev, N = 32)

# Specify model in BUGS language
cat(file = "squeezed_count_GLM.txt","
model {

  # Priors
  alpha ~ dnorm(0, 1.0E-06)
  beta ~ dnorm(0, 1.0E-06)

  # Likelihood
  for (i in 1:M){
    mu[i] <- alpha + beta * elev[i] # linear model for expected response
    logit(p[i]) <- logit(mu[i] / N) # express param as function of first moment
    for(j in 1:J){
      y[i,j] ~ dbin(p[i], N)
    }
  }
}
")

# Initial values
inits <- function() list(alpha = 1.7, beta = -1.2)          # works always
# inits <- function() list(alpha = runif(1), beta = runif(1)) # works sometimes

# Parameters monitored
params <- c("alpha", "beta", "mu")

# MCMC settings
ni <- 3000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call JAGS from R (ART <1 min)
library(jagsUI)
out7 <- jags(win.data, inits, params, "squeezed_count_GLM.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
# par(mfrow = c(4,2))  #  ~~~ replace with 'layout' argument
traceplot(out7, layout=c(4,2))
print(out7, 3)

