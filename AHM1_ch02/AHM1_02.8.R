#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 2. What are hierarchical models and how do we analyze them?
# =========================================================================

# 2.8 Assessment of model fit
# ===========================

# 2.8.1 Parametric bootstrapping example
# ------------------------------------------------------------------------

sim.data <- function(beta0 = -3, beta1 = 2, p = 0.6, x=NULL){
  # Function allows input of covariate "x", or simulates new

  M <- 100
  if(is.null(x))
     vegHt <- runif(M, 1, 3) # uniform from 1 to 3

  # Suppose that occupancy probability increases with vegHt
  # The relationship is described (default) by an intercept of -3 and
  #    a slope parameter of 2 on the logit scale
  # plogis is the inverse-logit (constrains us back to the [0-1] scale)
  psi <- plogis(beta0 + beta1*vegHt)

  # Now we simulated true presence/absence for 100 sites
  z <- rbinom(M, 1, psi)

  # Now generate observations
  J <- 3 # sample each site 3 times
  y <- rbinom(M,J,p*z)

  list(y=y, J=J, vegHt=vegHt)
}

# This is the negative log-likelihood based on the marginal distribution
# of y. It is the pmf of a zero-inflated binomial random variable.
#
negLogLikeocc <- function(beta, y, x, J) {
   beta0 <- beta[1]
   beta1 <- beta[2]
   p<- plogis(beta[3])
   psi <- plogis(beta0 + beta1*x)
   marg.likelihood <- dbinom(y, J, p) * psi + ifelse(y==0, 1, 0) * (1-psi)
   return(-sum(log(marg.likelihood)))
}

data <- sim.data()        # Generate a data set

# Let's minimize it
starting.values <- c(beta0=0, beta1=0, logitp=0)
opt.out <- optim(starting.values, negLogLikeocc, y=data$y, x=data$vegHt,
    J=data$J, hessian=TRUE)
(mles <- opt.out$par)

# Make a table with estimates, SEs, and 95% CI
mle.table <- data.frame(Est=mles,
                        SE = sqrt(diag(solve(opt.out$hessian))))
mle.table$lower <- mle.table$Est - 1.96*mle.table$SE
mle.table$upper <- mle.table$Est + 1.96*mle.table$SE
mle.table


# Define a fit statistic
fitstat <- function(y, Ey){
  sum((sqrt(y) - sqrt(Ey)))
}
# Compute it for the observed data
# ~~~ 3 lines of code added to ensure we are using output from sim.data(), see Errata 2021-10-09
y <- data$y
J <- data$J
vegHt <- data$vegHt
T.obs <- fitstat(y, J*plogis(mles[1] + mles[2]*vegHt)*plogis(mles[3]))

# Get bootstrap distribution of fit statistic
T.boot <- rep(NA, 100)
for(i in 1:100){
  # Simulate a new data set and extract the elements. Note we use
  # the previously simulated "vegHt" covariate
  data <- sim.data(beta0=mles[1],beta1=mles[2],p=plogis(mles[3]),x=vegHt)
  # Next we fit the model
  starting.values <- c(0,0,0)
  opt.out <- optim(starting.values, negLogLikeocc, y=data$y, x= data$vegHt, J=data$J, hessian=TRUE)
  (parms <- opt.out$par)
  # Obtain the fit statistic
  T.boot[i]<- fitstat(y, J*plogis(parms[1] + parms[2]*vegHt)*plogis(parms[3]) )
}

(T.obs)

summary(T.boot)


# 2.8.2 Bayesian p-value (no code)
# ------------------------------------------------------------------------

