#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 5. Fitting models using the Bayesian modeling software BUGS and JAGS
# =========================================================================

# ~~~~~ this section requires the following code from section 5.3 ~~~~~~~~~~
library(AHMbook)
# Generate data with data.fn from chapter 4
set.seed(24)
data <- data.fn(show.plot=FALSE)
attach(data)
# Summarize data by taking mean at each site and plot
Cmean <- apply(C, 1, mean)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 5.8 Fitting a model with non-standard likelihood using the zeros or the ones tricks
# ===================================================================================

# Package the data needed in a bundle
win.data <- list(Cmean1 = Cmean, Cmean2 = Cmean, zeros = rep(0, M),
    ones = rep(1, M), M = length(Cmean), elev = elev, forest = forest) # note 2 copies of response

# Write text file with model description in BUGS language
cat(file = "multiple_linear_regression_model.txt",
"model {

  # Priors
  for(k in 1:3){ # Loop over three ways to specify likelihood
    alpha0[k] ~ dnorm(0, 1.0E-06)           # Prior for intercept
    alpha1[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
    alpha2[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
    alpha3[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
    sd[k] ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
  }
  var1 <- pow(sd[1], 2)                      # Variance in zeros trick
  var2 <- pow(sd[2], 2)                      # Variance in ones trick
  tau <- pow(sd[3], -2)                      # Precision tau = 1/(sd^2)

  C1 <- 10000 # zeros trick: make large enough to ensure lam >= 0
  C2 <- 10000 # ones trick: make large enough to ensure p <= 1
  pi <- 3.1415926

  # Three variants of specification of the likelihood
  for (i in 1:M){
    # 'Zeros trick' for normal likelihood
    zeros[i] ~ dpois(phi[i])  # likelihood contribution is exp(-phi)
    #   negLL[i] <- log(sd[1]) + 0.5 * pow((Cmean1[i] - mu1[i]) / sd[1],2 )
    negLL[i] <- -log(sqrt(1/(2*pi*var1))) + pow(Cmean1[i]-mu1[i],2)/(2*var1)
    phi[i] <- negLL[i] + C1
    mu1[i] <- alpha0[1] + alpha1[1]*elev[i] + alpha2[1]*forest[i] + alpha3[1]*elev[i]*forest[i]

    # 'Ones trick' for normal likelihood
    ones[i] ~ dbern(p[i])  # likelihood contribution is p directly
    L[i] <- sqrt(1/(2*pi*var2)) * exp(-pow(Cmean1[i]-mu2[i],2)/(2*var2))
    p[i] <- L[i] / C2
    mu2[i] <- alpha0[2] + alpha1[2]*elev[i] + alpha2[2]*forest[i] + alpha3[2]*elev[i]*forest[i]

    # Standard distribution function for the normal
    Cmean2[i] ~ dnorm(mu3[i], tau)
    mu3[i] <- alpha0[3] + alpha1[3]*elev[i] + alpha2[3]*forest[i] + alpha3[3]*elev[i]*forest[i]
  }
}"
)

# Initial values
inits <- function() list(alpha0 = rnorm(3, 0, 10), alpha1 = rnorm(3,0,10),
    alpha2 = rnorm(3,0,10), alpha3 = rnorm(3,0,10))

# Parameters monitored (i.e., for which estimates are saved)
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd")

# MCMC settings
ni <- 1200   ;   nt <- 1   ;   nb <- 200   ;  nc <- 3    # For JAGS

# Call JAGS
library(jagsUI)
outX <- jags(win.data, inits, params, "multiple_linear_regression_model.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(outX)


# Define negative log-likelihood function
neglogLike <- function(param) {
  alpha0 = param[1]
  alpha1 = param[2]
  alpha2 = param[3]
  alpha3 = param[4]
  sigma = exp(param[5])    # Estimate sigma on log-scale
  mu = alpha0 + alpha1*elev + alpha2*forest + alpha3*elev*forest
  #   -sum(dnorm(Cmean, mean=mu, sd=sigma, log=TRUE))  # cheap quick way
  sum(-log(sqrt(1/(2*3.1415926*sigma^2))) + (Cmean-mu)^2/(2*sigma^2))
}

# Find parameter values that minimize function value
(fit <- optim(par = rep(0, 5), fn = neglogLike, method = "BFGS"))

exp(fit$par[5])            # Backtransform to get sigma

