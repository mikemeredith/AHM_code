#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7. Modeling abundance using multinomial N-mixture models
# =========================================================================

library(unmarked)
library(R2WinBUGS)
bd <- "C:/WinBUGS14"

# ~~~~~ this section uses data from section 7.5 ~~~~~~~~~~~~~
data(ovendata)
ovenFrame <- unmarkedFrameMPois(y = ovendata.list$data,
    siteCovs = as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 7.6 Bayesian Analysis in BUGS using the Conditional Multinomial (3-part model)
# ==============================================================================

# ----- pseudo-code, can't be run as-is ~~~~~~~~~~~~~~~
# Set-up data with a missing value for element "not captured"
# y[i,] <- c(y[i,],NA)

# Then, in BUGS, do this:
# y[i,] ~ dmulti(probs[i,], N[i])
# N[i] ~ dpois(lambda[i])
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Harvest the data and bundle it up for sending to BUGS
y <- as.matrix(getY(ovenFrame))
ncap <- apply(y, 1, sum)   # number of individuals removed per point
data <- list(y = y, M = nrow(y), n = ncap, X=as.matrix(siteCovs(ovenFrame)))
str(data)                  # Good practice to always inspect your BUGS data

# Write BUGS model
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- logit(p0)
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)

  for(i in 1:M){ # Loop over sites
    # Conditional multinomial cell probabilities
    pi[i,1] <- p[i]
    pi[i,2] <- p[i]*(1-p[i])
    pi[i,3] <- p[i]*(1-p[i])*(1-p[i])
    pi[i,4] <- p[i]*(1-p[i])*(1-p[i])*(1-p[i])
    pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3] + pi[i,4])
    pcap[i] <- 1 - pi0[i]
    for(j in 1:4){
      pic[i,j] <- pi[i,j] / pcap[i]
    }

    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]

    # Model specification, three parts:
    y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional
                                        #    cell probabilities
    n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the
                                        #    observed sample size
    N[i] ~ dpois(lambda[i])             # component 3 is the process model

    # log-linear model for abundance: UFC + TRBA + UFC:TRBA
    log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
  }
}
",fill=TRUE, file="model.txt")

# Initial values
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0 = runif(1), N = ncap+2)
}

# Parameters monitored

parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3")

# MCMC settings
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1

# Experience the power of BUGS and print posterior summary
# ~~~~~ WinBUGS fails with "Rejection1" appearing in the status bar ~~~~~~~~~~~~~
# ~~~~~ The JAGS code below works fine.
# library(R2WinBUGS)
# bd <- "c:/Program Files/WinBUGS14/" # ~~~~~ doesn't work
# out <- bugs(data, inits, parameters, "model.txt", n.thin = nt,
       # n.chains = nc, n.burnin = nb, n.iter = ni, debug = TRUE,
       # bugs.directory = bd)
# print(out, 3)
# To make trace plots using CODA functionality
# plot(as.mcmc.list(out), ask = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(jagsUI)
out <- jags(data, inits, parameters, "model.txt", n.thin=nt,
   n.chains = nc, n.burnin = nb, n.iter = ni)
print(out, 3)


# 7.6.1 Goodness-of-fit using Bayesian p-values
# ------------------------------------------------------------------------

# ~~~~~ The code supplied needs to be inserted in the BUGS model ~~~~~~~~~
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- logit(p0)
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)

  for(i in 1:M){ # Loop over sites
    # Conditional multinomial cell probabilities
    pi[i,1] <- p[i]
    pi[i,2] <- p[i]*(1-p[i])
    pi[i,3] <- p[i]*(1-p[i])*(1-p[i])
    pi[i,4] <- p[i]*(1-p[i])*(1-p[i])*(1-p[i])
    pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3] + pi[i,4])
    pcap[i] <- 1 - pi0[i]
    for(j in 1:4){
      pic[i,j] <- pi[i,j] / pcap[i]
    }

    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]

    # Model specification, three parts:
    y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional
                                        #    cell probabilities
    n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the
                                        #    observed sample size
    N[i] ~ dpois(lambda[i])             # component 3 is the process model

    # log-linear model for abundance: UFC + TRBA + UFC:TRBA
    log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
  }
  # ~~~ inserted code:
  for(i in 1:M){
    n.pred[i] ~ dbin(pcap[i], N[i])
    y.pred[i, 1:4] ~ dmulti(pic[i, 1:4], n[i]) #note this is cond'l on ncap[i]
    for(k in 1:4) {
      e1[i,k] <- pic[i,k]*n[i]
      resid1[i,k] <- pow(pow(y[i,k],0.5)-pow(e1[i,k],0.5),2)
      resid1.pred[i,k] <- pow(pow(y.pred[i,k],0.5)-pow(e1[i,k],0.5),2)
    }
    e2[i] <- pcap[i]*lambda[i]
    resid2[i] <- pow(pow(n[i],0.5)-pow(e2[i],0.5),2)
    resid2.pred[i] <- pow(pow(n.pred[i],0.5)-pow(e2[i],0.5),2)
  }
  fit1.data <- sum(resid1[,])       # fit statistic for observed data y
  fit1.pred <- sum(resid1.pred[,])

  fit2.data <- sum(resid2[])       # fit statistic for new data n
  fit2.pred <- sum(resid2.pred[])
}
",fill=TRUE, file="model_gof.txt")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


parameters <- c("N", "p0", "beta0", "beta1", "beta2", "beta3",
   "fit1.data", "fit1.pred", "fit2.data", "fit2.pred")

# ~~~~~ WinBUGS sometimes fails with "Rejection1" appearing in the status bar ~~~~~~~~~~~~~
# out <- bugs (data, inits, parameters, "model_gof.txt",n.thin=nt, n.chains=nc,
       # n.burnin=nb, n.iter=ni, debug=TRUE, bugs.dir = bd)
# ~~~~ use JAGS instead
out <- jags (data, inits, parameters, "model_gof.txt",n.thin=nt, n.chains=nc,
       n.burnin=nb, n.iter=ni, parallel=TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mean(out$sims.list$fit1.pred > out$sims.list$fit1.data)

mean(out$sims.list$fit2.pred > out$sims.list$fit2.data)


# 7.6.2 Model selection in BUGS
# ------------------------------------------------------------------------
# Write BUGS model
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.1)
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  beta3 ~ dnorm(0, 0.1)
  w1 ~ dbern(0.5)
  w2 ~ dbern(0.5)
  w3 ~ dbern(0.5)

  for(i in 1:M){
    # Conditional multinomial cell probabilities
    pi[i,1] <- p[i]
    pi[i,2] <- p[i] * (1-p[i])
    pi[i,3] <- p[i] * (1-p[i]) * (1-p[i])
    pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i])
    pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3]+ pi[i,4])
    pcap[i] <- 1 - pi0[i]
    for(j in 1:4){
      pic[i,j] <- pi[i,j] / pcap[i]
    }

    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]

    # Model specification, 3 parts:
    y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional
                                        #    cell probabilities
    n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the
                                        #    observed sample size
    N[i] ~ dpois(lambda[i])             # component 3 is the process model

    # log-linear model for abundance: effects of UFC and TRBA, with weights
    log(lambda[i])<- beta0 + w1*beta1*X[i,1] + w2*beta2*X[i,2] +
        w1*w2*w3*beta3*X[i,2]*X[i,1]
  }
}
",fill=TRUE,file="model.txt")


# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2",
   "beta3", "w1", "w2", "w3")

# Initial values
set.seed(2015)
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0 = runif(1), N = ncap+2,
         w1=1, w2=1, w3=1)
}
# Call WinBUGS from R and summarize marginal posteriors
out <- bugs(data, inits, parameters, "model.txt",
  n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni,
  # debug = TRUE, bugs.dir=bd)
  debug = FALSE, bugs.dir=bd)  # ~~~~ for autotesting
print(out, digits = 3)


# Extract the indicator variables for model selection
w1 <- out$sims.list$w1
w2 <- out$sims.list$w2

# Create a new w3 variable which takes on the value 1 if the
#     interaction is in the model, also requires w1 = 1 AND w2=1
w3 <- out$sims.list$w3 * w1 * w2

# Combine into a model indicator string and tabulate posterior frequencies
mod <- paste(w1, w2, w3, sep = "")
table(mod)


# 7.6.3 Poisson formulation of the multinomial mixture model
# ------------------------------------------------------------------------
# Specify model in BUGS language
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- logit(p0)
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)

  for(i in 1:M){
    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]
    # log-linear model for abundance: UFC + TRBA + UFC:TRBA
    log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]

    # Poisson parameter = multinomial cellprobs x expected abundance
    pi[i,1] <- p[i] * lambda[i]
    pi[i,2] <- p[i] * (1-p[i]) * lambda[i]
    pi[i,3] <- p[i] * (1-p[i]) * (1-p[i]) * lambda[i]
    pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i]) * lambda[i]

    for(j in 1:4){
      y[i,j] ~ dpois(pi[i,j])
    }
    # Generate predictions of N[i]
    N[i] ~ dpois(lambda[i])
  }
}
",fill=TRUE,file="modelP.txt")

# Bundle up the data and inits
data <- list(y = y, M = nrow(y), X = as.matrix(siteCovs(ovenFrame)))
inits <- function(){
  list (p0 = runif(1), alpha1=runif(1), beta0=runif(1), beta1=runif(1), beta2=runif(1), beta3=runif(1))
}

# Define parameters to save and MCMC settings
parameters <- c("p0", "alpha1", "beta0", "beta1", "beta2", "beta3", "N")
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1


# Call WinBUGS from R and summarize marginal posteriors
out <- bugs(data, inits, parameters, "modelP.txt",
  n.thin=nt, n.chains = nc, n.burnin = nb, n.iter = ni,
  # debug = TRUE, bugs.dir=bd)
  debug = FALSE, bugs.dir=bd) # ~~~~ for autotesting
print(out, 3)


# Specify model in BUGS language
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- logit(p0)
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)

  tau ~ dgamma(0.1,0.1)  # Excess-Poisson variation (precision)
  sigma <- sqrt(1 / tau)

  for(i in 1:M){
    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]
    # Normal random effects
    eta[i] ~ dnorm(0, tau)  # 'residuals' for extra-Poisson noise
    # log-linear model for abundance: UFC + TRBA + UFC:TRBA + eta
    log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1] + eta[i]
        # note 'extra-residual' for overdispersion

    # Poisson parameter = multinomial cellprobs x expected abundance
    pi[i,1] <- p[i] * lambda[i]
    pi[i,2] <- p[i] * (1-p[i]) * lambda[i]
    pi[i,3] <- p[i] * (1-p[i]) * (1-p[i]) * lambda[i]
    pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i]) * lambda[i]

    for(j in 1:4){
      y[i,j] ~ dpois(pi[i,j])
    }
    # Generate predictions of N[i]
    N[i] ~ dpois(lambda[i])
  }
}
",fill=TRUE,file="modelP.txt")

# Inits
inits <- function(){
  list (p0 = runif(1), alpha1=runif(1), beta0=runif(1), beta1=runif(1),
         beta2=runif(1), beta3=runif(1), tau = 1)  }

# Define parameters to save and MCMC settings
parameters <- c("p0", "alpha1", "beta0", "beta1", "beta2", "beta3", "sigma")
nc <- 3   ;   ni <- 32000   ;   nb <- 2000   ;   nt <- 1

# Call WinBUGS from R and summarize marginal posteriors
out <- bugs(data, inits, parameters, "modelP.txt",
  n.thin=nt, n.chains = nc, n.burnin = nb, n.iter = ni,
  # debug = TRUE, bugs.dir = bd)
  debug = FALSE, bugs.dir = bd) # ~~~~ for autotesting
print(out, digits=3)

plot(density(out$sims.list$sigma), frame = FALSE, main = "") # Fig. 7-4

