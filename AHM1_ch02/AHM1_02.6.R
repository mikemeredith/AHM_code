#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 2. What are hierarchical models and how do we analyze them?
# =========================================================================

# ~~~~~~~ This requires results from section 2.4 ~~~~~~~~~~~~~~~~~~~~
source("AHM1_02.4.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.6 Basic Markov chain Monte Carlo (MCMC)
# =========================================

# 2.6.1 Metropolis-Hastings algorithm (no code)

# 2.6.2 Illustration: using MH for a binomial model
# ------------------------------------------------------------------------
# Simulate data
set.seed(2016)
y <- rbinom(2, size=10, p = 0.5)

# Define the joint distribution (= likelihood) which we will maximize
jointdis <- function(data, K, p){
   prod(dbinom(data, size=K, p=p))
}

# Posterior is proportional to likelihood times prior
posterior <- function(p, data, K, a, b){
   prod(dbinom(data, size=K, p=p)) * dbeta(p, a, b)
}

# Do 100,000 MCMC iterations using Metropolis algorithm
# Assume vague prior which is beta(1,1) = Unif(0,1)
mcmc.iters <- 100000
out <- rep(NA, mcmc.iters)

# Starting value
p <- 0.2

# Begin the MCMC loop
for(i in 1:mcmc.iters){

   # Use a uniform candidate generator (not efficient)
   p.cand <- runif(1, 0, 1)

   # Alternative: random walk proposal
   # p.cand <- rnorm(1, p, .05)  # Need to reject if > 1 or < 0
   # if(p.cand < 0 | p.cand > 1 ) next

   r <- posterior(p=p.cand, y, K=10, a=1, b=1) / posterior(p=p, y, K=10, a=1, b=1)
   # Generate a uniform r.v. and compare with "r", this imposes the
   #    correct probability of acceptance
   if(runif(1) < r)
    p <- p.cand

   # Save the current value of p
   out[i] <- p
}


mean(out)

sd(out)

quantile(out, c(0.025, 0.975))

# Evaluate likelihood for a grid of values of p
p.grid <- seq(0.1, 0.9, , 200)
likelihood <- rep(NA, 200)

for(i in 1:200){
   likelihood[i] <- jointdis(y, K=10, p=p.grid[i])
}

op <- par(mfrow=c(2,1), mar = c(5,5,3,2))
plot(p.grid, likelihood, xlab="", ylab="Likelihood", xlim=c(0,1), ty = "l", main = "Likelihood function")
p.hat <- p.grid[likelihood == max(likelihood)]
abline(v = p.hat)
text(p.hat, 0.005, paste("MLE = ", round(p.hat, 3), sep= ""))

plot(density(out), xlim=c(0,1), main = "Posterior distribution", xlab = "p", ylab = "Posterior")
p.mean <- mean(out)
abline(v = p.mean)
text(p.mean, 0.5, paste("Post. mean = ", round(p.mean, 3),sep=" "))
par(op)


# 2.6.3 Metropolis algorithm for multi-parameter models
# ----------------------------------------------------------------------------------------
log.posterior <- function(beta0, beta1, z, vegHt){
# Note: "z" and "vegHt" must be input
   loglike <- -1 * negLogLike(c(beta0, beta1), z, vegHt)
   logprior <- dnorm(c(beta0, beta1), 0, 10, log=TRUE)
   return(loglike + logprior[1] + logprior[2])
}

niter <- 50000
out <- matrix(NA, niter, 2, dimnames = list(NULL, c("beta0", "beta1")))

# Initialize parameters
beta0 <- rnorm(1)
beta1 <- rnorm(1)

# Current value of the log(posterior)
logpost.curr <- log.posterior(beta0, beta1, z, vegHt)

# Run MCMC algorithm
for(i in 1:niter){
   if(i %% 1000 == 0)                     # report progress
      cat("iter", i, "\n")
      # Update intercept (beta0)
      # Propose candidate values of beta
      # If the proposal was not symmetric, would be Metrop-*Hastings*
      beta0.cand <- rnorm(1, beta0, 0.3) # 0.3 is tuning parameter
      # Evaluate the log(posterior)
      logpost.cand <- log.posterior(beta0.cand, beta1, z, vegHt)
      # Compute Metropolis acceptance probability, r
      r <- exp(logpost.cand - logpost.curr)
      # Keep candidate if it meets criterion (u < r)
      if(runif(1) < r){
         beta0 <- beta0.cand
         logpost.curr <- logpost.cand
      }

      # Update slope (beta1)
      beta1.cand <- rnorm(1, beta1, 0.3) # 0.3 is tuning parameter
      # Evaluate the log(posterior)
      logpost.cand <- log.posterior(beta0, beta1.cand, z, vegHt)

       # Compute Metropolis acceptance probability
       r <- exp(logpost.cand - logpost.curr)
       # Keep candidate if it meets criterion (u < r)
       if(runif(1) < r){
          beta1 <- beta1.cand
          logpost.curr <- logpost.cand
       }
   out[i,] <- c(beta0, beta1) # Save samples for iteration
}

# Plot
op <- par("mfrow", oma=c(0,0,0,0),mar=c(5,4,1,1))
layout(rbind(c(1,1),
             c(2,2),
             c(3,4)), respect = TRUE) # <- play with these settings

plot(out[,1], type="l", xlab="Iteration", ylab="beta0")
plot(out[,2], type="l", xlab="Iteration", ylab="beta1")
plot(density(out[,1]), xlab="beta0", main="")
plot(density(out[,2]), xlab="beta1", main="")
par(op)

# 2.6.4 Why we need to use MCMC for logistic regression (no code)
# 2.6.5 Gibbs sampling (no code)
# 2.6.6 Convergence and mixing of Markov chains (no code)
# 2.6.6.1 Slow mixing and thinning of Markov chains (no code)
# 2.6.6.2 Effective sample size and Monte Carlo error (no code)

# 2.6.7 Bayesian analysis of hierarchical models
# ------------------------------------------------------------------------

# The occupancy model
# '''''''''''''''''''

# Simulate the data set
set.seed(2014)
M <- 100                        # number of sites
vegHt <- runif(M, 1, 3)         # uniform from 1 to 3
psi <- plogis(-3 + 2*vegHt)     # occupancy probability
z <- rbinom(M, 1, psi)          # realised presence/absence
p <- 0.6                        # detection probability
J <- 3                          # sample each site 3 times
y <-rbinom(M, J, p*z)           # observed detection frequency

# Number of MCMC iterations to to
niter <- 50000

# Matrix to hold the simulated values
out <- matrix(NA, niter, 3, dimnames = list(NULL, c("beta0", "beta1", "p")))

# Initialize parameters, likelihood, and priors
starting.values <- c(beta0=0, beta1=0)
beta0 <- starting.values[1]
beta1 <- starting.values[2]
z <-ifelse(y>0, 1, 0)
p <- 0.2

# NOTE: using logistic reg. likelihood function here!
loglike <- -1*negLogLike(c(beta0, beta1), z, vegHt)
logprior <- dnorm(c(beta0, beta1), 0, 10, log=TRUE)

# Run MCMC algorithm
for(i in 1:niter) {
   if(i %% 1000 ==0 ) # report progress
   cat("iter", i, "\n")

   # PART 1 of algorithm -- same as before
   # Update intercept (beta0)
   # propose candidate values of beta
   beta0.cand <- rnorm(1, beta0, 0.3) # 0.3 is tuning parameter
   # evaluate likelihood and priors for candidates
   loglike.cand <- -1*negLogLike(c(beta0.cand, beta1), z, vegHt)
   logprior.cand <- dnorm(beta0.cand, 0, 10, log=TRUE)
   # Compute Metropolis acceptance probability (r)
   r <- exp((loglike.cand+logprior.cand) - (loglike + logprior[1]))
   # Keep candidate if it meets the criterion
   if(runif(1) < r){
      beta0 <- beta0.cand
      loglike <- loglike.cand
      logprior[1] <- logprior.cand
   }
   # Update slope (beta1)
   beta1.cand <- rnorm(1, beta1, 0.3) # 0.3 is tuning parameter
   # evaluate likelihood and priors for candidates
   loglike.cand <- -1*negLogLike(c(beta0,beta1.cand), z, vegHt)
   logprior.cand <- dnorm(beta1.cand, 0, 10, log=TRUE)
   # Compute Metropolis acceptance probability r
   r <- exp((loglike.cand+logprior.cand) - (loglike + logprior[2]))
   # Keep the candidates if they meet the criterion
   if(runif(1) < r) {
      beta1 <- beta1.cand
      loglike <- loglike.cand
      logprior[2] <- logprior.cand
   }

   # Part 2 of the algorithm
   # update z. Note we only need to update z if y=0.
   # The full conditional has known form

   psi <- plogis(beta0 + beta1 * vegHt)
   psi.cond <- dbinom(0,J,p) * psi /(dbinom(0, J, p) * psi + (1-psi))
   z[y==0] <- rbinom(sum(y==0), 1, psi.cond[y==0])
   loglike <- -1 * negLogLike(c(beta0, beta1), z, vegHt)

   # Part 3: update p
  ## The commented code will update p using Metropolis
  ## loglike.p <- sum(log(dbinom(y[z==1],J,p)))
  ## p.cand <- runif(1, 0, 1)
  ## loglike.p.cand <- sum(log(dbinom(y[z==1], J, p.cand)))
  ## if(runif(1) < exp(loglike.p.cand-loglike.p))
  ##    p<-p.cand
  ## This bit draws p directly from its full conditional
   p <- rbeta(1, 1+ sum(y), sum(z)*J +1 - sum(y) )

   # Save MCMC samples
   out[i,] <- c(beta0,beta1,p)
}

# Plot bivariate representation of joint posterior
pairs(out)

# Trace/history plots for each parameter (Fig. 2.8)
op <- par(mfrow=c(2,1))
plot(out[,1], type="l", xlab="Iteration", ylab="beta0")
abline(h = mean(out[,1]), col = "blue", lwd = 2)
abline(h = -3, col = "red", lwd = 2)
plot(out[,2], type="l", xlab="Iteration", ylab="beta1")
abline(h = mean(out[,2]), col = "blue", lwd = 2)
abline(h = 2, col = "red", lwd = 2)
par(op)
