#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 2. What are hierarchical models and how do we analyze them?
# =========================================================================

# 2.4 Classical inference based on likelihood
# ===========================================

# 2.4.1 The frequentist interpretation (no code)
# 2.4.2 Properties of MLEs (no code)
# 2.4.3 Delta approximation (no code)

# 2.4.4 Example: Classical inference for logistic regression
# ------------------------------------------------------------------------
# Simulate a covariate called vegHt for 100 sites
set.seed(2014)  # Set seed so we all get the same values of vegHt
M <- 100        # Number of sites surveyed
vegHt <- runif(M, 1, 3) # uniform from 1 to 3

# Suppose that occupancy probability increases with vegHt
# The relationship is described by an intercept of -3 and
#    a slope parameter of 2 on the logit scale
beta0 <- -3
beta1 <- 2
psi <- plogis(beta0 + beta1*vegHt) # apply inverse logit

# Now we go to 100 sites and observe presence or absence
z <- rbinom(M, 1, psi)

# Definition of negative log-likelihood.
negLogLike <- function(beta, y, x) {
    beta0 <- beta[1]
    beta1 <- beta[2]
    psi <- plogis(beta0 + beta1*x)
    likelihood <- psi^y * (1-psi)^(1-y) # same as next line:
#   likelihood <- dbinom(y, 1, psi)
    return(-sum(log(likelihood)))
}

# Look at (negative) log-likelihood for 2 parameter sets
negLogLike(c(0,0), y=z, x=vegHt)
negLogLike(c(-3,2), y=z, x=vegHt) # Lower is better!

# Let's minimize it formally by function minimisation
starting.values <- c(beta0=0, beta1=0)
opt.out <- optim(starting.values, negLogLike, y=z, x=vegHt, hessian=TRUE)
(mles <- opt.out$par)     # MLEs are pretty close to truth

# Alternative 1: Brute-force grid search for MLEs
mat <- as.matrix(expand.grid(seq(-10,10,0.1), seq(-10,10,0.1)))
                                               # above: Can vary resolution
nll <- array(NA, dim = nrow(mat))
for (i in 1:nrow(mat)){
  nll[i] <- negLogLike(mat[i,], y = z, x = vegHt)
}
which(nll == min(nll))
mat[which(nll == min(nll)),]

# Produce a likelihood surface, shown in Fig. 2-2.
library(raster)
r <- rasterFromXYZ(data.frame(x = mat[,1], y = mat[,2], z = nll))
mapPalette <- colorRampPalette(rev(c("grey", "yellow", "red")))
plot(r, col = mapPalette(100), main = "Negative log-likelihood",
       xlab = "Intercept (beta0)", ylab = "Slope (beta1)")
contour(r, add = TRUE, levels = seq(50, 2000, 100))

# Alternative 2: Use canned R function glm as a shortcut
(fm <- glm(z ~ vegHt, family = binomial)$coef)

# Add 3 sets of MLEs into plot
# 1. Add MLE from function minimisation
points(mles[1], mles[2], pch = 1, lwd = 2)
abline(mles[2],0)  # Put a line through the Slope value
lines(c(mles[1],mles[1]),c(-10,10))
# 2. Add MLE from grid search
points(mat[which(nll == min(nll)),1], mat[which(nll == min(nll)),2],
       pch = 1, lwd = 2)
# 3. Add MLE from glm function
points(fm[1], fm[2], pch = 1, lwd = 2)

# Note they are essentially all the same (as they should be)

Vc <- solve(opt.out$hessian)         # Get variance-covariance matrix
ASE <- sqrt(diag(Vc))                # Extract asymptotic SEs
print(ASE)

# Compare to SEs reported by glm() function (output thinned)
summary(glm(z ~ vegHt, family = binomial))

# Make a table with estimates, SEs, and 95% CI
mle.table <- data.frame(Est=mles,
                        ASE = sqrt(diag(solve(opt.out$hessian))))
mle.table$lower <- mle.table$Est - 1.96*mle.table$ASE
mle.table$upper <- mle.table$Est + 1.96*mle.table$ASE
mle.table

# Plot the actual and estimated response curves
plot(vegHt, z, xlab="Vegetation height", ylab="Occurrence probability")
plot(function(x) plogis(beta0 + beta1 * x), 1.1, 3, add=TRUE, lwd=2)
plot(function(x) plogis(mles[1] + mles[2] * x), 1.1, 3, add=TRUE,
     lwd=2, col="blue")
legend(1.1, 0.9, c("Actual", "Estimate"), col=c("black", "blue"), lty=1,
       lwd=2)


# 2.4.5 Bootstrapping
# ------------------------------------------------------------------------
nboot <- 1000   # Obtain 1000 bootstrap samples
boot.out <- matrix(NA, nrow=nboot, ncol=3)
dimnames(boot.out) <- list(NULL, c("beta0", "beta1", "psi.bar"))

for(i in 1:1000){
   # Simulate data
   psi <- plogis(mles[1] + mles[2] * vegHt)
   z <- rbinom(M, 1, psi)

   # Fit model
   tmp <- optim(mles, negLogLike, y=z, x=vegHt, hessian=TRUE)$par
   psi.mean <- plogis(tmp[1] + tmp[2] * mean(vegHt))
   boot.out[i,] <- c(tmp, psi.mean)
}

SE.boot <- sqrt(apply(boot.out, 2, var))  # Get bootstrap SE
names(SE.boot) <- c("beta0", "beta1", "psi.bar")

# 95% bootstrapped confidence intervals
apply(boot.out,2,quantile,c(0.025,0.975))

# Boostrap SEs
SE.boot

# Compare these with the ASEs for regression parameters
mle.table


# 2.4.6 Likelihood analysis of hierarchical models (no code)
# ------------------------------------------------------------------------

# 2.4.6.1 Discrete random variable
# ------------------------------------------------------------------------
set.seed(2014)
M <- 100                             # number of sites
vegHt <- runif(M, 1, 3)              # uniform from 1 to 3
psi <- plogis(beta0 + beta1 * vegHt) # occupancy probability
z <- rbinom(M, 1, psi)               # realised presence/absence
p <- 0.6                             # detection probability
J <- 3                               # sample each site 3 times
y <-rbinom(M, J, p*z)                # observed detection frequency

# Define negative log-likelihood.
negLogLikeocc <- function(beta, y, x, J) {
    beta0 <- beta[1]
    beta1 <- beta[2]
    p <- plogis(beta[3])
    psi <- plogis(beta0 + beta1*x)
    marg.likelihood <- dbinom(y, J, p)*psi + ifelse(y==0,1,0)*(1-psi)
    return(-sum(log(marg.likelihood)))
}

starting.values <- c(beta0=0, beta1=0,logitp=0)
(opt.out <- optim(starting.values, negLogLikeocc, y=y, x=vegHt,J=J,
                 hessian=TRUE))

sqrt(diag(solve(opt.out$hessian)))


# 2.4.6.2 A continuous latent variable
# ------------------------------------------------------------------------
# ~~~~ Following code not executable as is: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# marg <- rep(NA, J+1)
# for(j in 0:J){
# marg[j] <- integrate(
   # function(x){
      # dbinom(j, J, plogis(x)) * dnorm(x, mu, sigma)},
      # lower=-Inf,upper=Inf)$value
   # }
# }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# nx = encounter frequencies, number inds. encountered 1, 2, ..., 14 times
nx <- c(34, 16, 10, 4, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0)
nind <- sum(nx)      # Number of individuals observed
J <- 14              # Number of sample occasions

# Model Mh likelihood
Mhlik <- function(parms){
   mu <- parms[1]
   sigma <- exp(parms[2])
   # n0 = number of UNobserved individuals: N = nind + n0
   n0 <- exp(parms[3])

  # Compute the marginal probabilities for each possible value j=0,..,14
   marg <- rep(NA,J+1)
   for(j in 0:J){
      marg[j+1] <- integrate(
      function(x){dbinom(j, J, plogis(x)) * dnorm(x, mu, sigma)},
       lower=-Inf,upper=Inf)$value
   }

  # The negative log likelihood involves combinatorial terms computed
  # using lgamma()
  -1*(lgamma(n0 + nind + 1) - lgamma(n0 + 1) + sum(c(n0, nx) * log(marg)))
}
(tmp <- nlm(Mhlik, c(-1, 0, log(10)), hessian=TRUE))


(SE <- sqrt( (exp(tmp$estimate[3])^2)* diag(solve(tmp$hessian))[3] ) )


# 2.4.7 The R package ‘unmarked’ (no code)
# ------------------------------------------------------------------------
