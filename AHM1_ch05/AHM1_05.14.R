#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5. Fitting models using the Bayesian modeling software BUGS and JAGS
# =========================================================================

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14/"          # Place where your WinBUGS installed
library(jagsUI)

# ~~~~~ this section requires the following code from section 5.3 ~~~~~~~~~~
# Generate data with data.fn from chapter 4
set.seed(24)
data <- data.fn(show.plot=FALSE)
attach(data)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.14 Random-effects binomial GLM (binomial GLMM)
# ================================================


# Get detection/nondetection response
y <- C
y[y > 0] <- 1

# Bundle data
win.data <- list(y = y, M = nrow(y), J = ncol(y), elev = elev, forest = forest,
    elev.forest = elev * forest, wind = wind)
str(win.data)

# Specify model in BUGS language
cat(file = "RE.Bernoulli.txt","
model {

  # Priors
  mu.alpha0 <- logit(mean.theta)              # Random intercepts
  mean.theta ~ dunif(0,1)
  tau.alpha0 <- pow(sd.alpha0, -2)
  sd.alpha0 ~ dunif(0, 10)
  mu.alpha4 ~ dnorm(0, 0.001)                 # Random slope on wind
  tau.alpha4 <- pow(sd.alpha4, -2)
  sd.alpha4 ~ dunif(0, 10)
  for(k in 1:3){
    alpha[k] ~ dnorm(0, 0.001)               # Slopes
  }

  # Likelihood
  for (i in 1:M){
    alpha0[i] ~ dnorm(mu.alpha0, tau.alpha0) # Intercept random effects
    re00[i] <- alpha0[i] - mu.alpha0         # same zero-centered
    alpha4[i] ~ dnorm(mu.alpha4, tau.alpha4) # Slope random effects
    re04[i] <- alpha4[i] - mu.alpha4         # same zero-centered
    for(j in 1:J){
      y[i,j] ~ dbern(theta[i,j])
      logit(theta[i,j]) <- alpha0[i] + alpha[1] * elev[i] + alpha[2] * forest[i] +
          alpha[3] * elev.forest[i] + alpha4[i] * wind[i,j]
    }
  }
}")

# Other model run preparations
inits <- function() list(alpha0 = rnorm(M), alpha4 = rnorm(M))# Inits
params <- c("mu.alpha0", "sd.alpha0", "alpha0", "alpha", "mu.alpha4", "sd.alpha4",
    "alpha4", "re00", "re04")                        # Params
ni <- 30000 ; nt <- 25 ; nb <- 5000 ; nc <- 3                 # MCMC settings

# Call WinBUGS from R .... and crash !
# ~~~~~~~ ...and hence commented out! ~~~~~~~~~~~~~~~~~~~~~~~~
# out9 <- bugs(win.data, inits, params, "RE.Bernoulli.txt", n.chains = nc,
#   n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir,
#   working.directory = getwd())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Call JAGS from R (ART 2.5 min)
out9 <- jags(win.data, inits, params, "RE.Bernoulli.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out9, c("mu.alpha0", "sd.alpha0", "alpha[1:3]", "mu.alpha4", "sd.alpha4"),
    layout=c(2,2))
print(out9, 3)

yvec <- as.vector(y)            # Vector of M*J counts
elev.vec <- rep(elev, J)        # Vectorized elevation covariate
forest.vec <- rep(forest, J)    # Vectorized forest covariate
wind.vec <- as.vector(wind)     # Vectorized wind covariate
fac.site <- factor(rep(1:M, J)) # Site indicator (factor)
cbind(yvec, fac.site, elev.vec, forest.vec, wind.vec) # Look at data

# Fit same model using maximum likelihood
library(lme4)                   # Load package
summary(frem <- glmer(yvec ~ elev.vec*forest.vec + wind.vec + (wind.vec || fac.site),
    family = binomial))              # Fit model


# Compare Bayesian and non-Bayesian estimates
print(out9$summary[c(1:2, 270:274),c(1:3,7:9)], 4)

(re <- ranef(frem))                 # Print zero-centered random effects

op <- par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
pop.mean.int <- summary(frem)$coef[1,1]
pop.mean.slope <- summary(frem)$coef[4,1]
plot(sort(wind.vec), plogis(pop.mean.int + sort(wind.vec) * pop.mean.slope),
    type = "l", xlab = "Wind speed", ylab = "Prob. to count >0 great tits",
    lwd = 3, frame.plot = FALSE)
for(i in 1:267){
  lines(sort(wind.vec), plogis(pop.mean.int + re$fac.site[i,1] + sort(wind.vec) *
      (pop.mean.slope + re$fac.site[i,2])), lwd = 1, col = i)
}
title(main = "A", cex.main = 2)

# Compute expected detection/nondetection probability for a grid of elevation and
#   forest cover, at wind-speed = 0 (covariate average) and for hypermean of
#   intercepts alpha0
n.sims <- length(out9$sims.list$mu.alpha0)
elev.pred <- seq(-1, 1,,100)                       # Values of elevation
forest.pred <- seq(-1,1,,100)                      # Values of forest cover
pred.array <- array(NA, dim = c(100, 100, n.sims)) # Prediction array
for(i in 1:100){
  for(j in 1:100){
    pred.array[i,j,] <- plogis(out9$sims.list$mu.alpha0 +
        out9$sims.list$alpha[,1] * elev.pred[i] + out9$sims.list$alpha[,2] *
        forest.pred[j] + out9$sims.list$alpha[,3] * elev.pred[i] * forest.pred[j])
  }
}
pm.pred.array <- apply(pred.array, c(1,2), mean)   # Get posterior mean
psd.pred.array <- apply(pred.array, c(1,2), sd)    # Get posterior sd

mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=elev.pred, y= forest.pred, z=pm.pred.array, col = mapPalette(100),
    xlab = "Elevation", ylab = "Forest cover")
contour(x=elev.pred, y=forest.pred, z= pm.pred.array, add = TRUE, lwd = 1)
title(main = "B", cex.main = 2)

image(x=elev.pred, y= forest.pred, z=psd.pred.array, col = mapPalette(100),
    xlab = "Elevation", ylab = "Forest cover")
contour(x=elev.pred, y=forest.pred, z= psd.pred.array, add = TRUE, lwd = 1)
title(main = "C", cex.main = 2)
par(op)
