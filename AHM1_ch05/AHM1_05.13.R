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


# 5.13 Random-effects Poisson GLM (Poisson GLMM)
# ==============================================


# Bundle data
win.data <- list(C = C, M = nrow(C), J = ncol(C), elev = elev, forest = forest,
    elev.forest = elev * forest, wind = wind)

# Specify model in BUGS language
cat(file = "RE.Poisson.txt","
model {

# Priors
mu.alpha ~ dnorm(0, 0.001)                # Mean hyperparam
tau.alpha <- pow(sd.alpha, -2)
sd.alpha ~ dunif(0, 10)                   # sd hyperparam
for(k in 1:4){
   alpha[k] ~ dunif(-10, 10)              # Regression params
}

# Likelihood
for (i in 1:M){
   alpha0[i] ~ dnorm(mu.alpha, tau.alpha) # Random effects and hyperparams
   re0[i] <- alpha0[i] - mu.alpha         # zero-centered random effects
   for(j in 1:J){
      C[i,j] ~ dpois(lambda[i,j])
      log(lambda[i,j]) <- alpha0[i] + alpha[1] * elev[i] + alpha[2] * forest[i] +
          alpha[3] * elev.forest[i] + alpha[4] * wind[i,j]
   }
}
}")

# Other model run preparations
inits <- function() list(alpha0 = rnorm(M), alpha = rnorm(4)) # Inits
params <- c("mu.alpha", "sd.alpha", "alpha0", "alpha", "re0") # Params
ni <- 30000 ; nt <- 25 ; nb <- 5000 ; nc <- 3                 # MCMC settings

# Call WinBUGS or JAGS from R (ART 6-7 min) and summarize posteriors
out8 <- bugs(win.data, inits, params, "RE.Poisson.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)   # ~~~~ for autotesting

out8 <- jags(win.data, inits, params, "RE.Poisson.txt",
  # n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
# par(mfrow = c(3,2))  # ~~~ replace with 'layout'  argument
traceplot(out8, c("mu.alpha", "sd.alpha", "alpha[1:3]"), layout=c(3,2))

print(out8, 3)

Cvec <- as.vector(C)            # Vector of M*J counts
elev.vec <- rep(elev, J)        # Vectorized elevation covariate
forest.vec <- rep(forest, J)    # Vectorized forest covariate
wind.vec <- as.vector(wind)     # Vectorized wind covariate
fac.site <- factor(rep(1:M, J)) # Site indicator (factor)
cbind(Cvec, fac.site, elev.vec, forest.vec, wind.vec) # Look at data

# Fit same model using maximum likelihood (NOTE: glmer uses ML instead of REML)
library(lme4)
summary(fm <- glmer(Cvec ~ elev.vec*forest.vec + wind.vec + (1| fac.site),
    family = poisson))              # Fit model
ranef(fm)                       # Print zero-centered random effects


# Compare fixed-effects estimates (in spite of the confusing naming in glmer output),
#   Bayesian post. means and sd left, frequentist MLEs and SEs right
print(cbind(out8$summary[c(1:2, 270:273), 1:2], rbind(summary(fm)$coef[1,1:2],
    c(sqrt(summary(fm)$varcor$fac.site), NA), summary(fm)$coef[c(2,3,5,4),1:2])), 3)

# Compare graphically non-Bayesian and Bayesian random effects estimates
Freq.re <- ranef(fm)$fac.site[,1]         # Non-Bayesian estimates (MLEs)
Bayes.re <- out8$summary[274:540,]        # Bayesian estimates

op <- par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(Freq.re, Bayes.re[,1], xlab = "Non-Bayesian (glmer)",
    ylab = "Bayesian (BUGS)", xlim = c(-0.4, 0.42), ylim = c(-2, 2),
    frame.plot = FALSE, type = "n")
segments(Freq.re, Bayes.re[,3], Freq.re, Bayes.re[,7], col = "grey", lwd = 0.5)
abline(0, 1, lwd = 2)
points(Freq.re, Bayes.re[,1])
text(-0.38, 2, "A", cex=1.6)

wind.pred <- seq(-1, 1, , 1000)     # Covariate values for prediction
pred <- array(NA, dim = c(1000, 267))
for(i in 1:267){
  pred[,i]<- exp(out8$mean$alpha0[i] + out8$mean$alpha[1] * 0 + out8$mean$alpha[2] * 0 +
      out8$mean$alpha[3] * 0 + out8$mean$alpha[4] * wind.pred)    # Predictions for each site
}

matplot(wind.pred, pred, type = "l", lty = 1, col = "grey",
    xlab = "Wind speed", ylab = "Expected count", frame.plot = FALSE, ylim = c(0, 4))
lines(wind.pred, exp(out8$mean$mu.alpha + out8$mean$alpha[4] * wind.pred), col = "black", lwd = 3)
text(-0.9, 4, "B", cex=1.6)
par(op)
