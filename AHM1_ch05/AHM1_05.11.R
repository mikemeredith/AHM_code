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
# ~~~~~ and this bit from section 5.6 ~~~~~~~~~~
# Generate factor and plot raw data in boxplot as function of factor A
facFor <- as.numeric(forest < -0.5)         # Factor level 1
facFor[forest < 0 & forest > -0.5] <- 2     # Factor level 2
facFor[forest < 0.5 & forest > 0] <- 3      # Factor level 3
facFor[forest > 0.5] <- 4                   # Factor level 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.11 Binomial generalised linear model (binomial GLM, logistic regression)
# ==========================================================================


# Quantize counts from first survey and describe
y1 <- as.numeric(C[,1] > 0)  # Gets 1 if first count greater than zero
table(y1)

mean(N > 0)          # True occupancy
mean(y1)             # Observed occupancy after first survey

# Bundle data
win.data <- list(y1 = y1, M = length(y1), elev = elev, facFor = facFor)

# Specify model in BUGS language
cat(file = "Bernoulli_GLM.txt","
model {

  # Priors
  for(k in 1:4){
    alpha[k] <- logit(mean.psi[k])     # intercepts
    mean.psi[k] ~ dunif(0,1)
    beta[k] ~ dnorm(0, 1.0E-06)        # slopes
  }

  # Likelihood
  for (i in 1:M){
    y1[i] ~ dbern(theta[i])
    logit(theta[i]) <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
  }
}
")

# Initial values
inits <- function() list(mean.psi = runif(4), beta = rnorm(4,,3))   # Priors 2

# Parameters monitored
params <- c("mean.psi", "alpha", "beta", "theta")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS or JAGS from R (ART <1 min)
out6 <- bugs(win.data, inits, params, "Bernoulli_GLM.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)  # ~~~~ for autotesting

out6J <- jags(win.data, inits, params, "Bernoulli_GLM.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
# par(mfrow = c(4,2))  # ~~~ replaced with 'layout' argument
traceplot(out6J, c("alpha[1:4]", "beta[1:4]"), layout=c(4,2))
print(out6, 2)

# Compare with MLEs
summary(glm(y1 ~ factor(facFor)*elev-1-elev, family = binomial))


# Plot of observed response vs. two covariates
op <- par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
F1 <- facFor == 1 ; F2 <- facFor == 2 ; F3 <- facFor == 3 ; F4 <- facFor == 4
plot(jitter(y1,,0.05) ~ facFor, xlab = "Forest factor",
    ylab = "Observed occupancy probability", frame.plot = FALSE, ylim = c(0, 1.15))
lines(1:4, out6$summary[1:4,1], lwd = 2)
segments(1:4, out6$summary[1:4,3], 1:4, out6$summary[1:4,7])
text(1.15, 1.1, "A", cex=1.6)

plot(elev[F1], jitter(y1,,0.1)[F1], xlab = "Elevation", ylab = "", col = "red", frame.plot = FALSE)
points(elev[F2], jitter(y1,,0.05)[F2], col = "blue")
points(elev[F3], jitter(y1,,0.05)[F3], col = "green")
points(elev[F4], jitter(y1,,0.05)[F4], col = "grey")
lines(sort(elev[F1]), out6$mean$theta[F1][order(elev[F1])], col="red", lwd=2)
lines(sort(elev[F2]), out6$mean$theta[F2][order(elev[F2])], col="blue", lwd=2)
lines(sort(elev[F3]), out6$mean$theta[F3][order(elev[F3])], col="green", lwd=2)
lines(sort(elev[F4]), out6$mean$theta[F4][order(elev[F4])], col="grey", lwd=2)
text(-0.9, 1.1, "B", cex=1.6)
par(op)
