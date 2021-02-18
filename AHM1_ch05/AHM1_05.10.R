#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5. Fitting models using the Bayesian modeling software BUGS and JAGS
# =========================================================================

library(AHMbook)
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
# ~~~~~ and this from section 5.9 ~~~~~~~~~~
# Summarize data by taking max at each site
Cmax <- apply(C, 1, max)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 5.10 Goodness of fit assessment: posterior predictive checks and the parametric bootstrap
# =========================================================================================

# Bundle data
win.data <- list(Cmax = Cmax, M = length(Cmax), elev = elev, facFor = facFor, e = 0.0001)

# Specify model in BUGS language
cat(file = "Poisson_GLM.txt","
model {
  # Priors
  for(k in 1:4){
    alpha[k] ~ dnorm(0, 1.0E-06)
    beta[k] ~ dnorm(0, 1.0E-06)
  }

  # Likelihood and computations for posterior predictive check
  for (i in 1:M){
    Cmax[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]

    # Fit assessments: Chi-squared test statistic and posterior predictive check
    chi2[i] <- pow((Cmax[i]-lambda[i]),2) / (sqrt(lambda[i])+e)         # obs.
    Cmax.new[i] ~ dpois(lambda[i])      # Replicate (new) data set
    chi2.new[i] <- pow((Cmax.new[i]-lambda[i]),2) / (sqrt(lambda[i])+e) # exp.
  }
  # Add up discrepancy measures for entire data set
  fit <- sum(chi2[])                     # Omnibus test statistic actual data
  fit.new <- sum(chi2.new[])             # Omnibus test statistic replicate data

  # range of data as a second discrepancy measure
  obs.range <- max(Cmax[]) - min(Cmax[])
  exp.range <- max(Cmax.new[]) - min(Cmax.new[])
}
")

# Initial values
inits <- function() list(alpha = rnorm(4,,3), beta = rnorm(4,,3))

# Parameters monitored
params <- c("chi2", "fit", "fit.new", "obs.range", "exp.range")
params <- c("Cmax.new", "chi2.new", "chi2", "fit", "fit.new", "obs.range", "exp.range")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call JAGS from R and summarize posteriors
out5.1 <- jags(win.data, inits, params, "Poisson_GLM.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, modules = 'bugs')

print(out5.1, 2)

op <- par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(out5.1$sims.list$fit, out5.1$sims.list$fit.new, xlim = c(200, 1000),
    ylim = c(200, 1000), main = "", xlab = "Discrepancy observed data",
    ylab = "Discrepancy expected data", frame.plot = FALSE, cex = 1.5)
abline(0,1, lwd = 2)
text(240, 1000, "A", cex = 2)



(bpv <- mean(out5.1$sims.list$fit.new > out5.1$sims.list$fit))

plot(Cmax, out5.1$mean$chi2, xlab = "Observed data", ylab = "Chi2 contribution",
    frame.plot = FALSE, ylim = c(0, 70), cex = 1.5)
lines(0:30, sqrt(0:30), lwd = 2)
text(2, 70, "B", cex = 2)
par(op)

# Fit model to actual data and compute two fit statistics
fm <- glm(Cmax ~ as.factor(facFor)*elev, family = poisson) # Fit model
observed <- Cmax
expected <- predict(fm, type = "response")
plot(observed, expected)
abline(0,1)
chi2.obs <- sum((observed - expected)^2 / (expected + 0.0001)) # fit stat 1
range.obs <- diff(range(observed))                             # fit stat 2

# Generate reference distribution of fit stat under a fitting model
# simrep <- 100000      # Might want to try 1000 first, takes a looong time!
simrep <- 1000
chi2vec <- rangevec <- numeric(simrep)  # vectors for chi2 and maximum
for(i in 1:simrep){
  cat(paste(i, "\n"))
  Cmaxrep <- rpois(n = 267, lambda = expected)   # Generate replicate data set
  fmrep <- glm(Cmaxrep ~ as.factor(facFor)*elev, family = poisson) # Refit model
  expectednew <- predict(fmrep, type = "response")
  chi2vec[i] <- sum((Cmaxrep - expectednew)^2 / (expectednew + 0.0001))
  rangevec[i] <- diff(range(Cmaxrep))
}

# Summarize bootstrap results and compare with posterior predictive dist.
op <- par(mfrow = c(2, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(out5.1$sims.list$fit.new, col = "grey", main = "", breaks = 100,
    xlim = c(180, 900), freq = FALSE, ylim = c(0, 0.01))
abline(v = mean(out5.1$sims.list$fit), col = "red", lwd = 2)
text(200, 0.009, "A", cex = 2)
hist(out5.1$sims.list$exp.range, col = "grey", main = "", breaks = 50,
    xlim = c(10, 32), freq = FALSE, ylim = c(0, 0.40))
abline(v = mean(out5.1$sims.list$obs.range), col = "red", lwd = 2)
text(11, 0.38, "B", cex = 2)
hist(chi2vec, col = "grey", main = "", breaks = 100, xlim = c(180, 900),
    freq = FALSE, ylim = c(0, 0.02))
abline(v = chi2.obs, col = "red", lwd = 2)
text(200, 0.018, "C", cex = 2)
hist(rangevec, col = "grey", main = "", breaks = 50, xlim = c(10, 32), freq = FALSE)
abline(v = range.obs, col = "red", lwd = 2)
text(11, 0.38, "D", cex = 2)
par(op)

# Lack of fit ratio in PPD and parboot
mean(out5.1$sims.list$fit/out5.1$sims.list$fit.new) # ppc
mean(chi2.obs/chi2vec)                              # parboot

(pval1 <- 1-rank(c(chi2vec, chi2.obs))[simrep+1]/(simrep+1))
(pval2 <- 1-rank(c(rangevec, range.obs))[simrep+1]/(simrep+1))

