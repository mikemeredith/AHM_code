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

# Summarize data by taking mean at each site and plot
Cmean <- apply(C, 1, mean)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.6 Linear model with normal response (normal GLM): analysis of covariance (ANCOVA)
# -----------------------------------------------------------------------------------

# Generate factor and plot raw data in boxplot as function of factor A
facFor <- as.numeric(forest < -0.5)         # Factor level 1
facFor[forest < 0 & forest > -0.5] <- 2     # Factor level 2
facFor[forest < 0.5 & forest > 0] <- 3      # Factor level 3
facFor[forest > 0.5] <- 4                   # Factor level 4
table(facFor)                               # every site assigned a level OK

op <- par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(Cmean ~ factor(facFor), col = c("red", "blue", "green", "grey"),
    xlab = "Forest cover class", ylab = "Mean count of great tits",
    frame.plot = FALSE, ylim = c(0,20))
text(0.8, 20, "A", cex=1.6)


# Bundle data
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = elev, facFor = facFor)

# Specify model in BUGS language in effects parameterisation
cat(file = "ANCOVA1.txt","
model {

  # Priors
  alpha ~ dnorm(0, 1.0E-06)            # Prior for intercept = effect of level 1 of forest factor
  beta2 ~ dnorm(0, 1.0E-06)            # Prior for slope = effect of elevation for level 1 of forest factor
  beta1[1] <- 0                        # Set to zero effect of first level of facFor
  beta3[1] <- 0                        # Set to zero effect of first level of facFor of elevation
  for(k in 2:4){
    beta1[k] ~ dnorm(0, 1.0E-06)       # Prior for effects of factor facFor
    beta3[k] ~ dnorm(0, 1.0E-06)       # Prior for effects of factor facFor
  }
  tau <- pow(sd, -2)
  sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale

  # Likelihood
  for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)          # precision tau = 1 / variance
    mu[i] <- alpha + beta1[facFor[i]] + beta2 * elev[i] + beta3[facFor[i]] * elev[i]
  }
}
")

# Initial values
inits <- function() list(alpha = rnorm(1,,10), beta1 = c(NA, rnorm(3,,10)),
    beta2 = rnorm(1,,10), beta3 = c(NA, rnorm(3,,10)))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "sd")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS or JAGS from R (ART <1 min)
out3 <- bugs(win.data, inits, params, "ANCOVA1.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)  # ~~~~ for automated testing

out3J <- jags(win.data, inits, params, "ANCOVA1.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb)
traceplot(out3J)

# Fit model using least-squares (produces MLEs)
(fm <- summary(lm(Cmean ~ as.factor(facFor)*elev)))

# Summarize posteriors
print(out3, 3)


# Specify model in BUGS language
cat(file = "ANCOVA2.txt","
model {

  # Priors
  for(k in 1:4){
    alpha[k] ~ dnorm(0, 1.0E-06)       # Priors for intercepts
    beta[k] ~ dnorm(0, 1.0E-06)        # Priors for slopes
  }
  tau <- pow(sd, -2)
  sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale

  # Likelihood
  for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)          # precision tau = 1 / variance
    mu[i] <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
  }

  # Derived quantities: comparison of slopes (now you can forget the delta rule !)
  for(k in 1:4){
    diff.vs1[k] <- beta[k] - beta[1]    # Differences relative to beta[1]
    diff.vs2[k] <- beta[k] - beta[2]    # ... relative to beta[2]
    diff.vs3[k] <- beta[k] - beta[3]    # ... relative to beta[3]
    diff.vs4[k] <- beta[k] - beta[4]    # ... relative to beta[4]
  }
}
")

# Initial values
inits <- function() list(alpha = rnorm(4,,10), beta = rnorm(4,,10))

# Parameters monitored
params <- c("alpha", "beta", "sd", "diff.vs1", "diff.vs2", "diff.vs3", "diff.vs4")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS or JAGS from R (ART <1 min) and summarize posteriors
out4 <- bugs(win.data, inits, params, "ANCOVA2.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)  # ~~~~ for automated testing

system.time(out4J <- jags(win.data, inits, params, "ANCOVA2.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb))
traceplot(out4J)

print(out4, 2)

# Fit model using least-squares (produces MLEs)
(fm <- summary(lm(Cmean ~ as.factor(facFor)*elev-1-elev)))


plot(elev[facFor==1], Cmean[facFor==1], col = "red", ylim = c(0, 20),
    xlab = "Elevation", ylab = "", frame.plot = FALSE)
points(elev[facFor==2], Cmean[facFor==2], col = "blue")
points(elev[facFor==3], Cmean[facFor==3], col = "green")
points(elev[facFor==4], Cmean[facFor==4], col = "black")
abline(fm$coef[1,1], fm$coef[5,1], col = "red")
abline(fm$coef[2,1], fm$coef[6,1], col = "blue")
abline(fm$coef[3,1], fm$coef[7,1], col = "green")
abline(fm$coef[4,1], fm$coef[8,1], col = "black")
text(-0.8, 20, "B", cex=1.6)
par(op)

attach.bugs(out4)     # Allows to directly address the sims.list
str(diff.vs3)
op <- par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(diff.vs3[,1], col = "grey", breaks = 100, main = "", freq=FALSE, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-1.2, 0.8, "A", cex = 2)
hist(diff.vs3[,2], col = "grey", breaks = 100, main = "", freq=FALSE, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-1.4, 0.8, "B", cex = 2)
hist(diff.vs3[,4], col = "grey", breaks = 100, main = "", freq=FALSE, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-2.2, 0.8, "C", cex = 2)
par(op)

# Prob. difference greater than 1
mean(diff.vs3[,1] > 1)
mean(diff.vs3[,2] > 1)
mean(diff.vs3[,4] > 1)

# 5.7 Proportion of variance explained (R2)
# =========================================

cat(file = "Model0.txt","
model {
  # Priors
  mu ~ dnorm(0, 1.0E-06)
  tau <- pow(sd, -2)
  sd ~ dunif(0, 1000)
  # Likelihood
  for (i in 1:M){
    Cmean[i] ~ dnorm(mu, tau)
  }
}
")
inits <- function() list(mu = rnorm(1))
params <- c("mu", "sd")
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3
out0 <- jags(win.data, inits, params, "Model0.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb)

print(out0)

# Compute R2 from BUGS analysis
(total.var <- mean(out0$sims.list$sd^2))       # Total variance around the mean
(unexplained.var <- mean(out3$sims.list$sd^2)) # Not explained by the ANCOVA
(prop.explained <- (total.var - unexplained.var)/total.var)

