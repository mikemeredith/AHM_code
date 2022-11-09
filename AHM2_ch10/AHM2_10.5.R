#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10 : INTEGRATED MODELS FOR MULTIPLE TYPES OF DATA
# =========================================================
# Code from proofs dated 2020-08-19


# Approximate execution time for this code: 4 mins
# Run time with full number of iterations: 25 mins

library(AHMbook)
library(jagsUI)
library(berryFunctions)

# 10.5 Example 3: Combination of counts plus abundance-class data
# ======================================================================

library(AHMbook)

# Simulate two data sets from same process
# Constants for simulation function
M1 <- 500           # Number of sites with abundance-class surveys
J1 <- 5             # Number of surveys in abundance-class survey
M2 <- 100           # Number of sites with full-count surveys
J2 <- 2             # Number of surveys in full-count surveys
mean.lam <- 50      # Abundance intercept
beta.lam <- 0.5     # Coefficients of a site covariate in abundance
mean.p1 <- 0.4      # Detection intercept in abundance-class survey
mean.p2 <- 0.5      # Detection intercept in full-count survey
beta.p <- -1        # Coefficients of a site covariate in detection

# Simulate data set 1 for abundance-class surveys
set.seed(1)
str(data1 <- simNmix(nsite = M1, nvisit = J1, mean.lam = mean.lam,
    mean.p = mean.p1, beta2.lam = beta.lam, beta3.p = beta.p,
    show.plot = FALSE))

# Simulate data set 2 for full-count surveys
str(data2 <- simNmix(nsite = M2, nvisit = J2, mean.lam = mean.lam,
    mean.p = mean.p2, beta2.lam = beta.lam, beta3.p = beta.p,
    show.plot = FALSE))

library(berryFunctions) # Load package berryFunctions
breaks <- c(0, 10, 25, 50, 100, 200) # Up to and including
Cclass <- classify(c(data1$C), method = 'custom', breaks = breaks)$index
Aclass <- matrix(Cclass, byrow = FALSE, ncol = data1$nvisit)
head(data1$C)           # Look at the relationship between the raw data..
head(Aclass)            # ... and the abundance-class data

table(Aclass)           # Summary of the class data
# Aclass
#   1   2   3   4  5
# 687 897 543 348 25

breaks                  # Remember the class boundaries
# [1] 0 10 25 50 100 200

# ~~~ extra code for figure 10.3 ~~~~~~~~~~
# Visualize full and binned counts in first data set (Fig. 21-3)
op <- par(mfrow = c(1, 2), mar = c(6, 6, 5, 5), cex.lab = 1.5, cex.axis = 1.5,
    cex.main = 1.5)
plot(table(data1$C), type = 'h', lend = 'butt', lwd = 10, col = 'gray50',
    xlab = 'Counts', ylab = 'Number', frame = FALSE)
abline(v=breaks)
plot(data1$C, Aclass, xlab = 'Counts', ylab = 'Abundance class', frame = FALSE,
    pch = 16, col = 'gray50', cex = 2)
abline(v = breaks, lwd = 3)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get lower and upper class boundary for every binned count
n <- length(Cclass)
limits <- array(NA, c(n, 2), list(1:n, c('Lower', 'Upper')))
for(i in 1:n){
  limits[i, 1:2] <- c(breaks[Cclass[i]], breaks[Cclass[i]+1])
}
head(cbind('True count' = c(data1$C), 'Abundance class' = c(Aclass), limits))
#   True count Abundance class Lower Upper
# 1         15               2    10    25
# 2         21               2    10    25
# 3         44               3    25    50
# 4          8               1     0    10
# 5         41               3    25    50
# 6         98               4    50   100

# Vectorize the environmental covariate and get site index
X3vec <- rep(data1$site.cov[,3], 5)
sitevec <- rep(1:500, 5) # 500 sites with 5 reps each

# 10.5.1 Fitting a binomial N-mixture model with interval-censoring
#        to abundance-class data
# ------------------------------------------------------------------

# Response is simply a vector of ones !
y <- rep(1, 2500)

# Bundle data
str(bdata <- list(y = y, M = nrow(Aclass), X2 = data1$site.cov[,2],
    X3vec = X3vec, sitevec = sitevec, n = length(Cclass), limits = limits))
# List of 7
# $ y      : num [1:2500] 1 1 1 1 1 1 1 1 1 1 ...
# $ M      : int 500
# $ X2     : num [1:500] 0.217 0.753 0.632 0.653 -0.111 ...
# $ X3vec  : num [1:2500] 0.123 0.739 -0.467 1.82 -1.527 ...
# $ sitevec: int [1:2500] 1 2 3 4 5 6 7 8 9 10 ...
# $ n      : int 2500
# $ limits : num [1:2500, 1:2] 10 10 25 0 25 50 10 25 10 0 ...
# ..- attr(*, "dimnames") = List of 2
# .. ..$ : chr [1:2500] "1" "2" "3" "4" ...
# .. ..$ : chr [1:2] "Lower" "Upper"

# Specify model in BUGS language
cat(file = "model1.txt", "
model {
  # Priors
  alpha.p <- logit(mean.p)
  mean.p ~ dunif(0, 1)            # Detection intercept
  beta.p ~ dnorm(0, 0.1)          # Detection slope on X3
  alpha.lam ~ dnorm(0, 0.01)      # Abundance intercept
  mean.lam <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 0.1)        # Abundance slope on X2

  # Likelihood
  # Model for latent abundance
  for (i in 1:M){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha.lam + beta.lam * X2[i]
  }

  # Observation model for observed counts and for detection
  for (i in 1:n){
    y[i] ~ dinterval(C[i], limits[i,]) # specify interval censoring !
    C[i] ~ dbin(p[i], N[sitevec[i]])   # Count becomes estimated quantity
    logit(p[i]) <- alpha.p + beta.p * X3vec[i]
  }
  # Derived quantities
  Ntotal <- sum(N[])
}
")

# Initial values
# Need to give inits for both N and for C !
tmp <- matrix(limits[,2], ncol = 5)
Nst <- round(apply(tmp, 1, mean))+20
Cst <- limits[,1]+1
inits <- function() list(N = Nst, C = Cst, mean.p = 0.5,
    alpha.lam = log(mean(Nst)))

# Parameters monitored
params <- c("alpha.lam", "beta.lam", "mean.lam", "alpha.p", "beta.p",
    "mean.p", "Ntotal", "N", "C")

# MCMC settings
# na <- 1000 ; nc <- 3 ; ni <- 20000 ; nb <- 10000 ; nt <- 10
na <- 1000 ; nc <- 3 ; ni <- 2000 ; nb <- 1000 ; nt <- 1  # ~~~ for testing, 1 min

# Call JAGS (ART 13 min), gauge convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "model1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out1)
print(out1$summary[1:10,], 3) # not shown

# ~~~~ extra code for figure 10.4 ~~~~~~~~~~~~~~~~
# Graphs of estimated vs true N and same for C
op <- par(mfrow = c(1,2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5,
    cex.main = 1.5)
plot(data1$N, out1$mean$N, xlab = 'True N', ylab = 'Estimated N',
    frame = FALSE, pch = 16, cex = 1)
abline(0, 1, col = 'red', lwd = 3)
plot(data1$C, out1$mean$C, xlab = 'True C', ylab = 'Estimated C',
    frame = FALSE, pch = 16, cex = 1)
abline(0, 1, col = 'red', lwd = 3)
abline(v = breaks, lwd = 3)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Truth and estimates for some key estimands
# ~~~ extra code for the table ~~~~~~~~~~~~~~~~~~~~~~~
truth <- c('mean.p' = mean.p1, 'alpha1' = data1$beta3.p,
    'mean.lambda' = data1$mean.lam, 'beta1' = data1$beta2.lam,
    'Ntotal' = data1$Ntotal)
pm <- c(out1$mean$mean.p, out1$mean$beta.p, out1$mean$mean.lam,
    out1$mean$beta.lam, out1$mean$Ntotal)
LCI <- c(out1$q2.5$mean.p, out1$q2.5$beta.p, out1$q2.5$mean.lam,
    out1$q2.5$beta.lam, out1$q2.5$Ntotal)
UCI <- c(out1$q97.5$mean.p, out1$q97.5$beta.p, out1$q97.5$mean.lam,
    out1$q97.5$beta.lam, out1$q97.5$Ntotal)
cbind(truth, pm, LCI, UCI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               truth           pm            LCI           UCI
# mean.p          0.4     0.4377538     0.4031205     0.4749968
# alpha1         -1.0    -1.0433506    -1.1099861    -0.9804712
# mean.lambda    50.0    46.9308778    43.6895419    50.2425595
# beta1           0.5     0.5054756     0.4887845     0.5230504
# Ntotal      29955.0 28006.0973333 26153.8750000 29943.0750000


# 10.5.2 Fitting the integrated model
# -----------------------------------

# Bundle data
str(bdata <- list(y = y, M1 = nrow(Aclass), X21 = data1$site.cov[,2],
    X3vec1 = X3vec, sitevec = sitevec, n = length(Cclass), limits = limits,
    C2 = data2$C, M2 = nrow(data2$C), J2 = ncol(data2$C),
    X22 = data2$site.cov[,2], X32 = data2$site.cov[,3]))

# Specify model in BUGS language
cat(file = "model2.txt", "
model {
  # Priors
  # For data set 1 (binned counts)
  alpha.p1 <- logit(mean.p1)
  mean.p1 ~ dunif(0, 1)           # Detection intercept
  # For data set 2 (full counts)
  alpha.p2 <- logit(mean.p2)
  mean.p2 ~ dunif(0, 1)           # Detection intercept
  # Priors for the parameters shared in both data sets
  alpha.lam ~ dnorm(0, 0.01)      # Abundance intercept
  mean.lam <- exp(alpha.lam)
  beta.lam ~ dnorm(0, 0.1)        # Abundance slope on X2
  beta.p ~ dnorm(0, 0.1)          # Detection slope on X3

  # Likelihood
  # Model for latent abundance in both data sets
  # Note same parameter names means parameters are identical
  for (i in 1:M1){                # Data set 1
    N1[i] ~ dpois(lambda1[i])
    log(lambda1[i]) <- alpha.lam + beta.lam * X21[i]
  }
  for (i in 1:M2){                # Data set 2
    N2[i] ~ dpois(lambda2[i])
    log(lambda2[i]) <- alpha.lam + beta.lam * X22[i]
  }

  # Observation model for observed counts and for detection
  # Observation model for data set 1 (binned counts)
  for (i in 1:n){
    y[i] ~ dinterval(C1[i], limits[i,])       # interval censoring
    C1[i] ~ dbin(p1[i], N1[sitevec[i]])
    logit(p1[i]) <- alpha.p1 + beta.p * X3vec1[i]
  }
  # Observation model for data set 2 (full counts)
  for (i in 1:M2){ # Data set 2
    logit(p2[i]) <- alpha.p2 + beta.p * X32[i]
    for(j in 1:J2){
      C2[i,j] ~ dbin(p2[i], N2[i])
    }
  }

  # Derived quantities
  Ntotal1 <- sum(N1[])            # Total abundance in data set 1
  Ntotal2 <- sum(N2[])            # Total abundance in data set 2
  GTotalN <- Ntotal1 + Ntotal2    # Grand total at all 600 sites
}
")

# Initial values
tmp <- matrix(limits[,2], ncol = 5)
Nst1 <- round(apply(tmp, 1, mean))+20
Cst <- limits[,1]+1
Nst2 <- apply(bdata$C2, 1, max)+10
inits <- function() list(N1 = Nst1, C1 = Cst, N2 = Nst2)

# Parameters monitored
params <- c("mean.lam", "beta.lam", "mean.p1", "mean.p2", "beta.p",
    "Ntotal1", "Ntotal2", "GTotalN", "N1", "C1", "N2")

# MCMC settings
# na <- 1000 ; nc <- 3 ; ni <- 50000 ; nb <- 10000 ; nt <- 40
na <- 1000 ; nc <- 3 ; ni <- 5000 ; nb <- 1000 ; nt <- 4  # ~~~~ for testing, 2 mins

# Call JAGS (ART 22 min), gauge convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "model2.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out2)
jags.View(out2) # not shown

# Comparison of truth and estimates from simple model and IM
# ~~~~ extra code for the table ~~~~~~~~~~~~~~~~~~~~
truth <- c('mean.lam' = data1$mean.lam, 'beta.lam' = data1$beta2.lam,
    'mean.p1' = mean.p1, 'mean.p2' = mean.p2, 'beta.p' = data1$beta3.p,
    'Ntotal1' = data1$Ntotal, 'Ntotal2' = data2$Ntotal,
    'GNtotal' = data1$Ntotal + data2$Ntotal)
esti2 <- round(out2$summary[1:8, c(1,3,7)], 2)
esti1 <- array(NA, dim = c(8, 3))
esti1[1,] <- round(out1$summary[3, c(1,3,7)], 2)
esti1[2,] <- round(out1$summary[2, c(1,3,7)], 2)
esti1[3,] <- round(out1$summary[6, c(1,3,7)], 2)
esti1[5,] <- round(out1$summary[5, c(1,3,7)], 2)
esti1[6,] <- round(out1$summary[7, c(1,3,7)], 2)
colnames(esti1) <- c('mean(1)', '2.5% (1)', ' 97.5% (1)')
colnames(esti2) <- c('mean(2)', '2.5% (2)', ' 97.5% (2)')
cbind(truth, esti1, esti2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            truth     mean     2.5%    97.5%     mean     2.5%    97.5%
# mean.lam    50.0    46.93    43.69    50.24    48.74    46.06    51.99
# beta.lam     0.5     0.51     0.49     0.52     0.50     0.48     0.51
# mean.p1      0.4     0.44     0.40     0.47     0.42     0.39     0.45
# mean.p2      0.5       NA       NA       NA     0.52     0.48     0.56
# beta.p      -1.0    -1.04    -1.11    -0.98    -1.02    -1.08    -0.97
# Ntotal1  29955.0 28006.10 26153.88 29943.08 28921.67 27345.82 30844.10
# Ntotal2   5549.0       NA       NA       NA  5345.71  5056.93  5688.00
# GNtotal  35504.0       NA       NA       NA 34267.38 32432.00 36485.03
