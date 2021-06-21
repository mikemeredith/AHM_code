#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

library(AHMbook)
library(jagsUI)

# 4.6 Trend estimation with occupancy data
# ========================================

# Pick colonization, extinction; compute equilibrium occupancy
gam <- 0.1
eps <- 0.2
( psi.eq <- gam / (gam + eps) )
# [1] 0.3333333

# ~~~~ code for figure 4.7 is at the end of the script ~~~~

# Simulate a data set
set.seed(24)
data <- simDynocc(mean.psi1 = 0.8, range.phi = c(1-eps, 1-eps),
    range.gamma = c(gam, gam), range.p = c(0.4, 0.4))    # library(AHMbook)

# Stack detection histories
ystack <- array(NA, dim = c(data$nsites * data$nyears, data$nsurveys))
for(t in 1:data$nyears){
  ystack[((t-1)*data$nsites+1):(data$nsites*t),] <- data$y[,,t]
}

# Create year covariate and factor (both site covariates in unmarked)
year <- 1:data$nyears
yr <- rep(1:data$nyears, each = data$nsites)             # year as cont. cov.
yrfac <- as.factor(yr)                                   # year as a factor

# Format and summarize data
library(unmarked)
summary(umf <- unmarkedFrameOccu(y = ystack, siteCovs = data.frame(yr = yr,
    yrfac = yrfac)) )                                    # require(unmarked)

# (1) Analysis of stacked data: treating year as a factor
summary(fm1 <- occu(~1 ~yrfac-1, data = umf))
nd <- data.frame(yrfac = as.factor(1:10))
pred.yrfac <- predict(fm1, type = "state", newdata = nd)

# (2) Analysis of stacked data: fitting a trend of year
summary(fm2 <- occu(~1 ~yr, data = umf))
nd <- data.frame(yr = 1:10)
pred.yr <- predict(fm2, type = "state", newdata = nd)

# Create array to hold predictions
simrep <- 1000
bs.pred <- array(NA, dim = c(nrow(pred.yr), simrep))

# Nonparametric bootstrap for prediction of trend line (ART 2 min)
for(b in 1:simrep){
  cat(paste("\n** Bootstrap rep", b, "**"))
  # Sample sites with replacement (get their index)
  bs.samp <- sample(1:250, 250, replace = TRUE)
  # Repeat stacking with bootstrap sample
  ystack <- array(NA, dim = c(250 * 10, 3))
  for(t in 1:data$nyears){
    ystack[((t-1)*data$nsites+1):(data$nsites*t),] <- data$y[bs.samp,,t]
  }
  # Create unmarked data frame, fit model and predict trend line
  umf <- unmarkedFrameOccu(y = ystack, siteCovs = data.frame(yr = yr))
  fm <- occu(~1 ~ yr, data = umf, se = FALSE)
  tmp <- predict(fm, type = "state", newdata = data.frame(yr = 1:10))
  # Save param estimates
  bs.pred[,b] <- tmp[,1]
}

# Get bootstrap SE and CI for all annual predictions
se.bs <- apply(bs.pred, 1, sd)
ci.bs <- t(apply(bs.pred, 1, function(x)quantile(x, c(0.025, 0.975))))

# Compare asymptotic and bootstrapped SEs and CIs
round(cbind('ASE' = pred.yr[,2], 'Asymp_LCL' = pred.yr[,3],
    'Asymp_UCL' = pred.yr[,4], 'Bootstrapped SE' = se.bs,
    'Bootstrapped CI' = ci.bs), 3)
#         ASE Asymp_LCL Asymp_UCL Bootstrapped SE  2.5% 97.5%
# [1,]  0.023     0.632     0.722           0.030 0.620 0.735
# [2,]  0.021     0.591     0.673           0.028 0.579 0.686
# [3,]  0.019     0.548     0.622           0.026 0.535 0.635
# [4,]  0.017     0.503     0.568           0.024 0.490 0.583
# [5,]  0.015     0.456     0.515           0.023 0.443 0.531
# [6,]  0.014     0.408     0.464           0.022 0.393 0.477
# [7,]  0.015     0.358     0.417           0.023 0.342 0.432
# [8,]  0.016     0.310     0.373           0.024 0.292 0.386
# [9,]  0.017     0.264     0.332           0.025 0.247 0.347
# [10,] 0.018     0.222     0.295           0.026 0.206 0.308

# Bundle and summarize data
str(bdata <- list(y = data$y, nsites = dim(data$y)[1],
    nsurveys = dim(data$y)[2], nyears = dim(data$y)[3]))

# Specify model in BUGS language
cat(file = "occ.txt","
model {

  # Specify priors and specify model for occupancy
  for (t in 1:nyears){
    logit(psi[t]) <- alpha + beta.trend * (t-5.5) + eps.year[t]
    p[t] ~ dunif(0, 1)
    eps.year[t] ~ dnorm(0, tau.lpsi)
  }
  alpha <- logit(mean.psi)
  mean.psi ~ dunif(0,1)
  beta.trend ~ dnorm(0, 0.01)
  tau.lpsi <- pow(sd.lpsi, -2)
  sd.lpsi ~ dunif(0, 10)

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){
    for (t in 1:nyears){
      z[i,t] ~ dbern(psi[t])
      # Observation model
      for (j in 1:nsurveys){
        y[i,j,t] ~ dbern(z[i,t]*p[t])
      }
    }
  }

  # Derived parameters
  for (t in 1:nyears){
    n.occ[t] <- sum(z[1:nsites,t]) # Finite sample occupancy
    logit(psi.trend[t]) <- alpha + beta.trend * (t-5.5)
  }
}
")

# Initial values
inits <- function(){ list(z = apply(data$y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "psi.trend", "mean.psi", "alpha", "beta.trend",
    "sd.lpsi", "p", "n.occ")

# MCMC settings
na <- 1000 ; ni <- 6000 ; nt <- 1 ; nb <- 2000 ; nc <- 3

# Call JAGS (ART 2 min), check convergence and summarize posteriors
out <- jags(bdata, inits, params, "occ.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out)
print(out, dig = 2) # not shown

# ~~~~~~~~~~ code for figure 4.7 ~~~~~~~~~
op <- par(mfrow = c(1, 3), mar = c(5,5,5,2), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
# Plot A: True and observed proportion of occupied sites
plot(year, data$psi.fs, xlab = "Year", ylab = "Occupancy", col = "red", xlim = c(0,data$nyear+1),
    ylim = c(0,1), lwd = 2, frame.plot = F, type = 'b', main = 'True and observed occupancy',
    cex = 2, pch = 16)
points(year, data$psi.app, type = 'b', col = "black", pch = 16, cex = 2, lwd = 3)
legend('bottomleft', c('True', 'Observed'), col = c('red', 'black'), bty = 'n', cex = 1.2,
    lty = c(1,1,2), pch = c(16, 16, NA), lwd = 3)
text(0.3, 0.99, 'A', cex = 2)

# Plot B: True occupancy and estimates from stacked data-analyses in unmarked
plot(year, data$psi.fs, xlab = "Year", ylab = "Occupancy", col = "red", xlim = c(0,data$nyear+1),
    ylim = c(0,1), lwd = 2, frame.plot = F, type = 'b', main = 'unmarked', cex = 2, pch = 16)

off <- 0.3
points(year+off, pred.yrfac[,1], col = "blue", pch = 16, cex = 2)
segments(year+off, pred.yrfac[,3], year+off, pred.yrfac[,4], col = 'blue', lwd =2)
lines(year, pred.yr[,1], col = 'blue', lwd = 2)
matlines(year, pred.yr[,3:4], col = 'blue', lty = 1, lwd = 1)
legend('bottomleft', c('True', 'Estimated trend', 'Estimated year means'),
    col = c('red', 'blue', 'blue'), bty = 'n', cex = 1.2, lty = c(1,1,NA),
    pch = c(16, NA, 16), lwd = 3)
text(0.3, 0.99, 'B', cex = 2)

# Plot C: True occupancy and estimates from JAGS
off <- 0.3
plot(year, data$psi.fs, xlab = "Year", ylab = "Occupancy", col = "red", xlim = c(0,data$nyear+1),
    ylim = c(0,1), lwd = 2, frame.plot = F, type = 'b', main = 'JAGS', cex = 2, pch = 16)
points(year+off, out$mean$psi, col = "blue", pch = 16, cex = 2)
segments(year+off, out$q2.5$psi, year+off, out$q97.5$psi, col = 'blue', lwd = 2)
lines(year, out$mean$psi.trend, col = "blue", lwd = 3)
matlines(year, cbind(out$q2.5$psi.trend, out$q97.5$psi.trend), col = 'blue', lty = 1, lwd = 1)
legend('bottomleft', c('Truth', 'Estimated trend', 'Estimated year means'),
    col = c('red', 'blue', 'blue'), bty = 'n', cex = 1.2, lty = c(1, 1, NA),
    pch = c(16, NA, 16), lwd = 3)
text(0.3, 0.99, 'C', cex = 2)
par(op)
