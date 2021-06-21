#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

library(unmarked)

# 4.3 Simulation and analysis of the simplest dynamic occupancy model
# ===================================================================

# Choose sample sizes and prepare arrays for z and y
nsites <- 100            # Number of sites
nyears <- 12             # Number of years
nsurveys <- 2            # Number of presence/absence measurements
z <- array(NA, dim = c(nsites, nyears))           # latent presence/absence
y <- array(NA, dim = c(nsites, nsurveys, nyears)) # observed data

# Set parameter values as per above
psi1 <- 0.7                         # Prob. of initial occupancy or presence
phi <- 0.9                          # Persistence probability
gamma <- 0.05                       # Colonization probability
p <- 0.25                           # Probability of detection
(psi.eq <- gamma / (gamma+(1-phi))) # Equilibrium occupancy

# Generate initial presence/absence (i.e., the truth in year 1)
set.seed(1)               # So we all get same data set
z[,1] <- rbinom(n = nsites, size = 1, prob = psi1)
sum(z[,1]) / nsites       # True occupancy proportion in year 1
# [1] 0.68

# Generate presence/absence (i.e., the truth) in subsequent years
for(t in 2:nyears){
  exp.z <- z[,t-1] * phi + (1 - z[,t-1]) * gamma
  z[,t] <- rbinom(n = nsites, size = 1, prob = exp.z)
}
apply(z, 2, sum) / nsites # True occupancy proportions
# [1] 0.68 0.68 0.63 0.58 0.49 0.46 0.47 0.42 0.38 0.37 0.31 0.34

# Detection/nondetection data (i.e. presence/absence measurements)
for(t in 1:nyears){             # Loop over years 1 to 12
  for(j in 1:nsurveys){         # Loop over repeat visits 1 and 2
    y[,j,t] <- rbinom(n = nsites, size = 1, prob = z[,t] * p)
  }
}
y ; str(y)                      # Look at the data thus far simulated

# Generate missing values: create simple version of unbalanced data set
prob.missing <- 0.2           # Constant NA probability
y2 <- y                       # Duplicate balanced data set
for(i in 1:nsites){           # Loop over every datum in 3D array
  for(j in 1:nsurveys){
    for(t in 1:nyears){
      turnNA <- rbinom(1, 1, prob.missing)
      y2[i,j,t] <- ifelse(turnNA == 1, NA, y2[i,j,t])
    }
  }
}
y2 ; str(y2)                  # Look at the data now

table(nvisits <- apply(y2, c(1,3), function(x) sum(!is.na(x))))

# Compute true expected and realized occupancy (psi and psi.fs)
psi <- numeric(nyears) ; psi[1] <- psi1
for(t in 2:nyears){            # Compute true values of psi
  psi[t] <- psi[t-1] * phi + (1 - psi[t-1]) * gamma
}
psi.fs <- colSums(z) / 100     # True realized occupancy

# Compute observed occupancy proportion
zobs <- apply(y2, c(1,3), function(x) max(x, na.rm = TRUE))
zobs[zobs == '-Inf'] <- NA     # 13 site/years without visits
psi.obs <- apply(zobs, 2, sum, na.rm = TRUE) / apply(zobs, 2, function(x)
    sum(!is.na(x)))

# ~~~~~ code for figure 4.2 ~~~~~~~~~~~~
# Plot trajectories of psi, psi.fs, psi.eq and psi.obs
plot(1:nyears, psi, type = 'l', col = "red", lwd = 3, ylim = c(0,1), xlab = "Year",
    ylab = "Occupancy of Jura Asp Vipers", frame = F)
points(1:nyears, psi.fs, type = 'b', col = "red", pch = 16, cex = 2)
abline(h = gamma / (gamma+(1-phi)), lty = 2)
lines(1:nyears, psi.obs, col = "blue", lwd = 3)
legend('topright', c('True psi', 'True finite-sample psi', 'Equilibrium psi', 'Observed psi'),
    col = c('red', 'red', 'black', 'blue'), pch = c(NA, 16, NA, NA), lty = c(1,1,2,1),
    lwd = c(3,1,2,3), cex = 1.2, bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# fit the dynocc model in unmarked
# ''''''''''''''''''''''''''''''''
library(unmarked)
yy <- matrix(y2, nrow = nsites, ncol = nsurveys * nyears)
str(yy)
# int [1:100, 1:24] 0 NA NA 0 0 NA NA 0 1 0 ...

# Package and summarize data set
summary(umf <- unmarkedMultFrame(y = yy, numPrimary = nyears))
# unmarkedFrame Object

# 100 sites
# Maximum number of observations per site: 24
# Mean number of observations per site: 18.86
# Number of primary survey periods: 12
# Number of secondary survey periods: 2
# Sites with at least one detection: 70

# Tabulation of y observations:
# 0 1 <NA>
# 1651 235 514

# Fit dynamic occupancy model and look at estimates
summary(fm <- colext(psiformula = ~1, # First-year occupancy
    gammaformula = ~ 1,               # Colonization
    epsilonformula = ~ 1,             # Extinction
    pformula = ~ 1,                   # Detection
    data = umf))
# Call:
# colext(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1,
#   pformula = ~1, data = umf)

# Initial (logit-scale):
# Estimate SE z P(>|z|)
# 1.07 0.435 2.47 0.0135

# Colonization (logit-scale):
# Estimate SE z P(>|z|)
# -3.76 0.93 -4.04 5.3e-05

# Extinction (logit-scale):
# Estimate SE z P(>|z|)
# -2.4 0.304 -7.88 3.24e-15

# Detection (logit-scale):
# Estimate SE z P(>|z|)
# -1.16 0.127 -9.1 9.18e-20

# AIC: 1167.49
# Number of sites: 100
# optim convergence code: 0
# optim iterations: 47
# Bootstrap iterations: 0

# Backtransform estimates to probability scale
backTransform(fm, type = "psi")     # First-year occupancy
backTransform(fm, type = "col")     # Colonization probability
backTransform(fm, type = "ext")     # Extinction probability
backTransform(fm, type = "det")     # Detection probability
# Estimate     SE LinComb (Intercept)
#    0.745 0.0825    1.07           1

# Estimate     SE LinComb (Intercept)
#   0.0228 0.0207   -3.76           1

# Estimate     SE LinComb (Intercept)
#   0.0835 0.0233    -2.4           1

# Estimate     SE LinComb (Intercept)
#    0.239 0.0231   -1.16           1

# ~~~~~~~~~~~~~~~~~
# MLEs <- c(
# `psi(Int)` = backTransform(fm, type = "psi")@estimate, # First-year occupancy
# `col(Int)` = backTransform(fm, type = "col")@estimate, # Colonization probability
# `ext(Int)` = backTransform(fm, type = "ext")@estimate, # Extinction probability
# `det(Int)` = backTransform(fm, type = "det")@estimate) # Detection probability
# ~~~~~~~~~~~~~~~~~~~~~

# For Null model point estimates can simply do this
(MLEs <- plogis(coef(fm)))

# Get 95% CIs on probability scale and print them along with MLEs
( MLEandCI <- cbind(MLEs, rbind(plogis(confint(fm, type = "psi")),
    plogis(confint(fm, type = "col")),
    plogis(confint(fm, type = "ext")),
    plogis(confint(fm, type = "det"))) ))
#                MLEs       0.025     0.975
# psi(Int) 0.74517648 0.555131656 0.8726586
# col(Int) 0.02277496 0.003751754 0.1260508
# ext(Int) 0.08351979 0.047826505 0.1418819
# p(Int)   0.23907422 0.196680790 0.2873369

# fit the same model in BUGS
# ''''''''''''''''''''''''''
# Bundle and summarize data set
str(bdata <- list(y = y2, nsites = dim(y2)[1], nsurveys = dim(y2)[2],
    nyears = dim(y2)[3]))
# List of 4
# $ y        : int [1:100, 1:2, 1:12] 0 NA NA 0 0 NA NA 0 1 0 ...
# $ nsites   : int 100
# $ nsurveys : int 2
# $ nyears   : int 12

# Specify model in BUGS language
cat(file = "dynocc.txt","
model {

  # Priors
  psi1 ~ dunif(0, 1)
  phi ~ dunif(0, 1)
  gamma ~ dunif(0, 1)
  p ~ dunif(0, 1)

  # Likelihood
  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){ # Loop over nsites sites
    # Initial conditions of system
    z[i,1] ~ dbern(psi1) # Presence/absence at start of study
    # State transitions
    for (t in 2:nyears){ # Loop over nyears years
      z[i,t] ~ dbern(z[i,t-1] * phi + (1-z[i,t-1]) * gamma)
    }
  }

  # Observation model
  for (i in 1:nsites){
    for (j in 1:nsurveys){
      for (t in 1:nyears){
        y[i,j,t] ~ dbern(z[i,t] * p)
      }
    }
  }

  # Derived parameters
  eps <- 1 - phi # Extinction probability
  psi[1] <- psi1 # Population occupancy
  for (t in 2:nyears){
    psi[t] <- psi[t-1] * phi + (1-psi[t-1]) * gamma
  }
}
")

# Initial values
zst <- zobs # Take observed presence/absence as inits
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi1", "phi", "eps", "gamma", "p", "psi") # Could add 'z'

# MCMC settings
na <- 1000 ; ni <- 25000 ; nt <- 10 ; nb <- 5000 ; nc <- 3

# Call JAGS (ART 2 min), check convergence and summarize posteriors
library(jagsUI)
out1 <- jags(bdata, inits, params, "dynocc.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
print(out1, 3) # not shown

# ~~~~~~~ code for the table ~~~~~~~~~~~~~~~~~~
truth <- c("psi1" = psi1, "col" = gamma, "ext" = 1-phi, "det" = p)
Bayes <- out1$summary[c(1, 15, 14, 16), c(1,3,7)]
round(cbind(truth, MLEandCI, Bayes), 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      truth  MLEs 0.025 0.975  mean  2.5% 97.5%
# psi1  0.70 0.745 0.555 0.873 0.742 0.581 0.885
# col   0.05 0.023 0.004 0.126 0.086 0.050 0.139
# ext   0.10 0.084 0.048 0.142 0.032 0.002 0.079
# det   0.25 0.239 0.197 0.287 0.234 0.198 0.277


# ~~~~~~ code for figure 4.3 ~~~~~~~~~~~~~~~~~~~
# Plot true and estimated population occupancy probability
op <- par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
off <- 0.2       # Graphical offset
plot(1:nyears, psi, main = "True and estimated occupancy probability", xlab = "Year",
    ylab = "Occupancy", xlim = c(0.5, 12.5), ylim = c(0.3, 0.9), pch = 15, type = 'b',
    col = 'red', frame = F, cex = 2, lwd = 3)
points(1:nyears-off, fm@projected[2,,1], type = 'b', pch = 1, cex = 2)
points(1:nyears, out1$mean$psi, type = 'b', col = "black", pch = 16, cex = 2)
segments(1:nyears, out1$q2.5$psi, 1:nyears, out1$q97.5$psi, col = "black", lwd = 1)
legend(8, 0.9, c('Truth', 'MLEs', 'Posterior means'), col=c("red", "black", "black"),
    pch = c(15, 1, 16), cex = 1.5, bty = 'n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
