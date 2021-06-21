#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 20 mins
# Run time with the full number of iterations: 4.5 hrs

library(AHMbook)
library(unmarked)
options(stringsAsFactors=TRUE)

# 4.5 Simulation and analysis of a time-dependent data set with unmarked and BUGS
# ===============================================================================

set.seed(1)
str(data <- simDynocc(nsites = 250, nyears = 10, nsurveys = 3, mean.psi1 = 0.6,
    range.phi = c(0.5, 1), range.gamma = c(0, 0.5), range.p = c(0.1, 0.5)))
    # produces figure 4.4

library(unmarked)
str(yy <- matrix(data$y, data$nsites, data$nsurveys * data$nyears) )
# int [1:250, 1:30] 0 0 0 1 0 0 0 1 1 0 ...

# Create a matrix indicating the year each site is surveyed
year <- matrix(c('01','02','03','04','05','06','07','08','09','10'), nrow = nrow(yy),
    ncol = data$nyears, byrow=TRUE)
head(year) # not shown

xxp <- matrix(data$Xp, data$nsites, data$nsurveys * data$nyears) # Xp in 2D
summary( umf <- unmarkedMultFrame(y = yy, siteCovs = data.frame(Xpsi1 = data$Xpsi1),
    yearlySiteCovs = list(year = year, Xphi = data$Xphi, Xgamma = data$Xgamma),
    obsCovs = list(Xp = xxp), numPrimary = data$nyears) )
# unmarkedFrame Object

# 250 sites
# Maximum number of observations per site: 30
# Mean number of observations per site: 30
# Number of primary survey periods: 10
# Number of secondary survey periods: 3
# Sites with at least one detection: 224

# Tabulation of y observations:
# 0 1
# 6492 1008

# Site-level covariates:
# Xpsi1
# Min. :-1.94769
# 1st Qu.:-0.86551
# Median :-0.01341
# Mean : 0.03930
# 3rd Qu.: 0.92469
# Max. : 1.97074

# Observation-level covariates:
# Xp
# Min. :-1.999574
# 1st Qu.:-0.981910
# Median : 0.008498
# Mean : 0.009213
# 3rd Qu.: 1.003295
# Max. : 1.999421

# Yearly-site-level covariates:
#     year        Xphi               Xgamma
# 01 : 250   Min. :-1.99758     Min. :-1.999199
# 02 : 250   1st Qu.:-1.07270   1st Qu.:-0.983174
# 03 : 250   Median :-0.08815   Median :-0.008474
# 04 : 250   Mean :-0.03517     Mean : 0.008517
# 05 : 250   3rd Qu.: 1.01944   3rd Qu.: 1.092480
# 06 : 250   Max. : 1.99972     Max. : 1.996779
# (Other):1000

fm1 <- colext(psiformula = ~1, # First-year occupancy
    gammaformula = ~ 1,        # Colonization
    epsilonformula = ~ 1,      # Extinction
    pformula = ~ 1,            # Detection
    data = umf)
system.time(summary(fm2 <- colext(psiformula = ~1, gammaformula = ~year-1,
    epsilonformula = ~1, pformula = ~1, data = umf,
    control = list(trace = TRUE, REPORT = 5), se = TRUE)))
system.time(summary(fm3 <- colext(psiformula = ~1, gammaformula = ~1,
    epsilonformula = ~year-1, pformula = ~1, data = umf,
    control = list(trace = TRUE, REPORT = 5), se = TRUE)))
system.time(summary(fm4 <- colext(psiformula = ~1, gammaformula = ~1,
    epsilonformula = ~1, pformula = ~ year-1, data = umf,
    control = list(trace = TRUE, REPORT = 5), se = TRUE)))
system.time(summary(fm5 <- colext(psiformula = ~1, gammaformula = ~year-1,
    epsilonformula = ~year-1, pformula = ~1, data = umf,
    control = list(trace = TRUE, REPORT = 5), se = TRUE)))
system.time(summary(fm6 <- colext(psiformula = ~1, gammaformula = ~year-1,
    epsilonformula = ~1, pformula = ~year-1, data = umf,
    control = list(trace = TRUE, REPORT = 5, maxit = 500), se = TRUE)))
system.time(summary(fm7 <- colext(psiformula = ~1, gammaformula = ~1,
    epsilonformula = ~year-1, pformula = ~year-1, data = umf,
    control = list(trace = TRUE, REPORT = 5), se = TRUE)))
system.time(summary(fm8 <- colext(psiformula = ~1, gammaformula = ~year-1,
    epsilonformula = ~year-1, pformula = ~year-1, data = umf,
    control = list(trace = TRUE, REPORT = 5), se = TRUE)))

# Generate a fit list and compare the models using AIC (and check out some German)
models <- fitList(
    'psi(.)gam(.)eps(.)p(.)' = fm1,
    'psi(.)gam(Y)eps(.)p(.)' = fm2,
    'psi(.)gam(.)eps(Y)p(.)' = fm3,
    'psi(.)gam(.)eps(.)p(Y)' = fm4,
    'psi(.)gam(Y)eps(Y)p(.)' = fm5,
    'psi(.)gam(Y)eps(.)p(Y)' = fm6,
    'psi(.)gam(.)eps(Y)p(Y)' = fm7,
    'psi(.)gam(Y)eps(Y)p(Y)' = fm8)

(ms <- modSel(models))
#                        nPars     AIC  delta   AICwt cumltvWt
# psi(.)gam(Y)eps(Y)p(Y)    29 5365.38   0.00 1.0e+00        1
# psi(.)gam(Y)eps(.)p(Y)    21 5385.80  20.42 3.7e-05        1
# psi(.)gam(.)eps(Y)p(Y)    21 5388.12  22.73 1.2e-05        1
# psi(.)gam(Y)eps(Y)p(.)    20 5408.31  42.92 4.8e-10        1
# psi(.)gam(.)eps(.)p(Y)    13 5413.76  48.38 3.1e-11        1
# psi(.)gam(Y)eps(.)p(.)    12 5450.96  85.58 2.6e-19        1
# psi(.)gam(.)eps(Y)p(.)    12 5463.51  98.12 4.9e-22        1
# psi(.)gam(.)eps(.)p(.)     4 5496.68 131.29 3.1e-29        1

# Warnmeldung:
# In sqrt(diag(vcov(x, altNames = TRUE))) : NaNs wurden erzeugt

# Bundle and summarize data
str(bdata <- list(y = data$y, nsites = data$nsites, nsurveys = data$nsurveys,
    nyears = data$nyears))

# Specify model in BUGS language
cat(file = "dynocc.txt","         # overwrite previous file
  model {

  # Specify priors
  psi1 ~ dunif(0, 1)
  for (t in 1:(nyears-1)){
    phi[t] ~ dunif(0, 1)
    gamma[t] ~ dunif(0, 1)
    p[t] ~ dunif(0, 1)
  }
  p[nyears] ~ dunif(0, 1)

  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsites){
    z[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      z[i,t] ~ dbern(z[i,t-1]*phi[t-1] + (1-z[i,t-1])*gamma[t-1])
    }
  }

  # Observation model
  for (i in 1:nsites){
    for (j in 1:nsurveys){
      for (t in 1:nyears){
        y[i,j,t] ~ dbern(z[i,t]*p[t])
      }
    }
  }

  # Derived parameters
  # Sample and population occupancy, growth rate and turnover
  # Also, logit-scale params for direct comparison with unmarked
  lpsi1 <- logit(psi1)
  lp[1] <- logit(p[1])
  psi[1] <- psi1
  n.occ[1] <- sum(z[1:nsites,1])
  for (t in 2:nyears){
    psi[t] <- psi[t-1]*phi[t-1] + (1-psi[t-1])*gamma[t-1]
    n.occ[t] <- sum(z[1:nsites,t])
    growthr[t-1] <- psi[t]/psi[t-1]
    turnover[t-1] <- (1 - psi[t-1]) * gamma[t-1]/psi[t]
    lgamma[t-1] <- logit(gamma[t-1])
    leps[t-1] <- logit(1-phi[t-1])
    lp[t] <- logit(p[t])
  }
}
")

# Initial values
zst <- apply(data$y, c(1, 3), max)     # Obs. occurrence as inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover",
    "lpsi1", "lgamma", "leps", "lp")   # could add 'z'

# MCMC settings
na <- 1000 ; ni <- 20000 ; nt <- 10 ; nb <- 10000 ; nc <- 3

# Call JAGS (ART 3 min), check convergence and summarize posteriors
library("jagsUI")
out2 <- jags(bdata, inits, params, "dynocc.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out2)
print(out2, dig = 2) # not shown

# system.time(fm8 <- nonparboot(fm8, B = 500) ) # Takes about 3 hours
system.time(fm8 <- nonparboot(fm8, B = 50) )    # ~~~~ for testing
cbind(psi = data$mean.psi, smoothed = smoothed(fm8)[2,],
    SE = fm8@smoothed.mean.bsse[2,])            # Finite sample occupancy
cbind(psi = data$mean.psi, projected = projected(fm8)[2,],
    SE = fm8@projected.mean.bsse[2,])           # Population occupancy

round(cbind(psi = data$mean.psi, ML.estimates = projected(fm8)[2,],
    ML.ASE = fm8@projected.mean.bsse[2,], Bayesian.estimates = out2$summary[1:10,
    c(1:2)]), 4)
#           psi ML.estimates ML.ASE   mean     sd
# Year1  0.6000       0.6272 0.0536 0.5790 0.0561
# Year2  0.6192       0.6262 0.0406 0.6007 0.0420
# Year3  0.6371       0.5905 0.0890 0.5272 0.0748
# Year4  0.5209       0.4812 0.0827 0.5028 0.0712
# Year5  0.4750       0.3359 0.0603 0.3471 0.0528
# Year6  0.3698       0.3959 0.0383 0.4007 0.0386
# Year7  0.4939       0.4541 0.0635 0.4553 0.0601
# Year8  0.4244       0.4120 0.0576 0.4304 0.0580
# Year9  0.4245       0.4000 0.0343 0.3989 0.0333
# Year10 0.2922       0.3052 0.0459 0.3172 0.0479

# ~~~~~~~~ code for figure 4.5 ~~~~~~~~~~~~~~~~~~~~
off <- 0.05         # graphical offset
plot(1:data$nyears, data$psi[1,], type = "l", xlab = "Year", ylab = "Occupancy probability",
    col = "red", xlim = c(0,data$nyears+1), ylim = c(0,1), lwd = 2, lty = 1,
    frame.plot = FALSE, las = 1)
lines(1:data$nyears, data$psi.app, type = "l", col = "black", lwd = 2)
points(1:data$nyears-off, out2$mean$psi, type = "l", col = "blue", lwd = 2)
segments(1:data$nyears-off, out2$mean$psi-out2$sd$psi, 1:data$nyears-off,
    out2$mean$psi+out2$sd$psi, col = "blue", lwd = 1)
points(1:data$nyears+off, projected(fm8)[2,], type = "l", col = "green", lwd = 2)
segments(1: data$nyears+off, projected(fm8)[2,] - fm8@projected.mean.bsse[2,],
    1: data$nyears+off, projected(fm8)[2,] + fm8@projected.mean.bsse[2,], col = "green", lwd = 1)
legend('topright', legend = c('True psi', 'Obs. psi', 'Estimated psi (unmarked)',
    'Estimated psi (JAGS)'), col = c('red', 'black', 'green', 'blue'), lty = 1,
    lwd = 3, bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nd1 <- data.frame(year = c('01','02','03','04','05','06','07','08','09'))
nd2 <- data.frame(year = c('01','02','03','04','05','06','07','08','09','10'))

E.ext <- predict(fm8, type = 'ext', newdata = nd1)
E.col <- predict(fm8, type = 'col', newdata = nd1)
E.det <- predict(fm8, type = 'det', newdata = nd2)

# ~~~~~~~~ code for figure 4.6 ~~~~~~~~~~~~~~~~~
# Plot for extinction probability
op <- par(mfrow=c(3,1), mar=c(5, 4, 3, 3), cex = 0.8)
with(E.ext, { # simplify requesting columns of the data frame returned by predict
   plot(1:9, Predicted, pch=16, xaxt='n', xlab='',
       main=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylab = "", ylim=c(0,1), col="green")
   axis(1, at=1:9, labels=nd1$year[1:9])
   segments(1:9, lower, 1:9, upper, col="green")
   points((1:9)-0.1, 1-data$mean.phi, col="red", lwd = 1, pch=16)
   points((1:9)+0.1, 1-out2$summary[11:19,1], col="blue", lwd = 1, pch=16)
   segments(((1:9)+0.1), (1-out2$summary[11:19,3]), ((1:9)+0.1),
      (1-out2$summary[11:19,7]), col = "blue", lwd = 1)
})

# Plot for colonization probability
with(E.col, {
   plot(1:9, Predicted, pch=16, xaxt='n', xlab='',
       main=expression(paste('Colonization probability ( ', gamma, ' )')),
       ylab = "", ylim=c(0,1), col="green")
   axis(1, at=1:9, labels=nd1$year[1:9])
   segments(1:9, lower, 1:9, upper, col="green")
   points((1:9)-0.1, data$mean.gamma, col="red", lwd = 1, pch=16)
   points((1:9)+0.1, out2$summary[20:28,1], col="blue", lwd = 1, pch=16)
   segments(((1:9)+0.1), out2$summary[20:28,3], ((1:9)+0.1),
      out2$summary[20:28,7], col =   "blue", lwd = 1)
})

# Plot for detection probability: note 10 years
with(E.det, {
   plot(1:10, Predicted, pch=16, xaxt='n', xlab='Year',
       main=expression(paste('Detection probability ( ', p, ' )')),
       ylab = "", ylim=c(0,1), col="green")
   axis(1, at=1:10, labels=nd2$year[1:10])
   segments(1:10, lower, 1:10, upper, col="green")
   points((1:10)-0.1, data$mean.p, col="red", lwd = 1, pch=16)
   points((1:10)+0.1, out2$summary[29:38,1], col="blue", lwd = 1, pch=16)
   segments(((1:10)+0.1), out2$summary[29:38,3], ((1:10)+0.1),
      out2$summary[29:38,7], col = "blue", lwd = 1)
   legend(1, 1, c('Truth', 'MLEs', 'Posterior means'), col=c("red", "green",
      "blue"), pch=16, cex=0.8)
})
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~