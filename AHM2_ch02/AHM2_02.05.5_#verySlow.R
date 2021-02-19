#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# Approximate execution time for this code: 1.9 days
#   excluding the 2 blocks wrapped in if(FALSE){} which take
#   9 hrs and 5 days respectively.

library(jagsUI)
library(unmarked)
library(AHMbook)

# 2.5 Dynamic N-mixture model of Dail-Madsen
# ==========================================

# 2.5.5 Fitting some complex covariate effects in the dynamics
# ------------------------------------------------------------

# Test with same settings as before, that is, without covariates
set.seed(2017, kind = "L'Ec")
data <- simDM(nsites = 100, nsurveys = 3, nyears = 5, mean.lambda = 4,
    mean.phi = 0.8, mean.gamma.rel = 0.5, mean.p = 0.7, beta.lam = 0,
    beta.gamma = 0, beta.phi = 0, beta.p = 0)

# Fit model in unmarked with Dynamics = "autoreg" ART ~ 30 secs
summary(umf <- unmarkedFramePCO(y = data$yy, numPrimary = data$nyears))
system.time(summary(fm3 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1,
    p = ~1, data = umf, dynamics = "autoreg", K = 100,
    control = list(trace = TRUE, REPORT = 1))))  # 36 secs

# Check with higher value of K: find that K = 100 is OK. ART ~ 1.5 min
system.time(summary(fm3x <- pcountOpen(lam = ~1, gam = ~1, omega = ~1,
    p = ~1, data = umf, dynamics = "autoreg", K = 150,
    control = list(trace = TRUE, REPORT = 1)))) # not shown # 2 mins

# Backtransform parameter estimates and compare with truth
(lam <- exp(coef(fm3, type = "lambda")))  # Initial abundance
(gam <- exp(coef(fm3, type = "gamma")))   # Recruitment
(om <- plogis(coef(fm3, type = "omega"))) # Apparent survival
(p <- plogis(coef(fm3, type = "det")))    # Detection
str(data[4:7]) # Compare with The Truth

# lam(Int)
# 4.088409

# gamAR(Int)
# 0.5014593

# omega(Int)
# 0.8047996

# p(Int)
# 0.698113

# List of 4
# $ mean.lambda : num 4
# $ mean.gamma.rel : num 0.5
# $ mean.phi : num 0.8
# $ mean.p : num 0.7

# Execute function (plot output omitted)
set.seed(24, kind = "Mersenne-Twister")
str(data <- simDM(nsites = 50, nsurveys = 3, nyears = 5, mean.lambda = 4,
    mean.phi = 0.8, mean.gamma.rel = 0.5, mean.p = 0.7, beta.lam = 1,
    beta.gamma = 1, beta.phi = -1, beta.p = -1) ) # Details omitted

# Prepare data for unmarked
summary(umf <- unmarkedFramePCO(y = data$yy, numPrimary = data$nyears,
    siteCovs = data.frame(cov.lam = data$cov.lam, cov.gamma = data$cov.gamma,
    cov.phi = data$cov.phi), obsCovs = list(cov.p = data$ccov.p)))

# Fit Null model, switch off computation of SEs (ART 300 sec)
system.time(summary(fm4.1 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p= ~1,
    data = umf, dynamics = "autoreg", se = FALSE,
    control = list(trace = TRUE, REPORT = 1))))

# Add in site covariate on lambda (ART 120 sec, no SEs)
tmp <- coef(fm4.1)
inits <- c(tmp[1], 0, tmp[2:4]) # choose 0 for the new param in model
system.time(summary(fm4.2 <- pcountOpen(lam = ~ cov.lam, gam = ~1,
    omega = ~1, p = ~1, data = umf, dynamics = "autoreg", se = FALSE,
    control = list(trace = TRUE, REPORT = 1), starts = inits)))  # 1 min

# Add in observation covariate on p. ART ~ 2.5 mins
tmp <- coef(fm4.2)
inits <- c(tmp, 0)
system.time(summary(fm4.3 <- pcountOpen(lam = ~ cov.lam, gam = ~1,
    omega = ~1, p = ~ cov.p, data = umf, dynamics = "autoreg", se = FALSE,
    control = list(trace = TRUE, REPORT = 1), starts = inits) ))  # 80 secs

# Add in site covariate on gamma. ART ~ 14 hrs
tmp <- coef(fm4.3)
inits <- c(tmp[1:3], 0, tmp[4:6])
system.time(summary(fm4.4 <- pcountOpen(lam = ~ cov.lam,
    gam = ~cov.gamma, omega = ~1, p = ~ cov.p, data = umf, dynamics = "autoreg",
    se = FALSE, control = list(trace = TRUE, REPORT = 1), starts = inits) )) # 7.2 hrs
# ~~~ save the work so far ~~~~~
save.image("AHM2_02.05.5a.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add in final covariate, the site covariate on phi/omega. ART ~ 2.6 hrs
tmp <- coef(fm4.4)
inits <- c(tmp[1:5], 0, tmp[6:7])
system.time(summary(fm4.5 <- pcountOpen(lam = ~ cov.lam, gam = ~cov.gamma,
    omega = ~ cov.phi, p = ~ cov.p, data = umf,
    dynamics = "autoreg", se = FALSE, control = list(trace = TRUE, REPORT = 1),
    starts = inits) ))   # 9.6 hrs / 14 hrs
# ~~~ save the work so far ~~~~~
save.image("AHM2_02.05.5b.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Put the fitted models in a "fitList"
fms <- fitList(
"lam(.)gamma(.)phi(0)p(.)" = fm4.1,
"lam(cov)gamma(.)phi(0)p(.)" = fm4.2,
"lam(cov)gamma(.)phi(0)p(cov)" = fm4.3,
"lam(cov)gamma(cov)phi(0)p(cov)" = fm4.4,
"lam(cov)gamma(cov)phi(cov)p(cov)" = fm4.5)

# Rank them by AIC (singular Hessian message is due to SE=F option)
(ms <- modSel(fms))
#                                  nPars     AIC  delta    AICwt cumltvWt
# lam(cov)gamma(cov)phi(cov)p(cov)     8 2666.10   0.00  1.0e+00        1
# lam(cov)gamma(cov)phi(0)p(cov)       7 2688.62  22.52  1.3e-05        1
# lam(cov)gamma(.)phi(0)p(cov)         6 2813.65 147.55  9.1e-33        1
# lam(cov)gamma(.)phi(0)p(.)           5 3207.85 541.75 2.3e-118        1
# lam(.)gamma(.)phi(0)p(.)             4 3231.83 565.73 1.4e-123        1

# Coefficients for each model and Table With Everything (not printed)
print(coef(ms), digits=3)
toExport <- as(ms, "data.frame")
str(toExport) # output not printed

# ~~~~ These models can be skipped to save time, the results are not used subsequently ~~~~
# ~~~~ They take around 10 hours to run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(FALSE) {
  # Fit the basic model with no covariates
  system.time(summary(fm4.1.se <- pcountOpen(lam = ~1, gam = ~1,
      omega = ~1, p = ~1, data = umf, dynamics = "autoreg",
      control = list(trace = TRUE, REPORT = 1))))  # 3 mins

  # Add in site covariate on lambda (ART 179 sec, no SEs)
  tmp <- coef(fm4.1.se)
  inits <- c(tmp[1], 0, tmp[2:4])
  system.time(summary(fm4.2.se <- pcountOpen(lam = ~ cov.lam, gam = ~1,
      omega = ~1, p = ~1, data = umf, dynamics = "autoreg",
      control = list(trace = TRUE, REPORT = 1), starts = inits)))

  # Add in observation covariate on p
  tmp <- coef(fm4.2.se)
  inits <- c(tmp, 0)
  system.time(summary(fm4.3.se <- pcountOpen(lam = ~ cov.lam, gam = ~1,
      omega = ~1, p = ~ cov.p, data = umf, dynamics = "autoreg",
      control = list(trace = TRUE, REPORT = 1), starts = inits) ))

  # Add in site covariate on gamma. Starts need a bit of help....
  tmp <- coef(fm4.3.se)
  inits <- c(tmp[1:2], -0.6, 0.9, qlogis(0.8), tmp[5:6])
  system.time(summary(fm4.4.se <- pcountOpen(lam = ~ cov.lam,
      gam = ~ cov.gamma, omega = ~1, p= ~ cov.p, data = umf, dynamics = "autoreg",
      control = list(trace = TRUE, REPORT = 1), starts = inits) ))  # 6.6 hrs
  # ~~~ save the work so far ~~~~
  save.image("AHM2_02.05.5c.RData")
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}

# Add in final covariate, the site covariate on phi/omega
# tmp <- coef(fm4.4.se)
tmp <- coef(fm4.4)  # use this if you didn't do fm4.4.se
inits <- c(tmp[1:5], 0, tmp[6:7])
system.time(summary(fm4.5.se <- pcountOpen(lam = ~ cov.lam,
    gam = ~ cov.gamma, omega = ~ cov.phi, p= ~ cov.p, data = umf,
    dynamics = "autoreg", control = list(trace = TRUE, REPORT = 1), starts = inits) )) # 7.7 / 10.5 / 18 hrs
# ~~~ save the work so far ~~~~
save.image("AHM2_02.05.5d.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# To get a sense of the time comparison with and without SE computation, we ran model fm4.5 both
# ways (here note the use of estimates from model fm4.5 as starting values in both cases):
# ~~~ These take several days, do not include in routine tests ~~~~~~~~~~~~~~~~~~~~~~
if(FALSE) {
# Model 4.5 without SEs
tmp <- coef(fm4.4)
inits <- c(tmp[1:5], 0, tmp[6:7])
system.time(summary(fm4.5 <- pcountOpen(lam = ~ cov.lam,
    gam = ~ cov.gamma, omega = ~ cov.phi, p = ~ cov.p, data = umf,
    dynamics = "autoreg", K = max(data$yy) + 100, se = FALSE,
    control = list(trace = TRUE, REPORT = 1), starts = inits) ))

# user system elapsed
# 177435.19 1.32 177449.20 # 2.1 days

# Model 4.5 with SEs
system.time(summary(fm4.5.se <- pcountOpen(lam = ~ cov.lam, gam = ~ cov.gamma, omega = ~
cov.phi, p = ~ cov.p, data = umf, dynamics = "autoreg", K = max(data$yy) + 100, control =
list(trace = TRUE, REPORT = 1), starts = inits) ))

# user system elapsed
# 240161.03 1.86 240178.64 # 2.8 days
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fm4.5.se
# Abundance:
# Estimate SE z P(>|z|)
# (Intercept) 1.439 0.0754 19.09 2.99e-81
# cov.lam 0.632 0.1350 4.68 2.83e-06

# Recruitment:
# Estimate SE z P(>|z|)
# (Intercept) -0.552 0.125 -4.41 1.05e-05
# cov.gamma 0.880 0.104 8.48 2.25e-17

# Apparent Survival:
# Estimate SE z P(>|z|)
# (Intercept) 0.943 0.333 2.83 0.00470
# cov.phi -0.793 0.248 -3.20 0.00137

# Detection:
# Estimate SE z P(>|z|)
# (Intercept) 0.841 0.0807 10.4 2.12e-25
# cov.p -1.012 0.0668 -15.2 6.49e-52

# Running it in JAGS
# ''''''''''''''''''
# Bundle data
str(bdata <- list(C = data$y, nsites = dim(data$y)[1],
    nsurveys = dim(data$y)[3], nyears = dim(data$y)[2], cov.lam = data$cov.lam,
    cov.gamma = data$cov.gamma, cov.phi = data$cov.phi, cov.p = data$cov.p))

# Specify model in BUGS language
cat(file = "DM2.txt","
model {
  # Priors
  alpha.lam <- log(mean.lambda)
  mean.lambda ~ dunif(0, 100) # Population growth rate
  alpha.phi ~ dnorm(0,0.1) # Apparent survival (or omega)
  alpha.gamma ~ dnorm(0,0.1) # Per-capita recruitment rate
  alpha.p ~ dnorm(0,0.1) # Detection probability
  beta.lam ~ dnorm(0,0.1) # Coefs of 4 covariates
  beta.phi ~ dnorm(0,0.1)
  beta.gamma ~ dnorm(0,0.1)
  beta.p ~ dnorm(0,0.1)

  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha.lam + beta.lam * cov.lam[i]
    logit(phi[i]) <- alpha.phi + beta.phi*cov.phi[i]
    log(gamma[i]) <- alpha.gamma + beta.gamma*cov.gamma[i]

    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi[i], N[i,t])
      # R[i,t+1] ~ dpois(gamma[i]) # 'absolute' recruitment = 'constant'
      tmp[i,t] <- N[i,t] * gamma[i] # per-capita recruitment = 'autoreg'
      R[i,t+1] ~ dpois(tmp[i,t]) # per-capita recruitment = 'autoreg'
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    }

    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        logit(p[i,t,j]) <- alpha.p + beta.p*cov.p[i,t,j]
        C[i,t,j] ~ dbin(p[i,t,j], N[i,t])
      }
    }
  }
  # Derived quantities
  mean.phi <- ilogit(alpha.phi)
  mean.gamma <- exp(alpha.gamma)
  mean.p <- ilogit(alpha.p)
}
")

# Initial values that usually seem to work
R1 <- apply(data$y, c(1,2), max) + 10 # Use observed max. counts + 10
                                      # as inits for recruitment
R1[,1] <- NA
Nst <- apply(data$y, c(1,2), max)+2
Nst[,2:ncol(Nst)] <- NA
inits <- function(){ list(N = Nst, R = R1, beta.lam = 0, beta.phi = 0,
    beta.gamma = 0, beta.p = 0) }

# Parameters monitored
# could also monitor the latent variables: "N", "R", "S"
params <- c("mean.lambda", "mean.phi", "mean.gamma", "mean.p",
    "beta.lam", "beta.phi", "beta.gamma", "beta.p", "alpha.lam", "alpha.phi",
    "alpha.gamma", "alpha.p")

# MCMC settings
# na <- 1000 ; ni <- 500000 ; nt <- 250 ; nb <- 250000 ; nc <- 3  # 35 mins
na <- 1000 ; ni <- 50000 ; nt <- 25 ; nb <- 25000 ; nc <- 3 # ~~~ for testing, 4 mins

# Call JAGS (ART 33 min), check convergence and summarize posteriors
out4 <- jags(bdata, inits, params, "DM2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
op <- par(mfrow = c(2, 3)) ; traceplot(out4)
par(op)
print(out4, 2) # not printed

# Compare estimates with truth
cbind('truth' = unlist(data[4:11]), round(out4$summary[1:8, c(1,3,7)],3))
#                truth   mean   2.5%  97.5%
# mean.lambda      4.0  4.247  3.645  4.857
# mean.gamma.rel   0.5  0.706  0.513  0.837
# mean.phi         0.8  0.590  0.451  0.793
# mean.p           0.7  0.695  0.659  0.727
# beta.lam         1.0  0.632  0.364  0.891
# beta.gamma       1.0 -0.813 -1.357 -0.337
# beta.phi        -1.0  0.873  0.648  1.085
# beta.p          -1.0 -1.005 -1.140 -0.870

# ~~~ save the work so far ~~~~~
save.image("AHM2_02.05.5.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
