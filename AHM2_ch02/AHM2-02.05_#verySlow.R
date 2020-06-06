#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-01-09

library(jagsUI)
library(unmarked)
library(AHMbook)

# 2.5 Dynamic N-mixture model of Dail-Madsen
# ==========================================

# 2.5.1 A Dail-Madsen simulator
# -----------------------------

simDM0 <- function(nsites = 50, nsurveys = 3, nyears = 5, lambda = 4,
    phi = 0.8, gamma = 1.5, p = 0.7){
  ## Simulation for multiple-visit data (from pcountOpen help file)
  ## No covariates, constant time intervals between primary periods
  # nsites: Number of sites
  # nsurveys: Number of rep. (secondary) samples within period of closure
  # nyears: Number of primary samples (= period of closure):
  # years, seasons etc.
  # lambda: Initial expected abundance
  # phi, gamma: apparent survival and recruitment rates, respectively
  # p: detection probability
  y <- array(NA, dim = c(nsites, nyears, nsurveys))
  N <- matrix(NA, nsites, nyears)
  S <- R <- matrix(NA, nsites, nyears-1)
  N[,1] <- rpois(nsites, lambda) # Initial state
  for(t in 1:(nyears-1)) { # State dynamics
    S[,t] <- rbinom(nsites, N[,t], phi) # Survival process
    R[,t] <- rpois(nsites, gamma) # Recruitment process
    N[,t+1] <- S[,t] + R[,t]
  }
  for(j in 1:nsurveys){ # Observation process
    y[,,j] <- rbinom(nsites*nyears, N, p)
  }
  # Put observed data into two dimensions
  yy <- array(NA, dim = c(nsites, nsurveys*nyears))
  for(t in 1:nyears){
    yy[,(nsurveys * t-(nsurveys-1)):(nsurveys*t)] <- y[,t,]
  }
  return(list(nsites = nsites, nsurveys = nsurveys, nyears = nyears,
      lambda = lambda, phi = phi, gamma = gamma, p = p, N = N, S = S, R = R, y = y, yy = yy))
}

# Execute function
# set.seed(2017, kind = "L’Ecuyer")
set.seed(2017, kind = "L'Ecuyer")
str(data <- simDM0(nsites = 50, nsurveys = 3, nyears = 5, lambda = 4,
    phi = 0.8, gamma = 1.5, p = 0.7))


# 2.5.2 Fitting the DM model in BUGS
# ----------------------------------

# Bundle data set
str(bdata <- list(C = data$y, nsites = dim(data$y)[1], nsurveys = dim(data$y)[3],
  nyears = dim(data$y)[2]))
# List of 4
# $ C: int [1:50, 1:5, 1:3] 3 4 4 1 2 5 3 5 7 5 ...
# $ nsites: int 50
# $ nsurveys: int 3
# $ nyears: int 5

# Specify model in BUGS language
cat(file = "DM1.txt","
model {
  # Priors
  lambda ~ dunif(0, 100) # Population growth rate
  phi ~ dunif(0, 1) # apparent survival (omega in paper/unmarked)
  gamma ~ dunif(0, 5) # per-capita recruitment rate
  p ~ dunif(0, 1) # Detection probability
  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t]) # Survival process
      # R[i,t+1] ~ dpois(gamma) # 'absolute’ recruitment = 'constant’
      R[i,t+1] ~ dpois(N[i,t] * gamma) # per-capita recruitment = 'autoreg’
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    } # end t
    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,t,j] ~ dbin(p, N[i,t])
      } # end j
    } # end t
  } # end i
}
")

# Initial values
Nst <- apply(data$y, c(1,2), max) + 2
Nst[, 2:5] <- NA # cols 2:5 of N are deterministic, N <- S + R.
R1 <- apply(data$y, c(1,2), max) # Observed max. counts + 1 as inits
R1[,1] <- NA
inits <- function(){list( lambda = runif(1, 6, 16), phi = runif(1),
    gamma = runif(1), p = runif(1), N = Nst, R = R1 + 1 )}
# Parameters monitored
params <- c("lambda", "phi", "gamma", "p")
# MCMC settings
na <- 1000 ; ni <- 25000 ; nt <- 4 ; nb <- 5000 ; nc <- 3
# Call JAGS (ART 2 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "DM1.txt", n.adapt = na, n.chains = nc, n.thin = nt,
  n.iter = ni, n.burnin = nb, parallel = TRUE)

par(mfrow = c(2, 3)) ; traceplot(out1) ; par(mfrow = c(1, 1))
print(out1, 3)
# Per-capita recruitment parameterisation
# mean sd 2.5% 50% 97.5% overlap0 f Rhat n.eff
# lambda 4.168 0.315 3.575 4.159 4.813 FALSE 1 1.000 15000
# phi 0.863 0.034 0.789 0.865 0.923 FALSE 1 1.002 1982
# gamma 0.232 0.036 0.168 0.230 0.311 FALSE 1 1.002 1682
# p 0.696 0.018 0.659 0.696 0.729 FALSE 1 1.001 5292

# To choose the absolute recruitment parameterization, we edit the BUGS code above to fit the
# absolute (or “constant”) recruitment parameterization simply by commenting out the recruitment line
# for “per capita” and then uncommenting the line previous to that.
# ~~~~~~~~~~~~~~~~~~ here's the new JAGS code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat(file = "DM1b.txt","
model {
  # Priors
  lambda ~ dunif(0, 100) # Population growth rate
  phi ~ dunif(0, 1) # apparent survival (omega in paper/unmarked)
  gamma ~ dunif(0, 5) # per-capita recruitment rate
  p ~ dunif(0, 1) # Detection probability
  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t]) # Survival process
      R[i,t+1] ~ dpois(gamma) # 'absolute’ recruitment = 'constant’
      # R[i,t+1] ~ dpois(N[i,t] * gamma) # per-capita recruitment = 'autoreg’
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    } # end t
    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,t,j] ~ dbin(p, N[i,t])
      } # end j
    } # end t
  } # end i
}
")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Then you can execute the line of code
# that calls JAGS, which should produce results similar to the following:
# Call JAGS (ART 1.3 min), check convergence and summarize posteriors
out1b <- jags(bdata, inits, params, "DM1b.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2, 3)) ; traceplot(out1b) ; par(mfrow = c(1, 1))
print(out1b, 2)
# Absolute parameterisation
# mean sd 2.5% 50% 97.5% overlap0 f Rhat n.eff
# lambda 4.09 0.31 3.51 4.08 4.72 FALSE 1 1 8586
# phi 0.86 0.03 0.79 0.86 0.92 FALSE 1 1 748
# gamma 1.14 0.17 0.85 1.13 1.51 FALSE 1 1 1061
# p 0.70 0.02 0.66 0.70 0.73 FALSE 1 1 5387

# 2.5.3 The unmarked function pcountOpen
# --------------------------------------

# Prepare data
summary(umf <- unmarkedFramePCO(y = data$yy, numPrimary = data$nyears))
# Fit model, backtransform and compare with truth
# Dynamics [ constant (as in data simulation)
(fm1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K = max(data$yy) + 100,
    dynamics = "constant", control = list(trace=TRUE, REPORT=1)) )
# Abundance:
# Estimate SE z P(>|z|)
# 1.4 0.0755 18.5 1.66e-76
# Recruitment:
# Estimate SE z P(>|z|)
# 0.104 0.144 0.721 0.471
# Apparent Survival:
# Estimate SE z P(>|z|)
# 1.87 0.276 6.75 1.44e-11
# Detection:
# Estimate SE z P(>|z|)
# 0.87 0.0825 10.5 5.93e-26
# AIC: 2466.449
# Back-transformation of parameters (full output not printed)
(lam <- coef(backTransform(fm1, "lambda"))) # or
(om <- plogis(coef(fm1, type="omega"))) # Apparent survival !
(gam <- exp(coef(fm1, type="gamma")))
(p <- plogis(coef(fm1, type="det")))
# lam(Int)
# 4.046917
# omega(Int)
# 0.8660311
# gamConst(Int)
# 1.109454
# p(Int)
# 0.7046613
# Dynamics = "autoreg"
# (not consistent with data generating model)
(fm2 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p = ~1, data = umf,
    dynamics = "autoreg", K = max(data$yy) + 100,
    control = list(trace=TRUE, REPORT = 1)))
# Abundance (log-scale):
# Estimate SE z P(>|z|)
# 1.42 0.0757 18.7 3.3e-78
# Recruitment (log-scale):
# Estimate SE z P(>|z|)
# -1.48 0.157 -9.4 5.4e-21
# Apparent Survival (logit-scale):
# Estimate SE z P(>|z|)
# 1.87 0.297 6.29 3.1e-10
# Detection (logit-scale):
# Estimate SE z P(>|z|)
# 0.844 0.0842 10 1.19e-23
# AIC: 2490.501
# Backtransformation of parameters (output omitted)
(lam <- exp(coef(fm2, type = "lambda")))
(om <- plogis(coef(fm2, type = "omega")))
(gam <- exp(coef(fm2, type = "gamma")))
(p <- plogis(coef(fm2, type = "det")))


# 2.5.4 Nonrobust design data
# (no code)

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
# Check with higher value of K: find that K[100 is OK. ART ~ 1.5 min
system.time(summary(fm3x <- pcountOpen(lam = ~1, gam = ~1, omega = ~1,
    p = ~1, data = umf, dynamics = "autoreg", K = 150,
    control = list(trace = TRUE, REPORT = 1)))) # not shown # 2 mins
# Backtransform parameter estimates and compare with truth
(lam <- exp(coef(fm3, type = "lambda"))) # Initial abundance
(gam <- exp(coef(fm3, type = "gamma"))) # Recruitment
(om <- plogis(coef(fm3, type = "omega"))) # Apparent survival
(p <- plogis(coef(fm3, type = "det"))) # Detection
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
# Add in final covariate, the site covariate on phi/omega. ART ~ 2.6 hrs
tmp <- coef(fm4.4)
inits <- c(tmp[1:5], 0, tmp[6:7])
system.time(summary(fm4.5 <- pcountOpen(lam = ~ cov.lam, gam = ~cov.gamma,
    omega = ~ cov.phi, p = ~ cov.p, data = umf,
    dynamics = "autoreg", se = FALSE, control = list(trace = TRUE, REPORT = 1),
    starts = inits) ))   # 9.6 hrs / 14 hrs
# Put the fitted models in a "fitList"
fms <- fitList(
"lam(.)gamma(.)phi(0)p(.)" = fm4.1,
"lam(cov)gamma(.)phi(0)p(.)" = fm4.2,
"lam(cov)gamma(.)phi(0)p(cov)" = fm4.3,
"lam(cov)gamma(cov)phi(0)p(cov)" = fm4.4,
"lam(cov)gamma(cov)phi(cov)p(cov)" = fm4.5)
# Rank them by AIC (singular Hessian message is due to SE=F option)
(ms <- modSel(fms))
# nPars AIC delta AICwt cumltvWt
# lam(cov)gamma(cov)phi(cov)p(cov) 8 2666.10 0.00 1.0e+00 1.00
# lam(cov)gamma(cov)phi(0)p(cov) 7 2688.62 22.52 1.3e-05 1.00
# lam(cov)gamma(.)phi(0)p(cov) 6 2813.65 147.55 9.1e-33 1.00
# lam(cov)gamma(.)phi(0)p(.) 5 3207.85 541.75 2.3e-118 1.00
# lam(.)gamma(.)phi(0)p(.) 4 3231.83 565.73 1.4e-123 1.00
# Coefficients for each model and Table With Everything (not printed)
print(coef(ms), digits=3)
toExport <- as(ms, "data.frame")
str(toExport) # output not printed

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
# Add in final covariate, the site covariate on phi/omega
tmp <- coef(fm4.4.se)
inits <- c(tmp[1:5], 0, tmp[6:7])
system.time(summary(fm4.5.se <- pcountOpen(lam = ~ cov.lam,
    gam = ~ cov.gamma, omega = ~ cov.phi, p= ~ cov.p, data = umf,
    dynamics = "autoreg", control = list(trace = TRUE, REPORT = 1), starts = inits) )) # 7.7 / 10.5hrs

# To get a sense of the time comparison with and without SE computation, we ran model fm4.5 both
# ways (here note the use of estimates from model fm4.5 as starting values in both cases):
# ~~~ These take a week, do not include in routine tests ~~~~~~~~~~~~~~~~~~~~~~
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
      # R[i,t+1] ~ dpois(gamma[i]) # 'absolute’ recruitment = 'constant’
      tmp[i,t] <- N[i,t] * gamma[i] # per-capita recruitment = 'autoreg’
      R[i,t+1] ~ dpois(tmp[i,t]) # per-capita recruitment = 'autoreg’
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    } # end t
    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        logit(p[i,t,j]) <- alpha.p + beta.p*cov.p[i,t,j]
        C[i,t,j] ~ dbin(p[i,t,j], N[i,t])
      } # end j
    } # end t
  } # end i
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
# na <- 1000 ; ni <- 500000 ; nt <- 250 ; nb <- 250000 ; nc <- 3
na <- 1000 ; ni <- 50000 ; nt <- 25 ; nb <- 25000 ; nc <- 3 # ~~~ for testing, 4 mins
# Call JAGS (ART 66 min), check convergence and summarize posteriors
out4 <- jags(bdata, inits, params, "DM2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2, 3)) ; traceplot(out4) ; par(mfrow = c(1, 1))
print(out4, 2) # not printed
# Compare estimates with truth
cbind('truth' = unlist(data[4:11]), round(out4$summary[1:8, c(1,3,7)],3))
# truth mean 2.5% 97.5%
# mean.lambda 4.0 4.247 3.645 4.857
# mean.gamma.rel 0.5 0.706 0.513 0.837
# mean.phi 0.8 0.590 0.451 0.793
# mean.p 0.7 0.695 0.659 0.727
# beta.lam 1.0 0.632 0.364 0.891
# beta.gamma 1.0 -0.813 -1.357 -0.337
# beta.phi -1.0 0.873 0.648 1.085
# beta.p -1.0 -1.005 -1.140 -0.870

# 2.5.6 Case study of the Swiss MHB data for the green woodpecker:
# Bayesian analysis of the Dail-Madsen model for robust design data with BUGS
# ---------------------------------------------------------------------------
source(file="AHM2-02.02.R")

# Bundle data
str(bdata <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2], nyears = dim(C)[3],
    elev = elev, forest = forest, DATE = DATE, INT = INT))
# List of 8
# $ C : int [1:267, 1:3, 1:14] 0 3 0 0 0 0 0 0 0 0 ...
# $ nsites : int 267
# $ nsurveys : int 3
# $ nyears : int 14
# $ elev : num [1:267, 1] -1.206 -1.16 -0.184 -0.439 -0.126 ...
# $ forest : num [1:267, 1] -1.1529 -0.467 0.0023 -0.9002 -0.106 ...
# $ DATE : num [1:267, 1:3, 1:14] -1.09 -1.32 -1.23 -1.27 -1.36 ...
# $ INT : num [1:267, 1:3, 1:14] -0.531 -0.957 0.168 -0.451 -0.864 ...
# MCMC settings
na <- 1000 ; ni <- 150000 ; nt <- 10 ; nb <- 50000 ; nc <- 3
# Call JAGS (ART 88 min), check convergence and summarize posteriors
# out5 <- jags(bdata, inits, params, "DM3.txt", n.adapt = na,
    # n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3)) ; traceplot(out5) ; par(mfrow = c(1, 1))
# print(out5, 3)
# mean sd 2.5% 50% 97.5% overlap0 f Rhat n.eff
# lambda 0.599 0.068 0.472 0.597 0.740 FALSE 1 1.001 3482
# phi 0.821 0.012 0.797 0.821 0.845 FALSE 1 1.005 425
# gamma 0.378 0.017 0.345 0.378 0.412 FALSE 1 1.001 1582
# p 0.290 0.009 0.273 0.290 0.308 FALSE 1 1.004 575

# Bundle data
str(bdata <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
  nyears = dim(C)[3], elev = as.vector(elev),forest = as.vector(forest),
  DATE = DATE, length = peckers$route.length, INT = INT)) # note length added
# Specify model in BUGS language
cat(file = "DM4.txt","
model {
  # Priors
  alpha.lam ~ dnorm(0,0.1) # Abundance parameters
  beta.elev ~ dnorm(0,0.1)
  beta.elev2 ~ dnorm(0,0.1)
  beta.for ~ dnorm(0,0.1)
  beta.ilen ~ dunif(-10,0) # Enforce negative, see section 7.9 of AHM1
  phi ~ dunif(0, 1) # App. survival (omega in paper/unmarked)
  beta.gamma ~ dnorm(0,0.1) # Recruitment model parameter
  gamma0 <- exp(beta.gamma)
  alpha.p ~ dnorm(0,0.1) # Detection probability parameters
  beta.jul ~ dnorm(0,0.1)
  beta.jul2 ~ dnorm(0,0.1)
  beta.int ~ dnorm(0,0.1)
  beta.int2 ~ dnorm(0,0.1)
  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    log(lambda[i]) <- alpha.lam + beta.elev*elev[i] +
    beta.elev2*elev[i]*elev[i] + beta.for*forest[i] +
    beta.ilen*(1/length[i])
    N[i,1] ~ dpois(lambda[i])
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t])
      R[i,t+1] ~ dpois(gamma0) # constant recruitment
      ###R[i,t+1] ~ dpois(N[i,t] * gamma0) # per-capita recruitment
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    } # end t
    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        logit(p[i,j,t]) <- alpha.p + beta.jul*DATE[i,j,t] +
            beta.jul2*DATE[i,j,t]*DATE[i,j,t] + beta.int*INT[i,j,t] +
            beta.int2*INT[i,j,t]*INT[i,j,t]
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t])
      } # end j
    } # end t
  } # end i
  # Derived quantities
  Nbar <- mean(N)
  lambar <- mean(lambda)
} # end model
")
# Initial values, sometimes very challenging
Rst <- apply(C, c(1,3), max, na.rm = TRUE)
Rst[Rst == '-Inf'] <- 1
Rst[,1] <- NA
Nst <- array(NA, dim = dim(Rst))
# tmp <- apply(C[,1], 1, max, na.rm = TRUE)
tmp <- apply(apply(C, c(1,3), max, na.rm = TRUE),1,max, na.rm = TRUE)
tmp[tmp == '-Inf'] <- 2
Nst[,1] <- tmp
# Initial values
inits <- function(){list(R = Rst, N = Nst+1, alpha.lam = rnorm(1, -1,1),
    beta.elev = rnorm(1), beta.elev2 = rnorm(1), beta.for = rnorm(1), beta.ilen = runif(1,-1,0),
    phi = 0.5, alpha.gamma = rnorm(1), beta.gamma = rnorm(1), sigma.gamma = 0.5, alpha.p = rnorm(1),
    beta.jul = rnorm(1), beta.jul2 = rnorm(1), beta.int = rnorm(1), beta.int2 = rnorm(1))}
# Parameters monitored
params <-c("alpha.lam", "beta.elev", "beta.elev2", "beta.for", "beta.ilen", "phi",
    "alpha.gamma", "beta.gamma", "sigma.gamma", "gamma0", "alpha.p", "beta.jul",
    "beta.jul2", "beta.int", "beta.int2")
# MCMC settings
# na <- 1000 ; ni <- 70000 ; nt <- 4 ; nb <- 10000 ; nc <- 4
na <- 1000 ; ni <- 7000 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ for testing, 7 mins
# try different seeds in hopes of getting one that runs
set.seed(299, kind = "Mersenne-Twister")
# Call JAGS (ART 83 min), check convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "DM4.txt", n.adapt = na, n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2,3)) ; traceplot(out6) ; par(mfrow = c(1, 1))
print(out6, digits=2)
# mean sd 2.5% 50% 97.5% overlap0 f Rhat n.eff
# alpha.lam 0.92 0.44 0.04 0.92 1.75 FALSE 0.98 1.00 1120
# beta.elev -0.43 0.14 -0.70 -0.43 -0.16 FALSE 1.00 1.00 867
# beta.elev2 -0.40 0.18 -0.77 -0.39 -0.06 FALSE 0.99 1.00 2197
# beta.for 0.18 0.10 -0.01 0.18 0.37 TRUE 0.97 1.00 4082
# beta.ilen -5.15 2.00 -9.05 -5.14 -1.28 FALSE 1.00 1.00 1558
# phi 0.83 0.01 0.80 0.83 0.86 FALSE 1.00 1.01 273
# beta.gamma -0.99 0.05 -1.09 -0.99 -0.89 FALSE 1.00 1.00 715
# gamma0 0.37 0.02 0.34 0.37 0.41 FALSE 1.00 1.00 706
# alpha.p -0.99 0.05 -1.09 -0.99 -0.89 FALSE 1.00 1.01 243
# beta.jul -0.22 0.02 -0.27 -0.22 -0.18 FALSE 1.00 1.00 5005
# beta.jul2 0.04 0.02 0.00 0.04 0.08 TRUE 0.97 1.00 10157
# beta.int 0.24 0.04 0.16 0.24 0.31 FALSE 1.00 1.01 362
# beta.int2 -0.09 0.02 -0.13 -0.09 -0.05 FALSE 1.00 1.00 4514

# 2.5.7 Likelihood analysis of the Swiss woodpecker data in unmarked
# ------------------------------------------------------------------

# Data management
nyears <- dim(C)[3]
Cwide <- C[,,1]
date <- DATE[,,1]
dur <- DUR[,,1]
int <- INT[,,1]
for(i in 2:nyears){
  Cwide <- cbind(Cwide, C[,,i])
  date <- cbind(date, DATE[,,i])
  dur <- cbind(dur, DUR[,,i])
  int <- cbind(int, INT[,,i])
}
library(unmarked)
sitecovs <- data.frame(cbind(elev = as.numeric(elev),
    forest = as.numeric(forest), ilength = 1 / peckers$route.length))
obscovs <- list(date = date, dur = dur, int = int)
summary(umf <- unmarkedFramePCO(y = Cwide, numPrimary = nyears,
    siteCovs = sitecovs, obsCovs = obscovs))

# We set Kmax here and then vary it later to check stability of the MLEs
Kmax <- 150
# Fit models with autoregression formulation
# Go-drink-some-coffee warning: ART about 10 min for each
system.time(dm0 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p = ~1,
    data = umf, dynamics = "autoreg", K = Kmax ) )
dm1 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p = ~date, data = umf,
    dynamics = "autoreg", K = Kmax )
dm2 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p = ~ date + I(date^2),
    data = umf, dynamics = "autoreg", K = Kmax )
dm3 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p = ~date + I(date^2) +
    int, data = umf, dynamics = "autoreg", K = Kmax )
system.time(dm4 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p = ~date +
    I(date^2) + int + I(int^2), data = umf, dynamics = "autoreg", K = Kmax ))  # 80 secs
# Models with constant dynamics: ART 1 – 10 min
system.time(dm0b <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p = ~ 1,
    data = umf, dynamics = "constant", K = Kmax ) )
dm1b <- pcountOpen(lam = ~elev, gam = ~1, omega = ~1, p = ~date , data = umf,
    dynamics = "constant", K = Kmax)
dm2b <- pcountOpen(lam = ~elev + I(elev^2), gam = ~1, omega = ~1, p = ~date + I(date^2),
    data = umf, dynamics = "constant", K = Kmax)
dm3b <- pcountOpen(lam = ~elev + I(elev^2) + forest, gam = ~1, omega = ~1,
    p = ~date + I(date^2) + int, data = umf, dynamics = "constant", K = Kmax)
system.time(dm4b <- pcountOpen(lam = ~elev + I(elev^2) + forest + ilength, gam = ~1,
    omega = ~1, p = ~date + I(date^2) + int + I(int^2), data = umf,
    dynamics = "constant", K = Kmax) )  # 9 mins
# Construct fitList and produce AIC table
fl <- fitList(dm0 = dm0, dm1 = dm1, dm2 = dm2, dm3 = dm3, dm4 = dm4, dm0b = dm0b, dm1b = dm1b,
    dm2b = dm2b, dm3b = dm3b, dm4b = dm4b)
modSel(fl)
# nPars AIC delta AICwt cumltvWt
# dm4 8 15745.49 0.00 1.0e+00 1.00
# dm3 7 15757.21 11.72 2.8e-03 1.00
# dm2 6 15765.96 20.47 3.6e-05 1.00
# dm1 5 15767.22 21.73 1.9e-05 1.00
# ... truncated ...
# Now build-up some abundance models. ART 5 – 10 min
system.time(dm5 <- pcountOpen(lam = ~elev , gam = ~1, omega = ~1, p = ~date + I(date^2) + int +
    I(int^2), data = umf, dynamics = "constant", K = Kmax ))  # 5 mins
dm6 <- pcountOpen(lam = ~elev + I(elev^2) , gam = ~1, omega = ~1, p = ~date + I(date^2) + int +
    I(int^2), data = umf, dynamics = "constant", K = Kmax )
dm7 <- pcountOpen(lam = ~elev + I(elev^2) + forest , gam = ~1, omega = ~1, p = ~date +
    I(date^2) + int + I(int^2), data = umf, dynamics = "constant", K = Kmax)
system.time(dm8 <- pcountOpen(lam = ~elev + I(elev^2) + forest + ilength, gam = ~1,
    omega = ~1, p = ~date + I(date^2) + int + I(int^2), data = umf, dynamics = "constant", K = Kmax)) # 9 mins
# Organize the results into a fit list and produce AIC table
fl <- fitList(dm5, dm6, dm7, dm8)
modSel(fl)
# nPars AIC delta AICwt cumltvWt
# dm8 12 16101.60 0.00 0.80390 0.80
# dm7 11 16105.55 3.95 0.11159 0.92
# dm6 10 16106.12 4.52 0.08369 1.00
# dm5 9 16115.39 13.79 0.00081 1.00

dm8
# Abundance:
# Estimate SE z P(>|z|)
# (Intercept) 0.894 0.4781 1.87 0.06144
# elev -0.413 0.1368 -3.02 0.00252
# I(elev^2) -0.375 0.1786 -2.10 0.03589
# forest 0.183 0.0977 1.88 0.06050
# ilength -5.012 2.1848 -2.29 0.02178
# Recruitment:
# Estimate SE z P(>|z|)
# -0.978 0.0515 -19 1.27e-80
# Apparent Survival:
# Estimate SE z P(>|z|)
# 1.57 0.0965 16.3 2.22e-59
# Detection:
# Estimate SE z P(>|z|)
# (Intercept) -0.9814 0.0532 -18.45 4.89e-76
# date -0.2274 0.0247 -9.20 3.64e-20
# I(date^2) 0.0367 0.0213 1.73 8.41e-02
# int 0.2430 0.0382 6.36 2.01e-10
# I(int^2) -0.0879 0.0200 -4.39 1.13e-05
