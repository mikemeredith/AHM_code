#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 2 hrs

library(unmarked)
library(jagsUI)
library(AHMbook)

# ~~~ need to prepare data ~~~~
source(file="AHM2_02.02.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.5 Dynamic N-mixture model of Dail-Madsen
# ==========================================

# 2.5.6 Case study of the Swiss MHB data for the green woodpecker:
#   Bayesian analysis of the Dail-Madsen model for robust design data with BUGS
# -----------------------------------------------------------------------------

# Bundle data
str(bdata <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3], elev = elev, forest = forest, DATE = DATE, INT = INT))
# List of 8
# $ C       : int [1:267, 1:3, 1:14] 0 3 0 0 0 0 0 0 0 0 ...
# $ nsites  : int 267
# $ nsurveys: int 3
# $ nyears  : int 14
# $ elev    : num [1:267, 1] -1.206 -1.16 -0.184 -0.439 -0.126 ...
# $ forest  : num [1:267, 1] -1.1529 -0.467 0.0023 -0.9002 -0.106 ...
# $ DATE    : num [1:267, 1:3, 1:14] -1.09 -1.32 -1.23 -1.27 -1.36 ...
# $ INT     : num [1:267, 1:3, 1:14] -0.531 -0.957 0.168 -0.451 -0.864 ...

# Specify model in BUGS language
cat(file = "DM3.txt","
model {
  # Priors
  lambda ~ dunif(0, 100)    # Population growth rate
  phi ~ dunif(0, 1)         # apparent survival (or omega)
  gamma ~ dunif(0, 5)       # per-capita recruitment rate
  p ~ dunif(0, 1)           # Detection probability

  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t])
      R[i,t+1] ~ dpois(gamma)              # 'absolute' recruitment = 'constant'
      ###R[i,t+1] ~ dpois(N[i,t] * gamma)  # per-capita recr. = 'autoreg'
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    }

    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,j,t] ~ dbin(p, N[i,t])
      }
    }
  }
}
")

# Initial values
Rst <- apply(C, c(1,3), max, na.rm = TRUE)
Rst[Rst == '-Inf'] <- 1
Rst[,1] <- NA
N1 <- apply(apply(C, c(1,3), max, na.rm = TRUE), 1, max, na.rm = TRUE)
N1[N1 == '-Inf'] <- 2
Nst <- array(NA, dim = dim(Rst))
Nst[,1] <- N1
inits <- function(){list(lambda = runif(1, 1, 8), phi = runif(1),
    gamma = runif(1), p = runif(1), R = Rst+1, N = Nst+2)}

# Parameters monitored
params <- c("lambda", "phi", "gamma", "p")

# MCMC settings
# na <- 1000  ; ni <- 150000  ;  nt <- 10  ;  nb <- 50000  ;  nc <- 3
na <- 1000  ; ni <- 15000  ;  nt <- 1  ;  nb <- 5000  ;  nc <- 3  # ~~~~ for testing

# Call JAGS (ART 88 min), check convergence and summarize posteriors
out5 <- jags(bdata, inits, params, "DM3.txt", n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# par(mfrow = c(2, 3))  #  ~~~ replace with 'layout' argument
traceplot(out5, layout=c(2,3))
print(out5, 3)
              # mean      sd      2.5%       50%     97.5% overlap0 f  Rhat n.eff
# lambda       0.599   0.068     0.472     0.597     0.740    FALSE 1 1.001  3482
# phi          0.821   0.012     0.797     0.821     0.845    FALSE 1 1.005   425
# gamma        0.378   0.017     0.345     0.378     0.412    FALSE 1 1.001  1582
# p            0.290   0.009     0.273     0.290     0.308    FALSE 1 1.004   575


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
      R[i,t+1] ~ dpois(gamma0)             # constant recruitment
      ###R[i,t+1] ~ dpois(N[i,t] * gamma0) # per-capita recruitment
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    }

    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        logit(p[i,j,t]) <- alpha.p + beta.jul*DATE[i,j,t] +
            beta.jul2*DATE[i,j,t]*DATE[i,j,t] + beta.int*INT[i,j,t] +
            beta.int2*INT[i,j,t]*INT[i,j,t]
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t])
      }
    }
  }
  # Derived quantities
  Nbar <- mean(N)
  lambar <- mean(lambda)
}
")

# Initial values, sometimes very challenging
Rst <- apply(C, c(1,3), max, na.rm = TRUE)
Rst[Rst == '-Inf'] <- 1
Rst[,1] <- NA
Nst <- array(NA, dim = dim(Rst))
tmp <- apply(C, 1, max, na.rm = TRUE)
tmp[tmp == '-Inf'] <- 2
Nst[,1] <- tmp

# Initial values
inits <- function(){list(R = Rst, N = Nst+1, alpha.lam = rnorm(1, -1,1),
    beta.elev = rnorm(1), beta.elev2 = rnorm(1), beta.for = rnorm(1), beta.ilen = runif(1,-1,0),
    phi = 0.5, alpha.gamma = rnorm(1), beta.gamma = rnorm(1), sigma.gamma = 0.5, alpha.p = rnorm(1),
    beta.jul = rnorm(1), beta.jul2 = rnorm(1), beta.int = rnorm(1), beta.int2 = rnorm(1))}

# Parameters monitored
params <-c("alpha.lam", "beta.elev", "beta.elev2", "beta.for",
    "beta.ilen", "phi", "alpha.gamma", "beta.gamma", "sigma.gamma", "gamma0",
    "alpha.p", "beta.jul", "beta.jul2", "beta.int", "beta.int2")

# MCMC settings
# na <- 1000 ; ni <- 70000 ; nt <- 4 ; nb <- 10000 ; nc <- 4
na <- 1000 ; ni <- 7000 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ for testing, 7 mins

# try different seeds in hopes of getting one that runs
set.seed(299, kind = "Mersenne-Twister")
# Call JAGS (ART 83 min), check convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "DM4.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,3))  #  ~~~ replace with 'layout' argument
traceplot(out6, layout=c(2,3))
print(out6, digits=2)
#             mean   sd  2.5%   50% 97.5% overlap0    f Rhat n.eff
# alpha.lam   0.92 0.44  0.04  0.92  1.75    FALSE 0.98 1.00  1120
# beta.elev  -0.43 0.14 -0.70 -0.43 -0.16    FALSE 1.00 1.00   867
# beta.elev2 -0.40 0.18 -0.77 -0.39 -0.06    FALSE 0.99 1.00  2197
# beta.for    0.18 0.10 -0.01  0.18  0.37     TRUE 0.97 1.00  4082
# beta.ilen  -5.15 2.00 -9.05 -5.14 -1.28    FALSE 1.00 1.00  1558
# phi         0.83 0.01  0.80  0.83  0.86    FALSE 1.00 1.01   273
# beta.gamma -0.99 0.05 -1.09 -0.99 -0.89    FALSE 1.00 1.00   715
# gamma0      0.37 0.02  0.34  0.37  0.41    FALSE 1.00 1.00   706
# alpha.p    -0.99 0.05 -1.09 -0.99 -0.89    FALSE 1.00 1.01   243
# beta.jul   -0.22 0.02 -0.27 -0.22 -0.18    FALSE 1.00 1.00  5005
# beta.jul2   0.04 0.02  0.00  0.04  0.08     TRUE 0.97 1.00 10157
# beta.int    0.24 0.04  0.16  0.24  0.31    FALSE 1.00 1.01   362
# beta.int2  -0.09 0.02 -0.13 -0.09 -0.05    FALSE 1.00 1.00  4514

# ~~~ save the work so far ~~~~
save.image("AHM2_02.05.6.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ extra code for figure 2.10 ~~~~~~~~~~~~~~~~~
x.date <- seq(min(DATE,na.rm=TRUE), max(DATE,na.rm=TRUE),,100)

elev.grid <- seq(min(peckers$elev),max(peckers$elev),,100)
x.elev <- standardize2match(elev.grid, peckers$elev)

len <- peckers$route
len.grid <- seq(min(len),max(len),,100)

forest.grid <- seq(min(peckers$forest),max(peckers$forest),,100)
x.fore <- standardize2match(forest.grid, peckers$forest)

op <- par(mfrow=c(1, 3), mar=c(5,5,3,3),cex.lab = 1.6, cex.main = 1.2)
parms <- out6$summary[,"mean"]
lam.fit <- exp(cbind(rep(1,100), x.elev,x.elev^2, 0)%*%parms[c("alpha.lam","beta.elev","beta.elev2","beta.for")])
plot(elev.grid, lam.fit, xlab="Elevation (m)", ylab="E(N)", type="l",
    ylim=c(0,3), main="(A)", frame=FALSE, lwd=3)
lam.fit <- cbind(rep(1,100), 0, 0, x.fore)%*%parms[c("alpha.lam","beta.elev","beta.elev2","beta.for")]
plot(forest.grid, lam.fit, xlab="Forest cover (%)", ylab="E(N)", type="l",
    ylim=c(0,3), main="(B)", frame=FALSE,lwd=3)
plot(len.grid, exp(parms["beta.ilen"]*(1/len.grid)), type="l", xlab="Route length (km)",
    ylab="Quadrat saturation", lwd=3, main="(C)", frame=FALSE)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
fl <- fitList(dm0 = dm0, dm1 = dm1, dm2 = dm2, dm3 = dm3, dm4 = dm4, dm0b = dm0b,
    dm1b = dm1b, dm2b = dm2b, dm3b = dm3b, dm4b = dm4b)
modSel(fl)
#     nPars      AIC delta   AICwt cumltvWt
# dm4     8 15745.49  0.00 1.0e+00        1
# dm3     7 15757.21 11.72 2.8e-03        1
# dm2     6 15765.96 20.47 3.6e-05        1
# dm1     5 15767.22 21.73 1.9e-05        1
# ... truncated ...

# Now build-up some abundance models. ART 5 – 10 min
system.time(dm5 <- pcountOpen(lam = ~elev , gam = ~1, omega = ~1,
    p = ~date + I(date^2) + int + I(int^2), data = umf, dynamics = "constant",
    K = Kmax ))  # 5 mins
dm6 <- pcountOpen(lam = ~elev + I(elev^2) , gam = ~1, omega = ~1,
    p = ~date + I(date^2) + int + I(int^2), data = umf, dynamics = "constant",
    K = Kmax )
dm7 <- pcountOpen(lam = ~elev + I(elev^2) + forest , gam = ~1, omega = ~1,
    p = ~date + I(date^2) + int + I(int^2), data = umf, dynamics = "constant",
    K = Kmax)
system.time(dm8 <- pcountOpen(lam = ~elev + I(elev^2) + forest + ilength,
    gam = ~1, omega = ~1, p = ~date + I(date^2) + int + I(int^2), data = umf,
    dynamics = "constant", K = Kmax)) # 9 mins

# Organize the results into a fit list and produce AIC table
fl <- fitList(dm5, dm6, dm7, dm8)
modSel(fl)
#     nPars      AIC delta   AICwt cumltvWt
# dm8    12 16101.60  0.00 0.80390     0.80
# dm7    11 16105.55  3.95 0.11159     0.92
# dm6    10 16106.12  4.52 0.08369     1.00
# dm5     9 16115.39 13.79 0.00081     1.00

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

# ~~~ save the work so far ~~~~~
save.image("AHM2_02.05.7.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
