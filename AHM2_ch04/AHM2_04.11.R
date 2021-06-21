#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# Run time with the full number of iterations: 10 hrs

library(AHMbook)
library(jagsUI)

# 4.11 Accounting for preferential sampling in a bird population study
# ====================================================================

# Load the data set from AHMbook
data(FrenchPeregrines)
str(dat <- FrenchPeregrines)

# Extract data for modeling
ain <- which(dat$department == 'Ain')             # Sites in Dep. Ain
jura <- which(dat$department == 'Jura')           # Sites in Dep. Jura
doubs <- which(dat$department == 'Doubs')         # Sites in Dep. Doubs
ht <- as.numeric(dat$height)                      # Cliff height
y <- as.matrix(dat[,4:56])                        # Detection/Nondetection data
nsites <- nrow(y)
nyears <- ncol(y)
year <- 1964:2016

# Produce some summaries, including ratio estimator
n.occ.obs <- apply(y, 2, sum, na.rm = TRUE)          # Observed N pairs
n.visited <- apply(y, 2, function(x) sum(!is.na(x))) # Number of visited sites
n.ratio <- n.occ.obs / (n.visited / 284)

# ~~~~~~~~ code for figure 4.30 ~~~~~~~~~~
plot(year, n.visited, xlab = 'Year', ylab = 'Number', type = 'o',
    pch = 0, ylim = c(0, 284), frame = FALSE)
abline(h = nsites, lty = 2)
points(year, n.occ.obs, pch = 16, type = 'o')
points(year, n.ratio, pch = 1, type = 'o', col = 'blue')
legend('bottomright', c('Total number of available sites',
    'Number of sites visited',
    'Ratio estimator of population size', 'Number of pairs observed'),
    lty = c(2, 1,1,1), pch = c(NA, 0, 1, 16),
    col = c('black', 'black', 'blue', 'black'), bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute visitation data V
V <- y ; V[V == 0] <- 1 ; V[is.na(V)] <- 0
y[1:5, 1:5] # Look at parts of y and V
#      yr1964 yr1965 yr1966 yr1967 yr1968
# [1,]      1      1      1      1      1
# [2,]     NA     NA     NA     NA     NA
# [3,]     NA     NA     NA     NA     NA
# [4,]     NA     NA      0      0      0
# [5,]     NA     NA     NA     NA     NA

V[1:5, 1:5]
#      yr1964 yr1965 yr1966 yr1967 yr1968
# [1,]      1      1      1      1      1
# [2,]      0      0      0      0      0
# [3,]      0      0      0      0      0
# [4,]      0      0      1      1      1
# [5,]      0      0      0      0      0

# Compute prop. years with pairs (given surveyed) and of years with
# surveys and plot, both starting from first ever visit of a site
first.visit <- apply(!is.na(y), 1, function(x) which(x)[1])
prop.visit <- prop.pair <- numeric(length = nsites)
for(i in 1:nsites){
  prop.visit[i] <- mean(V[i, first.visit[i]:nyears],na.rm = TRUE)
  prop.pair[i] <- mean(y[i, first.visit[i]:nyears],na.rm = TRUE)
}

# Fig. 4.31
op <- par(mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
plot(prop.pair, prop.visit, xlab = 'Proportion of years occupied',
    ylab = 'Proportion of years visited', pch = 16, cex = 2, ylim = c(0,1),
    frame = FALSE, col = rgb(0,0,0,0.4))
par(op)

# Data bundle
str(bdata <- list(y = y, height = ht, nsites = nsites, nyears = nyears,
    ain = ain, jura = jura, doubs = doubs, year = ((1964:2016)-1990) / 26))
# List of 8
# $ y     : int [1:284, 1:53] 1 NA NA NA NA NA NA NA 1 1 ...
# $ height: num [1:284] 2 2 1 3 1 3 2 3 3 3 ...
# $ nsites: int 284
# $ nyears: int 53
# $ ain   : int [1:93] 1 2 3 4 5 6 7 8 9 10 ...
# $ jura  : int [1:89] 94 95 96 97 98 99 100 101 102 103 ...
# $ doubs : int [1:102] 183 184 185 186 187 188 189 190 191 192 ...
# $ year  : num [1:53] -1 -0.962 -0.923 -0.885 -0.846 ...

# Model 1: no preferential sampling (PS)
# Specify model in BUGS language
cat(file = "dynocc1.txt", "
model {

  # Priors and models for parameters
  psi1 ~ dbeta(1, 1)                         # Initial occupancy

  # Model for phi and gamma: cliff height + site + smooth random year effects
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      logit(phi[i,t]) <- lphi[i,t]
      lphi[i,t] <- lphi.site[i] + lphi.year[t]
      logit(gamma[i,t]) <- lgamma[i,t]
      lgamma[i, t] <- lgamma.site[i] + lgamma.year[t]
    }
    lphi.site[i] ~ dnorm(alpha.lphi[height[i]], tau.lphi.site)
    lgamma.site[i] ~ dnorm(alpha.lgamma[height[i]], tau.lgamma.site)
  }

  # Priors for phi and gamma intercepts
  for(k in 1:3){
    alpha.lphi[k] <- logit(initial.phi[k])
    initial.phi[k] ~ dbeta(1, 1)
    alpha.lgamma[k] <- logit(initial.gamma[k])
    initial.gamma[k] ~ dbeta(1, 1)
  }
  tau.lphi.site <- pow(sd.lphi.site, -2)
  sd.lphi.site ~ dunif(0, 2)
  tau.lgamma.site <- pow(sd.lgamma.site, -2)
  sd.lgamma.site ~ dunif(0, 2)

  # Priors for year effects on phi and gamma with rw smoothers
  lphi.year[1] <- 0 # Set to zero to avoid overparameterization
  lgamma.year[1] <- 0
  for (t in 2:(nyears-1)){
    lphi.year[t] ~ dnorm(lphi.year[t-1], tau.eps.lphi)
    lgamma.year[t] ~ dnorm(lgamma.year[t-1], tau.eps.lgamma)
  }
  tau.eps.lphi <- pow(sd.eps.lphi,-2) # Hyperpriors for variances
  sd.eps.lphi ~ dunif(0, 1)
  tau.eps.lgamma <- pow(sd.eps.lgamma,-2)
  sd.eps.lgamma ~ dunif(0, 1)

  # Ecological and observation submodels confounded (no p)
  for (i in 1:nsites){
    y[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      y[i,t] ~ dbern(y[i,t-1]*phi[i,t-1] + (1-y[i,t-1])*gamma[i,t-1])
    }
  }

  # Derived parameters
  # Population occupancy and population size
  psi[1] <- psi1
  n.occ[1] <- sum(y[1:nsites,1])
  for (t in 2:nyears){
    n.occ[t] <- sum(y[1:nsites,t]) # Number of occupied sites
  }

  # Year-specific average values of phi and gamma
  for(t in 1:(nyears-1)){
    mean.phi.year[t] <- mean(phi[,t])
    mean.gamma.year[t] <- mean(gamma[,t])
  }

  # Average gamma and phi per cliff height category in central year (1989)
  for(k in 1:3){
    logit(gamma.cliff[k]) <- alpha.lgamma[k] + lgamma.year[26]
    logit(phi.cliff[k]) <- alpha.lphi[k] + lphi.year[26]
  }

  # Population size in each French Jura Departement
  for (t in 1:nyears){
    n.ain[t] <- sum(y[ain, t])
    n.jura[t] <- sum(y[jura, t])
    n.doubs[t] <- sum(y[doubs, t])
  }
}
")

# Initial values
inits <- function(){ list(psi1 = runif(1))}

# Parameters monitored
params <- c("psi1", "alpha.lphi", "initial.phi", "alpha.lgamma",
    "initial.gamma", "lphi.year", "lgamma.year", "mean.phi.year",
    "mean.gamma.year", "sd.eps.lphi", "sd.eps.lgamma", "gamma.cliff",
    "phi.cliff", "n.occ", "n.ain", "n.jura", "n.doubs", "y")

# MCMC settings
# na <- 1000 ; ni <- 20000 ; nt <- 10 ; nb <- 10000 ; nc <- 3  # 35 mins
na <- 1000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ testing, 6 mins

# Call JAGS (ART 43 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "dynocc1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
summary(out1) ; jags.View(out1) # not shown

# Model 2 : With PS
# '''''''''''''''''
# Model 2: modeling preferential sampling (PS) by joint modeling of
# occupancy and visitation
# Specify model in BUGS language
cat(file = "dynocc2.txt", "
model {

  # Submodel 1: model for detection/nondetection data (y)
  # ----------------------------------------------------

  # all this first part is identical to dynocc1 above and omitted in the book

  # ~~~ inserted code copied from above ~~~~~~~~~~~~~~~~
  # Priors and models for parameters
  psi1 ~ dbeta(1, 1) # Initial occupancy

  # Model for phi and gamma: cliff height + site + smooth random year effects
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      logit(phi[i,t]) <- lphi[i,t]
      lphi[i,t] <- lphi.site[i] + lphi.year[t]
      logit(gamma[i,t]) <- lgamma[i,t]
      lgamma[i, t] <- lgamma.site[i] + lgamma.year[t]
    }
    lphi.site[i] ~ dnorm(alpha.lphi[height[i]], tau.lphi.site)
    lgamma.site[i] ~ dnorm(alpha.lgamma[height[i]], tau.lgamma.site)
  }

  # Priors for phi and gamma intercepts
  for(k in 1:3){
    alpha.lphi[k] <- logit(initial.phi[k])
    initial.phi[k] ~ dbeta(1, 1)
    alpha.lgamma[k] <- logit(initial.gamma[k])
    initial.gamma[k] ~ dbeta(1, 1)
  }
  tau.lphi.site <- pow(sd.lphi.site, -2)
  sd.lphi.site ~ dunif(0, 2)
  tau.lgamma.site <- pow(sd.lgamma.site, -2)
  sd.lgamma.site ~ dunif(0, 2)

  # Priors for year effects on phi and gamma with rw smoothers
  lphi.year[1] <- 0 # Set to zero to avoid overparameterization
  lgamma.year[1] <- 0
  for (t in 2:(nyears-1)){
    lphi.year[t] ~ dnorm(lphi.year[t-1], tau.eps.lphi)
    lgamma.year[t] ~ dnorm(lgamma.year[t-1], tau.eps.lgamma)
  }
  tau.eps.lphi <- pow(sd.eps.lphi,-2) # Hyperpriors for variances
  sd.eps.lphi ~ dunif(0, 1)
  tau.eps.lgamma <- pow(sd.eps.lgamma,-2)
  sd.eps.lgamma ~ dunif(0, 1)

  # Ecological and observation submodels confounded (no p)
  for (i in 1:nsites){
    y[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      y[i,t] ~ dbern(y[i,t-1]*phi[i,t-1] + (1-y[i,t-1])*gamma[i,t-1])
    }
  }

  # Derived parameters
  # Population occupancy and population size
  psi[1] <- psi1
  n.occ[1] <- sum(y[1:nsites,1])
  for (t in 2:nyears){
    n.occ[t] <- sum(y[1:nsites,t]) # Number of occupied sites
  }

  # Year-specific average values of phi and gamma
  for(t in 1:(nyears-1)){
    mean.phi.year[t] <- mean(phi[,t])
    mean.gamma.year[t] <- mean(gamma[,t])
  }

  # Average gamma and phi per cliff height category in central year (1989)
  for(k in 1:3){
    logit(gamma.cliff[k]) <- alpha.lgamma[k] + lgamma.year[26]
    logit(phi.cliff[k]) <- alpha.lphi[k] + lphi.year[26]
  }

  # Population size in each French Jura Departement
  for (t in 1:nyears){
    n.ain[t] <- sum(y[ain, t])
    n.jura[t] <- sum(y[jura, t])
    n.doubs[t] <- sum(y[doubs, t])
  }
  # ~~~~~~~ end of inserted code ~~~~~~~~~~~~

  # Submodel 2: model for whether a site is visited or not (V)
  # ----------------------------------------------------------
  # Priors and linear models
  # Including an effect of past occupation history in visitation from t=2
  for (i in 1:nsites){
    for (t in 1:nyears){
      theta[i,t] <- ilogit(alpha.visit + beta.visit[1] * year[t] +
      beta.visit[2] * pow(year[t],2) + beta.visit[3] * pow(year[t],3) +
      kappa.lphi * lphi.site[i] + kappa.lgamma * lgamma.site[i])
    }
  }
  alpha.visit <- logit(theta.int)
  theta.int ~ dbeta(1, 1)
  for(v in 1:3){ # Coefficients of time
    beta.visit[v] ~ dnorm(0, 0.1)
  }
  kappa.lphi ~ dnorm(0, 0.5)               # first coefficient for PS
  kappa.lgamma ~ dnorm(0,0.5)              # second coefficient for PS
  # curve(dnorm(x, 0, sqrt(1/0.5)), 0, 20) # howsit look like ?

  # Logistic regression for visits (V)
  for (i in 1:nsites){
    for (t in 1:nyears){
      V[i,t] ~ dbern(theta[i,t])
    }
  }
}
")

# Initial values
# initialize at solutions of model 1 without PS and with a guess of kappa
# (For this, you have to run model 1 beforehand)
tmp <- out1$mean
inits <- function(){list(psi1 = tmp$psi1, initial.phi = tmp$initial.phi,
    initial.gamma = tmp$initial.gamma, sd.eps.lphi = tmp$sd.eps.lphi,
    sd.eps.lgamma = tmp$sd.eps.lgamma, kappa.lphi = 3, kappa.lgamma = 3)}

# Parameters monitored
params <- c("psi1", "alpha.lphi", "initial.phi", "sd.lphi.site",
    "alpha.lgamma", "initial.gamma", "sd.lgamma.site", "lphi.year",
    "lgamma.year", "mean.phi.year", "mean.gamma.year", "sd.eps.lphi",
    "sd.eps.lgamma", "gamma.cliff", "phi.cliff", "alpha.visit", "beta.visit",
    "kappa.lphi", "kappa.lgamma", "n.occ", "n.ain", "n.jura", "n.doubs",
    "lphi.site", "lgamma.site", "y")

# MCMC settings
# na <- 1000 ; ni <- 300000 ; nt <- 150 ; nb <- 150000 ; nc <- 3  # 8.3 hrs
na <- 1000 ; ni <- 3000 ; nt <- 1 ; nb <- 1500 ; nc <- 3 # ~~~ for testing, 8 mins

# Call JAGS (ART 18 h), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "dynocc2.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out2, layout=c(2,2))
summary(out2) ; jags.View(out2) # not shown


# ~~~~ code for models 3 and 4 from MS dated 2019-08-28 ~~~~~~~~~~~~~~~~~~~~~~
# In model 3, we use y[i, t-1] as a predictor for theta[i,t], i.e.,
#   whether or not a pair was detected at that site during the previous year.

# Model 3
# Specify model in BUGS language
cat(file = "dynocc3.txt", "
model {

  # Submodel 1: model for detection/nondetection data (y)
  # ----------------------------------------------------
  # Priors
  psi1 ~ dbeta(1, 1)      # Initial occupancy

  # Model for phi and gamma: cliff height + site + smooth random year effects
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      logit(phi[i,t]) <- lphi[i,t]
      lphi[i,t] <- lphi.site[i] + lphi.year[t]
      logit(gamma[i,t]) <- lgamma[i,t]
      lgamma[i, t] <- lgamma.site[i] + lgamma.year[t]
    }
    lphi.site[i] ~ dnorm(alpha.lphi[height[i]], tau.lphi.site)
    lgamma.site[i] ~ dnorm(alpha.lgamma[height[i]], tau.lgamma.site)
  }
  # Priors for phi and gamma intercepts
  for(k in 1:3){
    alpha.lphi[k] <- logit(initial.phi[k])
    initial.phi[k] ~ dbeta(1, 1)
    alpha.lgamma[k] <- logit(initial.gamma[k])
    initial.gamma[k] ~ dbeta(1, 1)
  }
  tau.lphi.site <- pow(sd.lphi.site, -2)
  sd.lphi.site ~ dunif(0, 2)
  # sd.lphi.site ~ dnorm(0, 1)I(0,)
  tau.lgamma.site <- pow(sd.lgamma.site, -2)
  sd.lgamma.site ~ dunif(0, 2)
  # sd.lgamma.site ~ dnorm(0, 1)I(0,)

  # Priors for year effects on phi and gamma with rw smoothers
  lphi.year[1] <- 0     # Set to zero to avoid overparameterization
  lgamma.year[1] <- 0
  for (t in 2:(nyears-1)){
    lphi.year[t] ~ dnorm(lphi.year[t-1], tau.eps.lphi)
    lgamma.year[t] ~ dnorm(lgamma.year[t-1], tau.eps.lgamma)
  }
  tau.eps.lphi <- pow(sd.eps.lphi,-2)   # Hyperpriors for variances
  sd.eps.lphi ~ dunif(0, 1)
  tau.eps.lgamma <- pow(sd.eps.lgamma,-2)
  sd.eps.lgamma ~ dunif(0, 1)

  # Ecological and observation submodels confounded (no p)
  for (i in 1:nsites){
    y[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      y[i,t] ~ dbern(y[i,t-1]*phi[i,t-1] + (1-y[i,t-1])*gamma[i,t-1])
    }
  }

  # Derived parameters
  # Population occupancy and population size
  psi[1] <- psi1
  n.occ[1]<-sum(y[1:nsites,1])
  for (t in 2:nyears){
    n.occ[t] <- sum(y[1:nsites,t])     # Number of occupied sites
  }
  # Year-specific average values of phi and gamma
  for(t in 1:(nyears-1)){
    mean.phi.year[t] <- mean(phi[,t])
    mean.gamma.year[t] <- mean(gamma[,t])
  }
  # Average gamma and phi per cliff height category in middle year (1989)
  for(k in 1:3){
    logit(gamma.cliff[k]) <- alpha.lgamma[k] + lgamma.year[26]
    logit(phi.cliff[k]) <- alpha.lphi[k] + lphi.year[26]
  }

  # Population size in each departement
  for (t in 1:nyears){
    n.ain[t] <- sum(y[ain, t])
    n.jura[t] <- sum(y[jura, t])
    n.doubs[t] <- sum(y[doubs, t])
  }

  # Submodel 2: model for whether a site is visited or not (V)
  # ---------------------------------------------------------
  # Priors and linear models
  # No effect of past occupation history in visitation at t=1
  for (i in 1:nsites){
    theta[i,1] <- ilogit(alpha.visit + beta.visit[1] * year[1] + beta.visit[2] * pow(year[1],2) + beta.visit[3] * pow(year[1],3))
  }

  # Including an effect of past occupation history in visitation from t=2
  for (i in 1:nsites){
    for (t in 2:nyears){
      theta[i,t] <- ilogit(alpha.visit + beta.visit[1] * year[t] +
        beta.visit[2] * pow(year[t],2) + beta.visit[3] * pow(year[t],3) +
        kappa * y[i,t-1])
    }
  }
  alpha.visit <- logit(theta.int)
  theta.int ~ dbeta(1, 1)
  for(v in 1:3){            # Coefficients of time
    beta.visit[v] ~ dunif(-5, 5)
  }
  kappa ~ dunif(-10, 10)     # coefficient for preferential sampling

  # Logistic regression for visits (V)
  for (i in 1:nsites){
    for (t in 1:nyears){
      V[i,t] ~ dbern(theta[i,t])
    }
  }
}
")

# Initial values
# initialize at solutions of model 1 without PS and with a guess of kappa
tmp <- out1$mean
inits <- function(){list(psi1 = tmp$psi1, initial.phi = tmp$initial.phi,
    initial.gamma = tmp$initial.gamma, sd.eps.lphi = tmp$sd.eps.lphi,
    sd.eps.lgamma = tmp$sd.eps.lgamma, kappa = 4)}

# Parameters monitored
params <- c("psi1", "alpha.lphi", "initial.phi", "alpha.lgamma",
    "initial.gamma", "lphi.year", "lgamma.year", "mean.phi.year",
    "mean.gamma.year", "sd.eps.lphi", "sd.eps.lgamma", "gamma.cliff",
    "phi.cliff", "alpha.visit", "beta.visit", "kappa", "n.occ",
    "n.ain", "n.jura", "n.doubs", "y")

# MCMC settings
# na <- 1000  ;  ni <- 20000  ;  nt <- 30  ;  nb <- 10000  ;  nc <- 3
na <- 1000  ;  ni <- 2000  ;  nt <- 2  ;  nb <- 1000  ;  nc <- 3  # 11 mins

# Call JAGS (ART 94 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "dynocc3.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out3, layout=c(2,2))
summary(out3)   ;   jags.View(out3)   # not shown

# Finally, in model 4, we use as a predictor for theta[i, t] the mean number of
# years over the past 10 in which peregrines were detected at a site. We compute
# this in the model and for this need to supply the range of years in the data,
# which is necessarily fewer than 10 at the start of the study.

# Compute the year range for which site occupancy is averaged
# as a measure of site quality (in the model, year is indexed t)
startYear <- c(1, rep(1, 10), 2:43)
endYear <- c(1, 1:52)     # Fill in initial 1, but not used in BUGS code
lengthPeriod <- endYear-startYear+1
# Look at year range and duration
cbind(1964:2016, startYear, endYear, startYear+1963, endYear+1963, lengthPeriod)
#            startYear endYear           lengthPeriod
#  [1,] 1964         1       1 1964 1964            1  # not used in model
#  [2,] 1965         1       1 1964 1964            1
#  [3,] 1966         1       2 1964 1965            2
#  [4,] 1967         1       3 1964 1966            3
#  [5,] 1968         1       4 1964 1967            4
#  [6,] 1969         1       5 1964 1968            5
#  [7,] 1970         1       6 1964 1969            6
#  [8,] 1971         1       7 1964 1970            7
#  [9,] 1972         1       8 1964 1971            8
# [10,] 1973         1       9 1964 1972            9
# [11,] 1974         1      10 1964 1973           10
# [12,] 1975         2      11 1965 1974           10
# ......
# [52,] 2015        42      51 2005 2014           10
# [53,] 2016        43      52 2006 2015           10

# Data bundle
str(bdata <- list(y = y, height = ht, nsites = nsites, nyears = nyears,
    ain = ain, jura = jura, doubs = doubs, V = V, year = ((1964:2016)-1990)/26,
    startYear = startYear, endYear = endYear))


# Model 4
# Specify model in BUGS language
cat(file = "dynocc4.txt", "
model {

  # Submodel 1: model for detection/nondetection data (z)
  # ----------------------------------------------------
  # Priors
  psi1 ~ dbeta(1, 1)      # Initial occupancy

  # Model for phi and gamma: cliff height + site + smooth random year effects
  for (i in 1:nsites){
    for(t in 1:(nyears-1)){
      logit(phi[i,t]) <- lphi[i,t]
      lphi[i,t] <- lphi.site[i] + lphi.year[t]
      logit(gamma[i,t]) <- lgamma[i,t]
      lgamma[i, t] <- lgamma.site[i] + lgamma.year[t]
    }
    lphi.site[i] ~ dnorm(alpha.lphi[height[i]], tau.lphi.site)
    lgamma.site[i] ~ dnorm(alpha.lgamma[height[i]], tau.lgamma.site)
  }
  # Priors for phi and gamma intercepts
  for(k in 1:3){
    alpha.lphi[k] <- logit(initial.phi[k])
    initial.phi[k] ~ dbeta(1, 1)
    alpha.lgamma[k] <- logit(initial.gamma[k])
    initial.gamma[k] ~ dbeta(1, 1)
  }
  tau.lphi.site <- pow(sd.lphi.site, -2)
  sd.lphi.site ~ dunif(0, 2)
  # sd.lphi.site ~ dnorm(0, 1)I(0,)
  tau.lgamma.site <- pow(sd.lgamma.site, -2)
  sd.lgamma.site ~ dunif(0, 2)
  # sd.lgamma.site ~ dnorm(0, 1)I(0,)

  # Priors for year effects on phi and gamma with rw smoothers
  lphi.year[1] <- 0     # Set to zero to avoid overparameterization
  lgamma.year[1] <- 0
  for (t in 2:(nyears-1)){
    lphi.year[t] ~ dnorm(lphi.year[t-1], tau.eps.lphi)
    lgamma.year[t] ~ dnorm(lgamma.year[t-1], tau.eps.lgamma)
  }
  tau.eps.lphi <- pow(sd.eps.lphi,-2)   # Hyperpriors for variances
  sd.eps.lphi ~ dunif(0, 1)
  tau.eps.lgamma <- pow(sd.eps.lgamma,-2)
  sd.eps.lgamma ~ dunif(0, 1)

  # Ecological and observation submodels confounded (no p)
  for (i in 1:nsites){
    y[i,1] ~ dbern(psi1)
    for (t in 2:nyears){
      y[i,t] ~ dbern(y[i,t-1]*phi[i,t-1] + (1-y[i,t-1])*gamma[i,t-1])
    }
  }

  # Derived parameters
  # Population occupancy and population size
  psi[1] <- psi1
  n.occ[1]<-sum(y[1:nsites,1])
  for (t in 2:nyears){
    n.occ[t] <- sum(y[1:nsites,t])     # Number of occupied sites
  }
  # Year-specific average values of phi and gamma
  for(t in 1:(nyears-1)){
    mean.phi.year[t] <- mean(phi[,t])
    mean.gamma.year[t] <- mean(gamma[,t])
  }
  # Average gamma and phi per cliff height category in middle year (1989)
  for(k in 1:3){
    logit(gamma.cliff[k]) <- alpha.lgamma[k] + lgamma.year[26]
    logit(phi.cliff[k]) <- alpha.lphi[k] + lphi.year[26]
  }

  # Running past sum of years occupied = measure of current site quality, to be used in the model for visitation
  # Compute the proportion of years (during the last 10) in which a site was occupied (has to be a proportion to take account of the first 10 years)
  for (i in 1:nsites){
    for(t in 1:nyears){ # Value for year 1 is nonsense and not used below
      OccPast10[i,t] <- mean(y[i,startYear[t]:endYear[t]])
    }
  }

  # Population size in each departement
  for (t in 1:nyears){
    n.ain[t] <- sum(y[ain, t])
    n.jura[t] <- sum(y[jura, t])
    n.doubs[t] <- sum(y[doubs, t])
  }

  # Submodel 2: model for whether a site is visited or not (V)
  # ---------------------------------------------------------
  # Priors and linear models
  # No effect of past occupation history in visitation at t=1
  for (i in 1:nsites){
    theta[i,1] <- ilogit(alpha.visit + beta.visit[1] * year[1] + beta.visit[2] * pow(year[1],2) + beta.visit[3] * pow(year[1],3))
  }

  # Including an effect of past occupation history in visitation from t=2
  for (i in 1:nsites){
    for (t in 2:nyears){
      theta[i,t] <- ilogit(alpha.visit + beta.visit[1] * year[t] +
        beta.visit[2] * pow(year[t],2) + beta.visit[3] * pow(year[t],3) +
        kappa10 * OccPast10[i,t])
    }
  }
  alpha.visit <- logit(theta.int)
  theta.int ~ dbeta(1, 1)
  for(v in 1:3){            # Coefficients of time
    beta.visit[v] ~ dunif(-5, 5)
  }
  kappa10 ~ dunif(-10, 10)     # coefficient for preferential sampling

  # Logistic regression for visits (V)
  for (i in 1:nsites){
    for (t in 1:nyears){
      V[i,t] ~ dbern(theta[i,t])
    }
  }
}
")

# Initial values
# initialize at solutions of model 1 without PS and with guess of kappa10
tmp <- out1$mean
inits <- function(){list(psi1 = tmp$psi1, initial.phi = tmp$initial.phi,
    initial.gamma = tmp$initial.gamma, sd.eps.lphi = tmp$sd.eps.lphi,
    sd.eps.lgamma = tmp$sd.eps.lgamma, kappa10 = 4)}

# Parameters monitored
params <- c("psi1", "alpha.lphi", "initial.phi", "alpha.lgamma", "initial.gamma",
    "lphi.year", "lgamma.year", "mean.phi.year", "mean.gamma.year", "sd.eps.lphi",
    "sd.eps.lgamma", "gamma.cliff", "phi.cliff", "alpha.visit", "beta.visit",
    "kappa10", "n.occ", "n.ain", "n.jura", "n.doubs", "y")

# MCMC settings
# na <- 1000  ;  ni <- 20000  ;  nt <- 10  ;  nb <- 10000  ;  nc <- 3
na <- 1000  ;  ni <- 2000  ;  nt <- 1  ;  nb <- 1000  ;  nc <- 3  # ~~~ for testing, 17 mins

# Call JAGS (ART 156 min), check convergence and summarize posteriors
out4 <- jags(bdata, inits, params, "dynocc4.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out4, layout=c(2,2))
summary(out4)   ;   jags.View(out4)  # not shown


# ~~~~~~~~~ code for figures 4.32 ~~~~~~~~~~~~~
# Compare population size estimates (fig 4.32)
plot(year, n.occ.obs, xlab = 'Year', ylab = 'Population size',
    type = 'b', pch = 16, ylim = c(0, 250), frame = FALSE)
points(year, n.ratio, pch = 1, type = 'b', col = 'black')
points(year, out1$mean$n.occ, pch = '1', type = 'b')
points(year, out2$mean$n.occ, pch = '2', type = 'b')
points(year, out3$mean$n.occ, pch = '3', type = 'b')
points(year, out4$mean$n.occ, pch = '4', type = 'b')
legend('bottomright', c('Ratio estimator', 'Observed number of pairs'),
    lty = 1, pch = c(1, 16), col = 'black', bty = 'n')

# Compare estimates of persistence and colonization rates (Fig. 4-33)
phi.ma <- apply(rbind(out2$sims.list$mean.phi.year,
    out3$sims.list$mean.phi.year), 2, mean)   # Model average annual phi
gamma.ma <- apply(rbind(out2$sims.list$mean.gamma.year,
    out3$sims.list$mean.gamma.year), 2, mean) # Model average annual phi

plot(year[-53], out1$mean$mean.phi.year, xlab = 'Year', ylab = 'Probability',
    type = 'b', pch = 0, ylim = c(0, 1), frame = FALSE)
points(year[-53], phi.ma, pch = 15, type = 'b')
points(year[-53], out1$mean$mean.gamma.year,  pch = 1, type = 'b')
points(year[-53], gamma.ma,  pch = 16, type = 'b')
legend('bottomright', c('Persistence (without PS)', 'Persistence (with PS)',
    'Colonization (without PS)', 'Colonization (with PS)'), lty = 1,
    pch = c(0,15,1,16), bty = 'n')

# PS-model-average y estimates and plot them: Fig. 4–34
library(abind)
ysims <- abind(out2$sims.list$y, out3$sims.list$y, out4$sims.list$y, along = 1)
yhat <- apply(ysims, 2:3, mean)
mapPalette <- colorRampPalette(c("white", "black"))
image(x = 1964:2016, y = 1:nsites, z = t(yhat), col = mapPalette(10),
    axes = TRUE, xlab = "Year", ylab = "Site", main = '')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Population size comparisons in a table (not all printed out)
n.occ.ma <- colMeans(rbind(out2$sims.list$n.occ, out3$sims.list$n.occ,
    out3$sims.list$n.occ)) # Model average population size estimates from 3 PS models
print(cbind(n.visited, n.occ.obs, n.ratio, 'model1' = out1$mean$n.occ,
    'model2' = out2$mean$n.occ, 'model3' = out3$mean$n.occ, 'model4' = out4$mean$n.occ,
    'PS average' = n.occ.ma, 'missed' = n.occ.ma - n.occ.obs,
    'Prop.detected' = n.occ.obs / n.occ.ma), 2)
#        n.visited n.occ.obs n.ratio model1 model2 model3 model4 PS average missed Prop.detected
# yr1964       100        34      97     95     98     52     75         67   33.3          0.51
# yr1965       137        45      93     83     75     49     56         57   12.3          0.79
# yr1966       149        39      74     68     59     41     44         47    7.7          0.83
# yr1967       160        36      64     57     50     37     39         41    5.0          0.88
# yr1968       168        20      34     35     30     22     22         25    4.5          0.82
# yr1969       172        18      30     28     25     18     18         20    2.3          0.89
# yr1970       171        21      35     30     28     22     21         24    2.8          0.88
# yr2010       249       186     212    203    202    190    192        194    7.8          0.96
# yr2011       246       184     212    203    201    190    191        193    9.4          0.95
# yr2012       248       191     219    211    208    196    199        200    9.4          0.95
# yr2013       235       173     209    198    196    181    186        186   12.8          0.93
# yr2014       244       183     213    205    202    186    194        192    8.5          0.96
# yr2015       245       158     183    179    175    163    171        167    9.0          0.95
# yr2016       234       169     205    194    191    184    189        186   17.2          0.91


# PS-model-average cliff height effects in phi and gamma
# ~~~~ code to produce the table ~~~~~~~~~~~~~~~~~~~~~~
p.cl <- rbind(out2$sims.list$phi.cliff, out3$sims.list$ phi.cliff,
    out4$sims.list$ phi.cliff)
g.cl <- rbind(out2$sims.list$gamma.cliff, out3$sims.list$ gamma.cliff,
    out4$sims.list$ gamma.cliff)
cliff.tab <- cbind('post.mean phi' = apply(p.cl, 2, mean),
    'CRI' = t(apply(p.cl, 2, quantile, probs=c(0.025, 0.975))),
    'post.mean gamma' = apply(g.cl, 2, mean),
    'CRI' = t(apply(g.cl, 2, quantile, probs=c(0.025, 0.975))))
rownames(cliff.tab) <- c('low', 'medium', 'tall')
print(cliff.tab, 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model-averaged phi and gamma per cliff-height category
#        post.mean phi  2.5% 97.5% post.mean gamma   2.5% 97.5%
# low            0.897 0.842 0.937           0.162 0.0878 0.312
# medium         0.886 0.797 0.938           0.143 0.0597 0.342
# tall           0.943 0.914 0.966           0.305 0.1788 0.549
