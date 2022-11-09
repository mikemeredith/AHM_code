#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10 : INTEGRATED MODELS FOR MULTIPLE TYPES OF DATA
# =========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 25 mins
# Run time with the full number of iterations: 2.6 hrs

library(AHMbook)
library(jagsUI)

# 10.8 Example 6: A spatial integrated model combining spatial
#      capture-recapture data and counts or occupancy data
#      from unmarked individuals
# ====================================================

# 10.8.1 A brief introduction to spatial capture-recapture (no code)

# 10.8.2 Simulating spatial capture-recapture data, and counts of
#        unmarked individuals as a sort of “degraded” SCR data
# ----------------------------------------------------------------

library(AHMbook)

# Make a trapping grid of 64 hair snags
traplocs <- expand.grid(1:10, 1:10)
ntraps <- nrow(traplocs)

# Define state-space of point process by buffering traps
delta <- 2                          # Buffer width
Xl <- min(traplocs[,1] - delta)     # Lower x
Xu <- max(traplocs[,1] + delta)     # Upper x
Yl <- min(traplocs[,2] - delta)     # Lower y
Yu <- max(traplocs[,2] + delta)     # Upper y

# Distribute population of size N = 50 in the state space
N <- 50
K <- 4            # K = number of weeks of hair collection ("effort")

# Simulate activity centers
set.seed(121, kind = "Mersenne-Twister")    # Initialize RNGs
sx <- runif(N, Xl, Xu)                      # x coordinate
sy <- runif(N, Yl, Yu)                      # y coordinate
smat <- cbind(sx, sy)                       # Activity center x,y coordinates

# Only some traps collect SCR data ...
inner36 <- traplocs[,1] >= 3 & traplocs[,1] <= 8 & traplocs[,2] >= 3 &
    traplocs[,2] <= 8
# ... while others only collect count data
outer64 <- 1:ntraps
outer64 <- outer64[!inner36]

# Plot the state space, trap locations and individual ACs (Fig. 10.9)
plot(c(Xl,Xu), c(Yl,Yu), type = 'n', xlab = "x", ylab = "y", asp = 1)
points(smat, pch = '+', cex = 1.6)
points(traplocs, pch = 0)
points(traplocs[inner36,], pch=15)

# Parameters for the SCR model simulation
lam0 <- 0.6                       # Baseline encounter rate
sigma <- 0.6                      # Half-normal detection scale

# Generate the encounters of every individual in every trap
D <- e2dist(smat, traplocs)       # Ind/trap distance matrix
muy <- lam0 * exp(-(D*D) / (2 * sigma^2))
Y <- matrix(NA, nrow = N, ncol = ntraps)
for(i in 1:N){
  Y[i,] <- rpois(ntraps, K*muy[i,])
}

# Now take the SCR data from the inner 16 traps
Yscr <- Y[, inner36]
# Captured individuals appear in the data set
totalcaps <- apply(Yscr, 1, sum)
Yscr <- Yscr[totalcaps > 0,]
# Only count (i.e., 'degraded') data appear in the other data set
Yocc <- Y[, outer64]
# Total number of detections observed in outer traps
n <- apply(Yocc, 2, sum)

# Trap arrays for both data types
scrtraps <- traplocs[inner36,]
occtraps <- traplocs[outer64,]

# 10.8.3 Analysis of SCR data using an SCR model: the basic SCR0 model
#        (model 1)
# ---------------------------------------------------------------------

# Set up data augmentation for the encounter histories
M <- 150
Yaug <- matrix(0, M, dim(scrtraps)[1])
Yaug[1:dim(Yscr)[1],] <- Yscr

# Bundle and summarize data set
str(bdata <- list(y = Yaug, M = M, nsurveys = K, ntraps = nrow(scrtraps),
    Xl = Xl, Yl = Yl, Xu = Xu, Yu = Yu, X = as.matrix(scrtraps)),1)
# List of 9
# $ y       : num [1:150, 1:36] 0 0 1 0 0 0 0 0 0 0 ...
# $ M       : num 150
# $ nsurveys: num 4
# $ ntraps  : int 36
# $ Xl      : num -1
# $ Yl      : num -1
# $ Xu      : num 12
# $ Yu      : num 12
# $ X       : int [1:36, 1:2] 3 4 5 6 7 8 3 4 5 6 ...

# Specify simple SCR0 model for SCR data in BUGS language
cat(file = "SCR0.txt", "
model {
  # Priors
  lam0 ~ dunif(0,5)
  sigma ~ dunif(0, 10)
  psi ~ dunif(0,1)

  # 'Likelihood'
  for(i in 1:M){
    # Process model
    z[i] ~ dbern(psi) # 'Existence' of individual ...
    s[i,1] ~ dunif(Xl, Xu) # ... and its location: x ...
    s[i,2] ~ dunif(Yl, Yu) # .... and y coordinate of activity centre
    # Observation model
    for(j in 1:ntraps){
      d[i,j] <- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
      lambda[i,j] <- z[i]*lam0*exp(-(d[i,j]*d[i,j])/(2*sigma*sigma))
      y[i,j] ~ dpois(nsurveys * lambda[i,j])
    }
  }
  N <- sum(z[]) # Population size as a derived quantity
}
")

# Initial values
# inits for individual ACs: mean observed capture location
SinX <- runif(M, Xl, Xu)
SinY <- runif(M, Yl, Yu)
for(i in 1:dim(Yscr)[1]){
    SinX[i] <- sum(Yscr[i,]*scrtraps[,1])/(sum(Yscr[i,]))
    SinY[i] <- sum(Yscr[i,]*scrtraps[,2])/(sum(Yscr[i,]))
}
inits <- function() {list(z = c(rep(1, dim(Y)[1]), rbinom(M-dim(Y)[1], 1, 0.5)),
    psi = runif(1), s = cbind(SinX, SinY), lam0 = runif(1, 0.5, 1.5),
    sigma = runif(1, 0.5, 3)) }

# Params monitored
params <- c('psi', 'lam0', 'sigma', 'N')

# MCMC settings
# na <- 1000 ; ni <- 30000; nt <- 20 ; nb <- 10000 ; nc <- 3
na <- 1000 ; ni <- 3000; nt <- 2 ; nb <- 1000 ; nc <- 3  # ~~~ for testing, 3 mins

# Call JAGS (ART 15 min), assess convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "SCR0.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out1, layout=c(2,3))
print(out1, 3)
#         mean     sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
# psi    0.423  0.096  0.252  0.415  0.627    FALSE 1 1.002   818
# lam0   0.554  0.127  0.344  0.541  0.832    FALSE 1 1.000  3000
# sigma  0.573  0.060  0.472  0.568  0.706    FALSE 1 1.001  2163
# N     63.380 13.565 40.000 62.000 93.000    FALSE 1 1.002   853


# 10.8.4 Analysis of counts of unmarked individuals using an SCR
#        model: the Chandler-Royle model (model 2)
# ---------------------------------------------------------------

# Bundle and summarize data set
str(bdata <- list(n = n, M = M, nsurveys = K, ntraps = nrow(occtraps), Xl = Xl,
    Yl = Yl, Xu = Xu, Yu = Yu, X = as.matrix(occtraps)))
# List of 9
# $ n       : int [1:64] 2 0 0 5 2 7 3 1 0 0 ...
# $ M       : num 150
# $ nsurveys: num 4
# $ ntraps  : int 64
# $ Xl      : num -1
# $ Yl      : num -1
# $ Xu      : num 12
# $ Yu      : num 12
# $ X       : int [1:64, 1:2] 1 2 3 4 5 6 7 8 9 10 ...

# Specify basic Chandler-Royle (2013) model for counts in BUGS language
cat(file = "CRmodel.txt", "
model{
  # Priors
  lam0 ~ dunif(0,5)
  sigma ~ dunif(0, 10)
  psi ~ dunif(0,1)

  # 'Likelihood'
  for(i in 1:M) {
    # Process model
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl, Xu)
    s[i,2] ~ dunif(Yl, Yu)
  }
  # Observation model
  for(j in 1:ntraps) {
    for(i in 1:M) {
      d2[i,j] <- pow(s[i,1] - X[j,1],2) + pow(s[i,2] - X[j,2],2)
      lam[i,j] <- lam0*exp(-(d2[i,j])/(2*sigma*sigma))*z[i]
    }
    lambda[j] <- sum(lam[,j])     # Expected captures/occasion at trap j
    n[j] ~ dpois(nsurveys * lambda[j])
  }
  N <- sum(z[])                   # Population size as derived quantity
}
")

# Initial values
Sx <- runif(M, Xl, Xu)            # Inits for activity centers
Sy <- runif(M, Yl, Yu)
inits <- function() {list(z = as.vector(rep(1, M)), psi = runif(1),
    s = cbind(Sx, Sy), lam0 = runif(1, 0.5, 1.5), sigma = runif(1, 0.5, 3))
}

# Parameters monitored
params <- c('psi', 'lam0', 'sigma', 'N')

# MCMC settings
# na <- 1000 ; ni <- 20000; nt <- 3 ; nb <- 5000 ; nc <- 4
na <- 1000 ; ni <- 2000; nt <- 3 ; nb <- 500 ; nc <- 3  # ~~~ for testing, 5 mins

# Call JAGS (ART 49 min), assess convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "CRmodel.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out2, layout=c(2,3))
print(out2, 3)
#         mean     sd   2.5%    50%   97.5% overlap0 f  Rhat n.eff
# psi    0.622  0.236  0.113  0.643   0.979    FALSE 1 1.003   985
# lam0   0.638  0.221  0.263  0.618   1.119    FALSE 1 1.007   516
# sigma  0.505  0.215  0.303  0.442   1.192    FALSE 1 1.018   443
# N     93.634 35.409 17.000 97.000 148.000    FALSE 1 1.003   996


# ~~~ Extra code to run the model with M = 300 ~~~~~~~~
bdata300 <- list(n = n, M = 300, nsurveys = K, ntraps = nrow(occtraps), Xl = Xl,
    Yl = Yl, Xu = Xu, Yu = Yu, X = as.matrix(occtraps))
# Initial values
Sx <- runif(300, Xl, Xu) # Inits for activity centers
Sy <- runif(300, Yl, Yu)
inits <- function() {list(z = as.vector(rep(1, 300)), psi = runif(1),
    s = cbind(Sx, Sy), lam0 = runif(1, 0.5, 1.5), sigma = runif(1, 0.5, 3))
}
out2.300 <- jags(bdata300, inits, params, "CRmodel.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(out2.300, 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Results for M = 300:
#          mean     sd   2.5%     50%   97.5% overlap0 f  Rhat n.eff
# psi     0.495  0.233  0.086   0.466   0.955    FALSE 1 1.013   197
# lam0    0.634  0.223  0.257   0.612   1.137    FALSE 1 1.002  2009
# sigma   0.415  0.191  0.223   0.375   1.045    FALSE 1 1.044   165
# N     148.478 70.032 26.000 140.000 287.000    FALSE 1 1.013   194

# Remember what the true values were ...
(truth <- c('N' = N, 'lam0' = lam0, 'sigma' = sigma))
#    N lam0 sigma
# 50.0  0.6   0.6


# 10.8.5 Model 3: The integrated model for joint analysis of SCR and count data
# -----------------------------------------------------------------------------

# Bundle and summarize data set
str(bdata <- list(yscr = Yaug, n = n, M = M, nsurveys = K,
    n.scrtraps = nrow(scrtraps), n.occtraps = length(n),
    Xl = Xl, Yl = Yl, Xu = Xu, Yu = Yu,
    scrtraps = as.matrix(scrtraps), occtraps = as.matrix(occtraps)))
# List of 12
# $ yscr : int [1:150, 1:36] 0 0 1 0 0 0 0 0 0 0 ...
# $ n         : int [1:64] 2 0 0 5 2 7 3 1 0 0 ...
# $ M         : num 150
# $ nsurveys  : num 4
# $ n.scrtraps: int 36
# $ n.occtraps: int 64
# $ Xl        : num -1
# $ Yl        : num -1
# $ Xu        : num 12
# $ Yu        : num 12
# $ scrtraps  : int [1:36, 1:2] 3 4 5 6 7 8 3 4 5 6 ...
# $ occtraps  : int [1:64, 1:2] 1 2 3 4 5 6 7 8 9 10 ...

# Specify integrated Chandler-Royle+SCR model in BUGS language
cat(file = "IM_SCR.txt", "
model {
  # Priors
  psi ~ dunif(0, 1)
  lam0 ~ dunif(0, 5)
  sigma ~ dunif(0, 10)

  # 'Likelihood'
  for(i in 1:M){
    # Shared process model for both data sets
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl,Xu)
    s[i,2] ~ dunif(Yl,Yu)
    # Observation model for SCR data
    for(j in 1:n.scrtraps){
      d[i,j] <- pow(pow(s[i,1] - scrtraps[j,1],2) +
          pow(s[i,2] - scrtraps[j,2],2),0.5)
      lambda1[i,j] <- z[i]*lam0*exp(-(d[i,j]*d[i,j])/(2*sigma*sigma))
      yscr[i,j] ~ dpois(nsurveys * lambda1[i,j])  # Response 1: SCR data
    }
  }
  # Observation model for unmarked counts
  for(j in 1:n.occtraps) {
    for(i in 1:M){
      d2[i,j] <- pow(s[i,1] - occtraps[j,1],2) + pow(s[i,2] - occtraps[j,2],2)
      lam[i,j] <- z[i]*lam0*exp(-(d2[i,j])/(2*sigma*sigma))
    }
    lambda2[j] <- sum(lam[,j])          # Expected captures/occasion at trap t
    n[j] ~ dpois(nsurveys * lambda2[j]) # Response 2: unmarked counts
  }
  N <- sum(z[])                   # Population size as a derived quantity
}
")

# Inits for activity centers
SinX <- SinY <- rep(NA,M)
for(i in 1:dim(Yscr)[1]){
  SinX[i] <- sum(Yscr[i,]*scrtraps[,1])/(sum(Yscr[i,]))
  SinY[i] <- sum(Yscr[i,]*scrtraps[,2])/(sum(Yscr[i,]))
}
inits <- function() {list(z = c(rep(1, dim(Y)[1]),
    rbinom(M - dim(Y)[1], 1, 0.5)),
    psi = runif(1), s = cbind(SinX, SinY), sigma = runif(1, 0.5, 3)) }

# Parameters monitored
params <- c('psi', 'lam0', 'sigma', 'N')

# MCMC settings
# na <- 1000 ; ni <- 10000; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~ for testing, 6 mins

# Call JAGS (ART 36 min), assess convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "IM_SCR.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out3, layout=c(2,3))
print(out3, 3)
#         mean     sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
# psi    0.394  0.081  0.250  0.392  0.563    FALSE 1 1.001  1181
# lam0   0.559  0.103  0.374  0.553  0.780    FALSE 1 1.002   909
# sigma  0.580  0.060  0.480  0.574  0.717    FALSE 1 1.002   888
# N     58.856 10.582 40.000 58.000 81.000    FALSE 1 1.001  1795

# ~~~ extra code for figure 10.10 ~~~~~~~~~~~
op <- par(mfrow = c(1, 3))
hist(out2$sims.list$N, col = 'grey', xlab = 'N', ylab = 'Density',
    main = 'CR model, M = 150',
    freq = FALSE, axes = FALSE)
axis(1) ; axis(2)

hist(out2.300$sims.list$N, col = 'grey', xlab = 'N', ylab = 'Density',
    main = 'CR model, M = 300', freq = FALSE, axes = FALSE)
axis(1) ; axis(2)

hist(out3$sims.list$N, col = 'grey', xlab = 'N', ylab = 'Densiry',
    main = 'Integrated model, M = 150', freq = FALSE, axes = FALSE,
    xlim = c(0, 150))
axis(1) ; axis(2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 10.8.6 Model 4: The integrated model for joint analysis of SCR and
#        occupancy data
# -------------------------------------------------------------------

# Generate the encounters of every individual in every trap
set.seed(121, kind = "Mersenne-Twister") # Initialize RNGs
D <- e2dist(smat, traplocs) # Ind/trap distance matrix
muy <- lam0 * exp(-(D*D) / (2 * sigma^2))
Y <- array(NA, dim = c(N, ntraps, K) )
for(i in 1:N){
  for(k in 1:K){
    Y[i,,k] <- rpois(ntraps, muy[i,])
  }
}

# Now take the SCR data from the inner 36 traps
Yscr <- Y[, inner36, ]
# Captured individuals appear in the data set
totalcaps <- apply(Yscr, 1, sum)
nind <- sum(totalcaps>0)
Yscr <- Yscr[totalcaps > 0,,]

# Only count (i.e., 'degraded') data appear in the other data set
Yocc <- Y[, outer64, ]
# Total number of detections observed in outer traps
Yocc <- apply(Yocc, c(2,3), sum)

# Set up data augmentation for the encounter histories
M <- 150
Yaug <- array(0, dim = c(M, dim(scrtraps)[1],K))
Yaug[1:dim(Yscr)[1],,] <- Yscr

# Now convert to binary presence/absence (or occupancy) data:
Yaug[Yaug > 1] <- 1
Yocc[Yocc > 1] <- 1

# Bundle and summarize data set
str(bdata <- list(yscr = Yaug, yocc=Yocc, M = M, nsurveys = K,
    n.scrtraps = nrow(scrtraps), n.occtraps = length(n), Xl = Xl, Yl = Yl,
    Xu = Xu, Yu = Yu, scrtraps = as.matrix(scrtraps),
    occtraps = as.matrix(occtraps)) ,1)
# List of 12
# $ yscr      : num [1:150, 1:36, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
# $ yocc      : num [1:64, 1:4] 0 0 0 1 0 0 0 0 0 1 ...
# $ M         : num 150
# $ nsurveys  : num 4
# $ n.scrtraps: int 36
# $ n.occtraps: int 64
# $ Xl        : num -1
# $ Yl        : num -1
# $ Xu        : num 12
# $ Yu        : num 12
# $ scrtraps  : int [1:36, 1:2] 3 4 5 6 7 8 3 4 5 6 ...
# $ occtraps  : int [1:64, 1:2] 1 2 3 4 5 6 7 8 9 10 ...

# Specify integrated SCR/occupancy model in BUGS language
cat(file = "SCRocc.txt","
model{
  psi ~ dbeta(1, 1)
  lam0 ~ dunif(0, 5)
  lam0occ ~ dunif(0,5)
  sigma ~ dunif(0, 10)
  N <- sum(z[])
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl, Xu)
    s[i,2] ~ dunif(Yl, Yu)
    # Compute detection probability for SCR
    for(j in 1:n.scrtraps) {
      d2scr[i,j] <- (s[i,1]-scrtraps[j,1])^2 + (s[i,2]-scrtraps[j,2])^2
      lam[i,j] <- lam0*exp(-d2scr[i,j] / (2*sigma^2))
      # Note cloglog link here: Could use half-normal model too
      pscr[i,j] <- 1 - exp(-lam[i,j]) ####p0scr*exp(-d2scr[i,j] / (2*sigma^2))
    }
    # Compute detection probability for occupancy
    for(j in 1:n.occtraps) {
      d2occ[i,j] <- (s[i,1] - occtraps[j,1])^2 + (s[i,2] - occtraps[j,2])^2
      # Note: cloglog link here: could use half-normal model too
      # pocc[i,j] <- lam0occ*exp(-d2occ[i,j] / (2*sigma^2))
      lam2[i,j] <- lam0occ*exp(-d2occ[i,j] /(2*sigma^2))
      pocc[i,j] <- 1 - exp(-lam2[i,j])
    }
    for(j in 1:n.scrtraps) {
      for(k in 1:nsurveys) {
        # SCR encounters
        yscr[i,j,k] ~ dbern(pscr[i,j]*z[i])
      }
    }
    for(j in 1:n.occtraps){
      # for PA data compute probability of not captured
      pn[i,j] <- (1 - (pocc[i,j]*z[i]) )
    }
  }
  # Model for the presence-absence data
  for(j in 1:n.occtraps) {
    for(k in 1:nsurveys) {
      yocc[j,k] ~ dbern(1-prod(pn[,j])) #
    }
  }
}
")

# Inits
SinX <- runif(M, Xl, Xu)
SinY <- runif(M, Yl, Yu)
ncaps <- apply(Yscr, c(1,2), sum)
ncaps[ncaps>1] <- 1
for(i in 1:nind){
  SinX[i] <- sum(ncaps[i,]*scrtraps[,1])/(sum(ncaps[i,]))
  SinY[i] <- sum(ncaps[i,]*scrtraps[,2])/(sum(ncaps[i,]))
}
inits <- function() {list(z = c(rep(1, dim(Y)[1]), rbinom(M-dim(Y)[1], 1, 0.5)),
    psi = runif(1), s = cbind(SinX, SinY), sigma = runif(1, 0.5, 3),
    lam0occ = runif(1,0.2,0.5) ) }

# Parameters monitored
params <- c('psi', 'lam0', 'lam0occ', 'sigma', 'N')

# MCMC settings
# na <- 1000 ; ni <- 12000; nt <- 2 ; nb <- 2000 ; nc <- 3
na <- 1000 ; ni <- 1200; nt <- 1 ; nb <- 200 ; nc <- 3  # ~~~ for testing, 16 mins

# Call JAGS (ART 63 min), assess convergence and summarize posteriors
out4 <- jags(bdata, inits, params, "SCRocc.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out4, layout=c(2,3))
print(out4, 3)
#           mean     sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
# psi      0.406  0.077  0.270  0.402  0.571    FALSE 1 1.000 15000
# lam0     0.812  0.165  0.533  0.798  1.171    FALSE 1 1.000 15000
# lam0occ  0.573  0.143  0.342  0.559  0.893    FALSE 1 1.000  8749
# sigma    0.574  0.044  0.496  0.570  0.670    FALSE 1 1.001 10135
# N       60.731 10.036 44.000 60.000 83.000    FALSE 1 1.000 15000


# 10.8.7 Concluding comments on data integration in the context of
#   spatial capture-recapture models (no code)
