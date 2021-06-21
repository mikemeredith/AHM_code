#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using nsites and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 15 mins

library(AHMbook)
library(jagsUI)

# 2.9 The multistate Dail-Madsen model
# ====================================

# 2.9.1 Example 1: diseased frogs
# -------------------------------
# Specify model in BUGS language
cat(file = "frogs.txt", "
model {
  # Priors
  # ------- NOT Infected
  for(t in 1:(nyears-1)) {
    gammaN[t] <- mean.gammaN
    omegaN[t] <- mean.omegaN
    psi_NI[t] <- mean.psi_NI
  }
  alpha.lamN ~ dnorm(0,0.01)
  mean.pN ~ dunif(0,1)
  mean.gammaN ~ dnorm(0,0.01)
  mean.omegaN ~ dunif(0,1)
  mean.psi_NI ~ dunif(0,1)
  # ------- Infected
  for(t in 1:(nyears-1)) {
    gammaI[t] <- mean.gammaI
    omegaI[t] <- mean.omegaI
    psi_IN[t] <- mean.psi_IN
  }
  alpha.lamI ~ dnorm(0,0.01)
  mean.pI ~ dunif(0,1)
  mean.gammaI ~ dnorm(0,0.01)
  mean.omegaI ~ dunif(0,1)
  mean.psi_IN ~ dunif(0,1)
  #------ Detection
  for(t in 1:nyears) {
    pN[t] <- mean.pN
    pI[t] <- mean.pI
  }
  # Ecological state model
  #---- First year
  for(i in 1:nsites) {
    #------ Not infected
    NN[i, 1] ~ dpois(lambdaN[i])
    log(lambdaN[i]) <- alpha.lamN
    #----- Infected
    NI[i, 1] ~ dpois(lambdaI[i])
    log(lambdaI[i]) <- alpha.lamI
  } # end i
  #------ Second and subsequent years
  for(t in 2:nyears) {
    for(i in 1:nsites) {
    #------- NOT Infected
      SN[i,t] ~ dbin( omegaN[t-1], NN[i,t-1]) # Total survivors
      TN[i,t] ~ dbin( psi_NI[t-1], SN[i,t]) # Survive, become infected
      GN[i,t] ~ dpois(GaN[i, t])
      log(GaN[i, t]) <- gammaN[t-1]
      #------- Infected
      SI[i,t] ~ dbin(omegaI[t-1], NI[i,t-1] ) # Infecteds who survive
      TI[i,t] ~ dbin(psi_IN[t-1], SI[i,t] ) # Recover: transition to uninfected
      GI[i,t] ~ dpois(GaI[i, t])
      log(GaI[i, t]) <- gammaI[t-1]
      # Subsequent population size, add everything up and subtract out-transitions
      NN[i, t] <- SN[i,t] - TN[i,t] + GN[i,t] + TI[i,t]
      NI[i, t] <- SI[i,t] - TI[i,t] + GI[i,t] + TN[i,t]
    } # end i
  } # end t
  # Observation model
  for(i in 1:nsites) {
    for(j in 1:nsurveys) {
      for(t in 1:nyears) {
        yN[i, j, t] ~ dbin(pN[t], NN[i, t])
        yI[i, j, t] ~ dbin(pI[t], NI[i, t])
        # Predictions of new data
        y.newN[i, j, t] ~ dbin(pN[t], NN[i, t])
        y.newI[i, j, t] ~ dbin(pI[t], NI[i, t])
      } # end t
    } # end j
  } # end i
  # Bayesian p-value calculation
  for(t in 1:nyears) {
    for(i in 1:nsites) {
      evalN[i,t] <- pN[t] * NN[i,t]
      evalI[i,t] <- pI[t] * NI[i,t]
      for(j in 1:nsurveys) {
        EN[i,j,t] <- pow((yN[i,j,t] - evalN[i,t]),2) / (evalN[i,t] + 0.5)
        EI[i,j,t] <- pow((yI[i,j,t] - evalI[i,t]),2) / (evalI[i,t] + 0.5)
        E.newN[i,j,t] <- pow((y.newN[i,j,t] - evalN[i,t]),2) / (evalN[i,t] + 0.5)
        E.newI[i,j,t] <- pow((y.newI[i,j,t] - evalI[i,t]),2) / (evalI[i,t] + 0.5)
      } # end j
    } # end i
  } # end t
  fitN <- sum(EN[,,])
  fitN.new <- sum(E.newN[,,])
  fitI <- sum(EI[,,])
  fitI.new <- sum(E.newI[,,])
} # end model
")


# Simulate a data set (using function from book website)
set.seed(2019)
str(sodata <- simFrogDisease(nsites = 100, nyears = 3, nsurveys = 3,
    alpha.lam = 3, omega = c(0.9, 0.7), gamma = c(2,1),
    p = c(0.8, 0.8, 0.8), recovery = 0.1, infection = 0.1))

# Bundle the data
str(bdata <- list(yN = sodata$yN, yI = sodata$yI,
    nsites = dim(sodata$yN)[1], nsurveys = dim(sodata$yN)[2],
    nyears = dim(sodata$yN)[3]))
# List of 5
# $ yN       : int [1:100, 1:3, 1:3] 4 3 1 3 1 0 4 0 1 2 ...
# $ yI       : int [1:100, 1:3, 1:3] 2 3 2 5 1 3 2 3 0 1 ...
# $ nsites   : int 100
# $ nsurveys : int 3
# $ nyears   : int 3

# Initial values
inits <- function() {list(alpha.lamN = runif(1, 2, 3), mean.pN = runif(1, 0.9, 1),
    mean.omegaN = runif(1, 0.7, 1), mean.gammaN = runif(1, 2, 3),
    mean.psi_NI = runif(1, 0, 0.3), alpha.lamI = runif(1, 2, 3),
    mean.pI = runif(1, 0.9, 1), mean.omegaI = runif(1, 0.7, 1),
    mean.gammaI = runif(1, 2, 3), mean.psi_IN = runif(1, 0, 0.3))}

# Parameters monitored
params <- c("alpha.lamN", "alpha.lamI", "mean.pN", "mean.pI",
    "mean.omegaN",  "mean.omegaI", "mean.gammaN", "mean.gammaI",
    "mean.psi_NI", "mean.psi_IN", "fitN", "fitN.new", "fitI", "fitI.new")

# MCMC settings
# na <- 1000 ; ni <- 60000 ; nb <- 10000 ; nt <- 5 ; nc <- 10
na <- 1000 ; ni <- 60000 ; nb <- 10000 ; nt <- 5 ; nc <- 3  # ~~~ for testing

# Call JAGS (ART 14 min), gauge convergence and summarize posteriors
set.seed(127) # Cheating: use this seed to get good initial values.
(out9 <- jags(bdata, inits, params, "frogs.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE))
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out9, layout=c(2,2))
print(out9, dig = 2)
#               mean    sd   2.5%    50%  97.5% overlap0    f Rhat  n.eff
# alpha.lamN    1.16  0.06   1.04   1.16   1.27    FALSE 1.00 1.00  39674
# alpha.lamI    1.15  0.06   1.03   1.15   1.26    FALSE 1.00 1.00 100000
# mean.pN       0.81  0.01   0.78   0.81   0.83    FALSE 1.00 1.00  20958
# mean.pI       0.81  0.01   0.79   0.81   0.83    FALSE 1.00 1.00  34295
# mean.omegaN   0.87  0.05   0.77   0.87   0.96    FALSE 1.00 1.00   5198
# mean.omegaI   0.76  0.07   0.63   0.77   0.88    FALSE 1.00 1.01    578
# mean.gammaN   0.64  0.13   0.38   0.65   0.88    FALSE 1.00 1.00   1681
# mean.gammaI   0.16  0.18  -0.21   0.17   0.47     TRUE 0.82 1.00   1331
# mean.psi_NI   0.07  0.04   0.01   0.06   0.14    FALSE 1.00 1.00   1665
# mean.psi_IN   0.16  0.06   0.03   0.16   0.29    FALSE 1.00 1.01    829
# fitN        140.91  4.62 132.62 140.65 150.69    FALSE 1.00 1.00  12765
# fitN.new    146.90 10.90 126.58 146.58 169.22    FALSE 1.00 1.00  19525
# fitI        131.04  4.22 123.59 130.76 140.08    FALSE 1.00 1.00  16866
# fitI.new    138.13 10.47 118.57 137.80 159.49    FALSE 1.00 1.00  40418

op <- par(mfrow=c(1,2)) # Not shown
pp.check(out9, observed = 'fitN', simulated = 'fitN.new')
pp.check(out9, observed = 'fitI', simulated = 'fitI.new')
par(op)

# 2.9.2 Example 2: dusky salamanders -- a three-state model with two observable states
# ------------------------------------------------------------------------------------

data(duskySalamanders) # require(AHMbook)
str(duskySalamanders)
n <- duskySalamanders
nsites <- dim(n)[1]   # Total number of locations
nyears <- dim(n)[2]   # Total number of survey years
nObsAges <- dim(n)[3] # Number of observable stages
nsurveys <- dim(n)[4] # Number of replicate surveys

# Bundle data
str(bdata <- list(nsites = nsites, nyears = nyears, nsurveys = nsurveys, n = n))
# List of 4
# $ nsites   : int 21
# $ nyears   : int 7
# $ nsurveys : int 2
# $ n        : int [1:21, 1:7, 1:2, 1:2] 12 0 0 0 0 0 1 3 6 0 ...

# Specify model in BUGS language
cat(file = "ZipkinModel.txt", "
model {
  # Priors
  for(c in 1:3){ # Loop over 3 states
    lambda[c] ~ dunif(0, 50) # Initial abundance
    gamma[c] ~ dunif(0, 50) # Recruitment
  }
  phi[1] ~ dunif(0, 1) # Apparent survival
  phi[2] ~ dunif(0, 1)
  p[1] ~ dunif(0, 1) # Detection
  p[2] ~ dunif(0, 1)

  # Process model
  for(i in 1:nsites) { # Loop over sites
    #Initial abundance state
    N[i,1,1] ~ dpois(lambda[1]) # year 1 juveniles
    N[i,1,2] ~ dpois(lambda[2]) # year 2 juveniles
    N[i,1,3] ~ dpois(lambda[3]) # adults
    #Specify the model for years 2 through nYears
    for(t in 2:nyears) {
      # Number of survivors in each age class
      S[i,t,1] ~ dbin(phi[1], N[i,t-1,1])
      S[i,t,2] ~ dbin(phi[1], N[i,t-1,2])
      S[i,t,3] ~ dbin(phi[2], N[i,t-1,3])
      # Number of recruits (gamma[1]) and immigrants (gamma[2], gamma[3])
      G[i,t,1] ~ dpois(gamma[1]*N[i,t-1,3] + gamma[2])
      G[i,t,2] ~ dpois(gamma[2])
      G[i,t,3] ~ dpois(gamma[3])
      # Sum all stages to get total N at each site i in each year t
      N[i,t,1] <- G[i,t,1]
      N[i,t,2] <- S[i,t,1] + G[i,t,2]
      N[i,t,3] <- S[i,t,2] + S[i,t,3] + G[i,t,3]
    }

    # Observation model
    for(t in 1:nyears) {
      for(j in 1:nsurveys) {
        # Stages S1 and S2 are indistinguishable
        n[i,t,1,j] ~ dbin(p[1], (N[i,t,1]+N[i,t,2]))
        # Detection probability is the same for all adults
        n[i,t,2,j] ~ dbin(p[2], N[i,t,3])
      }
    }
  }

  # Derived quantities: Total N for each stage and year
  for (t in 1:nyears) {
    Ntotal[1,t] <- sum(N[,t,1])
    Ntotal[2,t] <- sum(N[,t,2])
    Ntotal[3,t] <- sum(N[,t,3])
  }
}
")

# Find initial values that work for JAGS ...
lamNew <- NA ; phiNew <- NA ; gammaNew <- NA ; pNew <- NA
lamNew[1] <- 5 * 10     # Initial population size (juveniles)
lamNew[2] <- 5 * 10     # Initial population size (adults)
phiNew[1] <- 0.90       # Survival rate (juveniles)
phiNew[2] <- 0.90       # Survival rate (adults)
gammaNew[1] <- 5 * 10   # Recruitment rate
gammaNew[2] <- 5 * 10   # Movement rate (juveniles)
gammaNew[3] <- 5 * 10   # Movement rate (adults)
pNew[1] <- 0.8          # Detection probability
pNew[2] <- 0.8

# Starting values for N[,1,] can be all-important !!!
N <- N1 <- array(NA, dim = c(nsites, nyears, 3) )
for(i in 1:nsites){
  for (t in 1:nyears){
    max.it1 <- max(n[i,t,1,]) + 5
    max.it2 <- max(n[i,t,2,]) + 5
    N[i,t,1] <- max.it1
    N[i,t,2] <- max.it1
    N[i,t,3] <- max.it2 + 2
  }
}
N1[,1,] <- N[,1,]

# Package all that into the Initial values function
inits <- function() list(phi = phiNew, gamma = runif(3,1,3), p = pNew, N = N1)

# Parameters monitored
params <- c("lambda", "phi", "gamma", "p", "Ntotal")

# MCMC settings
na <- 1000 ; ni <- 100000 ; nb <-50000 ; nt <- 2 ; nc <- 3

# Call JAGS (ART 7 min), check convergence and summarize posteriors
out10 <- jags(bdata, inits, params, "ZipkinModel.txt", n.adapt = na, n.thin = nt,
    n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out10)
print(out10, dig = 3)
#               mean     sd   2.5%    50%   97.5% overlap0 f  Rhat n.eff
# lambda[1]    1.961  1.545  0.098  1.629   5.908    FALSE 1 1.009   683
# lambda[2]    4.804  1.593  2.360  4.571   8.596    FALSE 1 1.003  1343
# lambda[3]    5.546  0.836  4.037  5.500   7.295    FALSE 1 1.004   874
# phi[1]       0.480  0.100  0.294  0.476   0.685    FALSE 1 1.003  1305
# phi[2]       0.701  0.067  0.559  0.706   0.822    FALSE 1 1.008   491
# gamma[1]     0.786  0.335  0.355  0.713   1.634    FALSE 1 1.006  1070
# gamma[2]     0.026  0.027  0.001  0.018   0.099    FALSE 1 1.000 17349
# gamma[3]     0.054  0.053  0.001  0.038   0.195    FALSE 1 1.002  4467
# p[1]         0.183  0.054  0.092  0.178   0.300    FALSE 1 1.003  1210
# p[2]         0.369  0.039  0.299  0.367   0.452    FALSE 1 1.008   464
# Ntotal[1,1] 40.160 31.792  2.000 33.000 121.000    FALSE 1 1.009   658
# Ntotal[2,1] 99.831 31.896 52.000 95.000 177.000    FALSE 1 1.003  1329
# ... output truncated ...
