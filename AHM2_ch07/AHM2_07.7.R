#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7 : MODELING FALSE POSITIVES
# ====================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 6 mins

library(jagsUI)
library(AHMbook)

# 7.7 False positives models in open populations
# ==============================================

# Simulation settings
set.seed(129)
nsites <- 200            # number of sites
nsurv1 <- 3              # number of occasions with Type 1 data
nsurv2 <- 4              # number of occasions with Type 2 data
Tsurv <- nsurv1 + nsurv2 # Total number of surveys
nyears <- 5              # number of years
psi <- 0.6               # initial expected proportion of sites occupied
gamma <- 0.1             # colonization rate
phi <- 0.88              # survival rate
p <- c(0.7, 0.5)         # detection prob of method 1 and method 2
fp <- 0.05               # false-positive error probability (p_10)
b <- 0.2                 # probability y is recorded as certain

# Simulate latent state matrix (now 2D) and the observations (now 3D)
z <- matrix(NA, nrow = nsites, ncol = nyears)
z[, 1] <- rbinom(nsites, 1, psi)
y <- array(NA, dim = c(nsites, Tsurv, nyears) )
for(t in 1:nyears){
  if(t > 1)
    z[,t] <- rbinom(nsites, 1, phi*z[,t-1] + gamma*(1-z[,t-1]) )
  for(i in 1:nsites){
    p1 <- p[1]*z[i,t]                     # certain detection (method 1)
    p2 <- p[2]*z[i,t] + fp*(1-z[i,t])     # uncertain detection (method 2)
    y[i, 1:3, t] <- rbinom(nsurv1, 1, p1) # simulate method 1 data
    y[i, 4:7, t] <- rbinom(nsurv2, 1, p2) # simulate method 2 data
    # now introduce certain observations: (type 3 data in unmarked)
    pr.certain <- z[i, t] * y[i, 4:7, t] * b
    y[i, 4:7, t] <- y[i, 4:7, t] + rbinom(4, 1, pr.certain)
  }
}

# Convert to categorical data so non-detection = 1, certain detection = 3
y[ , (nsurv1 + 1):(Tsurv), 1:nyears] <- 1 + y[ ,(nsurv1 + 1):(Tsurv), 1:nyears]

# Bundle data and summarize data bundle
str( bdata <- list(y = y, nsites = nrow(y), nsurv1 = nsurv1, nsurv2 = nsurv2,
    nyears = nyears) )
# List of 5
# $ y        : num [1:200, 1:7, 1:5] 0 1 1 0 1 0 1 0 0 1 ...
# $ nsites   : int 200
# $ nsurveys1: num 3
# $ nsurveys2: num 4
# $ nyears   : num 5

# Specify model in BUGS language
cat(file = "occufp3.txt","
model {

  # Priors
  psi ~ dunif(0, 1)
  fp ~ dunif(0, 1)
  b ~ dunif(0, 1)
  gamma ~ dunif(0, 1)
  phi ~ dunif(0, 1)
  alpha0 ~ dnorm(0,0.01)    # Logit-linear parameters for detection
  alpha1 ~ dnorm(0,0.01)    # Method effect

  # Likelihood, t = 1
  for (i in 1:nsites) {
    z[i,1] ~ dbern(psi)                          # State model
    for(j in 1:(nsurv1 + nsurv2)){
      probs[i,j,1,1,1] <- 1-fp                   # z = 0 obs probs
      probs[i,j,2,1,1] <- fp
      probs[i,j,3,1,1] <- 0
      probs[i,j,1,1,2] <- 1 - p[i,j,1]           # z = 1 obs probs
      probs[i,j,2,1,2] <- (1-b)*p[i,j,1]
      probs[i,j,3,1,2] <- p[i,j,1]*b
    }
    for(j in 1:nsurv1) {                         # Loop over replicate surveys
      logit(p[i,j,1]) <- alpha0
      y[i,j,1] ~ dbern(z[i,1]*p[i,j,1] )         # Type 1 data
    }
    for (j in (nsurv1+1):(nsurv1+nsurv2) ) {     # Replicate surveys
      logit(p[i,j,1]) <- alpha0 + alpha1
      y[i,j,1] ~ dcat(probs[i,j,1:3,1, z[i,1]+1 ] ) # Type 3 data
    }

    # Model components for t = 2, ...
    for(t in 2:nyears) {
      z[i,t] ~ dbern(z[i,t-1]*phi + (1-z[i,t-1])*gamma ) # State model t>1
      for(j in 1:(nsurv1+nsurv2)) {
        probs[i,j,1,t,1] <- 1-fp                 # z = 0 obs probs
        probs[i,j,2,t,1] <- fp
        probs[i,j,3,t,1] <- 0
        probs[i,j,1,t,2] <- 1-p[i,j,t]           # z = 1 obs probs
        probs[i,j,2,t,2] <- (1-b)*p[i,j,t]
        probs[i,j,3,t,2] <- p[i,j,t]*b
      }
      for(j in 1:nsurv1) {                       # Loop over replicate surveys
        logit(p[i,j,t]) <- alpha0
        y[i,j,t] ~ dbern(z[i,t]*p[i,j,t] )       # Type 1 data
      }
      for (j in (nsurv1+1):(nsurv1+nsurv2) ) {   # Loop over surveys
        logit(p[i,j,t]) <- alpha0 + alpha1
        y[i,j,t] ~ dcat(probs[i,j,1:3,t, z[i,t]+1 ] ) # Type 3 data
      }
    }
  }
}
")

# Initial values.
zst <- apply(y[, 1:(nsurv1+nsurv2), ], c(1, 3), max)
zst[zst>1] <- 1
inits <- function(){ list(z = zst, alpha0 = rnorm(1, -1, 1),
    alpha1 = rnorm(1, 0, 1), phi = 0.7, gamma = 0.1, fp = 0.05, b = 0.1)}

# Parameters monitored
params <- c("psi", "alpha0", "alpha1", "phi", "gamma", "fp", "b")

# MCMC settings
na <- 1000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3

# Call JAGS (ART 4 min), assess convergence and summarize posteriors
out5 <- jags(bdata, inits, params, "occufp3.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out5, layout=c(2,3))
print(out5, dig = 3)
#          mean    sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
# psi     0.634 0.034  0.567  0.635  0.697    FALSE 1 1.001  1081
# alpha0  0.860 0.054  0.756  0.860  0.969    FALSE 1 1.013   156
# alpha1 -0.892 0.068 -1.027 -0.892 -0.762    FALSE 1 1.012   177
# gamma   0.088 0.015  0.061  0.088  0.118    FALSE 1 1.002  1563
# phi     0.839 0.019  0.802  0.839  0.874    FALSE 1 1.001  1627
# fp      0.048 0.005  0.039  0.048  0.059    FALSE 1 1.000  3000
# b       0.210 0.012  0.186  0.209  0.235    FALSE 1 1.000  3000


# With water vole data
# ''''''''''''''''''''
data(waterVoles)   # This is from AHMbook
(wv <- waterVoles) # not shown
usites <- levels(wv$Patch)
nsites <- nlevels(wv$Patch)
yr <- as.numeric(as.character(wv$Year))

# Make two dummy occasions for Type = 1 and Type = 3, all NA for the
# water vole data
y <- array(NA, dim = c(nsites, 5, 3))
m <- 1
for(i in c(2009, 2010, 2011)){
  x <- wv[yr==i,]
  y[match(x[,1], usites), 2:4, m] <- as.matrix(x[,2:4])
  m <- m+1
}
nsurv1 <- 1       # NA data Type 1
nsurv2 <- 3       # real data Type 2 data
nsurv3 <- 1       # NA data Type 3
Tsurv <- nsurv1 + nsurv2 + nsurv3
nyears <- 3

# Convert data into a categorical form 1/2/3 instead of 0/1/2
y[, (nsurv1 + nsurv2 + 1):Tsurv, 1:nyears] <- 1 +
    y[, (nsurv1 + nsurv2 + 1):Tsurv, 1:nyears]

# Bundle data and summarize data bundle
str( bdata <- list(y = y, nsites = nrow(y), nsurv1 = nsurv1, nsurv2 = nsurv2,
    nsurv3 = nsurv3, nyears = nyears, Tsurv = Tsurv) )

# Specify model in BUGS language
cat(file = "occufp4.txt","
model {

  # Priors
  psi ~ dunif(0, 1)
  b ~ dunif(0,1)
  for(t in 1:nyears) {
    Nocc[t] <- sum(z[,t])
    fp[t] ~ dunif(0, 1)
    p[t] ~ dunif(0,1)
  }
  for(t in 1:(nyears-1)) {
    gamma[t] ~ dunif(0, 1)
    phi[t] ~ dunif(0, 1)
  }

  # Likelihood
  # Model for t = 1
  for (i in 1:nsites) {
    z[i,1] ~ dbern(psi)
    for(j in 1:nsurv1) { # Type 1 data
      y[i,j,1] ~ dbern(z[i,1]*p[1] )
    }
    for(j in (nsurv1+1):(nsurv1+nsurv2) ) { # Type 2 data
      y[i,j,1] ~ dbern(z[i,1]*p[1] + fp[1]*(1-z[i,1]) )
    }
    for(j in (nsurv1+nsurv2+1):Tsurv ) { # Type 3 data
      probs[i,j,1,1,1] <- 1-fp[1] # z = 0 obs probs
      probs[i,j,2,1,1] <- fp[1]
      probs[i,j,3,1,1] <- 0
      probs[i,j,1,1,2] <- 1-p[1] # z = 1 obs probs
      probs[i,j,2,1,2] <- (1-b)*p[1]
      probs[i,j,3,1,2] <- p[1]*b # Observation model
      y[i,j,1] ~ dcat(probs[i,j,1:3,1, z[i,1]+1 ] )
    }

    # Model for t>1
    for(t in 2:nyears) {
      z[i,t] ~ dbern(z[i,t-1]*phi[t-1] + (1-z[i,t-1])*gamma[t-1])
      for(j in 1:nsurv1) { # Type 1 data
      y[i,j,t] ~ dbern(z[i,t]*p[t] )
      }
        for(j in (nsurv1+1):(nsurv1+nsurv2 )) { # Type 2 data
        y[i,j,t] ~ dbern(z[i,t]*p[t] + fp[t]*(1-z[i,t]) )
      }
      for(j in (nsurv1+nsurv2+1):Tsurv ) { # Type 3 data
        probs[i,j,1,t,1] <- 1-fp[t] # z = 0 obs probs
        probs[i,j,2,t,1] <- fp[t]
        probs[i,j,3,t,1] <- 0
        probs[i,j,1,t,2] <- 1-p[t] # z = 1 obs probs
        probs[i,j,2,t,2] <- (1-b)*p[t]
        probs[i,j,3,t,2] <- p[t]*b # Observation model
        y[i,j,t] ~ dcat(probs[i,j,1:3,t, z[i,t]+1 ] )
      }
    }
  }
}
")

# Initial values
zst <- apply(y[,,], c(1,3), max, na.rm = TRUE)
zst[zst>1] <- 1
zst[zst == -Inf] <- 0
inits <- function(){ list(z = zst, p = c(0.8, 0.8, 0.8), phi = c(0.5, 0.9),
    gamma = c(0.1, 0.1), fp = c(0.05, 0.05, 0.05), b = 0.1 ) }

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "fp", "b", "Nocc")

# MCMC settings
na <- 1000 ; ni <- 5000; nt <- 1 ; nb <- 1000 ; nc <- 3

# Call JAGS (ART < 1 min), assess convergence and summarize posteriors
out6 <- jags(bdata, inits, params, "occufp4.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(out6, layout=c(2,3))
print(out6, 3)
#            mean    sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
# psi       0.346 0.058  0.240  0.343  0.464    FALSE 1 1.003   710
# p[1]      0.703 0.060  0.581  0.705  0.817    FALSE 1 1.003   650
# p[2]      0.898 0.038  0.813  0.901  0.963    FALSE 1 1.001  4481
# p[3]      0.938 0.033  0.865  0.941  0.991    FALSE 1 1.001  1466
# gamma[1]  0.234 0.063  0.118  0.231  0.365    FALSE 1 1.001  8185
# gamma[2]  0.274 0.085  0.122  0.270  0.454    FALSE 1 1.001 12000
# phi[1]    0.701 0.089  0.520  0.704  0.866    FALSE 1 1.001  1813
# phi[2]    0.802 0.084  0.622  0.809  0.942    FALSE 1 1.000 12000
# fp[1]     0.038 0.020  0.006  0.036  0.081    FALSE 1 1.003   702
# fp[2]     0.139 0.032  0.081  0.138  0.206    FALSE 1 1.001  2970
# fp[3]     0.347 0.053  0.241  0.348  0.452    FALSE 1 1.000 12000
# b         0.501 0.289  0.025  0.502  0.979    FALSE 1 1.000  4071
# Nocc[1]  39.154 4.278 32.000 39.000 48.000    FALSE 1 1.004   497
# Nocc[2]  44.627 3.561 38.000 45.000 52.000    FALSE 1 1.002  2163
# Nocc[3]  54.910 6.244 43.000 55.000 68.000    FALSE 1 1.001 12000

# ~~~~ code for figure 7.6 ~~~~~~~~~~
op <- par(mfrow = c(1,2), mar = c(5,4,4,3))
Nocc<- out6$sims.list$Nocc
colnames(Nocc)<- c("2009","2010","2011")
boxplot(Nocc, ylab="Number of occupied sites",xlab="Year",
    outline = FALSE, frame = FALSE)
pmat<- out6$sims.list$p
colnames(pmat)<- c("2009","2010","2011")
fpmat<- out6$sims.list$fp
colnames(fpmat)<- colnames(pmat)
boxplot(pmat, ylim=c(0,1), col="cyan", outline = FALSE, frame = FALSE,
    xlab="Year")
boxplot(fpmat,add=TRUE, col="red", outline = FALSE, frame = FALSE)
legend(2.5, 0.16, legend=c("detection, p","false positive, fp"),
    pch=20, col=c("cyan","red"))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
