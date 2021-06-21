# Applied hierarchical modeling in ecology - vol.2 - 2021
# Marc Kéry & J. Andy Royle
#
# Chapter 11 : SPATIALLY EXPLICIT DISTANCE SAMPLING ALONG TRANSECTS
# =================================================================
# Code from proofs dated 2020-08-19

print("Approximate execution time for this code: 40 mins")
# Run time with the full number of iterations: 60 mins

library(AHMbook)
library(jagsUI)

# 11.9 Analysis of wiggly transects
# =================================

library(raster)
library(AHMbook)
data(wigglyLine)

RNGversion("3.5.0")
set.seed(123, kind = "Mersenne-Twister") # This is important!
plot(wigglyLine[,1], wigglyLine[,2], type = "l", xlab = "Easting",
    ylab = "Northing", pch = " ", lwd = 2, cex.axis = 1.5, cex = 2,
    cex.lab = 1.5, asp = 1, frame = FALSE)
points <- SpatialPoints(wigglyLine)
sLine <- Line(points)
regpoints <- spsample(sLine, 100, type = "regular")
points(regpoints, col = "black", pch = 20, lwd = 2)  # produces Fig. 11.11

# ~~~ extra code for figure 11.13 ~~~~~~~~~~~~~~~~
# Set up grid of pixels
xlim <- c(-1, 4)
ylim <- c(-1, 5)
U <- expand.grid(seq(xlim[1],xlim[2],0.1), seq(ylim[1],ylim[2],0.1))
# Calculate total hazard of detection for each pixel
sigma <- 0.4
alpha0 <- -2
pmat <- rep(NA,nrow(U))
for (i in 1:nrow(U)) {
  dvec <- sqrt((U[i, 1] - wigglyLine[, 1])^2 + (U[i,  2] - wigglyLine[, 2])^2)
  loghaz <- alpha0 - (1/(2*sigma*sigma)) * dvec * dvec
  H <- sum(exp(loghaz))
  pmat[i] <- 1 - exp(-H)
}
# Do the plot
library(raster)
r <- rasterFromXYZ(cbind(U, pmat))
plot(r,  frame=FALSE, box=FALSE, col=topo.colors(20))
lines(wigglyLine, lwd=2)
# Add the 'equidetectable' points
points(c(2.508281, 3.211953), c(1.946258, 1.959366), pch=20)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 11.9.1 Data simulation
# ----------------------

set.seed(2014, kind = "Mersenne-Twister")
# We first use a uniform distribution here so that individuals are not
# distributed in relation to landscape structure
xlim = c(-1.0, 4.0)
ylim = c(-1.0, 5.0)
N <- 120

# Central locations or activity centers
sx <- runif(N, xlim[1], xlim[2])
sy <- runif(N, ylim[1], ylim[2])
plot(cbind(sx, sy), pch = " ", xlim = xlim, ylim = ylim, asp = 1, xlab = " ",
    ylab = " ", frame = FALSE)                # begin Fig. 11.12
lines(regpoints@coords, col = "black", lwd = 4)
points(sx, sy, pch = 3, col = "black", lwd = 1)
sigma.move <- 0.35    # This is either ind's moving around or measurement error
sigma <- 0.4
alpha0 <- 0.8
# Line chopped up into points -- "traps"
X <- regpoints@coords
J <- nrow(X)

# We sample this population K = 2 times along the same line
# (need not be the same line though)
nsurveys <- 2
U <- array(NA, dim = c(N, nsurveys, 2))
y <- pmat <- matrix(NA, nrow = N, ncol = nsurveys)

for (i in 1:N) {
  for (k in 1:nsurveys) {
    # Instaneous locations or movement outcomes
    U[i, k, ] <- c(rnorm(1, sx[i], sigma.move), rnorm(1, sy[i], sigma.move))
    dvec <- sqrt((U[i, k, 1] - X[, 1])^2 + (U[i, k, 2] - X[, 2])^2)
    loghaz <- alpha0 - (1/(2*sigma*sigma)) * dvec *dvec
    # Total hazard model here
    H <- sum(exp(loghaz))
    pmat[i, k] <- 1 - exp(-H)
    # Alternative: closest distance to the line
    # pmat[i,k] <- plogis(alpha0)*exp(- (min(dvec)^2)/(2*sigma^2) )
    y[i, k] <- rbinom(1, 1, pmat[i, k])
  }
}
# For data organization we separate the x- and y- coordinates
Ux <- U[, , 1]
Uy <- U[, , 2]
points(Ux, Uy, pch = 1, col = "black")        # Add them to plot Fig. 11.12
# If the individual was not detected, its location is missing
Ux[y == 0] <- NA
Uy[y == 0] <- NA
points(Ux, Uy, pch = 20, col = "black")       # Add detections to plot Fig. 11.12
# Retain individuals detected at least once
ncap <- apply(y, 1, sum)
y <- y[ncap > 0, ]
Ux <- Ux[ncap > 0, ]
Uy <- Uy[ncap > 0, ]

# 11.9.2 Running the models in JAGS
# ---------------------------------

# Data augmentation
M <- 400
nind <- nrow(y)
y <- rbind(y, matrix(0, nrow = (M - nrow(y)), ncol = ncol(y)))

# Ux and Uy get augmented with missing values. S inits as just described
Namat <- matrix(NA, nrow = (M - nind), ncol = ncol(y))
Ux <- rbind(Ux, Namat)
Uy <- rbind(Uy, Namat)
S <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
for (i in 1:nind) {
  S[i, ] <- c(mean(Ux[i, ], na.rm = TRUE), mean(Uy[i, ], na.rm = TRUE))
}

# For missing locations, use S[i, ] or close-by points to avoid inconsistent
# node crash
Ux.st <- Ux
Uy.st <- Uy
for (i in 1:M) {
  Ux.st[i, !is.na(Ux[i, ])] <- NA
  Uy.st[i, !is.na(Uy[i, ])] <- NA
  Ux.st[i, is.na(Ux[i, ])] <- S[i, 1]
  Uy.st[i, is.na(Uy[i, ])] <- S[i, 2]
}

# Bundle the data
str(bdata <- list(y = y, u = Ux, v = Uy, X = X, nsurveys = nsurveys, M = M,
    J = J, xlim = xlim, ylim = ylim) )
# List of 9
# $ y       : num [1:400, 1:2] 0 1 1 0 1 1 1 1 1 1 ...
# $ u       : num [1:400, 1:2] NA 2.3487 -0.0551 NA 1.525 ...
# $ v       : num [1:400, 1:2] NA 3.29 2.75 NA 1.99 ...
# $ X       : num [1:100, 1:2] 0.155 0.209 0.263 0.317 0.369 ...
# $ nsurveys: num 2
# $ M       : num 400
# $ J       : int 100
# $ xlim    : num [1:2] -1 4
# $ ylim    : num [1:2] -1 5

# Define in BUGS the SCR model with total hazard detection probability
cat(file = "wiggly1.txt","
model {

  # Priors
  alpha0 ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 10)
  sigma.move ~ dunif(0, 10)
  tau <- 1/(sigma.move*sigma.move)
  psi ~ dunif(0, 1)

  # Likelihood
  for(i in 1:M){ # Loop over individuals in augmented list
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    for(k in 1:nsurveys){ # Loop over temporal replicates
      u[i,k] ~ dnorm(s[i,1],tau)
      v[i,k] ~ dnorm(s[i,2],tau)
      for(j in 1:J){ # Loop over each point defining line segments
      d[i,k,j] <- pow(pow(u[i,k]-X[j,1],2) + pow(v[i,k]-X[j,2],2),0.5)
      h[i,k,j] <- exp(alpha0 - (1/(2*sigma*sigma))*d[i,k,j]*d[i,k,j])
      }
      H[i,k] <- sum(h[i,k,1:J])
      p[i,k] <- z[i]*(1-exp(-H[i,k]))
      y[i,k] ~ dbern(p[i,k])
    }
  }

  # Population size is a derived quantity
  N <- sum(z[])
}
")

# Define inits and parameters to monitor
inits <- function() { list(alpha0 = rnorm(1,-1,0.5), sigma = 0.5,
    sigma.move = 0.5, s = S, z = c(rep(1, nind), rep(0, M - nind)),
    u = Ux.st, v = Uy.st) }
params <- c("alpha0", "N", "psi", "sigma.move", "sigma")

# MCMC settings
# na <- 500 ; nt <- 1 ; nc <- 5 ; nb <- 500 ; ni <- 3000
na <- 500 ; nt <- 1 ; nc <- 3 ; nb <- 50 ; ni <- 300  # ~~~ for testing, 10 mins

# Run JAGS (ART 74 min), assess convergence and summarize posteriors
wiggly.Thaz.2reps <- jags(bdata, inits, params, "wiggly1.txt", n.adapt = na,
    n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(wiggly.Thaz.2reps, layout=c(2,3))
print(wiggly.Thaz.2reps, 3)
#              mean   sd  2.5%    50%  97.5% overlap0    f Rhat n.eff
# alpha0       0.87 1.12 -0.88   0.71   3.41     TRUE 0.77 1.04    83
# N          113.09 9.91 96.00 112.00 134.00    FALSE 1.00 1.01   343
# psi          0.28 0.03  0.22   0.28   0.35    FALSE 1.00 1.00   616
# sigma.move   0.35 0.03  0.31   0.35   0.41    FALSE 1.00 1.01   379
# sigma        0.40 0.05  0.32   0.40   0.51    FALSE 1.00 1.04    88


# Define in BUGS the SCR model with minimum distance detection model
cat(file = "wiggly2.txt","
model {

  # Priors
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, 10)
  sigma.move ~ dunif(0, 10)
  tau <- 1/(sigma.move*sigma.move)
  psi ~ dunif(0, 1)

  # Likelihood
  for(i in 1:M){ # Loop over individuals
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    for(k in 1:nsurveys){ # Loop over temporal replicates
      u[i,k] ~ dnorm(s[i,1], tau)
      v[i,k] ~ dnorm(s[i,2], tau)
      for(j in 1:J){ # Loop over line segments
        d[i,k,j] <- pow(pow(u[i,k]-X[j,1],2) + pow(v[i,k]-X[j,2],2),0.5)
      }
      # Compute closest point on the line and specify detection model
      closestD[i,k] <- min(d[i,k,1:J])
      p[i,k] <- z[i]*p0*exp(-(closestD[i,k]^2 )/ (2*sigma*sigma) )
      y[i,k] ~ dbern(p[i,k])
    }
  }

  # Derived quantity
  N <- sum(z[])
}
")

# Data, inits and parameters to save
str(bdata <- list(y = y, u = Ux, v = Uy, X = X, nsurveys = nsurveys, M = M,
    J = J, xlim = xlim, ylim = ylim) ) # Output omitted
inits <- function() { list(p0 = runif(1,0,1), sigma = 0.5, sigma.move = 0.5,
    s = S, z = c(rep(1, nind), rep(0, M - nind)), u = Ux.st, v = Uy.st)}
params <- c("p0", "N", "psi", "sigma.move", "sigma")

# MCMC settings
na <- 500 ; nt <- 1 ; nc <- 5 ; nb <- 500 ; ni <- 3500

# Run JAGS (ART 11 min), assess convergence and summarize posteriors
wiggly.closeD.2reps <- jags(bdata, inits, params, "wiggly2.txt", n.adapt = na,
    n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(wiggly.closeD.2reps, layout=c(2,3))
print(wiggly.closeD.2reps, 3)
#              mean    sd  2.5%    50%  97.5% overlap0 f Rhat n.eff
# p0           0.98  0.02  0.93   0.99   1.00    FALSE 1 1.00  7379
# N          115.86 11.88 95.00 115.00 141.00    FALSE 1 1.00  1919
# psi          0.29  0.04  0.22   0.29   0.37    FALSE 1 1.00  2768
# sigma.move   0.35  0.03  0.30   0.34   0.40    FALSE 1 1.02   182
# sigma        0.73  0.06  0.62   0.73   0.87    FALSE 1 1.00  5878

# ~~~~ rerun the data generation code with sigma.move = 0 ~~~~~~~
set.seed(2014, kind = "Mersenne-Twister")
# We first use a uniform distribution here so that individuals are not
# distributed in relation to landscape structure
xlim = c(-1.0, 4.0)
ylim = c(-1.0, 5.0)
N <- 120
# Central locations or activity centers
sx <- runif(N, xlim[1], xlim[2])
sy <- runif(N, ylim[1], ylim[2])

sigma.move <- 0                 # <---- changed this
alpha0 <- 0.8
# Line chopped up into points -- "traps"
X <- regpoints@coords
J <- nrow(X)
# We sample this population K = 2 times along the same line
# (need not be the same line though)
nsurveys <- 2
U <- array(NA, dim = c(N, nsurveys, 2))
y <- pmat <- matrix(NA, nrow = N, ncol = nsurveys)
for (i in 1:N) {
  for (k in 1:nsurveys) {
    # Instaneous locations or movement outcomes
    U[i, k, ] <- c(rnorm(1, sx[i], sigma.move), rnorm(1, sy[i], sigma.move))
    dvec <- sqrt((U[i, k, 1] - X[, 1])^2 + (U[i, k, 2] - X[, 2])^2)
    loghaz <- alpha0 - (1/(2*sigma*sigma)) * dvec *dvec
    # Total hazard model here
    H <- sum(exp(loghaz))
    pmat[i, k] <- 1 - exp(-H)
    # Alternative: closest distance to the line
    # pmat[i,k] <- plogis(alpha0)*exp(- (min(dvec)^2)/(2*sigma^2) )
    y[i, k] <- rbinom(1, 1, pmat[i, k])
  }
}
# For data organization we separate the x- and y- coordinates
Ux <- U[, , 1]
Uy <- U[, , 2]
# If the individual was not detected, its location is missing
Ux[y == 0] <- NA
Uy[y == 0] <- NA
# Retain individuals detected at least once
ncap <- apply(y, 1, sum)
y <- y[ncap > 0, ]
Ux <- Ux[ncap > 0, ]
Uy <- Uy[ncap > 0, ]

# Data augmentation
M <- 400
nind <- nrow(y)
y <- rbind(y, matrix(0, nrow = (M - nrow(y)), ncol = ncol(y)))
# Ux and Uy get augmented with missing values. S inits are random uniform or mean observed location
Namat <- matrix(NA, nrow = (M - nind), ncol = ncol(y))
Ux <- rbind(Ux, Namat)
Uy <- rbind(Uy, Namat)
S <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
for (i in 1:nind) {
  S[i, ] <- c(mean(Ux[i, ], na.rm = TRUE), mean(Uy[i, ], na.rm = TRUE))
}
#
# For missing locations, use S[i, ] or close-by points to avoid inconsistent node
Ux.st <- Ux
Uy.st <- Uy
for (i in 1:M) {
Ux.st[i, !is.na(Ux[i, ])] <- NA
Uy.st[i, !is.na(Uy[i, ])] <- NA
Ux.st[i, is.na(Ux[i, ])] <- S[i, 1]
Uy.st[i, is.na(Uy[i, ])] <- S[i, 2]
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Define in BUGS the MRDS model using minimum distance to the line
cat(file = "wiggly3.txt","
model {

  # Priors
  sigma ~ dunif(0, 10)
  p0 ~ dunif(0, 1)
  psi ~ dunif(0, 1)

  # Likelihood
  for(i in 1:M){                          # Loop over individuals
    z[i] ~ dbern(psi)
    # s = observed location here
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    for(k in 1:nsurveys){ # Loop over temporal replicates
      for(j in 1:J){ # Loop over each line segment
        d[i,k,j] <- pow(pow(s[i,1]-X[j,1], 2) + pow(s[i,2]-X[j,2], 2), 0.5)
      }
      # Compute closest point on the line and define detection model
      closestD[i,k] <- min(d[i,k,1:J])
      p[i,k] <- z[i]*p0*exp(-(closestD[i,k]^2 )/ (2*sigma*sigma) )
      y[i,k] ~ dbern(p[i,k])
    }
  }

  # Derived quantity
  N <- sum(z[])
}
")

# Compute average observed location
sbar <- cbind(apply(Ux, 1, mean, na.rm = TRUE), apply(Uy, 1, mean,
    na.rm = TRUE))
sbar[is.nan(sbar)] <- NA

# Data, inits and parameters
str(bdata <- list(y = y, s = sbar, X = X, nsurveys = nsurveys, M = M, J = J,
    xlim = xlim, ylim = ylim)) # Output not shown
inits <- function() { list( p0 = runif(1, 0, 1), sigma = 0.5,
    z = c(rep(1, nind), rep(0, M - nind)) )}
params <- c("p0", "N", "psi", "sigma")

# MCMC settings
na <- 500 ; nt <- 1 ; nc <- 5 ; nb <- 500 ; ni <- 1500

# Run JAGS (ART 2 min), assess convergence and summarize posteriors
wiggly.MRDS <- jags(bdata, inits, params, "wiggly3.txt", n.adapt = na,
    n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(wiggly.MRDS, layout=c(2,3))
print(wiggly.MRDS, 3)
#         mean    sd  2.5%    50%  97.5% overlap0 f Rhat n.eff
# p0      0.99  0.01  0.95   0.99   1.00    FALSE 1 1.01   753
# N     109.12 11.06 90.00 108.00 133.00    FALSE 1 1.02   238
# psi     0.27  0.04  0.21   0.27   0.35    FALSE 1 1.01   314
# sigma   0.85  0.09  0.71   0.84   1.04    FALSE 1 1.01   336


# ~~~~~~~~ extra code for the last model ~~~~~~~~~~~~

# Define in BUGS the single replicate model
cat(file = "wiggly4.txt","
model {

  # Priors
  alpha0 ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 10)
  psi ~ dunif(0,1)

  # Likelihood
  for(i in 1:M){ # Loop over individuals
    z[i] ~ dbern(psi)
    for(k in 1:nsurveys){ # Loop over temporal replicates
      u[i,k] ~ dunif(xlim[1], xlim[2])
      v[i,k] ~ dunif(ylim[1], ylim[2])
      for(j in 1:J){ # Loop over each point defining line segments
        d[i,k,j] <- pow(pow(u[i,k]-X[j,1],2) + pow(v[i,k]-X[j,2],2),0.5)
        h[i,k,j] <- exp(alpha0- (1/(2*sigma*sigma))*d[i,k,j]*d[i,k,j])
      }
      # Compute distance, hazard and detection probability
      closestD[i,k] <- min(d[i,k,1:J])
      H[i,k] <- sum(h[i,k,1:J])
      p[i,k] <- z[i]*(1-exp(-H[i,k])) # Total hazard
      # p[i,k] <- z[i]*exp( -( closestD[i,k]^2)/(2*sigma*sigma) )
      y[i,k] ~ dbern(p[i,k])
    }
  }

  # Derived quantity
  N <- sum(z[])
}
")

# Data, inits and parameters to save
str(bdata <- list(y = matrix(y[,1], ncol=1), u = matrix(Ux[,1], ncol=1),
    v = matrix(Uy[,1], ncol=1), X = X, nsurveys = 1, M = M, J = J,
    xlim = xlim, ylim = ylim) )
inits <- function() { list(sigma = 0.5, z = c(rep(1, nind),
    rep(0,M - nind)), u = Ux.st[,1,drop=FALSE], v = Uy.st[,1,drop=FALSE]) }
params <- c("alpha0", "N", "psi", "sigma")

# MCMC settings
na <- 500 ; nt <- 1 ; nc <- 5 ; nb <- 500 ; ni <- 3000  # 20 mins

# Run JAGS (ART 3 min), assess convergence and summarize posteriors
wiggly.Thaz.1rep <- jags(bdata, inits, params, "wiggly4.txt",
    n.adapt = na, n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni,
    parallel = TRUE)
# par(mfrow = c(2, 3))  # ~~~ replaced with 'layout' argument
traceplot(wiggly.Thaz.1rep, layout=c(2,3))
print(wiggly.Thaz.1rep, 2)
#          mean    sd   2.5%    50%  97.5% overlap0    f Rhat n.eff
# alpha0  -0.61  1.00  -2.47  -0.65   1.49     TRUE 0.74 1.01   490
# N      135.88 19.15 106.00 133.00 181.00    FALSE 1.00 1.01   740
# psi      0.34  0.05   0.25   0.34   0.46    FALSE 1.00 1.01  1068
# sigma    0.50  0.08   0.38   0.49   0.69    FALSE 1.00 1.01   555
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
