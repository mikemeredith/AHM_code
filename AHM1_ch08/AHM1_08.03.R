#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8. Modeling abundance using hierarchical distance sampling (HDS)
# =========================================================================

library(R2WinBUGS)
bd <- "C:/WinBUGS14"    # May have to adapt
library(AHMbook)

# 8.3 Bayesian conventional distance sampling
# ===========================================

# 8.3.1 Bayesian analysis of line transect data
# ------------------------------------------------------------------------
# Get data and do data-augmentation
# Observed distances (meters)
x <- c(71.93, 26.05, 58.47, 92.35, 163.83, 84.52, 163.83, 157.33,
    22.27, 72.11, 86.99, 50.8, 0, 73.14, 0, 128.56, 163.83, 71.85,
    30.47, 71.07, 150.96, 68.83, 90, 64.98, 165.69, 38.01, 378.21,
    78.15, 42.13, 0, 400, 175.39, 30.47, 35.07, 86.04, 31.69, 200,
    271.89, 26.05, 76.6, 41.04, 200, 86.04, 0, 93.97, 55.13, 10.46,
    84.52, 0, 77.65, 0, 96.42, 0, 64.28, 187.94, 0, 160.7, 150.45,
    63.6, 193.19, 106.07, 114.91, 143.39, 128.56, 245.75, 123.13,
    123.13, 153.21, 143.39, 34.2, 96.42, 259.81, 8.72)

B <- 500 # Strip half-width. Larger than max distance
nind <- length(x)

# Analysis of continuous data using data augmentation (DA)
nz <- 200 # Augment observed data with nz = 200 zeroes
y <- c(rep(1, nind), rep(0, nz)) # Augmented inds. have y=0 by definition
x <- c(x, rep(NA, nz)) # Value of distance are missing for the augmented

# Bundle and summarize data set
str( win.data <- list(nind=nind, nz=nz, x=x, y=y, B=B) )

# Save text file with BUGS model
cat("
model {

  # Priors
  sigma ~ dunif(0,1000)  # Half-normal scale
  psi ~ dunif(0,1)       # DA parameter

  # Likelihood
  for(i in 1:(nind+nz)){
    # Process model
    z[i] ~ dbern(psi)   # DA variables
    x[i] ~ dunif(0, B)  # Distribution of distances
    # Observation model
    logp[i] <- -((x[i]*x[i])/(2*sigma*sigma)) # Half-normal detection fct.
    p[i] <- exp(logp[i])
    mu[i] <- z[i] * p[i]
    y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
  }
  # Derived quantities
  N <- sum(z[1:(nind + nz)]) # Population size
  D <- N / 60                # Density, with A = 60 km^2 when B = 500
}
",fill=TRUE,file="model1.txt")

# Inits
zst <- y
inits <- function(){ list (psi=runif(1), z=zst, sigma=runif(1,40,200)) }

# Params to save
params <- c("N", "sigma", "D")

# Experience the raw power of BUGS and summarize marginal posteriors
library(R2WinBUGS)
# bd <- "c:/Program Files/WinBUGS14/"    # May have to adapt
out1 <- bugs(win.data, inits, params, "model1.txt",
  n.thin=2,n.chains=3, n.burnin=1000, n.iter=11000,
  # debug=TRUE, DIC=FALSE, bugs.dir=bd)
  debug=FALSE, DIC=FALSE, bugs.dir=bd) #~~~~~ for automated testing
print(out1, 3)


# Analysis of binned data using data augmentation
delta <- 50                # Width of distance bins
xg <- seq(0, B, delta)     # Make the interval cut points
dclass <- x %/% delta + 1  # Convert distances to distance category
nD <- length(xg) -1        # N intervals = length(xg) if max(x) = B

# Bundle data
# Note data changed to include dclass, nG, bin-width delta and midpt
midpt <- xg[-1] - delta/2  # Interval mid-points
str( win.data <- list (nind=nind, nz=nz, dclass=dclass, y=y, B=B,
   delta=delta, nD=nD, midpt=midpt) )   # Bundle and summarize

# BUGS model specification
cat("
model{
  # Priors
  psi ~ dunif(0, 1)
  sigma ~ dunif(0, 1000)

  # Likelihood
  # construct conditional detection probability and Pr(x) for each bin
  for(g in 1:nD){        # midpt = mid point of each cell
    log(p[g]) <- -midpt[g] * midpt[g] / (2 * sigma * sigma)
    pi[g] <- delta / B  # probability of x in each interval
  }

  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)             # model for individual covariates
    dclass[i] ~ dcat(pi[])        # population distribution of distance class
    mu[i] <- z[i] * p[dclass[i]]  # p depends on distance class
    y[i] ~ dbern(mu[i])
  }
  # Derived quantities: Population size and density
  N <- sum(z[])
  D <- N / 60
}
",fill=TRUE, file = "model2.txt")

# Inits function
zst <- y # DA variables start at observed value of y
inits <- function(){ list (psi=runif(1), z=zst, sigma=runif(1,40,200)) }

# Parameters to save
params <- c("N", "sigma", "D")

# Unleash WinBUGS and summarize posteriors
# bd <- "c:/Program Files/WinBUGS14/"
out2 <- bugs(win.data, inits, params, "model2.txt",
  n.thin=2, n.chains=3, n.burnin=1000, n.iter=11000,
  # debug=TRUE, DIC=FALSE, bugs.dir = bd)
  debug=FALSE, DIC=FALSE, bugs.dir = bd) #~~~~~ for automated testing
print(out2, 2)


# 8.3.2 Other formulations of the distance sampling model
# ------------------------------------------------------------------------

# 8.3.3 A treatise on the integration of mathematical functions in one dimension
# ------------------------------------------------------------------------
sigma <- 2                        # normal scale (standard deviation)
curve(exp(-x^2 / (2*sigma^2)), 0, 10, frame = F)

delta <- 1                        # bin width
mid <- seq(0.5, 9.5, delta)       # 10 rectangles
f.mid <- exp(-mid^2 / (2*sigma^2))
barplot(f.mid, add=T, space=0, col="grey", width=delta)
curve(exp(-x^2 / (2*sigma^2)), 0, 10, add = TRUE, col = "blue", lwd = 3)

# Integral done using the integrate function
integrate( function(x){ exp(-x^2/(2*sigma^2)) }, lower=0, upper=100)

# Summing up the 10 rectangular areas:
areas <- f.mid * delta
sum(areas)


# 8.3.4 Bayesian analysis of point transect data
# ------------------------------------------------------------------------
### Version 1: Point count data in BUGS (conditional likelihood)
# Simulate a data set and harvest the output
set.seed(1234)
tmp <- sim.pdata(N=200, sigma=1, keep.all=FALSE, B=3)
# ~~~~~~ remove objects that mask elements of temp ~~~~~~~
rm(B, sigma, y)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
attach(tmp)

# Chop the data into bins
delta <- 0.1           # width of distance bins for approximation
xg <- seq(0, B, delta) # Make the mid points and chop up the data
midpt <- xg[-1] - delta/2

# Convert distances to categorical distances (which bin?)
dclass <- d %/% delta + 1
nD <- length(midpt) # how many intervals
nind <- length(dclass)

# Bundle and summarize data set
str( win.data <- list(midpt=midpt, delta=delta, B=B, nind=nind,
    nD=nD, dclass=dclass) )


# BUGS model specification, conditional version
cat("
model{

  # Prior for single parameter
  sigma ~ dunif(0, 10)

  # Construct cell probabilities for nG cells (rectangle approximation)
  for(g in 1:nD){    # midpt[g] = midpoint of each distance band
    log(p[g]) <- -midpt[g] * midpt[g] / (2*sigma*sigma)
    pi[g] <- (( 2 * midpt[g] ) / (B*B)) * delta
    f[g] <- p[g] * pi[g]
    fc[g] <- f[g] / pcap
  }
  pcap <- sum(f[]) # capture prob. is the sum of all rectangular areas

  # Categorical observation model
  for(i in 1:nind){
    dclass[i] ~ dcat(fc[])
  }
  # Derived quantity: population size
  N <- nind / pcap
  D<- N/(3.141*B*B)
}
",fill=TRUE, file="model3.txt")

# Inits function
inits <- function(){list (sigma=runif(1, 1, 10)) }

# Params to save
params <- c("sigma", "N","D")

# MCMC settings
# ni <- 62000   ;   nb <- 2000   ;   nt   <-   2   ;   nc <- 3
ni <- 6200   ;   nb <- 200   ;   nt   <-   1   ;   nc <- 3  # ~~~~ for testing

# Run BUGS and summarize posteriors
# bd <- "c:/Program Files/WinBUGS14/"
out3 <- bugs(win.data, inits, params, "model3.txt",
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
  debug=FALSE, bugs.dir = bd)

## Version 2: point count data (full likelihood with data augmentation)
# Do data augmentation (for same simulated data set)
M <- 400
nz <- M - nind
y <- c(rep(1, nind), rep(0, nz))
dclass <- c(dclass, rep(NA, nz))

# Bundle and summarize data set
str( win.data <- list(midpt=midpt, delta=delta, B=B, nind=nind, nD=nD, dclass=dclass, y=y, nz=nz) )

# BUGS model
cat("
model{

  # Priors
  sigma ~ dunif(0, 10)
  psi ~ dunif(0, 1)

  # Construct cell probabilities for nG cells (rectangle approximation)
  for(g in 1:nD){           # midpt[g] = midpoint of each distance band
    log(p[g]) <- -midpt[g] * midpt[g] / (2*sigma*sigma)
    pi[g] <- ((2 * midpt[g]) / (B * B)) * delta
    pi.probs[g] <- pi[g] / norm
    f[g] <- p[g] * pi[g]
    fc[g] <- f[g] / pcap   # conditional probabilities
  }
  pcap <- sum(f[])# capture prob. is the sum of all rectangular areas
  norm <- sum(pi[])

  # Categorical observation model
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    dclass[i] ~ dcat(pi.probs[])
    mu[i] <- p[dclass[i]] * z[i]
    y[i] ~ dbern(mu[i])
  }

  # Derived quantity: population size
  N <- sum(z[])
  D<- N/(3.141*B*B)

}
",fill=TRUE,file="model4.txt")

# Inits
inits <- function(){list (sigma=runif(1,1,10), psi=runif(1) ) }

# Params to save
params <- c("sigma", "N","D","psi")

# MCMC settings
# ni <- 62000   ;   nb <- 2000   ;   nt   <-   2   ;   nc <- 3  # 1.17 days
ni <- 6200   ;   nb <- 200   ;   nt   <-   2   ;   nc <- 3  # ~~~~~ for testing

# Run BUGS and summarize posteriors
out4 <- bugs(win.data, inits, params, "model4.txt",
  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
  debug=FALSE, bugs.dir = bd)

# Compare posterior summaries
print(out3,2)

print(out4,2)

# ~~~~~ maybe save the output ~~~~~~~~~~~~~~~~~~~~~~~~
save(out1, out2, out3, out4, file="AHM1_08.03_WinBUGSoutput.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~