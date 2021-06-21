#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 3 : HIERARCHICAL MODELS OF SURVIVAL
# ===========================================
# Code from proofs dated 2020-08-18

library(AHMbook)
library(jagsUI)

# 3.2 Basic Cormack-Jolly-Seber Models
# ====================================

# 3.2.1 The basic CJS model as a state-space, or hierarchical, model
# ------------------------------------------------------------------
# (no code)

# 3.2.2 A simulation function for data under basic CJS models
# -----------------------------------------------------------

# Fig. 3.1 (of course, have to load AHMbook first)
set.seed(1)
str(data <- simCJS(n.occ = 6, n.marked = 20, phi = 0.7, p = 0.4,
    show.plot = TRUE))
# List of 9
# $ n.occ   : num 6
# $ n.marked: num [1:5] 20 20 20 20 20
# $ phi     : num [1:5] 0.7 0.7 0.7 0.7 0.7
# $ p       : num [1:5] 0.4 0.4 0.4 0.4 0.4
# $ z       : num [1:100, 1:6] 1 1 1 1 1 1 1 1 1 1 ...
# $ ch      : num [1:100, 1:6] 1 1 1 1 1 1 1 1 1 1 ...
# $ f       : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
# $ n.ind   : num 100
# $ n.alive : num [1:6] 20 33 44 53 51 42

# Some settings for the function
str(data <- simCJS()) # Implicit defaults
str(data <- simCJS(n.occ = 20, phi = rep(0.7, 19),
    p = rep(0.4, 19))) # More occasions
str(data <- simCJS(n.occ = 6,
    n.marked = c(10, 20, 30, 40, 50)))                # Different numbers of annual
                                                      # release totals
str(data <- simCJS(n.occ = 20, phi = 0.9, p = 0.4))   # More occasions and higher survival
str(data <- simCJS(n.marked = 100))                   # More individuals
str(data <- simCJS(phi = c(0.1, 0.3, 0.5, 0.7, 0.9))) # phi time-dep
str(data <- simCJS(p = c(0.1, 0.3, 0.5, 0.7, 0.9)))   # p time-dep

# How to simulate CJS data with random time effects for phi and p
mean.phi <- 0.8 # Mean of logit-normal survival
mean.p <- 0.4   # ..... recapture
sd.lphi <- 1    # SD of logit-normal survival
sd.lp <- 2      # ..... recapture
nyear <- 20

op <- par(mfrow = c(1, 2))
hist(plogis(rnorm(10^6, qlogis(mean.phi), sd.lphi)), col = 'gray',
    main = 'Distribution of phi', xlab = 'Apparent survival (phi)', freq = FALSE)
hist(plogis(rnorm(10^6, qlogis(mean.p), sd.lp)), col = 'gray',
    main = 'Distribution of p', xlab = 'Recapture (p)', freq = FALSE) # not shown
par(op)

phi <- plogis(rnorm(nyear-1, qlogis(mean.phi), sd.lphi))
p <- plogis(rnorm(nyear-1, qlogis(mean.p), sd.lp))
str(data <- simCJS(n.occ = nyear, n.marked = 20, phi = phi, p = p))

# 3.2.3 Bayesian and frequentist analysis of the simplest possible CJS model
# --------------------------------------------------------------------------

# Generate data set
set.seed(1)
str(data <- simCJS())

# A quickie with wiqid (Frequentist inference)
library(wiqid)
(mle <- survCJS(data$ch, model=list(phi~1, p~1), freq=1, ci = 0.95))
(mle <- survCJS(data$ch)) # Same

# Call: survCJS(DH = data$ch)

# Real values (duplicates omitted):
# est lowCI uppCI
# phi1 0.6173 0.4893 0.7310
# p1 0.4385 0.2947 0.5933

# AIC: 244.6032

# Bundle and summarize data set
str(bdata <- list(y = data$ch, f = data$f, n.ind = data$n.ind,
    n.occ = data$n.occ))
# List of 4
# $ y    : num [1:100, 1:6] 1 1 1 1 1 1 1 1 1 1 ...
# $ f    : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
# $ n.ind: num 100
# $ n.occ: num 6

# Specify model in BUGS language
cat(file = "cjs1.txt","
model {

  # Priors
  phi ~ dunif(0, 1) # Apparent survival
  p ~ dunif(0, 1) # Recapture

  # Likelihood
  for (i in 1:n.ind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occ){
      # State process: the latent alive/dead state
      z[i,t] ~ dbern(z[i,t-1] * phi)
      # Observation process: relates true state to observed state, y [ ch
      y[i,t] ~ dbern(z[i,t] * p)
    }
  }
}
")

# Initial values
inits <- function(){list(z = zinit(data$ch))}

# Parameters monitored
params <- c("phi", "p", "z")

# MCMC settings
na <- 5000 ; ni <- 120000 ; nt <- 10 ; nb <- 20000 ; nc <- 3

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "cjs1.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,3))  #  ~~~ replace with 'layout' argument
traceplot(out1, layout=c(2,3))
print(out1, 3)
#           mean    sd  2.5%   50% 97.5% overlap0 f Rhat n.eff
# phi      0.618 0.062 0.501 0.616 0.745    FALSE 1    1  4897
# p        0.443 0.076 0.302 0.441 0.597    FALSE 1    1  9372
# z[1,1]   1.000 0.000 1.000 1.000 1.000    FALSE 1   NA     1
# z[2,1]   1.000 0.000 1.000 1.000 1.000    FALSE 1   NA     1
# [ ... Output truncated ... ]
# z[99,6]  0.479 0.500 0.000 0.000 1.000     TRUE 1    1 27604
# z[100,6] 0.476 0.499 0.000 0.000 1.000     TRUE 1    1 30000

print(cbind(mle$real[c(1,6),], out1$summary[1:2, c(1, 3, 7)]), 3)
#        est lowCI uppCI  mean  2.5% 97.5%
# phi1 0.617 0.489 0.731 0.618 0.501 0.745
# p1   0.438 0.295 0.593 0.443 0.302 0.597

# Plot estimates of z matrix (Fig. 3.2)
mapPalette <- colorRampPalette(c("white", "black"))
image(x = 1:data$n.occ, y = 1:data$n.ind, z = t(out1$mean$z), col = mapPalette(10),
    axes = TRUE, xlab = "Year", ylab = "Individual")


# 3.2.4 Bayesian analysis of the fully time-dependent CJS model
# -------------------------------------------------------------

# Generate a slightly larger data set
set.seed(24)
str(data <- simCJS(n.occ = 6, n.marked = 100, phi = runif(5, 0.2, 0.8),
    p = runif(5, 0.2, 0.9)))

# Bundle and summarize data set
str(bdata <- list(y = data$ch, f = data$f, n.ind = data$n.ind,
    n.occ = data$n.occ))

# Specify model in BUGS language
cat(file = "cjs2.txt","
model {

  # Priors
  for(t in 1:(n.occ-1)){
    phi[t] ~ dunif(0, 1)   # Apparent survival
    p[t] ~ dunif(0, 1)     # Recapture
  }

  # Likelihood
  for (i in 1:n.ind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occ){
      # State process: the latent alive/dead state
      z[i,t] ~ dbern(z[i,t-1] * phi[t-1])
      # Observation process: relates true state to observed state, y [ ch
      y[i,t] ~ dbern(z[i,t] * p[t-1])
    } #t
  } #i
}
")

# Initial values
inits <- function(){list(z = zinit(data$ch))}

# Parameters monitored
params <- c("phi", "p", "z")

# MCMC settings
na <- 1000 ; ni <- 20000 ; nt <- 5 ; nb <- 10000 ; nc <- 3

# Call JAGS (ART 2 min), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "cjs2.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(2,3))  #  ~~~ replace with 'layout' argument
traceplot(out2, layout=c(2,3))
print(out2, 3) # not shown

print(cbind(truth=c(data$phi, data$p), out2$summary[1:10, c(1,3,7)]), 2)
#        truth mean 2.5% 97.5%
# phi[1]  0.38 0.31 0.20  0.46
# phi[2]  0.33 0.34 0.22  0.47
# phi[3]  0.62 0.61 0.48  0.76
# phi[4]  0.51 0.48 0.38  0.60
# phi[5]  0.60 0.49 0.21  0.95
# p[1]    0.84 0.76 0.48  0.96
# p[2]    0.40 0.36 0.20  0.54
# p[3]    0.73 0.61 0.46  0.75
# p[4]    0.76 0.86 0.69  0.97
# p[5]    0.38 0.54 0.22  0.97
