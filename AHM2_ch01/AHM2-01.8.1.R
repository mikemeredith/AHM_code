#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================

library(AHMbook)
library(jagsUI)

# 1.8 MODELING POPULATION DYNAMICS AT TWO TEMPORAL SCALES
# =======================================================

# 1.8.1 ANALYSIS OF SIMULATED DATA UNDER A GAUSSIAN PHENOMENOLOGICAL MODEL
# ------------------------------------------------------------------------

str(out <- simPH(npop = 18, nyear = 17, nrep = 10, date.range = 1:150,
    initial.lambda = 300, gamma.parms = c(0, 0.3), mu.range = c(50, 80),
    sigma.range = c(10, 20), p.range = c(0.4, 0.6), show.plot = TRUE))

# Execute function
str(out <- simPH()) # Implicit defaults
str(out <- simPH(show.plot = F)) # Skip plot browsing
str(out <- simPH(npop = 1)) # Only one population
str(out <- simPH(npop = 100)) # Many populations (only first 16 plotted)
str(out <- simPH(nyear = 2)) # Two years only
str(out <- simPH(nyear = 50)) # Fifty years
str(out <- simPH(date.range = 50:70)) # (Too) narrow survey season
str(out <- simPH(date.range = -100:200)) # Very long survey season
str(out <- simPH(initial.lambda = 10)) # Very small populations (and many extinctions)
str(out <- simPH(gamma.parms = c(0, 0))) # Stable population, no annual variation in gamma
str(out <- simPH(mu.range = c(50, 50))) # No variation in mu
str(out <- simPH(mu.range = c(0, 100))) # Lots of variation in mu
str(out <- simPH(sigma.range = c(10, 80))) # Lots of variation in sigma

# Create a data set
set.seed(1)
str(data <- simPH(npop = 18, nyear = 17, nrep = 10, date.range = 1:150,
    initial.lambda = 500, mu.range = c(50, 80), sigma.range = c(20, 30),
    p.range = c(1, 1)))
# Bundle and summarize data set: add pop and year info
str(bugs.data <- list(C = data$C, date = data$date, npop = data$npop,
    nyear = data$nyear, nsurvey = data$nrep, pi = pi) )
# List of 6
# $ C : int [1:18, 1:17, 1:10] 0 4 3 1 0 1 0 2 0 1 ...
# $ date : num [1:18, 1:17, 1:10] 6 14 16 27 8 33 9 9 9 2 ...
# $ npop : num 18
# $ nyear : num 17
# $ nsurvey: num 10
# $ pi : num 3.14

# Specify model in BUGS language
cat(file = "modelPH.txt","
model {
  # Priors
  # Top-level priors
  for(i in 1:npop){
    for(t in 1:nyear){
      mu[i,t] ~ dnorm(0, 0.0001)
      # curve(dnorm(x, 0, sqrt(1/ 0.0001)), -200, 200) # how's it look like ?
    }
  }
  expn1 ~ dunif(1, 2000)
  for(t in 1:(nyear-1)){
    gamma[t] ~ dunif(0.01, 10)
    sigma[t] ~ dunif(0.01, 50)
  }
  sigma[nyear] ~ dunif(0.01, 50)
  # Likelihood
  # Model for between-year dynamics
  for(i in 1:npop){
    # Initial conditions
    n1[i] ~ dpois(expn1)
    n[i,1] <- n1[i]
    # Autoregressive (Markovian) transitions from t to tD1
    for(t in 2:nyear){
      n[i,t] ~ dpois(n[i,(t-1)]*gamma[t-1])
    }
  }
  # Phenomenological within-season population model
  for(i in 1:npop){
    for(t in 1:nyear){
    for(k in 1:nsurvey){
      C[i,t,k] ~ dpois(lambda[i,t,k])
      lambda[i,t,k] <- n[i,t]*(1 / (sigma[t]*sqrt(2*pi)) )*exp( - pow((date[i,t,k] - mu
      [i,t]),2) / (2*pow(sigma[t], 2)) )
    }
    }
  }
}
")
# Initial values
nst <- 50 * apply(data$C, c(1,2), max, na.rm = TRUE)
gammast <- runif(data$nyear-1, 0.8, 1.1)
inits <- function() list(n = nst, gamma = gammast)
# Parameters monitored
params <- c("expn1", "n", "n1", "mu", "gamma", "sigma")
# MCMC settings
na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
# Call JAGS (ART 3 min), assess convergence and summarize posteriors
out12 <- jags(bugs.data, inits, params, "modelPH.txt", n.adapt = na, n.chains = nc,
  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(3,3)) ; traceplot(out12, c("expn1", "n1", "gamma", "sigma"))
summary(out12) ; View(out12) ; print(out12, dig = 3)
