#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

library(AHMbook)
library(R2WinBUGS)
bd <- "C:/WinBUGS14"  # location od the "WinBUGS14.exe" application

# 10.10 Multi-scale occupancy models
# ==================================


data <- sim3Occ()       # Execute function with default args
str(data)               # Summary of output with some comments added

# 'Null' model (model 1)
str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8,
    beta.Xpsi = 0, sd.logit.psi = 0, mean.theta = 0.6,
    theta.time.range = c(0, 0), beta.Xtheta = 0, sd.logit.theta = 0,
    mean.p = 0.4, p.time.range = c(0,0), beta.Xp = 0, sd.logit.p = 0))

# No covariate effects, no random variability (model 2)
str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8,
    beta.Xpsi = 0, sd.logit.psi = 0, mean.theta = 0.6,
    theta.time.range = c(-1, 1), beta.Xtheta = 0, sd.logit.theta = 0,
    mean.p = 0.4, p.time.range = c(-2,2), beta.Xp = 0, sd.logit.p = 0))

# All covariate effects, but no random variability (model 3)
str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8,
    beta.Xpsi = 1, sd.logit.psi = 0, mean.theta = 0.6,
    theta.time.range = c(-1, 1), beta.Xtheta = 1, sd.logit.theta = 0,
    mean.p = 0.4, p.time.range = c(-2,2), beta.Xp = -1, sd.logit.p = 0))

# Most complex model with all effects allowed for by sim function (model 4)
str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8,
    beta.Xpsi = 1, sd.logit.psi = 1, mean.theta = 0.6,
    theta.time.range = c(-1, 1), beta.Xtheta = 1, sd.logit.theta = 1,
    mean.p = 0.4, p.time.range = c(-2,2), beta.Xp = -1, sd.logit.p = 1))


set.seed(1)
str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8,
    beta.Xpsi = 1, mean.theta = 0.3, theta.time.range = c(-1, 1),
    beta.Xtheta = 1, mean.p = 0.2, p.time.range = c(-1,1), beta.Xp = -1))

# Look at data
str(data$z)           # True quadrat (pond) occurrence state
str(data$a)           # True subquadrat (water sample) occurrence state
str(data$y)           # Observed data
cbind("pond"=data$z, "sample 1"= data$a[,1], "sample 2"= data$a[,2],
    "sample 3"= data$a[,3], "sample 4"= data$a[,4], "sample 5"= data$a[,5])
which(data$z-apply(data$a, 1, max)==1) # Fungus present in pond, but not in examined samples


# Bundle and summarize data set
y <- data$y
str( win.data <- list(y = y, n.pond = dim(y)[1], n.samples = dim(y)[2],
    n.pcr = dim(y)[3], covA = data$covA, covB = data$covB, covC = data$covC) )

# Define model in BUGS language
sink("model.txt")
cat("
model {

  # Priors and model for params
  int.psi ~ dunif(0,1)         # Intercept of occupancy probability
  for(t in 1:n.samples){
    int.theta[t] ~ dunif(0,1) # Intercepts availability probability
  }
  for(t in 1:n.pcr){
    int.p[t] ~ dunif(0,1)     # Intercepts detection probability (1-PCR error)
  }
  beta.lpsi ~ dnorm(0, 0.1)    # Slopes of three covariates
  beta.ltheta ~ dnorm(0, 0.1)
  beta.lp ~ dnorm(0, 0.1)

  # 'Likelihood' (or basic model structure)
  for (i in 1:n.pond){
    # Occurrence in pond i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi) + beta.lpsi * covA[i]
    for (j in 1:n.samples){
      # Occurrence in sample j
      a[i,j] ~ dbern(mu.a[i,j])
      mu.a[i,j] <- z[i] * theta[i,j]
      logit(theta[i,j]) <- logit(int.theta[j]) + beta.ltheta * covB[i,j]
      for (k in 1:n.pcr){
        # PCR detection error process in sample k
        y[i,j,k] ~ dbern(mu.y[i,j,k])
        mu.y[i,j,k] <- a[i,j] * p[i,j,k]
        logit(p[i,j,k]) <- logit(int.p[k]) + beta.lp * covC[i,j,k]
      }
    }
    tmp[i] <- step(sum(a[i,])-0.1)
  }

  # Derived quantities
  sum.z <- sum(z[])   # Total # of occupied ponds in sample
  sum.a <- sum(tmp[]) # Total # of ponds with presence in <=1 of the 5 samples
} # end model
",fill=TRUE)
sink()

# Initial values
zst <- apply(y, 1, max)        # inits for presence (z)
ast <- apply(y, c(1,2), max)   # inits for availability (a)
inits <- function() list(z = zst, a = ast, int.psi = 0.5, beta.lpsi = 0)

# Parameters monitored
params <- c("int.psi", "int.theta", "int.p", "beta.lpsi", "beta.ltheta",
    "beta.lp", "sum.z", "sum.a")

# MCMC setting
ni <- 5000   ;   nt <- 2   ;   nb <- 1000   ;   nc <- 3

# Call WinBUGS and summarize posterior
out <- bugs(win.data, inits, params, "model.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.dir = bd) # bd="c:/WinBUGS14/"
  debug = FALSE, bugs.dir = bd) # ~~~~~ for autotesting
print(out, 3)
data$sum.z

