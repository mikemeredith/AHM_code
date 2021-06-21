#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 15 hrs

library(jagsUI)
library(AHMbook)

# ~~~~ need the Green Woodpecker data prepared in 2.2 ~~~~~~~~
source("AHM2_02.02.R")
# ~~~~ and this from 2.3.4 ~~~~~~~~~
nyears <- 14
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2.4 Modeling Temporary Emigration (TE) with a three-level N-mixture model
# =========================================================================

# 2.4.1 Doing it with bugs
# ------------------------

# Bundle data
str(bdata <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3], DATE = DATE, INT = INT))

# Specify model in BUGS language
cat(file = "Nmix3.txt","
model {

  # Priors
  lambda ~ dunif(0, 100)
  theta ~ dunif(0, 1)
  for (t in 1:nyears){
    alpha0[t] <- logit(mean.p[t])
    mean.p[t] ~ dunif(0, 1)
  }
  for (k in 1:3){
    alpha[k] ~ dnorm(0, 0.01)
  }

  # Likelihood
  # Ecological model for true abundance
  for (i in 1:nsites){                      # Loop over sites
    M[i] ~ dpois(lambda)                    # Super-population size M
    for(t in 1:nyears){                     # Loop over years
      N[i,t] ~ dbin(theta, M[i])            # 'Current' population size N
      for (j in 1:nsurveys){                # Loop over reps
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t])   # Actual counts
        logit(p[i,j,t]) <- alpha0[t] + alpha[1] * DATE[i,j,t] + alpha[2] *
        pow(DATE[i,j,t],2) + alpha[3] * INT[i,j,t]
      }
    }
  }

  # Derived quantity: Total M and total N across all surveyed sites
  totalM <- sum(M[])
  for (t in 1:nyears){
    totalN[t] <- sum(N[,t])
  }
}
")

# Initial values: Need for both N and M !
Nst <- apply(C, c(1,3), max, na.rm = TRUE)+1 # Inits for latent N
Nst[Nst == '-Inf'] <- 1
Mst <- apply(Nst, 1, max)                    # .. and for latent M
inits <- function() list(M = Mst, N =Nst, alpha = runif(3), lambda = runif(1))

# Parameters monitored
params <- c("alpha0", "alpha", "lambda", "theta", "totalM", "totalN")

# MCMC settings
# na <- 1000 ; ni <- 15000 ; nt <- 10 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1500 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ testing

# Call JAGS (ART 21 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "Nmix3.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out3)
print(out3, 3) # partially shown
#                mean     sd     2.5%      50%    97.5% overlap0     f  Rhat n.eff
# alpha0[1]    -2.068  0.115   -2.296   -2.068   -1.843    FALSE 1.000 1.001  2208
# alpha0[2]    -1.037  0.106   -1.245   -1.037   -0.836    FALSE 1.000 1.000  3000
# [... output truncated ...]
# alpha0[14]   -0.599  0.086   -0.766   -0.601   -0.427    FALSE 1.000 1.000  3000
# alpha[1]     -0.226  0.025   -0.275   -0.225   -0.177    FALSE 1.000 1.007  3000
# alpha[2]      0.052  0.024    0.005    0.053    0.098    FALSE 0.986 1.004  1176
# alpha[3]      0.210  0.036    0.137    0.211    0.279    FALSE 1.000 1.009   250
# lambda        4.325  0.182    3.989    4.320    4.692    FALSE 1.000 1.002   898
# theta         0.312  0.012    0.289    0.312    0.335    FALSE 1.000 1.004   548
# totalM     1153.890 34.931 1090.975 1152.000 1228.000    FALSE 1.000 1.005   465
# totalN[1]   343.237 19.170  305.975  343.000  381.025    FALSE 1.000 1.000  3000
# totalN[2]   333.756 18.101  300.000  333.000  371.000    FALSE 1.000 1.001  3000
# [... output truncated ...]

# Specify model in BUGS language
cat(file = "Nmix3b.txt","
model {

  # Priors
  lambda ~ dunif(0, 100)
  for (t in 1:nyears){
    theta[t] ~ dunif(0,1)
    alpha0[t] <- logit(mean.p[t])
    mean.p[t] ~ dunif(0, 1)
  }
  for (k in 1:3){
    alpha[k] ~ dnorm(0, 0.01)
  }

  # Likelihood
  # Ecological model for true abundance
  for (i in 1:nsites){                    # Loop over sites
    M[i] ~ dpois(lambda)                  # Super-population size M
    for(t in 1:nyears){                   # Loop over years
      N[i,t] ~ dbin(theta[t], M[i])       # 'Current' population size N
      for (j in 1:nsurveys){              # Loop over reps
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t]) # Actual counts
        logit(p[i,j,t]) <- alpha0[t] + alpha[1] * DATE[i,j,t] + alpha[2] *
        pow(DATE[i,j,t],2) + alpha[3] * INT[i,j,t]
      }
    }
  }

  # Derived quantity: Total M and total N across all surveyed sites
  totalM <- sum(M[])
  for (t in 1:nyears){
    totalN[t] <- sum(N[,t])
  }
}
")

# Initial values (omitted, same as previous model)

# Parameters monitored
params <- c("alpha0", "alpha", "lambda", "theta", "totalM", "totalN")

# MCMC settings
# na <- 1000 ; ni <- 15000 ; nt <- 10 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1500 ; nt <- 1; nb <- 500 ; nc <- 3  # ~~~ for testing

# Call JAGS (ART 21 min), check convergence and summarize posteriors
out3b <- jags(bdata, inits, params, "Nmix3b.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3, 3))  # ~~~ no longer needed
traceplot(out3b)
print(out3b, 3)
#                mean     sd     2.5%      50%    97.5% overlap0     f  Rhat n.eff
# alpha0[1]    -0.196  0.200   -0.585   -0.194    0.199     TRUE 0.840 1.001  1718
# ...[output truncated]...
# alpha0[14]   -1.278  0.196   -1.725   -1.262   -0.934    FALSE 1.000 1.010   602
# alpha[1]     -0.244  0.027   -0.295   -0.244   -0.192    FALSE 1.000 1.007   305
# alpha[2]      0.036  0.023   -0.010    0.036    0.081     TRUE 0.941 1.001  1801
# alpha[3]      0.221  0.035    0.153    0.221    0.292    FALSE 1.000 1.005   567
# lambda        4.286  0.184    3.933    4.281    4.661    FALSE 1.000 1.003   586
# theta[1]      0.078  0.010    0.059    0.078    0.100    FALSE 1.000 1.002  1560
# theta[2]      0.201  0.020    0.164    0.201    0.242    FALSE 1.000 1.002  1493
# ...[output truncated]...
# theta[12]     0.424  0.054    0.334    0.419    0.543    FALSE 1.000 1.008   387
# theta[13]     0.456  0.059    0.361    0.449    0.589    FALSE 1.000 1.002  1461
# theta[14]     0.569  0.086    0.439    0.556    0.781    FALSE 1.000 1.014   285
# totalM     1143.470 35.000 1077.000 1142.000 1215.000    FALSE 1.000 1.005   387
# totalN[1]    88.370  7.606   76.000   87.000  106.000    FALSE 1.000 1.002   861
# . . .

# 2.4.2 N-mixture models with TE in unmarked
# ------------------------------------------

# Wide format is required in unmarked
ywide <- matrix(C, nrow = dim(C)[1])
DATEwide <- matrix(DATE, nrow = dim(C)[1])
INTwide <- matrix(INT, nrow = dim(C)[1])
nyears <- dim(C)[3]
yearmat <- col( matrix(1, nrow = nrow(ywide), ncol = nyears) )

# Set everything up in an unmarkedFrame
library(unmarked)
summary(umf <- unmarkedFrameGPC(y = ywide,
    siteCovs = data.frame(elev = elev, forest = forest),
    numPrimary = nyears,
    yearlySiteCovs = data.frame(trend = c(t(yearmat - 7.5))),
    obsCovs = list(date = DATEwide, int = INTwide)) )

# Fit some models and produce some summaries
# Go-drink-coffee warning: ART 18-30 mins
summary(fm0p <- gpcount(~1, ~1, ~1, umf, K=100, mixture="P"))
summary(fm0nb <- gpcount(~1, ~1, ~1, umf, K=100, mixture="NB"))

# Mean super-population size on natural scale
backTransform(fm0p, type="lambda")

# Backtransformed linear combination(s) of Abundance estimate(s)

# Estimate    SE LinComb (Intercept)
#     4.28 0.176    1.45           1

# Transformation: exp
# Other parameters (not shown)
backTransform(fm0p, type = "phi") # not shown
backTransform(fm0p, type = "det") # not shown

# We check the value of K here by investigating whether AIC changes as we increase K:
# ART ~ 2 hours
summary(fm0nb.200 <- gpcount(~1, ~1, ~1, umf, K=200, mixture="NB"))
# (omitted)

fm0nb.200@AIC
# [1] 15638.88
fm0nb@AIC
# [1] 15692.69

rbind(coef(fm0nb), coef(fm0nb.200))
#      lambda(Int)  phi(Int)     p(Int) alpha(alpha)
# [1,]    2.749116 -2.243541 -0.8483574   -0.5088587
# [2,]    3.316829 -2.825893 -0.8891441   -0.5673470

# Run the models and summarize them
# Go-get-some-sleep warning: ART 8 hours
summary(fm9p <- gpcount(~elev + I(elev^2)+ forest, ~1,
    ~date + I(date^2) + int, umf, K = 100, mixture = "P"))
summary(fm9bp <- gpcount(~elev + I(elev^2)+ forest, ~trend,
    ~date + I(date^2) + int, umf, K = 100, mixture = "P"))
summary(fm9nb <- gpcount(~elev + I(elev^2)+ forest, ~1,
    ~date + I(date^2) + int, umf, K = 100, mixture = "NB") )
summary(fm9bnb <- gpcount(~elev + I(elev^2)+ forest, ~trend,
    ~date + I(date^2) + int, umf, K = 100, mixture = "NB") )
fl <- fitList(fm0p = fm0p, fm0nb = fm0nb, fm9p = fm9p, fm9bp = fm9bp,
    fm9nb = fm9nb, fm9bnb = fm9bnb)
modSel(fl)
#        nPars      AIC   delta    AICwt cumltvWt
# fm9bnb    11 15346.44    0.00  1.0e+00        1
# fm9nb     10 15471.02  124.58  8.9e-28        1
# fm0nb      4 15692.69  346.24  6.5e-76        1
# fm9bp     10 16078.24  731.80 1.2e-159        1
# fm9p       9 16222.69  876.25 5.3e-191        1
# fm0p       3 16692.70 1346.26 4.6e-293        1
