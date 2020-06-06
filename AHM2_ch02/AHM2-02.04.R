#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-01-09

library(jagsUI)
library(AHMbook)

# ~~~~ need the Green Woodpecker data prepared in 2.2 ~~~~~~~~
source("AHM2-02.02.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2.4 Modeling Temporary Emigration (TE) with a three-level N-mixture model
# =========================================================================

# 2.4.1 Doing it with bugs

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
  for (i in 1:nsites){ # Loop over sites
    M[i] ~ dpois(lambda) # Super-population size M
    for(t in 1:nyears){ # Loop over years
      N[i,t] ~ dbin(theta, M[i]) # 'Current’ population size N
      for (j in 1:nsurveys){ # Loop over reps
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t]) # Actual counts
        logit(p[i,j,t]) <- alpha0[t] + alpha[1] * DATE[i,j,t] + alpha[2] *
        pow(DATE[i,j,t],2) + alpha[3] * INT[i,j,t]
      } # end j
    } # end t
  } # end i
  # Derived quantity: Total M and total N across all surveyed sites
  totalM <- sum(M[])
  for (t in 1:nyears){
    totalN[t] <- sum(N[,t])
  } # end t
} # end model
")

# Initial values: Need for both N and M !
Nst <- apply(C, c(1,3), max, na.rm = TRUE)+1 # Inits for latent N
Nst[Nst == '-Inf'] <- 1
Mst <- apply(Nst, 1, max) # .. and for latent M
inits <- function() list(M = Mst, N =Nst, alpha = runif(3), lambda = runif(1))
# Parameters monitored
params <- c("alpha0", "alpha", "lambda", "theta", "totalM", "totalN")
# MCMC settings
na <- 1000 ; ni <- 15000 ; nt <- 10 ; nb <- 5000 ; nc <- 3
# Call JAGS (ART 21 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "Nmix3.txt", n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(3, 3)) ; traceplot(out3) ; par(mfrow = c(1, 1))
print(out3, 3) # partially shown
# mean sd 2.5% 50% 97.5% overlap0 f Rhat n.eff
# alpha0[1] -2.068 0.115 -2.296 -2.068 -1.843 FALSE 1.000 1.001 2208
# alpha0[2] -1.037 0.106 -1.245 -1.037 -0.836 FALSE 1.000 1.000 3000
# [...outputtruncated...]
# alpha0[14] -0.599 0.086 -0.766 -0.601 -0.427 FALSE 1.000 1.000 3000
# alpha[1] -0.226 0.025 -0.275 -0.225 -0.177 FALSE 1.000 1.007 3000
# alpha[2] 0.052 0.024 0.005 0.053 0.098 FALSE 0.986 1.004 1176
# alpha[3] 0.210 0.036 0.137 0.211 0.279 FALSE 1.000 1.009 250
# lambda 4.325 0.182 3.989 4.320 4.692 FALSE 1.000 1.002 898
# theta 0.312 0.012 0.289 0.312 0.335 FALSE 1.000 1.004 548
# totalM 1153.890 34.931 1090.975 1152.000 1228.000 FALSE 1.000 1.005 465
# totalN[1] 343.237 19.170 305.975 343.000 381.025 FALSE 1.000 1.000 3000
# totalN[2] 333.756 18.101 300.000 333.000 371.000 FALSE 1.000 1.001 3000
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
  for (i in 1:nsites){ # Loop over sites
    M[i] ~ dpois(lambda) # Super-population size M
    for(t in 1:nyears){ # Loop over years
      N[i,t] ~ dbin(theta[t], M[i]) # 'Current’ population size N
      for (j in 1:nsurveys){ # Loop over reps
        C[i,j,t] ~ dbin(p[i,j,t], N[i,t]) # Actual counts
        logit(p[i,j,t]) <- alpha0[t] + alpha[1] * DATE[i,j,t] + alpha[2] *
        pow(DATE[i,j,t],2) + alpha[3] * INT[i,j,t]
      } # end j
    } # end t
  } # end i
  # Derived quantity: Total M and total N across all surveyed sites
  totalM <- sum(M[])
  for (t in 1:nyears){
    totalN[t] <- sum(N[,t])
  } # end t
} # end model
")

# Initial values (omitted, same as previous model)
# Parameters monitored
params <- c("alpha0", "alpha", "lambda", "theta", "totalM", "totalN")
# MCMC settings
na <- 1000 ; ni <- 15000 ; nt <- 10 ; nb <- 5000 ; nc <- 3
# Call JAGS (ART 21 min), check convergence and summarize posteriors
out3b <- jags(bdata, inits, params, "Nmix3b.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(3, 3)) ; traceplot(out3b) ; par(mfrow = c(1, 1))
print(out3b, 3)
# mean sd 2.5% 50% 97.5% overlap0 f Rhat n.eff
# alpha0[1] -0.196 0.200 -0.585 -0.194 0.199 TRUE 0.840 1.001 1718
# ...[output truncated]...
# alpha0[14] -1.278 0.196 -1.725 -1.262 -0.934 FALSE 1.000 1.010 602
# alpha[1] -0.244 0.027 -0.295 -0.244 -0.192 FALSE 1.000 1.007 305
# alpha[2] 0.036 0.023 -0.010 0.036 0.081 TRUE 0.941 1.001 1801
# alpha[3] 0.221 0.035 0.153 0.221 0.292 FALSE 1.000 1.005 567
# lambda 4.286 0.184 3.933 4.281 4.661 FALSE 1.000 1.003 586
# theta[1] 0.078 0.010 0.059 0.078 0.100 FALSE 1.000 1.002 1560
# theta[2] 0.201 0.020 0.164 0.201 0.242 FALSE 1.000 1.002 1493
# ...[output truncated] ...
# theta[12] 0.424 0.054 0.334 0.419 0.543 FALSE 1.000 1.008 387
# theta[13] 0.456 0.059 0.361 0.449 0.589 FALSE 1.000 1.002 1461
# theta[14] 0.569 0.086 0.439 0.556 0.781 FALSE 1.000 1.014 285
# totalM 1143.470 35.000 1077.000 1142.000 1215.000 FALSE 1.000 1.005 387
# totalN[1] 88.370 7.606 76.000 87.000 106.000 FALSE 1.000 1.002 861
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
# Estimate SE LinComb (Intercept)
# 4.28 0.176 1.45 1
# Transformation: exp
# Other parameters (not shown)
backTransform(fm0p, type = "phi") # not shown
backTransform(fm0p, type = "det") # not shown

# ART ~ 2 hours
summary(fm0nb.200 <- gpcount(~1, ~1, ~1, umf, K=200, mixture="NB"))
# (omitted)
fm0nb.200@AIC
# [1] 15638.88
fm0nb@AIC
# [1] 15692.69
rbind(coef(fm0nb), coef(fm0nb.200))
# lambda(Int) phi(Int) p(Int) alpha(alpha)
# [1,] 2.749116 -2.243541 -0.8483574 -0.5088587
# [2,] 3.316829 -2.825893 -0.8891441 -0.5673470

# Run the models and summarize them
# Go-get-some-sleep warning: ART 8 hours
summary(fm9p <- gpcount(~elev + I(elev^2)+ forest, ~1, ~date + I(date^2) +
    int, umf, K = 100, mixture = "P"))
summary(fm9bp <- gpcount(~elev + I(elev^2)+ forest, ~trend, ~date +
    I(date^2) + int, umf, K = 100, mixture = "P"))
summary(fm9nb <- gpcount(~elev + I(elev^2)+ forest, ~1, ~date + I(date^2) +
    int, umf, K = 100, mixture = "NB") )
summary(fm9bnb <- gpcount(~elev + I(elev^2)+ forest, ~trend, ~date +
    I(date^2) + int, umf, K = 100, mixture = "NB") )
fl <- fitList(fm0p = fm0p, fm0nb = fm0nb, fm9p = fm9p, fm9bp = fm9bp,
    fm9nb = fm9nb, fm9bnb = fm9bnb)

modSel(fl)
# nPars AIC delta AICwt cumltvWt
# fm9bnb 11 15346.44 0.00 1.0e+00 1.00
# fm9nb 10 15471.02 124.58 8.9e-28 1.00
# fm0nb 4 15692.69 346.24 6.5e-76 1.00
# fm9bp 10 16078.24 731.80 1.2e-159 1.00
# fm9p 9 16222.69 876.25 5.3e-191 1.00
# fm0p 3 16692.70 1346.26 4.6e-293 1.00

# 2.4.3 Estimating density from point counts using a spatial TE model
# -------------------------------------------------------------------
# Grab data set from AHMbook package
data(cswa)
str(cswa)

# Prepare ’covs’ and use that instead of cswa$covs for the analysis
covs <- cswa$covs
covs$time1[is.na(covs$time1)] <- mean(covs$time1, na.rm=TRUE)
covs$date1[is.na(covs$date1)] <- mean(covs$date1, na.rm=TRUE)
# Define associated obsToY matrix
o2y <- diag(9)
o2y[upper.tri(o2y)] <- 1
o2y
# Create unmarkedFrame for the generalized multinomial mixture models:
library(unmarked)
umfcs <- unmarkedFrameGMM(y = matrix(cswa$counts, 43),
    siteCovs = covs[,c("plotArea", "patchArea", "woodHt", "woodCov")],
    yearlySiteCovs = list(
    time = covs[,c("time1", "time2", "time3")],
    date = covs[,c("date1", "date2", "date3")],
    obs = covs[,c("obs1", "obs2", "obs3")]),
    piFun = "instRemPiFun", obsToY = o2y, numPrimary = 3)
# Standardize the first 2 yearlySiteCovs (time, date) and summarize the unmarked data frame
ysc <- yearlySiteCovs(umfcs)
ysc[,1:2] <- scale(ysc[,1:2])
yearlySiteCovs(umfcs) <- ysc
summary(umfcs) # Not shown

# Create a list to hold the models, fit a basic model and summarize
cswamods <- list()
(cswamods$Null <- gmultmix(~offset(log(plotArea)), ~1, ~1, umfcs))
cswamods$Null
# Call:
# gmultmix(lambdaformula = ~offset(log(plotArea)), phiformula = ~1,
# pformula = ~1, data = umfcs)
# Abundance:
# Estimate SE z P(>|z|)
# 0.94 0.214 4.39 1.11e-05
# Availability:
# Estimate SE z P(>|z|)
# -0.302 0.334 -0.905 0.365
# Detection:
# Estimate SE z P(>|z|)
# -0.519 0.177 -2.94 0.00332
# AIC: 343.4377

(lam.cswa <- backTransform(cswamods$Null, type = "lambda"))
# Backtransformed linear combination(s) of Abundance estimate(s)
# Estimate SE LinComb (Intercept)
# 2.56 0.548 0.94 1
# Transformation: exp
(theta.cswa <- backTransform(cswamods$Null, type = "phi"))
# Backtransformed linear combination(s) of Availability estimate(s)
# Estimate SE LinComb (Intercept)
# 0.425 0.0816 -0.302 1
# Transformation: logistic
coef(lam.cswa) * coef(theta.cswa)
# [1] 1.088081

# Density estimate from spot mapping
(cswaDsm <- sum(cswa$spotMaps$CSWA) / sum(cswa$spotMaps$parea))
# [1] 1.109612
mean(cswa$spotMaps$CSWA / cswa$spotMaps$parea)
# [1] 1.345462

cswamods$Wh.. <- gmultmix(~offset(log(plotArea)) + woodHt, ~1, ~1, umfcs)
cswamods$Wh2.. <- gmultmix(~offset(log(plotArea)) + woodHt + I(woodHt^2), ~1,~1, umfcs)
cswamods$Wc.. <- gmultmix(~offset(log(plotArea)) + woodCov, ~1, ~1, umfcs)
cswamods$A.. <- gmultmix(~offset(log(plotArea)) + patchArea, ~1, ~1, umfcs)
cswamods$.A. <- gmultmix(~offset(log(plotArea)), ~patchArea, ~1, umfcs)
cswamods$A.A. <- gmultmix(~offset(log(plotArea)) + patchArea, ~patchArea, ~1, umfcs)
cswamods$..T <- gmultmix(~offset(log(plotArea)), ~1, ~time, umfcs)
cswamods$..D <- gmultmix(~offset(log(plotArea)), ~1, ~date, umfcs)
cswamods$..O <- gmultmix(~offset(log(plotArea)), ~1, ~obs, umfcs)
cswamods$..TxD <- gmultmix(~offset(log(plotArea)), ~1, ~time*date, umfcs)
cswamods$Wh2Wc.. <- gmultmix(~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, ~1, ~1, umfcs)
cswamods$Wh2WcA.. <- gmultmix(~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov + patchArea, ~1, ~1, umfcs)
cswamods$Wh2Wc.A. <- gmultmix(~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, ~patchArea, ~1, umfcs)
cswamods$Wh2Wc..T <- gmultmix(~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, ~1, ~time, umfcs)
cswamods$Wh2Wc..D <- gmultmix(~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, ~1, ~date, umfcs)
cswamods$Wh2Wc..Wh <- gmultmix(~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, ~1, ~woodHt, umfcs)
# Create a fitList and compute a model selection table
cswafits <- fitList(fits = cswamods)
modSel(cswafits)
# nPars AIC delta AICwt cumltvWt
# Wh2Wc..Wh 7 325.71 0.00 2.9e-01 0.29
# Wh2Wc.. 6 326.10 0.39 2.4e-01 0.53
# Wh2Wc..T 7 326.87 1.15 1.6e-01 0.69
# Wh2Wc.A. 7 327.50 1.79 1.2e-01 0.81
# Wh2WcA.. 7 327.99 2.28 9.2e-02 0.90
# Wh2Wc..D 7 328.10 2.39 8.7e-02 0.98
# Wh2.. 5 331.92 6.21 1.3e-02 1.00
# Wc.. 4 335.82 10.10 1.8e-03 1.00
# .. truncated ..

# The big Table of Everything
(msout <- as(modSel(cswafits), "data.frame"))
write.csv(msout, "msout.csv")

# Top model
cswamods$Wh2Wc..Wh
# Call:
# gmultmix(lambdaformula = ~offset(log(plotArea)) + woodHt + I(woodHt^2) +
# woodCov, phiformula = ~1, pformula = ~woodHt, data = umfcs)
# Abundance:
# Estimate SE z P(>|z|)
# (Intercept) -2.30 0.953 -2.41 0.01595
# woodHt 4.25 1.402 3.03 0.00243
# I(woodHt^2) -1.44 0.510 -2.83 0.00462
# woodCov 2.22 0.776 2.87 0.00416
# Availability:
# Estimate SE z P(>|z|)
# -0.902 0.527 -1.71 0.0872
# Detection:
# Estimate SE z P(>|z|)
# (Intercept) 0.198 0.495 0.40 0.689
# woodHt -0.558 0.375 -1.49 0.136
# AIC: 325.7116

# Compute N / A (birds / ha)
getD <- function(fit, covs) {
lam <- predict(fit, type="lambda")$Predicted
  theta <- matrix(predict(fit, type = "phi")$Predicted, length(lam))
  sum(lam) / sum(covs$plotArea / theta[,1])
}
# Density estimate
getD(cswamods$Wh2Wc..Wh, covs=covs)
# [1] 1.092638
# Bootstrap assessment of uncertainty. ART ~ 1.5 mins
# cswa.Dhat <- parboot(cswamods$Wh2Wc..Wh, statistic = getD, covs = covs, nsim = 500)
cswa.Dhat <- parboot(cswamods$Wh2Wc..Wh, statistic = getD, covs = covs,
    nsim = 500, ncores=3)
plot(cswa.Dhat); abline(v = 1.109612, col = 4) # Not shown
summary(cswa.Dhat@t.star) # Not shown
# Goodness-of-fit assessment (ART 6 min)
(gof <- parboot(cswamods$Wh2Wc..Wh, statistic = fitstats, nsim = 1000, ncores=3))
plot(gof) # Not shown
# Parametric Bootstrap Statistics:
# t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
# SSE 63.6 -7.532 13.36 0.692
# Chisq 374.5 15.820 113.75 0.299
# freemanTukey 66.4 0.902 7.11 0.445
# t_B quantiles:
# 0% 2.5% 25% 50% 75% 97.5% 100%
# SSE 37 49 62 70 80 101 139
# Chisq 195 253 307 341 387 533 2076
# freemanTukey 40 51 61 66 70 79 88
# t0 = Original statistic compuated from data
# t_B = Vector of bootstrap samples

## Create the unmarkedFrame
summary(umf <- unmarkedFramePCount(y = apply(cswa$count, c(1,3), sum),
    siteCovs = covs[,c("plotArea", "patchArea", "woodHt", "woodCov")],
    obsCovs = list(time = covs[,c("time1", "time2", "time3")],
    date = covs[,c("date1", "date2", "date3")],
    obs = covs[,c("obs1", "obs2", "obs3")] ) ) )
# Create a list to hold the models. Note default choice of K here
cswapcmods <- list()
(cswapcmods$Null <- pcount(~1 ~offset(log(plotArea)), umf))
# Call:
# pcount(formula = ~1 ~ offset(log(plotArea)), data = umfcspc)

# Abundance:
# Estimate SE z P(>|z|)
# 0.94 0.214 4.39 1.11e-05
# Detection:
# Estimate SE z P(>|z|)
# -0.319 0.332 -0.961 0.337
# AIC: 241.481
# Back-transform density and availability probability
(lam.cswapc <- backTransform(cswapcmods$Null, type = "state"))
# Backtransformed linear combination(s) of Abundance estimate(s)
# Estimate SE LinComb (Intercept)
# 2.56 0.548 0.94 1
(theta.cswapc <- backTransform(cswapcmods$Null, type = "det"))
# Backtransformed linear combination(s) of Detection estimate(s)
# Estimate SE LinComb (Intercept)
# 0.421 0.0808 -0.319 1
# Estimate of density from this model
coef(lam.cswapc) * coef(theta.cswapc)
# [1] 1.077969
sum(cswa$spotMaps$CSWA) / sum(cswa$spotMaps$parea)
# [1] 1.109612
mean(cswa$spotMaps$CSWA / cswa$spotMaps$parea)
# [1] 1.345462
# Now we fit a whole suite of models (ART really quick)
cswapcmods$Wh. <- pcount(~1 ~offset(log(plotArea)) + woodHt, umf)
cswapcmods$Wh2. <- pcount(~1 ~offset(log(plotArea)) + woodHt +
    I(woodHt^2), umf)
cswapcmods$Wc. <- pcount(~1 ~offset(log(plotArea)) + woodCov, umf)
cswapcmods$A. <- pcount(~1 ~offset(log(plotArea)) + patchArea, umf)
cswapcmods$.T <- pcount(~time ~offset(log(plotArea)), umf)
cswapcmods$.D <- pcount(~date ~offset(log(plotArea)), umf)
cswapcmods$.O <- pcount(~obs ~offset(log(plotArea)), umf)
cswapcmods$.TxD <- pcount(~time*date ~offset(log(plotArea)), umf)
# ~~~~~~~~~~~~~~~~~~~ these do not appear in the book ~~~~~~~~~~~~~~
cswapcmods$Wh2Wc. <- pcount(~1 ~offset(log(plotArea)) + woodHt + I(woodHt^2) +
    woodCov, umf)
cswapcmods$Wh2WcA. <- pcount(~1 ~offset(log(plotArea)) + woodHt + I(woodHt^2) +
    woodCov + patchArea, umf)
cswapcmods$Wh2Wc.T <- pcount(~time ~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umf)
cswapcmods$Wh2Wc.D <- pcount(~date ~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umf)
cswapcmods$Wh2Wc.Wh <- pcount(~woodHt ~ offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umf)
cswapcmods$Wh2Wc.A <- pcount(~patchArea ~ offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umf)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Organize the model fits into a fitList and produce AIC table
cswapcfits <- fitList(fits = cswapcmods)
(ms.cswapc <- modSel(cswapcfits))
# nPars AIC delta AICwt cumltvWt
# Wh2Wc. 5 224.15 0.00 3.0e-01 0.30
# Wh2Wc.Wh 6 224.36 0.21 2.7e-01 0.58
# Wh2Wc.A 6 225.54 1.39 1.5e-01 0.73
# Wh2Wc.T 6 225.77 1.62 1.3e-01 0.86
# Wh2WcA. 6 226.03 1.88 1.2e-01 0.98
# Wh2. 4 229.96 5.82 1.7e-02 1.00
# Wc. 3 233.86 9.71 2.4e-03 1.00
# ... truncated ...
# Create predictions of E(N) for the AIC-best model for each plot and
# compute the average density
X <- model.matrix(~woodHt + I(woodHt^2) + woodCov, covs)

lam <- exp(X %*% coef(cswapcmods$Wh2Wc., type = "state") +
    log(covs$plotArea))
sum(lam) / sum(covs$plotArea)
# [1] 3.750416
# Parametric bootstrap approach here
# Compute N / A (birds / ha)
getDpc <- function(fit, covs) {
  lam <- predict(fit, type="state")$Predicted
  sum(lam) / sum(covs$plotArea)
}
# Adjusted density
getDpc(cswapcmods$Wh2Wc., covs = covs)
# [1] 3.750416
(cswaDsm <- sum(cswa$spotMaps$CSWA) / sum(cswa$spotMaps$parea))
# [1] 1.109612
# Naive density
mean(cswa$spotMaps$CSWA / cswa$spotMaps$parea)
# [1] 1.345462
# Parametric bootstrap (ART 3 min)
( cswapc.Dhat <- parboot(cswapcmods$Wh2Wc., statistic = getDpc, nsim = 1000,
    covs = covs, ncores=3) )
  # Parametric Bootstrap Statistics:
  # t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
  # 1 3.75 -2.02 5.72 0.584
  # t_B quantiles:
  # 0% 2.5% 25% 50% 75% 97.5% 100%
  # [1,] 1.6 2.2 3.1 4.1 5.7 27 48
  # t0 = Original statistic compuated from data
  # t_B = Vector of bootstrap samples
# Get (parametric) bootstrap SE and CI
(SE <- sd(cswapc.Dhat@t.star))
(CI <- quantile(cswapc.Dhat@t.star, c(0.025, 0.975)))

# Removal models fit to each time period (ART quick)
# Occasion 1
o2y1 <- diag(3)
o2y1[upper.tri(o2y1)] <- 1
umfcs1 <- unmarkedFrameMPois(y = cswa$count[,,1],
    siteCovs = covs[,c("plotArea", "patchArea", "woodHt", "woodCov",
    "time1", "date1", "obs1")],
    piFun = "instRemPiFun", obsToY = o2y1)
sc1 <- siteCovs(umfcs1)
colnames(sc1)[5:7] <- c("time", "date", "obs")
sc1[, 5:6] <- scale(sc1[, 5:6])
siteCovs(umfcs1) <- sc1
cswamods1 <- list()
cswamods1$Null <- multinomPois(~1 ~offset(log(plotArea)), umfcs1)
cswamods1$Wh. <- multinomPois(~1 ~offset(log(plotArea)) + woodHt, umfcs1)
cswamods1$Wh2. <- multinomPois(~1 ~offset(log(plotArea)) + woodHt +
I(woodHt^2), umfcs1)
cswamods1$Wc. <- multinomPois(~1 ~offset(log(plotArea)) + woodCov, umfcs1)
# ~~~~~~~~~~~~~~~ these models not in book ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cswamods1$A. <- multinomPois(~1 ~offset(log(plotArea)) + patchArea, umfcs1)
cswamods1$.T <- multinomPois(~time ~offset(log(plotArea)), umfcs1)
cswamods1$.D <- multinomPois(~date ~offset(log(plotArea)), umfcs1)
cswamods1$.O <- multinomPois(~obs ~offset(log(plotArea)), umfcs1)
cswamods1$.TxD <- multinomPois(~time*date ~    offset(log(plotArea)), umfcs1)
cswamods1$Wh2Wc. <- multinomPois(~1 ~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umfcs1)
cswamods1$Wh2WcA. <- multinomPois(~1 ~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov + patchArea, umfcs1)
cswamods1$Wh2Wc.T <- multinomPois(~time ~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umfcs1)
cswamods1$Wh2Wc.D <- multinomPois(~date ~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umfcs1)
cswamods1$Wh2Wc.Wh <- multinomPois(~woodHt ~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umfcs1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cswafits1 <- fitList(fits = cswamods1)
(ms.cswa1 <- modSel(cswafits1))
# nPars AIC delta AICwt cumltvWt
# Wh2Wc.Wh 6 134.61 0.00 0.48533 0.49
# Wh2Wc. 5 136.84 2.23 0.15937 0.64
# Wh2Wc.D 6 137.59 2.98 0.10938 0.75
# Wh2Wc.T 6 138.01 3.39 0.08899 0.84
# Wh2. 4 138.18 3.57 0.08150 0.92
# Wh2WcA. 6 138.83 4.21 0.05904 0.98
# Wc. 3 143.00 8.39 0.00731 0.99
# Null 2 145.58 10.97 0.00202 0.99
# ... truncated ...

# Run a parametric bootstrap on density
sum(predict(cswamods1$Wh2Wc.Wh, type = "state")$Predicted) /
    (sum(cswa$spotMaps$parea) - 0.62) # darty had no point count data
# [1] 0.4809083
getDr1 <- function(fit, spotMaps) {
  D <- sum( predict(fit, type = "state")$Predicted)/(sum(spotMaps$parea) -0.62)
  return(D)
}
# Compute the density for replicate 1
getDr1(cswamods1$Wh2Wc.Wh, cswa$spotMaps)
# [1] 0.4809083
# Bootstrap the density estimate (ART 20 sec)
Dhat.rem1 <- parboot(cswamods1$Wh2Wc.Wh, statistic = getDr1, nsim = 1000,
    spotMaps = cswa$spotMaps, ncores=3)
# Summarize the bootstrap output
Dhat.rem1
# Parametric Bootstrap Statistics:
# t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
# 1 0.481 -9.54 159 0.461
# t_B quantiles:
# 0% 2.5% 25% 50% 75% 97.5% 100%
# [1,] 0.26 0.33 0.43 0.5 0.56 1.4 3937
# t0 = Original statistic compuated from data
# t_B = Vector of bootstrap samples


