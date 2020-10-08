#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

library(AHMbook)
library(unmarked)
library(jagsUI)

# 2.4 Modeling Temporary Emigration (TE) with a three-level N-mixture model
# =========================================================================

# 2.4.3 Estimating density from point counts using a spatial TE model
# -------------------------------------------------------------------
# Grab data set from AHMbook package
data(cswa)
str(cswa)

# Prepare 'covs' and use that instead of cswa$covs for the analysis
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

# Standardize the first 2 yearlySiteCovs (time, date) and summarize
#   the unmarked data frame
ysc <- yearlySiteCovs(umfcs)
ysc[,1:2] <- scale(ysc[,1:2])
yearlySiteCovs(umfcs) <- ysc
summary(umfcs) # Not shown

# Create a list to hold the models, fit a basic model and summarize
cswamods <- list()
(cswamods$Null <- gmultmix(~offset(log(plotArea)), ~1, ~1, umfcs))
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
#           nPars    AIC delta  AICwt cumltvWt
# Wh2Wc..Wh     7 325.71  0.00 0.2900     0.29
# Wh2Wc..       6 326.10  0.39 0.2400     0.53
# Wh2Wc..T      7 326.87  1.15 0.1600     0.69
# Wh2Wc.A.      7 327.50  1.79 0.1200     0.81
# Wh2WcA..      7 327.99  2.28 0.0920     0.90
# Wh2Wc..D      7 328.10  2.39 0.0870     0.98
# Wh2..         5 331.92  6.21 0.0130     1.00
# Wc..          4 335.82 10.10 0.0018     1.00
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
cswa.Dhat <- parboot(cswamods$Wh2Wc..Wh, statistic = getD, covs = covs,
    nsim = 500, ncores=3)
plot(cswa.Dhat); abline(v = 1.109612, col = 4) # Not shown
summary(cswa.Dhat@t.star)                      # Not shown

# Goodness-of-fit assessment (ART 6 min)
gof <- parboot(cswamods$Wh2Wc..Wh, statistic = fitstats, nsim = 1000, ncores=3)
plot(gof) # Not shown
gof

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
cswapcmods$Wh2Wc. <- pcount(~1 ~offset(log(plotArea)) + woodHt + I(woodHt^2) +
    woodCov, umf)
# ~~~~~~~~~~~~~~~~~~~ these do not appear in the book ~~~~~~~~~~~~~~
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
#          nPars    AIC delta  AICwt cumltvWt
# Wh2Wc.       5 224.15  0.00 0.3000     0.30
# Wh2Wc.Wh     6 224.36  0.21 0.2700     0.58
# Wh2Wc.A      6 225.54  1.39 0.1500     0.73
# Wh2Wc.T      6 225.77  1.62 0.1300     0.86
# Wh2WcA.      6 226.03  1.88 0.1200     0.98
# Wh2.         4 229.96  5.82 0.0170     1.00
# Wc.          3 233.86  9.71 0.0024     1.00
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

# Parametric bootstrap (ART 3 min)
cswapc.Dhat <- parboot(cswapcmods$Wh2Wc., statistic = getDpc, nsim = 1000,
    covs = covs, ncores=3)
cswapc.Dhat
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
# ''''''''''''''''''''''''''''''''''''''''''''''''''
# ~~~ See below for a loop to run all three occasions ~~~
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cswamods1$Wh2Wc.Wh <- multinomPois(~woodHt ~offset(log(plotArea)) + woodHt +
    I(woodHt^2) + woodCov, umfcs1)

cswafits1 <- fitList(fits = cswamods1)
(ms.cswa1 <- modSel(cswafits1))
#          nPars    AIC delta   AICwt cumltvWt
# Wh2Wc.Wh     6 134.61  0.00 0.48533     0.49
# Wh2Wc.       5 136.84  2.23 0.15937     0.64
# Wh2Wc.D      6 137.59  2.98 0.10938     0.75
# Wh2Wc.T      6 138.01  3.39 0.08899     0.84
# Wh2.         4 138.18  3.57 0.08150     0.92
# Wh2WcA.      6 138.83  4.21 0.05904     0.98
# Wc.          3 143.00  8.39 0.00731     0.99
# Null         2 145.58 10.97 0.00202     0.99
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

# t0 = Original statistic computed from data
# t_B = Vector of bootstrap samples

# ~~~ extra code to do all 3 occasions ~~~~~~~~~~~~~~~~~~~~~~~~
# Objects to hold output
bestfit <- density <- Dhat.rem1 <- vector("list", 3)
best <- character(3)

o2y1 <- diag(3)
o2y1[upper.tri(o2y1)] <- 1
getDr1 <- function(fit, spotMaps) {
  D <- sum( predict(fit, type = "state")$Predicted)/(sum(spotMaps$parea) -0.62)
  return(D)
}

for(k in 1:3) {
  umfcs1 <- unmarkedFrameMPois(y = cswa$count[,,k],
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

  cswafits1 <- fitList(fits = cswamods1)
  (ms.cswa1 <- modSel(cswafits1))

  # pull out AIC-best model
  ( best[k] <- ms.cswa1@Full$model[1] )
  bestfit[[k]] <- cswamods1[[best[k]]]

  # Run a parametric bootstrap on density
  sum(predict(bestfit[[k]] , type = "state")$Predicted) /
      (sum(cswa$spotMaps$parea) - 0.62) # darty had no point count data

  # Compute the density for replicate 1
  density[[k]] <- getDr1(bestfit[[k]] , cswa$spotMaps)

  # Bootstrap the density estimate (ART 20 sec)
  Dhat.rem1[[k]] <- parboot(bestfit[[k]] , statistic = getDr1, nsim = 1000,
    spotMaps = cswa$spotMaps, ncores=3)
}

# Look at results
best
density

Dhat.rem1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
