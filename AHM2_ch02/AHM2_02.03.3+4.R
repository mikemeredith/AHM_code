#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 10 hrs

library(jagsUI)
library(unmarked)
library(AHMbook)

# ~~~~ need the Green Woodpecker data prepared in 2.2 ~~~~~~~~
source("AHM2_02.02.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.3 Year-stratified N-mixture model
# ===================================

# 2.3.3 Fitting the year-stratified model in unmarked (“stacking” the data)
# -------------------------------------------------------------------------

# Stack the data
nyears <- 14               # Number of years
Cstacked <- NULL
DATEstacked <- NULL
INTstacked <- NULL
for(t in 1:nyears){         # Loop over years
  Cstacked <- rbind(Cstacked,C[,,t])
  DATEstacked <- rbind(DATEstacked, DATE[,,t])
  INTstacked <- rbind(INTstacked, INT[,,t])
}
fill.in <- is.na(INTstacked) & !is.na(DATEstacked)
INTstacked[fill.in] <- 0 # almost innocuous mean imputation

# Stretch out, i.e., replicate, the site covariates
elevstretched <- rep(elev, nyears)     # Scaled 'elev’ from before
foreststretched <- rep(forest, nyears) # ditto

# Create a year and a trend variable
trend <- rep(1:nyears, each=dim(C)[1]) - 7.5

# Create and summarize an unmarked data frame
library(unmarked)
summary(umf <- unmarkedFramePCount(y = Cstacked,
    siteCovs = data.frame(elev = elevstretched,
    forest = foreststretched, trend = trend),
    obsCovs = list(date = DATEstacked, int = INTstacked)))

# Fit some models with default Poisson mixture for abundance.
# ART for this whole block of 9 models: 5 mins
system.time(fm0 <- pcount(~1 ~ 1, umf, se = FALSE)) # 13 secs
fm1 <- pcount(~1 ~ trend, umf, se = FALSE)
fm2 <- pcount(~date ~ 1, umf, se = FALSE)
fm3 <- pcount(~1 ~ elev, umf, se = FALSE)
fm4 <- pcount(~1 ~ elev + forest, umf, se = FALSE)
fm5 <- pcount(~1 ~ elev + forest + trend, umf, se = FALSE)
fm6 <- pcount(~1 ~ elev + I(elev^2)+ forest + trend, umf, se = FALSE)
fm7 <- pcount(~date ~ elev + I(elev^2)+ forest + trend, umf, se = FALSE)
fm8 <- pcount(~date + I(date^2) ~ elev + I(elev^2)+ forest + trend, umf, se = FALSE)
system.time(fm9 <- pcount(~date + I(date^2) + int ~ elev + I(elev^2)+
    forest + trend, umf, se = TRUE)) # Later we want those SEs ... # 4 mins

# Negative binomial models: models with covariates and trend
# ART for these NB models: 8 mins
system.time(fm5nb <- pcount(~1 ~ elev + forest + trend, umf,
    mixture = "NB", se = FALSE))
fm6nb <- pcount(~1 ~ elev + I(elev^2)+ forest + trend, umf,
    mixture = "NB", se = FALSE)
fm7nb <- pcount(~date ~ elev + I(elev^2)+ forest + trend, umf,
    mixture = "NB", se = FALSE)
fm8nb <- pcount(~date + I(date^2) ~ elev + I(elev^2)+ forest + trend,
    umf, mixture = "NB", se = FALSE)
system.time(fm9nb <- pcount(~date + I(date^2) + int ~ elev + I(elev^2)+
forest + trend, umf, mixture = "NB", se = FALSE))  # 2.6 mins

# Organize models into a fitList and create a model selection table
fl <- fitList(fm0 = fm0, fm1 = fm1, fm2 = fm2, fm3 = fm3, fm4 = fm4, fm5 = fm5,
    fm6 = fm6, fm7 = fm7, fm8 = fm8, fm9 = fm9, fm5nb = fm5nb, fm6nb = fm6nb,
    fm7nb = fm7nb, fm8nb = fm8nb, fm9nb = fm9nb)
modSel(fl)
#       nPars      AIC   delta   AICwt cumltvWt
# fm9nb    10 16006.57    0.00 1.0e+00        1
# fm8nb     9 16041.75   35.18 2.3e-08        1
# fm7nb     8 16042.66   36.08 1.5e-08        1
# fm6nb     7 16118.14  111.56 5.9e-25        1
# fm5nb     6 16155.55  148.97 4.5e-33        1
# fm9       9 17731.46 1724.89 0.0e+00        1
# [... output truncated ...]

# K = 200, ART 250 secs
fm9nb.K200 <- pcount(~date + I(date^2) + int ~ elev + I(elev^2) + forest +
    trend, umf, K = 200, mixture = "NB", control = list(trace = TRUE,
    REPORT = 1, maxit = 500), se = FALSE) # Also monitor fitting criterion

# K = 400, ART 500 secs
fm9nb.K400 <- pcount(~date + I(date^2) + int ~ elev + I(elev^2) + forest +
    trend, umf, K = 400, mixture = "NB", control = list(trace = TRUE,
    REPORT = 1, maxit = 500), se = FALSE)

cbind('K = 118' = coef(fm9nb), 'K = 200' = coef(fm9nb.K200),
    'K = 400' =coef(fm9nb.K400))
#                    K = 118     K = 200     K = 400
# lam(Int)        1.14020215  1.14967788  1.14967788
# lam(elev)      -0.62345703 -0.62357562 -0.62357562
# lam(I(elev^2)) -0.25303162 -0.25303890 -0.25303890
# lam(forest)     0.18713541  0.18724224  0.18724224
# lam(trend)      0.04353930  0.04348226  0.04348226
# p(Int)         -1.83171876 -1.84261596 -1.84261597
# p(date)        -0.16271954 -0.16265435 -0.16265435
# p(I(date^2))    0.03958014  0.03936762  0.03936762
# p(int)          0.25505006  0.25531341  0.25531341
# alpha(alpha)   -0.77469815 -0.77644326 -0.77644326

rbind('AIC(K = 118)' = fm9nb@AIC, 'AIC(K = 200)' = fm9nb.K200@AIC,
    'AIC(K = 400)' = fm9nb.K400@AIC) # Compare AICs
#                  [,1]
# AIC(K = 118) 16006.57
# AIC(K = 200) 16006.49
# AIC(K = 400) 16006.49

system.time(fm9nb <- pcount(~date + I(date^2) + int ~ elev + I(elev^2) +
    forest + trend, umf, K = 200, mixture = "NB",
    control = list(trace = TRUE, REPORT = 1, maxit = 500), se = TRUE))  # 9 mins

# Bootstrap simulation with 100 replicates for the AIC-best NegBin model
# ART 2 hours
# system.time(pb1 <- parboot(fm9nb, fitstats, nsim = 100, report = 1))  # 20 mins
system.time(pb1 <- parboot(fm9nb, fitstats, nsim = 100, report = 1, ncores = 3))  # 20 mins

# Call: parboot(object = fm9nb, statistic = fitstats, nsim = 100, report = 1)

# Parametric Bootstrap Statistics:
#                 t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
# SSE          11319         -916.9              897       0.8416
# Chisq        23281         1857.8              792       0.0099
# freemanTukey  4660          -90.6              124       0.7624

# t_B quantiles:
#                 0% 2.5%    25%   50%   75% 97.5%  100%
# SSE           9834 10363 11683 12289 12745 13700 15569
# Chisq        19343 19872 20951 21408 21901 23002 23360
# freemanTukey  4445  4520  4678  4761  4842  4967  5061

# t0 = Original statistic compuated from data
# t_B = Vector of bootstrap samples

# Bootstrap simulation with 100 replicates for the AIC-best Poisson model
# system.time(pb1P <- parboot(fm9, fitstats, nsim = 100, report = 1))  # 393 secs
system.time(pb1P <- parboot(fm9, fitstats, nsim = 100, report = 1, ncores=3))  # 40 mins
# Call: parboot(object = fm9, statistic = fitstats, nsim = 100, report = 1)

# Parametric Bootstrap Statistics:
#                 t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
# SSE          11232           6334            124.6            0
# Chisq        23355          12923            161.3            0
# freemanTukey  4588           1055             38.4            0

# t_B quantiles:
#                 0%  2.5%   25%   50%   75% 97.5%  100%
# SSE           4614  4653  4821  4901  4979  5147  5183
# Chisq        10080 10120 10324 10455 10538 10717 10758
# freemanTukey  3443  3455  3513  3531  3559  3604  3634

# t0 = Original statistic compuated from data
# t_B = Vector of bootstrap samples

# AIC-best (NegBin) model
fm9nb
# Call:
#   pcount(formula = ~date + I(date^2) + int ~ elev + I(elev^2) +
#   forest + trend, data = umf, K = 200, mixture = "NB",
#   se = TRUE, control = list(trace = TRUE, REPORT = 1, maxit = 500))

# Abundance:

# Estimate SE z P(>|z|)
# (Intercept) 1.1497 0.13627 8.44 3.26e-17
# elev -0.6236 0.04065 -15.34 4.18e-53
# I(elev^2) -0.2530 0.04662 -5.43 5.72e-08
# forest 0.1872 0.03782 4.95 7.41e-07
# trend 0.0435 0.00788 5.52 3.49e-08

# Detection:

# Estimate SE z P(>|z|)
# (Intercept) -1.8426 0.1478 -12.47 1.09e-35
# date -0.1627 0.0231 -7.04 1.99e-12
# I(date^2) 0.0394 0.0200 1.97 4.86e-02
# int 0.2553 0.0444 5.75 8.74e-09

# Dispersion:

# Estimate SE z P(>|z|)
# -0.776 0.0503 -15.4 9.12e-54

# AIC: 16006.49

# AIC-best Poisson model
fm9
# Call:
# pcount(formula = ~date + I(date^2) + int ~ elev + I(elev^2) +
# forest + trend, data = umf, K = 100, se = TRUE)
# Abundance:
# Estimate SE z P(>|z|)
# (Intercept) 0.0903 0.03574 2.53 1.15e-02
# elev -0.5728 0.02736 -20.94 2.53e-97
# I(elev^2) -0.1654 0.03184 -5.20 2.05e-07
# forest 0.1324 0.01942 6.82 9.24e-12
# trend 0.0466 0.00458 10.19 2.30e-24
# Detection:
# Estimate SE z P(>|z|)
# (Intercept) -0.5032 0.0397 -12.68 7.91e-37
# date -0.2251 0.0263 -8.54 1.30e-17
# I(date^2) 0.0558 0.0233 2.40 1.66e-02
# int 0.3094 0.0343 9.02 1.96e-19
# AIC: 17731.46

# ~~~~~~ figure code not in the book ~~~~~~~~~~~~

# Predictions for every datum
pred.nb <- predict(fm9nb, type="state", newdata=umf) ###, append = TRUE)
pred.p <- predict(fm9, type="state", newdata=umf)
# Figure (not shown)
plot(pred.p[,1], pred.nb[,1], xlab="N, Poisson abundance",
    ylab="N, neg. binomial abundance", pch=20)
abline(0, 1, lwd=2,col="red")

# Figure 2.5
# ----------
# Comparative plots of abundance relationships with elevation, forest cover and
#   time (i.e., trend) and detection with date underboth NegBin and Poisson
#   models
cnb <- coef(fm9nb) # Get coefficients or NegBin model
cp <- coef(fm9)        # ... and for the Poisson model

op <- par(mfrow  = c(2, 2), mar = c(5,5,4,3), cex.lab = 1.5,
    cex.lab = 1.5, cex.main = 1.2)
curve(exp(cnb[1] + cnb[2] * (x - 1190) / 637 + (cnb[3] * (x - 1190) / 637)^2),
    200, 3000, xlab = 'Elevation (m)', ylab = 'Expected abundance E(N)',
    ylim = c(0, 10), lwd = 3, col = 'red', frame = FALSE,
    main = 'Abundance predictions from NB model (red) \n and from Poisson model (blue)')
curve(exp(cp[1] + cp[2] * (x - 1190) / 637 + (cp[3] * (x - 1190) / 637)^2),
    200, 3000, lwd = 3, col = 'blue', add = TRUE)
curve(exp(cnb[1] + cnb[4] * (x - 31.9) / 27.7), 0, 100,
    xlab = 'Forest cover (%)', ylab = 'Expected abundance E(N)',
    ylim = c(0, 5), lwd = 3, col = 'red', frame = FALSE,
    main = 'Abundance predictions from NB model (red) \n and from Poisson model (blue)')
curve(exp(cp[1] + cp[4] * (x - 31.9) / 27.7), 0, 100, lwd = 3, col = 'blue',
    add = TRUE)
curve(exp(cnb[1] + cnb[5] * (x - 2005 - 7.5)), 2004, 2017, xlab = 'Year',
    ylab = 'Expected abundance E(N)', ylim = c(0, 5), lwd = 3, col = 'red',
    frame = FALSE, main = 'Abundance predictions from NB model (red) \n and from Poisson model (blue)')
curve(exp(cp[1] + cp[5] * (x - 2005 - 7.5)), 2004, 2017, lwd = 3, col = 'blue',
    add = TRUE)
curve(plogis(cnb[6] + cnb[7] * standardize2match(x, dateso) +
    cnb[8] * standardize2match(x, dateso)^2), 100, 200, xlab = 'Julian Date',
    ylab = 'Expected detection E(p)', ylim = c(0, 0.5), lwd = 3, col = 'red',
    frame = FALSE, main = 'Detection predictions from NB model (red) \n and from Poisson model (blue)')
curve(plogis(cp[6] + cp[7] * standardize2match(x, dateso) + cp[8] *
    standardize2match(x, dateso)^2), 100, 200, lwd = 3, col = 'blue',
    add = TRUE)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2.3.4 Nonparametric bootstrapping to account for potential
# dependence in stacked data
# -------------------------------------------------------------------------

# Create an array to hold coefficient estimates
# simrep <- 100
simrep <- 10  # ~~~ for testing
npbs.esti <- array(NA, dim = c(length(coef(fm9nb)), simrep))
rownames(npbs.esti) <- names(coef(fm9nb))

# Other preps
nyears <- 14
trend <- rep(1:nyears, each = dim(C)[1]) - 7.5

# Launch bootstrap
for(b in 1:simrep){
  cat(paste("\n*** Bootstrap rep", b, "***\n"))
  # Get index of chosen sites
  bs.samp <- sample(1:267, 267, replace = TRUE)
  # Repeat data preparation with bootstrap sample
  Cstacked <- NULL ; DATEstacked <- NULL ; INTstacked <- NULL
  for(t in 1:nyears){ # Loop over years
    Cstacked <- rbind(Cstacked,C[bs.samp,,t])
    DATEstacked <- rbind(DATEstacked, DATE[bs.samp,,t])
    INTstacked <- rbind(INTstacked, INT[bs.samp,,t])
  }
  elevstretched <- rep(elev[bs.samp], nyears)
  foreststretched <- rep(forest[bs.samp], nyears)
  # Create new unmarked data frame and fit model to bootstrapped data set
  umf <- unmarkedFramePCount(y = Cstacked,
      siteCovs = data.frame(elev = elevstretched, forest = foreststretched,
      trend = trend), obsCovs = list(date = DATEstacked, int = INTstacked))
  fm <- pcount(~date + I(date^2) + int ~ elev + I(elev^2) +
      forest + trend, umf, K = 200, mixture = "NB", control = list(trace = TRUE,
      REPORT = 10, maxit = 100), se = FALSE) # No SEs needed
  # Save param estimates (could also make predictions and save them !)
  npbs.esti[,b] <- coef(fm)
}

# Get bootstrap SE, CI and CI length for all 10 parameters
se.bs <- apply(npbs.esti, 1, sd)
ci.bs <- t(apply(npbs.esti, 1, function(x)quantile(x, c(0.025, 0.975))))
cil.bs <- abs(ci.bs[,2] - ci.bs[,1]) # CI length

# Extract SE and CI and compute CI length
tmp <- summary(fm9nb)
se <- c(tmp$state$SE, tmp$det$SE, tmp$alpha$SE)
ci <- rbind(confint(fm9nb, type = 'state'), confint(fm9nb, type = 'det'),
confint(fm9nb, type = 'alpha'))
cil <- abs(ci[,2] - ci[,1])

# Compare nominal and bootstrapped SE/CI/CI length
print(cbind('Nominal SE' = se, ci), digits = 2)
print(cbind('Boot SE' = se.bs, ci.bs, 'se/se.bs ratio (%)' =
    round(100*(se/se.bs),0), 'cil/cil.bs ratio (%)' = round(100*(cil/cil.bs),0)),
    digits = 2)
#                Nominal SE    0.025  0.975
# lam(Int)           0.1363  0.88259  1.417
# lam(elev)          0.0407 -0.70325 -0.544
# lam(I(elev^2))     0.0466 -0.34442 -0.162
# lam(forest)        0.0378  0.11311  0.261
# lam(trend)         0.0079  0.02803  0.059
# p(Int)             0.1478 -2.13224 -1.553
# p(date)            0.0231 -0.20797 -0.117
# p(I(date^2))       0.0200  0.00025  0.078
# p(int)             0.0444  0.16834  0.342
# alpha(alpha)       0.0503 -0.87502 -0.678


#                Boot SE   2.5%  97.5% se/se.bs ratio (%) cil/cil.bs ratio (%)
# lam(Int)        0.2263  0.742  1.589                 60                   63
# lam(elev)       0.1037 -0.787 -0.397                 39                   41
# lam(I(elev^2))  0.1297 -0.506 -0.069                 36                   42
# lam(forest)     0.0817  0.026  0.326                 46                   50
# lam(trend)      0.0087  0.027  0.063                 90                   87
# p(Int)          0.1920 -2.243 -1.540                 77                   82
# p(date)         0.0362 -0.235 -0.099                 64                   66
# p(I(date^2))    0.0264 -0.011  0.091                 76                   76
# p(int)          0.0788  0.101  0.401                 56                   58
# alpha(alpha)    0.0845 -0.940 -0.614                 60                   60

# Compute average magnitude of nominal SE/CI length in % of bootstrapped
mean(100*(se/se.bs))
mean(100*(cil/cil.bs))
# [1] 60.41515
# [1] 62.59937


