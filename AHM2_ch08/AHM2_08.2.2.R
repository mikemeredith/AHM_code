#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8 : MODELING INTERACTIONS AMONG SPECIES
# ===============================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 6 mins

library(AHMbook)
library(unmarked)
library(abind)

# 8.2 Joint occupancy models for few species with symmetric interactions
# ======================================================================

# 8.2.2 Implementation of the Rota et al. model in unmarked
# ----------------------------------------------------------

# Load libraries and check out help texts of two main functions
library(AHMbook)
library(unmarked)
?occuMulti # Check out help text for two main functions
?unmarkedFrameOccuMulti

# Read in data set and bundle up detection/nondetection data
data(MesoCarnivores)
?MesoCarnivores
str(data <- MesoCarnivores)
str(ylist <- list(bobcat = data$bobcat, coyote = data$coyote,
    redfox = data$redfox))
# List of 3
# $ bobcat: int [1:1437, 1:3] 0 0 0 0 0 0 0 0 0 0 ...
# $ coyote: int [1:1437, 1:3] 0 0 0 0 0 0 0 0 0 0 ...
# $ redfox: int [1:1437, 1:3] 0 0 0 0 0 1 0 0 0 0 ...

# Get sample sizes
nsites <- 1437      # row dimension of above data
nsurveys <- 3       # column dimension ...
nspecies <- 3

# Prepare site covariates and look at them
Dist <- scale(data$sitecovs[,'Dist_5km'])
HDens <- scale(data$sitecovs[,'HDens_5km'])
table(Trail <- data$sitecovs[,'Trail'])
head(sitecovs <- data.frame(Dist = Dist, HDens = HDens, Trail = Trail))

# Invent some random observation covariate data for illustration
date <- matrix(rnorm(n = nsites * nsurveys), ncol = nsurveys)

# Create obsCovs
str(obscovs <- list(date = date))

# Build unmarked frame and summarize it
summary(umf <- unmarkedFrameOccuMulti(y = ylist, siteCovs = sitecovs,
    obsCovs = obscovs))

# unmarkedFrame Object

# 1437 sites
# 3 species: bobcat coyote redfox
# Maximum number of observations per site: 3
# Mean number of observations per site:
# bobcat: 3 coyote: 3 redfox: 3
# Sites with at least one detection:
# bobcat: 196 coyote: 401 redfox: 161
# Tabulation of y observations:
# bobcat:
# 0 1
# 4057 254
# coyote:
# 0 1
# 3685 626
# redfox:
# 0 1
# 4054 257
# ... [ rest of output truncated ] ....

umf@fDesign
# f1[bobcat] f2[coyote] f3[redfox] f4[bobcat:coyote]
# psi[111] 1 1 1 1
# psi[110] 1 1 0 1
# psi[101] 1 0 1 0
# psi[100] 1 0 0 0
# psi[011] 0 1 1 0
# psi[010] 0 1 0 0
# psi[001] 0 0 1 0
# psi[000] 0 0 0 0
# f5[bobcat:redfox] f6[coyote:redfox] f7[bobcat:coyote:redfox]
# psi[111] 1 1 1
# psi[110] 0 0 0
# psi[101] 1 0 0
# psi[100] 0 0 0
# psi[011] 0 1 0
# psi[010] 0 0 0
# psi[001] 0 0 0
# psi[000] 0 0 0

# View order of params, helpful for setting up formulas
colnames(umf@fDesign)
# [1] "f1[bobcat]"                "f2[coyote]"
# [3] "f3[redfox]"                "f4[bobcat:coyote]"
# [5] "f5[bobcat:redfox]"         "f6[coyote:redfox]"
# [7] "f7[bobcat:coyote:redfox]"

# Specify models for occupancy in terms of natural parameters f
occ_formulae1 <- c(
    # bobcat, coyote, red fox occupancy, respectively
    # Modeled as constant
    ' ~ 1',' ~ 1',' ~ 1',
    # Two-way interactions are fixed at 0
    # bobcat:coyote, bobcat:redfox, coyote:redfox
    rep(0,3),
    # Fix 3-way interaction at 0
    0
)

# Constant detection by species
det_formulae1 <- rep(' ~ 1', 3)

# Fit model 1
(fm1 <- occuMulti(det_formulae1, occ_formulae1, umf))

# Occupancy:
# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -1.171 0.1328 -8.82 1.13e-18
# [coyote] (Intercept) -0.629 0.0744 -8.46 2.79e-17
# [redfox] (Intercept) -1.846 0.0944 -19.56 3.58e-85

# Detection:
# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -1.104 0.1396 -7.91 2.59e-15
# [coyote] (Intercept) -0.333 0.0764 -4.36 1.30e-05
# [redfox] (Intercept) -0.252 0.1182 -2.13 3.29e-02

# AIC: 6710.658

# Specify models for occupancy in terms of natural parameters f
occ_formulae2 <- c(
    # bobcat, coyote, red fox occupancy, respectively
    # Modeled as a linear function of Dist
    rep(' ~ Dist', 3),
    # same: ' ~ Dist',' ~ Dist',' ~ Dist',
    # Two-way interactions are fixed at 0
    # bobcat:coyote, bobcat:redfox, coyote:redfox
    rep(0, 3),
    # Fix 3-way interaction also at 0
    0
)
# Detection by species
det_formulae2 <- rep('~as.factor(Trail)', 3)
# det_formulae2 <- rep('~as.factor(Trail) + date', 3)

# Fit model 2
(fm2 <- occuMulti(det_formulae2, occ_formulae2, umf,
    control = list(maxit = 500, trace = TRUE, REPORT = 1)))

# Occupancy:

# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -0.5542 0.1587 -3.4930 4.78e-04
# [bobcat] Dist -0.5224 0.1293 -4.0413 5.32e-05
# [coyote] (Intercept) 0.1917 0.1058 1.8116 7.01e-02
# [coyote] Dist -0.0069 0.0845 -0.0817 9.35e-01
# [redfox] (Intercept) -1.5358 0.1195 -12.8554 8.02e-38
# [redfox] Dist -0.2912 0.1123 -2.5926 9.53e-03

# Detection:

# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -2.58 0.1619 -15.91 5.54e-57
# [bobcat] as.factor(Trail)1 1.89 0.1714 11.02 3.01e-28
# [coyote] (Intercept) -1.96 0.0998 -19.65 5.63e-86
# [coyote] as.factor(Trail)1 2.17 0.1227 17.68 5.52e-70
# [redfox] (Intercept) -1.40 0.1836 -7.61 2.80e-14
# [redfox] as.factor(Trail)1 1.63 0.2186 7.44 1.00e-13

# AIC: 6257.525

# Specify models for occupancy in terms of natural parameters f
occ_formulae3 <- c(
    # bobcat, coyote, red fox occupancy, respectively
    # Modeled as constant
    rep(' ~ 1', 3),
    # same: ' ~ 1',' ~ 1',' ~ 1',
    # Two-way interactions are intercept-only
    # bobcat:coyote, bobcat:redfox, coyote:redfox
    rep('~1', 3),
    # same: '~1', '~1', '~1',
    # Still fix 3-way interaction at 0 (for an intercept model would write: '~1')
    0
)
# Detection by species
det_formulae3 <- rep('~1', 3)

# Fit model 3
(fm3 <- occuMulti(det_formulae3, occ_formulae3, umf,
    control = list(maxit = 500, trace = TRUE, REPORT = 1)))

# Occupancy:
# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -1.76 0.179 -9.81 1.03e-22
# [coyote] (Intercept) -1.30 0.137 -9.54 1.37e-21
# [redfox] (Intercept) -2.20 0.152 -14.44 2.81e-47
# [bobcat:coyote] (Intercept) 1.72 0.262 6.56 5.48e-11
# [bobcat:redfox] (Intercept) -1.38 0.377 -3.66 2.57e-04
# [coyote:redfox] (Intercept) 1.41 0.248 5.69 1.31e-08

# Detection:
# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -1.106 0.1398 -7.91 2.59e-15
# [coyote] (Intercept) -0.331 0.0761 -4.35 1.38e-05
# [redfox] (Intercept) -0.253 0.1183 -2.13 3.28e-02

# AIC: 6626.111

confint(fm3, type = 'state') # CIs for occupancy parameters
#                                       0.025      0.975
# psi([bobcat] (Intercept))        -2.1100157 -1.4072536
# psi([coyote] (Intercept))        -1.5717770 -1.0362043
# psi([redfox] (Intercept))        -2.4983062 -1.9012367
# psi([bobcat:coyote] (Intercept))  1.2042902  2.2311451
# psi([bobcat:redfox] (Intercept)) -2.1154563 -0.6387872
# psi([coyote:redfox] (Intercept))  0.9254514  1.8992954

# Specify models for occupancy in terms of natural parameters f
occ_formulae4 <- c(
  # bobcat, coyote, red fox occupancy, respectively
  # Main species effect modeled with a covariate
  rep(' ~ Dist', 3),
  # Two-way interactions are still intercept-only
  # bobcat:coyote, bobcat:redfox, coyote:redfox
  rep('~1', 3),
  # Still fix 3-way interaction at 0
  0
)
# Detection by species
det_formulae4 <- rep('~as.factor(Trail)', 3)

# Fit model 4
(fm4 <- occuMulti(det_formulae4, occ_formulae4, umf,
    control = list(maxit = 500, trace = TRUE, REPORT = 1)))

# Occupancy:
# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -1.010 0.2393 -4.22 2.42e-05
# [bobcat] Dist -0.637 0.1438 -4.43 9.51e-06
# [coyote] (Intercept) -0.529 0.1919 -2.75 5.89e-03
# [coyote] Dist 0.160 0.0986 1.62 1.05e-01
# [redfox] (Intercept) -1.762 0.2311 -7.62 2.47e-14
# [redfox] Dist -0.463 0.1271 -3.65 2.67e-04
# [bobcat:coyote] (Intercept) 1.304 0.3223 4.05 5.20e-05
# [bobcat:redfox] (Intercept) -1.881 0.3996 -4.71 2.52e-06
# [coyote:redfox] (Intercept) 1.224 0.3179 3.85 1.18e-04

# Detection:
# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -2.59 0.162 -15.97 2.11e-57
# [bobcat] as.factor(Trail)1 1.89 0.171 11.08 1.60e-28
# [coyote] (Intercept) -1.94 0.100 -19.31 4.48e-83
# [coyote] as.factor(Trail)1 2.15 0.123 17.43 4.63e-68
# [redfox] (Intercept) -1.38 0.184 -7.50 6.46e-14
# [redfox] as.factor(Trail)1 1.61 0.219 7.34 2.10e-13

# AIC: 6213.46

# Can look at 95% CIs
confint(fm4, type = 'state') # For occupancy
confint(fm4, type = 'det') # For detection

# Specify models for occupancy in terms of natural parameters f
occ_formulae5 <- c(
    # bobcat, coyote, red fox occupancy, respectively
    # Main species effects all modeled with two covariates
    rep(' ~ Dist', 3),
    # Two-way interactions are a function of housing density
    # bobcat:coyote, bobcat:redfox, coyote:redfox
    rep('~HDens', 3),
    # Still fix 3-way interaction at 0
    0
)
# Detection by species
det_formulae5 <- rep('~as.factor(Trail)', 3)

# Fit model 5
(fm5 <- occuMulti(det_formulae5, occ_formulae5, umf,
    control = list(maxit = 500, trace = TRUE, REPORT = 1)))

# Occupancy:
# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -1.154 0.261 -4.42 9.81e-06
# [bobcat] Dist -0.648 0.146 -4.43 9.28e-06
# [coyote] (Intercept) -0.569 0.207 -2.75 5.94e-03
# [coyote] Dist 0.203 0.106 1.91 5.64e-02
# [redfox] (Intercept) -1.844 0.246 -7.50 6.57e-14
# [redfox] Dist -0.335 0.128 -2.62 8.91e-03
# [bobcat:coyote] (Intercept) 1.046 0.374 2.80 5.17e-03
# [bobcat:coyote] HDens -2.239 0.654 -3.43 6.14e-04
# [bobcat:redfox] (Intercept) -1.545 0.427 -3.62 2.95e-04
# [bobcat:redfox] HDens 1.284 0.432 2.98 2.93e-03
# [coyote:redfox] (Intercept) 1.421 0.366 3.89 1.01e-04
# [coyote:redfox] HDens 1.277 0.386 3.31 9.42e-04

# Detection:
# Estimate SE z P(>|z|)
# [bobcat] (Intercept) -2.54 0.1632 -15.57 1.22e-54
# [bobcat] as.factor(Trail)1 1.83 0.1707 10.72 7.99e-27
# [coyote] (Intercept) -2.00 0.0993 -20.13 4.41e-90
# [coyote] as.factor(Trail)1 2.16 0.1211 17.86 2.47e-71
# [redfox] (Intercept) -1.56 0.1621 -9.64 5.51e-22
# [redfox] as.factor(Trail)1 1.82 0.2044 8.92 4.63e-19

# AIC: 6123.852

modSel(fl <- fitList(fm1, fm2, fm3, fm4, fm5))
#     nPars     AIC  delta    AICwt cumltvWt
# fm5    18 6123.85   0.00  1.0e+00        1
# fm4    15 6213.46  89.61  3.5e-20        1
# fm2    15 6257.57 133.72  9.2e-30        1
# fm3     9 6626.11 502.26 8.6e-110        1
# fm1     6 6710.66 586.81 3.8e-128        1

# Create prediction covariate for Dist and form predictions
r <- range(data$sitecovs$Dist)
x <- seq(r[1], r[2], length.out = 100)
x_scaled <- (x -mean(data$sitecovs$Dist)) / sd(data$sitecovs$Dist)
nd <- data.frame(Dist = x_scaled, HDens = 0)
pred.bobcat <- predict(fm5, type = 'state', species = 'bobcat', newdata = nd,
    se.fit = TRUE, nsims = 10^5)
pred.coyote <- predict(fm5, type = 'state', species = 'coyote', newdata = nd,
    se.fit = TRUE, nsims = 10^5)
pred.redfox <- predict(fm5, type = 'state', species = 'redfox', newdata = nd,
    se.fit = TRUE, nsims = 10^5)

# ~~~ extra code for figure 8.2 ~~~~~~~~
# Plot predictions of marginal occupancy
cols <- rainbow(3)
plot(x = x, y = pred.bobcat$Predicted, col = cols[1], type='l', lwd=5,
    ylim=c(0, 1), ylab='Marginal Occupancy',
    xlab='Disturbance in 5 km radius (Dist)', frame = FALSE, lty = 2)
lines(x = x, y = pred.coyote$Predicted, col = cols[2], lwd=5, lty = 1)
lines(x = x, y = pred.redfox$Predicted, col = cols[3], lwd=5, lty = 3)
polygon(c(x, rev(x)), c(pred.bobcat[,3], rev(pred.bobcat[,4])),
    border = NA, col=rgb(1,0,0,0.2))
polygon(c(x, rev(x)), c(pred.coyote[,3], rev(pred.coyote[,4])),
    border = NA, col=rgb(0,1,0,0.2))
polygon(c(x, rev(x)), c(pred.redfox[,3], rev(pred.redfox[,4])),
    border = NA, col=rgb(0,0,1,0.2))
legend('topleft', lty= c(1,2,3), lwd=3, cex = 1.3, col=cols[c(2,1,3)],
    legend=c('Coyote', 'Bobcat', 'Red fox'), bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Predict conditional occupancy and present in a table
nd <- data.frame(Dist = 0, HDens = 0)
bob_none <- predict(fm5, type = 'state', species = 'bobcat',
    cond = c('-coyote','-redfox'), newdata = nd, nsims = 10^5)
bob_coyote <- predict(fm5,type = 'state',species = 'bobcat',
    cond = c('coyote','-redfox'), newdata = nd, nsims = 10^5)
bob_redfox <- predict(fm5,type = 'state',species = 'bobcat',
    cond = c('-coyote','redfox'), newdata = nd, nsims = 10^5)
bob_both <- predict(fm5,type = 'state',species = 'bobcat',
    cond = c('coyote','redfox'), newdata = nd, nsims = 10^5)

round(occtab <- rbind('Neither' = bob_none, 'Coyote' = bob_coyote,
  'Red fox' = bob_redfox, 'Both' = bob_both), 3)
#         Predicted    SE lower upper
# Neither     0.240 0.047 0.159 0.345
# Coyote      0.473 0.109 0.268 0.688
# Red fox     0.063 0.033 0.025 0.152
# Both        0.161 0.089 0.054 0.397

# ~~~ extra code to plot these (not shown) ~~~~~~~~
plot_dat <- rbind(bob_none, bob_coyote, bob_redfox, bob_both)
plot(1:4, plot_dat$Predicted, xlim=c(0.5,4.5), ylim=c(0,0.8), xaxt='n',
  pch=19, cex=2, ylab='Bobcat occupancy', xlab='Other species present', frame = FALSE)
axis(1, at=c(1:4), labels=c('Neither', 'Coyote', 'Red fox', 'Both'))
for (i in 1:4){
  segments(i,plot_dat$lower[i],i,plot_dat$upper[i],lwd=2)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Predict interaction parameter (Fig. 8.3, plotting code not shown)
r <- range(data$sitecovs$HDens)
x <- seq(r[1], r[2], length.out = 100)
x_scaled <- (x - mean(data$sitecovs$HDens)) / sd(data$sitecovs$HDens)
nd <- data.frame(Dist = 0, HDens = x_scaled)
pr_nocoyote <- predict(fm5, type = 'state', species = 'bobcat',
    cond = '-coyote', newdata = nd, se.fit = TRUE, nsims = 10^6)
pr_coyote <- predict(fm5, type = 'state', species = 'bobcat',
    cond = 'coyote', newdata = nd, se.fit = TRUE, nsims = 10^6)

# ~~~ extra code for figure 8.3 ~~~~~~~~
plot(x=x, y=pr_nocoyote$Predicted, type='l', col='red', lwd=5,
    ylim=c(0,1), xlab='Housing density', ylab='Bobcat occupancy',
    frame = FALSE)
lines(x=x, y=pr_coyote$Predicted, col='blue', lwd=5)
polygon(c(x, rev(x)), c(pr_nocoyote[,3], rev(pr_nocoyote[,4])),
    border = NA, col=rgb(1,0,0,0.2)) #
polygon(c(x, rev(x)), c(pr_coyote[,3], rev(pr_coyote[,4])),
    border = NA, col=rgb(0,0,1,0.2)) #
legend('topleft', lwd=3, col=c('blue', 'red'), cex = 1.5,
    legend=c('Coyote present', 'Coyote absent'), bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ plot for bobcat occupancy with vs without red fox (not shown) ~~~
pr_nofox <- predict(fm5, type='state', species='bobcat', cond='-redfox',
    newdata=nd, se.fit=TRUE, nsims = 10^5)
pr_fox <- predict(fm5, type='state', species='bobcat', cond='redfox',
    newdata=nd, se.fit=TRUE, nsims = 10^5)
plot(x=x, y=pr_nofox$Predicted, type='l', col='red', lwd=5, cex = 1.5,
    ylim=c(0,0.5), xlab='Housing density', ylab='Bobcat occupancy', frame = FALSE)
lines(x=x, y=pr_fox$Predicted, col='blue', lwd=5)
polygon(c(x, rev(x)), c(pr_nofox[,3], rev(pr_nofox[,4])),
    border = NA, col=rgb(1,0,0,0.2)) #
polygon(c(x, rev(x)), c(pr_fox[,3], rev(pr_fox[,4])),
    border = NA, col=rgb(0,0,1,0.2)) #
legend(100, 0.5, lwd=3, col=c('red','blue'), cex = 1.5,
    legend=c('Red fox absent', 'Red fox present'), bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
