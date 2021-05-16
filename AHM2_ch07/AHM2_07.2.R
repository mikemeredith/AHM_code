#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7 : MODELING FALSE POSITIVES
# ====================================
# Code from proofs dated 2020-08-19


library(AHMbook)
library(unmarked)
library(AICcmodavg)

# 7.2 Basic occupancy models with false positives
# ===============================================

# 7.2.1 Modeling unclassified false positives in unmarked
# --------------------------------------------------------

# Simulation settings
set.seed(1)            # Initialize RNGs
nsites <- 200          # number of sites (i = 1, ..., nsites=M)
nsurveys <- 7          # number of visits (j = 1, ..., nsurveys=J)
psi <- 0.6             # expected occupancy probability
p <- 0.7               # detection probability (p_11)
fp <- 0.05             # false-positive error probability (p_10)

# Simulate occupancy states and encounter histories
z <- rbinom(nsites, 1, psi)                     # occupancy states
y <- matrix(NA, nrow = nsites, ncol = nsurveys) # empty matrix for detections
for(i in 1:nsites){
  pr_yequals1 <- p*z[i] + fp*(1 - z[i])         # p11 + p10
  y[i,] <- rbinom(nsurveys, 1, pr_yequals1)     # realized observations
}

# Number of false-positive detections per occasion
apply(y[z==0,]>0, 2, sum)
# [1] 8 4 3 4 1 8 7

# Number of false-negative detections per occasion
apply(y[z==1,]==0, 2, sum)
# [1] 39 32 33 39 44 34 42

type <- c(0, 7, 0)

# Build the unmarkedFrame
library(unmarked)
summary(umf <- unmarkedFrameOccuFP(y = y, type = type)) # not shown

largerp11 <- qlogis(c(0.5, 0.7, 0.1))
largerp10 <- qlogis(c(0.5, 0.1, 0.7))

(m1 <- occuFP(detformula = ~1,     # model for p_11
    FPformula = ~1,                # model for p_10
    stateformula = ~1,             # model for psi
    data = umf,                    # umarkedFrameOccuFP object
    starts = largerp11) )          # add p_10 < p_11 constraint

# Occupancy (logit-scale):
# Estimate   SE    z P(>|z|)
#    0.323 0.15 2.15  0.0316

# Detection (logit-scale):
# Estimate SE z P(>|z|)
# 0.762 0.083 9.18 4.32e-20

# false positive (logit-scale):
# Estimate SE z P(>|z|)
# -2.69 0.199 -13.5 1.2e-41

# AIC: 1550.633

( m2 <- occuFP(detformula = ~1,    # model for p_11
    FPformula = ~1,                # model for p_10
    stateformula = ~1,             # model for psi
    data = umf,                    # umarkedFrameOccuFP object
    starts = largerp10) )          # add p_11 < p_10 constraint

m1@AIC ; m2@AIC
# [1] 1550.633
# [1] 1550.633

cbind("m1" = plogis( coef(m1) ), "m2" = plogis( coef(m2)))
#                  m1         m2
# psi(Int) 0.57997672 0.42001094
# p(Int)   0.68185780 0.06352388
# fp(Int)  0.06352318 0.68186426


# 7.2.2  Case study: transience induced false positives in water voles
# ---------------------------------------------------------------------

# Look at the water vole data (in the AHMbook package)
data(waterVoles)
wv <- waterVoles
head(wv, 3)
#   Patch y1 y2 y3 Year
# 1 cla01  1  1  0 2009
# 2 cla02  1  1  0 2009
# 3 cla03  1  1  0 2009

# Make the false-positive umf (note not all sites surveyed in each year)
summary(wv.umf <- unmarkedFrameOccuFP(y = wv[,c("y1","y2","y3")],
    siteCovs = wv[,c("Year"), drop = FALSE], type = c(0, 3, 0)) ) # not shown

# 'Means' parameterization of these models
stvals <- list("null" = qlogis(c(0.5, 0.5, 0.5, 0.7, 0.1)),
    "p11.t" = qlogis(c(0.5, 0.5, 0.5, 0.7, 0.7, 0.7, 0.1)),
    "p10.t" = qlogis(c(0.5, 0.5, 0.5, 0.7, 0.1, 0.1, 0.1)),
    "both.t" = qlogis(c(0.5, 0.5, 0.5, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1)))
cand.mods <- list(
    "null" = occuFP(detformula = ~1, FPformula = ~1,
        stateformula = ~Year-1, data = wv.umf, starts = stvals$null),
    "p11.t" = occuFP(detformula = ~Year-1, FPformula = ~1,
        stateformula = ~Year-1, data = wv.umf, starts = stvals$p11.t),
    "p10.t" = occuFP(detformula = ~1, FPformula = ~Year-1,
        stateformula = ~Year-1, data = wv.umf, starts = stvals$p10.t),
    "both.t" = occuFP(detformula = ~Year-1, FPformula = ~Year-1,
        stateformula = ~Year-1, data = wv.umf, starts = stvals$both.t))

# ~~~ modSelFP is now deprecated, use the 'unmarked' functions instead ~~~
# (modTab <- modSelFP(cand.mods))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(modTab <- modSel(fitList(fits=cand.mods)))
#        nPars     AIC  dAIC AICwt cuWt
# both.t     9 1054.35  0.00  0.99 0.99
# p10.t      7 1063.29  8.94  0.01 1.00
# p11.t      7 1073.17 18.82  0.00 1.00
# null       5 1073.30 18.95  0.00 1.00

library(AICcmodavg)
aictab(cand.mods)
#        K    AICc Delta_AICc AICcWt Cum.Wt      LL
# both.t 9 1054.93       0.00   0.99   0.99 -518.18
# p10.t  7 1063.65       8.72   0.01   1.00 -524.65
# null   5 1073.49      18.56   0.00   1.00 -531.65
# p11.t  7 1073.53      18.60   0.00   1.00 -529.58

# Select the AIC-top model
topmod <- cand.mods$both.t

# Values for prediction
pred.df <- data.frame(Year = factor(c(2009, 2010, 2011)))

# Predict
p1 <- predict(topmod, type = "det", newdata = pred.df)
p2 <- predict(topmod, type = "fp", newdata = pred.df)
p3 <- predict(topmod, type = "state", newdata = pred.df)

# Create a data frame of predictions
preds <- rbind(p1, p2, p3)
preds$Type <- rep(factor(c("Detection", "False Positive", "Occupancy")), each=3)
preds$Year <- rep(c("2009", "2010", "2011"), times=3) # produce estimates for Fig. 7.2

# ~~~~~ extra code for Fig 7.2 ~~~~~~~~~~~~
years <- 2009:2011
op <- par(mfrow = c(1,3))
plot(years, p1[,1], xlab = "Year", ylab = "p", type = "b", pch = 20,
    cex = 3, ylim = c(0, 1), cex.lab = 2, frame = FALSE, xaxt="n")
segments(x0 = years, y0 = p1[,1] - 1.96*p1[,2], x1 = years,
    y1 = p1[,1] + 1.96*p1[,2])
title("Detection probability", cex.main = 2)
axis(1, at=years)
plot(years, p2[,1], xlab = "Year", ylab = "fp", type = "b", pch = 20,
    cex = 3, ylim = c(0,1), cex.lab = 2, frame = FALSE, xaxt="n")
segments(x0 = years, y0 = p2[,1] - 1.96*p2[,2], x1 = years,
    y1 = p2[,1] + 1.96*p2[,2])
title("False positive probability", cex.main = 2)
axis(1, at=years)
plot(years, p3[,1], xlab = "Year", ylab = "psi", type = "b", pch = 20,
    cex = 3, ylim = c(0,1), cex.lab = 2, frame = FALSE, xaxt="n")
segments(x0 = years, y0 = p3[,1] - 1.96*p3[,2], x1 = years,
    y1 = p3[,1] + 1.96*p3[,2])
title("Occupancy probability", cex.main = 2)
axis(1, at=years)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
