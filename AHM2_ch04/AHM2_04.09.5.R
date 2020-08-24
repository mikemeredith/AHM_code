#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from MS dated 2019-06-17

# Approximate run time for this script: 2 hrs
# Run time with the full number of iterations: 19 hrs

# library(AHMbook)
library(unmarked)
library(AICcmodavg)

# ~~~ load crossbill data from 4.9.1 ~~~~~~~~~~
source("AHM2_04.09.1_Crossbills.R")
# ~~~~~ and fitted model fm50 from 4.9.3 ~~~~~~~~~
load("AHM2_04.09.3_fm50.RData")
c.hat <- 2.102584
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 4.9 Analysis and mapping of crossbill distribution and range dynamics in Switzerland
# ====================================================================================

#4.9.5 Forming predictions in one or two dimensions
# -------------------------------------------------

# Let elevation only go till 2250 (no forest higher up)
n.pred.points <- 100                      # not too high
ep.original <- seq(min(cb$elev), 2250, length.out = n.pred.points)
ep <- (ep.original - mean.ele)/sd.ele
fp.original <- seq(min(cb$forest), max(cb$forest), length.out = n.pred.points)
fp <- (fp.original - mean.forest)/sd.forest
first.survey <- min(as.numeric(dates), na.rm = TRUE)
dp.original <- seq(from = first.survey, to = 105, length.out = n.pred.points)
dp <- (dp.original - mean.date) / sd.date

# 4.9.5.1 Predictions of Year Effects on Occupancy, Colonization, Extinction,
#         and Detection
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Get bootstrapped estimates of SE and CI for annual occupancy
# nboot <- 1000              # number of bootstrap samples
nboot <- 100                 # ~~~~ for testing
boot.psi.hat <- array(NA, dim = c(12, nboot))
for(i in 1:nboot){ # Start loop
  cat(paste("\n ** Nonparametric bootstrap rep", i, "**\n") )
  # Draw bootstrap sample
  bootsamp <- sample(1:267, replace = TRUE)
  # Create unmarked data frame
  umfboot <- unmarkedMultFrame(y = as.matrix(cb[bootsamp,6:41]),
      siteCovs = data.frame(elev = elev[bootsamp], forest = forest[bootsamp]),
      yearlySiteCovs = list(year = year), obsCovs = list(date = DATE[bootsamp,]),
      numPrimary = 12)
  # Fit model 50, use estimates from fm50 as inits, do not compute SEs
  inits <- coef(fm50)
  fmtmp <- colext(~ (elev + I(elev^2)) * (forest + I(forest^2)) -
      I(elev^2):I(forest^2),
      ~ (year-1) + (elev + I(elev^2)) * (forest + I(forest^2)),
      ~ (year-1) + (elev + I(elev^2)) * (forest + I(forest^2)) -
      I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
      ~ (year-1) + elev + forest + I(forest^2) + date + I(date^2) +
      elev:forest, umfboot, starts = inits,
      control = list(trace = T, REPORT = 25, maxit = 500), se = FALSE)
  # Compute and save estimates of annual occupancy
  # -> Here could bootstrap additional functions of model parameters
  boot.psi.hat[,i] <- projected(fmtmp)[2,]
}

SE.occ <- apply(boot.psi.hat, 1, sd)
CI.occ <- apply(boot.psi.hat, 1, function(x) quantile(x, prob = c(0.025, 0.975)))

# ~~~~~~~ maybe save this work ~~~~~~~~~~~~
save(boot.psi.hat, file="AHM2_04.09.5_boot.psi.hat.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(unmarked)
library(AICcmodavg)
# Choose covariate values for the prediction
nd <- data.frame(year=factor(c('2001', '2002', '2003', '2004', '2005', '2006',
    '2007', '2008', '2009', '2010', '2011')), elev = rep(0,11), forest = rep(0,11))

# Prediction of annual colonization with unmarked predict function
E.col.1 <- predict(fm50, type = 'col', newdata = nd)

# With modavgPred (can inflate variances)
E.col.2 <- modavgPred(cand.set = list(fm50), newdata = nd, conf.level = 0.95,
    parm.type = 'gamma', c.hat = c.hat) # Here for c.hat of 2.12

# Compare estimates for two predict methods
print(cbind(year = 2001:2011, as.matrix(E.col.1), E.col.2$matrix.output),2)
#       year Predicted     SE   lower upper mod.avg.pred uncond.se lower.CL upper.CL
# [1,]  2001   0.47686 0.0814 3.2e-01  0.63      0.47686     0.118  2.6e-01     0.70
# [2,]  2002   0.40818 0.0923 2.5e-01  0.59      0.40818     0.134  1.9e-01     0.67
# [3,]  2003   0.21682 0.0864 9.3e-02  0.43      0.21682     0.125  6.1e-02     0.54
# [4,]  2004   0.18211 0.0758 7.6e-02  0.38      0.18211     0.110  5.0e-02     0.49
# [5,]  2005   0.43028 0.1041 2.5e-01  0.63      0.43028     0.151  1.8e-01     0.72
# [6,]  2006   0.24293 0.1347 7.1e-02  0.57      0.24293     0.195  3.8e-02     0.72
# [7,]  2007   0.04486 0.1012 4.6e-04  0.83      0.04486     0.147  5.7e-05     0.97
# [8,]  2008   0.51520 0.0999 3.3e-01  0.70      0.51520     0.145  2.5e-01     0.77
# [9,]  2009   0.00016 0.0091 4.7e-52  1.00      0.00016     0.013  1.9e-73     1.00
# [10,] 2010   0.31265 0.1008 1.5e-01  0.53      0.31265     0.146  1.1e-01     0.63
# [11,] 2011   0.03052 0.0695 3.1e-04  0.76      0.03052     0.101  4.0e-05     0.96

# ~~~~~~~~~~~~ remaining code from the MS ~~~~~~~~~~~~~~~~~~~~

E.eps <- modavgPred(cand.set = list(fm50), newdata = nd, parm.type = 'epsilon',
    c.hat = c.hat)
nd_for_p <- data.frame(year=factor(2001:2012), elev = rep(0,12), forest = rep(0,12),
    date = rep(0,12)) # 12 years and date
E.p <- modavgPred(cand.set = list(fm50), newdata = nd_for_p, parm.type = 'detect',
    c.hat = c.hat)

# Figure 4.15:
# ''''''''''''
op <- par(mfrow=c(2,2), mar=c(5,5,1,2), cex.lab = 1.5, cex.axis = 1.3)
# Plot for occupancy probability
plot(1:12, projected(fm50)[2,], pch=1, xaxt='n', xlab='',
    ylab='Occupancy', xlim = c(0.5, 12.5), ylim=c(0,1), cex = 1.5, frame = FALSE)
axis(1, at=1:12, labels=nd$year[1:12])
segments(1:12, CI.occ[1,], 1:12, CI.occ[2,], lwd = 1)

# Plot for colonization probability
plot(1:11, E.col.2$matrix.output[,1], pch=1, xaxt='n', xlab='Year',
    ylab='Colonization', xlim = c(0.5, 11.5), ylim=c(0,1), cex = 1.5, frame = FALSE)
axis(1, at=1:11, labels=nd$year[1:11])
segments(1:11, E.col.2$matrix.output[,3], 1:11, E.col.2$matrix.output[,4], lwd = 1)

# Plot for extinction probability
plot(1:11, E.eps$matrix.output[,1], pch=1, xaxt='n', xlab='Year',
    ylab='Extinction', xlim = c(0.5, 11.5), ylim=c(0,1), cex = 1.5, frame = FALSE)
axis(1, at=1:11, labels=nd$year[1:11])
segments(1:11, E.eps$matrix.output[,3], 1:11, E.eps$matrix.output[,4], lwd = 1)

# Plot for detection probability
plot(1:12, E.p$matrix.output[,1], pch=1, xaxt='n', xlab='Year',
    ylab='Detection', xlim = c(0.5, 12.5), ylim=c(0,1), cex = 1.5, frame = FALSE)
axis(1, at=1:12, labels=nd_for_p$year)
segments(1:12, E.p$matrix.output[,3], 1:12, E.p$matrix.output[,4], lwd = 1)
par(op)

# 4.9.5.2 Predictions of a single continuous covariate
# ''''''''''''''''''''''''''''''''''''''''''''''''''''

nd <- data.frame(elev = ep, forest = 0)
E.psi1 <- modavgPred(cand.set = list(fm50), newdata = nd, parm.type = 'psi',
    c.hat = c.hat)

nd <- data.frame(year=factor('2011', levels = c('2001','2002','2003','2004',
    '2005','2006','2007','2008','2009','2010','2011')), elev = ep, forest = 0)
E.col <- modavgPred(cand.set = list(fm50), newdata = nd, parm.type = 'gamma',
    c.hat = c.hat)
E.eps <- modavgPred(cand.set = list(fm50), newdata = nd, parm.type = 'epsilon',
    c.hat = c.hat)
nd_for_p <- data.frame(year=factor('2011', levels = c('2001','2002','2003',
    '2004','2005','2006','2007','2008','2009','2010','2011','2012')),
    elev = ep, forest = 0, date = 0)
E.p <- modavgPred(cand.set = list(fm50), newdata = nd_for_p, parm.type = 'detect',
    c.hat = c.hat)

# Figure 4.16
# '''''''''''
op <- par(mfrow=c(2,2))
# Plot for occupancy probability
plot(ep.original, E.psi1$matrix.output[,1], xlab = 'Elevation',
    ylab='Initial occupancy', xlim = c(200, 2300),
    ylim=c(0,1), type = 'l', frame = FALSE, lwd = 2)
polygon(c(ep.original, rev(ep.original)), c(E.psi1$matrix.output[,3],
    rev(E.psi1$matrix.output[,4])), col = 'grey', border = 'grey')
lines(ep.original, E.psi1$matrix.output[,1], lwd = 2)

# Plot for colonization probability
plot(ep.original, E.col$matrix.output[,1], xlab = 'Elevation',
    ylab='Colonization', xlim = c(200, 2300),
    ylim=c(0,1), type = 'l', frame = FALSE, lwd = 2)
polygon(c(ep.original, rev(ep.original)), c(E.col$matrix.output[,3],
    rev(E.col$matrix.output[,4])), col = 'grey', border = 'grey')
lines(ep.original, E.col$matrix.output[,1], lwd = 2)

# Plot for extinction probability
plot(ep.original, E.eps$matrix.output[,1], xlab = 'Elevation',
    ylab='Extinction', xlim = c(200, 2300),
    ylim=c(0,1), type = 'l', frame = FALSE, lwd = 2)
polygon(c(ep.original, rev(ep.original)), c(E.eps$matrix.output[,3],
    rev(E.eps$matrix.output[,4])), col = 'grey', border = 'grey')
lines(ep.original, E.eps$matrix.output[,1], lwd = 2)

# Plot for detection probability
plot(ep.original, E.p$matrix.output[,1], xlab = 'Elevation',
    ylab='Detection', xlim = c(200, 2300),
    ylim=c(0,1), type = 'l', frame = FALSE, lwd = 2)
polygon(c(ep.original, rev(ep.original)), c(E.p$matrix.output[,3],
    rev(E.p$matrix.output[,4])), col = 'grey', border = 'grey')
lines(ep.original, E.p$matrix.output[,1], lwd = 2)
par(op)


# 4.9.5.3 Predictions for two continuous covariates simultaneously
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Figure 4.17
# Get predictions first and only once
tmp <- summary(fm50)
# Predict first-year (2001) occupancy for elevation and forest cover
psipar <- tmp$psi[,1]  ;  names(psipar) <- rownames(tmp$psi)  ; psipar
pred.matrix1 <- array(NA, dim = c(100, 100)) # Predict 100x100 matrix
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(elev=ep[i], forest=fp[j])
    pred <- predict(fm50, type="psi", newdata=newData)
    pred.matrix1[i, j] <- pred$Predicted
  }
}
# Predict first-year (2001) occupancy for elevation and forest cover
psipar <- tmp$psi[,1]  ;  names(psipar) <- rownames(tmp$psi)  ; psipar
pred.matrix2 <- array(NA, dim = c(100, 100)) # Predict 100x100 matrix
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(elev=ep[i], forest=fp[j])
    pred <- predict(fm50, type="psi", newdata=newData)
    pred.matrix2[i, j] <- pred$Predicted
  }
}
# Predict extinction (2001/2002) for elevation and forest cover
extpar <- tmp$ext[,1]  ;  names(extpar) <- rownames(tmp$ext)  ;  extpar
pred.matrix3 <- array(NA, dim = c(100, 100)) # Predict 100x100 matrix
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(year=factor("2001", levels = c('2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011')), elev = ep[i], forest = fp[j])
    pred <- predict(fm50, type="ext", newdata=newData)
    pred.matrix3[i, j] <- pred$Predicted
  }
}
# Predict detection (2001) for elevation and forest cover
# (NOTE there is 1 more level now in year covariate in newData)
detpar <- tmp$det[,1]  ;  names(detpar) <- rownames(tmp$det) ;  detpar
pred.matrix4 <- array(NA, dim = c(100, 100)) # Predict 100x100 matrix
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(year=factor('2001', levels = c('2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012')), elev = ep[i], forest = fp[j], date = 0)
    pred <- predict(fm50, type="det", newdata=newData)
    pred.matrix4[i, j] <- pred$Predicted
  }
}
# Predict detection (2001) for elevation and date
# (NOTE 1 more levels in year in newData)
pred.matrix5 <- array(NA, dim = c(100, 100)) # Predict 100x100 matrix
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(year=factor('2001', levels = c('2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012')), elev = ep[i], forest = 0, date = dp[j])
    pred <- predict(fm50, type="det", newdata=newData)
    pred.matrix5[i, j] <- pred$Predicted
  }
}
# Get prediction SE for detection (2001) for elevation and date
# Inflate SE for overdispersion now
se.pred.matrix <- array(NA, dim = c(100, 100)) # Predict 100x100 matrix
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(year=factor('2001', levels = c('2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012')), elev = ep[i], forest = 0, date = dp[j])
    pred <- modavgPred(cand.set = list(fm50), newdata = newData, conf.level = 0.95, parm.type = 'detect', c.hat = c.hat)
    se.pred.matrix[i, j] <- pred$matrix.output[2] # Get inflated SE
  }
}

mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
op <- par(mfrow = c(3,2))

# Predict first-year (2001) occupancy for elevation and forest cover
image(x=ep.original, y=fp.original, z=pred.matrix1, col = mapPalette(100),
    frame = FALSE, xlab = "Elevation [m]", ylab = "Forest cover [%]",
    main = "Occupancy prob. 2001")
points(cb$elev, cb$forest, pch="+", cex=1)
contour(x=ep.original, y=fp.original, z=pred.matrix1, add = TRUE,
    col = "blue")
box()

# Predict colonisation (2001/2002) for elevation and forest cover
image(x=ep.original, y=fp.original, z=pred.matrix2, col = mapPalette(100),
    frame = FALSE, xlab = "Elevation [m]", ylab = "Forest cover [%]",
    main="Colonization prob. 2001-2002")
points(cb$elev, cb$forest, pch="+", cex=1)
contour(x=ep.original, y=fp.original, z=pred.matrix2, add = TRUE,
    col = "blue")
box()

# Predict extinction (2001/2002) for elevation and forest cover
image(x=ep.original, y=fp.original, z=pred.matrix3, col = mapPalette(100),
    axes = FALSE, xlab = "Elevation [m]", ylab = "Forest cover [%]",
    main="Extinction prob. 2001-2002")
points(cb$elev, cb$forest, pch="+", cex=1)
contour(x=ep.original, y=fp.original, z=pred.matrix3, add = TRUE,
    col = "blue")
box()

# Predict detection (2001) for elevation and forest cover
image(x=ep.original, y=fp.original, z=pred.matrix4, col = mapPalette(100),
    frame = FALSE, xlab = "Elevation [m]", ylab = "Forest cover [%]",
    main="Detection prob. 2001")
points(cb$elev, cb$forest, pch="+", cex=1)
contour(x=ep.original, y=fp.original, z=pred.matrix4, add = TRUE,
    col = "blue")
box()

# Predict detection (2001) for elevation and date
image(x=ep.original, y=dp.original, z=pred.matrix5, col = mapPalette(100),
    frame = FALSE, xlab = "Elevation [m]", ylab = "Date (1 = April 1)", 
    main="Detection prob. 2001")
matpoints(cb$elev, as.matrix(cb[,43:45]), pch="+")
contour(x=ep.original, y=fp.original, z=pred.matrix5, add = TRUE, 
    col = "blue")
box()

# Get prediction SE for detection (2001) for elevation and date
image(x=ep.original, y=dp.original, z=se.pred.matrix, col = mapPalette(100),
    frame = FALSE, xlab = "Elevation [m]", ylab = "Date (1 = April 1)",
    main="SE of Detection prob. 2001")
matpoints(cb$elev, as.matrix(cb[,43:45]), pch="+")
contour(x=ep.original, y=fp.original, z=se.pred.matrix, add = TRUE, 
    col = "blue")
box()
par(op)
