#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18 and MS dated 2019-01-04

# Approximate run time for this script: 25 mins

# library(AHMbook)
library(unmarked)
library(sp)
library(raster)
# library(rgdal)  # ~~~~ not necessary ~~~~

# ~~~ load crossbill data from 4.9.1 ~~~~~~~~~~
source("AHM2_04.09.1_Crossbills.R")
# ~~~~ and model 50 and c.hat from 4.9.3 ~~~~~~
load("AHM2_04.09.3_fm50.RData")
c.hat <- 2.102584
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 4.9 Analysis and mapping of crossbill distribution and range dynamics in Switzerland
# ====================================================================================

# 4.9.6 Prediction in geographical space to produce species distribution maps
# ---------------------------------------------------------------------------

require(unmarked)
data(Switzerland) # Load Swiss landscape data in unmarked
str(ch <- Switzerland) #

# ~~~~ extra code for Figure 4.18 ~~~~~~~~~~~~~~~~
# Map elevation and forest cover
op <- par(mfrow = c(1, 2), mar = c(1,2,1,7), cex = 1.5)
# par(mfrow = c(2, 1), mar = c(2,2,2,5))
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = ch$elevation))
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "",
    zlim = c(0, 4500))
par(cex = 1.5)
r2 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = ch$forest))
mapPalette2 <- colorRampPalette(c("grey", "lightgreen", "darkgreen"))
plot(r2, col = mapPalette2(100), axes = FALSE, box = FALSE, main = "",
    zlim = c(0, 100))
par(op)
# ~~~~~~~~~~~~ end of extra code ~~~~~~~~~~~~~~~~

nyears <- 12

# Standardize using same values as for original data set in model fitting
EP <- (ch$elevation - mean.ele) / sd.ele
FP <- (ch$forest - mean.forest) / sd.forest

# Extract parameter estimates for each parameter type
tmp <- summary(fm50)
psipar <- tmp$psi[,1]         # Params in model for initial occupancy
names(psipar) <- rownames(tmp$psi)
colpar <- tmp$col[,1]         # Params in model for colonization
names(colpar) <- rownames(tmp$col)
extpar <- tmp$ext[,1]         # Params in model for extinction
names(extpar) <- rownames(tmp$ext)
detpar <- tmp$det[,1]         # Params in model for detection
names(detpar) <- rownames(tmp$det)

# Create arrays to contain predictions for each km2 in Switzerland
# First-year occupancy
pred.occ1 <- numeric(nrow(ch))
# Colonization and extinction prob. (note one fewer year)
pred.col <- pred.ext <- matrix(NA, nrow = nrow(ch), ncol = (nyears-1))
# Occupancy (all years) and detection prob.
pred.occ <- pred.det <- matrix(NA, nrow = nrow(ch), ncol = nyears)

# Compute predicted first-year occupancy probability
pred.occ1 <- plogis(psipar[1] + psipar["elev"]*EP +
    psipar["I(elev^2)"]*EP^2 + psipar["forest"]*FP +
    psipar["I(forest^2)"]*FP^2 + psipar["elev:forest"]*EP*FP +
    psipar["I(elev^2):forest"]*EP^2*FP)
pred.occ[,1] <- pred.occ1

# Compute predicted colonization and extinction probability and occupancy
# for later years (j>1)
for(j in 1:(nyears-1)){
  # Colonization
  pred.col[,j] <- plogis(colpar[j] + colpar["elev"]*EP +
      colpar["I(elev^2)"]*EP^2 + colpar["forest"]*FP +
      colpar["I(forest^2)"]*FP^2 + colpar["elev:forest"]*EP*FP +
      colpar["elev:I(forest^2)"]*EP*FP^2 +
      colpar["I(elev^2):forest"]*EP^2*FP +
      colpar["I(elev^2):I(forest^2)"]*EP^2*FP^2)
  # Extinction
  pred.ext[,j] <- plogis(extpar[j] + extpar["elev"]*EP +
      extpar["I(elev^2)"]*EP^2 + extpar["forest"]*FP +
      extpar["elev:forest"]*EP*FP + extpar["I(elev^2):forest"]*EP^2*FP)
  # Compute occupancy for later years recursively
  pred.occ[,(j+1)] <- pred.occ[,j] * (1-pred.ext[,j]) +
      (1-pred.occ[,j]) * pred.col[,j]
}

# Compute predicted detection probability (for average date)
for(j in 1:nyears){
  pred.det[,j] <- plogis(detpar[j] + detpar["elev"]*EP +
      detpar["forest"]*FP + detpar["I(forest^2)"]*FP^2 + detpar["date"] * 0 +
      detpar["I(date^2)"] * 0 + detpar["elev:forest"]*EP*FP)
}

# Obtain prediction uncertainty for Swiss colonization 2001-2002
library(AICcmodavg)
se.pred.col <- matrix(NA, nrow = nrow(ch), ncol = 1)
newData <- data.frame(year = factor('2001',
    levels = c('2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011')),
    elev = EP, forest = FP) # NOTE: Only 11 yearly intervals
system.time(
pred <- modavgPred(cand.set = list(fm50), newdata = newData,
    parm.type = 'gamma', c.hat = c.hat)) # Takes a while
se.pred <- pred$matrix.output[,2] # Grab inflated SE


# ~~~~~~~~~~ remaining code from MS dated 2019-01-04 ~~~~~~~~~~~~
all.pred <- cbind(pred.occ, pred.col, pred.ext, pred.det)
colnames(all.pred) <- c("Occupancy 2001", "Occupancy 2002", "Occupancy 2003", "Occupancy 2004", "Occupancy 2005", "Occupancy 2006", "Occupancy 2007", "Occupancy 2008", "Occupancy 2009", "Occupancy 2010", "Occupancy 2011", "Occupancy 2012",
  "Colonization 2001-2002", "Colonization 2002-2003", "Colonization 2003-2004", "Colonization 2004-2005","Colonization 2005-2006", "Colonization 2006-2007",
"Colonization 2007-2008","Colonization 2008-2009","Colonization 2009-2010", "Colonization 2010-2011","Colonization 2011-2012",
  "Extinction 2001-2002", "Extinction 2002-2003", "Extinction 2003-2004", "Extinction 2004-2005", "Extinction 2005-2006", "Extinction 2006-2007",
  "Extinction 2007-2008","Extinction 2008-2009","Extinction 2009-2010", "Extinction 2010-2011","Extinction 2011-2012",
  "Detection 2001", "Detection 2002", "Detection 2003", "Detection 2004", "Detection 2005", "Detection 2006", "Detection 2007", "Detection 2008","Detection 2009","Detection 2010","Detection 2011","Detection 2012")

# We can look at the annual averages (over Switzerland) of all parameters in that matrix. Again, we see a fair amount of temporal variation.

# Get annual average probabilities
annual.averages <- apply(all.pred, 2, mean)
print(annual.averages, 2)

        # Occupancy 2001         Occupancy 2002         Occupancy 2003
               # 0.21377                0.31155                0.36604
        # Occupancy 2004         Occupancy 2005         Occupancy 2006
               # 0.34339                0.29860                0.35713
        # Occupancy 2007         Occupancy 2008         Occupancy 2009
               # 0.37615                0.31723                0.43310
        # Occupancy 2010         Occupancy 2011         Occupancy 2012
               # 0.34673                0.35677                0.27867
# Colonization 2001-2002 Colonization 2002-2003 Colonization 2003-2004
               # 0.24024                0.20423                0.10661
# Colonization 2004-2005 Colonization 2005-2006 Colonization 2006-2007
               # 0.08939                0.21728                0.12082
# Colonization 2007-2008 Colonization 2008-2009 Colonization 2009-2010
               # 0.02228                0.26553                0.00015
# Colonization 2010-2011 Colonization 2011-2012   Extinction 2001-2002
               # 0.15488                0.01562                0.59230
  # Extinction 2002-2003   Extinction 2003-2004   Extinction 2004-2005
               # 0.39085                0.44038                0.53398
  # Extinction 2005-2006   Extinction 2006-2007   Extinction 2007-2008
               # 0.47358                0.31517                0.46910
  # Extinction 2008-2009   Extinction 2009-2010   Extinction 2010-2011
               # 0.04267                0.46304                0.46818
  # Extinction 2011-2012         Detection 2001         Detection 2002
               # 0.54178                0.29342                0.32279
        # Detection 2003         Detection 2004         Detection 2005
               # 0.38541                0.28391                0.44855
        # Detection 2006         Detection 2007         Detection 2008
               # 0.18446                0.30675                0.35432
        # Detection 2009         Detection 2010         Detection 2011
               # 0.16564                0.43442                0.23425
        # Detection 2012
               # 0.41922

# Compute range size and CI for every year
# ----------------------------------------
rs <- apply(all.pred[,1:nyears], 2, sum)
# print(cbind('Range size (km2)' = rs), 0)  # ~~~ digits=0 is invalid
round(cbind('Range size (km2)' = rs), 0)   # not shown

# Calculate bootstrap CIs for range size
# Prepare covariates and constants
# nyears <- 12
# EP <- (ch$elevation - mean.ele) / sd.ele
# FP <- (ch$forest - mean.forest) / sd.forest

# Create arrays to contain predictions for each km2 in Switzerland
# First-year occupancy
pred.occ1 <- numeric(nrow(ch))
# Colonization and extinction prob. (note one fewer year)
pred.col <- pred.ext <- matrix(NA, nrow = nrow(ch), ncol = (nyears-1))
# Occupancy (all years) and detection prob.
pred.occ <- pred.det <- matrix(NA, nrow = nrow(ch), ncol = nyears)

# Get bootstrapped estimates of SE and CI for range size (rs)
# nboot <- 1000           # number of bootstrap samples
nboot <- 10  # ~~~ for testing
boot.rs.hat <- array(NA, dim = c(12, nboot))
system.time(
for(i in 1:nboot){      # Start loop
  cat(paste("\n ** Nonparametric bootstrap rep", i, "**\n") )
  # Draw bootstrap sample
  bootsamp <- sample(1:267, replace = TRUE)
  # Create unmarked data frame
  umfboot <- unmarkedMultFrame(y= as.matrix(cb[bootsamp,6:41]), siteCovs=data.frame(elev = elev[bootsamp], forest = forest[bootsamp]), yearlySiteCovs=list(year=year), obsCovs=list(date=DATE[bootsamp,]), numPrimary=12)
  # Fit model 50, use estimates from fm50 as inits, do not compute SE's
  inits <- coef(fm50)
  fmtmp <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
     I(elev^2):I(forest^2),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
    I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
  ~ (year-1) + elev + forest + I(forest^2) + date + I(date^2) +
    elev:forest, umfboot, starts = inits,
    control=list(trace=T, REPORT=25, maxit=500), se = FALSE)

  # Take estimates and plug them in to estimate occupancy in every km2 of Switzerland and then sum those up, just as we did for the actual analysis
# Extract parameter estimates for each parameter type
tmp <- summary(fmtmp)
psipar <- tmp$psi[,1]             # Params in model for initial occupancy
names(psipar) <- rownames(tmp$psi)
colpar <- tmp$col[,1]             # Params in model for colonization
names(colpar) <- rownames(tmp$col)
extpar <- tmp$ext[,1]             # Params in model for extinction
names(extpar) <- rownames(tmp$ext)
detpar <- tmp$det[,1]             # Params in model for detection
names(detpar) <- rownames(tmp$det)

# Compute predicted first-year occupancy probability
pred.occ1 <- plogis(psipar[1] + psipar["elev"]*EP +
   psipar["I(elev^2)"]*EP^2 + psipar["forest"]*FP +
   psipar["I(forest^2)"]*FP^2 + psipar["elev:forest"]*EP*FP +
   psipar["I(elev^2):forest"]*EP^2*FP)
pred.occ[,1] <- pred.occ1

# Compute predicted colonization and extinction probability and occupancy for later years (j>1)
for(j in 1:(nyears-1)){
  # Colonisation
  pred.col[,j] <- plogis(colpar[j] + colpar["elev"]*EP +
    colpar["I(elev^2)"]*EP^2 + colpar["forest"]*FP +
    colpar["I(forest^2)"]*FP^2 + colpar["elev:forest"]*EP*FP +
    colpar["elev:I(forest^2)"]*EP*FP^2 +
    colpar["I(elev^2):forest"]*EP^2*FP +
    colpar["I(elev^2):I(forest^2)"]*EP^2*FP^2)

  # Extinction
  pred.ext[,j] <- plogis(extpar[j] + extpar["elev"]*EP +
    extpar["I(elev^2)"]*EP^2 + extpar["forest"]*FP +
    extpar["elev:forest"]*EP*FP + extpar["I(elev^2):forest"]*EP^2*FP)

  # Compute occupancy for later years recursively
  pred.occ[,(j+1)] <- pred.occ[,j] * (1-pred.ext[,j]) +
    (1-pred.occ[,j]) * pred.col[,j]
}

  # Add up range size (rs) and save
  boot.rs.hat[,i] <- apply(pred.occ, 2, sum)
} )  # 10 took 9 mins

# Get 'bootstrapped 95% CI'
CI.rs <- apply(boot.rs.hat, 1, function(x) quantile(x, prob = c(0.025, 0.975)))

# Plot range size with 95% CIs - figure 4.19
plot(2001:2012, rs, xlab = "Year", ylab = "Number of occupied 1km2 quadrats",
    pch = 16, ylim = c(5000, 20000), type = "b", frame = FALSE, cex = 2)
segments(2001:2012, CI.rs[1,], 2001:2012, CI.rs[2,])


map.fn <- function(data = all.pred, elev.limit = 2250){
  #   Produces Swiss maps of occupancy-related stuff.
  #   Needs Swiss landscape data (from 'unmarked') in object named 'ch'
  #   Input is matrix 'all.pred', which has the same number of rows as ch (i.e., 42275)

  # Take following by hand to produce Figure 4.20 in book
  # 2 and 9
  # 13 and 15
  # 24 and 25
  # 43 and 44

  # devAskNewPage(ask = TRUE)

  # Load packages
  require(sp)   ;   require(raster) # require(rgdal) # ~~~~ not necessary ~~~~

  # Define custom color palette
  mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

  # Now plot all parameters
  nparam <- ncol(data)
  for(i in 1:nparam) {
    # Select column in prediction array
    sel.param <- data[,i]
    param.name <- colnames(data)[i]
    cat("*** Map for:", param.name, "(", i,") ***\n")
    r <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = sel.param))
    # Create mask for elevation
    elevRaster <- rasterFromXYZ(cbind(ch$x, ch$y, ch$elevation))
    elevRaster[elevRaster > elev.limit]<-NA
    r <- mask(r, elevRaster)
    # Plot the map using custom color palette
    plot(r, col = mapPalette(100), axes = FALSE, box = FALSE,
        main = param.name, zlim = c(0, 1))
  }
}

oldask <- devAskNewPage(ask = dev.interactive(orNone=TRUE))

map.fn()  # Browse through all 46 estimated quantities # Figure 4.20

devAskNewPage(oldask)

# ~~~~~~~~ extra code for Figure 4.21 ~~~~~~~~~~~
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
r <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = se.pred))
elevRaster <- rasterFromXYZ(cbind(ch$x, ch$y, ch$elevation))
elevRaster[elevRaster > 2250] <- NA
r <- mask(r, elevRaster)
plot(r, col = mapPalette(100), axes=FALSE, box=FALSE, main="", zlim=c(0, 0.3))
points(cb$coordx, cb$coordy, pch = '+', cex = 0.8) # Add sample quads
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
