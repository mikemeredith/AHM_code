#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# 4.9 Analysis and mapping of crossbill distribution and range dynamics in Switzerland
# ====================================================================================

# 4.9.1 Data manipulations and creation of unmarked data frame
# ------------------------------------------------------------

# Read in data set from AHMbook
library(AHMbook)
data(crossbillAHM)
str(cb <- crossbillAHM)

# Extract response (detection/nondetection data) and survey dates
y <- as.matrix(cb[,6:41])      # Detection/nondetection data
dates <- as.matrix(cb[,42:77]) # Survey dates

# Standardize covariates for elevation, forest and survey date
mean.ele <- mean(cb$ele, na.rm=TRUE)
sd.ele <- sd(cb$ele, na.rm=TRUE)
elev <- (cb$ele - mean.ele) / sd.ele
mean.forest <- mean(cb$forest, na.rm=TRUE)
sd.forest <- sd(cb$forest, na.rm=TRUE)
forest <- (cb$forest - mean.forest) / sd.forest
mean.date <- mean(dates, na.rm=TRUE)
sd.date <- sd(c(dates), na.rm=TRUE)
DATE <- (dates - mean.date) / sd.date
DATE[is.na(DATE)] <- 0        # Mean-impute missing dates

# Generate unmarked data frame
library(unmarked)
year <- matrix(as.character(2001:2012), 267, 12, byrow = TRUE) # Year covar.
summary(umf <- unmarkedMultFrame(y = y, siteCovs = data.frame(elev, forest),
    yearlySiteCovs = list(year = year), obsCovs=list(date = DATE), numPrimary = 12) )
