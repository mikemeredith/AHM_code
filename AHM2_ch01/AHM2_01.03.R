#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

# 1.3 Crested tit count data from the Swiss MHB Breeding Bird Survey
# ==================================================================

# Get Crested Tit data and look at summaries
library(AHMbook)
data(crestedTit)
str(dat <- crestedTit)     # Marc prefers short names (comment from Mike)
C <- as.matrix(dat[,6:23]) # grab counts 1999:2016
year <- 1999:2016

# ~~~~ code to plot figure 1.3 ~~~~~~~~~~
matplot(year, t(C), main = "", type = "l", lty = 1, xlab = "Year",
    ylab = "Territory count", lwd = 3, cex.lab = 1.5, frame = FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Grab data for survey dates (for 2-3 surveys per year) and for duration
# Put into 3D array first, then summarize over reps within a year
nsite <- nrow(C)
nyear <- length(year)
datetmp <- as.matrix(dat[,24:77])
datefull <- array(datetmp, dim = c(nsite, 3, nyear))
durtmp <- as.matrix(dat[,78:131])
durfull <- array(durtmp, dim = c(nsite, 3, nyear))

# Get mean date of survey and mean survey duration for each site and year
date <- apply(datefull, c(1,3), mean, na.rm = TRUE)
dur <- apply(durfull, c(1,3), mean, na.rm = TRUE)
date[date == 'NaN'] <- NA
dur[dur == 'NaN'] <- NA
