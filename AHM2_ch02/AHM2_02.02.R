#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

# 2.2 Swiss MHB Data for the Green Woodpecker
# ===========================================

# Load data set and grab some data
library(AHMbook)
data("greenWoodpecker")
str(peckers <- greenWoodpecker)
counts <- as.matrix(peckers[,7:48])     # Counts 2004:2017
dateso <- as.matrix(peckers[,49:90])    # Survey Julian dates 2004:2017
timeso <- as.matrix(peckers[,91:132])   # Survey durations 2004:2017
into <- timeso/peckers[,'route.length'] # Survey intensity (min / km)

# Quick visualizations of counts, survey date, duration and intensityop
op <- par(mfrow = c(2, 2)) # plot not shown
plot(table(counts), frame = FALSE, main = 'MHB Green Woodpecker Counts',
    type = 'h', lend = 'butt', lwd = 10)
hist(dateso, main = 'Julian Dates (start on week-end around April 15)',
    breaks = 50, col = 'grey')
hist(timeso, main = 'Survey duration (min)', breaks = 50, col = 'grey')
hist(into, main = 'Survey intensity (min / km route length)',
    breaks = 50, col = 'grey')
par(op)

# Standardize survey dates, intensity and times ([ duration)
dates <- standardize(dateso)
int <- standardize(into)
times <- standardize(timeso)

# Put data into 3D arrays and summarize. This works because surveys
#   are grouped within years otherwise be careful!
C <- array(counts, dim=c(267, 3, 14))
DATE <- array(dates, dim=c(267, 3, 14))
DUR <- array(times, dim=c(267, 3, 14))
INT <- array(int, dim=c(267, 3, 14))
mean.C <- apply(C, c(1,3), mean, na.rm = TRUE)      # Mean count per site,year
annual.mean <- apply(mean.C, 2, mean, na.rm = TRUE)
annual.mean2 <- apply(C, 3, mean, na.rm = TRUE)     # Direct 3D average
site.mean <- apply(mean.C, 1, mean, na.rm = TRUE)
nsites.with.data <- apply(!is.na(mean.C), 2, sum, na.rm = TRUE)

cat("N sites with data per year:\n", nsites.with.data, "\n")
# N sites with data per year:
# 265 265 261 263 265 267 267 264 263 263 265 265 264 262 262

nsites.detected <- apply(mean.C > 0, 2, sum, na.rm = TRUE)
cat("N sites with Pecker detected per year:\n", nsites.detected, "\n")
# N sites with Pecker detected per year:
# 49 77 88 93 90 84 103 114 99 127 114 121 136 147

cat("Observed occupancy prob.:\n",
    round(nsites.detected/nsites.with.data,2), "\n")
# Observed occupancy prob.:
# 0.18 0.29 0.34 0.35 0.34 0.31 0.39 0.43 0.38 0.48 0.43 0.46 0.52 0.56

cat("Observed mean count per year:\n", round(annual.mean, 2), "\n")
# Observed mean count per year:
# 0.15 0.33 0.42 0.48 0.45 0.43 0.39 0.54 0.49 0.53 0.43 0.52 0.56 0.58

cat("Observed mean count (site):\n") ; summary(site.mean)
# Observed mean count (site):
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.02381 0.23077 0.44850 0.54304 3.90476

# Scale elevation and forest cover
elevo <- peckers$elev      # Elevation (original)
elev <- standardize(elevo)
foresto <- peckers$forest  # Forest cover (original)
forest <- standardize(foresto)

# Mean-impute DATE and INT
DATE[is.na(DATE)] <- 0     # Innocuous mean imputation
INT[is.na(INT)] <- 0       # Almost innocuous mean imputation
