#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 11. Hierarchical models for communities
# =========================================================================

library(AHMbook)

# 11.3 Metacommunity data from the Swiss breeding bird survey MHB
# ===============================================================

## Code modified to use the built-in data set MHB2014 instead of the file "MHB_2014.csv"
data(MHB2014)
?MHB2014
str(MHB2014)
# NB some of the data preprocessing on p.644 has already been done.

# Check the detection data in 3D array MHB2014$count: site x rep x species
( nsite <- nrow(MHB2014$sites) )    # number of sites in Swiss MHB
nrep <- 3                           # maximum number of replicate surveys per season
( nspec <- nrow(MHB2014$species) )  # 158 species occur in the 2014 data
dim(MHB2014$count) == c(nsite, nrep, nspec) # check

# Create the detection/nondetection (1/0) array
y <- MHB2014$count ; y[y > 1] <- 1  ## 'Y' replaced with 'y'
str(y)

# Check data for one species, here chaffinch, and pull them out from 3D array
(tmp <- y[, , "Common Chaffinch"])

# Frequency distribution of number of surveys actually carried out per site in 2014
# NB MHB2014$sites$nsurvey gives the number of surveys *planned*.
table(nsurveys <- apply(!is.na(y[,,1]), 1, sum))

# Which site has all NA data in 2014 ?
(NAsites <- which(nsurveys == 0) )

# Observed number of occupied sites
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
# For the 'all NA' site, max returns -Inf with a warning
tmp[tmp == -Inf] <- NA         # Change -Inf to NA
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))

# Plot species 'occurrence frequency' distribution (not shown)
plot(sort(obs.occ), xlab = "Species number", ylab = "Number of quads with detections")

# Drop data from species that were not observed in 2014
toss.out <- which(obs.occ == 0)
y <- y[,,-toss.out]
obs.occ <- obs.occ[-toss.out]

# Redefine nspec as the number of species observed in 2014: 145
( nspec <- dim(y)[3] )

str(y)

# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum))     # Compute and print sorted species counts

plot(table(C), xlim = c(0, 60), xlab = "Observed number of species",
    ylab = "Number of quadrats", frame = FALSE)
abline(v = mean(C, na.rm = TRUE), col = "blue", lwd = 3)

