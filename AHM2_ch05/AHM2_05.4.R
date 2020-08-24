#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5 : MODELING METACOMMUNITY DYNAMICS USING DYNAMIC COMMUNITY MODELS
# ==========================================================================
# Code from proofs dated 2020-08-19

library(AHMbook)

# 5.4 Data formatting for the DCM: fun with multidimensional arrays
# =================================================================

# Get a 4D array
library(AHMbook)
set.seed(123)
tmp <- simDCM(show.plot = FALSE)
BigArray <- tmp$y                # Grab the 4D detection/nondetection array
str(BigArray)                    # 100 sites x 3 visits x 10 years x 50 species

# (1) Shoot holes, i.e., randomly turn 25% of the values into NAs
length(BigArray)                 # 150k values
out <- sample(1:length(BigArray), length(BigArray)/4) # Sample of 25%
BigArray[out] <- NA              # Turn them into NAs
sum(is.na(BigArray))             # we now have 37500 NAs

# (2) Turn data into vector and create indexing factors for the array dimensions
df <- data.frame(
y = c(BigArray),
site = c(slice.index(BigArray, 1)),
visit = c(slice.index(BigArray, 2)),
year = c(slice.index(BigArray, 3)),
species = c(slice.index(BigArray, 4)) )
View(df) ; summary(df)           # Check the 'max.' for each column
str(df)

# (3) Toss out the rows with missing response
df <- df[!is.na(df$y),]          # Select rows with non-missing y only
sum(is.na(df$y))                 # Convince yourself NAs are gone
head(df)                         # Look at first 6 rows in data frame
#   y site visit year species
# 1 0    1     1    1       1
# 3 0    3     1    1       1
# 4 0    4     1    1       1
# 6 0    6     1    1       1
# 7 0    7     1    1       1
# 8 0    8     1    1       1

# (4) Format data in a spreadsheet format into a 4D array
# Determine required dimensions of 4D array
nsite <- length(unique(df$site))     # Number of sites
nvisit <- length(unique(df$visit))   # Number of surveys or visits
nyear <- length(unique(df$year))     # Number of years
nspec <- length(unique(df$species))  # Number of species

# Prepare array and pre-fill array with NAs
BigArray2 <- array(NA, dim = c(nsite, nvisit, nyear, nspec))

# Fill array with the detection/nondetection data
# Loop over all rows in the spreadsheet data and fill them in
# at the right place in the 4D array
for(i in 1:nrow(df)){
  BigArray2[df$site[i], df$visit[i], df$year[i], df$species[i]] <- df$y[i]
}

# Do quick checks ... look good
sum(df$y) ; sum(BigArray2, na.rm = TRUE)               # quick sum check
length(which(is.na(BigArray2)))                        # Same 37500 as before
all.equal(BigArray, BigArray2, check.attributes=FALSE) # BigArray has names

# Get a data set with detections only, no nondetections
str(BigArray <- tmp$y)                # Grab the 4D detection/nondetection array again
df <- data.frame(y = c(BigArray), site = c(slice.index(BigArray, 1)),
    visit = c(slice.index(BigArray, 2)), year = c(slice.index(BigArray, 3)),
    species = c(slice.index(BigArray, 4)) ) # Turn into a spreadsheet format again
df <- df[df$y == 1,]                  # Toss out all nondetection data

# Prepare array by pre-filling it with zeroes instead of NAs
BigArray3 <- array(0, dim = c(100, 3, 10, 50)) # known array dims !

# Fill array with the detection data
for(i in 1:nrow(df)){
  BigArray3[df$site[i], df$visit[i], df$year[i], df$species[i]] <- df$y[i]
}
sum(BigArray3) - nrow(df)             # quick check they're identical

# Reformat array so sites come last
dim(BigArray3)                                # 100 3 10 50
dim(BA.v2 <- aperm(BigArray3, c(2,3,4,1)))    # 3 10 50 100
