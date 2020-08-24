#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6 : MULTISTATE OCCUPANCY MODELS
# =======================================
# Code from proofs dated 2020-08-19

library(jagsUI)

# 6.4 Case study: Swiss eagle owls
# ================================

# 6.4.1 Summarizing the Swiss eagle owl data set and preparation for modeling
# ---------------------------------------------------------------------------

# Read in the data
library(AHMbook)
data(SwissEagleOwls)
str(dat <- SwissEagleOwls)

# Append an index for site/year in 'obs'
dat$obs$site.year <- paste(dat$obs$site_name, dat$obs$year, sep = '.')

# Append a binary version of the detection data
dat$obs$dnd <- as.numeric(dat$obs$y > 0) # Simple detection/nondetection
table(dat$obs$dnd)

# Create year covariate
year <- 2007:2016

# Look at frequency of four states in original detection data
table(dat$obs$y)
#    0    1   2   3
# 2414 1977 635 948

# Convert to three states: yms = 'y multi-state'
dat$obs$yms <- dat$obs$y               # Copy
dat$obs$yms[dat$obs$yms == 3] <- 2     # Lump pairs: 0,1,2
dat$obs$yms <- dat$obs$yms + 1         # Renumber observed states: 1,2,3
table(dat$obs$yms)                     # Look at new response data
#    1    2    3
# 2414 1977 1583

# Put detection data and survey dates into a 3D array
(nsites <- nrow(dat$sites))
(nyears <- length(unique(dat$obs$year)))
(nsurveys <- max(as.numeric(names(table(table(dat$obs$site.year))))))
sitelist <- sort(unique(dat$obs$site_name))

# Prepare three empty nsites x nreps x nyears arrays and then fill them all
yms <- date <- y <- array(NA, dim = c(nsites, nsurveys, nyears),
    dimnames = list(sitelist, 1:nsurveys, 1:nyears))
for(i in 1:nsites){
  for(t in 1:nyears){
    sel.site.year <- paste(sitelist[i], t+2006, sep = '.')
    tmp <- dat$obs[dat$obs$site.year == sel.site.year,]
    nr <- nrow(tmp)
    if(nr > 0){
      yms[i,1:nr,t] <- tmp$yms
      date[i,1:nr,t] <- tmp$jdate
      y[i,1:nr,t] <- tmp$dnd
    }
  }
}
sum(yms, na.rm = TRUE) ; sum(dat$obs$yms) # Sum checks ... all OK
sum(date, na.rm = TRUE) ; sum(dat$obs$jdate)
sum(y, na.rm = TRUE) ; sum(dat$obs$dnd)

# Look at parts of these data (first 10 sites, 20 surveys, 3 years)
yms[1:10,,1:3]                         # Multi-state detections
date[1:10,,1:3]                        # Survey dates: these are NOT ordered !
y[1:10,,1:3]                           # Binary detection/nondetections

# Compute number of sites visited at least once per season
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == '-Inf'] <- NA ; tmp[tmp == 0] <- 1
nvisit <- apply(tmp, 2, sum, na.rm = TRUE)

# Compute the proportion of visited sites in all 10 years
propvisit <- nvisit / nsites

# Compute the observed number of occupied sites in all 10 years
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == '-Inf'] <- NA
obsnocc <- apply(tmp, 2, sum, na.rm = TRUE)

# Compute the ratio estimator estimate of population size
nratio <- obsnocc / propvisit

# ~~~~~~~~~~ extra code for figure 6.4 ~~~~~~~~~~~~~~~~~~
# Plot number of visited sites and observed number of occupied sites
ylim <- range(c(obsnocc, nratio))
op <- par(mar = c(5,6,3,3), cex.lab = 1.5, cex.axis = 1.5)
plot(year, nvisit, type = 'b', pch = 1, main = '', xlab = 'Year',
    ylab = 'Number', col = 'black', las = 1, cex = 2, ylim = ylim, frame = FALSE)
points(year, obsnocc, type = 'b', pch = 16, col = 'black', cex = 2)
points(year, nratio, type = 'b', pch = 15, col = 'black', cex = 2)
legend('bottomright', c('Ratio estimator of population size',
    'Number of sites visited', 'Number of sites with Eagle owls detected'),
    pch = c(15, 1, 16), lwd = 1, cex = 1.5, bty = 'n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

table(nvisit.site.year <- apply(y, c(1,3), function(x) sum(!is.na(x))))
#    0   1   2   3   4  5  6  7  8  9 10
# 1457 414 189 146 107 82 57 46 41 25 26
# 11 12 13 14 15 16 17 18 19 20
# 19 15 13 10 10 10  7  6  6 54

# Summarize yms by site and year
tapply(dat$obs$yms, list(dat$obs$site_name, dat$obs$year), max)

# Summarize multi-state detections (yms) by site
# (example only shown for last year)
# Table shows number of observations per observation states for each territory
table(dat$obs$yms[dat$obs$year == 2016], dat$obs$site_name[dat$obs$year == 2016])
#   1 2 3 6 9 10 11 12 17 18 20 22 23 24 26 27 30 31 33
# 1 3 0 3 3 1  1  4  4  0  0  3  0  0  1  1  2  4  5  1
# 2 3 6 0 0 1  0  2  2  1  3  4  3  1  0  3  0  1  2  0
# 3 2 1 1 0 0  0  1  1  4  1  8  0  0  0  1  0  0  0  4
# [output truncated]

# Compute proportion of missing values in the 3D reponse array
(propNA <- sum(is.na(yms)) / prod(dim(yms)))
# [1] 0.8909854

# Compute nsurvey matrix: number of surveys per site/year
nsurveys <- array(1, dim = c(274, 10)) # this is nsites x nyears
for(i in 1:nsites){
  for(t in 1:nyears){
    tmp <- which(!is.na(yms[i,,t]))
    if(length(tmp) > 0){
      nsurveys[i,t] <- max(tmp)
    }
  }
}
head(nsurveys) # Look at matrix (site x year)
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    1    5    2    6    3    7    4   10   17     8
# [2,]    1    1    1    6    1    2    5    3    4     7
# [3,]    3    3   18    1    4   19   20   15   14     4
# [4,]    1    1    1    1    1    1    1    5    1     3
# [5,]    1    1    1    1    1    1    1    1    1     1
# [6,]    1   10    8    8    6    4    4    1    2     2
