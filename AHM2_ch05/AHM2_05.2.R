#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5 : MODELING METACOMMUNITY DYNAMICS USING DYNAMIC COMMUNITY MODELS
# ==========================================================================
# Code from proofs dated 2020-08-19

library(AHMbook)

# 5.2 A general simulation function for the DCM model
# ===================================================

# Explicit defaults (need to load AHMbook)
str(data <- simDCM(nspec = 50, nsites = 100, nsurveys = 3, nyears = 10,
    mean.psi1 = 0.4, sig.lpsi1 = 1, mu.beta.lpsi1 = 0, sig.beta.lpsi1 = 0,
    range.mean.phi = c(0.8, 0.8), sig.lphi = 1, mu.beta.lphi = 0,
    sig.beta.lphi = 0, range.mean.gamma = c(0.2, 0.2), sig.lgamma = 1,
    mu.beta.lgamma = 0, sig.beta.lgamma = 0, range.mean.p = c(0.5, 0.5),
    sig.lp = 1, mu.beta.lp = 0, sig.beta.lp = 0, range.beta1.survey = c(0, 0),
    range.beta2.survey = c(0, 0), trend.sd.site = c(0, 0),
    trend.sd.survey = c(0, 0), show.plot = TRUE) )

str(data <- simDCM(nspec = 200))   # More species (looks great)
str(data <- simDCM(nspec = 1))     # A single species (works !)
str(data <- simDCM(nsites = 267))  # More sites
str(data <- simDCM(nsites = 1))    # A single site
str(data <- simDCM(nsurveys = 10)) # More visits
str(data <- simDCM(nyears = 25))   # More years
str(data <- simDCM(nyears = 2))    # Just two years
try(data <- simDCM(nyears = 1))    # A single year ... this crashes

# No species heterogeneity in parameters of initial occupancy
str(data <- simDCM(sig.lpsi1 = 0, sig.beta.lpsi1 = 0))

# No species heterogeneity in parameters of persistence
str(data <- simDCM(sig.lphi = 0, sig.beta.lphi = 0))

# No species heterogeneity in parameters of colonization
str(data <- simDCM(sig.lgamma = 0, sig.beta.lgamma = 0))

# No species heterogeneity in parameters of detection
str(data <- simDCM(sig.lp = 0, sig.beta.lp = 0))

# No annual variation in rates phi, gamma and p
str(data <- simDCM(range.mean.phi = c(0.8, 0.8), range.mean.gamma = c(0.3, 0.3),
    range.mean.p = c(0.6, 0.6)))


set.seed(1)
dat <- simDCM(nspec = 200, nsites = 20, nsurveys = 2, nyears = 10,
    mean.psi1 = 0.1, sig.lpsi1 = 5,
    range.mean.phi = c(0.3, 0.3), sig.lphi = 5,
    range.mean.gamma = c(0.1, 0.1), sig.lgamma = 5,
    range.mean.p = c(0.1, 0.1), sig.lp = 5)

# ** Number of species ever occurring: 177
# ** Number of species ever detected: 121
# ** Average number of years of occurrence: 6.87
# ** Average number of years with detection: 3.695

# Pull out data from one year: you could now feed this data set
# into a static community occupancy model analysis (Chapter 11 in AHM1)
str(yyr1 <- dat$y[,,1,]) # Pull out year 1 as an example
# int [1:20, 1:2, 1:200] 0 0 0 0 0 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 3
# ..$ : chr [1:20] "Site1" "Site2" "Site3" "Site4" ...
# ..$ : chr [1:2] "Survey1" "Survey2"
# ..$ : chr [1:200] "Spec1" "Spec2" "Spec3" "Spec4" ...

# Pull out data for one species: you could now feed this data set
# into a single-species dynocc model analysis (Chapter 4)
str(ysp5 <- dat$y[,,,5]) # Pull out species 5 as an example
# int [1:20, 1:2, 1:10] 0 0 0 0 0 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 3
# ..$ : chr [1:20] "Site1" "Site2" "Site3" "Site4" ...
# ..$ : chr [1:2] "Survey1" "Survey2"
# ..$ : chr [1:10] "Year1" "Year2" "Year3" "Year4" ...
