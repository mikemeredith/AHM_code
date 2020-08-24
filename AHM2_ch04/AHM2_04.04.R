#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

library(AHMbook)

# 4.4 A general data simulation function for dynocc models
# ========================================================

library(AHMbook)
str(data <- simDynocc( # Explicit defaults
    nsites = 250, nyears = 10, nsurveys = 3, year.of.impact = NA,
    mean.psi1 = 0.4, beta.Xpsi1 = 0,
    range.phi = c(0.5, 1), impact.phi = 0, beta.Xphi = 0,
    range.gamma = c(0, 0.5), impact.gamma = 0, beta.Xgamma = 0,
    range.p = c(0.1, 0.9), beta.Xp = 0,
    range.beta1.survey = c(0, 0), range.beta2.survey = c(0, 0),
    trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0),
    trend.sd.site.survey = c(0, 0), show.plot = TRUE))

# All four parameters constant
str(data <- simDynocc(nsites = 250, nyears = 10, nsurveys = 3, mean.psi1 = 0.6,
    range.phi = c(0.7, 0.7), range.gamma = c(0.3, 0.3), range.p = c(0.5, 0.5)))

# Full time-dependence
str(data <- simDynocc(mean.psi1 = 0.6, range.phi = c(0.5, 0.8),
    range.gamma = c(0.1, 0.5), range.p = c(0.1, 0.9)))

# Constant intercepts, but covariates in all parameters
str(data <- simDynocc(mean.psi1 = 0.6, beta.Xpsi1 = 1,
    range.phi = c(0.6, 0.6), beta.Xphi = 2, range.gamma = c(0.3, 0.3),
    beta.Xgamma = 2, range.p = c(0.2, 0.2), beta.Xp = -2) )

# Full time-dependence and and effects of all covariates (incl. season)
str(data <- simDynocc(mean.psi1 = 0.6, beta.Xpsi1 = 1,
    range.phi = c(0.6, 1), beta.Xphi = 2, range.gamma = c(0, 0.2),
    beta.Xgamma = 2, range.p = c(0.1, 0.9), beta.Xp = -2,
    range.beta1.survey = c(2, 10), range.beta2.survey = c(-10, -20)) )

# No detection error (i.e., p = 1)
str(data <- simDynocc(range.p = c(1, 1)) )

# Can do a single site ....
str( data <- simDynocc(nsites = 1) )

# ... but must have at least two years
str(data <- simDynocc(nyears = 2) )

str(data <- simDynocc(nsurveys = 12, mean.psi1 = 0.6,
    range.phi = c(0.6, 0.6), range.gamma = c(0.3, 0.3),
    range.p = c(0.5, 0.5), range.beta1.survey = c(-0.3, 0.4),
    range.beta2.survey = c(0, -0.7)) )

# Add detection heterogeneity at the site level
str(data <- simDynocc(trend.sd.site = c(3, 3)) ) # No time trend
str(data <- simDynocc(trend.sd.site = c(1, 3)) ) # With time trend

# Add detection heterogeneity at the level of the survey
str(data <- simDynocc(trend.sd.survey = c(3, 3)) ) # No time trend
str(data <- simDynocc(trend.sd.survey = c(1, 3)) ) # With time trend

# Add detection heterogeneity at the level of the individual visit
str(data <- simDynocc(trend.sd.site.survey = c(3, 3)) ) # No trend
str(data <- simDynocc(trend.sd.site.survey = c(1, 3)) ) # With trend

str(data <- simDynocc(nsites = 250, nyears = 20, nsurveys = 3,
    year.of.impact = 10, impact.phi = 80, impact.gamma = 50) )
