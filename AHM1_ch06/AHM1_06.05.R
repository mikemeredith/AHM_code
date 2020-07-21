#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

library(AHMbook)
library(unmarked)

# 6.5 A very general data simulation function for N-mixture models: simNmix
# =========================================================================

# Execute function and inspect results
data <- simNmix()                   # Default arguments
data <- simNmix(show.plot = FALSE)  # Default args, no plots
set.seed(24)
str(data <- simNmix(nsite = 267, nvisit = 3, mean.theta = 1, mean.lam = 2,
    mean.p = 0.6, area = FALSE, beta1.theta = 0, beta2.theta = 0,
    beta3.theta = 0, beta2.lam = 0, beta3.lam = 0, beta4.lam = 0,
    beta3.p = 0, beta5.p = 0, beta6.p = 0, beta.p.survey = 0, beta.p.N = 0,
    sigma.lam = 0, dispersion = 10, sigma.p.site = 0, sigma.p.visit = 0,
    sigma.p.survey = 0, sigma.p.ind = 0, Neg.Bin = FALSE, open.N = FALSE,
    show.plot = TRUE))
    # All default args explicit


str(data <- simNmix())                  # Null data-generating model
str(data <- simNmix(mean.theta = 0.60)) # ZIP with 40% structural zeroes
str(data <- simNmix(sigma.lam = 1))     # Poisson-lognormal (PLN) mixture
str(data <- simNmix(Neg.Bin = TRUE))    # Negative-binomial mixture
str(data <- simNmix(mean.theta = 0.6, sigma.lam = 1))  # Zero-inflated PLN
str(data <- simNmix(mean.theta = 0.6, Neg.Bin = TRUE)) # Zero-infl. NegBin
str(data <- simNmix(mean.p = 1))        # Perfect detection (p = 1)
str(data <- simNmix(mean.theta = 0.6, mean.p = 1))     # ZIP with p = 1
str(data <- simNmix(sigma.lam = 1, mean.p = 1))        # PLN with p = 1


areas <- runif(267, 1, 2)                # Generate vector with site area
str(data <- simNmix(nsite = 267, area = areas)) # Sites with variable area
str(data <- simNmix(nvisit = 1))         # Only one visit
str(data <- simNmix(sigma.p.site = 1))   # Random site effects in p
str(data <- simNmix(sigma.p.visit = 1))  # Random visit (= time) effects in p
str(data <- simNmix(sigma.p.survey = 1)) # Random site-by-visit effects in p
str(data <- simNmix(sigma.p.ind = 1))    # Random individual effects in p
str(data <- simNmix(mean.theta = 0.5, beta1.theta = 1)) # Site cov 1 in suit.
str(data <- simNmix(beta2.lam = 1))      # Site covariate 2 in abundance process
str(data <- simNmix(beta3.p = 1))        # Site covariate 3 in detection process
str(data <- simNmix(beta.p.N = 1))       # Positive density-dep. in p
str(data <- simNmix(beta.p.N = -1))      # Negative density-dep. in p
# Same covariate in suitab. and abund. (see Phillips & Elith, Ecology, 2014 !)
str(data <- simNmix(mean.theta = 0.5, beta2.theta = 1, beta2.lam = -1))
# Same covariate in abundance and detection (see Kéry, Auk, 2008)
str(data <- simNmix(beta3.lam = 1, beta3.p = -1))
# Same covariate in all three levels of model (ouch !)
str(data <- simNmix(mean.theta = 0.5, beta3.theta = 1, beta3.lam = 1, beta3.p = -1))


# Use unmarked to fit some models to these data sets
cov <- data$site.cov
summary(umf <- unmarkedFramePCount(
   y=data$C, siteCovs= data.frame(cov1=cov[,1], cov2=cov[,2], cov3=cov[,3],
      cov4=cov[,4], cov5=cov[,5], cov6=cov[,6], area = data$area),
      obsCovs = list(survey.cov = data$survey.cov)))
summary(fm <- pcount(~1 ~1, umf))
summary(fm <- pcount(~1 ~1, umf, mixture = "ZIP"))
summary(fm <- pcount(~1 ~1, umf, mixture = "NB"))
summary(fm <- pcount(~cov1+cov2+cov3 ~ cov1+cov2+cov3, umf))

