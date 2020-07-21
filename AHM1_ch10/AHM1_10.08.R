#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

library(AHMbook)

# 10.8 Goodness of fit
# ====================

simOcc(M = 267, J = 3, mean.occupancy = 0.6, beta1 = 0, beta2 = 1, beta3 = 0,
    mean.detection = 0.3, time.effects = c(0, 0), alpha1 = 1, alpha2 = -1,
    alpha3 = 0, sd.lp = 0, b = 0)

