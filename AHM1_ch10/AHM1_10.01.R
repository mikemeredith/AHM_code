#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

# 10.1. Introduction to the modeling of occurrence, including species distributions
# ---------------------------------------------------------------------------------


alpha <- seq(-10, 10, by = 1)
curve(plogis(-10 +  1 * x), -2, 2, lwd = 3, ylim = c(0,1),
    xlab = "Covariate X", ylab = "Occupancy prob.", frame = FALSE, main = "",
    cex.lab = 1.5, cex.axis = 1.5)
for(i in 2:21){  curve(plogis(alpha[i] +  1 * x), -2, 2, lwd = 3, add = TRUE)  }
(min <- (plogis(-10 +  1 * 2) - plogis(-10 +  1 * -2)) )
(max <- (plogis(0 +  1 * 2) - plogis(0 +  1 * -2)) )

