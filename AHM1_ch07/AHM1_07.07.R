#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7. Modeling abundance using multinomial N-mixture models
# =========================================================================

library(unmarked)

# 7.7 Building custom multinomial models in unmarked
# ==================================================


# Removal model: capture probs for 5 sites, with 3 removal periods
(pRem <- matrix(0.5, nrow=5, ncol=3))

removalPiFun(pRem)   # Multinomial cell probabilities for each site

# Double observer model: capture probs for 5 sites, with 2 observers
(pDouble <- matrix(0.5, 5, 2))

doublePiFun(pDouble)  # Multinomial cell probabilities for each site


instRemPiFun <- function(p){
   M <- nrow(p)
   J <- ncol(p)
   pi <- matrix(NA, M, J)
   p[,1] <- pi[,1] <- 1 - (1 - p[,1])^2
   p[,2] <- 1 - (1 - p[,2])^3
   p[,3] <- 1 - (1 - p[,3])^5
   for(i in 2:J) {
      pi[,i] <- pi[, i - 1]/p[, i - 1] * (1 - p[, i - 1]) * p[, i]
   }
   return(pi)
}


instRemPiFun(pRem)

o2y <- matrix(1, 2, 3)
o2y

