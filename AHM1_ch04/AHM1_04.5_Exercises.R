#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 4. Introduction to data simulation
# =========================================================================

library(AHMbook)

# 4.5 Exercises
# =============

# Generate data set and fit Poisson and Bernoulli models
data <- data.fn(M = 267, J = 1, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 1, show.plot = FALSE)
summary(fmPois <- glm(C ~ elev*forest, family = poisson, data = data))
presence <- ifelse(data$C > 0, 1, 0)
summary(fmBern <- glm(presence ~ elev*forest, family = binomial(link = "cloglog"), data = data))

# Compare Poisson and Bernoulli estimates with truth
print(cbind(Truth = c(data$beta0, data$beta1, data$beta2, data$beta3), summary(fmPois)$coef[,1:2], summary(fmBern)$coef[,1:2]),3)

