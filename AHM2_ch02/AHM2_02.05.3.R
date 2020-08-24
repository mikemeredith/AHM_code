#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 2 : MODELING POPULATION DYNAMICS WITH COUNT DATA
# ========================================================
# Code from proofs dated 2020-08-18

library(jagsUI)
library(unmarked)
library(AHMbook)

# ~~~~~~~~~~ recreate the data set (see 2.5.1) ~~~~~~~~~~
set.seed(2017, kind = "L'Ecuyer")
str(data <- simDM0(nsites = 50, nsurveys = 3, nyears = 5, lambda = 4,
    phi = 0.8, gamma = 1.5, p = 0.7))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.5 Dynamic N-mixture model of Dail-Madsen
# ==========================================

# 2.5.3 The unmarked function pcountOpen
# --------------------------------------

# Prepare data
summary(umf <- unmarkedFramePCO(y = data$yy, numPrimary = data$nyears))

# Fit model, backtransform and compare with truth
# Dynamics = constant (as in data simulation)
(fm1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K = max(data$yy) + 100,
    dynamics = "constant", control = list(trace=TRUE, REPORT=1)) )

# Abundance:
# Estimate SE z P(>|z|)
# 1.4 0.0755 18.5 1.66e-76

# Recruitment:
# Estimate SE z P(>|z|)
# 0.104 0.144 0.721 0.471

# Apparent Survival:
# Estimate SE z P(>|z|)
# 1.87 0.276 6.75 1.44e-11

# Detection:
# Estimate SE z P(>|z|)
# 0.87 0.0825 10.5 5.93e-26

# AIC: 2466.449

# Back-transformation of parameters (full output not printed)
(lam <- coef(backTransform(fm1, "lambda"))) # or
(om <- plogis(coef(fm1, type="omega"))) # Apparent survival !
(gam <- exp(coef(fm1, type="gamma")))
(p <- plogis(coef(fm1, type="det")))
# lam(Int)
# 4.046917
# omega(Int)
# 0.8660311
# gamConst(Int)
# 1.109454
# p(Int)
# 0.7046613

# Dynamics = "autoreg"
# (not consistent with data generating model)
(fm2 <- pcountOpen(lam = ~1, gam = ~1, omega = ~1, p = ~1, data = umf,
    dynamics = "autoreg", K = max(data$yy) + 100,
    control = list(trace=TRUE, REPORT = 1)))

# Abundance (log-scale):
# Estimate SE z P(>|z|)
# 1.42 0.0757 18.7 3.3e-78

# Recruitment (log-scale):
# Estimate SE z P(>|z|)
# -1.48 0.157 -9.4 5.4e-21

# Apparent Survival (logit-scale):
# Estimate SE z P(>|z|)
# 1.87 0.297 6.29 3.1e-10

# Detection (logit-scale):
# Estimate SE z P(>|z|)
# 0.844 0.0842 10 1.19e-23

# AIC: 2490.501

# Backtransformation of parameters (output omitted)
(lam <- exp(coef(fm2, type = "lambda")))
(om <- plogis(coef(fm2, type = "omega")))
(gam <- exp(coef(fm2, type = "gamma")))
(p <- plogis(coef(fm2, type = "det")))
