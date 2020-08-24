#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7 : MODELING FALSE POSITIVES
# ====================================
# Code from proofs dated 2020-08-19

library(unmarked)

# 7.4 Modeling classified false positive detections:
#   'multi-state design' of Miller et al. (2011)
# ================================================

# 7.4.1 Modeling classified false positives in unmarked
# ------------------------------------------------------

# Set parameter values of the simulation
set.seed(129)       # RNG seed
nsites <- 200       # number of sites
nsurveys1 <- 3      # number of occasions with Type 1 data
nsurveys2 <- 4      # number of occasions with Type 2 data
psi <- 0.6          # expected proportion of are occupied
p <- c(0.7,0.5)     # detection prob of method 1 and method 2
fp <- 0.05          # false-positive error probability (p_10)
b <- 0.2            # probability y is recorded as certain

# Simulate the occupancy states and data
z <- rbinom(nsites, 1, psi)
y <- matrix(NA, nrow = nsites, ncol = nsurveys1 + nsurveys2)
for(i in 1:nsites){
  p1 <- p[1]*z[i]                      # certainly detection (method 1)
  p2 <- p[2]*z[i] + fp*(1-z[i])        # uncertainly detection (method 2)
  y[i,1:3] <- rbinom(nsurveys1, 1, p1) # simulate method 1 data
  y[i,4:7] <- rbinom(nsurveys2, 1, p2) # simulate method 2 data
  # Now introduce certain observations:
  pr.certain <- z[i] * y[i,4:7] * b
  y[i, 4:7] <- y[i, 4:7] + rbinom(4, 1, pr.certain)
}

# Make a covariate to distinguish between the two methods
Method <- matrix(c(rep("1", 3), rep("2", 4)), nrow = nsites,
    ncol = nsurveys1 + nsurveys2, byrow = TRUE)

type <- c(nsurveys1, 0, nsurveys2)

summary(umf2 <- unmarkedFrameOccuFP(y = y, obsCovs = list(Method = Method),
    type = type)) # not printed
(m2 <- occuFP(detformula = ~ -1 + Method, FPformula = ~ 1, Bformula = ~ 1,
    stateformula = ~ 1, data = umf2) )

# Occupancy:
# Estimate SE z P(>|z|)
# 0.553 0.148 3.73 0.000192

# Detection:
# Estimate SE z P(>|z|)
# Method1 1.1176 0.1225 9.12 7.37e-20
# Method2 0.0268 0.0891 0.30 7.64e-01

# false positive:
# Estimate SE z P(>|z|)
# -3.27 0.331 -9.88 5.26e-23

# Pcertain:
# Estimate SE z P(>|z|)
# -1.47 0.16 -9.2 3.75e-20

# AIC: 1734.899

# Coefficients on the link (= "beta") scale
coef(m2)
#  psi(Int)    p(Int) p(Method2)    fp(Int)     b(Int)
# 0.5529142 1.1177836 -1.0912711 -3.2700680 -1.4721436

# Coefficients on the probability (="real") scale
pred.df <- data.frame(Method = c("1", "2"))
round(rbind(
    "det" = predict(m2, type = 'det', newdata = pred.df),
    "fp" = predict(m2, type = 'fp', newdata = pred.df[1,,drop=FALSE]),
    "b" = predict(m2, type = 'b', newdata = pred.df[1,,drop=FALSE]),
    "state" = predict(m2, type = 'state', newdata = pred.df[1,,drop=FALSE])),3)
#       Predicted    SE lower upper
# det.1     0.754 0.023 0.706 0.795
# det.2     0.507 0.022 0.463 0.550
# fp        0.037 0.012 0.019 0.068
# b         0.187 0.024 0.144 0.239
# state     0.635 0.034 0.565 0.699


# 7.4.2. A general multi-type model with covariates
# -------------------------------------------------

# Simulation settings
set.seed(2019)           # RNG seed
nsites <- 200            # number of sites
nsurveys <- 7            # number of occasions
habitat <- rnorm(nsites) # Some (continuous) habitat descriptor

# Simulate the occupancy states and data
alpha0 <- 0        # Intercept...
alpha1 <- 1        # ... and slope of psi-habitat regression
psi <- plogis(alpha0 + alpha1*habitat) # Occupancy
z <- rbinom(nsites, 1, psi)            # Latent p/a states
y <- matrix(0,nsites, nsurveys)
p <- c(0.7, 0.5)   # method 2 will have a lower p
b <- 0.5           # probability that a observed positive is determined to be certain
fp <- 0.05         # False-positive prob.

# Simulate data of all 3 types. Note p differs between occ 1-2 and 3-7.
# False positives occur in occasions 3-7 but in occasion 7 there are some
# confirmed positives
for(i in 1:nsites){
  # Normal occupancy data
  y[i, 1:2] <- rbinom(2, 1, p[1]*z[i])
  # False-positives mixed in
  y[i, 3:6] <- rbinom(4, 1, p[2]*z[i] + fp*(1-z[i]))
  # Type 3 observations are occupancy data contaminated with false
  # positives but then we identify some of them as true
  y[i, 7] <- rbinom(1, 1, p[2]*z[i] + fp*(1-z[i]))
}

# Next we set some of the detections to confirmed positives
true.positives <- z==1 & y[,7]==1
confirmed <- (rbinom(nsites, 1, b) == 1) & true.positives
y[confirmed, 7] <- 2

# Make a covariate to distinguish between the two methods
Method <- matrix(c(rep("1", 2), rep("2", 5)), nrow = nsites, ncol = 7,
    byrow = TRUE)

# Type indicates a mix of all 3 data types
type <- c(2, 4, 1)

# Same covariate structure as before
siteCovs <- data.frame(habitat = habitat)
obsCovs <- list(Method = Method)
summary(umf1 <- unmarkedFrameOccuFP(y, siteCovs = siteCovs, obsCovs = obsCovs,
    type = type)) # not shown

# fp starting value should be small (-1 here).
# Note: last parameter in this model is "Pcertain"
( m3 <- occuFP(detformula = ~ -1 + Method, FPformula = ~1, Bformula = ~1,
    stateformula = ~ habitat, data = umf1, starts=c(0, 0, 0, 0, -1, 0)) )

# Occupancy:
# Estimate SE z P(>|z|)
# (Intercept) -0.0813 0.160 -0.509 6.11e-01
# habitat 0.9495 0.198 4.791 1.66e-06
#
# Detection:
# Estimate SE z P(>|z|)
# Method1 0.9751 0.1802 5.412 6.22e-08
# Method2 0.0155 0.0954 0.163 8.71e-01
#
# false positive:
# Estimate SE z P(>|z|)
# -2.86 0.207 -13.8 1.29e-43
#
# Pcertain:
# Estimate SE z P(>|z|)
# 0.0388 0.292 0.133 0.894

(m3b <- occuFP(detformula = ~-1 + Method, FPformula = ~1, Bformula = ~ habitat,
    stateformula = ~ habitat, data = umf1, starts = c(0, 0, 0, 0, - 1, 0, 0)) )
