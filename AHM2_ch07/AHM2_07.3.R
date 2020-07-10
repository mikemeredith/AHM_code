#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7 : MODELING FALSE POSITIVES
# ====================================
# Code from proofs dated 2020-07-08

library(unmarked)
library(AHMbook)
library(AICcmodavg)

# 7.3 Joint modeling of type 1 and type 2 data:
#  the 'multi-method design' of Miller et al. (2011)
# ==================================================

# 7.3.1 Multi-method models in unmarked
# -------------------------------------

set.seed(129) # RNG seed
nsites <- 200 # number of sites
nsurveys1 <- 3 # number of occasions with Type 1 data
nsurveys2 <- 4 # number of occasions with Type 2 data
psi <- 0.6 # expected proportion of sites occupied
p <- c(0.7,0.5) # detection prob of method 1 and method 2
fp <- 0.05 # false-positive error probability (p_10)\

# Simulate latent occupancy states and data
z <- rbinom(nsites, 1, psi)
y <- matrix(NA, nrow = nsites, ncol = nsurveys1 + nsurveys2)
for(i in 1:nsites){
  p1 <- p[1]*z[i] # certainly detection (method 1)
  p2 <- p[2]*z[i] + fp*(1-z[i]) # uncertainly detection (method 2)
  y[i,1:3] <- rbinom(nsurveys1, 1, p1) # simulate method 1 data
  y[i,4:7] <- rbinom(nsurveys2, 1, p2) # simulate method 2 data
}

# Make a covariate to distinguish between the two methods
Method <- matrix(c(rep("1", 3), rep("2", 4)), nrow = nsites,
    ncol = nsurveys1 + nsurveys2, byrow = TRUE)

type <- c(nsurveys1, nsurveys2, 0)

summary(umf1 <- unmarkedFrameOccuFP(y = y, obsCovs = list(Method = Method),
    type = type)) # not shown
( m1 <- occuFP(detformula = ~ Method, FPformula = ~1, stateformula = ~ 1,
    data = umf1) )
# Occupancy:
# Estimate SE z P(>|z|)
# 0.544 0.149 3.65 0.000267
# Detection:
# Estimate SE z P(>|z|)
# (Intercept) 0.901 0.119 7.59 3.09e-14
# Method2 -0.967 0.147 -6.59 4.34e-11
# false positive:
# Estimate SE z P(>|z|)
# -2.9 0.282 -10.3 7.31e-25

# Coefficients on the link (= "beta") scale
coef(m1)
# psi(Int) p(Int) p(Method2) fp(Int)
# 0.5442685 0.9008509 -0.9670559 -2.9046690

# Coefficients on the probability (="real") scale
pred.df <- data.frame(Method = c("1","2"))
round(rbind(
    "det" = predict(m1, type = 'det', newdata = pred.df),
    "fp" = predict(m1, type = 'fp', newdata = pred.df[1,,drop=F]),
    "state" = predict(m1, type = 'state', newdata = pred.df[1,,drop=F])),3)
# Predicted SE lower upper
# det.1 0.711 0.024 0.661 0.756
# det.2 0.483 0.022 0.440 0.527
# fp 0.052 0.014 0.031 0.087
# state 0.633 0.035 0.563 0.698


# 7.3.2 Case Study: Eurasian lynx in the Alps
# -------------------------------------------

data(EurasianLynx)
str(lynx <- EurasianLynx)
# 'data.frame': 43332 obs. of 10 variables:
# $ type : Factor w/ 2 levels "certain","uncertain": 1 1 1 1 1 1 1 ...
# $ site.nr : int 1 2 3 4 5 6 7 8 9 10 ...
# $ y.1 : int NA NA NA NA NA NA NA NA NA NA ...
# $ y.2 : int 0 0 0 0 0 0 0 0 0 0 ...
# $ y.3 : int 0 0 0 0 0 0 0 0 0 0 ...
# $ Year : int 1994 1994 1994 1994 1994 1994 1994 1994 1994 1994 ...
# $ Cntry : Factor w/ 2 levels "Italy","Switzerland": 1 1 1 1 1 1 1 1 1 1
# $ forest : num 78.3 34.1 40.8 68.2 73.5 ...
# $ xcoord : num 4110 4120 4130 4110 4120 4130 4140 4150 4160 4170 ...
# $ ycoord : num 2300 2300 2300 2310 2310 2310 2310 2310 2310 2310 ...

# Add the columns we need for analysis in unmarked
lynx$occ.1 <- 1
lynx$occ.2 <- 2
lynx$occ.3 <- 3
lynx$sYear <- standardize(lynx$Year)

# Extract the type 1 and type 2 data separately and bind them together
lynx1 <- lynx[lynx[,"type"] == "certain", ]
lynx2 <- lynx[lynx[,"type"] == "uncertain", ]
lynx <- cbind(lynx1[,c(2,3:5)], lynx2[,3:5] )
colnames(lynx) <- c("site.nr", "y.1", "y.2", "y.3", "y.4", "y.5", "y.6")
occ <- cbind(lynx1[, c("occ.1", "occ.2", "occ.3")],
lynx2[, c("occ.1", "occ.2", "occ.3")])
colnames(occ) <- c("occ.1", "occ.2", "occ.3", "occ.4", "occ.5", "occ.6")
lynx <- cbind(lynx,lynx1[, c("Year", "sYear", "Cntry")])

# Make the false-positive unmarkedFrame. Be sure to indicate type!
y <- lynx[,paste0("y.", 1:6)]
siteCovs <- lynx[, c("sYear", "Year", "Cntry")]
obsCovs <- list(occ = occ)
summary(lynx.umf <- unmarkedFrameOccuFP(y = y,
    siteCovs = siteCovs, obsCovs = obsCovs, type = c(3, 3, 0)) )

# Models to investigate trends in recording (takes approx. 15 mins)
cand.mods <- list(
    "p(c)fp(c)psi(c*t)" = occuFP(~Cntry, ~Cntry, ~1, ~Cntry*sYear, data = lynx.umf),
    "p(c*t)fp(c)psi(c*t)" = occuFP(~Cntry*sYear, ~Cntry, ~1, ~Cntry*sYear,
        data = lynx.umf),
    "p(c)fp(c*t)psi(c*t)" = occuFP(~Cntry, ~Cntry*sYear, ~1, ~Cntry*sYear,
        data = lynx.umf),
    "p(c*t)fp(c*t)psi(c*t)" = occuFP(~Cntry*sYear, ~Cntry*sYear, ~1, ~Cntry*sYear,
        data = lynx.umf),
    "p(c+t)fp(c+t)psi(c*t)" = occuFP(~Cntry + sYear, ~Cntry + sYear, ~1,
        ~Cntry*sYear, data = lynx.umf),
    "p(c+t)fp(c+t)psi(c+t)" = occuFP(~Cntry + sYear, ~Cntry + sYear, ~1,
        ~Cntry + sYear, data = lynx.umf),
    "p(c)fp(c)psi(c)" = occuFP(~Cntry, ~Cntry, ~1, ~Cntry, data = lynx.umf),
    "p(c*t)fp(c)psi(c)" = occuFP(~Cntry*sYear, ~Cntry, ~1, ~Cntry, data = lynx.umf),
    "p(c)fp(c*t)psi(c)" = occuFP(~Cntry, ~Cntry*sYear, ~1, ~Cntry, data = lynx.umf),
    "p(c*t)fp(c*t)psi(c)" = occuFP(~Cntry*sYear, ~Cntry*sYear, ~1, ~Cntry,
        data = lynx.umf),
    "p(c+t)fp(c+t)psi(c)" = occuFP(~Cntry + sYear, ~Cntry + sYear, ~1, ~Cntry,
        data = lynx.umf),
    "p(c)fp(c)psi(t)" = occuFP(~Cntry, ~Cntry, ~1, ~sYear, data = lynx.umf),
    "p(c*t)fp(c)psi(t)" = occuFP(~Cntry*sYear, ~Cntry, ~1, ~sYear, data = lynx.umf),
    "p(c)fp(c*t)psi(t)" = occuFP(~Cntry, ~Cntry*sYear, ~1, ~sYear, data = lynx.umf),
    "p(c*t)fp(c*t)psi(t)" = occuFP(~Cntry*sYear, ~Cntry*sYear, ~1, ~sYear,
        data = lynx.umf),
    "p(c+t)fp(c+t)psi(t)" = occuFP(~Cntry + sYear, ~Cntry + sYear, ~1, ~sYear,
        data = lynx.umf))

(modTab <- modSelFP(cand.mods))
# nPars AIC dAIC AICwt cuWt
# p(c+t)fp(c+t)psi(c*t) 10 24386.76 0.00 0.87 0.87
# p(c*t)fp(c*t)psi(c*t) 12 24390.48 3.72 0.13 1.00
# p(c)fp(c*t)psi(c*t) 10 24411.85 25.09 0.00 1.00
# [... output truncated ...]

# Select and summarize the AIC-top model
( topmod <- cand.mods$`p(c+t)fp(c+t)psi(c*t)` )
# Occupancy:
# Estimate SE z P(>|z|)
# (Intercept) -3.649 0.0935 -39.01 0.00e+00
# CntrySwitzerland 2.508 0.0991 25.31 2.61e-141
# sYear -0.314 0.0766 -4.10 4.22e-05
# CntrySwitzerland:sYear 0.584 0.0807 7.24 4.49e-13
# Detection:
# Estimate SE z P(>|z|)
# (Intercept) -1.319 0.0904 -14.60 3.01e-48
# CntrySwitzerland 0.176 0.0953 1.85 6.48e-02
# sYear 0.154 0.0310 4.98 6.24e-07
# false positive:
# Estimate SE z P(>|z|)
# (Intercept) -6.546 0.253 -25.868 1.54e-147
# CntrySwitzerland -6.744 9.342 -0.722 4.70e-01
# sYear -0.978 0.200 -4.885 1.04e-06


# ~~~ figure 4.7 ~~~~~~~

