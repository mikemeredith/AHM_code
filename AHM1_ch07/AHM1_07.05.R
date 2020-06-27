#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7. Modeling abundance using multinomial N-mixture models
# =========================================================================

library(AHMbook)
library(unmarked)

# 7.5 Example 1: Bird point counts based on removal sampling
# ==========================================================


data(ovendata)
ovendata.list$data[11:20,]   # Look at a snippet of data set

apply(ovendata.list$data,2,sum)  # Removals in occasion 1-4


# 7.5.1 Setting up the data for analysis
# ------------------------------------------------------------------------
library(unmarked)
data(ovendata)
ovenFrame <- unmarkedFrameMPois(y = ovendata.list$data,
    siteCovs = as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")


# 7.5.2 Fitting models using function multinomPois
# ------------------------------------------------------------------------
# Fit models: multinomPois order of formulas: detection, abundance
fm0 <- multinomPois(~ 1 ~ 1, ovenFrame)
fm1 <- multinomPois(~ 1 ~ ufc, ovenFrame)
fm2 <- multinomPois(~ 1 ~ trba, ovenFrame)
fm3 <- multinomPois(~ 1 ~ ufc + trba, ovenFrame)
fm4 <- multinomPois(~ 1 ~ ufc + trba + ufc:trba, ovenFrame)
fm5 <- multinomPois(~ ufc ~ ufc + trba, ovenFrame)
fm6 <- multinomPois(~ ufc ~ ufc + trba + ufc:trba, ovenFrame)

# Rank models by AIC
ms <- fitList(
    "lam(.)p(.)"                                = fm0,
    "lam(ufc)p(.)"                              = fm1,
    "lam(trba)p(.)"                             = fm2,
    "lam(ufc+trba)p(.)"                         = fm3,
    "lam(ufc+trba+ufc:trba)p(.)"                = fm4,
    "lam(ufc+trba)p(ufc)"                       = fm5,
    "lam(ufc+trba+ufc:trba)p(ufc)"              = fm6)

(ms1 <- modSel(ms))

# Table with everything you could possibly need
coef(ms1)[,1:4]  # Only first 4 columns shown

output <- as(ms1, "data.frame")

# 7.5.3 Fitting models using function gmultmix
# ------------------------------------------------------------------------
ovenFrame <- unmarkedFrameGMM(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
       numPrimary=1,type = "removal")

fm0 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
       data=ovenFrame)

# Fit Poisson models
fm1 <- gmultmix(~ ufc, ~ 1, ~  1, data = ovenFrame)
fm2 <- gmultmix(~ trba, ~ 1, ~ 1, data = ovenFrame)
fm3 <- gmultmix(~ ufc + trba, ~ 1, ~ 1, data = ovenFrame)
fm4 <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ 1, data = ovenFrame)
# Maybe p also depends on understory foliage?
fm5 <- gmultmix(~ ufc + trba, ~ 1, ~ ufc, data = ovenFrame)
fm6 <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ ufc, data = ovenFrame)

# Fit analogous NegBin models
fm0nb <- gmultmix(~ 1, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm1nb <- gmultmix(~ ufc, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm2nb <- gmultmix(~ trba, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm3nb <- gmultmix(~ ufc + trba , ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm4nb <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ 1, mixture = "NB",
       data = ovenFrame)
# maybe p also depends on understory foliage?
fm5nb <- gmultmix(~ ufc + trba, ~ 1, ~ ufc, mixture = "NB",
       data = ovenFrame)
fm6nb <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ ufc, mixture = "NB",
       data = ovenFrame)

# Rank models by AIC
gms <- fitList(
    "lam(.)p(.)"                                = fm0,
    "lam(ufc)p(.)"                              = fm1,
    "lam(trba)p(.)"                             = fm2,
    "lam(ufc+trba)p(.)"                         = fm3,
    "lam(ufc+trba+ufc:trba)p(.)"                = fm4,
    "lam(ufc+trba)p(ufc)"                       = fm5,
    "lam(ufc+trba+ufc:trba)p(ufc)"              = fm6,
    "NB,lam(.)p(.)"                             = fm0nb,
    "NB,lam(ufc)p(.)"                           = fm1nb,
    "NB,lam(trba)p(.)"                          = fm2nb,
    "NB,lam(ufc+trba)p(.)"                      = fm3nb,
    "NB,lam(ufc+trba+ufc:trba)p(.)"             = fm4nb,
    "NB,lam(ufc+trba)p(ufc)"                    = fm5nb,
    "NB,lam(ufc+trba+ufc:trba)p(ufc)"           = fm6nb)

(gms1 <- modSel(gms))

# Table with everything you could possibly need
output <- as(gms1, "data.frame")

# Summary results
gms1

fm2nb


print(coef(gms1), digits = 2)


# 7.5.4 Assessing model fit in unmarked
# ------------------------------------------------------------------------
set.seed(2015)
(gof <- parboot(fm2, fitstats, nsim = 1000, report = 1))

