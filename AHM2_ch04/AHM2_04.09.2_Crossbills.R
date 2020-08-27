#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 3.5 hrs
# Run time with the full number of iterations: 38.3 hrs

library(AHMbook)
library(unmarked)
library(AICcmodavg)

# ~~~ load crossbill data from 4.9.1 ~~~~~~~~~~
source("AHM2_04.09.1_Crossbills.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 4.9 Analysis and mapping of crossbill distribution and range dynamics in Switzerland
# ====================================================================================

# 4.9.2 Fitting a large, “global” dynamic occupancy model in unmarked
#   in stepwise manner as a guard against local optima
# -------------------------------------------------------------------

# ~~~~~ bonus code, not in the printed book, dated 2019-01-04 ~~~~~~~~~~~
# Function colext fits dynamic site-occupancy model of MacKenzie et al. (2003)
# See ?colext and vignette("colext")

### Investigate temporal variation first by fitting a Null model, the
#     fully time-dependent model and the intermediate models
# A Null model with constant params throughout
system.time(
   (fm0 <- colext(psiformula = ~ 1, gammaformula = ~ 1,
      epsilonformula = ~ 1, pformula = ~ 1, umf,
      control = list(trace=TRUE, REPORT=5)))
   )  # 2 secs
summary(fm0)
cbind(Npar = length(fm0@opt$par), NLL = fm0@negLogLike, AIC = fm0@AIC)

# Only colonisation varies
system.time(
   (fm1 <- colext(psiformula = ~ 1, gammaformula = ~ year-1,
      epsilonformula = ~ 1, pformula = ~ 1, umf,
      control=list(trace=TRUE, REPORT=5)))
   )  # 14 secs
summary(fm1)
cbind(Npar = length(fm1@opt$par), NLL = fm1@negLogLike, AIC = fm1@AIC)

# Only extinction varies
system.time(
   (fm2 <- colext(psiformula = ~ 1, gammaformula = ~ 1,
        epsilonformula = ~ year-1, pformula = ~ 1, umf,
        control=list(trace=TRUE, REPORT=5)))
   )  # 12 secs
summary(fm2)
cbind(Npar = length(fm2@opt$par), NLL = fm2@negLogLike, AIC = fm2@AIC)

# Only detection varies
system.time(
   (fm3 <- colext(psiformula = ~ 1, gammaformula = ~ 1,
        epsilonformula = ~ 1, pformula = ~ year-1, umf,
        control=list(trace=TRUE, REPORT=5)))
   )  # 15 secs
summary(fm3)
cbind(Npar = length(fm3@opt$par), NLL = fm3@negLogLike, AIC = fm3@AIC)

# Fully time-dependent model
system.time(
   (fm4 <- colext(~ 1, ~ year-1, ~year-1, ~year-1, umf,
        control=list(trace=TRUE, REPORT=5)))
   )  # 70 secs
summary(fm4)
cbind(Npar = length(fm4@opt$par), NLL = fm4@negLogLike, AIC = fm4@AIC)
cbind(AIC.constant.model = fm0@AIC, AIC.time_dep.model = fm4@AIC)

### Select between time constancy and different forms of time-dependence
models <- fitList(
  'psi(.)gam(.)eps(.)p(.)' = fm0,
  'psi(.)gam(Y)eps(.)p(.)' = fm1,
  'psi(.)gam(.)eps(Y)p(.)' = fm2,
  'psi(.)gam(.)eps(.)p(Y)' = fm3,
  'psi(.)gam(Y)eps(Y)p(Y)' = fm4)
(ms <- modSel(models))

                       # nPars     AIC  delta   AICwt cumltvWt
# psi(.)gam(Y)eps(Y)p(Y)    35 6806.25   0.00 1.0e+00     1.00
# psi(.)gam(.)eps(.)p(Y)    15 6837.96  31.71 1.3e-07     1.00
# psi(.)gam(Y)eps(.)p(.)    14 6907.15 100.90 1.2e-22     1.00
# psi(.)gam(.)eps(.)p(.)     4 6938.86 132.61 1.6e-29     1.00
# psi(.)gam(.)eps(Y)p(.)    14 6944.23 137.98 1.1e-30     1.00

# This is overwhelming evidence for (inter-)annual variation in all parameters that can vary at all in this model: colonization, extinction and detection probability.

# The most complex model that we want to fit is this (written in the usual R model definition language, in which we could also fit it in unmarked, and will below):

# psi1((elev + elev2) * (forest + forest2))
# gamma(year + (elev + elev2) * (forest + forest2))
# eps(year + (elev + elev2) * (forest + forest2))
# p(year + (elev + elev2) * (forest + forest2) + date + date2 + date:elev + date:elev2 + date2:elev + date2:elev2)

# Build up model for initial occupancy (psi1)
system.time(
   (fm01<- colext(~ elev,
    ~ (year-1),
    ~ (year-1),
    ~ (year-1), umf,
control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 20 secs
summary(fm01)

system.time(
   (fm02<- colext(~ elev + forest,
    ~ (year-1),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)
summary(fm02)

system.time(
   (fm03<- colext(~ elev + forest + I(elev^2),
    ~ (year-1),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 30 secs
summary(fm03)

system.time(
   (fm04<- colext(~ elev + forest + I(elev^2) + I(forest^2) ,
    ~ (year-1),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 35 secs
summary(fm04)

system.time(
   (fm05<- colext(~ elev + forest + I(elev^2) + I(forest^2) +
       elev:forest ,
    ~ (year-1),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 45 secs
summary(fm05)

system.time(
   (fm06<- colext(~ elev + forest + I(elev^2) + I(forest^2) +
       elev:forest + I(elev^2):forest ,
    ~ (year-1),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 60 secs
summary(fm06)

system.time(
   (fm07<- colext(~ elev + forest + I(elev^2) + I(forest^2) +
       elev:forest + I(elev^2):forest + elev:I(forest^2) ,
    ~ (year-1),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 33 secs
summary(fm07)

system.time(
   (fm08<- colext(~ elev + forest + I(elev^2) + I(forest^2) +
       elev:forest + I(elev^2):forest + elev:I(forest^2) +
       I(elev^2):I(forest^2) ,
    ~ (year-1),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 40 secs
summary(fm08)

# Check number of parameters and NLL for this series of models, fig 4.13 top left
plot(1:8, cbind(fm01@negLogLike, fm02@negLogLike, fm03@negLogLike, fm04@negLogLike,
    fm05@negLogLike, fm06@negLogLike, fm07@negLogLike, fm08@negLogLike),
    type = "b", pch = 16, ylab = "NLL", xlab = "Model Number",
    main = "Negative log-likelihood for series of models for initial occupancy (psi1)")

# Build up model for colonization rate (gamma)
system.time(
   (fm09<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev,
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 70 secs
summary(fm09)

system.time(
   (fm10<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest,
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 35 secs
summary(fm10)

system.time(
   (fm11<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest +I(elev^2),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 45 secs
summary(fm11)

system.time(
   (fm12<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest +I(elev^2) +I(forest^2),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 150 secs
summary(fm12)

system.time(
   (fm13<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest +I(elev^2) +I(forest^2) +
          elev:forest,
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 55 secs
summary(fm13)

system.time(
   (fm14<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest +I(elev^2) +I(forest^2) +
          elev:forest + I(elev^2):forest,
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 95 secs
summary(fm14)

system.time(
   (fm15<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest +I(elev^2) +I(forest^2) +
          elev:forest + I(elev^2):forest + elev:I(forest^2),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 95 secs
summary(fm15)

system.time(
   (fm16<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest +I(elev^2) +I(forest^2) +
          elev:forest + I(elev^2):forest + elev:I(forest^2) +
          I(elev^2):I(forest^2),
    ~ (year-1),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 45 secs
summary(fm16)

# Check number of parameters and NLL for this series of models, fig 4.13 top right
plot(9:16, cbind(fm09@negLogLike, fm10@negLogLike, fm11@negLogLike, fm12@negLogLike,
    fm13@negLogLike, fm14@negLogLike, fm15@negLogLike, fm16@negLogLike),
    type = "b", pch = 16, ylab = "NLL", xlab = "Model Number",
    main = "Negative log-likelihood for series of models for colonization (gamma)")

# This doesn't look good and shows clearly that the function minimisation algorithm for at least the models 12, 13 and 15 was caught in a place that does not represent the global optimum. Yet, we have jacked up maxit to 500 and function optim records successful convergence for all these model runs.

# Build up model for extinction rate (eps)
system.time(
   (fm17<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1)+ elev,
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 45 secs
summary(fm17)

system.time(
   (fm18<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1)+ elev + forest,
    ~ (year-1), umf,
control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 75 secs
summary(fm18)

system.time(
   (fm19<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1)+ elev + forest + I(elev^2),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 115 secs
summary(fm19)

system.time(
   (fm20<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1)+ elev + forest + I(elev^2) + I(forest^2),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 75 secs
summary(fm20)

system.time(
   (fm21<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1)+ elev + forest + I(elev^2) + I(forest^2) +
          elev:forest,
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 70 secs
summary(fm21)

system.time(
   (fm22<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) ,
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1)+ elev + forest + I(elev^2) + I(forest^2) +
          elev:forest + I(elev^2):forest,
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 45 secs
summary(fm22)

system.time(
   (fm23<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1)+ elev + forest + I(elev^2) + I(forest^2) +
          elev:forest + I(elev^2):forest + elev:I(forest^2),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 80 secs
summary(fm23)

system.time(
   (fm24<- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1)+ elev + forest + I(elev^2) + I(forest^2) +
          elev:forest + I(elev^2):forest + elev:I(forest^2) +
          I(elev^2):I(forest^2),
    ~ (year-1), umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 200 secs
summary(fm24)

# Check number of parameters and NLL for this series of models, fig 4.13 bottom left
plot(17:24, cbind(fm17@negLogLike, fm18@negLogLike, fm19@negLogLike, fm20@negLogLike,
    fm21@negLogLike, fm22@negLogLike, fm23@negLogLike, fm24@negLogLike),
    type = "b", pch = 16, ylab = "NLL", xlab = "Model Number",
    main = "Negative log-likelihood for series of models for extinction (eps)")

# Hence, at least the models 19 and 23 are the result of the algorithm finding a local rather than the global optimum. Again, optim says that these model runs have converged (and we set maxit = 500).

# Build up model for detection probability (p)
system.time(
   (fm25 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev, umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 50 secs
summary(fm25)


system.time(
   (fm26 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest, umf,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 60 secs
summary(fm26)

system.time(
   (fm27 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest + I(elev^2), umf,
control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 65 secs
summary(fm27)

system.time(
   (fm28 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest + I(elev^2) + I(forest^2),
umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 60 secs
summary(fm28)

system.time(
   (fm29 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest + I(elev^2) + I(forest^2) +
      elev:forest,
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 70 secs
summary(fm29)

system.time(
   (fm30 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest + I(elev^2) + I(forest^2) +
          elev:forest + I(elev^2):forest,
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 85 secs
summary(fm30)

system.time(
   (fm31 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest + I(elev^2) + I(forest^2) +
          elev:forest + I(elev^2):forest + elev:I(forest^2),
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 90 secs
summary(fm31)

system.time(
   (fm32 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + elev + forest + I(elev^2) + I(forest^2) +
          elev:forest + I(elev^2):forest + elev:I(forest^2) +
I(elev^2):I(forest^2),
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 90 secs
summary(fm32)

system.time(
   (fm33 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) + date,
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 90 secs
summary(fm33)

system.time(
   (fm34 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2),
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 80 secs
summary(fm34)

system.time(
   (fm35 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
          date + I(date^2) + date:elev,
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 85 secs
summary(fm35)

system.time(
   (fm36 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
          date + I(date^2) + date:elev + date:I(elev^2),
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 85 secs
summary(fm36)

system.time(
   (fm37 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
          date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev,
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)  # 90 secs
summary(fm37)

system.time(
   (fm38 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
          date + I(date^2) + date:elev + date:I(elev^2) +
          I(date^2):elev + I(date^2):I(elev^2),
      umf, control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
)
summary(fm38)


# Check number of parameters and NLL for this series of models, fig 4.13 bottom right
plot(25:38, cbind(fm25@negLogLike, fm26@negLogLike, fm27@negLogLike, fm28@negLogLike,
    fm29@negLogLike, fm30@negLogLike, fm31@negLogLike, fm32@negLogLike, fm33@negLogLike,
    fm34@negLogLike, fm35@negLogLike, fm36@negLogLike, fm37@negLogLike, fm38@negLogLike),
    type = "b", pch = 16, ylab = "NLL", xlab = "Model Number",
    main = "NLL for series of models for detection (p)")

# Now, the following models have not converged at their global minimum: 33 and 38, though, again, optim says that these model runs have converged.

# Look at NLLs for all 4 series (Fig 4.13)
oldpar <- par(mfrow = c(2,2), mar = c(5,4,4,2))

# Check number of parameters and NLL for this series of models
plot(1:8, cbind(fm01@negLogLike, fm02@negLogLike, fm03@negLogLike, fm04@negLogLike,
    fm05@negLogLike, fm06@negLogLike, fm07@negLogLike, fm08@negLogLike),
    type = "b", pch = 16, ylab = "NLL", xlab = "Model Number",
    main = "NLL for series of models for initial occupancy (psi1)")

# Check number of parameters and NLL for this series of models
plot(9:16, cbind(fm09@negLogLike, fm10@negLogLike, fm11@negLogLike, fm12@negLogLike,
    fm13@negLogLike, fm14@negLogLike, fm15@negLogLike, fm16@negLogLike),
    type = "b", pch = 16, ylab = "NLL", xlab = "Model Number",
    main = "NLL for series of models for colonization (gamma)")

# Check number of parameters and NLL for this series of models
plot(17:24, cbind(fm17@negLogLike, fm18@negLogLike, fm19@negLogLike, fm20@negLogLike,
    fm21@negLogLike, fm22@negLogLike, fm23@negLogLike, fm24@negLogLike),
    type = "b", pch = 16, ylab = "NLL", xlab = "Model Number",
    main = "NLL for series of models for extinction (eps)")

# Check number of parameters and NLL for this series of models
plot(25:38, cbind(fm25@negLogLike, fm26@negLogLike, fm27@negLogLike, fm28@negLogLike,
    fm29@negLogLike, fm30@negLogLike, fm31@negLogLike, fm32@negLogLike, fm33@negLogLike,
    fm34@negLogLike, fm35@negLogLike, fm36@negLogLike, fm37@negLogLike, fm38@negLogLike),
    type = "b", pch = 16, ylab = "NLL", xlab = "Model Number",
    main = "NLL for series of models for detection (p)")
par(oldpar)

# We then tried to find the global optimum for model 38 by initializing a search at the solutions of the neighbouring model 37 and also by starting the search repeatedly from various sets of random starting values in the vicinity of a previous set of solutions of model 38.

# Set up list for various runs of model fm38
fm38.list <- list()    # list to save results from 10 runs of fm38

# (1) Free search without provision of inits
summary(
  fm38.list[[1]] <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev + I(date^2):I(elev^2),
    umf, control=list(trace=TRUE, REPORT=5), se = FALSE))


# (2) Use estimates of model 37 as starting values
inits <- c(coef(fm37), 0) # extra parameter initialised at 0
summary(
   fm38.list[[2]] <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
   ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
   ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
   ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
       date + I(date^2) + date:elev + date:I(elev^2) +
       I(date^2):elev + I(date^2):I(elev^2), umf, starts = inits,
     control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))


# (3-10) Disturb the previous solutions a little bit
#    and truncate probable boundary estimates
inits <- function() {  # Write a function for producing inits
  tmp <- rnorm(n = length(coef(fm38.list[[2]])),
      # mean = coef(fm38A), sd = 0.2)  ########
      mean = coef(fm38.list[[2]]), sd = 0.2)
  tmp[tmp < -3] <- -3    # Truncate boundary estimates at -3
  tmp
}
cbind(coef(fm38.list[[2]]), inits()) # Compare solutions and inits
plot(coef(fm38.list[[2]]), inits())  # Plot both
abline(0,1)                          # Add equality line

# Run eight optimisations with different starting values
for (k in 3:10){
  set.seed(k)
  cat("\n\n Run Number", k-2, "\n\n")
  summary(
     fm38.list[[k]] <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
     ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
     ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
     ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
           date + I(date^2) + date:elev + date:I(elev^2) +
           I(date^2):elev + I(date^2):I(elev^2), umf, starts = inits(),
       control=list(trace=TRUE, REPORT=5, maxit = 500), se = FALSE))
}

# Look at NLL of different sets of estimates of model 38
tmptab <- c(sapply(fm38.list, logLik)*-1, -logLik(fm37))
names(tmptab) <- c("free search for model 38:",
    "using solutions from model 37 as inits:",
    paste("random around previous solution", 1:8),
    "model 37:")
# data.frame(NLL=tmptab)
print(data.frame(NLL=tmptab), digits = 10)  # table foot of p277
# ~~~~~~~~~~~ the following table is in the book ~~~~~~~~~~~~~~~~~~~~~~~

#                                                 NLL
# free search for model 38:               2996.728978
# using solutions from model 37 as inits: 2995.688954
# random around previous solution 1       2995.646178
# random around previous solution 2       2995.646266
# random around previous solution 3       2995.646051
# random around previous solution 4       2995.646011
# random around previous solution 5       2995.646451
# random around previous solution 6       2995.646279
# random around previous solution 7       2995.645961
# random around previous solution 8       2995.646141
# model 37:                               2995.707745

# We see again that the free search did not permit optim to find the global optimum here. Using random starting values around a previous solutions suggested NLL = 2995.646 to be the minimum NLL, and showed that fm38.list[[9]] (which is the 'random solution' 7) had the smallest value of the NLL. We therefore take its solutions to be the MLEs.

# Identify model run with smallest NLL
(best <- which.min(tmptab))
# random around previous solution 7
#                                 9

# ~~~~~~~~~ more bonus code ~~~~~~~~~~~~~~
# Fig. 4.14 shows that these 10 models fits produce virtually identical estimates for most parameters, and that discrepancies occurred only for a few parameters with extreme values.
# Plot of estimates from these 10 model fits
oldpar <- par(cex.lab = 1.5, cex.axis = 1.5)
plot(coef(fm38.list[[1]]), xlab = 'Parameter', ylab = 'Estimate', frame = FALSE,
    pch = 1, cex = 1.5)
abline(h = 0)
abline(v = 1:73, col = 'grey')
points(coef(fm38.list[[best]]), pch = 16, cex = 1.5)
for(m in 2:10){
  points(coef(fm38.list[[m]]), pch = 1, cex = 1.5)
  # browser()  # removed for testing
}
par(oldpar)

# So the 'best' solution is fm38.list[[9]]
# ~~~~~~ following code is in the book ~~~~~~~~~~~~~~~~~

# We refit this model and compute the SEs
inits <- coef(fm38.list[[best]])
summary(fm38X <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
        date + I(date^2) + date:elev + date:I(elev^2) +
        I(date^2):elev + I(date^2):I(elev^2), umf, starts = inits,
    control = list(maxit = 500), se = TRUE))

# Compute Chi-square test statistic for actual data by season
library(AICcmodavg)
mb.chisq(fm38X, print.table = TRUE)
# Generate reference distribution of test statistic under H0 (takes 7h)
# system.time( gof <- mb.gof.test(fm38X, print.table = FALSE, nsim = 1000,
#    plot.hist = TRUE, plot.seasons = TRUE, report = 1) )
# ~~~~ for testing ~~~~~~ 15 mins ~~~~
system.time( gof <- mb.gof.test(fm38X, print.table = FALSE, nsim = 10,
    plot.hist = TRUE, plot.seasons = TRUE, report = 1) )
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(c.hat <- gof$c.hat.est)
# [1] 2.102584  # ~~~ from the full run of 1000

AICc(fm38X, c.hat = c.hat) # QAICc of the global model
# [1] 3055.302

# ~~~ save output for use subsequently ~~~~~~~~~
save(fm38X, umf, gof, file="AHM2_04.09.2_output.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
