#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from MS dated 2019-06-17

library(unmarked)
library(AICcmodavg)

load("AHM2_04.09.3_fm50.RData")
c.hat <- 2.102584


# 4.9 Analysis and mapping of crossbill distribution and range dynamics in Switzerland
# ====================================================================================

# 4.9.4 Inference under the AIC-best model
# ----------------------------------------

# ~~~~ code to produce the table ~~~~~~~~~~~
# Produce a table with means, sds and CIs for all parameters
tmp <- summary(fm50)        # Print summary
tmp <- cbind(MLE = coef(fm50), SE = c(tmp[[1]][,2], tmp[[2]][,2], tmp[[3]][,2],
    tmp[[4]][,2]), rbind(confint(fm50, type = "psi"), confint(fm50, type = "col"),
    confint(fm50, type = "ext"), confint(fm50, type = "det")))
print(tmp, 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                               # MLE      SE      0.025     0.975
# psi(Int)                   -0.0793  0.3833   -0.83043   0.67191
# psi(elev)                   1.9507  0.4551    1.05876   2.84256
# psi(I(elev^2))             -1.2445  0.5098   -2.24373  -0.24534
# psi(forest)                -0.7112  0.5104   -1.71165   0.28924
# psi(I(forest^2))            0.1655  0.2613   -0.34668   0.67766
# psi(elev:forest)            0.1665  0.4511   -0.71761   1.05057
# psi(elev:I(forest^2))      -0.7809  0.3540   -1.47469  -0.08721
# psi(I(elev^2):forest)       1.8743  0.5393    0.81743   2.93127
# col(year2001)              -0.0926  0.3264   -0.73241   0.54716
# col(year2002)              -0.3715  0.3821   -1.12040   0.37739
# ....
# col(year2009)              -8.7195 55.8571 -118.19737 100.75829
# col(year2010)              -0.7878  0.4690   -1.70704   0.13149
# col(year2011)              -3.4583  2.3500   -8.06411   1.14756
# col(elev)                   0.5002  0.2060    0.09642   0.90405
# col(I(elev^2))             -0.2159  0.2914   -0.78699   0.35515
# col(forest)                -0.3927  0.2801   -0.94163   0.15632
# col(I(forest^2))            0.2008  0.2139   -0.21852   0.62003
# col(elev:forest)            0.3629  0.1822    0.00591   0.71994
# col(elev:I(forest^2))       0.4318  0.2128    0.01465   0.84887
# col(I(elev^2):forest)       0.5789  0.2884    0.01368   1.14412
# col(I(elev^2):I(forest^2)) -0.9509  0.2876   -1.51459  -0.38726
# ext(year2001)              -1.3218  0.5058   -2.31307  -0.33057
# ext(year2002)              -3.2547  0.7072   -4.64071  -1.86872
# ...
# ext(year2008)              -9.9830 61.7477 -131.00622 111.04026
# ext(year2009)              -2.4937  0.4330   -3.34226  -1.64512
# ext(year2010)              -2.4372  0.5266   -3.46928  -1.40514
# ext(year2011)              -1.7589  0.4253   -2.59239  -0.92533
# ext(elev)                  -1.9752  0.2615   -2.48781  -1.46257
# ext(I(elev^2))              1.1320  0.2998    0.54445   1.71959
# ext(forest)                 0.4192  0.2536   -0.07780   0.91622
# ext(elev:forest)            0.8572  0.2235    0.41918   1.29519
# ext(I(elev^2):forest)      -1.3221  0.2844   -1.87960  -0.76469
# p(year2001)                -0.3512  0.1977   -0.73858   0.03622
# p(year2002)                -0.2137  0.1623   -0.53188   0.10452
# ...
# p(year2011)                -0.7098  0.1449   -0.99376  -0.42582
# p(year2012)                 0.2540  0.1630   -0.06543   0.57334
# p(elev)                     0.4714  0.0774    0.31958   0.62316
# p(forest)                   0.6493  0.0790    0.49455   0.80403
# p(I(forest^2))             -0.2784  0.0488   -0.37412  -0.18274
# p(date)                     0.0816  0.0440   -0.00465   0.16781
# p(I(date^2))                0.1135  0.0383    0.03839   0.18862
# p(elev:forest)              0.4188  0.0674    0.28672   0.55079

# ~~~~~~~~~~~ bonus code ~~~~~~~~~~~~~~~~~~~~
# Estimates accounting for a c-hat = 2.10 (with 95% CIs)
tt <- list()
tt[[1]] <- modavg(list(fm50), parm = '(Intercept)', c.hat = c.hat,
    parm.type = 'psi', warn=FALSE)
tt[[2]] <- modavg(list(fm50), parm = 'elev', c.hat = c.hat,
    parm.type = 'psi', warn=FALSE)
tt[[3]] <- modavg(list(fm50), parm = 'I(elev^2)', c.hat = c.hat,
    parm.type = 'psi', warn=FALSE)
tt[[4]] <- modavg(list(fm50), parm = 'forest', c.hat = c.hat,
    parm.type = 'psi', warn=FALSE)
tt[[5]] <- modavg(list(fm50), parm = 'I(forest^2)', c.hat = c.hat,
    parm.type = 'psi', warn=FALSE)
tt[[6]] <- modavg(list(fm50), parm = 'elev:forest', c.hat = c.hat,
    parm.type = 'psi', warn=FALSE)
tt[[7]] <- modavg(list(fm50), parm = 'elev:I(forest^2)', c.hat = c.hat,
    parm.type = 'psi', warn=FALSE)
tt[[8]] <- modavg(list(fm50), parm = 'I(elev^2):forest', c.hat = c.hat,
    parm.type = 'psi', warn=FALSE)
inflated.uncertainty <- array(NA, dim = c(8, 3),
    dimnames = list(NULL, c('SE(infl)', '95%LCL(infl)', '95%UCL(infl)')))
for(i in 1:8){
   inflated.uncertainty[i, ] <- unlist(tt[[i]][c(4,6,7)])
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compare estimates side by side
print(cbind(tmp[1:8,], inflated.uncertainty), 2)
                         # MLE   SE 0.025  0.975 SE(infl) 95%LCL(infl) 95%UCL(infl)
# psi(Int)              -0.079 0.38 -0.83  0.672     0.56        -1.17         1.01
# psi(elev)              1.951 0.46  1.06  2.843     0.66         0.66         3.24
# psi(I(elev^2))        -1.245 0.51 -2.24 -0.245     0.74        -2.69         0.20
# psi(forest)           -0.711 0.51 -1.71  0.289     0.74        -2.16         0.74
# psi(I(forest^2))       0.165 0.26 -0.35  0.678     0.38        -0.58         0.91
# psi(elev:forest)       0.166 0.45 -0.72  1.051     0.65        -1.12         1.45
# psi(elev:I(forest^2)) -0.781 0.35 -1.47 -0.087     0.51        -1.79         0.22
# psi(I(elev^2):forest)  1.874 0.54  0.82  2.931     0.78         0.34         3.41

