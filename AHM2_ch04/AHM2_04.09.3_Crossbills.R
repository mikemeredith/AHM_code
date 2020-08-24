#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 35 mins
# Run time with the full number of iterations: 48 mins

library(unmarked)
library(AICcmodavg)

# ~~~ load output from 4.9.2 ~~~~~~~~~~
load("AHM2_04.09.2_output.RData")
c.hat <- 2.102584
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 4.9 Analysis and mapping of crossbill distribution and range dynamics in Switzerland
# ====================================================================================

# 4.9.3 Model selection by backwards elimination
# ----------------------------------------------

# ~~~~~ bonus code, not in the printed book, dated 2019-01-04 ~~~~~~~~~~~
# From the summary table of fm38X, the following term is among the least significant:  I(elev^2):I(date^2) in the detection model. In the entire backwards elimination process, we use estimates from a previous model as starting values to speed up convergence. I mark green terms that put a cap on potential model simplification in that part of the model, and yellow terms that are candidates for tossing out.

# Toss out I(elev^2):I(date^2) in detection
inits <- coef(fm38X)[-73]
system.time(
   (fm39 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
     I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2.7 mins
)
summary(fm39)
aictab(list('fm38X' = fm38X, 'fm39' = fm39), c.hat = c.hat)

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm39  73 3051.49        0.00    0.87   0.87 -1424.75
# fm38X 74 3055.30        3.81    0.13   1.00 -1424.74

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)),
    # gammaformula = ~(year - 1) + (elev + I(elev^2)) * (forest +
        # I(forest^2)), epsilonformula = ~(year - 1) + (elev +
        # I(elev^2)) * (forest + I(forest^2)), pformula = ~(year -
        # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) + date +
        # I(date^2) + date:elev + date:I(elev^2) + I(date^2):elev +
        # I(date^2):I(elev^2) - I(elev^2):I(date^2), data = umf,
    # starts = inits, se = TRUE, control = list(trace = TRUE, REPORT = 5,
        # maxit = 500))

# Initial (logit-scale):
                      # Estimate    SE      z   P(>|z|)
# (Intercept)             -0.120 0.384 -0.312 0.7552222
# elev                     1.907 0.457  4.175 0.0000298
# I(elev^2)               -1.151 0.555 -2.075 0.0380319
# forest                  -0.751 0.635 -1.183 0.2366552
# I(forest^2)              0.222 0.386  0.574 0.5656876
# elev:forest              0.153 0.508  0.302 0.7626682
# elev:I(forest^2)        -0.740 0.386 -1.917 0.0551952
# I(elev^2):forest         1.910 0.823  2.322 0.0202338
# I(elev^2):I(forest^2)   -0.113 0.589 -0.192 0.8475527  # toss?

# Colonization (logit-scale):
                      # Estimate     SE      z P(>|z|)
# year2001               -0.0453  0.328 -0.138 0.89024
# year2002               -0.3452  0.390 -0.885 0.37590
# year2003               -1.2499  0.510 -2.451 0.01425
# year2004               -1.5693  0.536 -2.926 0.00343
# year2005               -0.2510  0.429 -0.586 0.55816
# year2006               -1.1757  0.765 -1.538 0.12412
# year2007               -2.5192  1.739 -1.449 0.14743
# year2008                0.0755  0.406  0.186 0.85225
# year2009               -8.7232 44.436 -0.196 0.84437
# year2010               -0.7154  0.474 -1.509 0.13133
# year2011               -3.8713  3.735 -1.036 0.29997
# elev                    0.4486  0.216  2.077 0.03784
# I(elev^2)              -0.1948  0.302 -0.645 0.51892
# forest                 -0.3319  0.314 -1.057 0.29072
# I(forest^2)             0.1303  0.227  0.573 0.56643
# elev:forest             0.3611  0.205  1.765 0.07764
# elev:I(forest^2)        0.4608  0.227  2.026 0.04275
# I(elev^2):forest        0.5342  0.319  1.674 0.09406
# I(elev^2):I(forest^2)  -0.8998  0.307 -2.933 0.00336

# Extinction (logit-scale):
                      # Estimate     SE      z         P(>|z|)
# year2001               -1.1664  0.524 -2.226 0.0260313098214
# year2002               -3.1549  0.702 -4.492 0.0000070464130
# year2003               -2.6358  0.618 -4.262 0.0000202386750
# year2004               -1.6991  0.446 -3.812 0.0001379843516
# year2005               -2.3307  0.749 -3.110 0.0018688060733
# year2006               -3.7510  1.171 -3.204 0.0013575280507
# year2007               -2.3007  0.562 -4.092 0.0000428282988
# year2008               -9.9939 36.209 -0.276 0.7825420603872
# year2009               -2.3715  0.432 -5.491 0.0000000398658
# year2010               -2.3523  0.529 -4.449 0.0000086297425
# year2011               -1.6388  0.434 -3.778 0.0001583341548
# elev                   -2.1651  0.318 -6.798 0.0000000000106
# I(elev^2)               1.1319  0.350  3.235 0.0012148004279
# forest                  0.3881  0.559  0.694 0.4877960888034
# I(forest^2)            -0.0692  0.273 -0.254 0.7996850892375
# elev:forest             0.5499  0.381  1.445 0.1485761636166
# elev:I(forest^2)        0.3523  0.258  1.365 0.1721221817518
# I(elev^2):forest       -1.0257  0.530 -1.937 0.0527834104713
# I(elev^2):I(forest^2)  -0.1805  0.345 -0.523 0.6008275051518  # toss?

# Detection (logit-scale):
                      # Estimate     SE      z  P(>|z|)
# year2001               -0.3598 0.2018 -1.783 7.46e-02
# year2002               -0.2236 0.1777 -1.259 2.08e-01
# year2003                0.0863 0.1527  0.565 5.72e-01
# year2004               -0.4391 0.1585 -2.770 5.60e-03
# year2005                0.3885 0.1859  2.090 3.66e-02
# year2006               -1.0690 0.1858 -5.754 8.71e-09
# year2007               -0.2827 0.1824 -1.550 1.21e-01
# year2008               -0.0713 0.1645 -0.433 6.65e-01
# year2009               -1.1946 0.1367 -8.738 2.36e-18
# year2010                0.3145 0.1494  2.104 3.53e-02
# year2011               -0.7307 0.1549 -4.717 2.39e-06
# year2012                0.2423 0.1727  1.403 1.61e-01
# elev                    0.6150 0.1119  5.495 3.91e-08
# I(elev^2)              -0.0309 0.1247 -0.248 8.04e-01
# forest                  0.5506 0.1355  4.064 4.82e-05
# I(forest^2)            -0.2270 0.0765 -2.967 3.01e-03
# date                    0.1142 0.0582  1.962 4.97e-02
# I(date^2)               0.1473 0.0508  2.899 3.74e-03
# elev:forest             0.3823 0.1483  2.578 9.94e-03
# elev:I(forest^2)       -0.0224 0.0873 -0.256 7.98e-01
# I(elev^2):forest        0.1737 0.1837  0.946 3.44e-01
# I(elev^2):I(forest^2)  -0.1007 0.1114 -0.904 3.66e-01  # toss?
# elev:date              -0.0181 0.0827 -0.219 8.27e-01
# I(elev^2):date          0.0811 0.1100  0.737 4.61e-01
# elev:I(date^2)         -0.1437 0.0729 -1.971 4.87e-02

# AIC: 6135.333
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 116
# Bootstrap iterations: 0

# Toss out I(elev^2):I(forest^2) in occupancy
inits <- coef(fm39)[-9]
system.time(
   (fm40 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
     I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2.7 mins
)
summary(fm40)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)), pformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) + date +
    # I(date^2) + date:elev + date:I(elev^2) + I(date^2):elev +
    # I(date^2):I(elev^2) - I(elev^2):I(date^2), data = umf, starts = inits,
    # se = TRUE, control = list(trace = TRUE, REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)        -0.114 0.380 -0.299 0.7651232
# elev                1.921 0.453  4.237 0.0000226
# I(elev^2)          -1.193 0.511 -2.334 0.0195994
# forest             -0.681 0.507 -1.344 0.1790536
# I(forest^2)         0.167 0.259  0.647 0.5178555
# elev:forest         0.187 0.466  0.401 0.6884183
# elev:I(forest^2)   -0.768 0.356 -2.158 0.0309238
# I(elev^2):forest    1.794 0.541  3.317 0.0009098   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z P(>|z|)
# year2001               -0.0430  0.328 -0.131 0.89563
# year2002               -0.3445  0.389 -0.885 0.37595
# year2003               -1.2499  0.510 -2.450 0.01427
# year2004               -1.5709  0.537 -2.928 0.00341
# year2005               -0.2506  0.428 -0.585 0.55872
# year2006               -1.1750  0.764 -1.538 0.12412
# year2007               -2.5154  1.727 -1.457 0.14518
# year2008                0.0756  0.405  0.186 0.85216
# year2009               -8.7230 44.496 -0.196 0.84458
# year2010               -0.7152  0.474 -1.510 0.13112
# year2011               -3.8743  3.738 -1.037 0.29995
# elev                    0.4483  0.216  2.076 0.03793
# I(elev^2)              -0.1965  0.302 -0.651 0.51483
# forest                 -0.3336  0.313 -1.064 0.28722
# I(forest^2)             0.1306  0.227  0.575 0.56527
# elev:forest             0.3615  0.205  1.767 0.07721
# elev:I(forest^2)        0.4586  0.227  2.021 0.04330
# I(elev^2):forest        0.5360  0.319  1.682 0.09251
# I(elev^2):I(forest^2)  -0.8967  0.307 -2.925 0.00344   # Blocks things here

# Extinction (logit-scale):
                      # Estimate     SE      z          P(>|z|)
# year2001               -1.1622  0.522 -2.226 0.02599975269952
# year2002               -3.1539  0.701 -4.497 0.00000689233177
# year2003               -2.6368  0.618 -4.265 0.00001998978872
# year2004               -1.6986  0.445 -3.814 0.00013648590551
# year2005               -2.3314  0.750 -3.109 0.00187448445593
# year2006               -3.7502  1.169 -3.209 0.00133192343702
# year2007               -2.2997  0.562 -4.093 0.00004259788041
# year2008               -9.9944 36.043 -0.277 0.78155540664620
# year2009               -2.3716  0.432 -5.493 0.00000003943543
# year2010               -2.3515  0.528 -4.450 0.00000858497696
# year2011               -1.6377  0.433 -3.779 0.00015749975577
# elev                   -2.1671  0.318 -6.807 0.00000000000999
# I(elev^2)               1.1306  0.350  3.231 0.00123393683333
# forest                  0.3860  0.558  0.692 0.48874784741420
# I(forest^2)            -0.0683  0.272 -0.251 0.80180604155955
# elev:forest             0.5483  0.380  1.443 0.14914145157491
# elev:I(forest^2)        0.3536  0.258  1.372 0.17012054894632
# I(elev^2):forest       -1.0227  0.528 -1.937 0.05280289530641
# I(elev^2):I(forest^2)  -0.1813  0.345 -0.526 0.59875864546732  # toss?

# Detection (logit-scale):
                      # Estimate     SE      z  P(>|z|)
# year2001               -0.3578 0.2013 -1.777 7.56e-02
# year2002               -0.2233 0.1776 -1.257 2.09e-01
# year2003                0.0866 0.1526  0.568 5.70e-01
# year2004               -0.4393 0.1584 -2.773 5.56e-03
# year2005                0.3890 0.1858  2.094 3.63e-02
# year2006               -1.0683 0.1857 -5.752 8.84e-09
# year2007               -0.2827 0.1820 -1.553 1.20e-01
# year2008               -0.0708 0.1645 -0.430 6.67e-01
# year2009               -1.1952 0.1367 -8.744 2.26e-18
# year2010                0.3152 0.1494  2.109 3.49e-02
# year2011               -0.7308 0.1548 -4.720 2.36e-06
# year2012                0.2425 0.1726  1.405 1.60e-01
# elev                    0.6146 0.1119  5.494 3.94e-08
# I(elev^2)              -0.0306 0.1247 -0.246 8.06e-01
# forest                  0.5484 0.1348  4.067 4.76e-05
# I(forest^2)            -0.2257 0.0762 -2.961 3.06e-03
# date                    0.1141 0.0582  1.961 4.99e-02
# I(date^2)               0.1473 0.0508  2.898 3.75e-03
# elev:forest             0.3793 0.1477  2.568 1.02e-02
# elev:I(forest^2)       -0.0207 0.0870 -0.238 8.12e-01
# I(elev^2):forest        0.1785 0.1822  0.980 3.27e-01
# I(elev^2):I(forest^2)  -0.1036 0.1107 -0.936 3.49e-01
# elev:date              -0.0180 0.0827 -0.218 8.27e-01
# I(elev^2):date          0.0811 0.1100  0.737 4.61e-01
# elev:I(date^2)         -0.1437 0.0729 -1.971 4.87e-02   # Blocks things here

# AIC: 6133.371
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 132
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm40  72 3047.71        0.00    0.85   0.85 -1424.76
# fm39  73 3051.49        3.78    0.13   0.98 -1424.75
# fm38X 74 3055.30        7.59    0.02   1.00 -1424.74

# Toss out I(elev^2):I(forest^2) in extinction
inits <- coef(fm40)[-46]
system.time(
   (fm41 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
     I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 3 mins
)
summary(fm41)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2),
    # pformula = ~(year - 1) + (elev + I(elev^2)) * (forest + I(forest^2)) +
        # date + I(date^2) + date:elev + date:I(elev^2) + I(date^2):elev +
        # I(date^2):I(elev^2) - I(elev^2):I(date^2), data = umf,
    # starts = inits, se = TRUE, control = list(trace = TRUE, REPORT = 5,
        # maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)        -0.111 0.380 -0.291 0.7708162
# elev                1.915 0.452  4.234 0.0000229
# I(elev^2)          -1.197 0.510 -2.348 0.0188506
# forest             -0.685 0.507 -1.350 0.1768656
# I(forest^2)         0.168 0.259  0.646 0.5183727
# elev:forest         0.196 0.465  0.420 0.6741588
# elev:I(forest^2)   -0.772 0.356 -2.167 0.0302531
# I(elev^2):forest    1.803 0.540  3.338 0.0008440   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z P(>|z|)
# year2001               -0.0505  0.329 -0.154 0.87797
# year2002               -0.3640  0.392 -0.929 0.35312
# year2003               -1.2546  0.511 -2.454 0.01413
# year2004               -1.5854  0.536 -2.957 0.00310
# year2005               -0.2476  0.431 -0.574 0.56588
# year2006               -1.1850  0.775 -1.530 0.12612
# year2007               -2.7420  1.991 -1.377 0.16846
# year2008                0.0535  0.406  0.132 0.89506
# year2009               -8.7227 43.256 -0.202 0.84019
# year2010               -0.7370  0.472 -1.561 0.11863
# year2011               -3.7598  3.382 -1.112 0.26628
# elev                    0.4382  0.214  2.047 0.04070
# I(elev^2)              -0.2036  0.302 -0.674 0.50005
# forest                 -0.2913  0.304 -0.959 0.33777
# I(forest^2)             0.1109  0.225  0.493 0.62170
# elev:forest             0.3734  0.202  1.852 0.06398
# elev:I(forest^2)        0.4535  0.227  1.998 0.04573
# I(elev^2):forest        0.4921  0.307  1.602 0.10919
# I(elev^2):I(forest^2)  -0.8688  0.302 -2.875 0.00404   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE     z          P(>|z|)
# year2001           -1.208  0.525 -2.30 0.02134327394998
# year2002           -3.193  0.710 -4.49 0.00000696449348
# year2003           -2.690  0.615 -4.37 0.00001235537315
# year2004           -1.719  0.462 -3.72 0.00019716118957
# year2005           -2.394  0.735 -3.26 0.00111875788862
# year2006           -3.882  1.144 -3.39 0.00069232598372
# year2007           -2.335  0.569 -4.11 0.00004033508907
# year2008           -9.994 38.389 -0.26 0.79459732875215
# year2009           -2.393  0.432 -5.54 0.00000002981053
# year2010           -2.369  0.534 -4.44 0.00000913864424
# year2011           -1.675  0.438 -3.82 0.00013186977778
# elev               -2.212  0.312 -7.09 0.00000000000131
# I(elev^2)           1.076  0.348  3.09 0.00197987258941
# forest              0.587  0.430  1.37 0.17218401258271
# I(forest^2)        -0.181  0.175 -1.04 0.29941762466661
# elev:forest         0.675  0.304  2.22 0.02633503603338
# elev:I(forest^2)    0.288  0.225  1.28 0.20103832256360  # toss?
# I(elev^2):forest   -1.229  0.369 -3.33 0.00085537191641   # Blocks things here

# Detection (logit-scale):
                      # Estimate     SE      z                P(>|z|)
# year2001               -0.3603 0.2016 -1.787 0.07391502285595356880
# year2002               -0.2323 0.1763 -1.318 0.18751768342891211860
# year2003                0.0827 0.1523  0.543 0.58701140396006679101
# year2004               -0.4479 0.1562 -2.867 0.00414330913934032700
# year2005                0.3809 0.1880  2.026 0.04275393726377579168
# year2006               -1.0798 0.1824 -5.920 0.00000000322334745069
# year2007               -0.2997 0.1734 -1.728 0.08400509609313749904
# year2008               -0.0787 0.1634 -0.482 0.62996750989028049705
# year2009               -1.1971 0.1368 -8.752 0.00000000000000000209
# year2010                0.3118 0.1484  2.101 0.03561615926803920334
# year2011               -0.7317 0.1553 -4.713 0.00000244286323231848
# year2012                0.2273 0.1692  1.344 0.17901629020583559315
# elev                    0.6109 0.1110  5.501 0.00000003775691397685
# I(elev^2)              -0.0313 0.1243 -0.251 0.80154905379461094395
# forest                  0.5686 0.1286  4.421 0.00000984206695265386
# I(forest^2)            -0.2352 0.0738 -3.189 0.00142985589478902288
# date                    0.1146 0.0581  1.971 0.04876177633311713777
# I(date^2)               0.1481 0.0507  2.918 0.00352289665326563145
# elev:forest             0.3844 0.1466  2.622 0.00875225168411745275
# elev:I(forest^2)       -0.0223 0.0867 -0.258 0.79656412918095798048
# I(elev^2):forest        0.1618 0.1793  0.902 0.36700703937934031629
# I(elev^2):I(forest^2)  -0.0931 0.1090 -0.854 0.39299842353875114043  # toss?
# elev:date              -0.0186 0.0826 -0.225 0.82203895434382134866
# I(elev^2):date          0.0814 0.1099  0.741 0.45868089179333948469
# elev:I(date^2)         -0.1440 0.0728 -1.977 0.04806435180384371203   # Blocks things here

# AIC: 6131.644
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 182
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm41  71 3044.09        0.00    0.84   0.84 -1424.83
# fm40  72 3047.71        3.63    0.14   0.98 -1424.76
# fm39  73 3051.49        7.40    0.02   1.00 -1424.75
# fm38X 74 3055.30       11.21    0.00   1.00 -1424.74

# Toss out I(elev^2):I(forest^2) in detection
inits <- coef(fm41)[-67]
system.time(
   (fm42 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
     I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2.6 mins
)
summary(fm42)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41, 'fm42' = fm42), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2),
    # pformula = ~(year - 1) + (elev + I(elev^2)) * (forest + I(forest^2)) +
        # date + I(date^2) + date:elev + date:I(elev^2) + I(date^2):elev +
        # I(date^2):I(elev^2) - I(elev^2):I(date^2) - I(elev^2):I(forest^2),
    # data = umf, starts = inits, se = TRUE, control = list(trace = TRUE,
        # REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)       -0.0959 0.382 -0.251 0.8016320
# elev               1.9187 0.455  4.221 0.0000244
# I(elev^2)         -1.2168 0.510 -2.384 0.0171342
# forest            -0.6971 0.511 -1.363 0.1728227
# I(forest^2)        0.1667 0.261  0.640 0.5222487
# elev:forest        0.1529 0.456  0.335 0.7372582
# elev:I(forest^2)  -0.7569 0.355 -2.134 0.0328499
# I(elev^2):forest   1.8420 0.539  3.418 0.0006304   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z P(>|z|)
# year2001               -0.0563  0.330 -0.171 0.86436
# year2002               -0.3692  0.394 -0.938 0.34829
# year2003               -1.2557  0.511 -2.456 0.01403
# year2004               -1.5641  0.533 -2.935 0.00334
# year2005               -0.2387  0.431 -0.554 0.57954
# year2006               -1.1924  0.779 -1.531 0.12576
# year2007               -2.7222  1.967 -1.384 0.16629
# year2008                0.0553  0.405  0.137 0.89142
# year2009               -8.7226 43.703 -0.200 0.84180
# year2010               -0.7393  0.473 -1.562 0.11836
# year2011               -3.7098  3.291 -1.127 0.25963
# elev                    0.4325  0.214  2.017 0.04373
# I(elev^2)              -0.1882  0.304 -0.620 0.53543
# forest                 -0.3079  0.303 -1.015 0.30998
# I(forest^2)             0.1224  0.225  0.543 0.58702
# elev:forest             0.3379  0.195  1.736 0.08262
# elev:I(forest^2)        0.4970  0.219  2.268 0.02335
# I(elev^2):forest        0.5343  0.302  1.770 0.07667
# I(elev^2):I(forest^2)  -0.9253  0.293 -3.153 0.00162   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE     z           P(>|z|)
# year2001           -1.220  0.524 -2.33 0.019833539511548
# year2002           -3.222  0.724 -4.45 0.000008586275409
# year2003           -2.679  0.612 -4.38 0.000012006705814
# year2004           -1.724  0.461 -3.74 0.000186225715137
# year2005           -2.385  0.730 -3.27 0.001086550670460
# year2006           -3.878  1.138 -3.41 0.000653619863533
# year2007           -2.349  0.570 -4.12 0.000038081807671
# year2008           -9.994 41.648 -0.24 0.810359599468290
# year2009           -2.400  0.432 -5.56 0.000000026688000
# year2010           -2.377  0.535 -4.44 0.000008836734543
# year2011           -1.693  0.438 -3.87 0.000110524898578
# elev               -2.210  0.310 -7.14 0.000000000000952
# I(elev^2)           1.100  0.347  3.17 0.001521476367323
# forest              0.606  0.427  1.42 0.156295214304632
# I(forest^2)        -0.189  0.174 -1.08 0.279619235028127
# elev:forest         0.677  0.303  2.24 0.025141406021983
# elev:I(forest^2)    0.293  0.225  1.30 0.192354582120452
# I(elev^2):forest   -1.265  0.369 -3.42 0.000617262831241   # Blocks things here

# Detection (logit-scale):
                 # Estimate     SE      z                P(>|z|)
# year2001          -0.3576 0.2014 -1.776 0.07581363710000138534
# year2002          -0.2280 0.1756 -1.299 0.19409136257103529188
# year2003           0.0852 0.1524  0.559 0.57632654511724812352
# year2004          -0.4403 0.1561 -2.821 0.00478169038824988093
# year2005           0.3861 0.1875  2.059 0.03952063567940369143
# year2006          -1.0765 0.1829 -5.885 0.00000000397505886613
# year2007          -0.2926 0.1737 -1.685 0.09207008818663685312
# year2008          -0.0727 0.1632 -0.446 0.65588780108463951013
# year2009          -1.1923 0.1368 -8.718 0.00000000000000000283
# year2010           0.3182 0.1482  2.148 0.03170305923391564251
# year2011          -0.7255 0.1556 -4.663 0.00000311785443882654
# year2012           0.2275 0.1685  1.350 0.17697954778536170717
# elev               0.6178 0.1105  5.589 0.00000002289181166898
# I(elev^2)         -0.0609 0.1193 -0.510 0.60989839622643415851
# forest             0.6274 0.1083  5.793 0.00000000690720418530
# I(forest^2)       -0.2798 0.0521 -5.375 0.00000007647306199955
# date               0.1143 0.0582  1.963 0.04961977990682846923
# I(date^2)          0.1499 0.0507  2.954 0.00313383786453014105
# elev:forest        0.4593 0.1170  3.927 0.00008589734440792231
# elev:I(forest^2)  -0.0620 0.0729 -0.851 0.39481762549817717156
# I(elev^2):forest   0.0429 0.1131  0.380 0.70422060150623355668  # toss?
# elev:date         -0.0206 0.0826 -0.249 0.80333458642858224241
# I(elev^2):date     0.0858 0.1098  0.781 0.43457137317539323851
# elev:I(date^2)    -0.1475 0.0727 -2.029 0.04241755751404011066   # Blocks things here

# AIC: 6130.373
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 151
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm42  70 3040.72        0.00    0.82   0.82 -1425.00
# fm41  71 3044.09        3.37    0.15   0.97 -1424.83
# fm40  72 3047.71        6.99    0.02   1.00 -1424.76
# fm39  73 3051.49       10.77    0.00   1.00 -1424.75
# fm38X 74 3055.30       14.58    0.00   1.00 -1424.74

# Toss out I(elev^2):forest in detection
inits <- coef(fm42)[-66]
system.time(
   (fm43 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
     I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest,
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE)) # 2.7 mins
)
summary(fm43)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41, 'fm42' = fm42, 'fm43' = fm43), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2),
    # pformula = ~(year - 1) + (elev + I(elev^2)) * (forest + I(forest^2)) +
        # date + I(date^2) + date:elev + date:I(elev^2) + I(date^2):elev +
        # I(date^2):I(elev^2) - I(elev^2):I(date^2) - I(elev^2):I(forest^2) -
        # I(elev^2):forest, data = umf, starts = inits, se = TRUE,
    # control = list(trace = TRUE, REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z  P(>|z|)
# (Intercept)       -0.0873 0.383 -0.228 0.819616
# elev               1.9242 0.454  4.234 0.000023
# I(elev^2)         -1.2312 0.510 -2.413 0.015816
# forest            -0.7096 0.512 -1.385 0.166036
# I(forest^2)        0.1668 0.261  0.640 0.522343
# elev:forest        0.1388 0.453  0.306 0.759517
# elev:I(forest^2)  -0.7513 0.354 -2.122 0.033831
# I(elev^2):forest   1.8632 0.537  3.470 0.000520   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z P(>|z|)
# year2001               -0.0615  0.329 -0.187 0.85178
# year2002               -0.3763  0.394 -0.956 0.33905
# year2003               -1.2573  0.509 -2.469 0.01355
# year2004               -1.5513  0.528 -2.936 0.00333
# year2005               -0.2440  0.430 -0.568 0.57014
# year2006               -1.1809  0.768 -1.538 0.12402
# year2007               -2.7598  2.044 -1.350 0.17688
# year2008                0.0546  0.404  0.135 0.89251
# year2009               -8.7221 44.217 -0.197 0.84363
# year2010               -0.7543  0.467 -1.617 0.10592
# year2011               -3.4965  2.522 -1.386 0.16565
# elev                    0.4400  0.213  2.068 0.03861
# I(elev^2)              -0.1917  0.303 -0.634 0.52634
# forest                 -0.3067  0.304 -1.010 0.31245
# I(forest^2)             0.1276  0.224  0.569 0.56945
# elev:forest             0.3269  0.191  1.710 0.08719
# elev:I(forest^2)        0.5032  0.217  2.317 0.02052
# I(elev^2):forest        0.5403  0.300  1.798 0.07213
# I(elev^2):I(forest^2)  -0.9355  0.291 -3.217 0.00130   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE      z           P(>|z|)
# year2001           -1.243  0.522 -2.383 0.017191449998245
# year2002           -3.237  0.726 -4.459 0.000008252213814
# year2003           -2.683  0.612 -4.385 0.000011615761277
# year2004           -1.741  0.464 -3.753 0.000174467205284
# year2005           -2.378  0.722 -3.292 0.000995696591907
# year2006           -3.895  1.133 -3.437 0.000587552856548
# year2007           -2.366  0.570 -4.150 0.000033277658124
# year2008           -9.993 42.768 -0.234 0.815253630139557
# year2009           -2.409  0.430 -5.599 0.000000021557569
# year2010           -2.380  0.534 -4.458 0.000008290883865
# year2011           -1.710  0.437 -3.915 0.000090557373625
# elev               -2.210  0.309 -7.160 0.000000000000807
# I(elev^2)           1.116  0.346  3.228 0.001247955241152
# forest              0.635  0.422  1.504 0.132552561818856
# I(forest^2)        -0.192  0.174 -1.103 0.269829532688534
# elev:forest         0.689  0.302  2.282 0.022504955836342
# elev:I(forest^2)    0.288  0.224  1.281 0.200031595953416  # toss?
# I(elev^2):forest   -1.295  0.362 -3.575 0.000350504952661   # Blocks things here

# Detection (logit-scale):
                 # Estimate     SE      z                P(>|z|)
# year2001          -0.3666 0.2001 -1.833 0.06686948896021108202
# year2002          -0.2409 0.1719 -1.401 0.16108273013185864242
# year2003           0.0733 0.1490  0.492 0.62283126349619610806
# year2004          -0.4508 0.1539 -2.929 0.00340319013807123538
# year2005           0.3724 0.1847  2.017 0.04373258929387979776
# year2006          -1.0858 0.1813 -5.988 0.00000000212045222125
# year2007          -0.3047 0.1707 -1.785 0.07432555091161795191
# year2008          -0.0816 0.1615 -0.506 0.61318225856878605384
# year2009          -1.2039 0.1334 -9.025 0.00000000000000000018
# year2010           0.3098 0.1462  2.119 0.03406178924353796961
# year2011          -0.7335 0.1529 -4.796 0.00000161639386467754
# year2012           0.2170 0.1657  1.310 0.19028966440469297083
# elev               0.6112 0.1089  5.615 0.00000001964070020444
# I(elev^2)         -0.0374 0.1014 -0.369 0.71202745664088784583
# forest             0.6520 0.0867  7.523 0.00000000000005362559
# I(forest^2)       -0.2817 0.0518 -5.437 0.00000005418712208687
# date               0.1184 0.0573  2.067 0.03872942561923223975
# I(date^2)          0.1504 0.0508  2.961 0.00306618261847400007
# elev:forest        0.4818 0.1010  4.770 0.00000183997513844583
# elev:I(forest^2)  -0.0726 0.0675 -1.076 0.28204116144084806495  # toss?
# elev:date         -0.0211 0.0826 -0.255 0.79886790375801197683
# I(elev^2):date     0.0764 0.1074  0.712 0.47664054009585182792  # toss?
# elev:I(date^2)    -0.1460 0.0727 -2.009 0.04451771371781823583   # Blocks things here

# AIC: 6128.514
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 166
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm43  69 3037.11        0.00    0.83   0.83 -1425.04
# fm42  70 3040.72        3.61    0.14   0.97 -1425.00
# fm41  71 3044.09        6.98    0.03   1.00 -1424.83
# fm40  72 3047.71       10.61    0.00   1.00 -1424.76
# fm39  73 3051.49       14.38    0.00   1.00 -1424.75
# fm38X 74 3055.30       18.20    0.00   1.00 -1424.74

# Toss out I(elev^2):date in detection
inits <- coef(fm43)[-67]
system.time(
   (fm44 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
     I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date,
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2.5 mins
)
summary(fm44)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41,
    'fm42' = fm42, 'fm43' = fm43, 'fm44' = fm44), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2),
    # pformula = ~(year - 1) + (elev + I(elev^2)) * (forest + I(forest^2)) +
        # date + I(date^2) + date:elev + date:I(elev^2) + I(date^2):elev +
        # I(date^2):I(elev^2) - I(elev^2):I(date^2) - I(elev^2):I(forest^2) -
        # I(elev^2):forest - I(elev^2):date, data = umf, starts = inits,
    # se = TRUE, control = list(trace = TRUE, REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)        -0.094 0.381 -0.246 0.8053092
# elev                1.925 0.453  4.246 0.0000217
# I(elev^2)          -1.228 0.508 -2.415 0.0157358
# forest             -0.701 0.511 -1.373 0.1696383
# I(forest^2)         0.166 0.260  0.638 0.5235038
# elev:forest         0.147 0.452  0.325 0.7448983
# elev:I(forest^2)   -0.755 0.354 -2.133 0.0329446
# I(elev^2):forest    1.853 0.535  3.464 0.0005329   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z P(>|z|)
# year2001               -0.0489  0.328 -0.149 0.88146
# year2002               -0.3673  0.393 -0.936 0.34938
# year2003               -1.2499  0.508 -2.460 0.01389
# year2004               -1.5435  0.528 -2.923 0.00347
# year2005               -0.2438  0.430 -0.567 0.57098
# year2006               -1.1857  0.771 -1.539 0.12386
# year2007               -2.6642  1.864 -1.430 0.15282
# year2008                0.0532  0.404  0.132 0.89537
# year2009               -8.7219 44.064 -0.198 0.84309
# year2010               -0.7458  0.465 -1.604 0.10870
# year2011               -3.5208  2.569 -1.370 0.17058
# elev                    0.4441  0.212  2.091 0.03650
# I(elev^2)              -0.2087  0.300 -0.695 0.48686
# forest                 -0.3041  0.303 -1.002 0.31615
# I(forest^2)             0.1208  0.224  0.538 0.59037
# elev:forest             0.3267  0.190  1.715 0.08633
# elev:I(forest^2)        0.4998  0.217  2.305 0.02115
# I(elev^2):forest        0.5397  0.300  1.799 0.07202
# I(elev^2):I(forest^2)  -0.9236  0.290 -3.185 0.00145   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE      z          P(>|z|)
# year2001           -1.238  0.521 -2.378 0.01740267265844
# year2002           -3.225  0.722 -4.468 0.00000787786323
# year2003           -2.674  0.610 -4.383 0.00001168403754
# year2004           -1.740  0.464 -3.754 0.00017388682794
# year2005           -2.383  0.731 -3.261 0.00110963139489
# year2006           -3.872  1.123 -3.448 0.00056411653737
# year2007           -2.364  0.570 -4.148 0.00003358889138
# year2008           -9.993 42.079 -0.237 0.81228982262475
# year2009           -2.412  0.431 -5.589 0.00000002283287
# year2010           -2.380  0.535 -4.449 0.00000863546564
# year2011           -1.699  0.436 -3.897 0.00009755681894
# elev               -2.207  0.310 -7.125 0.00000000000104
# I(elev^2)           1.107  0.344  3.214 0.00130909773932
# forest              0.629  0.423  1.487 0.13689392249139
# I(forest^2)        -0.193  0.175 -1.106 0.26865542765377
# elev:forest         0.684  0.302  2.264 0.02354612053018
# elev:I(forest^2)    0.289  0.225  1.286 0.19857323839838  # toss?
# I(elev^2):forest   -1.282  0.360 -3.555 0.00037728288384   # Blocks things here

# Detection (logit-scale):
                 # Estimate     SE      z                 P(>|z|)
# year2001          -0.3656 0.2000 -1.828 0.067576422279007383742
# year2002          -0.2416 0.1718 -1.407 0.159521674180083061767
# year2003           0.0715 0.1488  0.481 0.630679243071087691774
# year2004          -0.4515 0.1540 -2.932 0.003364086156159064919
# year2005           0.3709 0.1844  2.011 0.044282562885703570521
# year2006          -1.0861 0.1813 -5.990 0.000000002101193897252
# year2007          -0.3044 0.1707 -1.783 0.074646212069622144836
# year2008          -0.0858 0.1613 -0.532 0.594900646453873749309
# year2009          -1.2034 0.1334 -9.020 0.000000000000000000188
# year2010           0.3069 0.1462  2.099 0.035815079108475832148
# year2011          -0.7397 0.1522 -4.859 0.000001179480834723715
# year2012           0.2156 0.1659  1.299 0.193953691084641804965
# elev               0.6081 0.1088  5.591 0.000000022568103728522
# I(elev^2)         -0.0198 0.0985 -0.201 0.840998506552927072022
# forest             0.6473 0.0865  7.481 0.000000000000073480308
# I(forest^2)       -0.2811 0.0518 -5.422 0.000000058777812765438
# date               0.1365 0.0514  2.655 0.007929168371991828293
# I(date^2)          0.1507 0.0508  2.964 0.003033372130728684713
# elev:forest        0.4753 0.1005  4.730 0.000002240927867279611
# elev:I(forest^2)  -0.0699 0.0673 -1.037 0.299577073721385478411  # toss?
# elev:date         -0.0155 0.0826 -0.188 0.850890845499991166356
# elev:I(date^2)    -0.1114 0.0540 -2.063 0.039146445127194548497

# AIC: 6127.023
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 141
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm44  68 3033.71        0.00    0.82   0.82 -1425.16
# fm43  69 3037.11        3.40    0.15   0.97 -1425.04
# fm42  70 3040.72        7.01    0.02   0.99 -1425.00
# fm41  71 3044.09       10.38    0.00   1.00 -1424.83
# fm40  72 3047.71       14.01    0.00   1.00 -1424.76
# fm39  73 3051.49       17.78    0.00   1.00 -1424.75
# fm38X 74 3055.30       21.60    0.00   1.00 -1424.74

# Toss out elev:I(forest^2) in detection
inits <- coef(fm44)[-65]
system.time(
   (fm45 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date -
      elev:I(forest^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2.5 mins
)
summary(fm45)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41,
    'fm42' = fm42, 'fm43' = fm43, 'fm44' = fm44, 'fm45' = fm45), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2),
    # pformula = ~(year - 1) + (elev + I(elev^2)) * (forest + I(forest^2)) +
        # date + I(date^2) + date:elev + date:I(elev^2) + I(date^2):elev +
        # I(date^2):I(elev^2) - I(elev^2):I(date^2) - I(elev^2):I(forest^2) -
        # I(elev^2):forest - I(elev^2):date - elev:I(forest^2),
    # data = umf, starts = inits, se = TRUE, control = list(trace = TRUE,
        # REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)        -0.092 0.383 -0.240 0.8100265
# elev                1.941 0.452  4.289 0.0000179
# I(elev^2)          -1.241 0.509 -2.441 0.0146646
# forest             -0.725 0.514 -1.412 0.1579277
# I(forest^2)         0.178 0.261  0.681 0.4955646
# elev:forest         0.196 0.450  0.436 0.6626068
# elev:I(forest^2)   -0.786 0.353 -2.229 0.0258468
# I(elev^2):forest    1.855 0.535  3.469 0.0005229

# Colonization (logit-scale):
                      # Estimate     SE      z P(>|z|)
# year2001               -0.0564  0.328 -0.172 0.86334
# year2002               -0.3870  0.394 -0.983 0.32555
# year2003               -1.2515  0.506 -2.473 0.01340
# year2004               -1.5301  0.522 -2.929 0.00340
# year2005               -0.2697  0.428 -0.630 0.52841
# year2006               -1.1501  0.736 -1.562 0.11831
# year2007               -2.7200  1.929 -1.410 0.15849
# year2008                0.0620  0.402  0.154 0.87757
# year2009               -8.7218 43.104 -0.202 0.83965
# year2010               -0.7451  0.465 -1.601 0.10933
# year2011               -3.4638  2.482 -1.395 0.16287
# elev                    0.4528  0.210  2.152 0.03136
# I(elev^2)              -0.2162  0.299 -0.724 0.46925
# forest                 -0.2838  0.305 -0.931 0.35175
# I(forest^2)             0.1129  0.224  0.505 0.61363
# elev:forest             0.3506  0.188  1.870 0.06155
# elev:I(forest^2)        0.4717  0.214  2.206 0.02738
# I(elev^2):forest        0.5147  0.299  1.724 0.08469
# I(elev^2):I(forest^2)  -0.8959  0.286 -3.131 0.00174

# Extinction (logit-scale):
                 # Estimate     SE      z           P(>|z|)
# year2001           -1.273  0.524 -2.430 0.015083230751263
# year2002           -3.258  0.721 -4.519 0.000006224807175
# year2003           -2.697  0.608 -4.440 0.000008996619528
# year2004           -1.765  0.468 -3.767 0.000165327947638
# year2005           -2.388  0.729 -3.273 0.001062856918601
# year2006           -3.881  1.111 -3.494 0.000475778990504
# year2007           -2.380  0.570 -4.176 0.000029614192439
# year2008           -9.992 43.221 -0.231 0.817170526667665
# year2009           -2.417  0.430 -5.620 0.000000019150754
# year2010           -2.386  0.537 -4.444 0.000008816574850
# year2011           -1.713  0.438 -3.913 0.000091093232302
# elev               -2.225  0.305 -7.284 0.000000000000323
# I(elev^2)           1.136  0.346  3.284 0.001023678956340
# forest              0.670  0.427  1.570 0.116300978814268
# I(forest^2)        -0.203  0.175 -1.160 0.245961278558443
# elev:forest         0.653  0.299  2.183 0.029022055208857
# elev:I(forest^2)    0.316  0.221  1.431 0.152572695962936  # toss?
# I(elev^2):forest   -1.314  0.362 -3.631 0.000282033317312

# Detection (logit-scale):
               # Estimate     SE      z                 P(>|z|)
# year2001        -0.3588 0.2006 -1.788 0.073695450500577283637
# year2002        -0.2434 0.1715 -1.420 0.155743981669096337450
# year2003         0.0733 0.1484  0.494 0.621394475235115262279
# year2004        -0.4535 0.1533 -2.959 0.003087619611562527830
# year2005         0.3671 0.1852  1.982 0.047446712582744056397
# year2006        -1.0812 0.1814 -5.960 0.000000002526310108337
# year2007        -0.3062 0.1696 -1.806 0.070957205824624086654
# year2008        -0.0825 0.1614 -0.511 0.609036222929274395632
# year2009        -1.2079 0.1334 -9.054 0.000000000000000000138
# year2010         0.3094 0.1458  2.123 0.033788594598392215906
# year2011        -0.7404 0.1532 -4.834 0.000001337513802822184
# year2012         0.2198 0.1663  1.322 0.186250551303934364089
# elev             0.5644 0.0995  5.672 0.000000014134417228935
# I(elev^2)       -0.0085 0.0978 -0.087 0.930681007717892305742
# forest           0.6731 0.0832  8.090 0.000000000000000597023
# I(forest^2)     -0.2933 0.0507 -5.786 0.000000007223521220065
# date             0.1373 0.0515  2.668 0.007619852826309091338
# I(date^2)        0.1516 0.0509  2.979 0.002894729949930842719
# elev:forest      0.4013 0.0709  5.662 0.000000014941232156517
# elev:date       -0.0163 0.0826 -0.197 0.843788863578160519552
# elev:I(date^2)  -0.1118 0.0541 -2.068 0.038661465090657709531

# AIC: 6126.104
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 145
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm45  67 3030.62        0.00    0.79   0.79 -1425.41
# fm44  68 3033.71        3.09    0.17   0.96 -1425.16
# fm43  69 3037.11        6.49    0.03   0.99 -1425.04
# fm42  70 3040.72       10.10    0.01   1.00 -1425.00
# fm41  71 3044.09       13.47    0.00   1.00 -1424.83
# fm40  72 3047.71       17.10    0.00   1.00 -1424.76
# fm39  73 3051.49       20.87    0.00   1.00 -1424.75
# fm38X 74 3055.30       24.69    0.00   1.00 -1424.74

# Toss out elev:I(forest^2) in extinction
inits <- coef(fm45)[-44]
system.time(
   (fm46 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2) - elev:I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date -
      elev:I(forest^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2.4 mins
)
summary(fm46)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41,
    'fm42' = fm42, 'fm43' = fm43, 'fm44' = fm44, 'fm45' = fm45, 'fm46' = fm46),
    c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2) -
    # elev:I(forest^2), pformula = ~(year - 1) + (elev + I(elev^2)) *
    # (forest + I(forest^2)) + date + I(date^2) + date:elev + date:I(elev^2) +
    # I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
    # I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date -
    # elev:I(forest^2), data = umf, starts = inits, se = TRUE,
    # control = list(trace = TRUE, REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)       -0.0878 0.384 -0.229 0.8191819
# elev               1.9311 0.452  4.271 0.0000194
# I(elev^2)         -1.2339 0.509 -2.425 0.0153016
# forest            -0.7249 0.514 -1.411 0.1583541
# I(forest^2)        0.1758 0.261  0.672 0.5013538
# elev:forest        0.2075 0.449  0.463 0.6437143
# elev:I(forest^2)  -0.7914 0.353 -2.243 0.0249066
# I(elev^2):forest   1.8540 0.535  3.465 0.0005312   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE       z P(>|z|)
# year2001               -0.0822  0.329 -0.2501 0.80251
# year2002               -0.4308  0.397 -1.0863 0.27734
# year2003               -1.2591  0.501 -2.5117 0.01202
# year2004               -1.5554  0.522 -2.9771 0.00291
# year2005               -0.2759  0.434 -0.6364 0.52452
# year2006               -1.2266  0.772 -1.5888 0.11211
# year2007               -2.7600  2.013 -1.3713 0.17027
# year2008                0.0242  0.407  0.0593 0.95271
# year2009               -8.7223 39.271 -0.2221 0.82423
# year2010               -0.7847  0.464 -1.6927 0.09052
# year2011               -3.4402  2.284 -1.5065 0.13194
# elev                    0.4911  0.210  2.3425 0.01916
# I(elev^2)              -0.1746  0.299 -0.5848 0.55870
# forest                 -0.2589  0.302 -0.8581 0.39085
# I(forest^2)             0.1184  0.222  0.5329 0.59408
# elev:forest             0.3908  0.185  2.1094 0.03491
# elev:I(forest^2)        0.4230  0.212  1.9959 0.04595
# I(elev^2):forest        0.4931  0.296  1.6677 0.09537
# I(elev^2):I(forest^2)  -0.9047  0.285 -3.1706 0.00152   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE      z           P(>|z|)
# year2001           -1.394  0.519 -2.687 0.007214644529635
# year2002           -3.310  0.735 -4.500 0.000006795368613
# year2003           -2.761  0.613 -4.505 0.000006648001148
# year2004           -1.866  0.484 -3.859 0.000113955248144
# year2005           -2.475  0.712 -3.475 0.000510372999729
# year2006           -4.046  1.201 -3.369 0.000753913214526
# year2007           -2.498  0.583 -4.288 0.000018061323416
# year2008           -9.989 66.099 -0.151 0.879881921115889
# year2009           -2.491  0.435 -5.721 0.000000010619107
# year2010           -2.464  0.537 -4.585 0.000004544270362
# year2011           -1.812  0.441 -4.105 0.000040415110956
# elev               -2.061  0.283 -7.293 0.000000000000303
# I(elev^2)           1.282  0.332  3.862 0.000112481455071
# forest              0.870  0.416  2.088 0.036762500298937
# I(forest^2)        -0.262  0.171 -1.525 0.127196921704892  # toss?
# elev:forest         0.926  0.248  3.730 0.000191662345016
# I(elev^2):forest   -1.509  0.338 -4.469 0.000007841975357   # Blocks things here

# Detection (logit-scale):
               # Estimate     SE       z                 P(>|z|)
# year2001       -0.36550 0.2011 -1.8178 0.069099052243009356777
# year2002       -0.25069 0.1706 -1.4694 0.141721715957158045374
# year2003        0.08040 0.1500  0.5360 0.591986317162055519780
# year2004       -0.45121 0.1546 -2.9178 0.003524783938836379381
# year2005        0.35589 0.1872  1.9012 0.057281199779643590664
# year2006       -1.08998 0.1805 -6.0379 0.000000001561620425943
# year2007       -0.30867 0.1724 -1.7902 0.073423301171369845242
# year2008       -0.09451 0.1612 -0.5864 0.557583053862158672942
# year2009       -1.21150 0.1337 -9.0646 0.000000000000000000125
# year2010        0.30487 0.1463  2.0841 0.037152154421218976099
# year2011       -0.74252 0.1518 -4.8911 0.000001002762840838343
# year2012        0.20101 0.1652  1.2164 0.223826497498105292383
# elev            0.57100 0.0995  5.7394 0.000000009499934856696
# I(elev^2)      -0.00754 0.0977 -0.0772 0.938465105971895585668
# forest          0.68385 0.0823  8.3111 0.000000000000000094853
# I(forest^2)    -0.29812 0.0504 -5.9109 0.000000003402537709987
# date            0.13805 0.0514  2.6845 0.007263211938880680986
# I(date^2)       0.15203 0.0509  2.9897 0.002792107736349363264
# elev:forest     0.39749 0.0708  5.6178 0.000000019338253978770
# elev:date      -0.01602 0.0826 -0.1939 0.846235928070424137104
# elev:I(date^2) -0.11208 0.0540 -2.0772 0.037784763301656948409   # Blocks things here

# AIC: 6126.118
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 163
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm46  66 3028.00        0.00    0.75   0.75 -1425.89
# fm45  67 3030.62        2.61    0.20   0.95 -1425.41
# fm44  68 3033.71        5.70    0.04   0.99 -1425.16
# fm43  69 3037.11        9.10    0.01   1.00 -1425.04
# fm42  70 3040.72       12.71    0.00   1.00 -1425.00
# fm41  71 3044.09       16.08    0.00   1.00 -1424.83
# fm40  72 3047.71       19.71    0.00   1.00 -1424.76
# fm39  73 3051.49       23.48    0.00   1.00 -1424.75
# fm38X 74 3055.30       27.30    0.00   1.00 -1424.74

# Try to toss out I(forest^2) in extinction
inits <- coef(fm46)[-42]
system.time(
   (fm47 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date -
      elev:I(forest^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2.4 mins
)
summary(fm47)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41,
    'fm42' = fm42, 'fm43' = fm43, 'fm44' = fm44, 'fm45' = fm45, 'fm46' = fm46,
    'fm47' = fm47), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2) -
    # elev:I(forest^2) - I(forest^2), pformula = ~(year - 1) +
    # (elev + I(elev^2)) * (forest + I(forest^2)) + date + I(date^2) +
    # date:elev + date:I(elev^2) + I(date^2):elev + I(date^2):I(elev^2) -
    # I(elev^2):I(date^2) - I(elev^2):I(forest^2) - I(elev^2):forest -
    # I(elev^2):date - elev:I(forest^2), data = umf, starts = inits,
    # se = TRUE, control = list(trace = TRUE, REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)       -0.0951 0.381 -0.250 0.8027439
# elev               1.9363 0.454  4.261 0.0000204
# I(elev^2)         -1.2242 0.510 -2.403 0.0162767
# forest            -0.6958 0.507 -1.372 0.1699840
# I(forest^2)        0.1626 0.260  0.624 0.5324325
# elev:forest        0.1687 0.450  0.374 0.7080779
# elev:I(forest^2)  -0.7770 0.354 -2.198 0.0279812
# I(elev^2):forest   1.8584 0.537  3.459 0.0005416   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z P(>|z|)
# year2001               -0.0852  0.326 -0.261 0.79382
# year2002               -0.3581  0.381 -0.940 0.34738
# year2003               -1.2720  0.511 -2.491 0.01273
# year2004               -1.5208  0.520 -2.926 0.00343
# year2005               -0.2623  0.425 -0.617 0.53743
# year2006               -1.1643  0.756 -1.541 0.12341
# year2007               -2.9921  2.325 -1.287 0.19806
# year2008                0.0633  0.400  0.158 0.87438
# year2009               -8.7211 51.394 -0.170 0.86525
# year2010               -0.7733  0.470 -1.646 0.09977
# year2011               -3.4418  2.374 -1.450 0.14713
# elev                    0.5019  0.208  2.410 0.01595
# I(elev^2)              -0.2259  0.291 -0.776 0.43803
# forest                 -0.3902  0.280 -1.392 0.16406
# I(forest^2)             0.1948  0.215  0.907 0.36418
# elev:forest             0.3688  0.184  2.006 0.04491
# elev:I(forest^2)        0.4273  0.213  2.003 0.04513
# I(elev^2):forest        0.5769  0.290  1.991 0.04652
# I(elev^2):I(forest^2)  -0.9392  0.288 -3.256 0.00113   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE      z           P(>|z|)
# year2001           -1.315  0.505 -2.603 0.009234078295600
# year2002           -3.244  0.705 -4.601 0.000004194979833
# year2003           -2.710  0.605 -4.482 0.000007384107179
# year2004           -1.831  0.462 -3.964 0.000073621696332
# year2005           -2.397  0.681 -3.522 0.000428936274346
# year2006           -4.159  1.130 -3.680 0.000232971211342
# year2007           -2.422  0.565 -4.288 0.000017991428265
# year2008           -9.986 57.475 -0.174 0.862069120260057
# year2009           -2.474  0.435 -5.691 0.000000012620461
# year2010           -2.429  0.530 -4.584 0.000004562921717
# year2011           -1.743  0.427 -4.082 0.000044697237011
# elev               -1.968  0.266 -7.403 0.000000000000134
# I(elev^2)           1.107  0.302  3.671 0.000241242635436
# forest              0.407  0.254  1.604 0.108804108735872
# elev:forest         0.857  0.223  3.838 0.000124226140721   # Blocks things here
# I(elev^2):forest   -1.293  0.282 -4.587 0.000004500569298   # Blocks things here

# Detection (logit-scale):
               # Estimate     SE      z  P(>|z|)
# year2001        -0.3622 0.2002 -1.809 7.04e-02
# year2002        -0.2314 0.1714 -1.350 1.77e-01
# year2003         0.0730 0.1488  0.491 6.23e-01
# year2004        -0.4432 0.1549 -2.861 4.22e-03
# year2005         0.3657 0.1849  1.977 4.80e-02
# year2006        -1.0797 0.1808 -5.973 2.33e-09
# year2007        -0.3182 0.1636 -1.945 5.17e-02
# year2008        -0.0807 0.1617 -0.499 6.18e-01
# year2009        -1.2071 0.1334 -9.052 1.40e-19
# year2010         0.3041 0.1472  2.066 3.88e-02
# year2011        -0.7363 0.1526 -4.826 1.39e-06
# year2012         0.2284 0.1661  1.375 1.69e-01
# elev             0.5626 0.0992  5.674 1.40e-08
# I(elev^2)       -0.0106 0.0979 -0.108 9.14e-01
# forest           0.6486 0.0794  8.172 3.02e-16
# I(forest^2)     -0.2784 0.0491 -5.670 1.43e-08
# date             0.1366 0.0514  2.655 7.93e-03
# I(date^2)        0.1504 0.0509  2.955 3.13e-03
# elev:forest      0.4090 0.0701  5.832 5.49e-09
# elev:date       -0.0145 0.0826 -0.176 8.61e-01
# elev:I(date^2)  -0.1112 0.0540 -2.061 3.93e-02  # toss?

# AIC: 6126.675
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 165
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm47  65 3025.69        0.00    0.70   0.70 -1426.50
# fm46  66 3028.00        2.32    0.22   0.92 -1425.89
# fm45  67 3030.62        4.93    0.06   0.98 -1425.41
# fm44  68 3033.71        8.02    0.01   1.00 -1425.16
# fm43  69 3037.11       11.42    0.00   1.00 -1425.04
# fm42  70 3040.72       15.03    0.00   1.00 -1425.00
# fm41  71 3044.09       18.40    0.00   1.00 -1424.83
# fm40  72 3047.71       22.03    0.00   1.00 -1424.76
# fm39  73 3051.49       25.80    0.00   1.00 -1424.75
# fm38X 74 3055.30       29.61    0.00   1.00 -1424.74

# QAICc continues to come down....

# Try to toss out elev:I(date^2) in detection
inits <- coef(fm47)[-64]
system.time(
   (fm48 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date -
      elev:I(forest^2) - I(date^2):elev,
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2.2 mins
)
summary(fm48)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41,
    'fm42' = fm42, 'fm43' = fm43, 'fm44' = fm44, 'fm45' = fm45, 'fm46' = fm46,
    'fm47' = fm47, 'fm48' = fm48), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2) -
    # elev:I(forest^2) - I(forest^2), pformula = ~(year - 1) +
    # (elev + I(elev^2)) * (forest + I(forest^2)) + date + I(date^2) +
    # date:elev + date:I(elev^2) + I(date^2):elev + I(date^2):I(elev^2) -
    # I(elev^2):I(date^2) - I(elev^2):I(forest^2) - I(elev^2):forest -
    # I(elev^2):date - elev:I(forest^2) - I(date^2):elev, data = umf,
    # starts = inits, se = TRUE, control = list(trace = TRUE, REPORT = 5,
        # maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)       -0.0837 0.383 -0.218 0.8270847
# elev               1.9446 0.457  4.259 0.0000206
# I(elev^2)         -1.2362 0.512 -2.414 0.0157725
# forest            -0.7070 0.510 -1.386 0.1656282
# I(forest^2)        0.1643 0.261  0.629 0.5294254
# elev:forest        0.1684 0.451  0.373 0.7089892
# elev:I(forest^2)  -0.7801 0.354 -2.203 0.0275747  # toss?
# I(elev^2):forest   1.8713 0.539  3.473 0.0005144   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z  P(>|z|)
# year2001               -0.0936  0.326 -0.287 0.774306
# year2002               -0.3677  0.382 -0.962 0.335856
# year2003               -1.2825  0.510 -2.515 0.011905
# year2004               -1.5032  0.510 -2.948 0.003200
# year2005               -0.2813  0.424 -0.663 0.507400
# year2006               -1.1352  0.733 -1.550 0.121217
# year2007               -3.0147  2.314 -1.303 0.192689
# year2008                0.0623  0.401  0.155 0.876451
# year2009               -8.7202 54.541 -0.160 0.872974
# year2010               -0.7853  0.470 -1.672 0.094620
# year2011               -3.4473  2.354 -1.464 0.143077
# elev                    0.4970  0.207  2.396 0.016564
# I(elev^2)              -0.2151  0.292 -0.737 0.461245
# forest                 -0.3922  0.280 -1.401 0.161178
# I(forest^2)             0.2002  0.214  0.937 0.349013
# elev:forest             0.3659  0.184  1.993 0.046314
# elev:I(forest^2)        0.4318  0.213  2.026 0.042736
# I(elev^2):forest        0.5757  0.289  1.990 0.046635
# I(elev^2):I(forest^2)  -0.9509  0.288 -3.304 0.000953   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE      z           P(>|z|)
# year2001           -1.317  0.507 -2.597 0.009403465618673
# year2002           -3.255  0.711 -4.578 0.000004691657785
# year2003           -2.723  0.607 -4.489 0.000007147280604
# year2004           -1.824  0.460 -3.968 0.000072344150220
# year2005           -2.390  0.669 -3.571 0.000356204914608
# year2006           -4.194  1.126 -3.724 0.000196069362362
# year2007           -2.423  0.566 -4.284 0.000018357507105
# year2008           -9.984 63.168 -0.158 0.874411302388129
# year2009           -2.488  0.435 -5.719 0.000000010731842
# year2010           -2.435  0.529 -4.604 0.000004150030181
# year2011           -1.753  0.427 -4.109 0.000039741730348
# elev               -1.973  0.263 -7.504 0.000000000000062
# I(elev^2)           1.123  0.305  3.687 0.000227309910074
# forest              0.417  0.254  1.645 0.099878452407150
# elev:forest         0.857  0.224  3.829 0.000128453832458
# I(elev^2):forest   -1.319  0.284 -4.640 0.000003476001683   # Blocks things here

# Detection (logit-scale):
            # Estimate     SE      z  P(>|z|)
# year2001    -0.35054 0.2007 -1.746 8.07e-02
# year2002    -0.21046 0.1713 -1.229 2.19e-01
# year2003     0.09577 0.1486  0.644 5.19e-01
# year2004    -0.42244 0.1541 -2.741 6.12e-03
# year2005     0.38631 0.1844  2.095 3.62e-02
# year2006    -1.04466 0.1796 -5.816 6.03e-09
# year2007    -0.29827 0.1616 -1.846 6.50e-02
# year2008    -0.05571 0.1618 -0.344 7.31e-01
# year2009    -1.18581 0.1331 -8.910 5.11e-19
# year2010     0.32182 0.1468  2.192 2.84e-02
# year2011    -0.71243 0.1522 -4.682 2.84e-06
# year2012     0.25260 0.1653  1.528 1.26e-01
# elev         0.47937 0.0903  5.311 1.09e-07
# I(elev^2)   -0.00342 0.0977 -0.035 9.72e-01  # toss?
# forest       0.64799 0.0793  8.169 3.11e-16
# I(forest^2) -0.27775 0.0492 -5.651 1.59e-08
# date         0.08378 0.0445  1.883 5.97e-02
# I(date^2)    0.12219 0.0491  2.487 1.29e-02
# elev:forest  0.41477 0.0702  5.907 3.48e-09
# elev:date   -0.02195 0.0829 -0.265 7.91e-01  # toss?

# AIC: 6128.881
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 139
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm48  64 3024.19        0.00    0.60   0.60 -1427.50
# fm47  65 3025.69        1.50    0.28   0.88 -1426.50
# fm46  66 3028.00        3.82    0.09   0.97 -1425.89
# fm45  67 3030.62        6.43    0.02   0.99 -1425.41
# fm44  68 3033.71        9.52    0.01   1.00 -1425.16
# fm43  69 3037.11       12.92    0.00   1.00 -1425.04
# fm42  70 3040.72       16.53    0.00   1.00 -1425.00
# fm41  71 3044.09       19.90    0.00   1.00 -1424.83
# fm40  72 3047.71       23.52    0.00   1.00 -1424.76
# fm39  73 3051.49       27.30    0.00   1.00 -1424.75
# fm38X 74 3055.30       31.11    0.00   1.00 -1424.74

# QAICc continues to come down....

# Try to toss out I(elev^2) in detection
inits <- coef(fm48)[-57]
system.time(
   (fm49 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date -
      elev:I(forest^2) - I(date^2):elev - I(elev^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 1.8 mins
)
summary(fm49)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41,
    'fm42' = fm42, 'fm43' = fm43, 'fm44' = fm44, 'fm45' = fm45, 'fm46' = fm46,
    'fm47' = fm47, 'fm48' = fm48, 'fm49' = fm49), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2) -
    # elev:I(forest^2) - I(forest^2), pformula = ~(year - 1) +
    # (elev + I(elev^2)) * (forest + I(forest^2)) + date + I(date^2) +
    # date:elev + date:I(elev^2) + I(date^2):elev + I(date^2):I(elev^2) -
    # I(elev^2):I(date^2) - I(elev^2):I(forest^2) - I(elev^2):forest -
    # I(elev^2):date - elev:I(forest^2) - I(date^2):elev - I(elev^2),
    # data = umf, starts = inits, se = TRUE, control = list(trace = TRUE,
        # REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)       -0.0839 0.383 -0.219 0.8264747
# elev               1.9449 0.455  4.276 0.0000191
# I(elev^2)         -1.2364 0.509 -2.427 0.0152081
# forest            -0.7070 0.510 -1.387 0.1653614
# I(forest^2)        0.1644 0.261  0.629 0.5290436
# elev:forest        0.1684 0.451  0.373 0.7088202
# elev:I(forest^2)  -0.7799 0.354 -2.203 0.0275753
# I(elev^2):forest   1.8711 0.539  3.474 0.0005124   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z  P(>|z|)
# year2001               -0.0935  0.326 -0.287 0.774222
# year2002               -0.3676  0.382 -0.963 0.335646
# year2003               -1.2826  0.510 -2.516 0.011880
# year2004               -1.5030  0.510 -2.948 0.003200
# year2005               -0.2815  0.424 -0.664 0.506704
# year2006               -1.1353  0.732 -1.552 0.120762
# year2007               -3.0148  2.303 -1.309 0.190597
# year2008                0.0624  0.400  0.156 0.876093
# year2009               -8.7202 54.849 -0.159 0.873679
# year2010               -0.7855  0.469 -1.674 0.094202
# year2011               -3.4474  2.334 -1.477 0.139657
# elev                    0.4976  0.206  2.413 0.015842
# I(elev^2)              -0.2156  0.291 -0.740 0.459364
# forest                 -0.3922  0.280 -1.402 0.161044
# I(forest^2)             0.2000  0.214  0.936 0.349444
# elev:forest             0.3653  0.182  2.002 0.045269
# elev:I(forest^2)        0.4321  0.213  2.029 0.042431
# I(elev^2):forest        0.5763  0.288  1.999 0.045645
# I(elev^2):I(forest^2)  -0.9509  0.288 -3.306 0.000947   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE      z            P(>|z|)
# year2001           -1.317  0.506 -2.604 0.0092056110225532
# year2002           -3.255  0.710 -4.582 0.0000045972186352
# year2003           -2.723  0.606 -4.494 0.0000069814427334
# year2004           -1.824  0.456 -3.996 0.0000643066794642
# year2005           -2.390  0.665 -3.595 0.0003244296574923
# year2006           -4.194  1.123 -3.736 0.0001871213046904
# year2007           -2.423  0.565 -4.286 0.0000181522324314
# year2008           -9.984 63.254 -0.158 0.8745813048576252
# year2009           -2.488  0.433 -5.747 0.0000000091071364
# year2010           -2.435  0.528 -4.613 0.0000039667257737
# year2011           -1.753  0.425 -4.120 0.0000378551355351
# elev               -1.973  0.262 -7.532 0.0000000000000501
# I(elev^2)           1.125  0.301  3.741 0.0001836260567897
# forest              0.417  0.253  1.646 0.0998330477638100
# elev:forest         0.856  0.224  3.831 0.0001274709056734
# I(elev^2):forest   -1.320  0.284 -4.643 0.0000034389083485   # Blocks things here

# Detection (logit-scale):
            # Estimate     SE      z  P(>|z|)
# year2001     -0.3511 0.1975 -1.778 7.54e-02
# year2002     -0.2121 0.1624 -1.306 1.92e-01
# year2003      0.0946 0.1436  0.659 5.10e-01
# year2004     -0.4230 0.1504 -2.812 4.92e-03
# year2005      0.3855 0.1788  2.156 3.11e-02
# year2006     -1.0463 0.1719 -6.086 1.16e-09
# year2007     -0.2994 0.1539 -1.946 5.17e-02
# year2008     -0.0564 0.1589 -0.355 7.22e-01
# year2009     -1.1867 0.1279 -9.280 1.69e-20
# year2010      0.3210 0.1422  2.257 2.40e-02
# year2011     -0.7135 0.1457 -4.898 9.66e-07
# year2012      0.2523 0.1631  1.547 1.22e-01
# elev          0.4780 0.0808  5.914 3.33e-09
# forest        0.6479 0.0791  8.187 2.68e-16
# I(forest^2)  -0.2779 0.0489 -5.686 1.30e-08
# date          0.0836 0.0445  1.879 6.03e-02
# I(date^2)     0.1223 0.0484  2.524 1.16e-02
# elev:forest   0.4154 0.0683  6.084 1.18e-09
# elev:date    -0.0230 0.0759 -0.303 7.62e-01  # toss?

# AIC: 6126.882
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 52
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm49  63 3020.73        0.00    0.77   0.77 -1427.50
# fm48  64 3024.19        3.46    0.14   0.91 -1427.50
# fm47  65 3025.69        4.96    0.06   0.97 -1426.50
# fm46  66 3028.00        7.28    0.02   0.99 -1425.89
# fm45  67 3030.62        9.89    0.01   1.00 -1425.41
# fm44  68 3033.71       12.98    0.00   1.00 -1425.16
# fm43  69 3037.11       16.38    0.00   1.00 -1425.04
# fm42  70 3040.72       19.99    0.00   1.00 -1425.00
# fm41  71 3044.09       23.36    0.00   1.00 -1424.83
# fm40  72 3047.71       26.99    0.00   1.00 -1424.76
# fm39  73 3051.49       30.76    0.00   1.00 -1424.75
# fm38X 74 3055.30       34.58    0.00   1.00 -1424.74

# QAICc continues to come down....

# Try to toss out elev:date in detection
inits <- coef(fm49)[-62]
system.time(
   (fm50 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date -
      elev:I(forest^2) - I(date^2):elev - I(elev^2) - elev:date,
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2 mins
)
summary(fm50)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41,
    'fm42' = fm42, 'fm43' = fm43, 'fm44' = fm44, 'fm45' = fm45, 'fm46' = fm46,
    'fm47' = fm47, 'fm48' = fm48, 'fm49' = fm49, 'fm50' = fm50), c.hat = c.hat)

# Call:
# colext(psiformula = ~(elev + I(elev^2)) * (forest + I(forest^2)) -
    # I(elev^2):I(forest^2), gammaformula = ~(year - 1) + (elev +
    # I(elev^2)) * (forest + I(forest^2)), epsilonformula = ~(year -
    # 1) + (elev + I(elev^2)) * (forest + I(forest^2)) - I(elev^2):I(forest^2) -
    # elev:I(forest^2) - I(forest^2), pformula = ~(year - 1) +
    # (elev + I(elev^2)) * (forest + I(forest^2)) + date + I(date^2) +
    # date:elev + date:I(elev^2) + I(date^2):elev + I(date^2):I(elev^2) -
    # I(elev^2):I(date^2) - I(elev^2):I(forest^2) - I(elev^2):forest -
    # I(elev^2):date - elev:I(forest^2) - I(date^2):elev - I(elev^2) -
    # elev:date, data = umf, starts = inits, se = TRUE, control = list(trace = TRUE,
    # REPORT = 5, maxit = 500))

# Initial (logit-scale):
                 # Estimate    SE      z   P(>|z|)
# (Intercept)       -0.0793 0.383 -0.207 0.8361599
# elev               1.9507 0.455  4.287 0.0000181
# I(elev^2)         -1.2445 0.510 -2.441 0.0146385
# forest            -0.7112 0.510 -1.393 0.1635247
# I(forest^2)        0.1655 0.261  0.633 0.5265466
# elev:forest        0.1665 0.451  0.369 0.7120662
# elev:I(forest^2)  -0.7809 0.354 -2.206 0.0273602
# I(elev^2):forest   1.8743 0.539  3.476 0.0005093   # Blocks things here

# Colonization (logit-scale):
                      # Estimate     SE      z  P(>|z|)
# year2001               -0.0926  0.326 -0.284 0.776591
# year2002               -0.3715  0.382 -0.972 0.330912
# year2003               -1.2843  0.509 -2.524 0.011609
# year2004               -1.5021  0.509 -2.953 0.003150
# year2005               -0.2807  0.425 -0.661 0.508650
# year2006               -1.1367  0.733 -1.552 0.120730
# year2007               -3.0584  2.363 -1.294 0.195499
# year2008                0.0608  0.400  0.152 0.879224
# year2009               -8.7195 55.857 -0.156 0.875951
# year2010               -0.7878  0.469 -1.680 0.093030
# year2011               -3.4583  2.350 -1.472 0.141121
# elev                    0.5002  0.206  2.428 0.015184
# I(elev^2)              -0.2159  0.291 -0.741 0.458651
# forest                 -0.3927  0.280 -1.402 0.160953
# I(forest^2)             0.2008  0.214  0.938 0.348005
# elev:forest             0.3629  0.182  1.992 0.046326
# elev:I(forest^2)        0.4318  0.213  2.029 0.042477
# I(elev^2):forest        0.5789  0.288  2.007 0.044709
# I(elev^2):I(forest^2)  -0.9509  0.288 -3.307 0.000945   # Blocks things here

# Extinction (logit-scale):
                 # Estimate     SE      z            P(>|z|)
# year2001           -1.322  0.506 -2.614 0.0089598376388647
# year2002           -3.255  0.707 -4.603 0.0000041731853998
# year2003           -2.729  0.606 -4.506 0.0000066051474440
# year2004           -1.832  0.458 -4.001 0.0000629802006562
# year2005           -2.398  0.662 -3.623 0.0002914545859385
# year2006           -4.198  1.118 -3.755 0.0001734660491413
# year2007           -2.430  0.565 -4.300 0.0000171120827765
# year2008           -9.983 61.748 -0.162 0.8715627422218372
# year2009           -2.494  0.433 -5.760 0.0000000084251477
# year2010           -2.437  0.527 -4.628 0.0000036848041496
# year2011           -1.759  0.425 -4.136 0.0000353733760410
# elev               -1.975  0.262 -7.552 0.0000000000000429
# I(elev^2)           1.132  0.300  3.776 0.0001592998002870
# forest              0.419  0.254  1.653 0.0982960408791153
# elev:forest         0.857  0.223  3.836 0.0001252123094098
# I(elev^2):forest   -1.322  0.284 -4.649 0.0000033432656187   # Blocks things here

# Detection (logit-scale):
            # Estimate     SE      z  P(>|z|)
# year2001     -0.3512 0.1977 -1.777 7.56e-02
# year2002     -0.2137 0.1623 -1.316 1.88e-01
# year2003      0.0972 0.1433  0.678 4.98e-01
# year2004     -0.4199 0.1501 -2.797 5.15e-03
# year2005      0.3857 0.1791  2.153 3.13e-02
# year2006     -1.0458 0.1718 -6.087 1.15e-09
# year2007     -0.2973 0.1532 -1.941 5.23e-02
# year2008     -0.0545 0.1586 -0.344 7.31e-01
# year2009     -1.1855 0.1278 -9.274 1.80e-20
# year2010      0.3231 0.1420  2.275 2.29e-02
# year2011     -0.7098 0.1449 -4.899 9.64e-07
# year2012      0.2540 0.1630  1.558 1.19e-01
# elev          0.4714 0.0774  6.086 1.16e-09
# forest        0.6493 0.0790  8.224 1.97e-16
# I(forest^2)  -0.2784 0.0488 -5.703 1.18e-08
# date          0.0816 0.0440  1.854 6.37e-02
# I(date^2)     0.1135 0.0383  2.962 3.06e-03  # toss?
# elev:forest   0.4188 0.0674  6.216 5.10e-10   # Blocks things here

# AIC: 6124.975
# Number of sites: 267
# optim convergence code: 0
# optim iterations: 98
# Bootstrap iterations: 0

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)

       # K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm50  62 3017.34        0.00    0.81   0.81 -1427.52
# fm49  63 3020.73        3.39    0.15   0.96 -1427.50
# fm48  64 3024.19        6.85    0.03   0.98 -1427.50
# fm47  65 3025.69        8.35    0.01   0.99 -1426.50
# fm46  66 3028.00       10.66    0.00   1.00 -1425.89
# fm45  67 3030.62       13.28    0.00   1.00 -1425.41
# fm44  68 3033.71       16.37    0.00   1.00 -1425.16
# fm43  69 3037.11       19.77    0.00   1.00 -1425.04
# fm42  70 3040.72       23.38    0.00   1.00 -1425.00
# fm41  71 3044.09       26.75    0.00   1.00 -1424.83
# fm40  72 3047.71       30.37    0.00   1.00 -1424.76
# fm39  73 3051.49       34.15    0.00   1.00 -1424.75
# fm38X 74 3055.30       37.96    0.00   1.00 -1424.74

# QAICc continues to come down....

# Try to toss out I(date^2) in detection
inits <- coef(fm50)[-60]
system.time(
   (fm51 <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
       I(elev^2):I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
      I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
    ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) +
      date + I(date^2) + date:elev + date:I(elev^2) +
      I(date^2):elev + I(date^2):I(elev^2) - I(elev^2):I(date^2) -
      I(elev^2):I(forest^2) - I(elev^2):forest - I(elev^2):date -
      elev:I(forest^2) - I(date^2):elev - I(elev^2) - elev:date -
      I(date^2),
      umf, starts = inits,
      control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))  # 2 mins
)
summary(fm51)
aictab(list('fm38X' = fm38X, 'fm39' = fm39, 'fm40' = fm40, 'fm41' = fm41,
    'fm42' = fm42, 'fm43' = fm43, 'fm44' = fm44, 'fm45' = fm45, 'fm46' = fm46,
    'fm47' = fm47, 'fm48' = fm48, 'fm49' = fm49, 'fm50' = fm50, 'fm51' = fm51),
    c.hat = c.hat)

# ~~~~~~~~ this table of results is in the book ~~~~~~~~~~~

# Model selection based on QAICc:
# (c-hat estimate = 2.102584)
#        K   QAICc Delta_QAICc QAICcWt Cum.Wt Quasi.LL
# fm50  62 3017.34        0.00    0.53   0.53 -1427.52
# fm51  61 3018.17        0.83    0.35   0.87 -1429.64
# fm49  63 3020.73        3.39    0.10   0.97 -1427.50
# fm48  64 3024.19        6.85    0.02   0.99 -1427.50
# fm47  65 3025.69        8.35    0.01   1.00 -1426.50
# [ ... list truncated ...]
# fm39  73 3051.49       34.15    0.00   1.00 -1424.75
# fm38X 74 3055.30       37.96    0.00   1.00 -1424.74

# ~~~~~~~~~~~~ more bonus code ~~~~~~~~~~~~~~~~~~~~~
# Now, QAICc comes up again. All the other critical terms in the model have p-values smaller than I(date^2), so it is unlikely that dropping them would lead to an improvement in QAICc. Hence, we assume that model 50 is the one with the best QAICc that can be identified by single-step changes in the model.

# Write the model formula more concisely
inits <- coef(fm50)
(fm50x <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
     I(elev^2):I(forest^2),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
    I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
  ~ (year-1) + elev + forest + I(forest^2) + date + I(date^2) +
    elev:forest, umf, starts = inits,
     control=list(trace=TRUE, REPORT=5, maxit = 500), se = TRUE))

summary(fm50x)
summary(lm(coef(fm50) ~coef(fm50x)))
# Seems virtually identical   ..... though why not exactly ???

# Rename 'original' fm50 - > fm50o
fm50o <- fm50

# Rename 'concise' fm50x - > fm50
fm50 <- fm50x

# Try model 50 once more without inits
inits <- NULL
 (fm50x <- colext(~ (elev + I(elev^2))* (forest + I(forest^2)) -
     I(elev^2):I(forest^2),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)),
  ~ (year-1) + (elev + I(elev^2))* (forest + I(forest^2)) -
    I(elev^2):I(forest^2) - elev:I(forest^2) - I(forest^2),
  ~ (year-1) + elev + forest + I(forest^2) + date + I(date^2) +
    elev:forest, umf, starts = inits,
     control=list(trace=TRUE, REPORT=25, maxit = 500), se = TRUE))

summary(fm50x)
summary(lm(coef(fm50) ~coef(fm50x)))
cbind('old' = coef(fm50), 'no inits' = coef(fm50x))
plot(coef(fm50), coef(fm50x), cex = 1.5, pch = 16) ; abline(0,1)

# ~~~ save output for use subsequently ~~~~~~~~~
save(fm50, umf, file = "AHM2_04.09.3_fm50.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
