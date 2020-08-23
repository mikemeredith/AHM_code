#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from MS dated 2018-11-30

# 1.7 “Demographic” state-space models for inference about relative abundance
# ===========================================================================

# 1.7.3 Comments on demographic state-space models
# -----------------------------------------------------------

# ~~~~ code for figure 1.14 ~~~~~~~~~~
# Compare GLM, GLMM, Gaussian and demographic SSM in terms of standardized
#    population trajectory

# Load saved output
load("AHM2_01.04_out1.RData")
load("AHM2_01.05.3_out4.RData")
load("AHM2_01.07.2extra_outs.RData")
load("AHM2_01.07.2_out11.RData")

# Standardize to estimate in 2007
tmp1 <- out1$sims.list$popindex / out1$sims.list$popindex[,9] # Model 1
tmp4 <- out4$sims.list$popindex / out4$sims.list$popindex[,9] # Model 4
tmp8 <- out8$sims.list$popindex / out8$sims.list$popindex[,9] # Model 8
tmp11 <- out11$sims.list$popindex / out11$sims.list$popindex[,9] # Model 11
# Get posterior means and CRIs
pm <- cbind(apply(tmp1, 2, mean), apply(tmp4, 2, mean), apply(tmp8, 2, mean),
    apply(tmp11, 2, mean))
fn1 <- function(x) quantile(x, 0.025)
fn2 <- function(x) quantile(x, 0.975)
LCL <- cbind(apply(tmp1, 2, fn1), apply(tmp4, 2, fn1), apply(tmp8, 2,fn1),
    apply(tmp11, 2, fn1))
UCL <- cbind(apply(tmp1, 2, fn2), apply(tmp4, 2, fn2), apply(tmp8, 2, fn2),
  apply(tmp11, 2, fn2))

ylim = c(0.5, 1.3)
off <- 0.05
year <- 1999:2016

plot(year-3*off, pm[,1], type = 'b', pch = 16, lty = 1, xlab = 'Year', col = 1,
    ylim = ylim, ylab = 'Standardized population level', las = 1, frame = FALSE)
points(year-1*off, pm[,2], type = 'b', pch = 16, lty = 1, col = 2)
points(year+1*off, pm[,3], type = 'b', pch = 16, lty = 1, col = 3)
points(year+3*off, pm[,4], type = 'b', pch = 16, lty = 1, col = 4)
abline(h = 1, col = 'grey')
abline(v = 2007, col = 'grey')
segments(year-3*off, LCL[,1], year-3*off, UCL[,1], col = 1)
segments(year-1*off, LCL[,2], year-1*off, UCL[,2], col = 2)
segments(year+1*off, LCL[,3], year+1*off, UCL[,3], col = 3)
segments(year+3*off, LCL[,4], year+3*off, UCL[,4], col = 4)
legend('bottomright', c('GLM (Model 1)', 'GLMM (Model 4)',
    'Gaussian SSM (Model 8)', 'Demographic SSM (Model 11)'),
    col=1:4, pch=16, bty = 'n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
