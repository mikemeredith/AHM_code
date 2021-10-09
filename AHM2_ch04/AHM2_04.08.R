#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4 : MODELING SPECIES DISTRIBUTION AND RANGE DYNAMICS, AND POPULATION
#             DYNAMICS USING DYNAMIC OCCUPANCY MODELS
# ============================================================================
# Code from proofs dated 2020-08-18

# Run time with the full number of iterations: 32 mins

library(AHMbook)
library(AICcmodavg)
library(jagsUI)


# 4.8 Goodness of fit
# ===================

set.seed(68)
str(data <- simDynocc(nsites = 250, nyears = 10, nsurveys = 3,
    mean.psi1 = 0.4, beta.Xpsi1 = 1,
    range.phi = c(0.2, 1), beta.Xphi = 1,
    range.gamma = c(0, 0.8), beta.Xgamma = -1,
    range.p = c(0.1, 0.9), beta.Xp = 2,
    range.beta1.survey = c(0, 3), range.beta2.survey = c(-3, 0),
    trend.sd.site = c(0, 2), trend.sd.survey = c(0, 2)) )

# 4.8.1 MacKenzie & Bailey goodness of fit test for dynocc models
# --------------------------------------------------------------

# Fit constant dynocc model
library(unmarked) # Load the package
yy <- matrix(data$y, nrow = data$nsites, ncol = data$nsurveys *
    data$nyears)
summary(umf <- unmarkedMultFrame(y = yy, numPrimary = data$nyears))
summary(fm <- colext(~1, ~ 1, ~ 1, ~ 1, data = umf))

# Compute Chi-square test statistic for actual data by season and
# generate reference distribution of test statistic under H0 (~1.1 h)

# ~~~ reduce for testing ~~~~~~
# system.time(gof <- mb.gof.test(fm, print.table = FALSE, nsim = 1000,
#     plot.hist = TRUE, plot.seasons = TRUE, report = 1) ) # load(AICcmodavg)
system.time(gof <- mb.gof.test(fm, print.table = FALSE, nsim = 10,
    plot.hist = TRUE, plot.seasons = TRUE, report = 1, ncores = 3) )
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gof
# Goodness-of-fit for dynamic occupancy model

# Number of seasons: 10

# Chi-square statistic:
# Season 1 Season 2 Season 3 Season 4 Season 5
#  18.5134  73.8895  20.3665  54.8196  13.8125
# Season 6 Season 7 Season 8 Season 9 Season 10
#  84.9100  67.5141  67.5075  96.9058   56.9778

# Total chi-square = 555.2166
# Number of bootstrap samples = 1000
# P-value = 0

# Quantiles of bootstrapped statistics:
# 0% 25% 50% 75% 100%
# 39 58 66 74 117

# Estimate of c-hat = 8.37


# 4.8.2 GoF tests based on Bayesian posterior predictive distributions
# --------------------------------------------------------------------
# ~~~~~~~~~~~~~ code below from MS dated 2019-01-04 ~~~~~~~~~~~~~~~~~~~~

# Bundle and summarize data set
str(bdata <- list(y = data$y, nsite = data$nsite, nyear = data$nyear,
    nsurvey = data$nsurvey, Xpsi1 = data$Xpsi1, Xphi = data$Xphi,
    Xgamma = data$Xgamma, Xp = data$Xp, e = 0.0001))

# Specify model in BUGS language
cat(file = "dynocc.gof.txt","
  model {

  # (1) Model fitting part of code
  # ------------------------------
  # Priors
  alpha.lpsi <- logit(mean.psi1)
  mean.psi1 ~ dunif(0,1)
  beta.lpsi ~ dunif(-5, 5)
  alpha.lphi <- logit(mean.phi)
  mean.phi ~ dunif(0,1)
  alpha.lgamma <- logit(mean.gamma)
  mean.gamma ~ dunif(0,1)
  alpha.lp <- logit(mean.p)
  mean.p ~ dunif(0,1)
  beta.lphi ~ dunif(-5, 5)
  beta.lgamma ~ dunif(-5, 5)
  beta.lp ~ dunif(-5, 5)

  # Ecological submodel
  for (i in 1:nsite){
    # Initial conditions
    z[i,1] ~ dbern(psi1[i])
    logit(psi1[i]) <- alpha.lpsi + beta.lpsi * Xpsi1[i]
    psi[i, 1] <- psi1[i]     # Population occupancy as derived quantity
    # State transitions
    for (t in 2:nyear){
      z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1] + (1-z[i,t-1]) * gamma[i,t-1])
      logit(phi[i, t-1]) <- alpha.lphi + beta.lphi * Xphi[i, t-1]
      logit(gamma[i, t-1]) <- alpha.lgamma + beta.lgamma * Xgamma[i, t-1]
      psi[i, t]<-psi[i,t-1]*phi[i,t-1]+(1-psi[i,t-1])*gamma[i,t-1] # Derv'd
    }
  }

  # Observation model
  for (i in 1:nsite){
    for (j in 1:nsurvey){
      for (t in 1:nyear){
        y[i,j,t] ~ dbern(z[i,t] * p[i,j,t])
        logit(p[i,j,t]) <- alpha.lp + beta.lp * Xp[i,j,t]
      }
    }
  }

  # (2) GoF computation part of code
  # (based on posterior predictive distributions)
  # --------------------------------------------
  # Draw a replicate data set under the fitted model
  for (i in 1:nsite){
    for (t in 1:nyear){
      for (j in 1:nsurvey){
        yrep[i,j,t] ~ dbern(z[i,t] * p[i,j,t])
      }
    }
  }

  # (2a) Computations for the GoF of the open part of the model
  # (based on number of state transitions)
  # ----------------------------------------------------------
  # Compute observed z matrix for observed and replicated data
  for (i in 1:nsite){
    for (t in 1:nyear){
      zobs[i,t] <- max(y[i,,t])       # For observed data
      zobsrep[i,t] <- max(yrep[i,,t]) # For replicated data
    }

  # Identify extinctions, persistence, colonization and non-colonizations
  for (t in 2:nyear){
    # ... for observed data
    ext[i,(t-1)] <- equals(zobs[i,t],0) * equals(zobs[i,t-1],1)
    nonext[i,(t-1)] <- equals(zobs[i,t],1) * equals(zobs[i,t-1],1)
    colo[i,(t-1)] <- equals(zobs[i,t],1) * equals(zobs[i,t-1],0)
    noncolo[i,(t-1)] <- equals(zobs[i,t],0) * equals(zobs[i,t-1],0)
    # ... for replicated data
    extrep[i,(t-1)] <- equals(zobsrep[i,t],0) * equals(zobsrep[i,t-1],1)
    nonextrep[i,(t-1)] <- equals(zobsrep[i,t],1) * equals(zobsrep[i,t-1],1)
    colorep[i,(t-1)] <- equals(zobsrep[i,t],1) * equals(zobsrep[i,t-1],0)
    noncolorep[i,(t-1)] <- equals(zobsrep[i,t],0)*equals(zobsrep[i,t-1],0)
    }
  }

  # Tally up number of transitions and put into a matrix for each year
  for(t in 1:(nyear-1)){
    # ... for observed data
    tm[1,1,t] <- sum(noncolo[,t]) # transition mat for obs. data
    tm[1,2,t] <- sum(colo[,t])
    tm[2,1,t] <- sum(ext[,t])
    tm[2,2,t] <- sum(nonext[,t])
    # ... for replicated data
    tmrep[1,1,t] <- sum(noncolorep[,t]) # transition mat for rep. data
    tmrep[1,2,t] <- sum(colorep[,t])
    tmrep[2,1,t] <- sum(extrep[,t])
    tmrep[2,2,t] <- sum(nonextrep[,t])
  }

  # Compute expected numbers of transitions under the model
  # Probability of each individual transition
  for(i in 1:nsite){
   for(t in 1:(nyear-1)){
    noncolo.exp[i,t] <- (1-psi[i,t]) * (1-gamma[i,t])
    colo.exp[i,t] <- (1-psi[i,t]) * gamma[i,t]
    ext.exp[i,t] <- psi[i,t] * (1-phi[i,t])
    nonext.exp[i,t] <- psi[i,t] * phi[i,t]
   }
  }
  # Sum up over sites to obtain the expected number of those transitions
  for(t in 1:(nyear-1)){
    Etm[1,1,t] <- sum(noncolo.exp[,t])
    Etm[1,2,t] <- sum(colo.exp[,t])
    Etm[2,1,t] <- sum(ext.exp[,t])
    Etm[2,2,t] <- sum(nonext.exp[,t])
  }

  # Compute Chi-square discrepancy ~~~ see Errata 2021-10-09
  for(t in 1:(nyear-1)){
    # ... for observed data
    x2Open[1,1,t] <- pow((tm[1,1,t] - Etm[1,1,t]), 2) / (Etm[1,1,t]+e)
    x2Open[1,2,t] <- pow((tm[1,2,t] - Etm[1,2,t]), 2) / (Etm[1,2,t]+e)
    x2Open[2,1,t] <- pow((tm[2,1,t] - Etm[2,1,t]), 2) / (Etm[2,1,t]+e)
    x2Open[2,2,t] <- pow((tm[2,2,t] - Etm[2,2,t]), 2) / (Etm[2,2,t]+e)
    # ... for replicated data
    x2repOpen[1,1,t] <- pow((tmrep[1,1,t]-Etm[1,1,t]),2)/(Etm[1,1,t]+e)
    x2repOpen[1,2,t] <- pow((tmrep[1,2,t]-Etm[1,2,t]),2)/(Etm[1,2,t]+e)
    x2repOpen[2,1,t] <- pow((tmrep[2,1,t]-Etm[2,1,t]),2)/(Etm[2,1,t]+e)
    x2repOpen[2,2,t] <- pow((tmrep[2,2,t]-Etm[2,2,t]),2)/(Etm[2,2,t]+e)
  }

  # Add up overall test statistic and compute fit stat ratio (open part)
  Chi2Open <- sum(x2Open[,,])       # Chisq. statistic for observed data
  Chi2repOpen <- sum(x2repOpen[,,]) # Chisq. statistic for replicated data
  Chi2ratioOpen <- Chi2Open / Chi2repOpen


  # (2b) Computations for the GoF of the closed part of the model
  # (based on the number of times detected, i.e., detection freqiencies)
  # --------------------------------------------------------------------
  # Compute detection frequencies for observed and replicated data
  for (i in 1:nsite){
    for (t in 1:nyear){
      # Det. frequencies for observed and replicated data
      detfreq[i,t] <- sum(y[i,,t])
      detfreqrep[i,t] <- sum(yrep[i,,t])

      # Expected detection frequencies under the model
      for (j in 1:nsurvey){
        tmp[i,j,t] <- z[i, t] * p[i,j,t]
      }
      E[i,t] <- sum(tmp[i,,t])     # Expected number of detections

      # Chi-square and Freeman-Tukey discrepancy measures
      # ..... for actual data set
      x2Closed[i,t] <- pow((detfreq[i,t] - E[i,t]),2) / (E[i,t]+e)
      ftClosed[i,t] <- pow((sqrt(detfreq[i,t]) - sqrt(E[i,t])),2)

      # ..... for replicated data set
      x2repClosed[i,t] <- pow((detfreqrep[i,t] - E[i,t]),2) / (E[i,t]+e)
      ftrepClosed[i,t] <- pow((sqrt(detfreqrep[i,t]) - sqrt(E[i,t])),2)
    }
  }

  # Add up Chi-square and FT discrepancies and compute fit stat ratio (closed part)
  Chi2Closed <- sum(x2Closed[,])
  FTClosed <- sum(ftClosed[,])
  Chi2repClosed <- sum(x2repClosed[,])
  FTrepClosed <- sum(ftrepClosed[,])
  Chi2ratioClosed <- Chi2Closed / Chi2repClosed
  FTratioClosed <- FTClosed / FTrepClosed
}
")

# Initial values
inits <- function(){list(z = array(1, dim=c(data$nsite, data$nyear)))}

# Parameters monitored
params <- c('mean.psi1', 'beta.lpsi', 'mean.phi', 'beta.lphi', 'mean.gamma',
    'beta.lgamma', 'mean.p', 'beta.lp', 'Chi2Open', 'Chi2repOpen',
    'Chi2ratioOpen', 'Chi2Closed', 'Chi2repClosed', 'Chi2ratioClosed',
    'FTClosed', 'FTrepClosed', 'FTratioClosed', 'tm', 'tmrep', 'Etm')

# MCMC settings
na <- 1000  ;  ni <- 2000  ;  nt <- 2  ;  nb <- 1000  ;  nc <- 3

# Call JAGS from R, check convergence and summarize posteriors
out.gof <- jags(bdata, inits, params, "dynocc.gof.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out.gof)
print(out.gof, 2)

# JAGS output for model 'dynocc.gof.txt', generated by jagsUI.
# Estimates based on 3 chains of 2000 iterations,
# adaptation = 1000 iterations (sufficient),
# burn-in = 1000 iterations and thin rate = 2,
# yielding 1500 total samples from the joint posterior.
# MCMC ran in parallel for 3.385 minutes at time 2018-12-06 15:03:11.

                   # mean    sd    2.5%     50%   97.5% overlap0 f Rhat n.eff
# mean.psi1          0.44  0.04    0.36    0.44    0.53    FALSE 1 1.00  1500
# beta.lpsi          1.11  0.18    0.76    1.11    1.48    FALSE 1 1.00   653
# mean.phi           0.45  0.02    0.41    0.45    0.49    FALSE 1 1.00  1500
# beta.lphi          0.73  0.08    0.58    0.72    0.88    FALSE 1 1.01   295
# mean.gamma         0.48  0.02    0.44    0.48    0.52    FALSE 1 1.00   569
# beta.lgamma       -0.76  0.08   -0.93   -0.76   -0.62    FALSE 1 1.00  1500
# mean.p             0.37  0.01    0.34    0.37    0.39    FALSE 1 1.00  1500
# beta.lp            1.28  0.05    1.19    1.28    1.38    FALSE 1 1.00   620
# Chi2Open        1142.14 81.12  995.07 1139.76 1301.98    FALSE 1 1.00   555
# Chi2repOpen      489.98 78.25  354.01  484.31  658.66    FALSE 1 1.00  1500
# Chi2ratioOpen      2.38  0.33    1.77    2.36    3.10    FALSE 1 1.00  1351
# Chi2Closed       608.89 35.01  543.96  607.09  683.72    FALSE 1 1.00  1500
# Chi2repClosed    542.39 31.90  486.15  541.06  608.17    FALSE 1 1.00   981
# Chi2ratioClosed    1.12  0.06    1.01    1.12    1.25    FALSE 1 1.00  1500
# FTClosed         273.26 18.02  241.01  272.26  312.18    FALSE 1 1.00  1098
# FTrepClosed      263.69 16.20  233.63  263.45  296.72    FALSE 1 1.00  1500
# FTratioClosed      1.04  0.06    0.92    1.03    1.17    FALSE 1 1.00   696
# tm[1,1,1]         72.00  0.00   72.00   72.00   72.00    FALSE 1   NA     1
# tm[2,1,1]         56.00  0.00   56.00   56.00   56.00    FALSE 1   NA     1
# tm[1,2,1]         83.00  0.00   83.00   83.00   83.00    FALSE 1   NA     1
      # [  ... output truncated ...  ]

# Summary of test results for the open and the closed parts
op <- par(mfrow = c(1,3))

# Plots of expected versus observed value of fit stats
# Open part
pl <- range(c(out.gof$sims.list$Chi2Open, out.gof$sims.list$Chi2repOpen))
plot(out.gof$sims.list$Chi2Open, out.gof$sims.list$Chi2repOpen,
    xlab = "Chi2 observed data", ylab = "Chi2 expected data",
    main = "Open part of model ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(450, 1430, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repOpen >
    out.gof$sims.list$Chi2Open), 2)), cex = 2)

# Closed part of model: Chi-squared
pl <- range(c(out.gof$sims.list$Chi2Closed, out.gof$sims.list$Chi2repClosed))
plot(out.gof$sims.list$Chi2Closed, out.gof$sims.list$Chi2repClosed,
    xlab = "Chi2 observed data", ylab = "Chi2 expected data",
    main = "Closed part of model (Chi-squared)", xlim = pl, ylim = pl,
    frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(525, 720, paste('Bpv = ', round(mean(out.gof$sims.list$Chi2repClosed >
    out.gof$sims.list$Chi2Closed), 2)), cex = 2)

# Closed part of model: Freeman-Tukey
pl <- range(c(out.gof$sims.list$FTClosed, out.gof$sims.list$FTrepClosed))
plot(out.gof$sims.list$FTClosed, out.gof$sims.list$FTrepClosed,
    xlab = "FT observed data", ylab = "FT expected data",
    main = "Closed part of model (Freeman-Tukey)", xlim = pl, ylim = pl,
    frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(240, 335, paste('Bpv = ', round(mean(out.gof$sims.list$FTrepClosed >
    out.gof$sims.list$FTClosed), 2)), cex = 2)
par(op)

# Inspect observed and expected numbers of transitions for each interval
obs.trans <- out.gof$mean$tm
rep.trans <- out.gof$mean$tmrep
exp.trans <- out.gof$mean$Etm
dimnames(obs.trans) <- dimnames(rep.trans) <- dimnames(exp.trans) <-
    list(c("From Non-occ", "From Occ"), c("To Non-occ", "To Occ"))

for(t in 1:(data$nyear-1)){
  cat("\n*** Interval", t, " ***\n")
  cat("* Observed transitions*\n")
  print(obs.trans[,,t])

  cat("\n* Replicate data transitions*\n")
  print(rep.trans[,,t])

  cat("\n* Expected transitions*\n")
  print(exp.trans[,,t])
}
# ~~~~~~~~~ following is from book ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# *** Interval 1 ***
# * Observed transitions*
# To Non-occ To Occ
# From Non-occ 72 83
# From Occ 56 39
# * Replicate data transitions*
# To Non-occ To Occ
# From Non-occ 83.29733 74.59733
# From Occ 57.13000 34.97533
# * Expected transitions*
# To Non-occ To Occ
# From Non-occ 68.87431 66.69349
# From Occ 61.81234 52.61985
# [ .... truncated output ...]
