#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# Approximate execution time for this code: 1.8 hrs
# Run time with the full number of iterations: 4.8 days

library(AHMbook)
library(unmarked)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14/"          # Place where your WinBUGS installed

# ~~~~~ this section needs 'umf' and 'fm5ZIP' from section 6.9 ~~~~~~~~~~~~~~~~~~~~~
data(SwissTits)
# Select Great tit and covariate data from 2013 and
#   drop 4 sites not surveyed in 2013
y0 <- SwissTits$counts[, , '2013', 'Great tit']
( NA.sites <- which(rowSums(is.na(y0)) == 3) ) # Unsurveyed sites
y <- y0[-NA.sites, ]                 # Drop them from the count data
tits <- SwissTits$sites[-NA.sites, ] # Also drop from the site covariates
elev.mean <- mean(tits$elev) ; elev.sd <- sd(tits$elev)
forest.mean <- mean(tits$forest) ; forest.sd <- sd(tits$forest)
# Get date and duration data for 2013, without the NA.sites rows:
date <- SwissTits$date[-NA.sites, , '2013']
dur <- SwissTits$dur[-NA.sites, , '2013']
time <- matrix(rep(as.character(1:3), nrow(y)), ncol = 3, byrow = TRUE)
umf <- unmarkedFramePCount(y = y,
  siteCovs=data.frame(elev=scale(tits[,"elev"]), forest=scale(tits[,"forest"]), iLength=1/tits[,"rlength"]),
  obsCovs=list(time = time, date = scale(date), dur = scale(dur)))

fm5ZIP <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
      - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
      - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
      - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
      - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2)
      - elev:I(date^2) - I(date^2):dur
      ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength
      - I(elev^2):forest - I(elev^2):I(forest^2),
      starts = c(2.67, -1.21, -0.27,  0.08, -0.17, -1.59, -0.03, -0.08, -0.14, -0.34, -0.29,  0.13,
        0.10,  0.03,  0.53,  0.12, -0.20, -0.55, -0.30,  0.03,  0.14,  0.28, -0.04, -0.09,  0),
      umf, control=list(trace=TRUE, REPORT=5), mixture = "ZIP")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 6.11 Bayesian modeling of Swiss Great tits with BUGS
# ====================================================


# 6.11.1 Bayesian fitting of the basic ZIP N-mixture model
# ------------------------------------------------------------------------
# Prepare data for BUGS data bundle
elev <- umf@siteCovs$elev   ;   elev2 <- elev^2
forest <- umf@siteCovs$forest   ;   forest2 <- forest^2
date <- matrix(umf@obsCovs$date, ncol = 3, byrow = TRUE)
dur <- matrix(umf@obsCovs$dur, ncol = 3, byrow = TRUE)
date[is.na(date)] <- 0   ;   date2 <- date^2
dur[is.na(dur)] <- 0   ;   dur2 <- dur^2
iRoute <- umf@siteCovs$iLength

# Design matrix for abundance model (no intercept)
lamDM <- model.matrix(~ elev + elev2 + forest + forest2 + elev:forest +
    elev:forest2 + iRoute)[,-1]


# Specify model in BUGS language
sink("ZIPNmix.txt")
cat("
model {

  # Specify priors
  # zero-inflation/suitability
  phi ~ dunif(0,1)          # proportion of suitable sites
  theta <- 1-phi            # zero-inflation (proportion of unsuitable)
  ltheta <- logit(theta)

  # abundance
  beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
  for(k in 1:7){            # Regression params in lambda
     beta[k] ~ dnorm(0, 1)
  }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ dunif(0, 2)      # site heterogeneity in lambda

  # detection
  for(j in 1:3){
     alpha0[j] <- logit(mean.p[j])
     mean.p[j] ~ dunif(0, 1)# p intercept for occasions 1-3
  }
  for(k in 1:13){           # Regression params in p
     alpha[k] ~ dnorm(0, 1)
  }
  tau.p.site <- pow(sd.p.site, -2)
  sd.p.site ~ dunif(0, 2)   # site heterogeneity in p
  tau.p.survey <- pow(sd.p.survey, -2)
  sd.p.survey ~ dunif(0, 2) # site-survey heterogeneity in p

  # ZIP model for abundance
  for (i in 1:nsite){
     a[i] ~ dbern(phi)
     eps.lam[i] ~ dnorm(0, tau.lam)       # Random site effects in log(abundance)
     loglam[i] <- beta0 + inprod(beta[], lamDM[i,]) + eps.lam[i] * hlam.on
     loglam.lim[i] <- min(250, max(-250, loglam[i]))  # 'Stabilize' log
     lam[i] <- exp(loglam.lim[i])
     mu.poisson[i] <- a[i] * lam[i]
     N[i] ~ dpois(mu.poisson[i])
  }

  # Measurement error model
  for (i in 1:nsite){
    eps.p.site[i] ~ dnorm(0, tau.p.site) # Random site effects in logit(p)
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # 'Stabilize' logit
      lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
        alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
        alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
        alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
        alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
        alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
        alpha[13] * date[i,j] * dur2[i,j] +
        eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
        eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
     }
  }
  # Posterior predictive distributions of chi2 discrepancy
  for (i in 1:nsite) {
    for (j in 1:nrep) {
      y.sim[i,j] ~ dbin(p[i,j], N[i]) # Create new data set under model
      e.count[i,j] <- N[i] * p[i,j]   # Expected datum
      # Chi-square discrepancy for the actual data
      chi2.actual[i,j] <- pow((y[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
      # Chi-square discrepancy for the simulated ('perfect') data
      chi2.sim[i,j] <- pow((y.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
      # Add small value e to denominator to avoid division by zero
    }
  }
  # Add up individual chi2 values for overall fit statistic
  fit.actual <- sum(chi2.actual[,])  # Fit statistic for actual data set
  fit.sim <- sum(chi2.sim[,])        # Fit statistic for a fitting model
  bpv <- step(fit.sim-fit.actual)    # Bayesian p-value
  c.hat <- fit.actual/fit.sim        # c-hat estimate

  # Derived parameters: Total abundance at 263 sampled sites
  Ntotal263 <- sum(N[])
}
",fill = TRUE)
sink()

# Initial values
Nst <- apply(y, 1, max, na.rm = T) + 1
Nst[is.na(Nst)] <- round(mean(y, na.rm = TRUE))
Nst[Nst == "-Inf"] <- round(mean(y, na.rm = TRUE))
inits <- function(){ list(N = Nst, beta0 = 0, mean.p = rep(0.5,3),
    beta = runif(7, 0,0), alpha = runif(13, 0,0))}

# Parameters monitored
params <- c("theta", "ltheta", "phi", "beta0", "beta", "sd.lam", "alpha0",
    "mean.p", "alpha", "sd.p.site", "sd.p.survey", "fit.actual", "fit.sim",
    "bpv", "c.hat", "Ntotal263")

# Bundle data and choose to fit simple ZIP model (model 1)
win.data1 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 0, hp.site.on = 0,
   hp.survey.on = 0)

# MCMC settings
# ni <- 50000    ;    nt <- 4    ;    nb <- 10000    ;    nc <- 3
ni <- 5000    ;    nt <- 1    ;    nb <- 1000    ;    nc <- 3  # ~~~~ for testing

# Call WinBUGS from R (ART 93 min) and summarize posteriors
out1 <- bugs(win.data1, inits, params, "ZIPNmix.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)  #  ~~~ for autotesting
print(out1, dig = 3)


# Compare MLEs and Bayesian posterior means (order first): table and graph
tmp <- summary(fm5ZIP)
ord.MLE <- rbind(tmp$psi[,1:2], tmp$state[,1:2], tmp$det[c(7:9, 1:6, 10:16),1:2])
ord.Bayes <- out1$summary[-c(1,3,12,16:18,32:39), 1:2]
cbind(ord.MLE, ord.Bayes)
op <- par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(ord.MLE[,1], ylim = c(-3,3), pch = 16, col = "black", main = "",
    frame = FALSE, xlab = "Parameters (zero-inflation, abundance, detection)",
    ylab = "Parameter estimate (+/- 1 SE)", cex = 1.5)
segments(1:25, ord.MLE[,1]-ord.MLE[,2], 1:25, ord.MLE[,1]+ord.MLE[,2], lwd = 2)
abline(h = 0)
abline(v = c(1.5, 9.5), col = "grey")
points((1:25)+0.3, ord.Bayes[,1], pch = 16, col = "blue", cex = 1.5)
segments(1:25+0.3, ord.Bayes[,1]-ord.Bayes[,2], 1:25+0.3,
    ord.Bayes[,1]+ord.Bayes[,2], col = "blue", lwd = 2)
par(op)

library(unmarked)
data(Switzerland)             # Load Swiss landscape data in unmarked
CH <- Switzerland


ELEV <- (CH$elev-elev.mean)/elev.sd
ELEV2 <- ELEV^2
FOREST <- (CH$forest-forest.mean)/forest.sd
FOREST2 <- FOREST^2
CHdata <- cbind(elev = ELEV, elev2 = ELEV2, forest = FOREST, forest2 = FOREST^2,
    iRoute = rep(0, length(CH$elev)), elev.forest = ELEV * FOREST,
    elev.forest2 = ELEV * FOREST2)
str(CHdata)                  # This is a design matrix

MCMCout <- out1              # Choose results output from model 1
(nsamp <- length(MCMCout$sims.list$theta))  # how many MCMC samples do we have ?

# Subsample
sub.sample.size <- 3000      # choose sample of 3000
selection <- sort(sample(1:nsamp, sub.sample.size))

# Array to hold predictions for every Swiss 1km2 quadrat
lamPred <- array(NA, dim =c(length(CH[,1]), sub.sample.size))

# Fill the array
for(i in 1:sub.sample.size){
  MCMCstep <- selection[i]
  lamPred[,i] <- (1-MCMCout$sims.list$theta[MCMCstep]) *
      exp(MCMCout$sims.list$beta0[MCMCstep] +
      CHdata %*% MCMCout$sims.list$beta[MCMCstep,1:7])
}

# Get posterior means for every quadrat and check if sensible
meanlam <- apply(lamPred, 1, mean)   # Get posterior mean
max(meanlam)                         # Check maximum
sum(meanlam > 100)                   # Are any predictions >100 ?
plot(CH$elev, meanlam, ylim = c(0, max(meanlam))) ; abline(v=2250, col="red")
plot(CH$elev, meanlam, ylim = c(0,100)) ; abline(v = 2250, col = "red", lwd = 2)
meanlam[meanlam > 100] <- 100        # censor all predictions
lamPred[lamPred > 100] <- 100

# Produce map of posterior mean of lambda
library(raster)
# library(rgdal)  # ~~~~ not necessary ~~~~
r <- rasterFromXYZ(data.frame(x = CH$x, y = CH$y, z = meanlam))
elevation <- rasterFromXYZ(cbind(CH$x, CH$y,CH$elevation))
elevation[elevation > 2250] <- NA
r <- mask(r, elevation)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
op <- par(mfrow = c(1,2), mar = c(5,5,1,5))
plot(r, col = mapPalette(100), axes = FALSE, box = FALSE, main ="")
# ~~~~ shape files not available ~~~~~~~~~~~~~~
# lakes <- readOGR(".", "lakes")
# rivers <- readOGR(".", "rivers")
# border <- readOGR(".", "border")
# plot(rivers, col = "dodgerblue", add = TRUE)
# plot(border, col = "transparent", lwd = 1.5, add = TRUE)
# plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

elev.class <- 100*(CH$elev %/% 100 + 1)      # elevation class of each km2
tmp <- aggregate(lamPred, by = list(elev.class), FUN = sum)
N.elev <- as.matrix(tmp[,-1])              # Posterior sample of Ntotal per band
band <- tmp[,1]                            # elevation band (in m)
meanN <- apply(N.elev, 1, mean)
barplot(meanN, col = "grey", horiz = T, xlab = "Number of Great tit territories",
    ylab = "Elevation band (100m)", xlim = c(0, 200000))
axis(2, at = 1:length(band), labels = band)
par(op)

# Posterior distribution of total number of great tit territories in 2013
keep <- which((CH$water < 50) & (CH$elev < 2251))
Ntot <- apply(lamPred[keep,], 2, sum)
hist(Ntot, breaks = 100, col = "grey",
    main = "Posterior of national population size")

# Point estimate and 95% CRI
mean(Ntot)
quantile(Ntot, prob = c(0.025, 0.975))


# 6.11.2 Adding random effects in BUGS
# ------------------------------------------------------------------------

# 6.11.2.1 Accounting for overdispersion at multiple scales
# ------------------------------------------------------------------------
# MCMC settings
# ni <- 10^6    ;    nt <- 80    ;    nb <- 200000    ;    nc <- 3
ni <- 12000    ;    nt <- 1     ;    nb <- 2000    ;    nc <- 3  # ~~~~ for testing

# Bundle data and select model 2
win.data2 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 0,
   hp.survey.on = 0)

# Call WinBUGS from R (ART 4050 min ~ 3 days) and summarize posteriors
out2 <- bugs(win.data2, inits, params, "ZIPNmix.txt", n.chains = nc, n.thin = nt,
    n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
# ~~~ After running for 3 days, you should probably save the result ~~~~~~~~~
save(out2, file="AHM1_06.11_out2.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(out2, dig = 3)


# Bundle data and select model 3
win.data3 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 0, hp.site.on = 1,
   hp.survey.on = 0)

# Call WinBUGS from R (ART 4200 min) and summarize posteriors
out3 <- bugs(win.data3, inits, params, "ZIPNmix.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
# ~~~ After running for 3 days, you should probably save the result ~~~~~~~~~~~
save(out3, file="AHM1_06.11_out3.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(out3, dig = 3)


# Bundle data and select model 4
win.data4 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 0, hp.site.on = 0,
   hp.survey.on = 1)

# Call WinBUGS from R (ART 4020 min) and summarize posteriors
out4 <- bugs(win.data4, inits, params, "ZIPNmix.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
# ~~~ After running for 3 days, you should probably save the result ~~~~~~~~~
save(out4, file="AHM1_06.11_out4.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(out4, dig = 3)


# Bundle data and select model 5
win.data5 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 1,
   hp.survey.on = 0)

# Call WinBUGS from R (ART 4250 min) and summarize posteriors
out5 <- bugs(win.data5, inits, params, "ZIPNmix.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
# ~~~ you should probably save the result ~~~~~~~~~~~~~~~~~
save(out5, file="AHM1_06.11_out5.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(out5, dig = 3)


# Bundle data and select model 6
win.data6 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 0,
   hp.survey.on = 1)

# Call WinBUGS from R (ART 4230 min) and summarize posteriors
out6 <- bugs(win.data6, inits, params, "ZIPNmix.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
# ~~~ After running for 3 days, you should probably save the result:
save(out6, file="AHM1_06.11_out6.RData")
print(out6, dig = 3)


# Bundle data and select model 7
win.data7 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
   date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 1,
   hp.survey.on = 1)

# Call WinBUGS from R (ART 4625 min) and summarize posteriors
out7 <- bugs(win.data7, inits, params, "ZIPNmix.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir)
# ~~~ you should probably save the result ~~~~~~~~~~~~
save(out7, file="AHM1_06.11_out7.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(out7, dig = 3)

# Look at posteriors for random effects
MCMCout <- out7       # Choose which model you want to plot random effects
op <- par(mfrow = c(1,3))
hist(MCMCout$sims.list$sd.lam, breaks = 60, col = "grey")
hist(MCMCout$sims.list$sd.p.site, breaks = 60, col = "grey")
hist(MCMCout$sims.list$sd.p.survey, breaks = 60, col = "grey")
par(op)

# 6.11.2.2 Linear modeling of a variance in the N-mixture model
# ------------------------------------------------------------------------
# Bundle and summarize data set
str( win.data8 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
   lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
   date2 = date2, dur2 = dur2) )

# Specify model in BUGS language
sink("Nmix.special.txt")
cat("
model {

  # Specify priors
  # abundance
  beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
  for(k in 1:7){            # Regression params in lambda
    beta[k] ~ dnorm(0, 1)
  }
  # Model for unexplained variance in lambda among sites
  for (i in 1:nsite){
    tau.lam[i] <- 1/var.lam[i]
    log(var.lam[i]) <- alpha.var.lam + beta.var.lam * elev[i]
  }
  # Priors for intercept and slope of linear model for variance
  alpha.var.lam ~ dunif(-1, 1)
  beta.var.lam ~ dunif(0, 3)

  # detection
  for(j in 1:3){
    alpha0[j] <- logit(mean.p[j])
    mean.p[j] ~ dunif(0, 1)# p intercept for occasions 1-3
  }
  for(k in 1:13){           # Regression params in p
    alpha[k] ~ dnorm(0, 1)
  }
  tau.p.survey <- pow(sd.p.survey, -2)
  sd.p.survey ~ dunif(0, 1) # site-survey heterogeneity in p

  # Poisson-lognormal model for abundance
  for (i in 1:nsite){
    eps.lam[i] ~ dnorm(0, tau.lam[i]) # Random site effects in log(abundance)
    loglam[i] <- beta0 + inprod(beta[], lamDM[i,]) + eps.lam[i]
    loglam.lim[i] <- min(250, max(-250, loglam[i]))  # 'Stabilize' log
    mu.poisson[i] <- exp(loglam.lim[i])
    N[i] ~ dpois(mu.poisson[i])
  }

  # Binomial measurement error model with extra-binomial dispersion
  for (i in 1:nsite){
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # 'Stabilize' logit
      lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
        alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
        alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
        alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
        alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
        alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
        alpha[13] * date[i,j] * dur2[i,j] + eps.p.survey[i,j]
        eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
  }
}
",fill = TRUE)
sink()


# Initial values
Nst <- apply(y, 1, max, na.rm = T) + 1
Nst[is.na(Nst)] <- round(mean(y, na.rm = TRUE))
Nst[Nst == "-Inf"] <- round(mean(y, na.rm = TRUE))
inits <- function(){ list(N = Nst, beta0 = 0, mean.p = rep(0.5,3),
    beta = runif(7, 0,0), alpha = runif(13, 0,0), alpha.var.lam = 0,
    beta.var.lam = 1.5, sd.p.survey = 0.3)}

# Parameters monitored
params <- c("beta0", "beta", "alpha.var.lam", "beta.var.lam", "alpha0",
    "mean.p", "alpha", "sd.p.survey")

# MCMC settings
# ni <- 180000    ;    nt <- 100    ;    nb <- 10000    ;    nc <- 3
ni <- 18000    ;    nt <- 10    ;    nb <- 1000    ;    nc <- 3  # ~~~~ for testing

# Call WinBUGS from R (ART 374 min) and summarize posteriors
out8 <- bugs(win.data8, inits, params, "Nmix.special.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir)  # ~~~ for autotesting
# ~~~ you should probably save the result ~~~~~~~~~~~~~~~
save(out8, file="AHM1_06.11_out8.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(out8, dig = 3)


print(cbind(out6$summary[c(4:11,13:31,33),1:2], out8$summary[c(1:8,11:30),1:2]),2)


# Predict site variance as a function of elevation
# Get posterior distribution of predictions first
orig.elev.pred <- seq(200, 2250, 50)
elev.pred <- (orig.elev.pred - mean(tits$elev)) / sd(tits$elev)
(n.mcmc <- length(out8$sims.list$alpha.var.lam))  # how many MCMC samples ?
post.sd.lam <- array(NA, dim = c(length(elev.pred), n.mcmc))
for(i in 1:length(elev.pred)){
  post.sd.lam[i,] <- sqrt(exp(out8$sims.list$alpha.var.lam +
      out8$sims.list$beta.var.lam * elev.pred[i]))
}

# Plot posterior mean and a sample of 500 regression lines from posterior
show <- sample(1:n.mcmc, 500)
matplot(orig.elev.pred, post.sd.lam[,show], xlab = "Elevation (m)",
    ylab = " sd.lam", type = "l", lty = 1, lwd = 1, col = "grey",
    frame = FALSE, ylim = c(0, 6))
lines(orig.elev.pred, apply(post.sd.lam, 1, mean), lwd = 3, col = "blue")

