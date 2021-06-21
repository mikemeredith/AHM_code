#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

# Approximate run time for this script: 3.5 hrs
# Run time with the full number of iterations: 40 hrs

library(AHMbook)
library(jagsUI)

# 1.8 Modeling population dynamics at two temporal scales
# =======================================================

# 1.8.2 Modeling the two-scale dynamics of British butterflies
# ------------------------------------------------------------

# Load Marbled White data and do some data summaries (fig not shown)
library(AHMbook)
data(UKmarbledWhite)
str(white <- UKmarbledWhite)
op <- par(mfrow = c(2,3))
tmp <- tapply(white$C, list(white$year, white$site), function(x) length(x))
plot(table(c(tmp)), main = "Counts per site and year")
plot(table(white$site), main = "Counts per site")
plot(table(white$north), main = "Counts per latitude class")
plot(table(white$date), main = "Counts per date")
plot(table(white$year), main = "Counts per year")
plot(table(white$C), main = "Frequency distribution of counts")
par(op)

# Variables for use later
year <- white$year - 1990

# Get latitude for each site (n [ 80) and scale it
latitude <- tapply(white$north, white$site, mean)
mean.lat <- mean(latitude)
sd.lat <- sd(latitude)
lat <- (latitude - mean.lat ) / sd.lat

# Get standardized date
jdate <- white$date - 150 # 1 = Day 151
date <- (jdate/ 10) - mean((jdate/ 10)) # Express in units of 10 days

# Data bundle
str(bdata <- list(C = white$C, pop = white$site, year = year,
    date = date, lat = lat, lat2 = lat^2, npop = max(white$site), nyear = max(year),
    nobs = length(year), pi = pi) )
# List of 10
# $ C    : int [1:9651] 15 59 53 4 4 16 25 63 74 28 ...
# $ pop  : num [1:9651] 1 1 1 1 1 1 1 1 1 1 ...
# $ year : num [1:9651] 1 1 1 1 1 2 2 2 2 2 ...
# $ date : num [1:9651] -0.546 0.654 1.254 2.954 3.454 ...
# $ lat  : num [1:80(1d)] -0.415 -0.646 1.004 0.961 1.112 ...
# $ lat2 : num [1:80(1d)] 0.172 0.418 1.009 0.924 1.236 ...
# $ npop : num 80
# $ nyear: num 25
# $ nobs : int 9651
# $ pi   : num 3.14

# Specify model in BUGS language
cat(file = "modelPH.txt","
model {

  ### Priors and linear models for parameters
  # Initial abundance
  for(i in 1:npop){
    expn1[i] <- exp(loglam[i])
    loglam[i] ~ dnorm(mu.llam[i], tau.llam)
    mu.llam[i] <- alpha.llam + beta.llam[1]*lat[i] + beta.llam[2]*lat2[i]
  }
  alpha.llam <- log(mean.lambda)
  mean.lambda ~ dunif(1, 300)
  beta.llam[1] ~ dnorm(0, 0.1)
  beta.llam[2] ~ dnorm(0, 0.1)
  tau.llam <- pow(sd.llam, -2)
  sd.llam ~ dunif(0.001, 5)

  # Interannual population growth rate
  for(t in 1:(nyear-1)){
    for(i in 1:npop){
      gamma[i, t] <- exp(loggam[i, t])
      loggam[i,t] <- alpha0.lgam + alpha.lgam.site[i] + alpha.lgam.time[t] +
          beta.lgam[1] * (t-13) + beta.lgam[2]*lat[i] + beta.lgam[3]*lat2[i] +
          beta.lgam[4]*lat[i] * (t-13) + beta.lgam[5]*lat2[i]*(t-13)
    }
  }
  alpha0.lgam <- log(mean.gamma)
  mean.gamma ~ dunif(0.001, 3)
  for(i in 1:npop){
    alpha.lgam.site[i] ~ dnorm(0, tau.lgam.site)
  }
  for(t in 1:(nyear-1)){
    alpha.lgam.time[t] ~ dnorm(0, tau.lgam.time)
  }
  for(v in 1:5){
    beta.lgam[v] ~ dnorm(0, 1)
  }
  tau.lgam.site <- pow(sd.lgam.site, -2)
  sd.lgam.site ~ dunif(0.001, 1)
  tau.lgam.time <- pow(sd.lgam.time, -2)
  sd.lgam.time ~ dunif(0.001, 1)

  # Mean (or peak) activity period
  for(t in 1:nyear){
    for(i in 1:npop){ # Peak activity period
      mu[i, t] <- alpha0.mu + alpha.mu.site[i] + alpha.mu.time[t] + beta.mu[1] *
      (t-13) + beta.mu[2]*lat[i] + beta.mu[3]*lat2[i] + beta.mu[4]*lat[i]*(t-13) +
      beta.mu[5]*lat2[i]*(t-13)
    }
  }
  alpha0.mu <- mean.mu
  mean.mu ~ dunif(-3, 3)
  for(i in 1:npop){
    alpha.mu.site[i] ~ dnorm(0, tau.mu.site)
  }
  for(t in 1:nyear){
    alpha.mu.time[t] ~ dnorm(0, tau.mu.time)
  }
  for(v in 1:5){
    beta.mu[v] ~ dnorm(0, 1)
  }
  tau.mu.site <- pow(sd.mu.site, -2)
  sd.mu.site ~ dunif(0.001, 1)
  tau.mu.time <- pow(sd.mu.time, -2)
  sd.mu.time ~ dunif(0.001, 1)

  # Half-length of activity period
  for(t in 1:nyear){
    for(i in 1:npop){ # Length of activity period
      sigma[i, t] <- exp(lsig[i,t]) # Width
      lsig[i,t] <- alpha0.lsig + alpha.lsig.site[i] + alpha.lsig.time[t] +
          beta.lsig[1] * (t-13) + beta.lsig[2]*lat[i] + beta.lsig[3]*lat2[i] +
          beta.lsig[4]*lat[i] * (t-13) + beta.lsig[5]*lat2[i]*(t-13)
    }
  }
  alpha0.lsig <- log(mean.sigma)
  mean.sigma ~ dunif(0.001, 3)
  for(i in 1:npop){
    alpha.lsig.site[i] ~ dnorm(0, tau.lsig.site)
  }
  for(t in 1:nyear){
    alpha.lsig.time[t] ~ dnorm(0, tau.lsig.time)
  }
  for(v in 1:5){
    beta.lsig[v] ~ dnorm(0, 1)
  }
  tau.lsig.site <- pow(sd.lsig.site, -2)
  sd.lsig.site ~ dunif(0.001, 1)
  tau.lsig.time <- pow(sd.lsig.time, -2)
  sd.lsig.time ~ dunif(0.001, 1)

  # 'Likelihood'
  # Model for between-year dynamics
  for(i in 1:npop){
    # Initial year
    n[i,1] ~ dpois(expn1[i])
    n1[i] <- n[i,1]
    # Autoregressive transitions from t to t+1
    for(t in 2:nyear){
      n[i,t] ~ dpois(gamma[i, t-1]*n[i,(t-1)])
    }
  }

  # Phenomenological within-season population model
  for(i in 1:nobs){
    C[i] ~ dpois(lambda[i])
    lambda[i] <- n[pop[i],year[i]]*(1 / (sigma[pop[i],year[i]]*sqrt(2*pi)) )*exp( -
        pow((date[i] - mu[pop[i], year[i]]),2) / (2*pow(sigma[pop[i],year[i]], 2)) )
  }
}
")

# Initial values
nst <- tapply(bdata$C, list(bdata$pop, bdata$year), max, na.rm = TRUE)
nst[is.na(nst)] <- round(mean(nst, na.rm = TRUE))
inits <- function() list(n = 2*nst)

# Parameters monitored
# Choose if want both hyperparams and latent variables
params <- c("mu.llam", "mean.lambda", "alpha.llam", "beta.llam",
    "sd.llam", "alpha0.lgam", "mean.gamma", "alpha.lgam", "beta.lgam",
    "sd.lgam.site", "sd.lgam.time", "mean.mu", "alpha0.mu", "beta.mu",
    "sd.mu.site", "sd.mu.time", "mean.sigma", "mu.lsig", "alpha0.lsig",
    "beta.lsig", "sd.lsig.site", "sd.lsig.time", "expn1", "n1", "n", "alpha.lgam.site", "alpha.lgam.time", "gamma", "alpha.mu.site",
    "alpha.mu.time", "mu", "alpha.lsig.site", "alpha.lsig.time", "sigma")

# MCMC settings
# na <- 10000 ; ni <- 150000 ; nt <- 50 ; nb <- 100000 ; nc <- 3
na <- 1000 ; ni <- 15000 ; nt <- 5 ; nb <- 10000 ; nc <- 3 # ~~~~ for testing

# Call JAGS (ART 33 hours)
out13 <- jags(bdata, inits, params, "modelPH.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
summary(out13) ; jags.View(out13)

# Convergence check
# par(mfrow = c(3,3), mar = c(3,3,3,2))  # ~~~ no longer needed
traceplot(out13) # All params

# Posterior summary only for parameters with Rhat > 1.1
print(out13$summary[which(out13$summary[,8] > 1.1), -c(4:6)], 3)

# Select some key parameters to summarize their posteriors
dim(out13$summary)
cbind(1:2000, out13$summary[1:2000, -c(4:6)])
sel.params <- c(81:112)
print(out13$summary[sel.params,-c(4:6)], 3)

# ~~~~~~~~~~ extra code for figure 1.18 ~~~~~~~~~~~~~~~~~~~
op <- par(mfrow = c(2, 2), mar = c(5,5,4,3), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

# Observed data: mean counts per year
Years <- 1991:2015
# 'Mean observed counts per year'
plot(Years, tapply(white$C, white$year, mean), cex = 2, pch = 16, type = 'b',
    main ='', xlab = 'Year', ylab = 'Mean count', ylim = c(16, 45), frame = FALSE)
text(1992, 0.95*45, '(a)', cex = 2)

# Mean estimated n per year
# Compute posterior of average n per site and plot average n per year
# 'Mean estimated n per year'
str(tmp <- apply(out13$sims.list$n, c(1,3), mean))
tmpm <- apply(tmp, 2, mean)
tmpCI <- apply(tmp, 2, function(x) quantile(x, c(0.025, 0.975)))
plot(1991:2015, tmpm, type = 'b', lty = 1, pch = 16, cex = 2, xlab = 'Year',
    ylab = 'Mean n', main ='', ylim = c(50, 220), frame = FALSE)
segments(1991:2015, tmpCI[1,], 1991:2015, tmpCI[2,])
text(1992, 0.95*220, '(b)', cex = 2)

# Estimated n per site and year
# 'Estimated n per year and site'
matplot(1991:2015, t(out13$mean$n), type = 'l', lty = 1, lwd = 2,
    xlab = 'Year', ylab = 'Estimated n', main ='', frame = FALSE)
text(1992, 0.95*1220, '(c)', cex = 2)

# Mean estimated n per site
# Compute mean and sd (over years) of estimate n at each site
# 'Mean estimated n and SD per site'
tmpm <- apply(out13$mean$n, 1, mean)
ooo <- order(tmpm)
tmpSD <- apply(out13$mean$n, 1, mean)[ooo]
plot(1:80, tmpm[ooo], pch = 16, cex = 1.5, lty = 1, lwd = 2,
    xlab = 'Site Number', ylab = 'Mean n', main ='', frame = FALSE,
    ylim = c(0, 1200))
segments(1:80,tmpm[ooo]-tmpSD, 1:80, tmpm[ooo]+tmpSD)
text(5, 0.95*1220, '(d)', cex = 2)
par(op)

# ~~~~~~~~~~~ code for figure 1.19 ~~~~~~~~~~~
op <- par(mfrow = c(3, 1), mar = c(5,5,4,3), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

# Estimated gamma per site and year
range(out13$mean$gamma)
matplot(1991:2014, t(out13$mean$gamma), type = 'l', lty = 1, lwd = 2,
    xlab = 'Year', ylab = 'gamma', main = '', frame = FALSE, ylim = c(0.5, 2.2))
abline(h = 1, lwd = 1)
text(1992.5, 0.95*2.1, '(a)', cex = 2)

# Estimated mu per site and year
range(back.mu <- (out13$mean$mu + mean((jdate / 10 )) ) * 10)
matplot(1991:2015, t(back.mu), type = 'l', lty = 1, lwd = 2, xlab = 'Year',
    ylab = 'mu', main = '', frame = FALSE, ylim = c(60, 110))
text(1992.5, 0.95*100, '(b)', cex = 2)

# Estimated sigma per site and year
range(back.sigma <- out13$mean$sigma * 10)
matplot(1991:2015, t(back.sigma), type = 'l', lty = 1, lwd = 2, xlab = 'Year',
    ylab = 'sigma', main = '', frame = FALSE, ylim = c(6, 22))
text(1992.5, 0.95*22, '(c)', cex = 2)
par(op)

# ~~~~~~~~~~~~ code for figure 1.20 ~~~~~~~~~~~~~~~~
op <- par(mfrow = c(2, 2), mar = c(5,5,4,3), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

# Annual mean of mu and sigma, backtransformed into (Julian) days
range(back.mean.mu <- ((out13$mean$alpha0.mu  + out13$mean$alpha.mu.time) +
    mean((jdate / 10 )) ) * 10)
range(back.mean.sigma <- (exp(out13$mean$alpha0.lsig  + out13$mean$alpha.lsig.time) * 10))

# Actual mu and sigma per site and year, against year
range(back.mu <- (out13$mean$mu + mean((jdate / 10 )) ) * 10)
range(back.sigma <- (out13$mean$sigma * 10))

# All at once
curve(dnorm(x, mean = back.mean.mu[1], sd = back.mean.sigma[1]), 30, 130,
    frame = FALSE, main = '', xlab = 'Julian Date', ylab = '',
    ylim = c(0, 0.06), col = 'blue', lwd = 3)
for(t in 2:25){
  curve(dnorm(x, mean = back.mean.mu[t], sd = back.mean.sigma[t]), col = t,
      lwd = 3, add = TRUE)
}
text(40, 0.95*0.06, '(a)', cex = 2)

# Now produce plots of the phenology for all sites and years
# Also add annual mean into the bundle
letters <- c('(b)', '(c)', '(d)')
sel.years <- seq(5, 25, 10)
for(k in 1:3){
  t <- sel.years[k]
  curve(dnorm(x, mean = back.mu[1,t], sd = back.sigma[1,t]), 30, 130,
      frame = FALSE, main = '', xlab = 'Julian Date', ylab = '', ylim = c(0, 0.08))
for(i in 2:80){
    curve(dnorm(x, mean = back.mu[i,t], sd = back.sigma[i,t]), add = TRUE)
  }
  curve(dnorm(x, mean = back.mean.mu[t], sd = back.mean.sigma[t]), col = 'blue',
      lwd = 3, add = TRUE)
text(40, 0.95*0.08, letters[k], cex = 2)
}
par(op)

# ~~~~~~~~~~~~~ code for figure 1.21 ~~~~~~~~~~~~~~~~
# Three-dimensional predictions
# Create covariate values for prediction and scale as in real analysis
# Latitude
lat.original <- seq(77000, 304000, length.out =100)
mean.lat <- mean(latitude)  # Remember mean and sd of lat in 80 sites
sd.lat <- sd(latitude)
lat <- (lat.original - mean.lat ) / sd.lat

# Time (years)
yr.original <- seq(1991, 2015, length.out = 100)
yr <- seq(1, 25, length.out = 100) - 13
tmp <- out13$mean

# Predictions for gamma along two covariate gradients
predmat1 <- array(NA, dim = c(100, 100))
dimnames(predmat1) <- list(yr.original, lat.original)
for(i in 1:100){                   # i is for year (time)
  for(j in 1:100){                 # j is for latitude (space)
    predmat1[i,j] <- exp(tmp$alpha0.lgam + tmp$beta.lgam[1] * yr[i] +
        tmp$beta.lgam[2] * lat[j] + tmp$beta.lgam[3] * lat[j]^2 +
        tmp$beta.lgam[4] * lat[j] * yr[i] + tmp$beta.lgam[5] * lat[j]^2 * yr[i])
   }
}

# Predictions for mu along two covariate gradients (backtransformed)
predmat2 <- array(NA, dim = c(100, 100))
dimnames(predmat2) <- list(yr.original, lat.original)
for(i in 1:100){                   # i is for year (time)
  for(j in 1:100){                 # j is for latitude (space)
    predmat2[i,j] <- tmp$alpha0.mu + tmp$beta.mu[1] * yr[i] +
        tmp$beta.mu[2] * lat[j] + tmp$beta.mu[3] * lat[j]^2 +
        tmp$beta.mu[4] * lat[j] * yr[i] + tmp$beta.mu[5] * lat[j]^2 * yr[i]
  }
}
predmat2 <- (predmat2 + mean((jdate / 10 )) ) * 10 # Backtransform

# Predictions for sigma along two covariate gradients (backtransformed)
predmat3 <- array(NA, dim = c(100, 100))
dimnames(predmat3) <- list(yr.original, lat.original)
for(i in 1:100){                   # i is for year (time)
  for(j in 1:100){                 # j is for latitude (space)
    predmat3[i,j] <- 10 * exp(tmp$alpha0.lsig + tmp$beta.lsig[1] * yr[i] +
        tmp$beta.lsig[2] * lat[j] + tmp$beta.lsig[3] * lat[j]^2 +
        tmp$beta.lsig[4] * lat[j] * yr[i] + tmp$beta.lsig[5] * lat[j]^2 * yr[i])
   }
}

# All three sets of predictions in a single graph
op <- par(mfrow = c(1, 3), mar = c(5,5,5,2), cex.lab = 2, cex.axis = 2, cex.main = 2)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

# gamma
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x = lat.original/1000, y = yr.original, z = predmat1,
    col = mapPalette(100), axes = FALSE, xlab = "Latitude", ylab = "Year",
    main = "gamma")
contour(x = lat.original/1000, y = yr.original, z = predmat1, add = TRUE,
    col = "blue", labcex = 1.2, lwd = 1.5)
axis(1, at = seq(min(lat.original/1000), max(lat.original/1000), by = 100))
axis(2, at = seq(1990, 2015, by = 5))
box()

# mu
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x = lat.original/1000, y = yr.original, z = predmat2,
    col = mapPalette(100), axes = FALSE, xlab = "Latitude",
    ylab = "", main = "mu")
contour(x = lat.original/1000, y = yr.original, z = predmat2, add = TRUE,
    col = "blue", labcex = 1.2, lwd = 1.5)
axis(1, at = seq(min(lat.original/1000), max(lat.original/1000), by = 100))
axis(2, at = seq(1990, 2015, by = 5))
box()

# sigma
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x = lat.original/1000, y = yr.original, z = predmat3,
    col = mapPalette(100), axes = FALSE, xlab = "Latitude",
    ylab = "", main = "sigma")
contour(x = lat.original/1000, y = yr.original, z = predmat3, add = TRUE,
    col = "blue", labcex = 1.2, lwd = 1.5)
axis(1, at = seq(min(lat.original/1000), max(lat.original/1000), by = 100))
axis(2, at = seq(1990, 2015, by = 5))
box()
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
