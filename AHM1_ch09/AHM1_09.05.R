#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9. Advanced Hierarchical Distance Sampling
# =========================================================================

# Approximate execution time for this code: 1.9 hrs

library(AHMbook)
library(unmarked)

# 9.5 Open HDS models: temporary emigration
# =========================================

# 9.5.1 Data and model structure (no code)
# 9.5.2 Cautionary note on temporary emigration processes (no code)


# 9.5.3 Modeling temporary emigration with distance sampling in unmarked using the function gdistsamp
# ------------------------------------------------------------------------
# Load the wagtail data, investigate NA patterns in detection data
data("wagtail")
str(wagtail)
Y <- wagtail$Y
table(n.missing <- rowSums(is.na(Y))) # Frequency distribution of number of NAs per site
n.missing
keep <- which(n.missing == 0)     # Sites with complete distance data
Y <- Y[keep,]                     # restrict analysis to those

# Harvest other data for sites with complete distance data
potato <- wagtail$potato[keep]   ;   grass <- wagtail$grass[keep]
lscale <- wagtail$lscale[keep]   ;   hour <- wagtail$hour[keep,]
date <- wagtail$date[keep,]   ;   rep <- wagtail$rep[keep,]
breaks <- wagtail$breaks

# Look at the distance data
str(Y)
tmp <- apply(Y, 2, sum, na.rm = TRUE)
matplot(1:6, t(matrix(tmp, nrow = 4, byrow= TRUE)), type = "b", ylim = c(0, 90),
    xlab = "Distance class", ylab = "Number of wagtails", frame = FALSE,
    lwd = 3, lty = 1)

# Standardize all continuous covariates
mn.potato <- mean(potato)   ;   sd.potato <- sd(potato)
mn.grass <- mean(grass)   ;   sd.grass <- sd(grass)
mn.lscale <- mean(lscale)   ;   sd.lscale <- sd(lscale)
mn.date <- mean(date)    ;   sd.date <- sd(c(date))
mn.hour <- mean(hour)   ;   sd.hour <- sd(c(hour))
POTATO <- (potato - mn.potato) / sd.potato
GRASS <- (grass - mn.grass) / sd.grass
LSCALE <- (lscale - mn.lscale) / sd.lscale
DATE <- (date - mn.date) / sd.date
HOUR <- (hour - mn.hour) / sd.hour

# Package into unmarked GDS data frame and inspect the data
umf <- unmarkedFrameGDS(y = Y[,1:24], survey="point", unitsIn="m",
   dist.breaks=breaks, numPrimary = 4,
   siteCovs = data.frame(POTATO, GRASS, LSCALE),
   yearlySiteCovs=list(rep = rep, DATE = DATE, HOUR = HOUR))
str(umf)
summary(umf)


# Model fitting: Null models fm0
# exponential detection function
summary(fm0.exp <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
   keyfun = "exp", output = "density", unitsOut = "ha",
   mixture = "P", K = 100, se = TRUE, data = umf) )

# hazard detection function
summary(fm0.haz <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
   keyfun = "haz", output = "density", unitsOut = "ha",
   mixture = "P", K = 100, se = TRUE, data = umf ) )


# half-normal detection function
summary(fm0.hn <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
   keyfun = "halfnorm", output = "density", unitsOut = "ha",
   mixture = "P", K = 100, se = TRUE, data = umf,control=list(trace=TRUE, REPORT=1)) )

# Compare AIC scores for 3 detection functions
rbind('AIC exp' = fm0.exp@AIC, 'AIC haz' = fm0.haz@AIC, 'AIC hn' = fm0.hn@AIC)

backTransform(fm0.haz, type="lambda")

backTransform(fm0.haz, type="det")

# plot(1:300, gxhaz(1:300, shape = exp(5.13), scale=1.59), frame = FALSE, type = "l", xlab = "Distance Wagtail-Observer (metres)", ylab = "Detection probability", lwd=3)
plot(1:300, gxhaz(1:300, shape = exp(5.13), scale=exp(1.59)), frame = FALSE,
    type = "l", xlab = "Distance Wagtail-Observer (metres)",
    ylab = "Detection probability", lwd=3)
# See errata, 8 Nov 2018

# Model with time-dependent phi
fm1 <- gdistsamp(lambdaformula = ~1, phiformula = ~rep-1,
   pformula = ~1, keyfun = "haz", output = "density", unitsOut = "ha",
   mixture = "P", K = 100, se = TRUE, data = umf)

# Compare AIC for models with phi constant and phi time-dependent
rbind('AIC phi constant' = fm0.haz@AIC, 'AIC phi time-dep' = fm1@AIC)


# Add covariates on lambda:  2-phase fitting to assist convergence using
#   first K = 20
summary( fm2.init <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
   phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
   unitsOut = "ha",   control=list(trace=TRUE, REPORT=1),
   mixture = "P", K = 20, se = TRUE, data = umf))

starts <- coef(fm2.init)
summary( fm2  <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
   phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
   unitsOut = "ha",  starts=starts, control=list(trace=TRUE,
   REPORT=1), mixture = "P", K = 100, se = TRUE, data = umf))


# Add covariates on lambda: optimisation with default, no inits
summary( fm2a <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
   phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
   unitsOut = "ha", control=list(trace=TRUE, REPORT=1),
   mixture = "P", K = 100, se = TRUE, data = umf))


# Estimates of availability probabilities
plogis(fm2@estimates@estimates$phi@estimates)

# Models with time-dependent phi, AIC-best key functions and
#    three covariates on phi. Use previous estimates as starting values.
(tmp <- coef(fm2))
starts <- c(tmp[1:4], tmp[5:8], 0,0,0, tmp[9], tmp[10])

summary(fm3 <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
   phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE, pformula = ~1,
   keyfun = "haz", output = "density", unitsOut = "ha",
   mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
   se = TRUE, data = umf, starts = starts))

# Models with time-dependent phi, AIC-best key function, 3 covariates on phi
#     and, in addition, date and hour on detection
# linear effects on detection
tmp <- fm3@estimates@estimates
starts <- c(tmp$lambda@estimates, tmp$phi@estimates, tmp$det@estimates,
    0, 0, tmp$scale@estimates)
summary(fm4A <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
   phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE, pformula = ~ DATE + HOUR,
   keyfun = "haz", output = "density", unitsOut = "ha",
   mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
   se = TRUE, data = umf, starts = starts) )

# quadratic effects on detection
tmp <- fm4A@estimates@estimates
p.start <- tmp$det@estimates
p.start <- c(p.start[1:2], 0, p.start[3], 0)
starts <- c(tmp$lambda@estimates, tmp$phi@estimates, p.start, tmp$scale@estimates)

summary(fm4B <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
   phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE,
   pformula = ~ DATE + I(DATE^2) + HOUR + I(HOUR^2),
   keyfun = "haz", output = "density", unitsOut = "ha",
   mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
   se = TRUE, data = umf, starts = starts) )

starts <- coef(fm4B)[-16]   # Drop coef for HOUR^2
summary(fm4C <- gdistsamp(~ POTATO+GRASS+LSCALE,
~(rep-1)+ POTATO+GRASS+LSCALE, ~ DATE + I(DATE^2) + HOUR,
   keyfun = "haz", output = "density", unitsOut = "ha",
   mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
   se = TRUE, data = umf, starts = starts) )

starts <- c(coef(fm4C), 0)
summary(fm5 <- gdistsamp(~ POTATO+GRASS+LSCALE,
   ~(rep-1)+ POTATO+GRASS+LSCALE, ~ DATE + I(DATE^2) + HOUR,
   keyfun = "haz", output = "density", unitsOut = "ha",
   mixture = "NB", K = 100, control=list(trace=TRUE, REPORT=1),
   se = TRUE, data = umf , starts = starts) )

# Now we create a model selection table of these various models
modSel(fitList(fm0.haz, fm1, fm2, fm3, fm4A, fm4B, fm4C, fm5) )

summary(fm5)

# Bootstrap Goodness-of-fit assessment: ART ~ 20 hours
set.seed(1234)
# (pb <- parboot(fm5, fitstats, nsim=100, report=5))
(pb <- parboot(fm5, fitstats, nsim=10, report=5))  # ~~~~ reduce the number for testing

# Compute magnitude of "overdispersion" c.hat as ratio of observed to expected
#    chisquare test statistic
(c.hat <- pb@t0[2] / mean(pb@t.star[,2]))  # c-hat as ratio of observed/expected


# Predictions of lambda for POTATO, GRASS and LSCALE
newdat1 <- data.frame(POTATO=0, GRASS=0, LSCALE = seq(-1.8,4.33,,100))
newdat2 <- data.frame(POTATO=seq(-0.75,3,,100), GRASS=0, LSCALE = 0)
newdat3 <- data.frame(POTATO=0, GRASS=seq(-0.4, 3.6,,100), LSCALE = 0)
pred1 <- predict(fm5, type="lambda", newdata=newdat1, append = TRUE)
pred2 <- predict(fm5, type="lambda", newdata=newdat2, append = TRUE)
pred3 <- predict(fm5, type="lambda", newdata=newdat3, append = TRUE)

# Predictions of phi for POTATO, GRASS and LSCALE and for rep = 1
newdat4 <- data.frame(rep = factor('1', levels = c('1','2','3','4')),
    POTATO=0, GRASS=0, LSCALE = seq(-1.8,4.33,,100))
newdat5 <- data.frame(rep = factor('1', levels = c('1','2','3','4')),
    POTATO=seq(-0.75,3,,100), GRASS=0, LSCALE = 0)
newdat6 <- data.frame(rep = factor('1', levels = c('1','2','3','4')),
    POTATO=0, GRASS=seq(-0.4, 3.6,,100), LSCALE = 0)
pred4 <- predict(fm5, type="phi", newdata=newdat4, append = TRUE)
pred5 <- predict(fm5, type="phi", newdata=newdat5, append = TRUE)
pred6 <- predict(fm5, type="phi", newdata=newdat6, append = TRUE)

# Predictions of detection function sigma for DATE and HOUR
newdat7 <- data.frame(DATE = seq(-1.51,1.69,,100), HOUR = 0)
newdat8 <- data.frame(DATE=0, HOUR = seq(-1.92,3.1,,100))
pred7 <- predict(fm5, type="det", newdata=newdat7, append = TRUE)
pred8 <- predict(fm5, type="det", newdata=newdat8, append = TRUE)

op <- par(mfrow = c(1,3), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
plot(newdat1$LSCALE, pred1[,1], xlab="Standardized covariate",
    ylab="Density (birds/ha)", lwd=3,type="l", frame = FALSE)
lines(newdat2$POTATO, pred2[,1], lwd=3, col="red")
lines(newdat3$GRASS, pred3[,1], lwd=3, col="blue")
legend(-1.6, 1.65, c("LSCALE", "POTATO", "GRASS"),
    col=c("black", "red", "blue"), lty=1, lwd=3, cex=1.2)

plot(newdat4$LSCALE, pred4[,1], xlab="Standardized covariate",
    ylab="Availability (phi)", lwd=3,type="l", frame = FALSE)
lines(newdat5$POTATO, pred5[,1], lwd=3, col="red")
lines(newdat6$GRASS, pred6[,1], lwd=3, col="blue")
legend(2, 0.65, c("LSCALE", "POTATO", "GRASS"),
    col=c("black", "red", "blue"), lty=1, lwd=3, cex=1.2)

plot(newdat7$DATE, pred7[,1], xlab="Standardized covariate",
    ylab="Detection function (sigma)", lwd=3,type="l", frame = FALSE,
    ylim = c(100, 200))
lines(newdat8$HOUR, pred8[,1], lwd=3, col="red")
legend(0.5, 140, c("DATE", "HOUR"), col=c("black", "red"),
    lty=1, lwd=3, cex=1.2)
par(op)

# 9.5.4 Fitting temporary emigration HDS models in BUGS
# ------------------------------------------------------------------------

# 9.5.4.1 Simulating a temporary emigration system
# ------------------------------------------------------------------------
simHDSopen(type="line", nsites = 100, mean.lam = 2, beta.lam = 0, mean.sig = 1,
    beta.sig = 0, B = 3, discard0=TRUE, nreps=2, phi=0.7, nyears=5,
    beta.trend = 0)

# Obtain a temporary emigration data set
set.seed(1234)
str(tmp <- simHDSopen("point", nreps=7, nyears=5, nsites=100)  )
attach(tmp)


apply(tmp$M.true,2,sum)

# Define distance class information
delta <- 0.5
nD <- B%/%delta                 # Number of distance classes
midpt <- seq(delta/2, B, delta) # Mid-point of distance intervals

# Create the 4-d array
y4d <- array(0,dim=c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- tmp$data[[yr]][[rep]]
    site <- data[,1]
    dclass <- data[,"d"]%/%delta + 1
    ndclass <- B%/%delta
    dclass <- factor(dclass, levels= 1:ndclass)
# ~~~~~ this cannot work ~~~~~~~~~~
    # y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass)
# ~~~~~ use this instead ~~~~~~~~~~~
    ttt <- table(site, dclass)
    siteID <- as.numeric(rownames(ttt))
    y4d[siteID,1:nD,rep,yr] <- ttt

  }
}

y3d <- y4d[,,,1]


# Bundle and summarize the data set
nobs <- apply(y3d, c(1,3), sum)  # Total detections per site and occasion
str( data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta,
    habitat=habitat, B=B, nobs = nobs) )


# Define model in BUGS
cat("
model {
  # Prior distributions
  beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
  mean.lam <- exp(beta0)
  beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
  phi ~ dunif(0,1)        # Probability of availability
  sigma ~ dunif(0.01,5)   # Distance function parameter

  # Detection probs for each distance interval and related things
  for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal
    f[b] <- (2*midpt[b]*delta)/(B*B)    # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
  }
  cellprobs[nD+1]<- 1-sum(cellprobs[1:nD])

  for (s in 1:nsites) {
    for (k in 1:K) {
      pdet[s,k] <- sum(cellprobs[1:nD])   # Distance class probabilities
      pmarg[s,k] <- pdet[s,k]*phi         # Marginal probability

      # Model part 4: distance class frequencies
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k])
      # Model part 3: total number of detections:
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])
      # nobs[s,k] ~ dbin(pdet[s,k], Navail[s,k]) # Alternative formulation
      # Model part 2: Availability. Not used in this model but simulated.
      Navail[s,k] ~ dbin(phi, M[s])
    }  # end k loop
    # Model part 1: Abundance model
    M[s] ~ dpois(lambda[s])
    log(lambda[s]) <- beta0 + beta1*habitat[s]
  }  # End s loop

  # Derived quantities
  Mtot <- sum(M[])
  for(k in 1:K){
    Ntot[k]<- sum(Navail[,k])
  }
} # End model
",file="model.txt")

# Assemble the initial values and parameters to save for JAGS
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max)  +2
inits <- function(){
  list(M=Mst, sigma = 1.0, phi=0.9, beta0=log(2), beta1=0.5)
}
params <- c("sigma", "phi", "beta0", "mean.lam", "beta1", "Mtot", "Ntot")

# MCMC settings
# ni <- 60000   ;   nb <- 10000   ;   nt <- 5   ;   nc <- 3
ni <- 6000   ;   nb <- 1000   ;   nt <- 1   ;   nc <- 3  # ~~~~~ reduce for testing

# Run WinBUGS or JAGS
# library("R2WinBUGS")
library("jagsUI")  # JAGS works but WinBUGS does not!
# bd <- "c:/WinBUGS14/"
# out1 <- bugs(data, inits, parameters, "model.txt", n.thin=nthin,
#       n.chains=nc, n.burnin=nb,n.iter=ni,debug=TRUE, bugs.dir = bd)
# We get this error: vector valued relation y3d must involve consecutive
# elements of variable

# Run JAGS: This fails quite often ('invalid parent node'), just keep trying
outTE1 <- jags(data, inits, params, "model.txt",
  n.thin=nt,n.chains=nc, n.burnin=nb,n.iter=ni, parallel = TRUE)
traceplot(outTE1)   ;    print(outTE1, 3)            # ART 4 min


# Put true values into a vector
truth <- c(tmp$parms[c(1:3,5)], Mtot = sum(tmp$M[,1]),
    Ntot = (apply(tmp$Na.real[,,1],2,sum)))

# Get posterior means and 2.5% and 97.5% percentiles (95% CRI)
post <- outTE1$summary[c("mean.lam", "beta1", "sigma", "phi", "Mtot",
    "Ntot[1]", "Ntot[2]", "Ntot[3]" ,"Ntot[4]", "Ntot[5]", "Ntot[6]",
    "Ntot[7]"), c(1,3,7)]

# Table compares truth with posterior mean and 95% CRI from JAGS
cbind(truth, posterior = round(post, 3))


# 9.5.4.2 Bayesian analysis of the Wagtail data
# ------------------------------------------------------------------------
y3d <- array(NA,dim=c(nrow(Y), 6, 4) )          # Create 3d array
y3d[,,1] <- Y[,1:6]  ;  y3d[,,2] <- Y[,7:12]    # Fill the array
y3d[,,3] <- Y[,13:18]  ;  y3d[,,4] <- Y[,19:24]

K <- 4                          # Number of primary occasions
nsites <- nrow(Y)               # Number of sites
nD <- 6                         # Number of distance classes
midpt <- seq(25,275,50)         # Class midpoint distance
delta <- 50                     # Class width
B <- 300                        # Maximum distance
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion

# Bundle and summarize data set
area <- pi*(300^2)/10000
str(data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta,
   B=B, nobs=nobs, area=area))

# Write out the BUGS model file
cat("
model {

  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)

  # Availability parameter
  phi ~ dunif(0,1)

  # Detection parameter
  sigma ~ dunif(0,500)

  # Multinomial cell probabilities
  for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma)  # Half-normal model
    f[b] <- (2*midpt[b]*delta)/(B*B) # Scaled radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
  }
  cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])

  for (s in 1:nsites) {
    for (k in 1:K) {
      # Conditional 4-part version of the model
      pdet[s,k] <- sum(cellprobs[1:nD])
      pmarg[s,k] <- pdet[s,k]*phi
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k]) # Part 4: distance
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: number of detected individuals
      Navail[s,k] ~ dbin(phi,M[s])        # Part 2: Number of available individuals
    }  # end k loop

  M[s] ~ dpois(lambda[s])    #  Part 1: Abundance model
  log(lambda[s]) <- beta0    #  Habitat variables would go here
  }  # end s loop

  # Derived quantities
  for(k in 1:K){
    Davail[k] <- phi*exp(beta0)/area
  }
  Mtotal <- sum(M[])
  Dtotal<- exp(beta0)/area
} # end model
",fill=TRUE,file="wagtail.txt")

# Inits
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max,na.rm=TRUE) + 2
inits <- function() list(M=Mst, sigma = 100.0)

# Parameters to save
params <- c("sigma", "phi", "beta0", "beta1", "Mtotal", "Davail", "Dtotal")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Run JAGS (ART 3 min)
library("jagsUI")
wag1 <- jags(data, inits, params, "wagtail.txt", n.thin=nt,
    n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(wag1)
summary(wag1)

exp(-2.06)


# Bundle and summmarize data set for BUGS
rep <- matrix(as.numeric(rep), ncol=4)
area <- pi*(300^2)/10000
str(data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt = midpt,
    delta=delta, B=B, nobs=nobs, POTATO=POTATO, GRASS=GRASS, LSCALE=LSCALE,
    rep=rep, DATE=DATE, HOUR=HOUR, area=area))

# Define model in BUGS
cat("
model {

  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)

  # Availability parameters
  phi0 ~ dunif(0,1)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4){
    gamma1[k] ~ dunif(0, 1) # Availability effects of surveys 1 - 4
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k]))
  }
  gamma2 ~ dnorm(0, 0.01)
  gamma3 ~ dnorm(0, 0.01)
  gamma4 ~ dnorm(0, 0.01)

  # Detection parameters
  sigma0 ~ dunif(0.1,500)   # Intercept
  alpha2 ~ dnorm(0, 0.01)   # effect of DATE (linear)
  alpha3 ~ dnorm(0, 0.01)   # effect of DATE (squared)
  alpha4 ~ dnorm(0, 0.01)   # effect of HOUR
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)

  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*POTATO[s] + gamma3*GRASS[s] + gamma4*LSCALE[s]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0) + alpha2*DATE[s,k] + alpha3*pow(HOUR[s,k],2) + alpha4*HOUR[s,k]
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])

      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop

     M[s] ~ dnegbin(prob[s], r)
     prob[s] <- r/(r+lambda[s])
     # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
     log(lambda[s]) <- beta0 + beta1*POTATO[s] + beta2*GRASS[s] + beta3*LSCALE[s]
  }  # end s loop

  # Derived quantities
  for(k in 1:K){
    Davail[k] <- mean(phi[,k])*exp(beta0)/area
  }
  Mtotal <- sum(M[])
  Dtotal <- exp(beta0)/area
} # End model
",fill=TRUE,file="wagtail2.txt")

# Inits
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max,na.rm=TRUE) + 2
inits <- function() list(M=Mst, sigma0 = 100, alpha2=0, alpha3=0, alpha4=0,
    gamma2=0, gamma3=0, gamma4=0, beta1=0,beta2=0,beta3=0, r = 1)

# Parameters to save
params <- c("r","sigma0", "beta0", "beta1", "beta2", "beta3", "Mtotal",
    "alpha2", "alpha3",  "alpha4", "theta", "Dtotal", "Davail", "phi0",
    "gamma1", "gamma2", "gamma3", "gamma4" ,"logit.gamma1")

# MCMC settings
ni <- 32000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 5

# Run JAGS (ART 79 min), check convergence and summarize posteriors
wag2 <- jags(data, inits, params, "wagtail2.txt", n.thin=nt,n.chains=nc,
    n.burnin=nb,n.iter=ni, parallel = TRUE)

# Compare posterior means to MLEs obtained from unmarked
mle <- coef(fm5)
mle[12] <- exp(mle[12]) # convert to distance units
mle[16] <- exp(mle[16]) # back-transform the hazard parameter
mle[17] <- exp(mle[17]) # back-transform the NB dispersion parameter
bayes <- wag2$summary[,1]
bayes <- c(bayes[3:6],bayes[c(25:28,22:24)],bayes[c(2,8:10)],bayes[11],bayes[1])
bayes[1] <- log(exp(bayes[1])/area)# Convert from N per site to log(density) per ha

round( cbind(mle,bayes), 3)


# 9.5.4.3 Robust Design: Replicates within and among years (no code)
# ------------------------------------------------------------------------

