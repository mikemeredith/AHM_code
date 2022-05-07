#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7. Modeling abundance using multinomial N-mixture models
# =========================================================================

# Approximate execution time for this code: 40 mins

library(unmarked)
library(jagsUI)

# 7.9 Example 3: Jays in the Swiss MHB
# ====================================


data(jay)        # Load data
str(jay)         # Inspect data list
dim(jay$caphist) # Look at detection history data

jay$caphist[1:4,]

# 7.9.1 Setting up the data and preparing for analysis
# ------------------------------------------------------------------------
crPiFun <- function(p) {
   p1 <- p[,1] # Extract the columns of the p matrix, one for
   p2 <- p[,2] #   each of J = 3 sample occasions
   p3 <- p[,3]
   cbind(      # define multinomial cell probabilities:
      "100" = p1 * (1-p2) * (1-p3),
      "010" = (1-p1) * p2 * (1-p3),
      "001" = (1-p1) * (1-p2) * p3,
      "110" = p1 * p2 * (1-p3),
      "101" = p1 * (1-p2) * p3,
      "011" = (1-p1) * p2 * p3,
      "111" = p1 * p2 * p3,
      "10x" = p1*(1-p2),
      "01x" = (1-p1)*p2,
      "11x" = p1*p2)
}


o2y <- matrix(1, 3, 10)

# Grab the data objects
covinfo <- jay$covinfo
gridinfo <- jay$gridinfo
sitecovs <- jay$sitecovs
caphist <- jay$caphist

# Get observation covariates to use in model
# Day of year, sample intensity and survey duration.
# Standardize them.
day <- scale(covinfo[,c("date1", "date2", "date3")])
dur <- as.matrix(covinfo[,c("dur1", "dur2", "dur3")])
dur[is.na(dur)] <- mean(dur, na.rm=TRUE) # Pad the 6 missing values
intensity <- dur / sitecovs[,"length"] # Sample rate = duration/length
dur <- scale(dur)
intensity <- scale(intensity)


reps <- apply(!is.na(day), 1, sum)
day[reps==2,3] <- 0
dur[reps==2,3] <- 0
intensity[reps==2, 3] <- 0

# Store the observation covariates in a list
obscovs <- list(intensity = intensity, dur = dur, day = day)

# Standardize site covariates
sitecovs[,"elev"] <- scale(sitecovs[,"elev"])
sitecovs[,"forest"] <- scale(sitecovs[,"forest"])
# NOTE: length is NOT standardized, but over-written with its reciprocal
sitecovs[,"iLength"] <- 1 / sitecovs[,"length"]

# Create unmarkedFrame (need crPiFun above and unmarked loaded)
caphist <- as.matrix(caphist)
mhb.umf <- unmarkedFrameMPois(y=caphist, siteCovs=as.data.frame(sitecovs),
   obsCovs=obscovs, obsToY=o2y, piFun="crPiFun")



# 7.9.2. Fitting some models
# ------------------------------------------------------------------------
# Fit a series of models
fm1 <- multinomPois(~1 ~1, mhb.umf)
fm2 <- multinomPois(~day ~1, mhb.umf)
fm3 <- multinomPois(~day + I(day^2) ~1, mhb.umf)
fm4 <- multinomPois(~intensity ~1, mhb.umf)
fm5 <- multinomPois(~intensity + I(intensity^2) ~1, mhb.umf)
fm6 <- multinomPois(~day + intensity ~1, mhb.umf)
fm7 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~1, mhb.umf)

# Assemble the models into a fitList and rank them by AIC
mspart1 <- fitList(
 "lam(.)p(.)" = fm1,
 "lam(.)p(day)" = fm2,
 "lam(.)p(day+day^2)" = fm3,
 "lam(.)p(intensity)" = fm4,
 "lam(.)p(intensity + intensity^2)" = fm5,
 "lam(.)p(day + rate)" = fm6,
 "lam(.)p(data + day^2 + intensity + intensity^2)" = fm7)

(mspart1 <- modSel(mspart1))

# Fit series of models
fm7 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~1, mhb.umf)
fm8 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~elev, mhb.umf)
fm9 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest, mhb.umf)
fm10 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~iLength, mhb.umf)
fm11 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + elev, mhb.umf)
fm12 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + iLength, mhb.umf)
fm13 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~elev + iLength, mhb.umf)
fm14 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + elev + iLength, mhb.umf)
fm15 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~elev + I(elev^2), mhb.umf)
fm16 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + elev + I(elev^2), mhb.umf)
fm17 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + elev + I(elev^2) + iLength, mhb.umf)


# Assemble the models into a fitList
mspart2 <- fitList(
    "lam(.)p(best)"                               = fm7,
    "lam(elev)p(best)"                            = fm8,
    "lam(forest)p(best)"                          = fm9,
    "lam(length)p(best)"                          = fm10,
    "lam(forest + elev)p(best)"                   = fm11,
    "lam(forest + iLength)p(best)"                = fm12,
    "lam(elev + iLength)p(best)"                  = fm13,
    "lam(forest + elev + length)p(best)"          = fm14,
    "lam(elev + elev^2)p(best)"                   = fm15,
    "lam(forest + elev + elev^2)p(best)"          = fm16,
    "lam(forest + elev + elev^2 + iLength)p(best)"= fm17)

# Rank them by AIC
(mspart2 <- modSel(mspart2))

fm17

mhb.umf2 <- unmarkedFrameGMM(y=caphist, numPrimary = 1,
    siteCovs=as.data.frame(sitecovs), obsCovs=obscovs, obsToY=o2y, piFun="crPiFun")

fm1NB <- gmultmix(~1, ~1, ~1, mix = "NB", data = mhb.umf2)

fm17P <- gmultmix(~forest + elev + I(elev^2) + iLength, ~1, ~day + I(day^2) +
    intensity + I(intensity^2), mix = "P",   data = mhb.umf2)

fm17NB <- gmultmix(~forest + elev + I(elev^2) + iLength, ~1, ~day + I(day^2) +
    intensity + I(intensity^2), mix = "NB",   data = mhb.umf2)


fm17P    # AIC-best Poisson mixture model


fm17NB  # NegBin version of AIC-best Poisson mixture model



# 7.9.3 Analysis of model fit
# ------------------------------------------------------------------------
# Define new fitstats function
fitstats2 <- function(fm) {
   observed <- getY(fm@data)
   expected <- fitted(fm)
   resids <- residuals(fm)
   n.obs<- apply(observed,1,sum,na.rm=TRUE)
   n.pred<- apply(expected,1,sum,na.rm=TRUE)
   sse <- sum(resids^2,na.rm=TRUE)
   chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
   freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
   freeTuke.n<- sum((sqrt(n.obs)-sqrt(n.pred))^2,na.rm=TRUE)
   sse.n <- sum( (n.obs -n.pred)^2,na.rm=TRUE)
   chisq.n <- sum((n.obs - n.pred)^2 / expected,na.rm=TRUE)

   out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke,
      SSE.n = sse.n, Chisq.n = chisq.n, freemanTukey.n=freeTuke.n)
   return(out)
}


(pb.mhb <- parboot(fm17, fitstats2, nsim=1000, report=1))


(pb.mhbNB <- parboot(fm17NB, fitstats2, nsim=1000, report=1))

# Compute c-hat
pb.mhbNB@t0[2] / mean(pb.mhbNB@t.star[,2])

n.obs <- apply(caphist, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(fm17P), 1, sum, na.rm=TRUE)
n.predNB <- apply(fitted(fm17NB), 1, sum, na.rm=TRUE)
plot(n.obs, n.predP, frame = FALSE)   # Fig. 7-7
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df = 4), type = "l", lwd = 2, col = "blue")
points(smooth.spline(n.predNB ~ n.obs, df=4), type= "l", lwd = 2, col = "red")



# 7.9.4 Summary analyses of jay models
# ------------------------------------------------------------------------
range(siteCovs(mhb.umf)$elev)

elev.mean <- attr(siteCovs(mhb.umf)$elev, "scaled:center")
elev.sd <- attr(siteCovs(mhb.umf)$elev, "scaled:scale")
elev.orig <- elev.sd*seq(-1.5, 2.42,,500)  + elev.mean


# Remember length = 0 is saturation sampling because length = 1/L
newL <- data.frame(elev = seq(-1.5,2.42,,500), elev.orig,
    forest = 0, iLength = 1/5.1)           # 'Low' prediction
newH <- data.frame(elev = seq(-1.5,2.42,,500), elev.orig,
    forest = 0, iLength = 0)           # 'High' prediction
predL <- predict(fm17NB, type="lambda", newdata=newL, appendData=TRUE)
predH <- predict(fm17NB, type="lambda", newdata=newH, appendData=TRUE)
head(cbind(low = predL[,1:2], high = predH[,1:2]))

plot(Predicted ~ elev.orig, predL, type="l", lwd = 3, xlab="Elevation",
 ylab="Expected # territories", ylim=c(0, 13), frame=FALSE, col = "red")
points(Predicted ~ elev.orig, predH, type="l", lwd = 3, col = "blue")
matlines(elev.orig, predL[,3:4], lty = 1, lwd = 1, col = "red")
matlines(elev.orig, predH[,3:4], lty = 1, lwd = 1, col = "blue")

b <- coef(fm17NB)[3]
c <- coef(fm17NB)[4]
elev.opt <- -b / (2*c)

(elev.opt <- elev.opt*elev.sd + elev.mean)

require(AICcmodavg)
model.list <- list(fm17NB) # candidate model list with single model
model.names <- c("AIC-best model")

# Compute model-averaged predictions of abundance for values of elevation,
#  with uncertainty (SE, CIs) adjusted for overdispersion (c.hat), with latter
#  estimated from bootstrapped Chisquare
pred.c.hatL <- modavgPred(cand.set = model.list, modnames = model.names,
    newdata = newL, parm.type = "lambda", type = "response", c.hat = 1.11)

# Compare predictions and SE without and with c.hat adjustment
head(cbind(predL[1:2], pred.c.hatL[2:3]), 10) ## see errata


rlength <- seq(1, 30, 0.01)         # Vary route length from 1 to 30 kms
newData <- data.frame(elev=0, forest=0, iLength=1/rlength)
pred <- predict(fm17NB, type="lambda", newdata=newData, c.hat = 1.11)
op <- par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.3)
plot(rlength, pred[[1]], type = "l", lwd = 3, col = "blue", frame = FALSE,
    xlab = "Transect length (km)", ylab = "Exposed population (lambda)",
    ylim = c(0, 16), axes = FALSE)
axis(1, at = seq(2,30,2))       ;      axis(2)
abline(v = c(1.2, 5.135, 9.4), lwd = 2)
matlines(rlength, cbind(pred[[1]]-pred[[2]], pred[[1]]+pred[[2]]),
    type = "l", lty = 1, lwd = 2, col = "gray")
sat.pred <- predict(fm17NB, type="lambda",
    newdata= data.frame(elev=0, forest=0, iLength=0), c.hat = 1.11)
abline(as.numeric(sat.pred[1]),0, lwd = 2, lty = 2)
par(op)


pred[round(rlength,2)==1.2,]/as.numeric(sat.pred[1])

pred[round(rlength,2)==5.14,]/as.numeric(sat.pred[1])

pred[round(rlength,2)==9.4,]/as.numeric(sat.pred[1])


# 7.9.5 Spatial prediction
# ------------------------------------------------------------------------
library(raster)
# library(rgdal)  # ~~~~ not necessary ~~~~

# Swiss landscape data and shape files
data(Switzerland)         # Load Swiss landscape data from unmarked
CH <- Switzerland         # this is for 'confoederatio helvetica'
head(CH)
gelev <- CH[,"elevation"] # Median elevation of quadrat
gforest <- CH[,"forest"]
grid <- CH[,c("x", "y")]
# ~~~~~ these shape files are not available ~~~~~~~~~~~~
# lakes <- readOGR(".", "lakes")
# rivers <- readOGR(".", "rivers")
# border <- readOGR(".", "border")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Draw two maps of Swiss elevation and forest cover (Fig. 7-10)
op <- par(mfrow = c(1,2), mar = c(1,2,3,5))
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))
mapPalette2 <- colorRampPalette(c("grey", "lightgreen", "darkgreen"))
r1 <- rasterFromXYZ(cbind(x = CH$x, y = CH$y, z = CH$elevation))
r2 <- rasterFromXYZ(cbind(x = CH$x, y = CH$y, z = CH$forest))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Elevation (m a.s.l.)", zlim = c(0, 4000))
# plot(rivers, col = "dodgerblue", add = TRUE)
# plot(border, col = "transparent", lwd = 1.5, add = TRUE)
# plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)

plot(r2, col = mapPalette2(100), axes = FALSE, box = FALSE, main = "Forest cover (%)", zlim = c(0, 100))
# plot(rivers, col = "dodgerblue", add = TRUE)
# plot(border, col = "transparent", lwd = 1.5, add = TRUE)
# plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
par(op)

# Standardize elevation for all grid cells using the mean at sample plots
elev.mean <- attr(siteCovs(mhb.umf)$elev, "scaled:center")
elev.sd <- attr(siteCovs(mhb.umf)$elev, "scaled:scale")
gelev <- (gelev - elev.mean) / elev.sd

# Standardize forest cover also using the mean at sample plots
forest.mean <- attr(siteCovs(mhb.umf)$forest, "scaled:center")
forest.sd <- attr(siteCovs(mhb.umf)$forest, "scaled:scale")
gforest <- (gforest - forest.mean) / forest.sd


# Form predictions for Swiss landscape (slooow!)
newL <- data.frame(elev=gelev, forest=gforest, iLength=1/5.1)
newH <- data.frame(elev=gelev, forest=gforest, iLength=0)
pred.mhb.NB.Low <- predict(fm17NB, type="lambda", newdata=newL, appendData=T)
pred.mhb.NB.High <- predict(fm17NB, type="lambda", newdata=newH, appendData=T)

# Create rasters and mask for elevation (mask areas > 2250 m)
r1 <- rasterFromXYZ(cbind(x = CH$x, y = CH$y, z = pred.mhb.NB.Low[,1]))
r2 <- rasterFromXYZ(cbind(x = CH$x, y = CH$y, z = pred.mhb.NB.Low[,2]))
elev <- rasterFromXYZ(cbind(CH$x, CH$y, gelev))
elev[elev > 2250] <- NA
r1 <- mask(r1, elev)
r2 <- mask(r2, elev)

# Draw maps of jay density and standard error of density (Fig. 7–11)
op <- par(mfrow = c(1,2), mar = c(1,2,3,5))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE,
    main = "Density of European Jay", zlim = c(0, 10))
# plot(rivers, col = "dodgerblue", add = TRUE)
# plot(border, col = "transparent", lwd = 1.5, add = TRUE)
# plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
plot(r2, col = mapPalette1(100), axes = FALSE, box = FALSE,
    main = "Standard errors of density", zlim = c(0, 1.5))
# plot(rivers, col = "dodgerblue", add = TRUE)
# plot(border, col = "transparent", lwd = 1.5, add = TRUE)
# plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
par(op)


# 7.9.6 Population size of jays in Switzerland
# ------------------------------------------------------------------------
out <- which(CH$water > 50 | CH$elevation > 2250)
(N.jay <- sum(pred.mhb.NB.Low[-out,1]))      # 'low prediction'

(N.jay <- sum(pred.mhb.NB.High[-out,1]))     # 'high prediction'


Nhat <- function(fm) {
   betavec <- coef(fm)[1:5]
   Xg <- cbind(rep(1, length(gforest)), gforest, gelev, gelev*gelev, 1/5.1)
   predLow <- as.numeric(exp(Xg%*%(betavec)))
   predHigh <- as.numeric(exp(Xg[,-5]%*%(betavec[-5])))
   out <- which(CH$water > 50 | CH$elevation > 2250)
   Nlow <- sum(predLow[-out])
   Nhigh <- sum(predHigh[-out])
   return(c(Nlow = Nlow, Nhigh = Nhigh))
}

set.seed(2015)
# ~~~~~~ parboot now runs in parallel by default, but Nhat doesn't work ~~~~~~~~
# (pb.N <- parboot(fm17NB, Nhat, nsim=100, report=1))
(pb.N <- parboot(fm17NB, Nhat, nsim=100, report=1, parallel=FALSE))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 7.9.7 Bayesian Analysis of the MHB Data
# ------------------------------------------------------------------------
# Extract data and do data augmentation up to M = 400
y <- as.matrix(getY(mhb.umf))

# Now  we have to stretch out the encounter frequencies into individuals...
# There were 439 unique individuals observed during the survey
eh <- unlist(dimnames(y)[2])
ehid <- col(y)[y>0 & !is.na(y)]  # Column ids, index to encounter history
eh <- eh[ehid]
siteid <- row(y)[y>0 & !is.na(y)] # Site ids
y <- y[y > 0 & !is.na(y)]   # Positive counts
eh <- rep(eh, y)
siteid <- rep(siteid, y)

eh.mat <- matrix(NA,nrow=length(eh),ncol=3)
for(i in 1:length(eh)){
  eh.mat[i,]<- as.numeric(unlist(strsplit(eh[i],split="")))
}


# Define some things and do the data augmentation
nsites = nrow(sitecovs)
nind <- nrow(eh.mat)
M <- 800
y <- rbind(eh.mat, matrix(0, nrow=(M-nind), ncol=3))

# Augment site ID
site <- c(siteid, rep(NA, M-nind))
sitecovs <- siteCovs(mhb.umf)

obscovs <- obsCovs(mhb.umf)
Intensity <- matrix(obscovs[,"intensity"], nrow=nsites, ncol=3,byrow=TRUE)

# Bundle data for BUGS
data <- list(y = y, J = 3, M = M , nsites=nsites, X = as.matrix(sitecovs),
    Intensity= Intensity, group=site)
str(data)

# Specify model in BUGS language
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0 / (1-p0))
  alpha1 ~ dnorm(0, 0.01)

  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  beta2 ~ dnorm(0,0.01)
  beta3 ~ dnorm(0,0.01)
  psi <- sum(lambda[]) / M   # psi is a derived parameter

  # Model for abundance: lambda depends on Elev, Length, Forest
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * X[s,1] + (beta2/X[s,2]) + beta3*X[s,3]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables

    # Observation model: p depends on Intensity
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * Intensity[group[i],j]

      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }
  }
}
",fill=TRUE,file="model.txt")

# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2",
    "beta3", "psi")

# Initial values
zst <- c(rep(1,M-100), rep(0,100)) ## see errata
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0=runif(1),
    beta1=rnorm(1), z= zst ) }

# MCMC settings
ni <- 11000   ;   nb <- 1000   ;   nt <- 4   ;   nc <- 3

# Call JAGS from R and summarize marginal posteriors

out <- jags(data, inits, parameters, "model.txt",
  # n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni)
  n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel=TRUE)  # ~~~ for testing

print(out, digits = 2)

fm17P    # AIC-best Poisson mixture model

