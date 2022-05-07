#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# Approximate execution time for this code: 2.25 hrs
# Run time with the full number of iterations: 31 hrs

library(AHMbook)
library(unmarked)

# 6.9 Abundance mapping of Swiss Great tits with unmarked
# =======================================================


# 6.9.1 Set up of the analysis
# ------------------------------------------------------------------------
## Code modified to use the SwissTits data set included in the AHMbook package
data(SwissTits)
?SwissTits
str(SwissTits)
SwissTits$species  # Available species

# Select Great tit and covariate data from 2013 and
#   drop 4 sites not surveyed in 2013
y0 <- SwissTits$counts[, , '2013', 'Great tit']
( NA.sites <- which(rowSums(is.na(y0)) == 3) ) # Unsurveyed sites
y <- y0[-NA.sites, ]                 # Drop them from the count data
tits <- SwissTits$sites[-NA.sites, ] # Also drop from the site covariates
str(y)  # Check the matrix of count data
# Get date and duration data for 2013, without the NA.sites rows:
date <- SwissTits$date[-NA.sites, , '2013']
dur <- SwissTits$dur[-NA.sites, , '2013']

# Plot observed data: counts vs survey date (Fig. 6-9)
matplot(t(date), t(y), type = "l", lwd = 3, lty = 1, frame = FALSE,
    xlab = "Survey data (1 = April 1)", ylab = "Count of Great Tits")

# Load unmarked, create unmarked data frame and inspect result
library(unmarked)
time <- matrix(rep(as.character(1:3), nrow(y)), ncol = 3, byrow = TRUE)
umf <- unmarkedFramePCount(y = y,
  siteCovs=data.frame(elev=scale(tits[,"elev"]), forest=scale(tits[,"forest"]),
      iLength=1/tits[,"rlength"]),
  obsCovs=list(time = time, date = scale(date), dur = scale(dur)))
summary(umf)                            # Summarize unmarked data frame
summary(apply(y, 1, max, na.rm = TRUE)) # Summarize max counts


# 6.9.2 Model fitting
# ------------------------------------------------------------------------
fm1 <- pcount(~ (elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
    ~ (elev+I(elev^2)) * (forest+I(forest^2)) + iLength,
    umf, control=list(trace=TRUE, REPORT=5))
summary(fm1)   ;   fm1@AIC


fm1.K500 <- pcount(fm1@formula, umf, control=list(trace=T, REPORT=5), K = 500)
summary(fm1.K500)   ;   fm1.K500@AIC


fm2 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
      - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
      - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
      - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
      ~ (elev+I(elev^2)) * (forest+I(forest^2))
      + iLength, starts = coef(fm1)[1:31],
      umf, control=list(trace=TRUE, REPORT=5))
summary(fm2)                      # AIC = 3695.792


fm3 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
      - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
      - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
      - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
      - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2)
      ~ (elev+I(elev^2)) * (forest+I(forest^2))
      + iLength, starts = coef(fm2)[-c(23, 27, 31)],
      umf, control=list(trace=TRUE, REPORT=5))
summary(fm3)                      # AIC = 3691.184


fm4 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
      - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
      - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
      - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
      - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2)
      - elev:I(date^2) - I(date^2):dur
      ~ (elev+I(elev^2)) * (forest+I(forest^2))
      + iLength, starts = coef(fm3)[-c(21, 28)],
      umf, control=list(trace=TRUE, REPORT=5))
summary(fm4)                      # AIC = 3687.565


fm5 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
      - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
      - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
      - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
      - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2)
      - elev:I(date^2) - I(date^2):dur
      ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength
      - I(elev^2):forest - I(elev^2):I(forest^2),
      starts = coef(fm4)[-c(9:10)],
      umf, control=list(trace=TRUE, REPORT=5))
summary(fm5)                      # AIC = 3686.094


# Negative binomial (NB) mixture
fm5NB <- pcount(fm5@formula, starts = c(coef(fm5),0),
      umf, control=list(trace=TRUE, REPORT=5), mixture = "NB")
summary(fm5NB)                      # AIC = 3181.046

# Zero-inflated Poisson (ZIP) mixture
fm5ZIP <- pcount(fm5@formula, starts = c(coef(fm5),0),
      umf, control=list(trace=TRUE, REPORT=5), mixture = "ZIP")
summary(fm5ZIP)                      # AIC = 3636.058


cbind(rbind("Poisson" = exp(coef(fm5)[1]), 
        "NegBin" = exp(coef(fm5NB)[1]), 
        "ZIP" = exp(coef(fm5ZIP)[1])), 
    rbind(plogis(coef(fm5)[15:17]), 
        plogis(coef(fm5NB)[15:17]), 
        plogis(coef(fm5ZIP)[15:17])))

# 6.9.3 Model criticism and goodness of fit
# ------------------------------------------------------------------------
library(AICcmodavg)
# ~~~ need to limit the cores used by Nmix.gof.test if other software is running
# ncores <- parallel::detectCores() - 1 # the default, ok if nothing else is running
ncores <- 3 # ~~~ for testing

system.time(gof.P <- Nmix.gof.test(fm5, nsim=100, ncores=ncores))  # 18 mins with 3 cores << 65 min
# ~~~~~~~~~~ changes in AICcmodavg::Nmix.gof.test wef 2.2-0 (25 February 2019) ~~~~~~~
# Nmix.gof.test now has a 'parallel' argument with default TRUE. The fm5NB and fm5ZIP models
#  using formula=fm5@formula will not work in parallel. With parallel=FALSE they work:

# system.time(gof.NB <- Nmix.gof.test(fm5NB, nsim=100, parallel=FALSE))   # 131 min
# system.time(gof.ZIP <- Nmix.gof.test(fm5ZIP, nsim=100, parallel=FALSE)) # 69 min

# Alternatively, rerun fm5NB and fm5ZIP with the formula specified in full:
fm5NB <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1 -
      elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur -
      elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2) -
      I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2) -
      I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2) -
      elev:I(date^2) - I(date^2):dur
      ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength -
      I(elev^2):forest - I(elev^2):I(forest^2), starts = c(coef(fm5),0),
      umf, control=list(trace=TRUE, REPORT=5), mixture = "NB")

fm5ZIP <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1 -
      elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur -
      elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2) -
      I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2) -
      I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2) -
      elev:I(date^2) - I(date^2):dur
      ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength -
      I(elev^2):forest - I(elev^2):I(forest^2), starts = c(coef(fm5),0),
      umf, control=list(trace=TRUE, REPORT=5), mixture = "ZIP")
system.time(gof.NB <- Nmix.gof.test(fm5NB, nsim=100, parallel=TRUE, ncores=ncores))   # 38 min << 131 min
system.time(gof.ZIP <- Nmix.gof.test(fm5ZIP, nsim=100, parallel=TRUE, ncores=ncores)) # 21 min << 69 min
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gof.P   ;   gof.NB   ;   gof.ZIP                        # print results


# Look at data, fitted values and residuals and produce plots
print(cbind(y, fitted(fm5), residuals(fm5)), 2)  # For Poisson model
plot_Nmix_resi(fm5, fm5NB, fm5ZIP)               # Produces Fig. 6–10 ## function renamed

plot_Nmix_resi <- function(fmP, fmNB, fmZIP){  # ~~~~~ this function is in AHMbook package
  # Function does diagnostic plots for one Nmix model fitted with all three
  #   mixture distributions currently available in unmarked:
  #   Poisson, negative binomial and zero-inflated Poisson
  # For each, fitted values vs. observed data and
  #   residuals vs. fitted values are plotted.
  library(unmarked)

  # Plot fitted vs. observed data
  op <- par(mfrow = c(2,3), mar = c(4,4,2,2), cex = 1.2)
  tmp1 <- range(c(fitted(fmP), fitted(fmNB), fitted(fmZIP)), na.rm = T)
  limits1 = round(c(tmp1[1], tmp1[2]))
  tmp2 <- range(c(residuals(fmP), residuals(fmNB), residuals(fmZIP)), na.rm = T)
  limits2 = round(c(tmp2[1], tmp2[2]))

  plot(fitted(fmP)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (P)",
      frame = FALSE, ylim = limits1)
  abline(0,1, lwd = 3 )
  abline(lm(c(fitted(fmP))~ c(fmP@data@y)), col = "blue", lwd = 3)
  plot(fitted(fmNB)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (NB)",
      frame = FALSE, ylim = limits1)
  abline(0,1, lwd = 3)
  abline(lm(c(fitted(fmNB))~ c(fmP@data@y)), col = "blue", lwd = 3)
  plot(fitted(fmZIP)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (ZIP)",
      frame = FALSE, ylim = limits1)
  abline(0,1, lwd = 3)
  abline(lm(c(fitted(fmZIP)) ~ c(fmP@data@y)), col = "blue", lwd = 3)

  # Plot residuals vs. fitted values
  plot(residuals(fmP)~ fitted(fmP), xlab = "Fitted values (P)", ylab = "Residuals",
      frame = FALSE, xlim = limits1, ylim = limits2)
  abline(h = 0, lwd = 2)
  abline(lm(c(residuals(fmP)) ~ c(fitted(fmP))), col = "blue", lwd = 3)
  plot(residuals(fmNB)~ fitted(fmNB), xlab = "Fitted values (NB)", ylab = "Residuals",
      frame = FALSE, xlim = limits1, ylim = limits2)
  abline(h = 0, lwd = 2)
  abline(lm(c(residuals(fmNB)) ~ c(fitted(fmNB))), col = "blue", lwd = 3)
  plot(residuals(fmZIP)~ fitted(fmZIP), xlab = "Fitted values (ZIP)", ylab = "Residuals",
      frame = FALSE, xlim = limits1, ylim = limits2)
  abline(h = 0, lwd = 2)
  abline(lm(c(residuals(fmZIP)) ~ c(fitted(fmZIP))), col = "blue", lwd = 3)
  par(op)
}


# Compute RMSE for all three models
(RMSEP <- sqrt(mean((y - fitted(fm5))^2, na.rm = TRUE)))      # Poisson
(RMSENB <- sqrt(mean((y - fitted(fm5NB))^2, na.rm = TRUE)))   # NB
(RMSEZIP <- sqrt(mean((y - fitted(fm5ZIP))^2, na.rm = TRUE))) # ZIP


map.Nmix.resi <- function(fm, x = tits$coordx, y = tits$coordy){
  # Function produces a map of the mean residuals from an N-mixture model
  #    object named fm, which was fit by function pcount in unmarked
  # Function arguments are the fitted model object and the x and y coordinates
  #    of every site
  library(sp)
  mean.resi <- apply(residuals(fm), 1, mean, na.rm = TRUE)
  mean.resi[mean.resi == "NaN"] <- mean(mean.resi, na.rm = TRUE)
  spdata <- data.frame(residuals = mean.resi, x = x, y = y)
  coordinates(spdata) <- c("x", "y")
  plot(bubble(spdata, "residuals", col = c("blue", "red"),
      main = paste("Average residuals of fitted N-mixture model")))
}

map.Nmix.resi(fm5, x = tits$coordx, y = tits$coordy)    # Map of average residuals for Poisson model ## function defaults changed in AHMbook 0.1.3
map.Nmix.resi(fm5NB, x = tits$coordx, y = tits$coordy)  # Map of average residuals for NB model
map.Nmix.resi(fm5ZIP, x = tits$coordx, y = tits$coordy) # Map of average residuals for ZIP model


binFittedP <- fitted(fm5) %/% 2.5                     # Bin fitted values
mMeanP <- tapply(fitted(fm5)^2, binFittedP, mean, na.rm = TRUE)   # Mean mean
mVarP <- tapply(residuals(fm5)^2, binFittedP, mean, na.rm = TRUE) # Mean variance
nsampleP <- table(binFittedP)                         # Sample size
binFittedNB <- fitted(fm5NB) %/% 2.5
mMeanNB <- tapply(fitted(fm5NB)^2, binFittedNB, mean, na.rm = TRUE)
mVarNB <- tapply(residuals(fm5NB)^2, binFittedNB, mean, na.rm = TRUE)
nsampleNB <- table(binFittedNB)
plot(mMeanP, mVarP, xlab = "Binned mean fitted response",
    ylab = "Mean variance of response", frame = FALSE, cex = log(nsampleP),
    col = "grey", pch = 16, cex.lab = 1.5)
points(mMeanNB, mVarNB, cex = log(nsampleNB), col = "green", pch = 16)


# 6.9.4 Analysis of results
# ------------------------------------------------------------------------
# Two new data sets for prediction: for lambda and for p (200 data points each)
lamNewData <- data.frame(elev = (seq(200, 2250,,200) - mean(tits[,"elev"]))/ sd(tits[,"elev"]),
    forest = 0, iLength = 1/5.1)
pNewData <- data.frame(elev = 0, time = factor("2", levels = c("1", "2", "3")),
    dur = 0, date = (seq(1,90,,200) - mean(date, na.rm = T))/sd(date, na.rm = T))

# Predictions for lambda and p, with SE and CIs, but no overdispersion
predict(fm5, type="state", newdata=lamNewData)    # Poisson model
predict(fm5, type="det", newdata=pNewData)
predict(fm5NB, type="state", newdata=lamNewData)  # NegBin model
predict(fm5NB, type="det", newdata=pNewData)
predict(fm5ZIP, type="state", newdata=lamNewData) # ZIP
predict(fm5ZIP, type="det", newdata=pNewData)


# Predictions for lambda only, incl. SE and CIs, with overdispersion
predictSE(fm5ZIP, newdata=lamNewData, print.matrix = TRUE, type="response",
    parm.type = "lambda", c.hat = 2.47)


# Predictions for lambda and p, incl. SE and with overdispersion, but no CIs
# Poisson model, for natural and for link scale of both lambda and p
modavgPred(cand.set = list(fm5), newdata=lamNewData, parm.type = "lambda",
    type = "response", c.hat = 3.82)
modavgPred(cand.set = list(fm5), newdata=lamNewData, parm.type = "lambda",
    type = "link", c.hat = 3.82)    # Could be used to get 95% CIs
modavgPred(cand.set = list(fm5), newdata=pNewData, parm.type = "detect",
    type = "response", c.hat = 3.82)
modavgPred(cand.set = list(fm5), newdata=pNewData, parm.type = "detect",
    type = "link", c.hat = 3.82)      # Could be used to get 95% CIs

# NegBin model, for natural and for link scale of both lambda and p
modavgPred(cand.set = list(fm5NB), newdata=lamNewData, parm.type = "lambda",
    type = "response", c.hat = 1.79)
modavgPred(cand.set = list(fm5NB), newdata=lamNewData, parm.type = "lambda",
    type = "link", c.hat = 1.79)    # Could be used to get 95% CIs
modavgPred(cand.set = list(fm5NB), newdata=pNewData, parm.type = "detect",
    type = "response", c.hat = 1.79)
modavgPred(cand.set = list(fm5NB), newdata=pNewData, parm.type = "detect",
    type = "link", c.hat = 1.79)         # Could be used to get 95% CIs

# ZIP model, for natural and for link scale of both lambda and p
modavgPred(cand.set = list(fm5ZIP), newdata=lamNewData, parm.type = "lambda",
    type = "response", c.hat = 2.47)
modavgPred(cand.set = list(fm5ZIP), newdata=lamNewData, parm.type = "lambda",
    type = "link", c.hat = 2.47)     # Not yet implemented (May 2015)
modavgPred(cand.set = list(fm5ZIP), newdata=pNewData, parm.type = "detect",
    type = "response", c.hat = 2.47)
modavgPred(cand.set = list(fm5ZIP), newdata=pNewData, parm.type = "detect",
    type = "link", c.hat = 2.47)         # this works, so we could get 95% CIs


rlength <- seq(1, 30, 0.01)         # Vary route length from 1 to 30 kms
newData <- data.frame(elev=0, forest=0, iLength=1/rlength)
pred <- predictSE(fm5ZIP, parm.type="lambda", newdata=newData, c.hat = 2.47)
op <- par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.3)
plot(rlength, pred[[1]], type = "l", lwd = 3, col = "blue", frame = FALSE,
    xlab = "Transect length (km)", ylab = "Exposed population (lambda)",
    ylim = c(0, 16), axes = FALSE)
axis(1, at = seq(2,30,2))       ;      axis(2)
abline(v = c(1.2, 5.135, 9.4), lwd = 2)
matlines(rlength, cbind(pred[[1]]-pred[[2]], pred[[1]]+pred[[2]]),
    type = "l", lty = 1, lwd = 2, col = "gray")


sat.pred <- predictSE(fm5ZIP, parm.type="lambda",
    newdata= data.frame(elev=0, forest=0, iLength=0), c.hat = 2.47)
abline(h = sat.pred$fit, lwd = 2, lty = 2)
par(op)

# Inspect the numbers
print(cbind("Route length" = rlength, "Exp. pop" = pred[[1]],
    "Rel. exp. pop" = pred[[1]] / sat.pred$fit), 3)


# Create covariate vectors for prediction and standardise as in analysis
ep.orig <- seq(200, 2250, length.out=100)  # Elevation between 200 and 2250 m
(elev.mean <- mean(tits[,"elev"]))
(elev.sd <- sd(tits[,"elev"]))
ep <- (ep.orig - elev.mean) / elev.sd      # Standardized for prediction
fp.orig <- seq(0, 100, length.out=100)     # Forest cover between 0 and 100%
(forest.mean <- mean(tits[,"forest"]))
(forest.sd <- sd(tits[,"forest"]))
fp <- (fp.orig - forest.mean) / forest.sd  # Standardised for prediction
date.orig <- seq(1, 90, length.out=100)  # Survey date from 1 April - 1 July
(date.mean <- mean(date, na.rm = TRUE))
(date.sd <- sd(date, na.rm = TRUE))
datep <- (date.orig - date.mean) / date.sd # Standardised for prediction
dur.orig <- seq(90, 420, length.out=100)   # Survey duration from 90 - 420 min
(dur.mean <- mean(dur, na.rm = TRUE))
(dur.sd <- sd(dur, na.rm = TRUE))
durp <- (dur.orig - dur.mean) / dur.sd     # Standardised for prediction

# ~~~~~~~ modavgPred now returns a list with 7 elements; for the plots we need [2:3] ~~~~~
# Do predictions along single covariate gradient
newData1 <- data.frame(elev=ep, forest=0, iLength=0, date=0, dur=0,
    time = factor("2", levels = c("1", "2", "3")))
pred1 <- predictSE(fm5ZIP, newdata=newData1, c.hat = 2.47)
# pred2 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData1, parm.type = "detect", type = "response", c.hat = 2.47)
pred2 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData1,
    parm.type = "detect", type = "response", c.hat = 2.47)[2:3] # ~~~~~
newData3 <- data.frame(elev=0, forest=fp, iLength=0, date=0, dur=0,
    time = factor("2", levels = c("1", "2", "3")))
pred3 <- predictSE(fm5ZIP, newdata=newData3, c.hat = 2.47)
newData4 <- data.frame(elev=0, forest=0, iLength=0, date=datep, dur=0,
    time = factor("2", levels = c("1", "2", "3")))
pred4 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData4,
    parm.type = "detect", type = "response", c.hat = 2.47)[2:3] # ~~~~~
newData5 <- data.frame(elev=0, forest=0, iLength=0, date=0, dur=durp,
    time = factor("2", levels = c("1", "2", "3")))
pred5 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData5,
    parm.type = "detect", type = "response", c.hat = 2.47)[2:3] # ~~~~~
newData6 <- data.frame(elev=0, forest=0, iLength=0, date=0, dur=0,
    time = c("1", "2", "3"), stringsAsFactors=TRUE) # ~~~~ change in default for stringsAsFactors
pred6 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData6,
    parm.type = "detect", type = "response", c.hat = 2.47)[2:3] # ~~~~~

# Plot these predictions along single covariate gradient
op <- par(mfrow = c(3,2), mar = c(5,5,3,2), cex.lab = 1.3, cex.axis = 1.3)
plot(ep.orig, pred1[[1]], type = "l", lwd = 2, col = "blue",
    xlab = "Elevation (m)", ylab = "Expected abundance", las = 1,
    ylim = c(0,50), frame = FALSE)
matlines(ep.orig, cbind(pred1[[1]]-pred1[[2]], pred1[[1]]+pred1[[2]]),
    type = "l", lty = 1, lwd = 1, col = "gray")
plot(fp.orig, pred3[[1]], type = "l", lwd = 2, col = "blue",
    xlab = "Forest cover (%)", ylab = "Expected abundance", las = 1, ylim = c(0, 18), frame = FALSE)
matlines(fp.orig, cbind(pred3[[1]]-pred3[[2]], pred3[[1]]+pred3[[2]]),
    type = "l", lty = 1, lwd = 1, col = "gray")
plot(ep.orig, pred2[[1]], type = "l", lwd = 2, col = "blue",
    xlab = "Elevation (m)", ylab = "Expected detection",
    las = 1, ylim = c(0,1), frame = FALSE)
matlines(ep.orig, cbind(pred2[[1]]-pred2[[2]], pred2[[1]]+pred2[[2]]),
    type = "l", lty = 1, lwd = 1, col = "gray")
plot(date.orig, pred4[[1]], type = "l", lwd = 2, col = "blue",
    xlab = "Survey date (1 = April 1)", ylab = "Expected detection",
    las = 1, ylim = c(0,1), frame = FALSE)
matlines(date.orig, cbind(pred4[[1]]-pred4[[2]], pred4[[1]]+pred4[[2]]),
    type = "l", lty = 1, lwd = 1, col = "gray")
plot(dur.orig, pred5[[1]], type = "l", lwd = 2, col = "blue",
    xlab = "Survey duration (min)", ylab = "Expected detection",
    las = 1, ylim = c(0,1), frame = FALSE)
matlines(dur.orig, cbind(pred5[[1]]-pred5[[2]], pred5[[1]]+pred5[[2]]),
    type = "l", lty = 1, lwd = 1, col = "gray")
barplot(pred6[[1]], names.arg = c("1", "2", "3"), ylim = c(0,1),
    ylab = "Expected detection", xlab = "Survey Number")
segments(c(0.7,1.9, 3.1), pred6[[1]]-pred6[[2]], c(0.7,1.9, 3.1),
    pred6[[1]]+pred6[[2]], lwd = 2)
par(op)

# Make predictions along two covariate gradients
# (1) Expected abundance (lambda) for forest and elevation
pred.matrix1 <- array(NA, dim = c(100, 100))
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(x=0, y=0, elev=ep[i], forest=fp[j], iLength=0)
    pred.matrix1[i,j] <- predict(fm5ZIP,type="state", newdata=newData)[1,1]
  }
}

# (2) Expected detection (p) for elevation and survey date
pred.matrix2 <- array(NA, dim = c(100, 100))
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(elev=ep[i], date=datep[j], dur=0,
        time = factor("2", levels = c("1", "2", "3")))
    pred.matrix2[i,j] <- predict(fm5ZIP, type="det", newdata=newData)[1,1]
  }
}

# (3) Expected detection (p) for elevation and survey duration
pred.matrix3 <- array(NA, dim = c(100, 100))
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(elev=ep[i], date=0, dur=durp[j],
        time = factor("2", levels = c("1", "2", "3")))
    pred.matrix3[i,j] <- predict(fm5ZIP, type="det", newdata=newData)[1,1]
  }
}

# Plot these prediction matrices
op <- par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

image(x=ep.orig, y=fp.orig, z=pred.matrix1, col = mapPalette(100), axes = FALSE,
    xlab = "Elevation (m)", ylab = "Forest cover (%)")
contour(x=ep.orig, y=fp.orig, z=pred.matrix1, add = T, col = "blue",
    labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250))
axis(2, at = seq(0, 100, by = 10))
box()
title(main = "(A)", font.main = 2)
points(tits$elev, tits$forest, pch="+", cex=1.5)

image(x=ep.orig, y=date.orig, z=pred.matrix2, col = mapPalette(100), axes = FALSE,
    xlab = "Elevation (m)", ylab = "Date (1 = April 1)")
contour(x=ep.orig, y=date.orig, z=pred.matrix2, add = T, col = "blue",
    labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250))
axis(2, at = seq(10, 120, by = 10))
box()
title(main = "(B)", font.main =2)
matpoints(tits$elev, date, pch="+", cex=1.5)

image(x=ep.orig, y=dur.orig, z=pred.matrix3, col = mapPalette(100), axes = FALSE,
    xlab = "Elevation (m)", ylab = "Duration (min)")
contour(x=ep.orig, y=dur.orig, z=pred.matrix3, add = T, col = "blue",
    labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250))
axis(2, at = seq(90, 420, by = 20))
box()
title(main = "(C)", font.main = 2)
matpoints(tits$elev, dur, pch="+", cex=1.5)
par(op)

data(Switzerland)             # Load Swiss landscape data in unmarked
CH <- Switzerland

# Predictions for lambda, with overdispersion
newData <- data.frame(elev = (CH$elev-elev.mean)/elev.sd,
    forest = (CH$forest-forest.mean)/forest.sd, iLength = 0, date=0, dur=0,
    time = factor("2", levels = c("1", "2", "3")))
predCH <- predictSE(fm5ZIP, newdata=newData, print.matrix = TRUE,
    type="response", parm.type = "lambda", c.hat = 2.47)


max(predCH[,1])                 # Look at the max prediction  --- 43.8
sum(predCH[,1] > 60)            # How many are > 60 ?  --- none
plot(CH$elev, predCH[,1])       # Relationship with elevation
predCH[1][predCH[1] > 60] <- 60 # Censor freak predicions (not req'd here)

# Prepare Swiss coordinates and produce map
library(raster)
# library(rgdal)  # ~~~~ not necessary

# Define a new dataframe with coordinates and outcome to be plotted
PARAM1 <- data.frame(x = CH$x, y = CH$y, z = predCH[,1])

# Convert the DataFrame into a raster object
r1 <- rasterFromXYZ(PARAM1)

# Create mask for elevation (mask areas > 2250 m)
elev <- rasterFromXYZ(cbind(CH$x, CH$y,CH$elevation))
elev[elev > 2250] <- NA
r1 <- mask(r1, elev)

# Create custom color palette
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

# Map expected abundance of great tits in Switzerland in 2013
op <- par(mfrow = c(1,2), mar = c(1,1,2,4))
plot(r1, col = mapPalette(100), axes = FALSE, box = FALSE, main ="")
# ~~~~~ these shape files were not distributed ~~~~~~~~~~~~~~
# lakes <- readOGR(".", "lakes")
# rivers <- readOGR(".", "rivers")
# border <- readOGR(".", "border")
# plot(rivers, col = "dodgerblue", add = TRUE)
# plot(border, col = "transparent", lwd = 1.5, add = TRUE)
# plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prepare raster with prediction SEs
r2 <- rasterFromXYZ(data.frame(x = CH$x, y = CH$y, z = predCH[,2]))
elev <- rasterFromXYZ(cbind(CH$x, CH$y,CH$elevation))
elev[elev > 2250] <- NA
r2 <- mask(r2, elev)

# Map prediction SEs of expected abundance
plot(r2, col = mapPalette(100), axes = FALSE, box = FALSE, main ="")
# ~~~~~ these shape files were not distributed ~~~~~~~~~~~~~~
# lakes <- readOGR(".", "lakes")
# rivers <- readOGR(".", "rivers")
# border <- readOGR(".", "border")
# plot(rivers, col = "dodgerblue", add = TRUE)
# plot(border, col = "transparent", lwd = 1.5, add = TRUE)
# plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(op)

# Predictions for p with overdispersion (very slow!)
newData <- data.frame(elev = (CH$elev-elev.mean)/elev.sd, date=0, dur=0,
    time = factor("2", levels = c("1", "2", "3")))
predCHp <- modavgPred(cand.set = list(fm5ZIP), newdata = newData,
    parm.type = "detect", type = "response", c.hat = 2.47)

(N <- sum(predCH[-which(CH$elev > 2250),1]))    # Nat'l population size

out <- which(CH$water > 50 | CH$elev > 2250)
(N <- sum(predCH[-out,1]))  # 'Terrestrial' population size < 2250 m elevation


# Remind ourselves of the relevant coefficients in the model
cbind(coef(fm5ZIP)[25])     # Zero-inflation model
cbind(coef(fm5ZIP)[1:8])    # Abundance model


pelev <- (CH$elev - elev.mean)/elev.sd
pforest <- (CH$forest - forest.mean)/forest.sd


Nhat <- function(fm = fm5ZIP, iLength = 0, area = 1) {
  betavec <- coef(fm)[1:8]
  DM <- cbind(rep(1,length(pelev)), pelev, pelev^2, pforest,
  pforest^2, rep(iLength,length(pelev)), pelev*pforest, pelev*pforest^2)
  pred <- exp(DM %*% betavec) * (1-plogis(coef(fm)[25])) * (1/area)
  pred2 <- pred
  N1 <- sum(pred[-which(CH$elev > 2250),])# Drop quads > 2250 m
  pred2[pred2 > 100] <- 100      # Censor freak-high estimates
  N2 <- sum(pred2[-which(CH$elev > 2250),]) # Estimate with freaks censored
  out <- which(CH$water > 50 | CH$elev > 2250)
  N3 <- sum(pred2[-out,])  # Estimate excluding water bodies
  return(c(N1 = N1, N2 = N2, N3 = N3))
}


(Nest <- Nhat(fm5ZIP))            # For default saturation density
(Nest <- Nhat(fm5ZIP, 1/10))      # For 10 km routes

# Launch the bootstrap (takes about 30h)
# ~~~~~~~~~~ changes in unmarked::parboot wef 0.12-0  ~~~~~~~
# parboot now has a 'parallel' argument with default TRUE. Unfortunately the 'Nhat' function
#  defined above won't work  in parallel, so need to specify parallel = FALSE.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# system.time(pb.N <- parboot(fm5ZIP, Nhat, nsim = 2500, report=5, parallel=FALSE))
system.time(pb.N <- parboot(fm5ZIP, Nhat, nsim = 25, report=5, parallel=FALSE)) # ~~~ use for testing, 15 mins

bs <- pb.N@t.star[,3]             # Extract the bootstrapped vals of N3

# Sample statistics
summary(bs)
quantile(bs, prob = c(0.025, 0.975))  # Get 95% CI

# Plot
hist(bs, breaks = 100, col = "grey", xlab = "National population size of great tits",
    cex = 1.5, main = "", xlim = c(400000, 1000000))
abline(v = Nest[3], lwd = 3)
abline(v = mean(bs), lwd = 3, col = "blue")
abline(v = quantile(bs, prob = c(0.025, 0.975)), lwd = 3, col = "blue", lty = 2)



(CI1 <- quantile(bs, prob = c(0.025, 0.975)))     # Percentile-based CI
(CI2 <- c(mean(bs-2*sd(bs)), mean(bs+2*sd(bs))))  # Normal approx.-based CI


c.hat <- gof.ZIP$c.hat.est
(CI <- c(mean(bs-2*sqrt(c.hat)*sd(bs)), mean(bs+2*sqrt(c.hat)*sd(bs))))


# 6.9.5 Conclusions on the analysis with unmarked (no code)
# ------------------------------------------------------------------------

# 6.10 The issue of space, or: what is your effective sample area ?
# =================================================================

AHR <- seq(0.001, 100, 0.1)         # Area home range (in km2): 0.1-100
RHR <- sqrt(AHR/pi)                 # Translate to home range radius (in km)
ESA <- (1 + 2*RHR)^2                # Eff. sample area for 1km2 nominal area
op <- par(mfrow = c(1,2), mar = c(5,5,3,2), cex.lab = 1.3, cex.axis = 1.3)
plot(AHR, ESA, xlab = "Home range size (km2)", ylab = "Effective sample area (km2)",
    type = "l", lwd = 3, frame = FALSE)
abline(h = 1, col = "red", lwd = 3) # Nominal sample area
abline(h = 0, col = "grey", lwd = 3)


GT.HR <- seq(0.003, 0.1,,1000)     # Great tit home ranges in km2
GT.rad <- sqrt(GT.HR/pi)           # Great tit home range radius (in km)
ESA.GT <- (1 + 2*GT.rad)^2         # Effective sample area for Great tit
NTOT <- numeric(length(ESA.GT))    # Adjusted national total population
for(i in 1:length(NTOT)){
   NTOT[i] <- Nhat(fm5ZIP, 0, ESA.GT[i])[3]
}
plot(100*GT.HR, NTOT, xlab = "Great tit home-range (ha)",
    ylab = "Adjusted national total (pairs)", type = "l", lwd = 3, frame = FALSE)
par(op)

