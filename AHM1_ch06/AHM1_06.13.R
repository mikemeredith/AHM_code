#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# Approximate execution time for this code: 35 mins

library(AHMbook)
library(unmarked)

# 6.13 The Royle-Nichols model and other non-standard N-mixture models
# ====================================================================


# 6.13.1 The Royle-Nichols or Poisson/Bernoulli N-mixture model
# ------------------------------------------------------------------------
playRN <- function(M = 267, J = 3, mean.abundance = 1, mean.detection = 0.3){
  # Function generates replicated count data under the Nmix model of Royle (2004),
  #   then 'degrades' the data to detection/nondetection and fits the RN model
  #   (Royle & Nichols 2003) using unmarked and estimates site-specific abundance.
  #   Requires function simNmix and package unmarked.
  #
  devAskNewPage(ask = FALSE)
  #
  # Simulate Nmix data under a range of abundance levels
  data <- simNmix(nsite = M, nvisit = J, mean.lam = mean.abundance, mean.p = mean.detection, beta2.lam = 1, beta3.p = -1, beta.p.survey = -1, show.plot = FALSE)
  # Turn counts into detection/nondetection data
  y <- data$C          # Copy counts C into y
  y[y>0] <- 1          # Turn counts >0 into 1
  # Load unmarked, format data and summarize
  library(unmarked)
  umf <- unmarkedFrameOccu(y=y, siteCovs= data.frame(cov2 = data$site.cov[,2], cov3 = data$site.cov[,3]), obsCovs = list(obscov = data$survey.cov))
  # Fit data-generating model
  fm <- occuRN(~cov3+obscov ~cov2, data=umf)
  # Estimate local abundance N and plot against true N (known in simulation)
  Nest <- bup(ranef(fm, K = ), "mean")
  op <- par(mfrow = c(1,1))
  plot(data$N, Nest, xlab = "True local abundance", ylab = "Estimated local abundance", frame = F)
  abline(0,1, lwd = 3)                              # 1:1 line
  abline(lm(Nest ~ data$N), col = "blue", lwd = 3)  # Regression
  par(op)
  slope <- coef(lm(Nest ~ data$N))[2]               # Is 1 if model perfect
  return(list(nsite = M, nvisit = J, coef = coef(fm), slope = slope))
}

# Execute the function using various settings
playRN(M = 100, J = 3, mean.abundance = 0.1)  # Increasing abundance
playRN(M = 100, J = 3, mean.abundance = 1)
playRN(M = 100, J = 3, mean.abundance = 5)
playRN(M = 100, J = 3, mean.detection = 0.3)  # Increasing detection
playRN(M = 100, J = 3, mean.detection = 0.5)
playRN(M = 100, J = 3, mean.detection = 0.7)
playRN(M = 100, J = 20)                       # More visits
playRN(M = 1000, J = 3)                       # More sites


# Run simulation
lam <- c(0.1, 0.5, 1, 2.5, 5, 10)    # 6 levels of mean abundance
simrep <- 100
results <- array(NA, dim = c(length(lam), simrep, 6))
for(i in 1:6){
  for(j in 1:simrep){
    cat(paste("\n *** lambda level", lam[i], ", simrep", j, "***\n"))
    tmp <- playRN(mean.abundance = lam[i])
    results[i,j,1:5] <- tmp$coef   # Coefficients of RN model
    results[i,j,6] <- tmp$slope    # Slope of regression of Nest on Ntrue
  }
}

# Summary of results for abundance (Fig. 6-23)
op <- par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
boxplot(t(exp(results[,,1])), names = as.character(lam), outline = FALSE,
    frame = FALSE, col = "grey", xlab = "Mean abundance (lambda)",
    ylab = "lambda intercept", ylim = c(0,10))
points(1:6, lam, pch = "*", col = "red", cex = 3)
boxplot(t(results[,,2]), names = as.character(lam), outline = FALSE,
    frame = FALSE, col = "grey", xlab = "Mean abundance (lambda)",
    ylab = "lambda slope", ylim = c(0.5, 1.5))
abline(h = 1, col = "red", lwd = 3)
par(op)

# Load data on Swiss tits
## Code modified to use the SwissTits data set included in the AHMbook package
data(SwissTits)
?SwissTits
str(SwissTits)

# Select the count data for 2013 (all species)
y0 <- SwissTits$count[, , '2013', ]
str(y0)
# We keep the sites with count data, remove those with 3 NAs
# See which sites have counts in 2013 for (say) Great tits:
keep <- which(rowSums(is.na(y0[, , "Great tit"])) != 3)
length(keep)
y3D <- y0[keep, , ]

# Get covariate data (site and observational) and drop unsurveyed sites
elev <- SwissTits$sites$ele[keep]
route <- SwissTits$sites$rlength[keep]
forest <- SwissTits$sites$forest[keep]
date <- SwissTits$date[keep, , '2013']  # Survey date
dur <- SwissTits$dur[keep, , '2013']    # Survey duration

# 'Degrade' counts to mere detection/nondetection data
# y3DRN <- y ; y3DRN[y3DRN > 1] <- 1  # ~~~~~~~~~ this is wrong
y3DRN <- y3D ; y3DRN[y3DRN > 1] <- 1  # Overwrite any count >1 with 1 (for RN model)

( spec.names <- paste0(SwissTits$species$name, "s") )

library(unmarked)

# Loop over 6 species of tits
op <- par(mfrow = c(2,2), mar = c(5,4,3,1))
oldask <- devAskNewPage(ask=dev.interactive(orNone=TRUE)) #~~~~ to replace browser call
for(k in 1:6){
  cat("\n*** Analysis for ", spec.names[k], "***\n")
  # Plot observed data: counts vs survey date
  matplot(t(date), t(y3D[,,k]), type = "l", lwd = 3, lty = 1, frame = FALSE,
      xlab = "Survey date (1 = April 1)", ylab = "Observed counts",
      main = paste("Counts of", spec.names[k], "as a function of survey date"))

  # Fit standard Nmix model (Nmix1)
  time <- matrix(rep(as.character(1:3), 263), ncol = 3, byrow = T)
  summary(umf1 <- unmarkedFramePCount(y = y3D[,,k], siteCovs=data.frame(elev=scale(elev),
      forest=scale(forest), iLength=1/route), obsCovs=list(time = time, date = scale(date),
      dur = scale(dur))) )
  Nmix1 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
        ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength,
        umf1, control=list(trace=TRUE, REPORT=5, maxit = 250))
  (tmp1 <- summary(Nmix1))

  # Fit RN model (Nmix2)
  summary(umf2 <- unmarkedFrameOccu(y = y3DRN[,,k], siteCovs=data.frame(elev=scale(elev),
      forest=scale(forest), iLength=1/route), obsCovs=list(time = time, date = scale(date), dur = scale(dur))))
  # Use solutions from Nmix1 as inits for Nmix2
  Nmix2 <- occuRN(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
        ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength,
        umf2, control=list(trace=TRUE, REPORT=5), starts = coef(Nmix1))
  (tmp2 <- summary(Nmix2))

  # Compare estimates under both models
  # Table with MLEs and SEs
  print(cbind(rbind(tmp1$state[,1:2], tmp1$det[,1:2]),
      rbind(tmp2$state[,1:2], tmp2$det[,1:2])))

  # Plot of all RN estimates versus all Nmix estimates
  plot(coef(Nmix1), coef(Nmix2), xlab = "Coefficients Nmix",
      ylab = "Coefficients RN", main = spec.names[k])
  abline(0,1, lwd = 2)
  abline(h = 0, lwd = 1, col = "grey")
  abline(v = 0, lwd = 1, col = "grey")
  abline(lm(coef(Nmix2) ~ coef(Nmix1)), lwd = 2, col = "blue")
  # browser() #~~~~~ incompatible with automated checking

  # Overall discrepancy measure (for state model only): slope of regression and r2
  print(slope <- coef(lm(coef(Nmix2)[1:10] ~ coef(Nmix1)[1:10]))[2]) # Slope
  print(r <- cor(coef(Nmix2)[1:10], coef(Nmix1)[1:10]))       # Correlation
}
devAskNewPage(oldask) #~~~~ clean up
par(op)

# 6.13.2 The Poisson/Poisson N-mixture model (no code)
# ------------------------------------------------------------------------

