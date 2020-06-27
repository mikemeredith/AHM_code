#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7. Modeling abundance using multinomial N-mixture models
# =========================================================================

library(unmarked)
library(AHMbook)

# ~~~~~~~~~~ section 7.8.5 uses the first data set from section 7.3 ~~~~~~~~~~~
set.seed(2014)
data <- simNmix(mean.lam = exp(1), beta3.lam = 1, mean.p = plogis(0),
       sigma.p.visit = 1, show.plot=FALSE)
str(data)
# Get detection history frequencies for each site (for exactly 3 surveys)
dhfreq <- array(NA, dim = c(data$nsite, 7),
  dimnames = list(NULL, c("100", "010", "001", "110", "101", "011", "111")))
for(i in 1:data$nsite){
  dhfreq[i,] <- table(factor(paste(data$DH[i,1,], data$DH[i,2,],
  data$DH[i,3,], sep = ""),
  levels = c("100", "010", "001", "110", "101", "011", "111")))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 7.8 Spatially Stratified Capture-Recapture Models
# =================================================

# 7.8.5 Example 2: analysis of data simulated with the simNmix data function
# --------------------------------------------------------------------------
# Fit model in unmarked
library(unmarked)
time <- matrix(as.character(1:3), data$nsite, 3, byrow = TRUE)

# Define pifun for J=3 occasion capture-recapture protocol
crPiFun <- function(p) {
   p1 <- p[,1] # Extract the columns of the p matrix, one for
   p2 <- p[,2] # each of J = 3 sample occasions
   p3 <- p[,3]
   cbind( # define multinomial cell probabilities:
   "100" = p1 * (1-p2) * (1-p3),
   "010" = (1-p1) * p2 * (1-p3),
   "001" = (1-p1) * (1-p2) * p3,
   "110" = p1 * p2 * (1-p3),
   "101" = p1 * (1-p2) * p3,
   "011" = (1-p1) * p2 * p3,
   "111" = p1 * p2 * p3)
}

# Define mapping function for missing values
o2y <- matrix(1, 3, 7)

# Create unmarked frame and fit couple of models
umf <- unmarkedFrameMPois(y = dhfreq, siteCovs = data.frame(cov3 = data$site.cov[,3]),
  obsCovs = list(time = time), obsToY = o2y, piFun = "crPiFun")
fm1 <- multinomPois(~1 ~1, umf)    # detection model before abundance model
fm2 <- multinomPois(~time-1 ~1, umf)
fm3 <- multinomPois(~1 ~cov3, umf)
fm4 <- multinomPois(~time-1 ~cov3, umf)

# Assemble the models into a fitList and rank using AIC
ms <- fitList(
    "lam(.)p(.)" = fm1,
    "lam(.)p(time)" = fm2,
    "lam(cov3)p(.)" = fm3,
    "lam(cov3)p(time)" = fm4)

(AICtable <- modSel(ms))

summary(fm4)
p.true <- qlogis(data$p[min(which(data$N>0)),,1])
tmp <- cbind(rbind(lam0 = log(data$mean.lam), beta3 = data$beta3.lam,
    logit.p1 = p.true[1], logit.p2 = p.true[2], logit.p3 = p.true[3]), coef(fm4))
colnames(tmp) <- c("Truth", "MLEs")
tmp


# No temporal variation in p and effect of site-covariate on lambda
set.seed(24)
data <- simNmix(mean.lam = exp(1), beta2.lam = 1, beta3.lam = 1,
    mean.p = plogis(0.2), beta3.p = -1, beta.p.survey = -1)
str(data$DH)

# Get occasions with first detection of each individual
f <- apply(data$DH, c(1,3), function(x) min(which(x != 0)))
head(f)   ;   str(f)

# Produce removal counts (for any number of occasions)
y <- array(NA, dim = c(data$nsite, data$nvisit),
    dimnames = list(NULL, as.factor(1:data$nvisit)))
for(i in 1:data$nsite){
  y[i,] <- table(factor(f[i,], levels = as.character(1:data$nvisit)))
}
y                # Look at removal data set

# Fit models in unmarked
summary(umf <- unmarkedFrameMPois(y = y,
    siteCovs = data.frame(cov2 = data$site.cov[,2], cov3 = data$site.cov[,3]),
    obsCovs = list(obs.cov = data$survey.cov),
    type = "removal"))  # Create and look at um data frame

fm1 <- multinomPois(~1 ~1, umf)    # Detection model before abundance model
(fm4 <- multinomPois(~cov3 + obs.cov ~cov2 + cov3, umf))

print(c(log(data$mean.lam), data$beta2.lam, data$beta3.lam,
       # logit(data$mean.p), data$beta3.p, data$beta.p.survey))
       qlogis(data$mean.p), data$beta3.p, data$beta.p.survey))


# Simulate detection frequency data  from back in section 7.3.1
set.seed(24)
data <- simNmix(mean.lam = exp(0), beta2.lam = 1, beta3.lam =-1,
   beta4.lam = 0.2, dispersion = 1, mean.p = plogis(0),
   beta3.p = 1, sigma.p.visit = 1, Neg.Bin = TRUE)
dhfreq <- array(NA, dim = c(data$nsite, 7),
  dimnames = list(NULL, c("100", "010", "001", "110", "101", "011", "111")))

for(i in 1:data$nsite){
  dhfreq[i,] <- table(factor(paste(data$DH[i,1,], data$DH[i,2,],
      data$DH[i,3,], sep = ""),
      levels = c("100", "010", "001", "110", "101", "011", "111")))
}
dhfreq                     # Look at resulting data set

# Bundle data in unmarked frame
time <- matrix(as.character(1:3), data$nsite, 3, byrow = TRUE)
summary(umf <- unmarkedFrameGMM(y = dhfreq, numPrimary = 1,
       siteCovs = data.frame(cov2 = data$site.cov[,2],
       cov3 = data$site.cov[,3], cov4 = data$site.cov[,4]),
       obsCovs = list(time = time), obsToY = o2y, piFun = "crPiFun"))

# Fit a couple of models, first for detection
fm1 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
    mix = "NB", data = umf)
fm2 <- gmultmix(~1, ~1, ~time-1, mix = "NB", data = umf)
fm3 <- gmultmix(~1, ~1, ~time-1+cov3, mix = "NB", data = umf)

# ... then for abundance,
fm4 <- gmultmix(~cov2, ~1, ~1, mix = "NB", data = umf)
fm5 <- gmultmix(~cov2+cov3, ~1, ~1, mix = "NB", data = umf)
fm6 <- gmultmix(~cov2+cov3+cov4, ~1, ~1, mix = "NB", data = umf)

# ... and the data-generating model
fm7 <- gmultmix(~cov2+cov3+cov4, ~1, ~time-1+cov3, mix = "NB", data = umf)

# Compare models with AIC
ms <- fitList(
    "lam(.)p(.)" = fm1,
    "lam(.)p(time)" = fm2,
    "lam(.)p(time+cov3)" = fm3,
    "lam(cov2)p(.)" = fm4,
    "lam(cov2+cov3)p(.)" = fm5,
    "lam(cov2+cov3+cov4)p(.)" = fm6,
    "lam(cov2+cov3+cov4)p(time+cov3)" = fm7)

(AICtable <- modSel(ms))

# Compare data-generation truth with estimates
p.true <- qlogis(data$mean.p) + data$eta.p.visit
Truth <- rbind(lam0 = log(data$mean.lam), beta2.lam = data$beta2.lam,
    beta3.lam = data$beta3.lam, beta4.lam = data$beta4.lam,
    logit.p1 = p.true[1], logit.p2 = p.true[2], logit.p3 = p.true[3],
    beta3.p = data$beta3.p, log.dispersion = log(data$dispersion))
MLEs <- coef(fm7)
print(cbind(Truth, MLEs), 3)

