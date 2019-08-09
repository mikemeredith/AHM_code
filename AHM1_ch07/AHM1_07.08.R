#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
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


crPiFun <- function(p) {
   p1 <- p[,1]
   p2 <- p[,2]
   p3 <- p[,3]
   cbind("001" = (1 - p1) * (1 - p2) *      p3,
         "010" = (1 - p1) *      p2  * (1 - p3),
         "011" = (1 - p1) *      p2  *      p3,
         "100" =      p1  * (1 - p2) * (1 - p3),
         "101" =      p1  * (1 - p2) *      p3,
         "110" =      p1  *      p2  * (1 - p3),
         "111" =      p1  *      p2  *      p3)
}


p <- matrix(0.4, 2, 3)
crPiFun(p)


# To compute pi0 we do this:
(pi0 <- 1 - rowSums(crPiFun(p)))


# 7.8.1 Example 2: Fitting Models M0, Mt, and Mx to Chandler’s Flycatcher data
# ------------------------------------------------------------------------
alfl <- read.csv(system.file("csv", "alfl.csv", package="unmarked"))
head(alfl, 5)

alfl.covs <- read.csv(system.file("csv", "alflCovs.csv",package="unmarked"),
        row.names=1)
head(alfl.covs)

alfl$captureHistory <- paste(alfl$interval1, alfl$interval2, alfl$interval3,
        sep="")
alfl$captureHistory <- factor(alfl$captureHistory,
     levels=c("001", "010", "011", "100", "101", "110", "111"))
alfl$id <- factor(alfl$id, levels=rownames(alfl.covs))

alfl.v1 <- alfl[alfl$survey==1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
head(alfl.H1, 5)

intervalMat <- matrix(c('1','2','3'), 50, 3, byrow=TRUE)
class(alfl.H1) <- "matrix"
o2y <- matrix(1, 3, 7)
umf.cr1 <- unmarkedFrameMPois(y=alfl.H1,
   siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
   obsCovs=list(interval=intervalMat), obsToY=o2y, piFun="crPiFun")
summary(umf.cr1)

M0 <- multinomPois(~ 1 ~ 1, umf.cr1)
Mt <- multinomPois(~ interval - 1 ~ 1, umf.cr1)
Mx <- multinomPois(~ time.1 ~ 1, umf.cr1)

(M0.woody <- multinomPois(~ 1 ~ woody, umf.cr1))

fl <- modSel(fitList(M0, Mt, Mx, M0.woody))

nd <- data.frame(woody=seq(0, 0.8, length=50))
E.abundance <- predict(M0.woody, type="state", newdata=nd, appendData=TRUE)
plot(Predicted ~ woody, E.abundance, type="l", ylim=c(0, 6),
      ylab="Alder flycatchers / plot", xlab="Woody vegetation cover",
      frame = F)
 lines(lower ~ woody, E.abundance, col=gray(0.7))
 lines(upper ~ woody, E.abundance, col=gray(0.7))


backTransform(M0.woody, type="det")

 round(getP(M0.woody), 2)[1,]


# 7.8.2 Models with behavioral response
# ------------------------------------------------------------------------
crPiFun.Mb <- function(p) {
 pNaive <- p[,1]
 pWise <- p[,3]
 cbind("001" = (1 - pNaive) * (1 - pNaive) *      pNaive,
       "010" = (1 - pNaive) *      pNaive  * (1 - pWise),
       "011" = (1 - pNaive) *      pNaive  *      pWise,
       "100" =      pNaive  * (1 - pWise)  * (1 - pWise),
       "101" =      pNaive  * (1 - pWise)  *      pWise,
       "110" =      pNaive  *      pWise   * (1 - pWise),
       "111" =      pNaive  *      pWise   *      pWise)
}

behavior <- matrix(c('Naive', 'Naive', 'Wise'), 50, 3, byrow=TRUE)
umf.cr1Mb <- unmarkedFrameMPois(y=alfl.H1,
   siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
   obsCovs=list(behavior=behavior),  obsToY=o2y, piFun="crPiFun.Mb")
M0 <- multinomPois(~1 ~1, umf.cr1Mb)
(Mb <- multinomPois(~behavior-1 ~1, umf.cr1Mb))


# 7.8.3 Models with individual heterogeneity model
# ------------------------------------------------------------------------
parID <- matrix(c('p','sig','sig'), 50, 3, byrow=TRUE)
umf.cr2 <- unmarkedFrameMPois(y=alfl.H1,
        siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
        obsCovs=list(parID=parID), obsToY=o2y, piFun="MhPiFun")



# 7.8.4 Bayesian analysis using data augmentation (DA): heterogeneity models in BUGS
# ----------------------------------------------------------------------------------
# Extract data and do data augmentation up to M = 400
y <- as.matrix(alfl[,c("interval1","interval2","interval3")] )
nind <- nrow(y)
M <- 400
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=3))

# Site ID
# Make site ID into an integer: This only works ok here because sites are in
#  alphabetical order in the data set!
site <- as.numeric(alfl$id)
site <- c(site, rep(NA, M-nind))

# Next we extract the covariates and standardize them
sitecovs <- scale(as.matrix(alfl.covs[,c("woody", "struct", "time.1")]))

# Bundle data for BUGS
data <- list(y = y, J = 3, M = M , nsites = 50, X = sitecovs, group=site)
str(data)

# Specify model in BUGS language
cat("
model {

# Prior distributions
p0 ~ dunif(0,1)
alpha0 <- log(p0 / (1-p0))    # same as logit(p0)
alpha1 ~ dnorm(0, 0.01)
alpha2 ~ dnorm(0,0.01)
beta0 ~ dnorm(0,0.01)
beta1 ~ dnorm(0,0.01)
psi <- sum(lambda[]) / M   # psi is a derived parameter

# log-linear model for abundance: lambda depends on WOODY
for(s in 1:nsites){
  log(lambda[s]) <- beta0 + beta1 * X[s,1]
  probs[s] <- lambda[s] / sum(lambda[])
}

# Model for individual encounter histories
for(i in 1:M){
  group[i] ~ dcat(probs[])  # Group == site membership
  z[i] ~ dbern(psi)         # Data augmentation variables

  # Observation model: p depends on 2 covariates: STRUCT + TIME
  for(j in 1:J){
    logit(p[i,j]) <- alpha0 + alpha1 * X[group[i],2] + alpha2*X[group[i],3]
    pz[i,j] <- p[i,j] * z[i]
    y[i,j] ~ dbern(pz[i,j])
  }
}
}
",fill=TRUE,file="model.txt")

# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "alpha2", "beta0", "beta1", "psi")

# Initial values
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), alpha2 = rnorm(1),beta0=runif(1),
    beta1=rnorm(1), z= c( rep(1,100), rep(0, 300))) }

# MCMC settings
ni <- 11000   ;   nb <- 1000   ;   nt <- 4   ;   nc <- 3

# Call JAGS from R and summarize marginal posteriors
library("jagsUI")
out <- jags(data, inits, parameters, "model.txt", n.thin = nt,
   n.chains = nc, n.burnin = nb, n.iter = ni)
print(out, digits = 3)


# Specify model in BUGS language
cat("
model {

# Prior distributions
p0 ~ dunif(0,1)
alpha0 <- log(p0/(1-p0))
alpha1 ~ dnorm(0, 0.01)
alpha2 ~ dnorm(0, 0.01)
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)
psi <- sum(lambda[])/M
tau ~ dgamma(0.1,0.1) # New parameter, precision of ind. random effects
sigma <- 1/sqrt(tau)

# log-linear model for abundance: lambda depends on WOODY
for(s in 1:nsites){
  log(lambda[s])<- beta0 + beta1*X[s,1]
  probs[s]<- lambda[s]/sum(lambda[])
}

# Model for individual encounter histories
for(i in 1:M){
  eta[i] ~ dnorm(alpha0, tau)  # Individual random effect
  group[i] ~ dcat(probs[])  # Group == site membership
  z[i] ~ dbern(psi)         # Data augmentation variables
  # Observation model: p depends on STRUCT + TIME + ind. heterogeneity
  for(j in 1:J){
    logit(p[i,j]) <- alpha1 * X[group[i],2] + alpha2*X[group[i],3] + eta[i]
    pz[i,j] <- p[i,j] * z[i]
    y[i,j] ~ dbern(pz[i,j])
  }
 }
}
",fill=TRUE,file="model.txt")

# Parameters monitored: add sigma
parameters <- c("p0", "alpha0", "alpha1", "alpha2", "beta0", "beta1", "psi", "sigma")

# Initial values: add tau
inits <- function(){
list (p0 = runif(1), alpha1 = runif(1), alpha2 = rnorm(1),beta0=runif(1),
      beta1=rnorm(1),z= c( rep(1,100), rep(0, 300)), tau = 1)  }

# MCMC settings (others as before)
ni <- 50000   ;   nb <- 10000   ;   nt <- 2   ;   nt <- 10

# Call JAGS from R and summarize marginal posteriors
out <- jagsUI(data, inits, parameters, "model.txt", n.thin = nt,
   n.chains = nc, n.burnin = nb, n.iter = ni)
print(out, digits = 3)


# 7.8.5 Example 2: analysis of data simulated with the simNmix data function
# --------------------------------------------------------------------------
# Fit model in unmarked
library(unmarked)
time <- matrix(as.character(1:3), data$nsite, 3, byrow = T)

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
tmp <- cbind(rbind(lam0 = log(data$mean.lam), beta3 = data$beta3.lam, logit.p1 = p.true[1], logit.p2 = p.true[2], logit.p3 = p.true[3]), coef(fm4))
colnames(tmp) <- c("Truth", "MLEs")
tmp


# No temporal variation in p and effect of site-covariate on lambda
set.seed(24)
data <- simNmix(mean.lam = exp(1), beta2.lam = 1, beta3.lam = 1, mean.p = plogis(0.2), beta3.p = -1, beta.p.survey = -1)
str(data$DH)

# Get occasions with first detection of each individual
f <- apply(data$DH, c(1,3), function(x) min(which(x != 0)))
head(f)   ;   str(f)

# Produce removal counts (for any number of occasions)
y <- array(NA, dim = c(data$nsite, data$nvisit), dimnames = list(NULL, as.factor(1:data$nvisit)))
for(i in 1:data$nsite){
   y[i,] <- table(factor(f[i,], levels = as.character(1:data$nvisit)))
}
y                # Look at removal data set

# Fit models in unmarked
summary(umf <- unmarkedFrameMPois(y = y, siteCovs = data.frame(cov2 = data$site.cov[,2], cov3 = data$site.cov[,3]), obsCovs = list(obs.cov = data$survey.cov), type = "removal"))  # Create and look at um data frame

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
time <- matrix(as.character(1:3), data$nsite, 3, byrow = T)
summary(umf <- unmarkedFrameGMM(y = dhfreq, numPrimary = 1,
       siteCovs = data.frame(cov2 = data$site.cov[,2],
       cov3 = data$site.cov[,3], cov4 = data$site.cov[,4]),
       obsCovs = list(time = time), obsToY = o2y, piFun = "crPiFun"))

# Fit a couple of models, first for detection
fm1 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1, mix = "NB", data = umf)
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
Truth <- rbind(lam0 = log(data$mean.lam), beta2.lam = data$beta2.lam, beta3.lam = data$beta3.lam, beta4.lam = data$beta4.lam, logit.p1 = p.true[1], logit.p2 = p.true[2], logit.p3 = p.true[3], beta3.p = data$beta3.p, log.dispersion = log(data$dispersion))
MLEs <- coef(fm7)
print(cbind(Truth, MLEs), 3)

