#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 11. Hierarchical models for communities
# =========================================================================

# 11.9 Species richness maps and species accumulation curves
# ------------------------------------------------------------------------


# Get Swiss landscape data and standardise covariates as for model 10
library(unmarked)
data(Switzerland)
ch <- Switzerland
ELE <- (ch$elevation - mean.ele) / sd.ele
FOREST <- (ch$forest - mean.forest) / sd.forest

nsamp <- nrow(all10)            # 1200   ..... far too many
nkm2 <- length(ch[[1]])         # 42275, that's a LOT!
select.samp <- sort(sample(1:nsamp, 50)) # Chose random sample of 50
nsamp <- length(select.samp)    # new sample size 50

# Create posterior predictive distribution for Z for Swiss landscape
str( zCH <- array(NA, dim = c(nkm2, 215, nsamp)) ) # BIG array !
W <- all10[,1722:1936]          # Grab MCMC samples from w
LPSI <- all10[,1507:1721]       # Grab MCMC samples from logit(psi)
BETALPSI1 <- all10[,646:860]    # Grab MCMC samples from betalpsi1
BETALPSI2 <- all10[,861:1075]   # Grab MCMC samples from betalpsi2
BETALPSI3 <- all10[,1076:1290]  # Grab MCMC samples from betalpsi3
for(i in 1:nkm2){               # takes about 5 mins !
   cat(paste("\nQuadrat", i, "\n"))
   for(u in 1:length(select.samp)){
      psi <- W[select.samp[u],] * plogis(LPSI[select.samp[u],] +
         BETALPSI1[select.samp[u],] * ELE[i] +
         BETALPSI2[select.samp[u],] * ELE[i]^2 +
         BETALPSI3[select.samp[u],] * FOREST[i] )
      zCH[i,,u] <- rbinom(215, 1, psi)
   }
}

# Compute posterior distribution of species richness by collapsing z array
SR <- apply(zCH, c(1,3), sum)   # posterior distribution
pmSR <- apply(SR, 1, mean)      # posterior mean
sdSR <- apply(SR, 1, sd)        # posterior standard deviation


library(raster)
library(rgdal)
par(mfrow = c(1,2), mar = c(2,2,3,5))
# Posterior mean map
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pmSR))
elev <- rasterFromXYZ(cbind(ch$x, ch$y,ch$elevation))
elev[elev > 2250] <- NA         # Mask areas > 2250 m a.s.l.
r1 <- mask(r1, elev)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette(100), axes = F, box = FALSE, main ="")
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)

# Posterior standard deviation map
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = sdSR))
elev <- rasterFromXYZ(cbind(ch$x, ch$y,ch$elevation))
elev[elev > 2250] <- NA         # Mask areas > 2250 m a.s.l.
r1 <- mask(r1, elev)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette(100), axes = F, box = FALSE, main ="")
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)


# Get 3,000 posterior samples of omega, and the mean and sd hyperparameters
omega <- out101$sims.list$omega
mu.lpsi <- out101$sims.list$mu.lpsi
str( sd.lpsi <- out101$sims.list$sd.lpsi )    # Confirms we have 3,000 draws

# compute posterior predictions of species occurrence probabilities
nsites <- 100
ndraws <- length(omega)
Nmax <- 215
psi <- matrix(NA, nrow=ndraws, ncol=Nmax)
for (i in 1:ndraws) {
   w <- rbinom(215, 1, omega[i])
   psi[i,] <- w * plogis(rnorm(Nmax, mean = mu.lpsi[i], sd=sd.lpsi[i]))
}

# compute posterior predictions of species presence at each site
z <- array(NA, dim=c(ndraws, Nmax, nsites))
for (i in 1:ndraws) {
  for (j in 1:Nmax) {
    z[i,j, ] <- rbinom(nsites, size=1, prob=psi[i,j])
  }
}

# compute posterior predictions of cumulative number of species present
Ntot <- matrix(NA, nrow=ndraws, ncol=nsites)
for (i in 1:ndraws) {
  for (j in 1:nsites) {
    zsum <- rep(NA, Nmax)
    if (j>1) {
      zsum <- apply(z[i, , 1:j], 1, sum)
    }
    else {
      zsum <- z[i, , 1]
    }
    Ntot[i,j] <- sum(zsum>0)
  }
}                        # takes about 4 min

# compute summary stats of species accumulation curve
nSpeciesPresent <- matrix(NA, nrow=3, ncol=nsites)
for (j in 1:nsites) {
  x <- Ntot[,j]
  nSpeciesPresent[1, j] <- mean(x)
  nSpeciesPresent[2:3, j] <- quantile(x, probs=c(0.05, 0.95))
}

# Plot species accumulation curve
ylim = c(min(nSpeciesPresent[2,]), max(nSpeciesPresent[3,]))
plot(1:nsites, nSpeciesPresent[1,], pch=16, ylim=ylim, type="b",
xlab="Number of sample locations", ylab="Number of occurring species",
las=1, cex.axis=1.2, cex.lab=1.5, cex=1.2, frame = F)
segments(1:nsites, nSpeciesPresent[2,], 1:nsites, nSpeciesPresent[3,])

