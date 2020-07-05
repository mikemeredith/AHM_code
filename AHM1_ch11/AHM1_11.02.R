#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 11. Hierarchical models for communities
# =========================================================================

# Approximate execution time for this code: 6 mins
# Run time with the full number of iterations: 30 mins

library(AHMbook)

# 11.2 Simulation of a metacommunity
# ==================================


simComm(type="det/nondet", nsite=30, nrep=3, nspec=100,
    mean.psi=0.25, sig.lpsi=1, mu.beta.lpsi=0, sig.beta.lpsi=0,
    mean.lambda=2, sig.loglam=1, mu.beta.loglam=1, sig.beta.loglam=1,
    mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE)


# Execute function with default arguments
set.seed(1234)
data <- simComm(type="det/nondet", nsite=30, nrep=3, nspec=100,
    mean.psi=0.25, sig.lpsi=1, mu.beta.lpsi=0, sig.beta.lpsi=0,
    mean.lambda=2, sig.loglam=1, mu.beta.loglam=1, sig.beta.loglam=1,
    mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE)
# data <- simComm() # same

str(data)

# Some possibly interesting settings of the function
data <- simComm(nsite = 267, nspec = 190, mean.psi = 0.25, sig.lpsi = 2,
    mean.p = 0.12, sig.lp = 2) # similar to Swiss MHB
data <- simComm(mean.psi = 1)         # all species occur at every site
data <- simComm(mean.p = 1)           # no measurement error (perfect detection)

# Effect of spatial sample size (nsite) on species richness in sample (Ntotal.fs)
data <- simComm(nsite=50, nspec = 200) # 1-3 are usually missed in sample
data <- simComm(nsite=30, nspec = 200) # 4-6 usually missed
data <- simComm(nsite=10, nspec = 200) # around 30 typically missed

# Check for frequentist characteristics of such statistics
temp <- rep(NA, 100)
for(i in 1:100){
   cat("\nSimrep", i)
   temp[i] <- simComm(nsite=10, nspec = 200, show.plot = FALSE)$Ntotal.fs
}
hist(200-temp, breaks = 30,
    main = "Number of species in the metacommunity \nthat do not occur in the 10 sampled sites", col = "gold")

# Simulation 1: effects of psi and sd(logit(psi)) on number of species actually occurring in the 50 sampled sites
# simrep <- 50                   # Run 50 simulation reps
simrep <- 10                   # ~~~~ reduce number for testing
mpsi <- seq(0.01, 0.25,,10)
slpsi <- seq(0.1, 5,,10)
results1 <- array(NA, dim = c(10, 10, simrep))
for(i in 1:10){      # Loop over levels of factor mean.psi (mpsi)
  for(j in 1:10){    # Loop over levels of factor sig.lpsi (slpsi)
    for(k in 1:simrep){
      cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
      tmp <-  simComm(nsite=50, nspec = 200, show.plot = FALSE, mean.psi = mpsi[i],
        sig.lpsi = slpsi[j])
      results1[i,j,k] <- tmp$Ntotal.fs
    }
  }
}


# Simulation 2: effects of p and sd(logit(p)) on the proportion of the species occurring in the 50 sampled sites that are detected at least once
# simrep <- 50         # Run 50 simulation reps again
simrep <- 10         # Reduce number for testing
mp <- seq(0.01, 0.25,,10)
slp <- seq(0.1, 5,,10)
results2 <- array(NA, dim = c(10, 10, simrep, 2))
for(i in 1:10){      # Loop over levels of factor mean.p (mp)
  for(j in 1:10){    # Loop over levels of factor sig.lp (slp)
    for(k in 1:simrep){
      cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
      tmp <-  simComm(nsite=50, nspec = 200, show.plot = FALSE, mean.p = mp[i],
        sig.lp = slp[j])
      results2[i,j,k,] <- c(tmp$Ntotal.fs, tmp$Ntotal.obs)
    }
  }
}


# Plot these two prediction matrices
op <- par(mfrow = c(1, 2), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

# Plot proportion of species occurring in sampled sites (Fig. 11-4 left)
z1 <- apply(results1/200, c(1,2), mean) # Prop species occurring
image(x=mpsi, y=slpsi, z=z1, col = mapPalette(100), axes = TRUE,
    xlab = "Occupancy probability (psi)",
    ylab = "Among-species variability in psi")
contour(x=mpsi, y=slpsi, z=z1, add = TRUE, col = "blue", labcex = 1.5, lwd = 1.5)

# Plot proportion of species detected in sampled sites (Fig. 11-4 right)
z2 <- apply(results2[,,,2] / results2[,,,1], c(1,2), mean)
image(x=mp, y=slp, z=z2, col = mapPalette(100), axes = TRUE,
    xlab = "Detection probability (p)",
    ylab = "Among-species variability in p")
contour(x=mp, y=slp, z=z2, add = TRUE, col = "blue", labcex = 1.5, lwd = 1.5)
par(op)
