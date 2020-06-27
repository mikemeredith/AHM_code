#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 1 Distribution, abundance and species richness in ecology
# =========================================================================

library(AHMbook)

# 1.1 Point processes, distribution, abundance and species richness
# =================================================================

sim.fn(quad.size = 10, cell.size = 1, intensity = 1)


set.seed(82)
tmp <- sim.fn(quad.size = 16, cell.size = 2, intensity = 0.5)


# Effect of grain size of study on abundance and occupancy (intensity constant)
tmp <- sim.fn(quad.size = 10, cell.size = 1, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 2, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 5, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 10, intensity = 0.5)


# Effect of intensity of point pattern (intensity) on abundance and occupancy
tmp <- sim.fn(intensity = 0.1) # choose default quad.size = 10, cell.size = 1
tmp <- sim.fn(intensity = 1)
tmp <- sim.fn(intensity = 5)
tmp <- sim.fn(intensity = 10)


simrep <- 100                 # Run 100 simulation reps ## See errata
grain <- c(0.1,0.2,0.25,0.5,1,2) # values will be fed into 'cell.size' argument
int <- seq(0.1, 3,,6)         # values will be fed into 'lambda' argument
n.levels <- length(grain)     # number of factor levels in simulation
results <- array(NA, dim = c(n.levels, n.levels, 2, simrep)) # 4-D array !
for(i in 1:n.levels){         # Loop over levels of factor grain
  for(j in 1:n.levels){       # Loop over levels of factor intensity
    for(k in 1:simrep){
      cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
      tmp <- sim.fn(cell.size = grain[i], intensity = int[j], show.plot = FALSE)
      results[i,j,1:2,k] <- c(mean(tmp$N), tmp$psi)
    }
  }
}


# Plot these two prediction matrices (NOT IN BOOK)
op <- par(mfrow = c(2, 2), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
# Plot mean abundance in sampled quadrats
z1 <- apply(results[,,1,], c(1,2), mean)   # mean abundance
image(x=grain, y=int, z=z1, col = mapPalette(100), axes = TRUE,
    xlab = "Grain size (cell.size)", ylab = "Intensity of PPP")
contour(x=grain, y=int, z=z1, add = TRUE, col = "blue", labcex = 1.5, lwd = 1.5)
# Plot mean occupancy in sampled quadrats
z2 <- apply(results[,,2,], c(1,2), mean)   # mean occupancy
image(x=grain, y=int, z=z2, col = mapPalette(100), axes = TRUE,
    xlab = "Grain size (cell.size)", ylab = "Intensity of PPP")
contour(x=grain, y=int, z=z2, add = TRUE, col = "blue", labcex = 1.5, lwd = 1.5)
# Plot relationship between occupancy and abundance for whole range of abundance
plot(results[,,1,], results[,,2,], xlab = "Mean abundance", ylab = "Occupancy", frame = FALSE)
lines(smooth.spline(results[,,2,] ~ results[,,1,], df = 4), lwd = 3, col = "blue")
#abline(0, 1, lwd = 3)
# ... and only for very small abundance
keep <- results[,,1,] < 0.25
plot(results[,,1,][keep], results[,,2,][keep], xlab = "Mean abundance", ylab = "Occupancy", frame = FALSE)
abline(0, 1, lwd = 3)
lines(smooth.spline(results[,,2,][keep] ~ results[,,1,][keep], df = 4), lwd = 3, col = "blue")
par(op)
