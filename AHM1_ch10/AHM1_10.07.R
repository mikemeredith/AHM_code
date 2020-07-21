#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

# Approximate execution time for this code: 10 mins

library(AHMbook)

# 10.7 Study design, and bias and precision of site-occupancy estimators
# ======================================================================


# Do simulation with 1000 reps
simreps <- 1000
library(unmarked)

# Define arrays to hold the results
p <- array(dim = c(simreps,3,3))
estimates <- array(dim = c(2,simreps,3,3))

# Choose number and levels of simulation factors
nsites <- c(20, 120, 250)     # number of sites
nsurveys <- c(2, 5, 10)       # number of repeat surveys

# Start simulation
system.time(                  # time whole thing
for(j in 1:3) {               # Loop j over site factor
  for(k in 1:3) {             # Loop k over survey factor
    for(i in 1:simreps){      # Loop i over simreps
      # Counter
      cat("** nsites", j, "nsurveys", k, "simrep", i, "***\n")
      # Generate a data set: pick p and use p in simOcc()
      det.prob <- runif(1, 0.01, 0.99)
      data <- simOcc(M = nsites[j], J = nsurveys[k], mean.occupancy = 0.5,
         beta1 = 0, beta2 = 0, beta3 = 0, mean.detection = det.prob,
         time.effects = c(0, 0), alpha1 = 0, alpha2 = 0, alpha3 = 0,
         sd.lp = 0, b = 0, show.plot = F)
      # Fit model
      umf <- unmarkedFrameOccu(y = data$y)
      tmp <- occu(~1 ~1, umf, se = FALSE)   # Only get MLEs, not SEs
      # Save results (p and MLEs)
      p[i,j,k] <- data$mean.det
      estimates[,i,j,k] <- coef(tmp)
    }
  }
}
)

# Plot results
op <- par(mfrow = c(3,3), mar = c(4,5,3,1), cex.main = 1.2)
for(j in 1:3){
   for(k in 1:3){
      lab <- paste(nsites[j],"sites,", nsurveys[k],"surveys")
      plot(p[,j,k], plogis(estimates[1,,j,k]), xlab = "Detection prob.",
         ylab = "Occupancy prob.", main = lab, ylim = c(0,1))
         abline(h = 0.5, col = "red", lwd = 2)
         lines(smooth.spline(plogis(estimates[1,,j,k])~ p[,j,k], df = 5),
         col = "blue", lwd = 2)
   }
}
par(op)
