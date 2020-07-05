#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# Approximate execution time for this code: 12 mins
# Run time with the full number of iterations: 5 hrs

library(AHMbook)
library(unmarked)

# 6.6 Study design, and bias and precision of the N-mixture estimator
# ===================================================================


# Define simulation settings and arrays for sim results
# simreps <- 1000                    # Simulate and analyse 1000 data sets
simreps <- 50                      # ~~~ smaller number for testing, ca. 15 mins
nsites <- c(20, 120, 250)          # Levels for nsites factor
nreps <- c(2, 5, 10)               # Levels of nrep factor
estimates <- array(NA, dim = c(2, simreps, 3, 3))

# Fill p with random numbers between 0.01 and 0.99
p <- array(runif(n=simreps*3*3, 0.01, 0.99), dim = c(simreps, 3, 3))

# Launch simulation (takes about 6.3 hours with simreps = 1000)
for(s in 1:3){                     # Loop over levels of nsites factor
  for(r in 1:3){                   # Loop over levels of nreps factor
    for(i in 1:simreps){           # Simulate and analyse 1000 data sets
      cat("*** Simrep number", i, "***\n")
      data <- simNmix(nsite=nsites[s], nvisit=nreps[r], mean.lam = 5,
          mean.p=p[i,s,r], show.plot = FALSE)        # Generate data set
      umf <- unmarkedFramePCount(y = data$C) # Bundle data for unmarked
      fm <- pcount(~1 ~1, umf)               # Fit model
      estimates[,i,s,r] <- coef(fm)          # Save estimates
    }
  }
}

# Visualisation
op <- par(mfrow = c(3,3), mar = c(4.5,4.5,2,2), cex.lab = 1.5, cex.axis = 1.3)
for(s in 1:3){                     # Loop over nsites
  for(r in 1:3){                   # Loop over nreps
    plot(p[,s,r], exp(estimates[1,,s,r]), xlab = "Detection probability",
      ylab = "lambda_hat", main = "", ylim = c(0, 75), frame = FALSE)
    text(0.75, 60, paste("M = ", nsites[s], ", J = ", nreps[r], sep = ""),
      cex = 1.5)
    abline(h = 5, col = "red", lwd = 2)
    lines(smooth.spline(exp(estimates[1,,s,r])~p[,s,r]), col="blue", lwd=2)
  }
}
par(op)
