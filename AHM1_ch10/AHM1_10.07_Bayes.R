#   Applied hierarchical modeling in ecology, vol 1 (2016)
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

# Approximate execution time for this code: < 10 mins
# With the full number of iterations, 4.6 hrs

library(AHMbook)
library(jagsUI)

# 10.7 Study design, and bias and precision of site-occupancy estimators
# ======================================================================

# Supplementary code to run the simulation using JAGS
#     (compiled by Mike Meredith)

# Create the model used in Section 10.3
sink("model.txt")
cat("
model {
  # Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  # Likelihood
  for (i in 1:M) {    # Loop over sites
    z[i] ~ dbern(psi)         # State model
    for (j in 1:J) { # Loop over replicate surveys
      y[i,j] ~ dbern(z[i]*p)  # Observation model (only JAGS !)
      # y[i,j] ~ dbern(mu[i])  # For WinBUGS define 'straw man'
    }
    # mu[i] <- z[i]*p          # Only WinBUGS
  }
}
",fill = TRUE)
sink()

# Do simulation with 1000 reps
simreps <- 10   # ~~~ for testing, 3 mins
# simreps <- 1000  # 4.6 hrs

# Define arrays to hold the results
p <- array(dim = c(simreps,3,3))
estimatesB <- array(dim = c(2,simreps,3,3))

# Choose number and levels of simulation factors
nsites <- c(20, 120, 250)     # number of sites
nsurveys <- c(2, 5, 10)       # number of repeat surveys

params <- c("psi", "p")
ni <- 5000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3

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
         sd.lp = 0, b = 0, show.plot = FALSE)
      # Fit model
      win.data <- list(y = data$y, M = nrow(data$y), J = ncol(data$y))
      zst <- apply(data$y, 1, max)
      inits <- function(){list(z = zst)}
      out <- jags(win.data, inits, params, "model.txt",
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
          DIC=FALSE, parallel=TRUE, verbose=FALSE)
      # Save results (p and posterior means)
      p[i,j,k] <- data$mean.det
      estimatesB[,i,j,k] <- unlist(out$mean)
    }
  }
}
)
save(estimatesB, file="AHM1_10.07_Bayes_estimatesB.RData")

# Plot results - figure 10.8
op <- par(mfrow = c(3,3), mar = c(4,5,3,1), cex.main = 1.2)
for(j in 1:3){
   for(k in 1:3){
      lab <- paste(nsites[j],"sites,", nsurveys[k],"surveys")
      plot(p[,j,k], estimatesB[1,,j,k], xlab = "Detection prob.",
         ylab = "Occupancy prob.", main = lab, ylim = c(0,1))
      abline(h = 0.5, col = "red", lwd = 2)
      lines(smooth.spline(estimatesB[1,,j,k] ~ p[,j,k], df = 5),
         col = "blue", lwd = 2)
   }
}
par(op)
