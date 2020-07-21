#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 4. Introduction to data simulation
# =========================================================================

# 4.3 Packaging everything in a function
# ======================================
# Function definition with set of default values
data.fn <- function(M = 267, J = 3, mean.lambda = 2, beta1 = -2,
    beta2 = 2, beta3 = 1, mean.detection = 0.3, alpha1 = 1, alpha2 = -3,
    alpha3 = 0, show.plot = TRUE){
#
# Function to simulate point counts replicated at M sites during J occasions.
# Population closure is assumed for each site.
# Expected abundance may be affected by elevation (elev),
# forest cover (forest) and their interaction.
# Expected detection probability may be affected by elevation,
# wind speed (wind) and their interaction.
# Function arguments:
#     M: Number of spatial replicates (sites)
#     J: Number of temporal replicates (occasions)
#     mean.lambda: Mean abundance at value 0 of abundance covariates
#     beta1: Main effect of elevation on abundance
#     beta2: Main effect of forest cover on abundance
#     beta3: Interaction effect on abundance of elevation and forest cover
#     mean.detection: Mean detection prob. at value 0 of detection covariates
#     alpha1: Main effect of elevation on detection probability
#     alpha2: Main effect of wind speed on detection probability
#     alpha3: Interaction effect on detection of elevation and wind speed
#     show.plot: if TRUE, plots of the data will be displayed;
#        set to FALSE if you are running simulations.

# Create covariates
elev <- runif(n = M, -1, 1)                         # Scaled elevation
forest <- runif(n = M, -1, 1)                       # Scaled forest cover
wind <- array(runif(n = M*J, -1, 1), dim = c(M, J)) # Scaled wind speed

# Model for abundance
beta0 <- log(mean.lambda)               # Mean abundance on link scale
lambda <- exp(beta0 + beta1*elev + beta2*forest + beta3*elev*forest)
N <- rpois(n = M, lambda = lambda)      # Realised abundance
Ntotal <- sum(N)                        # Total abundance (all sites)
psi.true <- mean(N>0)                   # True occupancy in sample

# Plots
if(show.plot){
  op <- par(mfrow = c(2, 2), cex.main = 1)
  # devAskNewPage(ask = TRUE)
  devAskNewPage(ask = dev.interactive(orNone=TRUE))  # ~~~~~ for autotesting
  curve(exp(beta0 + beta1*x), -1, 1, col = "red",
      main = "Relationship lambda-elevation \nat average forest cover",
      frame.plot = FALSE, xlab = "Scaled elevation")
  plot(elev, lambda, xlab = "Scaled elevation",
      main = "Relationship lambda-elevation \nat observed forest cover",
      frame.plot = FALSE)
  curve(exp(beta0 + beta2*x), -1, 1, col = "red",
      main = "Relationship lambda-forest \ncover at average elevation",
      xlab = "Scaled forest cover", frame.plot = FALSE)
  plot(forest, lambda, xlab = "Scaled forest cover",
      main = "Relationship lambda-forest cover \nat observed elevation",
      frame.plot = FALSE)
  par(op)
}

# Model for observations
alpha0 <- qlogis(mean.detection)        # mean detection on link scale
p <- plogis(alpha0 + alpha1*elev + alpha2*wind + alpha3*elev*wind)
C <- matrix(NA, nrow = M, ncol = J)     # Prepare matrix for counts
for (i in 1:J){                         # Generate counts by survey
  C[,i] <- rbinom(n = M, size = N, prob = p[,i])
}
summaxC <- sum(apply(C,1,max))          # Sum of max counts (all sites)
psi.obs <- mean(apply(C,1,max)>0)       # Observed occupancy in sample

# More plots
if(show.plot){
  op <- par(mfrow = c(2, 2))
  curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red",
      main = "Relationship p-elevation \nat average wind speed",
      xlab = "Scaled elevation", frame.plot = FALSE)
  matplot(elev, p, xlab = "Scaled elevation",
      main = "Relationship p-elevation\n at observed wind speed",
      pch = "*", frame.plot = FALSE)
  curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red",
      main = "Relationship p-wind speed \n at average elevation",
      xlab = "Scaled wind speed", frame.plot = FALSE)
  matplot(wind, p, xlab = "Scaled wind speed",
      main = "Relationship p-wind speed \nat observed elevation",
      pch = "*", frame.plot = FALSE)

  matplot(elev, C, xlab = "Scaled elevation",
      main = "Relationship counts and elevation",
      pch = "*", frame.plot = FALSE)
  matplot(forest, C, xlab = "Scaled forest cover",
      main = "Relationship counts and forest cover",
      pch = "*", frame.plot = FALSE)
  matplot(wind, C, xlab = "Scaled wind speed",
      main = "Relationship counts and wind speed", pch = "*", frame.plot = FALSE)
  desc <- paste('Counts at', M, 'sites during', J, 'surveys')
  hist(C, main = desc, breaks = 50, col = "grey")
  par(op)
}

# Output
return(list(M = M, J = J, mean.lambda = mean.lambda, beta0 = beta0, beta1 = beta1,
    beta2 = beta2, beta3 = beta3, mean.detection = mean.detection, alpha0 = alpha0,
    alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, elev = elev, forest = forest,
    wind = wind, lambda = lambda, N = N, p = p, C = C, Ntotal = Ntotal,
    psi.true = psi.true, summaxC = summaxC, psi.obs = psi.obs))
}


data.fn()                  # Execute function with default arguments
data.fn(show.plot = FALSE) # same, without plots
data.fn(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1,
    mean.detection = 0.3, alpha1 = 1, alpha2 = -3, alpha3 = 0) # Explicit defaults
data <- data.fn()          # Assign results to an object called 'data'

simrep <- 10000
NTOTAL <- SUMMAXC <- numeric(simrep)
for(i in 1:simrep){
   data <- data.fn(show.plot = FALSE)
   NTOTAL[i] <- data$Ntotal
   SUMMAXC[i] <- data$summaxC
}
plot(sort(NTOTAL), ylim = c(min(SUMMAXC), max(NTOTAL)), ylab = "",
    xlab = "Simulation", col = "red", frame = FALSE)
points(SUMMAXC[order(NTOTAL)], col = "blue")

data.fn(J = 2)              # Only 2 surveys
data.fn(J = 1)              # No temporal replicate
data.fn(M = 1, J = 100)     # No spatial replicates, but 100 counts
data.fn(alpha3 = 1)         # With interaction elev-wind on p ## see errata
data.fn(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1,
    mean.detection = 1)     # No obs. process (i.e., p = 1, perfect detection)
data.fn(mean.lambda = 50)   # Really common species
data.fn(mean.lambda = 0.05) # Really rare species

set.seed(24)
data <- data.fn()    # Default arguments
str(data)            # Look at the object

attach(data)         # Make objects inside of 'data' accessible directly

detach(data)         # Make clean up
