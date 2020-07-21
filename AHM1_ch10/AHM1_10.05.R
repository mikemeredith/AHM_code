#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10. Modeling static occurrence and species distributions using
#             site-occupancy models
# =========================================================================

library(AHMbook)

# 10.5 A general data-simulation function for static occupancy models: simOcc()
# =============================================================================


simOcc(M = 267, J = 3, mean.occupancy = 0.6, beta1 = -2, beta2 = 2, beta3 = 1,
    mean.detection = 0.3, time.effects = c(-1, 1), alpha1 = -1, alpha2 = -3,
    alpha3 = 0, sd.lp = 0.5, b = 2, show.plot = TRUE)

simOcc()                  # Execute function with default arguments
simOcc(show.plot = FALSE) #    same, without plots
simOcc(M = 267, J = 3, mean.occupancy = 0.6, beta1 = -2, beta2 = 2, beta3 = 1,
    mean.detection = 0.3, time.effects = c(-1, 1), alpha1 = -1, alpha2 = -3,
    alpha3 = 0, sd.lp = 0.5, b = 2, show.plot = TRUE) # Explicit defaults

# Create a 'fix' data set and look at what we created
set.seed(24)
data <- simOcc()          # Assign results to an object called 'data'
str(data)


# Simplest possible occupancy model, with constant occupancy and detection
tmp <- simOcc(mean.occ=0.6, beta1=0, beta2=0, beta3=0, mean.det=0.3,
    time.effects=c(0, 0), alpha1=0, alpha2=0, alpha3=0, sd.lp=0, b=0)
str(tmp)                 # give overview of results

# psi = 1 (i.e., species occurs at every site)
tmp <- simOcc(mean.occ=1)   ;   str(tmp)

# p = 1 (i.e., species is always detected when it occurs)
tmp <- simOcc(mean.det=1)   ;   str(tmp)

# Other potentially interesting settings include these:
simOcc(J = 2)                 # Only 2 surveys
simOcc(M = 1, J = 100)        # No spatial replicates, but 100 measurements
simOcc(beta3 = 1)             # Including interaction elev-wind on p
simOcc(mean.occ = 0.96)       # A really common species
simOcc(mean.occ = 0.05)       # A really rare species
simOcc(mean.det = 0.96)       # A really easy species
simOcc(mean.det = 0.05)       # A really hard species
simOcc(mean.det = 0)          # The dreaded invisible species
simOcc(alpha1=-2, beta1=2)    # Opposing effects of elev on psi and p
simOcc(J = 10, time.effects = c(-5, 5)) # Huge time effects on p
simOcc(sd.lp = 10)            # Huge (random) site effects on p
simOcc(J = 10, b = 0)         # No behavioural response in p
simOcc(J = 10, b = 2)         # Trap happiness
simOcc(J = 10, b = -2)        # Trap shyness

