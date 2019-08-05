#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# 6.10 The issue of space, or: what is your effective sample area ?
# ------------------------------------------------------------------------


AHR <- seq(0.001, 100, 0.1)         # Area home range (in km2): 0.1-100
RHR <- sqrt(AHR/pi)                 # Translate to home range radius (in km)
ESA <- (1 + 2*RHR)^2                # Eff. sample area for 1km2 nominal area
par(mfrow = c(1,2), mar = c(5,5,3,2), cex.lab = 1.3, cex.axis = 1.3)
plot(AHR, ESA, xlab = "Home range size (km2)", ylab = "Effective sample area (km2)", type = "l", lwd = 3, frame = F)
abline(h = 1, col = "red", lwd = 3) # Nominal sample area
abline(h = 0, col = "grey", lwd = 3)


GT.HR <- seq(0.003, 0.1,,1000)     # Great tit home ranges in km2
GT.rad <- sqrt(GT.HR/pi)           # Great tit home range radius (in km)
ESA.GT <- (1 + 2*GT.rad)^2         # Effective sample area for Great tit
NTOT <- numeric(length(ESA.GT))    # Adjusted national total population
for(i in 1:length(NTOT)){
   NTOT[i] <- Nhat(fm5ZIP, 0, ESA.GT[i])[3]
}
plot(100*GT.HR, NTOT, xlab = "Great tit home-range (ha)", ylab = "Adjusted national total (pairs)", type = "l", lwd = 3, frame = F)

