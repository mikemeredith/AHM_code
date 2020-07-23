#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10 : INTEGRATED MODELS FOR MULTIPLE TYPES OF DATA
# =========================================================
# Code from proofs dated 2020-07-23


library(AHMbook)

# 10.2 A simulation game to improve your intuition about point, abundance,
#      and occurrence patterns
# ========================================================================

# Function call with explicit default arguments(requires AHMbook)
str(dat <- simPPe(lscape.size = 150, buffer.width = 25, variance.X = 1,
    theta.X = 10, M = 250, beta = 1, quads.along.side = 6))

# Produce Fig. 10.2
set.seed(117, sample.kind="Rounding")
str(dat <- simPPe(lscape.size = 200, buffer.width = 25, variance.X = 1,
    theta.X = 70, M = 200, beta = 1, quads.along.side = 6))

# Smaller study area, fewer individuals (M)
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 1,
    theta.X = 10, M = 50, beta = 1, quads.along.side = 6))

# Stronger habitat heterogeneity (variance.X): more aggregation
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 10,
    theta.X = 10, M = 50, beta = 1, quads.along.side = 6))

# Longer habitat gradient (theta.X)
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 1,
    theta.X = 250, M = 250, beta = 1, quads.along.side = 6))

# No habitat variability (variance.X): homogeneous point process
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 0,
    theta.X = 10, M = 100, beta = 1, quads.along.side = 6))

# No habitat preference (beta): homogeneous point process
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 1,
    theta.X = 10, M = 100, beta = 0, quads.along.side = 6))

# Habitat heterogeneity at very small scale (theta.X) -> (almost)
# homogeneous point process (in spite of strong habitat preference)
str(dat <- simPPe(lscape.size = 1000, buffer.width = 20, variance.X = 1,
    theta.X = 0.001, M = 250, beta = 1, quads.along.side = 6))

str(simPPe(M = 1)) # This often produces no point at all
str(simPPe(M = 10))
str(simPPe(M = 100))
str(simPPe(M = 1000))

str(simPPe(M = 20, quads.along.side = 50)) # Lots of small sites
str(simPPe(M = 20, quads.along.side = 10))
str(simPPe(M = 20, quads.along.side = 5))
str(simPPe(M = 20, quads.along.side = 1)) # study area is one single site
