#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 4. Introduction to data simulation
# =========================================================================


# 4.2 Generation of a typical point count data set
# ================================================

# 4.2.1 Initial steps: sample size and covariate values
# ------------------------------------------------------------------------
M <- 267         # Number of spatial replicates (sites)
J <- 3           # Number of temporal replicates (counts)

set.seed(24)     # Can choose seed of your choice

elev <- runif(n = M, -1, 1)             # Scaled elevation of a site
forest <- runif(n = M, -1, 1)           # Scaled forest cover at each site
wind <- array(runif(n = M*J, -1, 1), dim = c(M, J)) # Scaled wind speed

# 4.2.2 Simulating the ecological process and its outcome: great tit abundance
# ------------------------------------------------------------------------
mean.lambda <- 2               # Mean expected abundance of great tits
beta0 <- log(mean.lambda)      # Same on log scale (= log-scale intercept)
beta1 <- -2                    # Effect (slope) of elevation
beta2 <- 2                     # Effect (slope) of forest cover
beta3 <- 1                     # Interaction effect (slope) of elev and forest

log.lambda <- beta0 + beta1 * elev + beta2 * forest + beta3 * elev * forest
lambda <- exp(log.lambda)      # Inverse link transformation

op <- par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
curve(exp(beta0 + beta1*x), -1, 1, col = "red", frame.plot = FALSE,
    ylim = c(0, 18), xlab = "Elevation", ylab = "lambda", lwd = 2)
text(-0.9, 17, "A", cex = 1.5)
plot(elev, lambda, frame.plot = FALSE, ylim = c(0, 38), xlab = "Elevation", ylab = "")
text(-0.9, 36, "B", cex = 1.5)
curve(exp(beta0 + beta2*x), -1, 1, col = "red", frame.plot = FALSE,
    ylim = c(0, 18), xlab = "Forest cover", ylab = "lambda", lwd = 2)
text(-0.9, 17, "C", cex = 1.5)
plot(forest, lambda, frame.plot = FALSE, ylim = c(0, 38), xlab = "Forest cover", ylab = "")
text(-0.9, 36, "D", cex = 1.5)
par(op)

# Compute expected abundance for a grid of elevation and forest cover
cov1 <- seq(-1, 1,,100)                       # Values for elevation
cov2 <- seq(-1,1,,100)                        # Values for forest cover
lambda.matrix <- array(NA, dim = c(100, 100)) # Prediction matrix,
# for every combination of values of elevation and forest cover
for(i in 1:100){
  for(j in 1:100){
    lambda.matrix[i, j] <- exp(beta0 + beta1 * cov1[i] + beta2 * cov2[j] +
        beta3 * cov1[i] * cov2[j])
  }
}

op <- par(mfrow = c(1, 2), mar = c(5,4,3,2), cex.main = 1.6)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x = cov1, y = cov2, z = lambda.matrix, col = mapPalette(100),
    xlab = "Elevation", ylab = "Forest cover", cex.lab = 1.2)
contour(x = cov1, y = cov2, z = lambda.matrix, add = TRUE, lwd = 1)
title(main = "A")
matpoints(elev, forest, pch="+", cex=0.8)
par(op)

N <- rpois(n = M, lambda = lambda)   # Realised abundance
sum(N)                               # Total population size at M sites
table(N)                             # Frequency distribution of tit abundance

# 4.2.3 Simulating the observation process and its outcome: point counts of great tits
# ------------------------------------------------------------------------
mean.detection <- 0.3            # Mean expected detection
alpha0 <- qlogis(mean.detection) # same on logit scale (intercept)
alpha1 <- 1                      # Effect (slope) of elevation
alpha2 <- -3                     # Effect (slope) of wind speed
alpha3 <- 0                      # Interaction effect (slope) of elev and wind

logit.p <- alpha0 + alpha1 * elev + alpha2 * wind + alpha3 * elev * wind
p <- plogis(logit.p)             # Inverse link transform
mean(p)                          # average per-individual p is 0.38

op <- par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", frame.plot = FALSE,
    ylim = c(0, 1.1), xlab = "Elevation", ylab = "p", lwd = 2)
text(-0.9, 1.05, "A", cex = 1.5)
matplot(elev, p, pch = "*", frame.plot = FALSE, ylim = c(0, 1.1),
    xlab = "Elevation", ylab = "")
text(-0.9, 1.05, "B", cex = 1.5)
curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", frame.plot = FALSE,
    ylim = c(0, 1.1), xlab = "Wind speed", ylab = "p", lwd = 2)
text(-0.9, 1.05, "C", cex = 1.5)
matplot(wind, p, pch = "*", frame.plot = FALSE, ylim = c(0, 1.1),
    xlab = "Wind speed", ylab = "p")
text(-0.9, 1.05, "D", cex = 1.5)
par(op)

# Compute expected detection probability for a grid of elevation and wind speed
cov1 <- seq(-1, 1,,100)                  # Values of elevation
cov2 <- seq(-1,1,,100)                   # Values of wind speed
p.matrix <- array(NA, dim = c(100, 100)) # Prediction matrix which combines every value in cov 1 with every other in cov2
for(i in 1:100){
  for(j in 1:100){
    p.matrix[i, j] <- plogis(alpha0 + alpha1 * cov1[i] + alpha2 * cov2[j] + alpha3 * cov1[i] * cov2[j])
  }
}
image(x = cov1, y = cov2, z = p.matrix, col = mapPalette(100),
    xlab = "Elevation", ylab = "Wind speed", cex.lab = 1.2)
contour(x = cov1, y = cov2, z = p.matrix, add = TRUE, lwd = 1)
title(main = "B")
matpoints(elev, wind, pch="+", cex=0.7, col = "black")

C <- matrix(NA, nrow = M, ncol = J)      # Prepare array for counts
for (i in 1:J){                          # Generate counts
  C[,i] <- rbinom(n = M, size = N, prob = p[,i])
}

head(cbind("True N" = N, "1st count" = C[,1], "2nd count" = C[,2],
    "3rd count" = C[,3]), 10)                    # First 10 rows (= sites)

table(C)

op <- par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
matplot(elev, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38),
    xlab = "Elevation", ylab = "Count (C)")
text(-0.9, 36, "A", cex = 1.5)
matplot(forest, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38),
    xlab = "Forest cover", ylab = "Count (C)")
text(-0.9, 36, "B", cex = 1.5)
matplot(wind, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38),
    xlab = "Wind speed", ylab = "Count (C)")
text(-0.9, 36, "C", cex = 1.5)
hist(C, breaks = 50, col = "grey", ylim = c(0, 460), main = "", xlab = "Count (C)")
text(3, 450, "D", cex = 1.5)
par(op)

sum(N)                   # True total abundance (all sites)
sum(apply(C, 1, max))    # 'Observed' total abundance (all sites)

sum(N>0)                 # True number of occupied sites
sum(apply(C, 1, max)>0)  # 'Observed' number of occupied sites
