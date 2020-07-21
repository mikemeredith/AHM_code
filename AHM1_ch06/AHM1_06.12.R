#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# Approximate execution time for this code: 20 mins

library(AHMbook)

# 6.12 Time-for-space substitution
# ================================

simpleNmix(nyear = 12, nrep = 4, beta0 = 2, beta1 = 0.1, alpha0 = 0.5,
    alpha1 = -0.1, alpha2 = 1)

# Define function to simulate such data
simpleNmix <- function(nyear = 12, nrep = 4, beta0 = 2, beta1 = 0.1,
    alpha0 = 0.5, alpha1 = -0.1, alpha2 = 1){
  # Simple function simulates data under binomial N-mixture model where you have
  # a single site that is survyed over 'nyear' primary sampling periods
  # ('seasons', 'years'), within which there are 'nrep' secondary samples each
  # alpha0, alpha1 are the logit-linear coefficients of detection (p) on Time
  #    and on a survey-specific covariate such as temperature (temp).
  # beta0 and beta1 are the log-linear coefficients of expected abundance
  #   (lambda) on Time.

  Time <- 1:nyear
  temp <- matrix(runif(nyear*nrep, -2, 2), ncol = nrep)
  N <- rpois(n = nyear, lambda = exp(beta0 + beta1 * Time))
  C <- array(NA, dim = c(nyear, nrep))
  p <- plogis(alpha0 + alpha1*Time + alpha2*temp)
  for(j in 1:nrep){
     C[,j] <- rbinom(n = nyear, size = N, prob =p[,j])
  }
  op <- par(mfrow = c(3, 2))
  curve(exp(beta0 + beta1 * x), 1, nyear, main = "Expected abundance (lambda) over time",
      frame = FALSE, lwd = 2, ylab = "lambda", xlab = "Time")
  plot(Time, N, main = "Realized abundance (N) over time", frame = FALSE)
  curve(plogis(alpha0 +alpha1 * x), 1, nyear, main = "p over time", frame = FALSE,
      lwd = 2, xlab = "Time", ylab = "p (at averate temp)")
  matplot(Time, C, main = "Counts (C) over time", frame = FALSE)
  curve(plogis(alpha0 + alpha2 * x), -2, 2, main = "p vs. Temperature",
      frame = FALSE, lwd = 2, xlab = "Temperature", ylab = "p (at start of study)")
  matplot(temp, C, main = "Counts (C) over time", frame = FALSE)
  par(op)

  return(list(nyear=nyear, nrep=nrep, beta0=beta0, beta1=beta1, alpha0=alpha0,
      alpha1=alpha1, alpha2=alpha2, N=N, C=C, Time=Time, temp = temp, p = p))
}


library(unmarked)
simrep <- 2500                  # Number of simreps # ca 25 mins
results <- array(NA, dim = c(simrep, 8)) # Array for results
for(i in 1:simrep){
  cat("Simrep", i, "\n")
  data <- simpleNmix(nyear = 12, nrep = 4) # Simulate a data set
  umf <- unmarkedFramePCount(y = data$C,
      siteCovs = data.frame(Time = data$Time), obsCov = list(temp = data$temp))
  fm1 <- pcount(~Time+temp ~Time, data = umf)
  fm2 <- glm(c(data$C)~rep(data$Time,data$nrep)+c(data$temp),family='poisson')
  results[i, 1:5] <- coef(fm1)
  results[i, 6:8] <- coef(fm2)
}
colnames(results) <- c(names(coef(fm1)), names(coef(fm2)))


op <- par(mfrow = c(2,2), mar = c(4,4,3,2))
hist(results[,1], breaks = 100, col = "grey", main = "lambda(Int) Nmix", xlim = c(0, 4))
abline(v = data$beta0, col = "red", lwd = 3)
abline(v = mean(results[,1]), col = "blue", lwd = 3)

hist(results[,2], breaks = 100, col = "grey", main = "lambda(Slope Time) Nmix",
    xlim = c(-0.1, 0.3))
abline(v = data$beta1, col = "red", lwd = 3)
abline(v = mean(results[,2]), col = "blue", lwd = 3)

hist(results[,6], breaks = 100, col = "grey", main = "lambda(Int) GLM",
    xlim = c(0, 4))
abline(v = data$beta0, col = "red", lwd = 3)
abline(v = mean(results[,6]), col = "blue", lwd = 3)

hist(results[,7], breaks = 100, col = "grey", main = "lambda(Slope Time) GLM",
    xlim = c(-0.1, 0.3))
abline(v = data$beta1, col = "red", lwd = 3)
abline(v = mean(results[,7]), col = "blue", lwd = 3)
par(op)
