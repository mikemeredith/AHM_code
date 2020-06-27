#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 7. Modeling abundance using multinomial N-mixture models
# =========================================================================

library(AHMbook)

# 7.3 Simulating multinomial observations in R
# ============================================

rmultinom(10, 5, c(0.1, 0.2, 0.7))

(rmn <- sample(1:3, 5, replace=TRUE, prob = c(0.1, 0.2, 0.7)))

set.seed(2015)                         # Initialize RNG

# Simulate covariate values and local population size for each point
x <- rnorm(100)
N <- rpois(100, lambda=exp(-1 + 1*x) ) # Intercept and slope equal to 1
table(N)                               # Summarize

# Define detection probabilities (p) for both observers
p1 <- 0.8
p2 <- 0.6

# Construct the multinomial cell probabilities (pi)
cellprobs <- c(p1*p2, p1*(1-p2), (1-p1)*p2, (1-p1)*(1-p2))

# Create a matrix to hold the data
y <- matrix(NA, nrow=100, ncol=4)
dimnames(y) <- list(1:100, c("11", "10", "01", "00"))

# Loop over sites and generate data with function rmultinom()
for(i in 1:100){
   y[i,] <- rmultinom(1, N[i], cellprobs)
}

# Remove 4th column ("not detected") and summarize results
y <- y[,-4]
apply(y, 2, sum)


# Generate specific pseudo-random data set
set.seed(2014)
data <- simNmix(mean.lam = exp(1), beta3.lam = 1, mean.p = plogis(0),
       sigma.p.visit = 1, show.plot=FALSE)
str(data$DH)

# View detection histories for site with max abundance (here, N = 30)
t(data$DH[min(which(data$N == max(dim(data$DH)[3]))),,])


# Get detection history frequencies for each site (for exactly 3 surveys)
dhfreq <- array(NA, dim = c(data$nsite, 7),
  dimnames = list(NULL, c("100", "010", "001", "110", "101", "011", "111")))
for(i in 1:data$nsite){
  dhfreq[i,] <- table(factor(paste(data$DH[i,1,], data$DH[i,2,],
      data$DH[i,3,], sep = ""),
      levels = c("100", "010", "001", "110", "101", "011", "111")))
}

head(dhfreq)      # Data for first 6 sites

# Get occasions with first detection of each individual
f <- apply(data$DH, c(1,3), function(x) min(which(x != 0)))
head(f)   ;   str(f)   ;   table(f)    # Inspect result

# Produce removal counts
y <- array(NA, dim = c(data$nsite, data$nvisit), dimnames = list(NULL,
       as.factor(1:data$nvisit)))
for(i in 1:data$nsite){
   y[i,] <- table(factor(f[i,], levels = as.character(1:data$nvisit)))
}
head(y)      # Data for first 6 sites

set.seed(24)
data <- simNmix(mean.lam = exp(0), beta2.lam = 1, beta3.lam =-1,
   beta4.lam = 0.2, dispersion = 1, mean.p = plogis(0),
   beta3.p = 1, sigma.p.visit = 1, Neg.Bin = TRUE)

