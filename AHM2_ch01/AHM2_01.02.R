#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
# Chapter 1 : RELATIVE ABUNDANCE MODELS FOR POPULATION DYNAMICS
# =============================================================
# Code from proofs dated 2020-08-18

library(AHMbook)

# 1.2 Risks in the naive interpretation of relative abundance
# ===========================================================

# Choose constants
set.seed(1)     # 'freeze' RNGs
M <- 250        # Number of sites
T <- 20         # Number of years
lambda <- 100   # Expected abundance at t = 1
gamma <- 1.02   # Population growth rate
p <- 0.6        # Detection probability

# Create array for true abundance and for counts
N <- C <- array(NA, dim = c(M, T))

# Simulate initial conditions of system: true state at t=1
N[,1] <- rpois(M, lambda)
table(N[,1]) # Summarize

# Simulate later true states
for(t in 2:T){
  N[,t] <- rpois(M, N[,t-1] * gamma)
}

# Simulate binomial observation process and generate actual counts
for(t in 1:T){
  C[,t] <- rbinom(M, N[,t], p)
}

op <- par(mfrow = c(1, 3))  # not shown
hist(N, breaks = 100, main = 'N', col = 'grey')
hist(C, breaks = 100, main = 'C', col = 'grey')
plot(N, C, xlab = 'True N', ylab = 'Observed C', frame = FALSE)
abline(0,1)
lm(c(C) ~ c(N)) # Check slope corresponds to p .... OK !
par(op)

op <- par(mfrow = c(2, 2))  # not shown
ylim <- range(N)
matplot(t(N), type = 'l', lty = 1, main = 'Trajectories of true N', frame = FALSE,
    ylim = ylim)
matplot(t(C), type = 'l', lty = 1, main = 'Trajectories of observed C', frame = FALSE,
    ylim = ylim)
plot(table(N[,1]), main = 'Initial N', frame = FALSE)
plot(table(N[,T]), main = 'Final N', frame = FALSE)
par(op)

library(AHMbook)
str(tmp <- simNpC())
str(tmp <- simNpC(T = 20, expN = c(100, 75), dp = c(0.5, 0.5))) # Explicit defaults
str(tmp <- simNpC(T = 20, expN = c(100, 75), dp = c(1, 1)))     # p = 1

# Simulate data for Fig. 1.1
set.seed(1)
# Declining population
str(tmp1 <- simNpC(T = 20, expN = c(100, 75), dp = c(1, 0.5)))    # p declining
str(tmp2 <- simNpC(T = 20, expN = c(100, 75), dp = c(0.5, 0.5)))  # p stable
str(tmp3 <- simNpC(T = 20, expN = c(100, 75), dp = c(0.5, 1)))    # p increasing
# Stable population
str(tmp4 <- simNpC(T = 20, expN = c(75, 75), dp = c(1, 0.5)))     # p declining
str(tmp5 <- simNpC(T = 20, expN = c(75, 75), dp = c(0.5, 0.5)))   # p stable
str(tmp6 <- simNpC(T = 20, expN = c(75, 75), dp = c(0.5, 1)))     # p increasing
# Increasing population
str(tmp7 <- simNpC(T = 20, expN = c(75, 100), dp = c(1, 0.5)))    # p declining
str(tmp8 <- simNpC(T = 20, expN = c(75, 100), dp = c(0.5, 0.5)))  # p stable
str(tmp9 <- simNpC(T = 20, expN = c(75, 100), dp = c(0.5, 1)))    # p increasing

