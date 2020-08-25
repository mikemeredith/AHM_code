#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 10 : INTEGRATED MODELS FOR MULTIPLE TYPES OF DATA
# =========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 12 mins
# Run time with the full number of iterations: 90 mins

library(AHMbook)
library(jagsUI)

# 10.7 Example 5: A simple integrated population model
# ====================================================

library(AHMbook)
data("greenWoodpecker") # See Chapter 2 for data formatting etc
str(C <- array(as.matrix(greenWoodpecker[,7:48]), dim=c(267, 3, 14)))

# 10.7.1 Fitting the DM model to the count data alone
# ---------------------------------------------------

# Bundle and summarize data
str(bdata <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3]))

# Specify model in BUGS language
cat(file = "DM.txt","
model {

  # Priors
  lambda ~ dunif(0, 5)            # Initial abundance
  phi ~ dunif(0, 1)               # Apparent survival
  gamma ~ dunif(0, 3)             # Recruitment rate (absolute)
  p ~ dunif(0, 1)                 # Detection probability

  # Likelihood for the counts (= DM model)
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t])
      R[i,t+1] ~ dpois(gamma)     # 'absolute' recruitment = 'constant'
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    }

    # Observation process for the counts
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,j,t] ~ dbin(p, N[i,t])
      }
    }
  }
}
")

# Initial values function
Rst <- apply(C, c(1,3), max, na.rm = TRUE)
Rst[Rst == '-Inf'] <- 1
Rst[,1] <- NA
N1 <- apply(apply(C, c(1,3), max, na.rm = TRUE),1,max, na.rm = TRUE)
N1[N1 == '-Inf'] <- 2
Nst <- array(NA, dim = dim(Rst))
Nst[,1] <- N1
inits <- function(){list(lambda = runif(1, 1, 2), phi = runif(1),
    gamma = runif(1), p = runif(1), R = Rst+1, N = Nst+2)}

# Parameters monitored
params <- c("lambda", "phi", "gamma", "p")

# MCMC settings
# na <- 1000 ; ni <- 30000 ; nt <- 10 ; nb <- 20000 ; nc <- 3
na <- 1000 ; ni <- 3000 ; nt <- 1 ; nb <- 2000 ; nc <- 3  # ~~~ for testing, 2 mins

# Call JAGS (ART 22 min), check convergence and summarize posteriors
(out1 <- jags(bdata, inits, params, "DM.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE))
#         mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# lambda 0.594 0.069 0.467 0.591 0.735    FALSE 1 1.010   218
# phi    0.821 0.012 0.795 0.821 0.843    FALSE 1 1.030    87
# gamma  0.378 0.017 0.346 0.378 0.411    FALSE 1 1.004   804
# p      0.291 0.009 0.274 0.291 0.309    FALSE 1 1.032    80


# 10.7.2 Fitting the CJS model to the capture-recapture data alone
# ----------------------------------------------------------------

# Simulate five data sets with different sample sizes
set.seed(24)
str(data1 <- simCJS(n.occ = 6, n.marked = 20, phi = 0.58, p = 0.4,
    show.plot = FALSE)) # results in 100 individuals over 6 years
# ~~~ additional code ~~~~~~~~~~~
str(data2 <- simCJS(n.occ = 6, n.marked = 100, phi = 0.58, p = 0.4,
    show.plot = FALSE)) # ... 500 individuals ...
str(data3 <- simCJS(n.occ = 6, n.marked = 200, phi = 0.58, p = 0.4,
    show.plot = FALSE)) # ... 1000 individuals ...
str(data4 <- simCJS(n.occ = 6, n.marked = 1000, phi = 0.58, p = 0.4,
    show.plot = FALSE)) # ...5000 individuals ...
str(data5 <- simCJS(n.occ = 6, n.marked = 2000, phi = 0.58, p = 0.4,
    show.plot = FALSE)) # 10000 individuals over 6 years
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Turn capture-histories into m-arrays
(m_array1 <- ch2marray(data1$ch))
(r1 <- rowSums(m_array1))
# ~~~ additional code ~~~~~~~~~~~
(m_array2 <- ch2marray(data2$ch))
(r2 <- rowSums(m_array2))
(m_array3 <- ch2marray(data3$ch))
(r3 <- rowSums(m_array3))
(m_array4 <- ch2marray(data4$ch))
(r4 <- rowSums(m_array4))
(m_array5 <- ch2marray(data5$ch))
(r5 <- rowSums(m_array5))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make five data bundles and summarize each
str(bdata1 <- list(m_array = m_array1, r = r1, n.occ = data1$n.occ))
# ~~~ additional code ~~~~~~~~~~~
str(bdata2 <- list(m_array = m_array2, r = r2, n.occ = data2$n.occ))
str(bdata3 <- list(m_array = m_array3, r = r3, n.occ = data3$n.occ))
str(bdata4 <- list(m_array = m_array4, r = r4, n.occ = data4$n.occ))
str(bdata5 <- list(m_array = m_array5, r = r5, n.occ = data5$n.occ))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Specify CJS model in the m-array format in the BUGS language
cat(file = "CJS.txt"," # This is the m-array format
model {

  # Priors for modeled parameters
  for(t in 1:(n.occ-1)){
    phi.lk[t] <- phi        # .lk stands for 'in the likelihood'
    p.lk[t] <- p
  }
  phi ~ dunif(0, 1)         # Apparent survival
  p ~ dunif(0, 1)           # Recapture probability

  # Multinomial likelihood for the m-array format of the capture-recapture
  # data (= CJS model)
  for (t in 1:(n.occ-1)){
    m_array[t,1:n.occ] ~ dmulti(pr[t, ], r[t])
  }
  # Define the cell probabilities of the m-array
  for (t in 1:(n.occ-1)){
    q[t] <- 1 - p.lk[t]
    pr[t,t] <- phi.lk[t]*p.lk[t]
    for (j in (t+1):(n.occ-1)){
      pr[t,j] <- prod(phi.lk[t:j])*prod(q[t:(j-1)])*p.lk[j]
    }
    for (j in 1:(t-1)){
      pr[t,j] <- 0
    }
  }
  for (t in 1:(n.occ-1)){
    pr[t,n.occ] <- 1 - sum(pr[t,1:(n.occ-1)])
  }
}
")

# Run all five analyses
inits <- function(){list(phi = runif(1), p = runif(1))}
params <- c("phi", "p")
na <- 1000 ; ni <- 150000 ; nt <- 100 ; nb <- 50000 ; nc <- 3
out2.1 <- jags(bdata1, inits, params, "CJS.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# ~~~ code for the other four cases inserted ~~~~
out2.2 <- jags(bdata2, inits, params, "CJS.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
out2.3 <- jags(bdata3, inits, params, "CJS.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
out2.4 <- jags(bdata4, inits, params, "CJS.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
out2.5 <- jags(bdata5, inits, params, "CJS.txt", n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Table of estimates of phi for 5 sample sizes
# ~~~ extra code ~~~
data.set <- paste(c(100, 500, 1000, 5000, 10000), 'marked birds')
phi.esti <- rbind(out2.1$summary[1,c(1,3,7)], out2.2$summary[1,c(1,3,7)],
    out2.3$summary[1,c(1,3,7)], out2.4$summary[1,c(1,3,7)],
    out2.5$summary[1,c(1,3,7)])
rownames(phi.esti) <- data.set
print(phi.esti, 3) # Remember phi simulated as 0.58
# ~~~~~~~~~~~~~~~~~~~~
#                     mean  2.5% 97.5%
# 100 marked birds   0.623 0.564 0.684
# 1000 marked birds  0.572 0.529 0.617
# 5000 marked birds  0.574 0.554 0.594
# 10000 marked birds 0.590 0.576 0.604


# 10.7.3 Fitting the integrated model
# -----------------------------------

# Bundle and summarize data set for five CJS data sample sizes
str(bdata1 <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3], # Data for DM model
    m_array = m_array1, r = r1, n.occ = data1$n.occ)) # Data for CJS model
# List of 7
# $ C       : int [1:267, 1:3, 1:14] 0 3 0 0 0 0 0 0 0 0 ...
# $ nsites  : int 267
# $ nsurveys: int 3
# $ nyears  : int 14
# $ m_array : num [1:5, 1:6] 459 0 0 0 0 162 584 0 0 0 ...
# $ r       : num [1:5] 2000 2459 2746 2918 3045
# $ n.occ   : num 6
# ~~~~ code for the other four cases inserted ~~~~~~~~~~~
str(bdata2 <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3], m_array = m_array2, r = r2, n.occ = data2$n.occ))
str(bdata3 <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3], m_array = m_array3, r = r3, n.occ = data3$n.occ))
str(bdata4 <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3], m_array = m_array4, r = r4, n.occ = data4$n.occ))
str(bdata5 <- list(C = C, nsites = dim(C)[1], nsurveys = dim(C)[2],
    nyears = dim(C)[3], m_array = m_array5, r = r5, n.occ = data5$n.occ))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Specify model in BUGS language
cat(file = "IPM.txt","
model {

  # Priors for jointly modeled parameters
  lambda ~ dunif(0, 10)           # Initial abundance
  phiJoint ~ dunif(0, 1)          # Apparent survival
  gamma ~ dunif(0, 3)             # Recruitment rate
  pDM ~ dunif(0, 1)               # Detection probability in the count data
  pCJS ~ dunif(0, 1)              # Detection probability in the cap-recap data

  # (1) Likelihood for the counts (= DM model as before)

  # Make one parameter identical in the two submodels
  phiDM <- phiJoint

  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phiDM, N[i,t])
      R[i,t+1] ~ dpois(gamma)             # 'absolute' recruitment = 'constant'
      # R[i,t+1] ~ dpois(N[i,t] * gamma)  # per-capita recr. = 'autoreg'
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    }

    # Observation process for the counts
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,j,t] ~ dbin(pDM, N[i,t])
      }
    }
  }

  # (2) Multinomial likelihood for the m-array format of the capture-
  # recapture data (= CJS model)

  # Constraints for CJS parameters
  for(t in 1:(n.occ-1)){
    phiCJS[t] <- phiJoint         # Survival constant and same as in DM model
    pCJS_t[t] <- pCJS             # Recapture constant
  }

  # Multinomial likelihood for the m-array format of the capture-recapture
  # data (= CJS model)
  for (t in 1:(n.occ-1)){
    m_array[t,1:n.occ] ~ dmulti(pr[t, ], r[t])
  }
  # Define the cell probabilities of the m-array
  for (t in 1:(n.occ-1)){
    q[t] <- 1 - pCJS_t[t]         # Probability of non-recapture
    pr[t,t] <- phiCJS[t]*pCJS_t[t]
    for (j in (t+1):(n.occ-1)){
      pr[t,j] <- prod(phiCJS[t:j])*prod(q[t:(j-1)])*pCJS_t[j]
    }
    for (j in 1:(t-1)){
      pr[t,j] <- 0
    }
  }
  for (t in 1:(n.occ-1)){
    pr[t,n.occ] <- 1-sum(pr[t,1:(n.occ-1)])
  }
}
")

# Initial values function
Rst <- apply(C, c(1,3), max, na.rm = TRUE)
Rst[Rst == '-Inf'] <- 1
Rst[,1] <- NA
N1 <- apply(apply(C, c(1,3), max, na.rm = TRUE),1,max, na.rm = TRUE)
N1[N1 == '-Inf'] <- 2
Nst <- array(NA, dim = dim(Rst))
Nst[,1] <- N1
inits <- function(){list(lambda = runif(1, 1, 2), phiJoint = runif(1),
    gamma = runif(1), pDM = runif(1), pCJS = runif(1), R = Rst+1, N = Nst+2)}

# Parameters monitored
params <- c("lambda", "phiJoint", "gamma", "pDM", "pCJS")

# MCMC settings
# na <- 1000 ; ni <- 30000 ; nt <- 10 ; nb <- 20000 ; nc <- 3
na <- 1000 ; ni <- 3000 ; nt <- 1 ; nb <- 2000 ; nc <- 3  # ~~~ for testing, 2 mins each

# Call JAGS, check convergence and summarize posteriors
out3.1 <- jags(bdata1, inits, params, "IPM.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# ~~~~ code for the other four cases inserted ~~~~~
out3.2 <- jags(bdata2, inits, params, "IPM.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
out3.3 <- jags(bdata3, inits, params, "IPM.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
out3.4 <- jags(bdata4, inits, params, "IPM.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
out3.5 <- jags(bdata5, inits, params, "IPM.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ extra code for figure 10.7 ~~~~~~~~~~~
off <- 0.05
op <- par(mar = c(8,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot((1:6)-off, c(out1$mean$phi, out3.1$mean$phiJoint, out3.2$mean$phiJoint,
    out3.3$mean$phiJoint, out3.4$mean$phiJoint, out3.5$mean$phiJoint),
    ylim = c(0.5, 0.9), xlab = '', ylab = 'Estimates of phi', main = '',
    pch = c(16, rep(15,5)), cex = 2, frame = FALSE, las = 1, axes = FALSE)
segments(1-off, out1$q2.5$phi, 1-off, out1$q97.5$phi)
segments(2-off, out3.1$q2.5$phiJoint, 2-off, out3.1$q97.5$phiJoint)
segments(3-off, out3.2$q2.5$phiJoint, 3-off, out3.2$q97.5$phiJoint)
segments(4-off, out3.3$q2.5$phiJoint, 4-off, out3.3$q97.5$phiJoint)
segments(5-off, out3.4$q2.5$phiJoint, 5-off, out3.4$q97.5$phiJoint)
segments(6-off, out3.5$q2.5$phiJoint, 6-off, out3.5$q97.5$phiJoint)
points((2:6)+off, c(out2.1$mean$phi, out2.2$mean$phi, out2.3$mean$phi,
    out2.4$mean$phi, out2.5$mean$phi), pch = 1, cex = 2)
segments(2+off, out2.1$q2.5$phi, 2+off, out2.1$q97.5$phi)
segments(3+off, out2.2$q2.5$phi, 3+off, out2.2$q97.5$phi)
segments(4+off, out2.3$q2.5$phi, 4+off, out2.3$q97.5$phi)
segments(5+off, out2.4$q2.5$phi, 5+off, out2.4$q97.5$phi)
segments(6+off, out2.5$q2.5$phi, 6+off, out2.5$q97.5$phi)
axis(1, at = 1:6, c('n = 0', 'n = 100', 'n = 500', 'n = 1000', 'n = 5000', 'n = 10000'), las = 2)
abline(h = 0.58, col = 'red')
axis(2)
legend('topright', c('DM','CJS', 'IPM'), pch = c(16,1,15), cex = 1.5, bty = 'n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ and for figure 10.8 ~~~~~~~~~~~~~
off <- 0.05
op <- par(mar = c(8,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot((1:6)-off, c(out1$mean$phi, out3.1$mean$phiJoint, out3.2$mean$phiJoint,
    out3.3$mean$phiJoint, out3.4$mean$phiJoint, out3.5$mean$phiJoint),
    ylim = c(0.5, 0.9), xlab = '', ylab = 'Estimates of phi', main = '',
    pch = c(16, rep(15,5)), cex = 2, frame = FALSE, las = 1, axes = FALSE)
segments(1-off, out1$q2.5$phi, 1-off, out1$q97.5$phi)
segments(2-off, out3.1$q2.5$phiJoint, 2-off, out3.1$q97.5$phiJoint)
segments(3-off, out3.2$q2.5$phiJoint, 3-off, out3.2$q97.5$phiJoint)
segments(4-off, out3.3$q2.5$phiJoint, 4-off, out3.3$q97.5$phiJoint)
segments(5-off, out3.4$q2.5$phiJoint, 5-off, out3.4$q97.5$phiJoint)
segments(6-off, out3.5$q2.5$phiJoint, 6-off, out3.5$q97.5$phiJoint)
points((2:6)+off, c(out2.1$mean$phi, out2.2$mean$phi, out2.3$mean$phi,
    out2.4$mean$phi, out2.5$mean$phi), pch = 1, cex = 2)
segments(2+off, out2.1$q2.5$phi, 2+off, out2.1$q97.5$phi)
segments(3+off, out2.2$q2.5$phi, 3+off, out2.2$q97.5$phi)
segments(4+off, out2.3$q2.5$phi, 4+off, out2.3$q97.5$phi)
segments(5+off, out2.4$q2.5$phi, 5+off, out2.4$q97.5$phi)
segments(6+off, out2.5$q2.5$phi, 6+off, out2.5$q97.5$phi)
axis(1, at = 1:6, c('n = 0', 'n = 100', 'n = 500', 'n = 1000', 'n = 5000',
    'n = 10000'), las = 2)
abline(h = 0.58, col = 'red')
axis(2)
legend('topright', c('DM','CJS', 'IPM'), pch = c(16,1,15), cex = 1.5, bty = 'n')
par(op)

# Plot lambda and gamma estimates as a function of the amount of CJS data
op <- par(mfrow = c(1, 2), mar = c(8,5,4,3))
# lambda
plot(1:6, c(out1$mean$lambda, out3.1$mean$lambda, out3.2$mean$lambda,
    out3.3$mean$lambda, out3.4$mean$lambda, out3.5$mean$lambda),
    ylim = c(0.2, 0.8), xlab = '', ylab = 'lambda', pch = c(16, rep(15,5)),
    cex = 2, frame = FALSE, las = 1, axes = FALSE)
segments(1, out1$q2.5$lambda, 1, out1$q97.5$lambda)
segments(2, out3.1$q2.5$lambda, 2, out3.1$q97.5$lambda)
segments(3, out3.2$q2.5$lambda, 3, out3.2$q97.5$lambda)
segments(4, out3.3$q2.5$lambda, 4, out3.3$q97.5$lambda)
segments(5, out3.4$q2.5$lambda, 5, out3.4$q97.5$lambda)
segments(6, out3.5$q2.5$lambda, 6, out3.5$q97.5$lambda)
axis(1, at = 1:6, c('n = 0', 'n = 100', 'n = 500', 'n = 1000',
    'n = 5000', 'n = 10000'), las = 2)
axis(2)

# gamma
plot(1:6, c(out1$mean$gamma, out3.1$mean$gamma, out3.2$mean$gamma,
    out3.3$mean$gamma, out3.4$mean$gamma, out3.5$mean$gamma),
    ylim = c(0.2, 0.6), xlab = '', ylab = 'gamma', pch = c(16, rep(15,5)),
    cex = 2, frame = FALSE, las = 1, axes = FALSE)
segments(1, out1$q2.5$gamma, 1, out1$q97.5$gamma)
segments(2, out3.1$q2.5$gamma, 2, out3.1$q97.5$gamma)
segments(3, out3.2$q2.5$gamma, 3, out3.2$q97.5$gamma)
segments(4, out3.3$q2.5$gamma, 4, out3.3$q97.5$gamma)
segments(5, out3.4$q2.5$gamma, 5, out3.4$q97.5$gamma)
segments(6, out3.5$q2.5$gamma, 6, out3.5$q97.5$gamma)
axis(1, at = 1:6, c('n = 0', 'n = 100', 'n = 500', 'n = 1000',
    'n = 5000', 'n = 10000'), las = 2)
axis(2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
