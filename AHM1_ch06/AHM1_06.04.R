#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

library(jagsUI)

# 6.4 A slightly more complex N-mixture model with covariates
# ===========================================================

# Choose sample sizes and prepare obs. data array y
set.seed(1)                   # So we all get same data set
M <- 100                      # Number of sites
J <- 3                        # Number of repeated abundance measurements
C <- matrix(NA, nrow = M, ncol = J) # to contain the observed data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for abundance model and compute lambda
beta0 <- 0                    # Log-scale intercept
beta1 <- 2                    # Log-scale slope for vegHt
lambda <- exp(beta0 + beta1 * vegHt) # Expected abundance
plot(vegHt, lambda, type = "l", lwd = 3)  # Expected abundance

# Draw local abundance and look at data so far
N <- rpois(M, lambda)
points(vegHt, N)              # Add realized abundance to plot
table(N)

# Plot the true system state (Fig. 6–2, left)
op <- par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
plot(vegHt, N, xlab="Vegetation height", ylab="True abundance (N)",
    frame = FALSE, cex = 1.5)
lines(seq(-1,1,,100), exp(beta0 + beta1* seq(-1,1,,100)), lwd=3, col = "red")


# Create a covariate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2                        # Logit-scale intercept
alpha1 <- -3                        # Logit-scale slope for wind
p <- plogis(alpha0 + alpha1 * wind) # Detection probability
#plot(p ~ wind, ylim = c(0,1))       # Look at relationship

# Take J = 3 abundance measurements at each site
for(j in 1:J) {
    C[,j] <- rbinom(M, N, p[,j])
}

# Plot observed data and effect of wind on det. probability (Fig. 6–2, middle)
plot(wind, C/max(C), xlab="Wind", ylab="Scaled counts: C/max(C)",
    frame = FALSE, cex = 1.5)
lines(seq(-1,1,,100), plogis(alpha0 + alpha1*seq(-1,1,,100)), lwd=3, col="red")


# Expected (lambda) and realized abundance (N) and measurements (C)
cbind(lambda=round(lambda,2), N=N, C1=C[,1], C2=C[,2], C3=C[,3])


# Create factors
time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
hab <- c(rep("A", 33), rep("B", 33), rep("C", 34))  # assumes M = 100


# Load unmarked, format data in unmarked data frame and summarize
library(unmarked)
umf <- unmarkedFramePCount(
   y=C,                                            # Counts matrix
   siteCovs= data.frame(vegHt = vegHt, hab = hab), # Site covariates
   obsCovs = list(time = time, wind = wind))       # Observation covs
summary(umf)


# Fit model and extract estimates
# linear model for p follows first tilde, then comes linear model for lambda
summary(fm.Nmix1 <- pcount(~wind ~vegHt, data=umf, control=list(trace=T, REPORT=1)))


fm.Nmix2 <- pcount(~wind ~vegHt, data=umf, mixture="NB",
    control=list(trace=TRUE, REPORT=5))
fm.Nmix3 <- pcount(~wind ~vegHt, data=umf, mixture="ZIP",
    control=list(trace=TRUE, REPORT=5))
cbind(AIC.P=fm.Nmix1@AIC, AIC.NB=fm.Nmix2@AIC, AIC.ZIP=fm.Nmix3@AIC)


# Predictions of lambda for specified values of vegHt, say 0.2 and 2.1
newdat <- data.frame(vegHt=c(0.2, 1.1))
predict(fm.Nmix1, type="state", newdata=newdat, append = T)


# ... or of p for values of wind of -1 to 1
newdat <- data.frame(wind=seq(-1, 1, , 5))
predict(fm.Nmix1, type="det", newdata=newdat, append = T)


# Predict lambda and detection for actual data set
(lambda.hat <- predict(fm.Nmix1, type="state"))     # lambda at every site
(p.hat <- predict(fm.Nmix1, type="det"))            # p during every survey


# Predict lambda and detection as function of covs
newdat <- data.frame(vegHt=seq(-1, 1, 0.01))
pred.lam <- predict(fm.Nmix1, type="state", newdata=newdat)
newdat <- data.frame(wind=seq(-1, 1, 0.1))
pred.det <- predict(fm.Nmix1, type="det", newdata=newdat)


# Fit detection-naive GLM to counts and plot comparison (Fig. 6–2, right)
summary(fm.glm <- glm(c(C) ~ rep(vegHt, 3), family=poisson)) # p-naive  model
matplot(vegHt, C, xlab="Vegetation height", ylab="Counts", frame = FALSE,
    cex = 1.5, pch = 1, col = "black")
lines(seq(-1,1,,100), exp(beta0 + beta1* seq(-1,1,,100)), lwd=3, col = "red")
curve(exp(coef(fm.glm)[1]+coef(fm.glm)[2]*x), -1, 1, type ="l", lwd=3, add=TRUE)
lines(vegHt, predict(fm.Nmix1, type="state")[,1], col = "blue", lwd = 3)
legend(-1, 7, c("Truth", "'Poisson GLM' with p", "Poisson GLM without p"),
    col=c("red", "blue", "black"), lty = 1, lwd=3, cex = 1.2)
par(op)

ranef(fm.Nmix1)

# calculate lambda.hat: exp(a0 + a1*vegHt)
lambda.hat <- predict(fm.Nmix1, type="state")[,1]

# calculate p.hat: plogis(b0 + b1*wind)
p.hat <- matrix(predict(fm.Nmix1, type="det")[,1], ncol=ncol(C), byrow=TRUE)

Ngrid <- 0:(100+max(umf@y, na.rm = TRUE))
posterior <- matrix(NA, nrow=nrow(C), ncol=length(Ngrid))
bup2 <- array(NA, dim = M)

for(i in 1:nrow(C)){      # Loop over sites
  # Compute prior using MLE
  gN <- dpois(Ngrid, lambda.hat[i])
  gN <- gN/sum(gN)

  # Compute likelihood for each possible value of N
  fy <- rep(NA, length(Ngrid))
  for(j in 1:length(Ngrid)){
     fy[j]<- prod(dbinom(C[i,], Ngrid[j], p.hat[i,]))
  }

  # Compute marginal of y. for denominator of Bayes rule
  qy <- sum(fy * gN)

  # Posterior
  posterior[i,] <- fy * gN / qy

  # N can't be less than max(C)
  if(max(C[i,] > 0))
    posterior[i,0:max(C[i,])]<- 0

  # Compute posterior mean (BUP)
   bup2[i] <- sum(posterior[i,] * Ngrid)
}

# Compare BUPS with true N and counts for first and last 5 sites
(bup1 <- bup(ranef(fm.Nmix1)))
cbind(N=N,count1=C[,1],count2=C[,2],count3=C[,3],BUP1=bup1,BUP2=bup2)[c(1:5, 91:95),]


plot(ranef(fm.Nmix1), xlim = c(0,12))[sort(sample(1:100, 12))]


# Main-effects ANCOVA: additive effects of factor and covariate
summary(fm.Nmix2 <- pcount(~ wind+time-1 ~ vegHt+hab-1, data=umf))

# Interaction-effects ANCOVA: multiplicative effects of factor and covariate
summary(fm.Nmix3 <- pcount(~ wind*time-1-wind ~ vegHt*hab-1-vegHt, data=umf))

# Get predictions for factor levels at average values of covariates
newdat <- data.frame(vegHt=0, hab = c("A", "B", "C"))
predict(fm.Nmix2, type="state", newdata=newdat, appendData = T)


newdat <- data.frame(time = c("1", "2", "3"), wind = 0)
predict(fm.Nmix3, type="det", newdata=newdat, appendData = T)


newdat <- data.frame(vegHt=seq(0, 2 ,by = 0.1),
    hab = factor("A", levels = c("A", "B", "C")))
predict(fm.Nmix2, type="state", newdata=newdat, appendData = TRUE)


LRT(fm.Nmix3, fm.Nmix1)


# Bundle data
win.data <- list(C = C, M = nrow(C), J = ncol(C), wind = wind, vegHt = vegHt,
    hab = as.numeric(factor(hab)), XvegHt = seq(-1, 1,, 100),
    Xwind = seq(-1, 1,,100) )
str(win.data)

# Specify model in BUGS language
cat(file = "model2.txt", "
model {
  # Priors
  for(k in 1:3){                # Loop over 3 levels of hab or time factors
    alpha0[k] ~ dunif(-10, 10) # Detection intercepts
    alpha1[k] ~ dunif(-10, 10) # Detection slopes
    beta0[k] ~ dunif(-10, 10)  # Abundance intercepts
    beta1[k] ~ dunif(-10, 10)  # Abundance slopes
  }

  # Likelihood
  # Ecological model for true abundance
  for (i in 1:M){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0[hab[i]] + beta1[hab[i]] * vegHt[i]
    # Some intermediate derived quantities
    critical[i] <- step(2-N[i])# yields 1 whenever N is 2 or less
    z[i] <- step(N[i]-0.5)     # Indicator for occupied site
    # Observation model for replicated counts
    for (j in 1:J){
      C[i,j] ~ dbin(p[i,j], N[i])
      logit(p[i,j]) <- alpha0[j] + alpha1[j] * wind[i,j]
    }
  }

  # Derived quantities
  Nocc <- sum(z[])         # Number of occupied sites among sample of M
  Ntotal <- sum(N[])       # Total population size at M sites combined
  Nhab[1] <- sum(N[1:33])  # Total abundance for sites in hab A
  Nhab[2] <- sum(N[34:66]) # Total abundance for sites in hab B
  Nhab[3] <- sum(N[67:100])# Total abundance for sites in hab C
  for(k in 1:100){         # Predictions of lambda and p ...
    for(level in 1:3){    #    ... for each level of hab and time factors
      lam.pred[k, level] <- exp(beta0[level] + beta1[level] * XvegHt[k])
      logit(p.pred[k, level]) <- alpha0[level] + alpha1[level] * Xwind[k]
    }
  }
  N.critical <- sum(critical[]) # Number of populations with critical size
}")

# Initial values
Nst <- apply(C, 1, max)+1   # Important to give good inits for latent N
inits <- function() list(N = Nst, alpha0 = rnorm(3), alpha1 = rnorm(3),
    beta0 = rnorm(3), beta1 = rnorm(3))

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "Nocc", "Ntotal", "Nhab",
    "N.critical", "lam.pred", "p.pred") # could also estimate N, bayesian counterpart to BUPs before: simply add "N" to the list

# MCMC settings
nc <- 3   ;   ni <- 22000   ;   nb <- 2000   ;   nt <- 10

# Call JAGS, time run (ART 1 min) and summarize posteriors
system.time(out <- jags(win.data, inits, params, "model2.txt", n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
)
traceplot(out, param = c('alpha0', 'alpha1', 'beta0', 'beta1', 'Nocc', 'Ntotal',
    'Nhab', 'N.critical'))
print(out, 2)


plot(table(out$sims.list$N.critical), xlab="Number of populations with critical size",
    ylab="Frequency", frame = FALSE)
abline(v = 74.5, col = "red", lwd = 3)


(metapop.extinction.risk <- mean(out$sims.list$N.critical > 74))


op <- par(mfrow = c(1,2), mar = c(5,5,3,2), cex.axis = 1.5, cex.lab = 1.5)
X <- seq(-1, 1,, 100)
plot(X, out$summary[219:318,1], xlab = "Vegetation Height",
    ylab = "Expected abundance (lambda)", ylim = c(0, 11),
    frame = FALSE, type = "l")
polygon(c(X, rev(X)), c(out$summary[219:318,3], rev(out$summary[219:318,7])),
    col = "gray", border = FALSE)
lines(X, out$summary[219:318,1], lty = 1, lwd = 3, col = "blue")
plot(X, out$summary[519:618,1], xlab = "Wind speed",
    ylab = "Detection probability (p)", ylim = c(0, 1),
    frame = FALSE, type = "l")
polygon(c(X, rev(X)), c(out$summary[519:618,3], rev(out$summary[519:618,7])),
    col = "gray", border = FALSE)
lines(X, out$summary[519:618,1], lty = 1, lwd = 3, col = "blue")
par(op)
