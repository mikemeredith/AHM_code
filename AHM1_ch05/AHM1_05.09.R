#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 5. Fitting models using the Bayesian modeling software BUGS and JAGS
# =========================================================================

# 5.9 Poisson generalized linear model (Poisson GLM)
# --------------------------------------------------


# Summarize data by taking max at each site
Cmax <- apply(C, 1, max)
table(Cmax)

# Bundle data
win.data <- list(Cmax = Cmax, M = length(Cmax), elev = elev, facFor = facFor, e = 0.0001)

# Specify model in BUGS language
cat(file = "Poisson_GLM.txt","
model {

# Priors
for(k in 1:4){
   alpha[k] ~ dnorm(0, 1.0E-06)       # Prior for intercepts
   beta[k] ~ dnorm(0, 1.0E-06)        # Prior for slopes
}

# Likelihood
for (i in 1:M){
   Cmax[i] ~ dpois(lambda[i])         # note no variance parameter
   log(lambda[i]) <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
   resi[i] <- (Cmax[i]-lambda[i]) / (sqrt(lambda[i])+e)   # Pearson resi
}
}
")


# Initial values
inits <- function() list(alpha = rnorm(4,,3), beta = rnorm(4,,3))

# Parameters monitored
params <- c("alpha", "beta", "lambda", "resi")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS or JAGS from R and summarize posteriors
out5 <- bugs(win.data, inits, params, "Poisson_GLM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

system.time(out5J <- jags(win.data, inits, params, "Poisson_GLM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb))
par(mfrow = c(4,2))    ;    traceplot(out5J, c("alpha[1:4]", "beta[1:4]"))
print(out5J, 3)

par(mfrow = c(1, 3), mar = c(5,5,3,2), cex = 1.3, cex.lab = 1.5, cex.axis = 1.5)
hist(out5$summary[276:542, 1], xlab = "Pearson residuals", col = "grey", breaks = 50, main = "", freq = F, xlim = c(-5, 5), ylim = c(0, 0.57))
abline(v = 0, col = "red", lwd = 2)
text(-4.7, 0.54, "A", cex = 1.5)

plot(1:267, out5$summary[276:542, 1], main = "", xlab = "Order of data", ylab = "Pearson residual", frame.plot = F)
abline(h = 0, col = "red", lwd = 2)
text(8, 4, "B", cex = 1.5)

plot(out5$summary[9:275, 1],out5$summary[276:542, 1], main = "", xlab = "Predicted values", ylab = "Pearson residual", frame.plot = F, xlim = c(-1, 14))
abline(h = 0, col = "red", lwd = 2)
text(-0.5, 4, "C", cex = 1.5)


summary(glm(Cmax ~ factor(facFor)*elev-1-elev, family = poisson))


lambda2 <- array(dim = c(15000, 267))
for(j in 1:267){                            # Loop over sites
   lambda2[,j] <- exp(out5$sims.list$alpha[,facFor[j]] + out5$sims.list$beta[,facFor[j]] * elev[j]) # linear regression/backtransform
}
plot(out5$sims.list$lambda ~ lambda2, pch = ".")  # Check the two are identical
lm(c(out5$sims.list$lambda) ~ c(lambda2))

sorted.ele1 <- sort(elev[facFor == 1])
sorted.y1 <- out5$summary[9:275,][facFor == 1,][order(elev[facFor == 1]),]

# Plot A
par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(elev[facFor == 1], jitter(Cmax[facFor ==1]), ylab = "Maximum count", xlab = "Elevation (scaled)", frame.plot=F), ylim = c(0, 6))
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2) # Post. mean
lines(sorted.ele1, sorted.y1[,3], col = "grey", lwd = 2) # Lower 95% CL
lines(sorted.ele1, sorted.y1[,7], col = "grey", lwd = 2) # Upper 95% CL
text(-0.8, 6, "A", cex = 2)

# Plot B
plot(sorted.ele1, sorted.y1[,1], type='n', xlab = "Elevation (scaled)", ylab = "", frame.plot = F, ylim = c(0, 6))
polygon(c(sorted.ele1, rev(sorted.ele1)), c(sorted.y1[,3], rev(sorted.y1[,7])), col='grey', border=NA)
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2)
text(-0.8, 6, "B", cex = 2)

# Plot C
elev.pred <- seq(-1,1, length.out = 200)  # Cov. for which to predict lambda
n.pred <- 50                             # Number of prediction profiles
pred.matrix <- array(NA, dim = c(length(elev.pred), n.pred))
for(j in 1:n.pred){
   sel <- sample(1:length(out5$sims.list$alpha[,1]),1) # Choose one post. draw
   pred.matrix[,j] <- exp(out5$sims.list$alpha[sel,1] + out5$sims.list$beta[sel,1] * elev.pred)
}
plot(sorted.ele1, sorted.y1[,1], type='n', xlab = "Elevation (scaled)", ylab = "", frame.plot = F, ylim = c(0, 6))
matlines(elev.pred, pred.matrix, col = "grey", lty = 1, lwd = 1)
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2)
text(-0.8, 6, "C", cex = 2)

