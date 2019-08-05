#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 5. Fitting models using the Bayesian modeling software BUGS and JAGS
# =========================================================================

# 5.4 The R package rjags
# ------------------------------------------------------------------------


library(rjags)
load.module("glm")     # Careful with that package, see JAGS discussion list
load.module("dic")

# Have to explicitly list the deviance if want samples
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "deviance")

# Adaptative phase to maximize MCMC efficiency
model <- jags.model(file = "multiple_linear_regression_model.txt", data = win.data, inits = inits, n.chains = nc, n.adapt = 1000)

# Burnin
update(model, nb)

# Generate posterior samples
samples <- coda.samples(model = model, variable.names = params, n.iter = ni - nb, thin = nt)

# Get the summary statistics for the posterior samples
summfit <- summary(samples)
print(summfit, 2)

# Traceplots and posterior densities
plot(samples[,1:4])

# Compute the Brooks-Gelman-Rubin statistic (R-hat)
gelman.diag(samples)

# Compute the effective sample size
effectiveSize(samples)


# Secondary burnin can be applied (e.g. another 500 samples tossed out)
#samples <- window(samples, start = nb + 500 + 1, end = ni)

# More samples can be drawn (starting where the old chains stopped, not starting from 0)
newsamples <- coda.samples(model = model, variable.names = params, n.iter = 1500, thin = nt)

# Combine the new samples with the old ones (ugly but works)
mc1 <- as.mcmc(rbind(samples[[1]], newsamples[[1]]))
mc2 <- as.mcmc(rbind(samples[[2]], newsamples[[2]]))
mc3 <- as.mcmc(rbind(samples[[3]], newsamples[[3]]))
allsamples <- as.mcmc.list(list(mc1, mc2, mc3))

# Mean deviance
Dbar <- summfit$statistics["deviance","Mean"]

# Variance of the deviance
varD <- summfit$statistics["deviance","SD"]^2

# Compute pD and DIC (according to A. Gelman, implemented in R2jags)
pD <- varD/2
DIC <- Dbar + pD

# Another DIC computation (according to M. Plummer). DIC = Penalized deviance
(dic.pD <- dic.samples(model, 2000, "pD"))

