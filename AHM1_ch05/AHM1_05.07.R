#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 5. Fitting models using the Bayesian modeling software BUGS and JAGS
# =========================================================================

# 5.7 Proportion of variance explained (R2)
# =========================================

cat(file = "Model0.txt","
model {
# Priors
mu ~ dnorm(0, 1.0E-06)
tau <- pow(sd, -2)
sd ~ dunif(0, 1000)
# Likelihood
for (i in 1:M){
   Cmean[i] ~ dnorm(mu, tau)
}
}
")
inits <- function() list(mu = rnorm(1))
params <- c("mu", "sd")
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3
out0 <- jags(win.data, inits, params, "Model0.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(out0)

# Compute R2 from BUGS analysis
(total.var <- mean(out0$sims.list$sd^2))       # Total variance around the mean
(unexplained.var <- mean(out3$sims.list$sd^2)) # Not explained by the ANCOVA
(prop.explained <- (total.var - unexplained.var)/total.var)

