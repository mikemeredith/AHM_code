#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 5. Fitting models using the Bayesian modeling software BUGS and JAGS
# =========================================================================

library(AHMbook)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14/"          # Place where your WinBUGS installed
library(jagsUI)

# ~~~~~ this section requires the following code from section 5.3 ~~~~~~~~~~
# Generate data with data.fn from chapter 4
set.seed(24)
data <- data.fn(show.plot=FALSE)
attach(data)

# Summarize data by taking mean at each site and plot
Cmean <- apply(C, 1, mean)

# Write text file with model description in BUGS language
cat(file = "multiple_linear_regression_model.txt",
"   # --- Code in BUGS language starts with this quotation mark ---
model {

  # Priors
  alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
  alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
  alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
  alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
  tau <- pow(sd, -2)                   # Precision tau = 1/(sd^2)
  sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale

  # Likelihood
  for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)      # dispersion tau is precision (1/variance)
    mu[i] <- alpha0 + alpha1*elev[i] + alpha2*forest[i] + alpha3*elev[i]*forest[i]
  }

  # Derived quantities
  for (i in 1:M){
    resi[i] <- Cmean[i] - mu[i]
  }
}"#  --- Code in BUGS language ends on this line ---
)

# Initial values (have to give for at least some estimands)
inits <- function() list(alpha0 = rnorm(1,0,10), alpha1 = rnorm(1,0,10),
    alpha2 = rnorm(1,0,10), alpha3 = rnorm(1,0,10))

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.5 Missing values in a Bayesian analysis
# =========================================


# 5.5.1 Some responses missing
# ------------------------------------------------------------------------
# Copy mean counts and turn first 10 into NAs
Cm <- Cmean         # Copy Cmean into Cm
Cm[1:10] <- NA      # turn first 10 into missing

# Bundle data (inside BUGS use Cm for Cmean)
win.data <- list(Cmean = Cm, M = length(Cm), elev = elev, forest = forest)

# Parameters monitored (i.e., for which estimates are saved)
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "Cmean", "mu")

# … or this to get a subset of the parameters
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "Cmean[1:10]", "mu[1:10]")

# Call WinBUGS or JAGS from R (ART <1 min) and summarize posteriors
out1.1 <- bugs(win.data, inits, params, "multiple_linear_regression_model.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd(), DIC = FALSE)
  debug = FALSE, bugs.directory = bugs.dir, DIC = FALSE) # ~~~ for autotesting

out1.1 <- jags(win.data, inits, params, "multiple_linear_regression_model.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(out1.1, 2)

print(cbind(Truth=Cmean[1:10], out1.1$summary[6:15,c(1:3,7)],
    out1.1$summary[16:25, c(1:3,7)]),3)


# 5.5.2 All responses missing
# ------------------------------------------------------------------------
# Bundle data: simply drop the response from list
win.data <- list(M = length(Cm), elev = elev, forest = forest)

# Alternatively, add all-NA data vector
win.data$Cmean <- as.numeric(rep(NA, length(Cmean)))
str(win.data)  # Cmean is numeric

# Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "mu[1:2]")

# Call WinBUGS from R (ART <1 min) and summarize posteriors
out1.2 <- bugs(win.data, inits, params, "multiple_linear_regression_model.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd(), DIC = FALSE)
  debug = FALSE, bugs.directory = bugs.dir,DIC = FALSE) # ~~~ for autotesting

print(out1.2, 2)

op <- par(mfrow = c(3, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(out1.2$sims.list$alpha0, breaks = 100, col = "grey", main = "")
text(-3500, 600, "A", cex = 2)
hist(out1.2$sims.list$alpha1, breaks = 100, col = "grey", main = "", ylab = "")
text(-3500, 580, "B", cex = 2)
hist(out1.2$sims.list$alpha2, breaks = 100, col = "grey", main = "")
text(-3700, 580, "C", cex = 2)
hist(out1.2$sims.list$sd, breaks = 100, col = "grey", main = "",
    ylim = c(0, 230), ylab = "")
text(30, 220, "D", cex = 2)
hist(out1.2$sims.list$mu[,1], breaks = 100, col = "grey", main = "")
text(-4000, 480, "E", cex = 2)
hist(out1.2$sims.list$mu[,2], breaks = 100, col = "grey", main = "", ylab = "")
text(-4000, 480, "F", cex = 2)
par(op)


# 5.5.3 Missing values in a covariate
# ------------------------------------------------------------------------
# Shoot 'holes' in the covariate data
ele <- elev          # copy of elevation covariate
ele[1:10] <- NA      # missing values in covariate elevation

# Bundle data: feed new 'ele' into 'elev' covariate inside of BUGS model
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = ele, forest = forest)

# Specify model in BUGS language
cat(file = "missing_cov_imputation_model_1.txt","
model {

  # Priors
  alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
  alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
  alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
  alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
  tau <- pow(sd, -2)
  sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale

  # Likelihood
  for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)      # precision tau = 1 / variance
    mu[i] <- alpha0 + alpha1 * elev[i] + alpha2 * forest[i] + alpha3 * elev[i] * forest[i]
  }

  # Model for missing covariates
  for (i in 1:M){
    elev[i] ~ dnorm(0, 0.001)
  }
}")

# Initial values
inits <- function() list(alpha0 = rnorm(1,,10), alpha1 = rnorm(1,,10),
    alpha2 = rnorm(1,,10), alpha3 = rnorm(1,,10))

# Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "elev")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS from R (ART <1 min)
out1.3 <- bugs(win.data, inits, params, "missing_cov_imputation_model_1.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir) # ~~~~ for autotesting


# Specify model in BUGS language
cat(file = "missing_cov_imputation_model_2.txt","
model {

  # Priors
  alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
  alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
  alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
  alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
  tau <- pow(sd, -2)
  sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale

  # Likelihood
  for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)      # precision tau = 1 / variance
    mu[i] <- alpha0 + alpha1 * elev[i] + alpha2 * forest[i] + alpha3 * elev[i] * forest[i]
  }

  # Covariate mean as a model for missing covariates
  for (i in 1:M){
    elev[i] ~ dnorm(mu.elev, tau.elev)    # Assume elevation normally distributed
  }
  mu.elev ~ dnorm(0, 0.0001)
  tau.elev <- pow(sd.elev, -2)
  sd.elev ~ dunif(0, 100)
}")

# Initial values
inits <- function() list(alpha0 = rnorm(1,,10), alpha1 = rnorm(1,,10),
    alpha2 = rnorm(1,,10), alpha3 = rnorm(1,,10))

# Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "elev", "mu.elev", "sd.elev")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS from R (ART <1 min)
out1.4 <- bugs(win.data, inits, params, "missing_cov_imputation_model_2.txt",
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
  # debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
  debug = FALSE, bugs.directory = bugs.dir) # ~~~ for autotesting

op <- par(cex = 1.5, lwd = 2)
plot(elev[1:10]-0.01, out1.3$summary[6:15,1], ylim = c(-10, 10), col = "red",
    xlab = "True value of covariate", ylab = "Estimate (with 95% CRI)",
    frame.plot =FALSE)
segments(elev[1:10]-0.01, out1.3$summary[6:15,3], elev[1:10]-0.01,
    out1.3$summary[6:15,7], col = "red")
points(elev[1:10]+0.01, out1.4$summary[6:15,1], ylim = c(-3, 3), col = "blue")
segments(elev[1:10]+0.01, out1.4$summary[6:15,3], elev[1:10]+0.01,
    out1.4$summary[6:15,7], col = "blue")
abline(0,1)
par(op)
