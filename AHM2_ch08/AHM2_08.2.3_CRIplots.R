#   Applied hierarchical modeling in ecology, vol.2, 2021
#   Marc Kéry & J. Andy Royle
#
# Chapter 8 : MODELING INTERACTIONS AMONG SPECIES
# ===============================================
# This script is by Mike Meredith.

library(jagsUI)

# 8.2 Joint occupancy models for few species with symmetric interactions
# ======================================================================

# 8.2.3 Implementation in JAGS
# ----------------------------

# Figures 8.4 to 8.7 with CRIs
# ''''''''''''''''''''''''''''
# In the book (p. 479) K&R say "We work with the posterior means only; if you want prediction SEs and CRIs, you must work with the full list of posterior draws instead, and most quantities in the code that follows here will need one added dimension."
# We need to derive MCMC chains for each of the quantities to be plotted, and then summarize to get means and CRIs as the last step.

# Load the output from the main script
load(file="AHM2_08.2.3.Rdata")

# Number of draws in the MCMC chains
ndraws <- out1$mcmc.info$n.samples

# Number of points to plot
npoints <- 100  # 100 enough for a decent plot, increase for publication-quality graphics


# Occupancy as a function of Disturbance
# (at average Housing density)
# ......................................

# Create prediction covariates
(r <- range(data$sitecovs[,'Dist_5km']))
(mnx <- mean(data$sitecovs[,'Dist_5km']))
(sdx <- sd(data$sitecovs[,'Dist_5km']))
Dist.pred.orig <- seq(r[1], r[2], length.out = npoints) # Create new data
Dist.pred <- (Dist.pred.orig-mnx) / sdx                 # ... and scale as the real data
HDens.pred <- rep(0, npoints)

# Create model matrices for the prediction covariates
head(psi_cov1 <- cbind(1, Dist.pred))    # We vary this covariate ...
head(psi_cov2 <- cbind(1, HDens.pred))   # .. while this stays constant


### Assemble the linear predictors
# The matrix multiplication is different to the book because (1) the coefficients are now MCMC chains in a 3000 x 2 matrix (instead of a length 2 vector), and (2) the output should be MCMC chains as columns in a 3000-row matrix (instead of a vector of point estimates).

# ...for states B, C, and F
psiB <- out1$sims.list$betaB %*% t(psi_cov1)  # produces 3000 x 100 matrix
psiC <- out1$sims.list$betaC %*% t(psi_cov1)
psiF <- out1$sims.list$betaF %*% t(psi_cov1)

# ...for states BC, CF, and BF (in that order)
psiBC <- psiB + psiC + out1$sims.list$betaBC %*% t(psi_cov2)
psiCF <- psiC + psiF + out1$sims.list$betaCF %*% t(psi_cov2)
psiBF <- psiB + psiF + out1$sims.list$betaBF %*% t(psi_cov2)

# ...for state BCF
psiBCF <- psiB + psiC + psiF + psiBC + psiCF + psiBF

### Assemble the state probability array
# Compute latent state vector (lsv)
lsv <- array(NA, dim=c(ndraws, npoints, ncat))  # 3000 x 100 x 8
lsv[, , 1] <- 1              # none of the three species
lsv[, , 2] <- exp( psiB )    # only B
lsv[, , 3] <- exp( psiC )    # only C
lsv[, , 4] <- exp( psiF )    # only F
lsv[, , 5] <- exp( psiBC)    # both B and C
lsv[, , 6] <- exp( psiCF)    # both C and F
lsv[, , 7] <- exp( psiBF)    # both B and F
lsv[, , 8] <- exp( psiBCF )  # all three species

# Probability of each state as MCMC chains
# denominator = sum  of lsv for each draw x point
denom <- apply(lsv, 1:2, sum)
lspMC <- sweep(lsv, 1:2, denom, "/")
dimnames(lspMC)[[3]] <- c('U', 'B', 'C', 'F', 'BC', 'CF', 'BF', 'BCF')
str(lspMC)

# Fill in summary values for 'lsp' for Fig. 8.4
lsp <- array(NA, dim=c(npoints, ncat, 4))
dimnames(lsp) <- list(NULL, c('U', 'B', 'C', 'F', 'BC', 'CF', 'BF', 'BCF'),
    c("mean", "SD", "lower", "upper"))
lsp[, , 1] <- apply(lspMC, 2:3, mean)
lsp[, , 2] <- apply(lspMC, 2:3, sd)
lsp[, , 3] <- apply(lspMC, 2:3, quantile, probs=0.025)
lsp[, , 4] <- apply(lspMC, 2:3, quantile, probs=0.975)

# Plot the latent state probabilities as a function of covariate Dist (Fig. 8.4)
matplot(Dist.pred.orig, lsp[,,'mean'], type = 'l', lty = 1, lwd = 5, col = 1:8,
    frame = FALSE, xlab = 'Disturbance', ylab = 'Probability of state',
    ylim = c(0, max(lsp[,,-2])), las = 1)
legend('topleft', lwd = 3, lty = 1, col = 1:8, colnames(lsp), bty = 'n',
    cex = 1.2, ncol = 2)
# Add CrIs
for(i in 1:8) {
  polygon(c(Dist.pred.orig, rev(Dist.pred.orig)), c(lsp[,i,'lower'], rev(lsp[,i,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}

# ... or do separate plots
op <- par(mfrow=c(2,4))
for(i in 1:8) {
  plot(Dist.pred.orig, lsp[,i,'mean'], type = 'n',
      frame = FALSE, xlab = 'Disturbance', ylab = 'Probability of state',
      ylim = c(0, max(lsp[,,-2])), las = 1)
  polygon(c(Dist.pred.orig, rev(Dist.pred.orig)), c(lsp[,i,'lower'], rev(lsp[,i,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
  lines(Dist.pred.orig, lsp[,i,'mean'], lwd = 2, col = i)
  legend('topleft', legend=colnames(lsp)[i], lwd = 2, col = i, bty = 'n')
}
par(op)

### Compute marginal occupancy for all three species
# Get the MCMC chains
marginalMC <- array(NA, dim=c(ndraws, npoints, 3))
dimnames(marginalMC)[[3]] <- c("B", "C", "F")
marginalMC[,,'B'] <- apply(lspMC[, ,c('B', 'BC', 'BF', 'BCF')], 1:2, sum) # Bobcat
marginalMC[,,'C'] <- apply(lspMC[, ,c('C', 'BC', 'CF', 'BCF')], 1:2, sum) # Coyote
marginalMC[,,'F'] <- apply(lspMC[, ,c('F', 'CF', 'BF', 'BCF')], 1:2, sum) # Fox

# Fill in summary values for 'marginal'
marginal <- array(NA, dim=c(npoints, 3, 4))
dimnames(marginal) <- list(NULL, c("bobcat", "coyote", "fox"),
    c("mean", "SD", "lower", "upper"))
marginal[, , 1] <- apply(marginalMC, 2:3, mean)
marginal[, , 2] <- apply(marginalMC, 2:3, sd)
marginal[, , 3] <- apply(marginalMC, 2:3, quantile, probs=0.025)
marginal[, , 4] <- apply(marginalMC, 2:3, quantile, probs=0.975)

# Plot marginal occupancy for all three species (not shown)
matplot(Dist.pred.orig, marginal[,,'mean'],
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE, xlab = 'Disturbance',
    ylab = 'Marginal occupancy probability', ylim = c(0, max(marginal[,,'upper'])))
legend('top', c('Bobcat', 'Coyote', 'Red Fox'), lwd = 3, lty = 1:3, col = 1:3,
    horiz = TRUE, bty = 'n')
# Add CrIs
for(i in 1:3) {
  polygon(c(Dist.pred.orig, rev(Dist.pred.orig)), c(marginal[,i,'lower'], rev(marginal[,i,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}

### Compute conditional occupancy for all species
# get MCMC chains
conditionalMC <- array(NA, dim=c(ndraws, npoints, 9))
dimnames(conditionalMC)[[3]] <- c("B.alone", "B.given.coyote", "B.given.fox",
                                  "C.alone", "C.given.bobcat", "C.given.fox",
                                  "F.alone", "F.given.bobcat", "F.given.coyote")

conditionalMC[,,'B.alone'] <- lspMC[, , 'B'] / apply(lspMC[, , c('U', 'B')], 1:2, sum)
conditionalMC[,,'B.given.coyote'] <- apply(lspMC[, , c('BC', 'BCF')], 1:2, sum) / marginalMC[,,'C']
conditionalMC[,,'B.given.fox'] <- apply(lspMC[, , c('BF', 'BCF')], 1:2, sum) / marginalMC[,,'F']

conditionalMC[,,'C.alone'] <- lspMC[, , 'C'] / apply(lspMC[, , c('U', 'C')], 1:2, sum)
conditionalMC[,,'C.given.bobcat'] <- apply(lspMC[, , c('BC', 'BCF')], 1:2, sum) / marginalMC[,,'B']
conditionalMC[,,'C.given.fox'] <- apply(lspMC[, , c('CF', 'BCF')], 1:2, sum) / marginalMC[,,'F']

conditionalMC[,,'F.alone'] <- lspMC[, , 'F'] / apply(lspMC[, , c('U', 'F')], 1:2, sum)
conditionalMC[,,'F.given.bobcat'] <- apply(lspMC[, , c('BF', 'BCF')], 1:2, sum) / marginalMC[,,'B']
conditionalMC[,,'F.given.coyote'] <- apply(lspMC[, , c('CF', 'BCF')], 1:2, sum) / marginalMC[,,'C']

# Fill in summary values for Fig. 8.5
conditional <- array(NA, dim=c(npoints, 9, 4))
dimnames(conditional) <- list(NULL, c("B.alone", "B.given.coyote", "B.given.fox",
                                      "C.alone", "C.given.bobcat", "C.given.fox",
                                      "F.alone", "F.given.bobcat", "F.given.coyote"),
                              c("mean", "SD", "lower", "upper"))
conditional[, , 1] <- apply(conditionalMC, 2:3, mean)
conditional[, , 2] <- apply(conditionalMC, 2:3, sd)
conditional[, , 3] <- apply(conditionalMC, 2:3, quantile, probs=0.025)
conditional[, , 4] <- apply(conditionalMC, 2:3, quantile, probs=0.975)

# Plot conditional probabilities of occupancy (Fig. 8.5)
op <- par(mfrow = c(1, 3))
matplot(Dist.pred.orig, conditional[, 1:3, 'mean'],
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Disturbance', ylab = 'Occupancy probability',
    ylim = c(0, 0.6), main = 'Bobcat')
legend('topright', c('alone', 'given Coyote presence', 'given Fox presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
for(i in 1:3) {
  polygon(c(Dist.pred.orig, rev(Dist.pred.orig)),
      c(conditional[,i,'lower'], rev(conditional[,i,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}
matplot(Dist.pred.orig, conditional[, 4:6, 'mean'],
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Disturbance', ylab = 'Occupancy probability',
    ylim = c(0, 1), main = 'Coyote')
legend('bottomright', c('alone', 'given Bobcat presence', 'given Fox presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
for(i in 1:3) {
  polygon(c(Dist.pred.orig, rev(Dist.pred.orig)),
      c(conditional[,i+3,'lower'], rev(conditional[,i+3,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}
matplot(Dist.pred.orig, conditional[, 7:9, 'mean'],
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Disturbance', ylab = 'Occupancy probability',
    ylim = c(0, 0.6), main = 'Red Fox')
legend('topright', c('alone', 'given Bobcat presence', 'given Coyote presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
for(i in 1:3) {
  polygon(c(Dist.pred.orig, rev(Dist.pred.orig)),
      c(conditional[,i+6,'lower'], rev(conditional[,i+6,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}
par(op)


# Occupancy as a function of Housing density
# (at average Disturbance)
# ......................................

# Create prediction covariates
(r <- range(data$sitecovs[,'HDens_5km']))
mnx <- mean(data$sitecovs[,'HDens_5km'])
sdx <- sd(data$sitecovs[,'HDens_5km'])
HDens.pred.orig <- seq(r[1], 120, length.out = npoints)
HDens.pred <- (HDens.pred.orig-mnx) / sdx
Dist.pred <- rep(0, npoints)

# Create model matrices for the prediction covariates
head(psi_cov1 <- cbind(1, Dist.pred))    # Now this covariate is fixed ...
head(psi_cov2 <- cbind(1, HDens.pred))   # .. while this varies

## The following code is the same as above, only the model matrices differ ##

### Assemble the linear predictors
# ...for states B, C, and F
psiB <- out1$sims.list$betaB %*% t(psi_cov1)  # 3000 x 100
psiC <- out1$sims.list$betaC %*% t(psi_cov1)
psiF <- out1$sims.list$betaF %*% t(psi_cov1)

# ...for states BC, CF, and BF (in that order)
psiBC <- psiB + psiC + out1$sims.list$betaBC %*% t(psi_cov2)
psiCF <- psiC + psiF + out1$sims.list$betaCF %*% t(psi_cov2)
psiBF <- psiB + psiF + out1$sims.list$betaBF %*% t(psi_cov2)

# ...for state BCF
psiBCF <- psiB + psiC + psiF + psiBC + psiCF + psiBF

### Assemble the state probability array
# Compute latent state vector (lsv)
lsv <- array(NA, dim=c(ndraws, npoints, ncat))  # 3000 x 100 x 8
lsv[, , 1] <- 1              # none of the three species
lsv[, , 2] <- exp( psiB )    # only B
lsv[, , 3] <- exp( psiC )    # only C
lsv[, , 4] <- exp( psiF )    # only F
lsv[, , 5] <- exp( psiBC)    # both B and C
lsv[, , 6] <- exp( psiCF)    # both C and F
lsv[, , 7] <- exp( psiBF)    # both B and F
lsv[, , 8] <- exp( psiBCF )  # all three species

# Probability of each state as MCMC chains
# denominator = sum  of lsv for each draw x point
denom <- apply(lsv, 1:2, sum)
lspMC <- sweep(lsv, 1:2, denom, "/")
dimnames(lspMC)[[3]] <- c('U', 'B', 'C', 'F', 'BC', 'CF', 'BF', 'BCF')
str(lspMC)

# Fill in summary values for 'lsp' for Fig. 8.6 left
lsp <- array(NA, dim=c(npoints, ncat, 4))
dimnames(lsp) <- list(NULL, c('U', 'B', 'C', 'F', 'BC', 'CF', 'BF', 'BCF'),
    c("mean", "SD", "lower", "upper"))
lsp[, , 1] <- apply(lspMC, 2:3, mean)
lsp[, , 2] <- apply(lspMC, 2:3, sd)
lsp[, , 3] <- apply(lspMC, 2:3, quantile, probs=0.025)
lsp[, , 4] <- apply(lspMC, 2:3, quantile, probs=0.975)

### Compute marginal occupancy for all three species
# Get the MCMC chains
marginalMC <- array(NA, dim=c(ndraws, npoints, 3))
dimnames(marginalMC)[[3]] <- c("B", "C", "F")
marginalMC[,,'B'] <- apply(lspMC[, ,c('B', 'BC', 'BF', 'BCF')], 1:2, sum) # Bobcat
marginalMC[,,'C'] <- apply(lspMC[, ,c('C', 'BC', 'CF', 'BCF')], 1:2, sum) # Coyote
marginalMC[,,'F'] <- apply(lspMC[, ,c('F', 'CF', 'BF', 'BCF')], 1:2, sum) # Fox

# Fill in summary values for 'marginal'
marginal <- array(NA, dim=c(npoints, 3, 4))  # this figure not shown
dimnames(marginal) <- list(NULL, c("bobcat", "coyote", "fox"),
    c("mean", "SD", "lower", "upper"))
marginal[, , 1] <- apply(marginalMC, 2:3, mean)
marginal[, , 2] <- apply(marginalMC, 2:3, sd)
marginal[, , 3] <- apply(marginalMC, 2:3, quantile, probs=0.025)
marginal[, , 4] <- apply(marginalMC, 2:3, quantile, probs=0.975)

# Plot the lsp and marginal occupancy of each species as a function of
#   covariate HDens while keeping Dist constant (Fig. 8.6)
op <- par(mfrow = c(1, 2))
matplot(HDens.pred.orig, lsp[,,'mean'], type = 'l', lty = 1, lwd = 5, col = 1:8,
    frame = FALSE, xlab = 'Housing density', ylab = 'Probability',
    ylim = c(0, 1.05), las = 1, main = 'Probability of each state')
legend('right', lwd = 3, lty = 1, col = 1:8, colnames(lsp), bty = 'n')
for(i in 1:8) {
  polygon(c(HDens.pred.orig, rev(HDens.pred.orig)), c(lsp[,i,'lower'], rev(lsp[,i,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}

matplot(HDens.pred.orig, marginal[,,'mean'],
    type = 'l', lty = 1:3, lwd = 5, col = 1:3, frame = FALSE,
    xlab = 'Housing density', ylab = 'Probability', ylim = c(0, 1),
    main = 'Marginal occupancy probability')
legend('right', c('Bobcat', 'Coyote', 'Red Fox'), lwd = 3, lty = 1:3,
    col = 1:3, bty = 'n')
for(i in 1:3) {
  polygon(c(HDens.pred.orig, rev(HDens.pred.orig)), c(marginal[,i,'lower'], rev(marginal[,i,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}
par(op)

### Compute conditional occupancy for all species
# get MCMC chains
conditionalMC <- array(NA, dim=c(ndraws, npoints, 9))
dimnames(conditionalMC)[[3]] <- c("B.alone", "B.given.coyote", "B.given.fox",
                                  "C.alone", "C.given.bobcat", "C.given.fox",
                                  "F.alone", "F.given.bobcat", "F.given.coyote")
conditionalMC[,,'B.alone'] <- lspMC[, , 'B'] / apply(lspMC[, , c('U', 'B')], 1:2, sum)
conditionalMC[,,'B.given.coyote'] <- apply(lspMC[, , c('BC', 'BCF')], 1:2, sum) / marginalMC[,,'C']
conditionalMC[,,'B.given.fox'] <- apply(lspMC[, , c('BF', 'BCF')], 1:2, sum) / marginalMC[,,'F']

conditionalMC[,,'C.alone'] <- lspMC[, , 'C'] / apply(lspMC[, , c('U', 'C')], 1:2, sum)
conditionalMC[,,'C.given.bobcat'] <- apply(lspMC[, , c('BC', 'BCF')], 1:2, sum) / marginalMC[,,'B']
conditionalMC[,,'C.given.fox'] <- apply(lspMC[, , c('CF', 'BCF')], 1:2, sum) / marginalMC[,,'F']

conditionalMC[,,'F.alone'] <- lspMC[, , 'F'] / apply(lspMC[, , c('U', 'F')], 1:2, sum)
conditionalMC[,,'F.given.bobcat'] <- apply(lspMC[, , c('BF', 'BCF')], 1:2, sum) / marginalMC[,,'B']
conditionalMC[,,'F.given.coyote'] <- apply(lspMC[, , c('CF', 'BCF')], 1:2, sum) / marginalMC[,,'C']

# Fill in summary values for Fig. 8.7
conditional <- array(NA, dim=c(npoints, 9, 4))
dimnames(conditional) <- list(NULL, c("B.alone", "B.given.coyote", "B.given.fox",
                                      "C.alone", "C.given.bobcat", "C.given.fox",
                                      "F.alone", "F.given.bobcat", "F.given.coyote"),
                              c("mean", "SD", "lower", "upper"))
conditional[, , 1] <- apply(conditionalMC, 2:3, mean)
conditional[, , 2] <- apply(conditionalMC, 2:3, sd)
conditional[, , 3] <- apply(conditionalMC, 2:3, quantile, probs=0.025)
conditional[, , 4] <- apply(conditionalMC, 2:3, quantile, probs=0.975)

# Plot conditional probabilities of occupancy (Fig. 8.7)
op <- par(mfrow = c(1, 3))
matplot(HDens.pred.orig, conditional[, 1:3, 'mean'],
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Housing density', ylab = 'Occupancy probability',
    ylim = c(0, 1), main = 'Bobcat')
legend('topright', c('alone', 'given Coyote presence', 'given Fox presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
for(i in 1:3) {
  polygon(c(HDens.pred.orig, rev(HDens.pred.orig)),
      c(conditional[,i,'lower'], rev(conditional[,i,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}
matplot(HDens.pred.orig, conditional[, 4:6, 'mean'],
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Housing density', ylab = 'Occupancy probability',
    ylim = c(0, 1), main = 'Coyote')
legend('right', c('alone', 'given Bobcat presence', 'given Fox presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
for(i in 1:3) {
  polygon(c(HDens.pred.orig, rev(HDens.pred.orig)),
      c(conditional[,i+3,'lower'], rev(conditional[,i+3,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}
matplot(HDens.pred.orig, conditional[, 7:9, 'mean'],
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Housing density', ylab = 'Occupancy probability',
    ylim = c(0, 1), main = 'Red Fox')
legend('right', c('alone', 'given Bobcat presence', 'given Coyote presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
for(i in 1:3) {
  polygon(c(HDens.pred.orig, rev(HDens.pred.orig)),
      c(conditional[,i+6,'lower'], rev(conditional[,i+6,'upper'])),
      border=NA, col=adjustcolor(i, alpha=0.2))
}
par(op)
