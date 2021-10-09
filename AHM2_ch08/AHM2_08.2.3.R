#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 2: Dynamic and Advanced models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8 : MODELING INTERACTIONS AMONG SPECIES
# ===============================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 40 mins
# Run time with the full number of iterations: 5 hrs

library(AHMbook)
library(abind)
library(jagsUI)

# 8.2 Joint occupancy models for few species with symmetric interactions
# ======================================================================

# 8.2.3 Implementation in JAGS
# ----------------------------

# Get data set and put into right format
library(AHMbook)
library(abind)
data(MesoCarnivores)
str(data <- MesoCarnivores)
str(ylist <- list(bobcat = data$bobcat, coyote = data$coyote,
    redfox = data$redfox))

# Put species detections in an array (site x rep x species)
y <- abind(ylist, along = 3)

# Get sample sizes
nsites <- dim(y)[1]       # 1437 sites
nsurveys <- dim(y)[2]     # 3 surveys
nspec <- dim(y)[3]        # 3 species
ncat <- 2^nspec           # 8 possible community states

# Prepare site covariates for occupancy and look at them
Dist <- scale(data$sitecovs[,'Dist_5km'])
HDens <- scale(data$sitecovs[,'HDens_5km'])
summary(data$sitecovs[,1]) ; summary(data$sitecovs[,2])

# Prepare one site covariate to serve as an observational covariate for detection
table(Trail <- data$sitecovs[,'Trail'])    # summarize
Trail <- matrix(cbind(Trail, Trail, Trail), ncol = 3)

# Condense multi-species detection array to be site x survey
ycat <- apply(y, c(1,2), paste, collapse = "")
ycat[ycat == "000"] <- 1    # Unoccupied (abbreviated 'U')
ycat[ycat == "100"] <- 2    # Only bobcat ('B') detected
ycat[ycat == "010"] <- 3    # Only coyote ('C') detected
ycat[ycat == "001"] <- 4    # Only red fox ('F') detected
ycat[ycat == "110"] <- 5    # Bobcat and coyote ('BC')
ycat[ycat == "011"] <- 6    # Coyote and fox ('CF')
ycat[ycat == "101"] <- 7    # Bobcat and fox ('BF')
ycat[ycat == "111"] <- 8    # All three species ('BCF')

# Convert each column to a numeric: this is our response variable
ycat <- apply(ycat, 2, as.numeric)

# Model matrix for first order occupancy (psi): main effects
psi_cov <- matrix(NA, ncol = 2, nrow = nsites)
psi_cov[,1] <- 1    # Intercept
psi_cov[,2] <- Dist # Slope of Dist

# Model matrix for second order occupancy (psi): 2-way interactions
psi_inxs_cov <- matrix(NA, ncol = 2, nrow = nsites)
psi_inxs_cov[,1] <- 1     # Intercept
psi_inxs_cov[,2] <- HDens # Slope of HDens

# Model matrix for first order detection (rho): main effects
rho_cov <- array(NA, dim = c(nsites, nsurveys, 2))
rho_cov[,,1] <- 1     # Intercept
rho_cov[,,2] <- Trail # Site off/on trail

# Model matrix for second order detection (rho): 2-way interactions
rho_inxs_cov <- rep(1, nsites) # Intercept

# Bundle and summarize data set
str(bdata <- list(y = ycat, psi_cov = psi_cov, psi_inxs_cov = psi_inxs_cov,
    rho_cov = rho_cov, rho_inxs_cov = rho_inxs_cov, nsites = nsites,
    nsurveys = nsurveys, nfirst_order_psi = ncol(psi_cov),
    nsecond_order_psi = ncol(psi_inxs_cov), nfirst_order_rho = dim(rho_cov)[3],
    nsecond_order_rho = 1, ncat = ncat) )
# List of 12
# $ y                : num [1:1437, 1:3] 1 1 1 1 1 4 1 1 1 1 ...
# $ psi_cov          : num [1:1437, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
# $ psi_inxs_cov     : num [1:1437, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
# $ rho_cov          : num [1:1437, 1:3, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
# $ rho_inxs_cov     : num [1:1437] 1 1 1 1 1 1 1 1 1 1 ...
# $ nsites           : int 1437
# $ nsurveys         : int 3
# $ nfirst_order_psi : int 2
# $ nsecond_order_psi: int 2
# $ nfirst_order_rho : int 2
# $ nsecond_order_rho: num 1
# $ ncat             : num 8

# Specify model in BUGS language
cat(file = 'static_categorical.txt', "
model{

  # --- Priors ---
  # First order psi
  betaB[1] <- logit(mean.psiB) # fo occupancy intercepts
  betaC[1] <- logit(mean.psiC)
  betaF[1] <- logit(mean.psiF)
  mean.psiB[1] ~ dunif(0, 1)
  mean.psiC[1] ~ dunif(0, 1)
  mean.psiF[1] ~ dunif(0, 1)
  for(fo_psi in 2:nfirst_order_psi){ # fo occupancy slopes
    betaB[fo_psi] ~ dnorm(0, 0.1)
    betaC[fo_psi] ~ dnorm(0, 0.1)
    betaF[fo_psi] ~ dnorm(0, 0.1)
  }
  # Second order psi priors
  for(so_psi in 1:nsecond_order_psi){
    betaBC[so_psi] ~ dnorm(0, 0.1)
    betaBF[so_psi] ~ dnorm(0, 0.1)
    betaCF[so_psi] ~ dnorm(0, 0.1)
  }
  # First order detection priors (rho)
  alphaB[1] <- logit(mean.pB) # fo detection intercepts
  alphaC[1] <- logit(mean.pC)
  alphaF[1] <- logit(mean.pF)
  mean.pB ~ dunif(0, 1)
  mean.pC ~ dunif(0, 1)
  mean.pF ~ dunif(0, 1)
  for(fo_rho in 2:nfirst_order_rho){ # fo detection slopes
    alphaB[fo_rho] ~ dnorm(0, 0.1)
    alphaC[fo_rho] ~ dnorm(0, 0.1)
    alphaF[fo_rho] ~ dnorm(0, 0.1)
  }
  # Second order detection priors (rho)
  # none in this model

  # --- 'Likelihood' ---
  # (1) Basic hierarchical model: states and observations
  # Latent state model
  for(i in 1:nsites) {
    z[i] ~ dcat(lsv[i, ( 1:ncat )] )
  }
  # Detection model
  for(i in 1:nsites) {
    for(j in 1:nsurveys) {
      y[i, j] ~ dcat(rdm[i, j, ( 1:ncat ) , z[i] ] )
    }
  }

  # (2) Define the latent state vector and the observation matrices
  for( i in 1:nsites ) {
    # Latent state probabilities in latent state vector (lsv)
    # Probabilities for each state
    lsv[i, 1] <- 1 #------------------------------------------| U
    lsv[i, 2] <- exp( psiB[i] ) #-----------------------------| B
    lsv[i, 3] <- exp( psiC[i] ) #-----------------------------| C
    lsv[i, 4] <- exp( psiF[i] ) #-----------------------------| F
    lsv[i, 5] <- exp( psiBC[i] ) #----------------------------| BC
    lsv[i, 6] <- exp( psiCF[i] ) #----------------------------| CF
    lsv[i, 7] <- exp( psiBF[i] ) #----------------------------| BF
    lsv[i, 8] <- exp( psiBCF[i] ) #--------------------------| BCF
    for(j in 1:nsurveys){
      # Detection matrix (OS = observed state, TS = true state)
      # rdm = rho detection matrix.
      # OS along rows, TS along columns
      # True state = U
      rdm[i, j, 1, 1] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 1] <- 0 #----------------------------------| OS = B
      rdm[i, j, 3, 1] <- 0 #----------------------------------| OS = C
      rdm[i, j, 4, 1] <- 0 #----------------------------------| OS = F
      rdm[i, j, 5, 1] <- 0 #----------------------------------| OS = BC
      rdm[i, j, 6, 1] <- 0 #----------------------------------| OS = CF
      rdm[i, j, 7, 1] <- 0 #----------------------------------| OS = BF
      rdm[i, j, 8, 1] <- 0 #----------------------------------| OS = BCF
      # True state = B
      rdm[i, j, 1, 2] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 2] <- exp( rhoB[i, j] ) #------------------| OS = B
      rdm[i, j, 3, 2] <- 0 #----------------------------------| OS = C
      rdm[i, j, 4, 2] <- 0 #----------------------------------| OS = F
      rdm[i, j, 5, 2] <- 0 #----------------------------------| OS = BC
      rdm[i, j, 6, 2] <- 0 #----------------------------------| OS = CF
      rdm[i, j, 7, 2] <- 0 #----------------------------------| OS = BF
      rdm[i, j, 8, 2] <- 0 #----------------------------------| OS = BCF
      # True state = C
      rdm[i, j, 1, 3] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 3] <- 0 #----------------------------------| OS = B
      rdm[i, j, 3, 3] <- exp( rhoC[i, j] ) #------------------| OS = C
      rdm[i, j, 4, 3] <- 0 #----------------------------------| OS = F
      rdm[i, j, 5, 3] <- 0 #----------------------------------| OS = BC
      rdm[i, j, 6, 3] <- 0 #----------------------------------| OS = CF
      rdm[i, j, 7, 3] <- 0 #----------------------------------| OS = BF
      rdm[i, j, 8, 3] <- 0 #----------------------------------| OS = BCF
      # True state = F
      rdm[i, j, 1, 4] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 4] <- 0 #----------------------------------| OS = B
      rdm[i, j, 3, 4] <- 0 #----------------------------------| OS = C
      rdm[i, j, 4, 4] <- exp( rhoF[i, j] ) #------------------| OS = F
      rdm[i, j, 5, 4] <- 0 #----------------------------------| OS = BC
      rdm[i, j, 6, 4] <- 0 #----------------------------------| OS = CF
      rdm[i, j, 7, 4] <- 0 #----------------------------------| OS = BF
      rdm[i, j, 8, 4] <- 0 #----------------------------------| OS = BCF
      # True state = BC
      rdm[i, j, 1, 5] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 5] <- exp( rhoBC[i, j] ) #-----------------| OS = B
      rdm[i, j, 3, 5] <- exp( rhoCB[i, j] ) #-----------------| OS = C
      rdm[i, j, 4, 5] <- 0 #----------------------------------| OS = F
      rdm[i, j, 5, 5] <- exp( rhoBC[i, j] + rhoCB[i, j]) #----| OS = BC
      rdm[i, j, 6, 5] <- 0 #----------------------------------| OS = CF
      rdm[i, j, 7, 5] <- 0 #----------------------------------| OS = BF
      rdm[i, j, 8, 5] <- 0 #----------------------------------| OS = BCF
      # True state = CF
      rdm[i, j, 1, 6] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 6] <- 0 #----------------------------------| OS = B
      rdm[i, j, 3, 6] <- exp( rhoCF[i, j] ) #-----------------| OS = C
      rdm[i, j, 4, 6] <- exp( rhoFC[i, j] ) #-----------------| OS = F
      rdm[i, j, 5, 6] <- 0 #----------------------------------| OS = BC
      rdm[i, j, 6, 6] <- exp( rhoCF[i, j] + rhoFC[i, j] ) #---| OS = CF
      rdm[i, j, 7, 6] <- 0 #----------------------------------| OS = BF
      rdm[i, j, 8, 6] <- 0 #----------------------------------| OS = BCF
      # True state = BF
      rdm[i, j, 1, 7] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 7] <- exp( rhoBF[i, j] ) #-----------------| OS = B
      rdm[i, j, 3, 7] <- 0 #----------------------------------| OS = C
      rdm[i, j, 4, 7] <- exp( rhoFB[i, j] ) #-----------------| OS = F
      rdm[i, j, 5, 7] <- 0 #----------------------------------| OS = BC
      rdm[i, j, 6, 7] <- 0 #----------------------------------| OS = CF
      rdm[i, j, 7, 7] <- exp( rhoBF[i, j] + rhoFB[i, j] ) #---| OS = BF
      rdm[i, j, 8, 7] <- 0 #----------------------------------| OS = BCF
      # True state = BCF
      rdm[i, j, 1, 8] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 8] <- exp( rhoBCF[i, j] ) #----------------| OS = B
      rdm[i, j, 3, 8] <- exp( rhoCBF[i, j] ) #----------------| OS = C
      rdm[i, j, 4, 8] <- exp( rhoFBC[i, j] ) #----------------| OS = F
      rdm[i, j, 5, 8] <- exp( rhoBCF[i, j] + rhoCBF[i, j] ) #-| OS = BC
      rdm[i, j, 6, 8] <- exp( rhoCBF[i, j] + rhoFBC[i, j] ) #-| OS = CF
      rdm[i, j, 7, 8] <- exp( rhoBCF[i, j] + rhoFBC[i, j] ) #-| OS = BF
      rdm[i, j, 8, 8] <- exp( rhoBCF[i, j] + rhoCBF[i, j] +
          rhoFBC[i, j] ) #------------------------------------| OS = BCF
    }

    # (3) Specify linear models for the parameters in lsv and rdm
    # Linear models for the occupancy parameters
    # ...for states B, C, and F
    psiB[i] <- inprod( betaB, psi_cov[i, ] )
    psiC[i] <- inprod( betaC, psi_cov[i, ] )
    psiF[i] <- inprod( betaF, psi_cov[i, ] )
    # ...for states BC, CF, and BF (in that order)
    psiBC[i] <- psiB[i] + psiC[i] + inprod( betaBC, psi_inxs_cov[i, ] )
    psiCF[i] <- psiC[i] + psiF[i] + inprod( betaCF, psi_inxs_cov[i, ] )
    psiBF[i] <- psiB[i] + psiF[i] + inprod( betaBF, psi_inxs_cov[i, ] )
    # ...for state BCF
    psiBCF[i] <- psiB[i] + psiC[i] + psiF[i] +
        inprod( betaBC, psi_inxs_cov[i, ] ) +
        inprod( betaCF, psi_inxs_cov[i, ] ) +
        inprod( betaBF, psi_inxs_cov[i, ] )

    # Linear models for the detection parameters
    # => Here we could specify detection interactions as well
    for(j in 1:nsurveys){
      # Baseline detection linear predictors
      # do not incorporate interactions.
      rhoB[i, j] <- inprod( alphaB, rho_cov[i, j, ] )
      rhoC[i, j] <- inprod( alphaC, rho_cov[i, j, ] )
      rhoF[i, j] <- inprod( alphaF, rho_cov[i, j, ] )
      # Asymmetric interactions between all 3 species
      rhoBC[i, j] <- rhoB[i, j]
      rhoBF[i, j] <- rhoB[i, j]
      rhoCB[i, j] <- rhoC[i, j]
      rhoCF[i, j] <- rhoC[i, j]
      rhoFB[i, j] <- rhoF[i, j]
      rhoFC[i, j] <- rhoF[i, j]
      # Asymmetric interactions when all 3 species are present
      rhoBCF[i, j] <- rhoB[i, j]
      rhoCBF[i, j] <- rhoC[i, j]
      rhoFBC[i, j] <- rhoF[i, j]
    }
  }
}
")

# get the maximum possible state across all 3 potential surveys at a site
# returns a site x species matrix
zinit <- apply(y, c(1,3), sum, na.rm = TRUE)
zinit[zinit>1] <- 1 # make binary

# convert to a category
zcat <- apply(zinit, 1, paste, collapse = '')
zcat[zcat == "000"] <- 1 # nobody there
zcat[zcat == "100"] <- 2 # only bobcat
zcat[zcat == "010"] <- 3 # only coyote
zcat[zcat == "001"] <- 4 # only fox
zcat[zcat == "110"] <- 5 # bobcat and coyote
zcat[zcat == "011"] <- 6 # coyote and fox
zcat[zcat == "101"] <- 7 # bobcat and fox
zcat[zcat == "111"] <- 8 # all three

# make numeric again
zcat <- as.numeric(zcat)

# Inits function
inits <- function() list(z = zcat)

# Parameters monitored
params <- c('betaB', 'betaC', 'betaF', 'betaBC', 'betaBF', 'betaCF',
    'alphaB', 'alphaC', 'alphaF', 'mean.psiB', 'mean.psiC', 'mean.psiF',
    'mean.pB', 'mean.pC', 'mean.pF', 'z')

# MCMC settings
# na <- 10000 ; nc <- 3 ; ni <- 50000 ; nb <- 30000 ; nt <- 20
na <- 1000 ; nc <- 3 ; ni <- 5000 ; nb <- 3000 ; nt <- 2  # ~~~ for testing, 40 mins

# Call JAGS (ART 335 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, 'static_categorical.txt', n.chains = nc,
    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
which(out1$summary[,8] > 1.1)
print(out1$summary[1:24, -c(4:6)], 3)
#              mean     sd     2.5%   97.5% Rhat n.eff overlap0     f
# betaB[1]  -1.1412 0.2712 -1.69251 -0.6280 1.01   339        0 1.000
# betaB[2]  -0.6594 0.1482 -0.96154 -0.3856 1.01   406        0 1.000
# betaC[1]  -0.5343 0.2144 -0.98863 -0.1479 1.03    74        0 0.998
# betaC[2]   0.2107 0.1129  0.00607  0.4518 1.01   315        0 0.979
# betaF[1]  -1.8587 0.2668 -2.39865 -1.3824 1.00  1296        0 1.000
# betaF[2]  -0.3743 0.1399 -0.66442 -0.1159 1.00  1084        0 0.999
# betaBC[1]  0.8798 0.4150  0.05126  1.6751 1.06    40        0 0.984
# betaBC[2] -2.8701 0.9379 -4.96531 -1.3914 1.03    81        0 1.000
# betaBF[1] -1.5007 0.4597 -2.42223 -0.6409 1.01   171        0 1.000
# betaBF[2]  1.7285 0.9461  0.46911  4.0528 1.03    68        0 0.997
# betaCF[1]  1.5942 0.4749  0.80321  2.6352 1.01   395        0 1.000
# betaCF[2]  1.9149 0.9643  0.82842  4.3680 1.04    70        0 1.000
# alphaB[1] -2.5295 0.1645 -2.85495 -2.2177 1.00  2034        0 1.000
# alphaB[2]  1.8092 0.1685  1.48266  2.1413 1.00  1414        0 1.000
# alphaC[1] -2.0232 0.1043 -2.22755 -1.8163 1.01   349        0 1.000
# alphaC[2]  2.1582 0.1195  1.92025  2.3862 1.00  3000        0 1.000
# alphaF[1] -1.6145 0.1739 -1.94695 -1.2800 1.01   275        0 1.000
# alphaF[2]  1.8453 0.2080  1.42525  2.2442 1.00  1203        0 1.000
# mean.psiB  0.2455 0.0495  0.15545  0.3480 1.01   373        0 1.000
# mean.psiC  0.3709 0.0490  0.27118  0.4631 1.03    75        0 1.000
# mean.psiF  0.1378 0.0306  0.08328  0.2006 1.00  1235        0 1.000
# mean.pB    0.0746 0.0113  0.05443  0.0982 1.00  2056        0 1.000
# mean.pC    0.1172 0.0108  0.09730  0.1399 1.01   356        0 1.000
# mean.pF    0.1674 0.0242  0.12489  0.2176 1.01   281        0 1.000

# ~~~ The plots in the book do not show CRIs. ~~~
# ~~~ For alternative code with CRIs, see file AHM2_08.2.3_CRIplots.R ~~~
# ~~~ save the workspace so far ~~~
save(list=ls(), file="AHM2_08.2.3.Rdata")


# Grab posterior means of parameters
(tmp <- out1$mean[1:6])

# Create prediction covariates for occupancy as a function of Disturbance
# (at average Housing density)
(r <- range(data$sitecovs[,'Dist_5km']))
(mnx <- mean(data$sitecovs[,'Dist_5km']))
(sdx <- sd(data$sitecovs[,'Dist_5km']))
Dist.pred.orig <- seq(r[1], r[2], length.out = 1000) # Create new data
Dist.pred <- (Dist.pred.orig-mnx) / sdx      # ... and scale as the real data
HDens.pred <- rep(0, 1000)

# Create model matrices for the prediction covariates
head(psi_cov1 <- cbind(1, Dist.pred))    # We vary this covariate ...
head(psi_cov2 <- cbind(1, HDens.pred))   # .. while this stays constant

# Assemble the linear predictors
# ...for states B, C, and F
psiB <- as.numeric( psi_cov1 %*% tmp$betaB )
psiC <- as.numeric( psi_cov1 %*% tmp$betaC )
psiF <- as.numeric( psi_cov1 %*% tmp$betaF )

# ...for states BC, CF, and BF (in that order)
psiBC <- psiB + psiC + as.numeric( psi_cov2 %*% tmp$betaBC )
psiCF <- psiC + psiF + as.numeric( psi_cov2 %*% tmp$betaCF )
psiBF <- psiB + psiF + as.numeric( psi_cov2 %*% tmp$betaBF )

# ...for state ABC
psiBCF <- psiB + psiC + psiF + psiBC + psiCF + psiBF

# Assemble the state probability vector for every DIST covariate value
lsv <- lsp <- matrix(NA, ncol = ncat, nrow = 1000 )
colnames(lsv) <- colnames(lsp) <- c('U', 'B', 'C', 'F', 'BC', 'CF', 'BF', 'BCF')

# Compute latent state vector (lsv) and latent state probabilities (lsp)
lsv[ , 1] <- 1              # none of the three species
lsv[ , 2] <- exp( psiB )    # only B
lsv[ , 3] <- exp( psiC )    # only C
lsv[ , 4] <- exp( psiF )    # only F
lsv[ , 5] <- exp( psiBC)    # both B and C
lsv[ , 6] <- exp( psiCF)    # both C and F
lsv[ , 7] <- exp( psiBF)    # both B and F
lsv[ , 8] <- exp( psiBCF )  # all three species

# Probability of each state
for(i in 1:nrow(lsv)){
  lsp[i, ] <- lsv[i, ( 1:ncat )] / sum(lsv[i, ( 1:ncat )])
}

# Plot the lsp as a function of covariate Dist (Fig. 8.4)
matplot(Dist.pred.orig, lsp, type = 'l', lty = 1, lwd = 5, col = 1:8,
    frame = FALSE, xlab = 'Disturbance', ylab = 'Probability of state',
    ylim = c(0, 0.65), las = 1)
legend('topleft', lwd = 3, lty = 1, col = 1:8, colnames(lsp), bty = 'n',
    cex = 1.2, ncol = 2)

# Compute and plot marginal occupancy for all three species (not shown)
bobcat.marginal <- rowSums(lsp[,c(2,5,7,8)])
coyote.marginal <- rowSums(lsp[,c(3,5,6,8)])
fox.marginal <- rowSums(lsp[,c(4,6,7,8)])
matplot(Dist.pred.orig, cbind(bobcat.marginal, coyote.marginal, fox.marginal),
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE, xlab = 'Disturbance',
    ylab = 'Marginal occupancy probability', ylim = c(0, 0.7))
legend('top', c('Bobcat', 'Coyote', 'Red Fox'), lwd = 3, lty = 1:3, col = 1:3,
    horiz = TRUE, bty = 'n')

# Compute conditional occupancy for all species
B.alone <- lsp[, 'B'] / rowSums(lsp[, c('U', 'B')])
B.given.coyote <- rowSums(lsp[, c('BC', 'BCF')]) /
    rowSums(lsp[, c('C', 'BC', 'CF', 'BCF')])
B.given.fox <- rowSums(lsp[, c('BF', 'BCF')]) /
    rowSums(lsp[, c('F', 'CF', 'BF', 'BCF')])
C.alone <- lsp[, 'C'] / rowSums(lsp[, c('U', 'C')])
C.given.bobcat <- rowSums(lsp[, c('BC', 'BCF')]) /
    rowSums(lsp[, c('B', 'BC', 'BF', 'BCF')])
C.given.fox <- rowSums(lsp[, c('CF', 'BCF')]) /
    rowSums(lsp[, c('F', 'CF', 'BF', 'BCF')])
F.alone <- lsp[, 'F'] / rowSums(lsp[, c('U', 'F')])
F.given.bobcat <- rowSums(lsp[, c('BF', 'BCF')]) /
    rowSums(lsp[, c('B', 'BC', 'BF', 'BCF')])
F.given.coyote <- rowSums(lsp[, c('CF', 'BCF')]) /
    rowSums(lsp[, c('C', 'BC', 'CF', 'BCF')])

# Plot them (Fig. 8.5) -> see website
# ~~~ extra code for figure 8.5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
op <- par(mfrow = c(1, 3), mar = c(6,6,6,3), cex.lab = 2, cex.axis = 2,
    cex.main = 2)
matplot(Dist.pred.orig, cbind(B.alone, B.given.coyote, B.given.fox),
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Disturbance', ylab = 'Occupancy probability',
    ylim = c(0, 0.7), main = 'Bobcat')
legend('topright', c('alone', 'given Coyote presence', 'given Fox presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
matplot(Dist.pred.orig, cbind(C.alone, C.given.bobcat, C.given.fox),
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Disturbance', ylab = 'Occupancy probability',
    ylim = c(0, 1), main = 'Coyote')
legend('bottomright', c('alone', 'given Bobcat presence', 'given Fox presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
matplot(Dist.pred.orig, cbind(F.alone, F.given.bobcat, F.given.coyote),
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Disturbance', ylab = 'Occupancy probability',
    ylim = c(0, 0.7), main = 'Red Fox')
legend('topright', c('alone', 'given Bobcat presence', 'given Coyote presence'),
    lwd = 3, lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ extra code for figures 8.6 and 8.7 ~~~~~~
# Create prediction covariates for occupancy as a function of Housing density
#   (at average Disturbance) - figure 8.6
(r <- range(data$sitecovs[,'HDens_5km']))
mnx <- mean(data$sitecovs[,'HDens_5km'])
sdx <- sd(data$sitecovs[,'HDens_5km'])
HDens.pred.orig <- seq(r[1], 120, length.out = 1000)
HDens.pred <- (HDens.pred.orig-mnx) / sdx
Dist.pred <- rep(0, 1000)

# Create new model matrices for the prediction covariates
head(psi_cov1 <- cbind(1, Dist.pred))  # Now we keep this constant
head(psi_cov2 <- cbind(1, HDens.pred)) # But we vary this

# Assemble the linear predictors
# ...for states B, C, and F
psiB  <- as.numeric( psi_cov1 %*% tmp$betaB )
psiC  <- as.numeric( psi_cov1 %*% tmp$betaC )
psiF  <- as.numeric( psi_cov1 %*% tmp$betaF )
# ...for states BC, CF, and BF (in that order)
psiBC <-  psiB + psiC +  as.numeric( psi_cov2 %*% tmp$betaBC )
psiCF <-  psiC + psiF +  as.numeric( psi_cov2 %*% tmp$betaCF )
psiBF <-  psiB + psiF +  as.numeric( psi_cov2 %*% tmp$betaBF )
# ...for state ABC
psiBCF <-  psiB + psiC + psiF + psiBC + psiCF + psiBF

# Assemble the state probability vector for every DIST covariate value
lsv <- lsp <-  matrix(NA, ncol = ncat, nrow = 1000 )
colnames(lsv) <- colnames(lsp) <- c('U', 'B', 'C', 'F', 'BC', 'CF', 'BF', 'BCF')
# Compute latent state vector (lsv) and latent state probabilities (lsp)
lsv[ , 1] <- 1                     # none of the three species
lsv[ , 2] <- exp( psiB )           # only B
lsv[ , 3] <- exp( psiC )           # only C
lsv[ , 4] <- exp( psiF )           # only F
lsv[ , 5] <- exp( psiBC)           # both B and C
lsv[ , 6] <- exp( psiCF)           # both C and F
lsv[ , 7] <- exp( psiBF)           # both B and F
lsv[ , 8] <- exp( psiBCF )         # all three species
# Probability of each state
for(i in 1:nrow(lsv)){
  lsp[i, ] <- lsv[i, ( 1:ncat )] / sum(lsv[i, ( 1:ncat )])
}

# Compute and plot marginal occupancy for all three species
bobcat.marginal <- rowSums(lsp[,c(2,5,7,8)])
coyote.marginal <- rowSums(lsp[,c(3,5,6,8)])
fox.marginal <- rowSums(lsp[,c(4,6,7,8)])

# Plot the lsp and marginal occupancy of each species as a function of covariate HDens while keeping Dist constant (Fig. 8.6)
op <- par(mfrow = c(1, 2))
matplot(HDens.pred.orig, lsp, type = 'l', lty = 1, lwd = 5, col = 1:8,
    frame = FALSE, xlab = 'Housing density', ylab = 'Probability',
    ylim = c(0, 1.05), las = 1, main = 'Probability of each state')
legend(70, 0.6, lwd = 3, lty = 1, col = 1:8, colnames(lsp), bty = 'n',
    cex = 1.1, ncol = 2)
matplot(HDens.pred.orig, cbind(bobcat.marginal, coyote.marginal, fox.marginal),
    type = 'l', lty = 1:3, lwd = 5, col = 1:3, frame = FALSE,
    xlab = 'Housing density', ylab = 'Probability', ylim = c(0, 1),
    main = 'Marginal occupancy probability')
legend(90, 0.6, c('Bobcat', 'Coyote', 'Red Fox'), lwd = 3, lty = 1:3,
    col = 1:3, bty = 'n')
par(op)

# Compute conditional occupancy for all species, now with Housing density
B.alone <- lsp[, 'B'] / rowSums(lsp[, c('U', 'B')])
B.given.coyote <- rowSums(lsp[, c('BC', 'BCF')]) /
    rowSums(lsp[, c('C', 'BC', 'CF', 'BCF')])
B.given.fox <- rowSums(lsp[, c('BF', 'BCF')]) /
    rowSums(lsp[, c('F', 'CF', 'BF', 'BCF')])
C.alone <- lsp[, 'C'] / rowSums(lsp[, c('U', 'C')])
C.given.bobcat <- rowSums(lsp[, c('BC', 'BCF')]) /
    rowSums(lsp[, c('B', 'BC', 'BF', 'BCF')])
C.given.fox <- rowSums(lsp[, c('CF', 'BCF')]) /
    rowSums(lsp[, c('F', 'CF', 'BF', 'BCF')])
F.alone <- lsp[, 'F'] / rowSums(lsp[, c('U', 'F')])
F.given.bobcat <- rowSums(lsp[, c('BF', 'BCF')]) /
    rowSums(lsp[, c('B', 'BC', 'BF', 'BCF')])
F.given.coyote <- rowSums(lsp[, c('CF', 'BCF')]) /
    rowSums(lsp[, c('C', 'BC', 'CF', 'BCF')])

# Plot them (Fig. 8.7)
op <- par(mfrow = c(1, 3), mar = c(6,6,6,3), cex.lab = 2, cex.axis = 2,
    cex.main = 2)
matplot(HDens.pred.orig, cbind(B.alone, B.given.coyote, B.given.fox),
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Housing density', ylab = 'Occupancy probability',
    ylim = c(0, 1), main = 'Bobcat')
legend(60, 0.7, c('alone', 'given Coyote', 'given Fox'), lwd = 3,
    lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)

matplot(HDens.pred.orig, cbind(C.alone, C.given.bobcat, C.given.fox),
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Housing density', ylab = 'Occupancy probability',
    ylim = c(0, 1), main = 'Coyote')
legend(60, 0.7, c('alone', 'given Bobcat', 'given Fox'), lwd = 3,
    lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)

matplot(HDens.pred.orig, cbind(F.alone, F.given.bobcat, F.given.coyote),
    type = 'l', lty = 1:3, lwd = 3, col = 1:3, frame = FALSE,
    xlab = 'Housing density', ylab = 'Occupancy probability',
    ylim = c(0, 1), main = 'Red Fox')
legend(60, 0.45, c('alone', 'given Bobcat', 'given Coyote'), lwd = 3,
    lty = 1:3, col = 1:3, bty = 'n', cex = 1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~ extra code to get detection probabilities ~~~~~~~~~~~~~
# Look at predictions of detection model
# Set up structures to hold results
nsims <- out1$mcmc.info$n.samples
Trail.sims <- array(NA, dim = c(3, 2, nsims))
dimnames(Trail.sims) <- list(c('Bobcat', 'Coyote', 'Red Fox'),
  c('off trail', 'on trail'), NULL)
tmp <- out1$sims.list

# Compute the posteriors for predictions
# For Bobcat
Trail.sims[1,1,] <- plogis(tmp$alphaB[,1])
Trail.sims[1,2,] <- plogis(tmp$alphaB[,2])
# For Coyote
Trail.sims[2,1,] <- plogis(tmp$alphaC[,1])
Trail.sims[2,2,] <- plogis(tmp$alphaC[,2])
# For Red Fox
Trail.sims[3,1,] <- plogis(tmp$alphaF[,1])
Trail.sims[3,2,] <- plogis(tmp$alphaF[,2])

# Get posterior mean and sd
round(pm <- apply(Trail.sims, c(1,2), mean), 2)
round(CRI <- apply(Trail.sims, c(1,2),
    quantile, probs=c(0.025, 0.975)), 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Posterior means of p
#         off trail on trail
# Bobcat       0.07     0.86
# Coyote       0.12     0.90
# Red Fox      0.17     0.86

# 95% CRIs
# , , off trail
#       Bobcat Coyote Red Fox
# 2.5%    0.05   0.10    0.12
# 97.5%   0.10   0.14    0.22
# , , on trail
#       Bobcat Coyote Red Fox
# 2.5%    0.81   0.87    0.81
# 97.5%   0.89   0.92    0.90

# Compute and plot most likely states for a sample of the sites (figure 8.8)
nsims <- out1$mcmc.info$n.samples
state.prob <- most.likely <- array(NA, dim = c(1437, 8))
colnames(state.prob) <- colnames(most.likely) <- c("U", "B", "C", "F",
    "BC", "CF", "BF", "BCF") # label
for(i in 1:1437){
  state.prob[i,] <- tabulate(out1$sims.list$z[,i], nbins = 8)/nsims
  tmp <- tabulate(out1$sims.list$z[,i], nbins = 8)
  most.likely[i,] <- floor(tmp / max(tmp))
}
head(state.prob) # Gives probability that a site is in a certain state
head(most.likely) # Just picks the most likely state
op <- par(mar = c(5,5,6,2), cex.main = 2) # Plot 8.8
image(x = 1:500, y = 1:8, z = most.likely[1:500,], col = c('white', 'black'),
    xlab = 'Sites 1-500', ylab = 'State', axes = FALSE)
axis(1)
axis(2, at = 1:8, labels = c("U", "B", "C", "F", "BC", "CF", "BF", "BCF"), las = 1)
par(op)
