#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 9. Advanced Hierarchical Distance Sampling
# =========================================================================

# Approximate execution time for this code: 25 mins
# Run time with the full number of iterations: 3.7 hrs

library(AHMbook)
library(unmarked)
library(raster)
library(plotrix)

# ~~~~~~~ changes to RNG defaults ~~~~~~~~~~~~~~~~~~~~~~~~
# Use the old default random number generator to get the printed numbers
RNGversion("3.2.0")

# 9.8 Spatial Distance Sampling: Modelling within-unit variation in density
# =========================================================================

# 9.8.1 Distance sampling with location of encounter
# ------------------------------------------------------------------------
# Simulate a data set and harvest the output
set.seed(1234)
str(tmp <- sim.pdata(N=200,sigma=1,keep.all=FALSE,B=3))

# Harvest some data objects
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
nind <- length(d)

# Data augmentation
M <- 400                       # max of 400 individuals
nz <- M-nind                   # augment by nz individuals
y <- c(rep(1,nind), rep(0,nz)) # augmented data augmentation variable
u <- cbind(u1,u2)               # Locations of individuals detected
u <- rbind(u, matrix(NA, nrow=nz, ncol=2))


# Bundle and summarize the data set
str(data <- list (B=B, nind=nind, u=u, y=y, nz=nz))


# Write out the BUGS model file
cat("
model{

  # Priors
  sigma ~ dunif(0,10)
  psi ~ dunif(0,1)

  # Categorical observation model
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(0, 2*B)  # Here is the uniformity assumption made explicit
    u[i,2] ~ dunif(0, 2*B)
    # Compute distance as a derived quantity
    d[i] <- pow( pow( u[i,1]-B,2) + pow(u[i,2]-B,2), 0.5) # Pythagoras
    p[i] <- exp(-d[i]*d[i] / (2*sigma*sigma))
    mu[i] <- p[i] * z[i]
    y[i] ~ dbern(mu[i])
  }
  # Other derived quantity
  N <- sum(z[])
  D <- N / (B*B)
}
",fill=TRUE,file="model1.txt")


# Load libraries and specify MCMC settings
library("R2WinBUGS")   ;   library(jagsUI)
ni <- 22000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3

# Inits and parameters
inits <- function(){
  list(sigma=runif(1,1,10), psi=runif(1),z = c(rep(1,nind),rep(0,nz)) ) }
params <- c("sigma", "N", "psi")

# Execute jags and summarize the posterior distributions
# out1 <- jags (data, inits, parameters, "model1.txt", n.thin=nthin,  # ~~~ oops!
out1 <- jags (data, inits, params, "model1.txt", n.thin=nthin,
   # n.chains=nc, n.burnin=nb,n.iter=ni, parallel = FALSE)
   n.chains=nc, n.burnin=nb,n.iter=ni, parallel = TRUE)  # ~~~~~ for autotesting
# par(mfrow = c(2,2))  #  ~~~ replace with 'layout' argument
traceplot(out1, layout=c(2,2))
print(out1, 2)

# 9.8.2 The line transect case (no code)

# 9.8.3 Modelling spatial covariates
# ------------------------------------------------------------------------
# Simulator function for spatial distance sampling data
sim.spatialDS <-
function(N=1000, beta = 1, sigma=1, keep.all=FALSE, B=B, model="halfnorm"){
  # Function simulates coordinates of individuals on a square
  # Square is [0,2B] x [0,2B], with a count location on the point (B, B)
  #   N: total population size in the square
  #   beta: coefficient of SOEMTHING on spatial covariate x
  #   sigma: scale of half-normal detection function
  #   B: circle radius
  #   keep.all: return the data for y=0 individuals or not
  library(raster)      # Load required packages
  library(plotrix)

  # Create coordinates for 30 x 30 grid
  delta <- (2*B-0)/30                # '2D bin width'
  grx <- seq(delta/2, 2*B - delta/2, delta) # mid-point coordinates
  gr <- expand.grid(grx,grx)         # Create grid coordinates

  # Create spatially correlated covariate x and plot it
  V <- exp(-e2dist(gr,gr)/1)
  x <- t(chol(V))%*%rnorm(900)
  op <- par(mar=c(3,3,3,6))
  image(rasterFromXYZ(cbind(gr,x)), col=topo.colors(10))
  draw.circle(3, 3, B)
  points(3, 3, pch="+", cex=3)
  image_scale(x, col=topo.colors(10))  ## function renamed

  # Simulate point locations as function of habitat covariate x
  probs <- exp(beta*x)/sum(exp(beta*x)) # probability of point in pixel (sum = 1)
  pixel.id <- sample(1:900, N, replace=TRUE, prob=probs)
  # could simulate randomly within the pixel but it won't matter so place centrally
  u1 <- gr[pixel.id,1]
  u2 <- gr[pixel.id,2]
  points(u1, u2, pch=20, col='black', cex = 0.8)  # plot points
  title("This is so cool !")         # express your appreciation of all this

  d <- sqrt((u1 - B)^2 + (u2-B)^2)   # distance to center point of square
  #plot(u1, u2, pch = 1, main = "Point transect")
  N.real <- sum(d<= B)               # Population size inside of count circle

  # Can only count individuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them)
  # p <- ifelse(d< B, 1, 0) * exp(-d*d/(2*(sigma^2)))
  # We do away with the circle constraint here.
  if(model=="hazard")
     p <- 1-exp(-exp(-d*d/(2*sigma*sigma)))
  if(model=="halfnorm")
     p <- exp(-d*d/(2*sigma*sigma))
  # Now we decide whether each individual is detected or not
  y <- rbinom(N, 1, p)                                           # detected or not
  points(u1[d<= B], u2[d<= B], pch = 16, col = "black", cex = 1) # not detected
  points(u1[y==1], u2[y==1], pch = 16, col = "red", cex = 1)     # detected
  par(op)

  # Put all of the data in a matrix
  if(!keep.all){
     u1 <- u1[y==1]
     u2 <- u2[y==1]
     d <- d[y==1]   }
  # Output
  return(list(model=model, N=N, beta=beta, B=B, u1=u1, u2=u2, d=d, y=y,
      N.real=N.real, Habitat=x, grid=gr))
}

# Generate one data set and harvest the output
set.seed(1234)
str(tmp <- sim.spatialDS(N=200, beta=1, sigma=1.5, keep.all=FALSE, B=3)) # Fig. 9-7


# Harvest data
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
Habitat <- as.vector(tmp$Habitat)
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
nind <- length(d)
G <- nrow(Habgrid)

# Do data augmentation, including for pixel ID
M <- 400
nz <- M-nind
pixel <- rep(NA, M)   # We use discrete "pixel ID" here instead of "s"
y <- c(rep(1,nind), rep(0,nz))

# Pick some starting values and figure out the pixel of each observation
s <- cbind(u1,u2)
s <- rbind(s, matrix(NA,nrow=nz,ncol=2))
D <- e2dist(s[1:nind,], Habgrid)
for(i in 1:nind){
   pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Bundle and summarize the data for BUGS
str(data <- list (B=B, nind=nind, y=y, nz=nz, Habitat=Habitat,
  # Habgrid=Habgrid, G=G, pixel=pixel))
  Habgrid=as.matrix(Habgrid), G=G, pixel=pixel)) # ~~~~ jagsUI no longer converts data internally


# Write BUGS model
cat("
model{

  # Prior distributions
  sigma ~ dunif(0,10)
  psi ~ dunif(0,1)
  beta ~ dnorm(0,0.01)

  for(g in 1:G){   # g is the pixel index, there are G total pixels
    probs.num[g] <- exp(beta*Habitat[g])
    probs[g] <- probs.num[g]/sum(probs.num[])
  }

  # Models for DA variables and location (pixel)
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    pixel[i] ~ dcat(probs[])
    s[i,1:2] <- Habgrid[pixel[i],]   # location = derived quantity
    # compute distance = derived quantity
    d[i] <- pow(   pow( s[i,1]-B,2) + pow(s[i,2]-B,2), 0.5)
    p[i] <- exp(-d[i]*d[i]/(2*sigma*sigma))  # Half-normal detetion function
    mu[i]<- p[i]*z[i]
    y[i] ~ dbern(mu[i])                      # Observation model
  }
  # Derived parameters
  N <- sum(z[])                      # N is a derived parameter
  D <- N/9                           # area = 9 ha
}
",fill=TRUE, file="spatialDS.txt")


# Load libraries and specify MCMC settings
library("R2WinBUGS")
library(jagsUI)
# ni <- 12000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3  # 2 hrs
ni <- 1200   ;   nb <- 200   ;   nthin <- 1   ;   nc <- 3  # ~~~~~~~~ for testing

# Create inits and define parameters to monitor
inits <- function(){  list (sigma=runif(1,1,10),psi=runif(1),
       z = c(rep(1,nind), rep(0,nz)) ) }
params <- c("sigma", "N", "psi", "beta", "D")

# Run JAGS, check convergence and summarize posteriors
out2 <- jags (data, inits, params, "spatialDS.txt", n.thin=nthin,
   n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
# par(mfrow = c(2,3))  #  ~~~ replace with 'layout' argument
traceplot(out2, layout=c(2,3))
print(out2, 2)

# Add pixel in order to make a density map
params <- c("sigma", "N", "psi", "beta", "D", "pixel")

# Run JAGS, check convergence and summarize posteriors
out2 <- jags (data, inits, params, "spatialDS.txt", n.thin=nthin,
   # n.chains=nc, n.burnin=nb, n.iter=ni, parallel = FALSE)
   n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)  # ~~~~ speeds up testing

# Plot density maps
library(raster)
op <- par(mfrow=c(1,2))
pixel <- out2$sims.list$pixel
post <- table(pixel)/nrow(pixel)   # Average number of locations in each pixel
prior.mean <- mean(out2$sims.list$beta)*as.vector(Habitat)
prior.mean <- mean(out2$sims.list$psi)*M*exp(prior.mean)/sum(exp(prior.mean))
plot(rast.data <- rasterFromXYZ(cbind(Habgrid,prior.mean)), axes=FALSE,
       col=topo.colors(10) )
title("Prior mean density (estimated)")
plot(rast.post <- rasterFromXYZ(cbind(Habgrid,as.vector(post))),axes=FALSE,
       col=topo.colors(10) )
title("Posterior mean density")
par(op)

# 9.8.4 Spatial HDS models in unmarked using the pcount function
# ------------------------------------------------------------------------
## see errata; here we use the "logit" detection function; the results will be different
##  but not 'nonsense'.
# ~~~~ Use the sim.spatialDS function in the AHMbook package for this. ~~~~~~~
rm(sim.spatialDS)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate a data set, N = 600 for the population size
set.seed(1234)
tmp <-sim.spatialDS(N=600, sigma=1.5, keep.all=FALSE, B=3, model= "logit")

# Harvest stuff
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
Habitat <- as.vector(tmp$Habitat)
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
nind <- length(d)
G <- nrow(Habgrid)

# Find which pixel each observation belongs to
s <- cbind(u1,u2)
D <- e2dist(s[1:nind,], Habgrid)
pixel <-rep(NA,nind)
for(i in 1:nind){
   pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Create a vector of counts in each pixel and pad it with zeros
pixel.count <- rep(0, G)
names(pixel.count) <- 1:G
pixel.count[names(table(pixel))] <- table(pixel)
# Create a covariate: distance between observer and pixel center
dist <- sqrt( (Habgrid[,1]-3)^2 + (Habgrid[,2]-3)^2  )
# Construct an unmarkedFrame
umf <- unmarkedFramePCount(y=matrix(pixel.count,ncol=1),
   siteCovs=data.frame(dist=dist,Habitat=Habitat))
summary(umf)


# Fit an N-mixture model with no intercept and distance squared using
#   the function pcount.spHDS now in package unmarked

(fm1 <- pcount.spHDS(~ -1 + I(dist^2) ~ Habitat, umf, K=20))

lam <- exp( coef(fm1)[1] + coef(fm1)[2]*Habitat )
pred <- predict(fm1, type='state')
sum(lam)
sum(pred[,1])  # Same


# 9.8.5 Hierarchical spatial distance sampling
# ------------------------------------------------------------------------
sim.spatialHDS(lam0 = 4, sigma = 1.5, B = 3, nsites = 100)
 # lam0 = expected population size per site
 # nsites = number of point count locations
 # B = count radius. Function simulates coordinates of individuals on a square
 #       [0,2*B] x[0,2*B], with a count location on the point (B,B)
 # sigma = scale of half-normal detection function

library(raster)

# Simulate a data set and harvest the output
set.seed(1234)
str(tmp <-sim.spatialHDS(lam0 = 3, sigma = 1.5, B = 3, nsites = 100))

# Process the simulated data set
data <- tmp$data
# To make it a ‘real’ data set:
data <- data[!is.na(data[,2]),] # Get rid of the 0 sites
data <- data[data[,"y"]==1,]    # Only keep detected individuals

# Now zero-pad the observed counts
nsites <- tmp$nsites
nobs <- rep(0, nsites)
names(nobs) <- 1:nsites
nobs[names(table(data[,1]))] <- table(data[,1])

# Extract data elements that we need
site <- data[,"site"]
s <- data[,c("u1","u2")]
B <- tmp$B
Habitat <- (tmp$Habitat) # Raster values
Habgrid <- tmp$grid   # Raster coordinates
nind <- nrow(data)
G <- nrow(Habgrid)

# We have to convert observed locations to pixels
pixel <- rep(NA,nrow(data))
D <- e2dist(s[1:nind,], Habgrid)
for(i in 1:nind){
   pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Do data augmentation of data from each site “S-fold” DA. Three objects need
# to have DA applied to them: Ymat, umat, pixmat
Msite <- 2*max(nobs)  # Perhaps use a larger value
Ymat <- matrix(0,nrow=Msite,ncol=nsites)
umat <- array(NA, dim=c(Msite,nsites,2))
pixmat <- matrix(NA,nrow=Msite,ncol=nsites)
for(i in 1:nsites){
  if(nobs[i]==0) next
    Ymat[1:nobs[i],i]<- data[data[,1]==i,"y"]
    umat[1:nobs[i],i,1:2]<- data[data[,1]==i,c("u1","u2")]
    pixmat[1:nobs[i],i]<- pixel[data[,1]==i]
}

# Bundle the data for BUGS
str(data <- list (y=Ymat, pixel=pixmat, Habitat=Habitat,
   # Habgrid=Habgrid, G = G,
   Habgrid=as.matrix(Habgrid), G = G, # ~~~ jagsUI no longer converts input internally
   nsites=nsites, M = Msite, B = B))

# Write out the BUGS model file
cat("
model{

  # Prior distributions
  sigma ~ dunif(0,10)
  beta1 ~ dnorm(0,0.01)
  beta0 ~ dnorm(0,0.01)
  lam0 <- exp(beta0)*G           # Baseline lambda in terms of E(N) per sample unit

  # For each site, construct the DA parameter as a function of lambda
  for(s in 1:nsites){
    lamT[s] <- sum(lambda[,s])   # total abundance at a site
    psi[s] <- lamT[s]/M
      for(g in 1:G){             # g is the pixel index, there are G total pixels
      lambda[g,s] <- exp(beta0 + beta1*Habitat[g,s])
      probs[g,s] <- lambda[g,s]/sum(lambda[,s])
    }
  }

  # DA variables and spatial location variables:
  for(s in 1:nsites){
    for(i in 1:M){
      z[i,s] ~ dbern(psi[s])
      pixel[i,s] ~ dcat(probs[,s])
      u[i,s,1:2] <- Habgrid[pixel[i,s],]   # location = derived quantity
      # distance = derived quantity
      d[i,s] <- pow(   pow( u[i,s,1]-B,2) + pow(u[i,s,2]-B,2), 0.5)
      p[i,s] <- exp(-d[i,s]*d[i,s]/(2*sigma*sigma))  # Half-normal model
      mu[i,s] <- p[i,s]*z[i,s]
      y[i,s] ~ dbern(mu[i,s])    # Observation model
    }
    # Derived parameters
    N[s]<- sum(z[,s])            # Site specific abundance
  }
  Ntotal <- sum(N[])             # Total across all sites
  D <- Ntotal/(9*nsites)         # Density: point area = 9 ha
}
",fill=TRUE, file="spatialHDS.txt")


# Inits and parameters saved
zst <- Ymat
inits <- function(){ list (sigma=1.5, z = zst, beta0 = -5, beta1=1 ) }
params <- c("sigma", "Ntotal", "beta1", "beta0", "D", "lam0")
# ~~~~~ to do the plots below, you need to add "N" to params ~~~~~~~~~~~~~~~~~
params <- c("sigma", "Ntotal", "beta1", "beta0", "D", "lam0", "N")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MCMC settings
# ni <- 12000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3
ni <- 1200   ;   nb <- 200   ;   nthin <- 1   ;   nc <- 3  # ~~~ for testing

# Call JAGS, check convergence and summarize the results
out3 <- jags(data, inits, params, "spatialHDS.txt", n.thin=nthin,
    n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
# par(mfrow = c(2,3))  #  ~~~ replace with 'layout' argument
traceplot(out3, layout=c(2,3))
print(out3, 3)

Ntrue <- tmp$N
Nhat <- out3$summary[7:106,1]
plot(Ntrue, Nhat, xlab="Local population size, N, for each site",
    ylab="Posterior mean N for each site", pch=20)
abline(0, 1, lwd=2)

# ~~~~~ these were long runs, so maybe save ~~~~~~~~~~
save(out1, out2, out3, file="AHM1_09.08_JAGSoutput.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
