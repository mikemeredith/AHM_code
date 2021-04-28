#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 9. Advanced Hierarchical Distance Sampling
# =========================================================================

# 9.2 Distance sampling with clusters, groups, or other individual covariates
# ------------------------------------------------------------------------


# 9.2.1 Simulating HDS data with group size
# ------------------------------------------------------------------------
# Function to simulate data under HDS protocol with groups
simHDSg <- function(type = "line", nsites = 100, lambda.group = 0.75,
    alpha0 = 0, alpha1 = 0.5, beta0 = 1,    beta1 = 0.5, B = 4,
    discard0 = TRUE){
  #
  # Function simulates hierarchical distance sampling (HDS) data for groups under
  #   either a line (type = "line") or a point (type = "point") transect protocol
  #   and using a half-normal detection function (Buckland et al. 2001).
  #   Other function arguments:
  #     nsites: Number of sites (spatial replication)
  #     lambda.group: Poisson mean of group size
  #     alpha0, alpha1: intercept and slope of log-linear model relating sigma of
  #        half-normal detection function to group size
  #     beta0, beta1: intercept and slope of log-linear model relating the Poisson
  #        mean of the number of groups per unit area to habitat
  #     B: strip half width
  #
  # Get covariates
  habitat <- rnorm(nsites)              # Simulated covariate

  # Simulate abundance model for groups (Poisson GLM for N)
  lambda <- exp(beta0 + beta1*habitat)  # Density of groups per "square"
  N <- rpois(nsites, lambda)            # site-specific number of groups
  N.true <- N                           # for point: inside of B

  # Simulate observation model
  data <- groupsize <- NULL

  for(i in 1:nsites){
    if(N[i]==0){
      data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # save site, y=1, u, v, d
      next
  }

  if(type=="line"){
    # Simulation of distances, uniformly, for each individual in the population
    d <- runif(N[i], 0, B)
    gs <- rpois(N[i],lambda.group) +1  # Observable group sizes >= 1
    groupsize<-c(groupsize,gs)
    sigma.vec <- exp(alpha0 + alpha1*(gs-1))  # Subtract 1 for interpretation
    # Detection probability for each group
    p <- exp(-d*d/(2*(sigma.vec^2)))
    # Determine if individuals are captured or not
    y <- rbinom(N[i], 1, p)
    u1 <- u2 <- rep(NA,N[i])
    # Subset to "captured" individuals only
    d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
  }

  if(type=="point"){
    # Simulation of data on a circle of radius B (algorithm of Wallin)
    angle <- runif(N[i], 0, 360)
    r2 <- runif(N[i], 0, 1)
    r <-  B*sqrt(r2)
    u1 <-  r*cos(angle)  + B
    u2 <-  r*sin(angle)  + B

    d <- sqrt((u1 - B)^2 + (u2-B)^2)
    N.true[i] <- sum(d<= B)    # Population size inside of count circle, should be N[i] here.
    gs <- rpois(N[i], lambda.group) + 1
    groupsize <-c(groupsize,gs)
    sigma.vec <- exp(alpha0 + alpha1*(gs-1))
    # For counting individuals on a circle so we truncate p here
    p <- ifelse(d<(B), 1, 0)*exp(-d*d/(2*(sigma.vec^2)))
    y <- rbinom(N[i], 1, p)
    # Subset to "captured" individuals only
    d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
  }
   # Now compile things into a matrix and insert NA if no individuals were
   # captured at site i. Coordinates (u,v) are preserved.
   if(sum(y) > 0)
     data <- rbind(data,cbind(rep(i, sum(y)), y, u1, u2, d, gs))
   else
     data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # make a row of missing data
   }
  # Subset to sites at which individuals were captured. You may or may not
  #  do this depending on how the model is formulated so be careful.
  if(discard0)
      data <- data[!is.na(data[,2]),]

  # Visualisation
  if(type=="line"){       # For line transect
    op <- par(mfrow = c(1, 3))
    hist(data[,"d"], col = "lightblue", breaks = 20, main =
      "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = FALSE)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency",
        frame = FALSE)
    par(op)
  }

  if(type=="point"){       # For point transect
    op <- par(mfrow = c(2,2))
    plot(data[,"u1"], data[,"u2"], pch = 16, main =
        "Located groups in point transects", xlim = c(0, 2*B),
        ylim = c(0, 2*B), col = data[,1], asp = 1)
    points(B, B, pch = "+", cex = 3)
    library(plotrix)
    draw.circle(B, B, B)
    hist(data[,"d"], col = "lightblue", breaks = 20, main =
        "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = FALSE)
    plot(table(data[,"gs"]), main = "Observed group sizes",
        ylab = "Frequency", frame = FALSE)
    par(op)
  }

  # Output
  list(type = type, nsites = nsites, lambda.group = lambda.group,
      alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1,
      B = B, data=data, habitat=habitat, N = N, N.true = N.true,
      groupsize=groupsize)
}


data <- simHDSg(type = "line")     # Defaults for line transect data
data <- simHDSg(type = "point")    # Default for point transect data
data <- simHDSg(lambda.group = 5)  # Much larger groups
data <- simHDSg(lambda.group = 5, alpha1 = 0) # No effect of groups size on p


# 9.2.2 Analysis in BUGS
# ------------------------------------------------------------------------
set.seed(1234)                 # we all create same data set
temp <- simHDSg(type="line")   # Execute function
data <- temp$data              # harvest data
B <- temp$B                    # Get strip half width
habitat <- temp$habitat        # habitat covariate
nsites <- temp$nsites          # Number of spatial replicates
groupsize <- data[,"gs"] -1    # Input groupsize-1 as data

M <- 400                        # Size of augmented data set is M
nz <- M-nrow(data)              # Number of "pseudo-groups" added
y <- c(data[,2],rep(0,nz))      # Indicator of capture (== 1 for all obs. groups)
nind <- nrow(data)              # Number of observed groups
site <- c(data[,1], rep(NA,nz)) # Site they belong to is unknown
d <- c(data[,5], rep(NA,nz))    # Their distance data are missing ...
groupsize <- c(groupsize, rep(NA,nz)) # .... as is their size
zst <- y                        # Starting values for data augmentation variable

# Bundle data and produce summary
str(bugs.data <- list (y=y, B=B, nind=nind, nsites=nsites, d=d, habitat=habitat,
   site=site, nz=nz, groupsize=groupsize))


# Define model in BUGS langauge
cat("
model{

  # Prior distributions for model parameters
  alpha0 ~ dunif(-10,10)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-10,10)
  beta1 ~ dunif(-10,10)
  lambda.group ~ dgamma(0.1, 0.1)
  # psi is a derived parameter
  psi <- sum(lambda[])/(nind+nz)

  # Individual level model: observations and process
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)                   # Data augmentation variables
    d[i] ~ dunif(0, B)                  # Distance is uniformly distributed
    groupsize[i] ~ dpois(lambda.group)  # Group size is Poisson

    log(sigma[i]) <- alpha0 +  alpha1*groupsize[i]
    mu[i] <- z[i]*exp(-d[i]*d[i]/(2*sigma[i]*sigma[i])) #p dep on dist class
    # here using the half normal detection function (Buckland et al. 2001)
    y[i] ~ dbern(mu[i])

    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
    zg[i]<- z[i]*(groupsize[i] + 1)      # Number of individuals in that group
  }

  for(s in 1:nsites){
    # Model for population size of groups
    N[s] ~ dpois(lambda[s])
    log(lambda[s])<- beta0 + beta1*habitat[s]
    site.probs[s]<- lambda[s]/sum(lambda[])
  }

  # Derived quantities
  G <- sum(z[])        # Total number of groups
  Ntotal <- sum(zg[])  # Total population size (all groups combined)
}
",fill=TRUE, file="model1.txt")

# Load some libraries, define MCMC settings, inits function and parameters to save
library("R2WinBUGS")
library("jagsUI")  #
ni <- 6000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3
inits <- function(){list(alpha0=0, alpha1=0.5, beta0=0, beta1=0, z=zst)}
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "G",
   "lambda.group")

# Call JAGS (ART 1.4 min), check convergence and summarize posterior distributions
out1 <- jags(bugs.data, inits, params, "model1.txt", n.thin=nt,
   n.chains=nc, n.burnin=nb,n.iter=ni)
traceplot(out1)   ;   print(out1, 3)


# 9.2.3 Imperfect observation of cluster size (no code)
# ------------------------------------------------------------------------

