#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 11. Hierarchical models for communities
# =========================================================================

library(AHMbook)
data(MHB2014)

# ~~~~ For this we need the 'all10' matrix from section 7 ~~~~~~~~~~~~~
load("AHM1_11.07_all10.RData")
# ~~~~~~ this section requires the data prepared in section 11.3 ~~~~~~~~~~
source("AHM1_11.03.R")
# ~~~~ also this from section 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nz <- 150
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 11.8 Inferences based on the estimated Z matrix: similarity among sites and species
# ===================================================================================


# Plug MCMC samples for full z matrix into 3D array
str(all10)
nsite <- 267
nspec <- 215
nsamp <- dim(all10)[1]        # 1200 MCMC samples
# ~~~~~ the order of the columns in all10 is different ~~~~~~~~~~~~~~~~~~~~~
# The elements of the Z matrix are now in columns 1945 to 59349
# Safer to use parameter names and grep instead of column numbers.
nms <- colnames(all10)
is.z <- grep("^z", nms)
range(is.z)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- array(NA, dim = c(nsite, nspec, nsamp))
Jacc <- array(NA, dim = c(nsite, nspec, nsamp))
for(j in 1:nsamp){    # Fill z matrix by column (default)
   cat(paste("\nMCMC sample", j, "\n"))
   # z[,,j] <- all10[j, 1937:59341]
   z[,,j] <- all10[j, is.z]
}

# Restrict computations to observed species
zobs <- z[,1:145,]      # Species 1 to 145

# Compute Jaccard index for sites and for species
Jsite <- array(NA, dim = c(nsite, nsamp))
Jspec <- array(NA, dim = c(145, nsamp))


# Choose reference site and species for Jaccard indices
ref.site <- 1         # Just choose first site
ref.species <- 13     # European Sparrowhawk (check object 'obs.occ')

# Get posterior distributions for Jsite and Jspec (for references)
for(k in 1:nsamp){
  for(i in 1:nsite){ # Jaccard index for sites (in terms of shared species)
    Jsite[i,k] <- sum(zobs[ref.site,,k] * zobs[i,,k]) /
      (sum(zobs[ref.site,,k]) + sum(zobs[i,,k]) -
       sum(zobs[ref.site,,k] * zobs[i,,k]))
  }
  # for(i in 1:(nspec-nz)){ # Jacc. index for species (in terms of shared sites)
  for(i in 1:145){ # Jacc. index for species (in terms of shared sites) # ~~~~ 145 species only
    Jspec[i,k] <- sum(zobs[,ref.species,k] * zobs[,i,k]) /
      (sum(zobs[,ref.species,k]) + sum(zobs[,i,k]) -
      sum(zobs[,ref.species,k] * zobs[,i,k]))
  }
}
# NA's arise when a site has no species or a species no sites

# Get posterior means, standard deviations and 95% CRI
# Jaccard index for sites, compared to reference site 1
pm <- apply(Jsite, 1, mean, na.rm = TRUE)  # Post. mean of Jsite wrt. site 1
psd <- apply(Jsite, 1, sd, na.rm = TRUE)   # Post. SD of Jsite wrt. site 1
cri <- apply(Jsite, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm = TRUE)) # CRI
cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])


# Make a map of Jaccard site indices (Fig. 11-26)
x <- 3        # proportional size of plotting symbol
plot(MHB2014$sites$coordx, MHB2014$sites$coordy, xlab = "x coordinate",
    ylab = "y coordinate", cex = x*pm, asp = 1, pch = 16)
points(MHB2014$sites$coordx[which(pm == 1)], MHB2014$sites$coordy[which(pm == 1)],
    cex = x*pm, col = "red", pch = 16)


# Jaccard index for species, compared to reference species
# (species 13, European Sparrowhawk)
pm <- apply(Jspec, 1, mean, na.rm = TRUE)  # Post. mean of Jspec wrt. species 1
psd <- apply(Jspec, 1, sd, na.rm = TRUE)   # Post. mean of Jspec wrt. species 1
cri <- apply(Jspec, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm = TRUE)) # CRI
tmp <- cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])
rownames(tmp) <- names(obs.occ)
# print(tmp])          # print in systematic order
print(tmp)          # print in systematic order
print(tmp[rev(order(tmp[,1])),]) # print in order of decreasing Jacc. values
plot(1:145, tmp[rev(order(tmp[,1])),1])   # can also plot

