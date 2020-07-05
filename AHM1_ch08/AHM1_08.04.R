#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 8. Modeling abundance using hierarchical distance sampling (HDS)
# =========================================================================

# Approximate execution time for this code: 12 - 15 mins

library(AHMbook)
library(unmarked)

# 8.4 Hierarchical Distance Sampling
# ==================================


# 8.4.1 HDS data structure and model (no code)
# 8.4.2 HDS in unmarked (no code)

# 8.4.3 Example: Estimating the global population size of the Island Scrub Jay
# ------------------------------------------------------------------------
# Load, view and format the ISSJ data
library(unmarked)
data(issj)
round(head(issj), 2)

# Package things up into an unmarkedFrame
covs <- issj[,c("elevation", "forest", "chaparral")]
area <- pi*300^2 / 100^2             # Area in ha
jayumf <- unmarkedFrameDS(y=as.matrix(issj[,1:3]),
   siteCovs=data.frame(covs, area),
   dist.breaks=c(0, 100, 200, 300),
   unitsIn="m", survey="point")


# Fit model 1
(fm1 <- distsamp(~chaparral ~chaparral + elevation + offset(log(area)),
    jayumf, keyfun="halfnorm", output="abund"))

# Fit model 2
(fm2 <- distsamp(~1 ~chaparral + elevation + offset(log(area)),
    jayumf, keyfun="halfnorm", output="abund"))


# (pb <- parboot(fm1, fitstats, nsim=1000, report=5))
(pb <- parboot(fm1, fitstats, nsim=100, report=5))  # ~~~~~ for testing
(c.hat <- pb@t0[2] / mean(pb@t.star[,2]))  # c-hat as ratio of observed
                           # and mean of expected value of Chi2 (under H0)
                           # (see, e.g., Johnson et al., Biometrics, 2010)

residuals(fm1)             # Can inspect residuals
# ~~~~~ only plot if using screen device ~~~~~
if(dev.interactive(orNone=TRUE))
  plot(pb)                   # Not shown
print(pb)


# Standardize the covariates
sc <- siteCovs(jayumf)
sc.s <- scale(sc)
sc.s[,"area"] <- pi*300^2 / 10000  # Don't standardize area
siteCovs(jayumf) <- sc.s
summary(jayumf)


# Fit a bunch of models and produce a model selection table.
fall <- list()   # make a list to store the models

# With the offset output=abund is the same as output = density
fall$Null <- distsamp(~1 ~offset(log(area)), jayumf, output="abund")
fall$Chap. <- distsamp(~1 ~chaparral + offset(log(area)), jayumf,
    output="abund")
fall$Chap2. <- distsamp(~1 ~chaparral+I(chaparral^2)+offset(log(area)),
    jayumf, output="abund")
fall$Elev. <- distsamp(~1 ~ elevation+offset(log(area)), jayumf,
    output="abund")
fall$Elev2. <- distsamp(~1 ~ elevation+I(elevation^2)+offset(log(area)),
    jayumf, output="abund")
fall$Forest. <- distsamp(~1 ~forest+offset(log(area)), jayumf,
    output="abund")
fall$Forest2. <- distsamp(~1 ~forest+I(forest^2)+offset(log(area)),
    jayumf, output="abund")
fall$.Forest <- distsamp(~forest ~offset(log(area)), jayumf,
    output="abund")
fall$.Chap <- distsamp(~chaparral ~offset(log(area)), jayumf,
    output="abund")
fall$C2E. <- distsamp(~1 ~ chaparral + I(chaparral^2) + elevation +
    offset(log(area)),jayumf, output="abund")
fall$C2F2. <- distsamp(~1 ~chaparral + I(chaparral^2) + forest +
    I(forest^2)+offset(log(area)), jayumf,  output="abund")
fall$C2E.F <- distsamp(~forest ~chaparral+I(chaparral^2)+elevation+
    offset(log(area)), jayumf, output="abund")
fall$C2E.C <- distsamp(~chaparral ~chaparral + I(chaparral^2) + elevation +
    offset(log(area)), jayumf, output="abund")

# Create a fitList and a model selection table
(msFall <- modSel(fitList(fits=fall)))


# Check out the best model
fall$C2E.C


# Check out the goodness-of-fit of this model
(pb.try2 <- parboot(fall$C2E.C, fitstats, nsim=1000, report=5))

# Express the magnitude of lack of fit by an overdispersion factor
(c.hat <- pb.try2@t0[2] / mean(pb.try2@t.star[,2]))  #    Chisq

covs <- issj[,c("elevation", "forest", "chaparral")]
area <- pi*300^2 / 100^2             # Area in ha
jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
   siteCovs=data.frame(covs, area), numPrimary=1,
   dist.breaks=c(0, 100, 200, 300),
   unitsIn="m", survey="point")
sc <- siteCovs(jayumf)
sc.s <- scale(sc)
sc.s[,"area"] <- pi*300^2 / 10000  # Don't standardize area
siteCovs(jayumf) <- sc.s
summary(jayumf)


# Fit the model using gdistsamp and look at the fit summary
(nb.C2E.C <- gdistsamp( ~chaparral + I(chaparral^2) + elevation +
   offset(log(area)), ~1, ~chaparral, data =jayumf, output="abund",
   mixture="NB", K = 150))

gdistsamp(lambdaformula = ~chaparral + I(chaparral^2) + elevation +
    offset(log(area)), phiformula = ~1, pformula = ~chaparral,
    data = jayumf, output = "abund", mixture = "NB", K = 150)

(pb.try3 <- parboot(nb.C2E.C, fitstats, nsim=1000, report=5))

(c.hat <- pb.try3@t0[2] / mean(pb.try3@t.star[,2]))  #


# *Expected* population size for the sample points
getN <- function(fm, newdata=NULL)
   sum(predict(fm, type="lambda", newdata=newdata)[,1])
getN(nb.C2E.C)

# This does the same thing as the following three commands
X <- model.matrix(~chaparral+I(chaparral^2)+elevation+log(offset(area)),
        siteCovs(jayumf))
head(X) # The design matrix

# Prediction of total expected population size at the sample points
sum(exp(X %*% c(coef(nb.C2E.C, type="lambda"), 1)))

# Empirical Bayes estimates of posterior distribution:
# Pr(N=x | y, lambda, sigma) for x=0,1,...,K
re.jay <- ranef(nb.C2E.C, K = 150)

# *Realized* population size
sum(bup(re.jay, "mean"))


summary(jayumf) # Note the range of chaparral which we need to know

# Create a new data frame with area 28.27 ha, the area of a 300 m circle
chap.orig <- seq(0, 1, 0.01)    # Values from 0 to 1 prop. chaparral
chap.pred <- (chap.orig - mean(issj$chaparral)) / sd(issj$chaparral)
newdat <- data.frame(chaparral = chap.pred, elevation = 0, area=28.27)

# Expected values of N for covariate values in "newdat"
E.N <- predict(fall$C2E.C, type="state", newdata=newdat, appendData=TRUE)
head(E.N)

# Make a plot of the response curve for the grid of chaparral values
plot(chap.orig, E.N[,"Predicted"], xlab="Proportion chaparral",
    ylab="Predicted jay abundance", type="l", ylim = c(0, 20),
    frame = FALSE, lwd = 2)
matlines(chap.orig, E.N[,3:4], lty = 1, col = "grey", lwd = 1)

attributes(sc.s) # means are "scaled:center". SDs are "scaled:scale"

cruz.s <- cruz   # Created a new data set for the scaled variables
cruz.s$elevation <- (cruz$elevation-202)/125
cruz.s$chaparral <- (cruz$chaparral-0.270)/0.234
cruz.s$area <- (300*300)/10000 # The grid cells are 300x300m=9ha
EN <- predict(nb.C2E.C, type="lambda", newdata=cruz.s)

# Total population size (by summing predictions for all pixels)
getN(nb.C2E.C, newdata=cruz.s)

# Parametric bootstrap for CI
# A much faster function could be written to doing the sum
set.seed(2015)
(EN.B <- parboot(nb.C2E.C, stat=getN, nsim=1000, report=5))


library(raster)
cruz.raster <- stack(rasterFromXYZ(cruz.s[,c("x","y","elevation")]),
   rasterFromXYZ(cruz.s[,c("x","y","chaparral")]),
   rasterFromXYZ(cruz.s[,c("x","y","area")]))
names(cruz.raster) # These should match the names in the formula

plot(cruz.raster)                      #  not shown
# Elevation map on the original scale (not shown)
plot(cruz.raster[["elevation"]]*125 + 202, col=topo.colors(20),
    main="Elevation (in feet) and Survey Locations", asp = 1)
points(issj[,c("x","y")], cex=0.8, pch = 16)


EN.raster <- predict(nb.C2E.C, type="lambda", newdata=cruz.raster)
plot(EN.raster, col = topo.colors(20), asp = 1)   # See Fig. 8-8

