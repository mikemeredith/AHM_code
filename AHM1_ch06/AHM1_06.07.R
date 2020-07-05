#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# Run time with the full number of iterations: 1.8 hrs

library(AHMbook)
library(unmarked)

# 6.7 Study of some assumption violations using function simNmix
# ==============================================================


# simreps <- 1000                 # Number of data sets created/analysed
                                  # ~~~~ 1000 takes around 90 mins
simreps <- 100                    # ~~~~ reduced number for testing
MLE <- array(dim = c(5, simreps)) # Array to hold MLEs

for(i in 1:simreps){              # Create and analyse 1000 data sets
  cat("*** Simrep number", i, "***\n")
  # Create data set with some extra (here: open populations)
  data <- simNmix(mean.lam=exp(1), mean.p=0.5, beta2.lam=1,
    beta3.p=1, beta.p.survey=1, open.N=TRUE, show.plot=FALSE)
  # Analyse data set with standard model (here: assuming closure)
  umf <- unmarkedFramePCount(y=data$C, siteCovs =
    data.frame(cov2=data$site.cov[,2], cov3=data$site.cov[,3]),
    obsCovs = list(survey.cov = data$survey.cov))
  fm <- pcount(~cov3+survey.cov ~cov2, umf, se = FALSE)
  # Save MLEs
  MLE[,i] <- coef(fm)
}

