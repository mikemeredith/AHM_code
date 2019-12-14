

# Install or update packages needed to run the code in this repository
# ====================================================================
# All of these are available from CRAN

# I suggest first running
update.packages(ask='graphics',checkBuilt=TRUE)
# to ensure everything is up to date, including dependencies.

needed <- c("AHMbook", "unmarked", "AICcmodavg", "sp", "rgdal", "plotrix", "raster",
  "lme4", "R2WinBUGS", "R2OpenBUGS", "jagsUI", "denstrip", "rjags", "coda", "devtools")
got <- rownames(installed.packages())

( notgot <- needed[!needed %in% got] )

install.packages(notgot, dependencies=TRUE)

# 'devel' versions of packages
# ----------------------------
# If you want to try out devel versions of packages from Github, install these
#  AFTER 'update.packages' as that will "downdate" to the latest CRAN version.
# For example:
# devtools::install_github("mikemeredith/AHMbook")
# packageVersion("AHMbook")
# news(package="AHMbook")

