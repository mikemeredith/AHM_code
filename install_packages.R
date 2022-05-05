

# Install or update packages needed to run the code in this repository
# ====================================================================
# All of these are available from CRAN

# I suggest first running
update.packages(ask='graphics',checkBuilt=TRUE)
# to ensure everything is up to date, including dependencies.

needed <- c("AHMbook", "AICcmodavg", "berryFunctions", "coda", "corrplot", "denstrip",
    "devtools", "doParallel", "fields", "foreach", "jagsUI", "lme4", "mcmcOutput",
    "nimble", "plotrix", "R2OpenBUGS", "R2WinBUGS", "raster", "rjags",
    "sp", "unmarked", "wiqid")  # package 'rgdal' is no longer needed
got <- rownames(installed.packages())

( notgot <- needed[!needed %in% got] )

install.packages(notgot, dependencies=TRUE)

# 'devel' versions of packages
# ----------------------------
# If you want to try out devel versions of packages from GitHub, install these
#  AFTER 'update.packages' as that will "downdate" to the latest CRAN version.
# For example:
# remotes::install_github("mikemeredith/AHMbook")
# packageVersion("AHMbook")
# remotes::install_github("rbchan/unmarked")
# remotes::install_github("kenkellner/unmarked")
# packageVersion("unmarked")
# remotes::install_github("kenkellner/jagsUI")
# packageVersion("jagsUI")

