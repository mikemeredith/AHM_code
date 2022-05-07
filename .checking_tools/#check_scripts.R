
# Source this script to check all the scripts in the current working directory
#   and its subdirectories.
# Files checked, errors, and times > 20 secs are listed in "#check_<date>.log".

# Finally check the log file for ERRORs, WARNINGs, NOTES and report at end of run
# (added 2021-08-15)

# Objects with initial '.' are not removed by rm(list = ls()), so persist during
#   running of individual scripts.

# -------- utility functions ----------------------------
# Function to convert secs to mins/hrs/days
.hms <- function(time) {
  if(time > 60*60*36)
    return(paste(round(time/3600/24, 2), "days"))
  if(time > 60*90)
    return(paste(round(time/3600, 2), "hrs"))
  if(time > 90)
    return(paste(round(time/60, 2), "mins"))
  return(paste(round(time, 2), "secs"))
}

# Function to compare current plotting 'par's with a default
.comparePars <- function(default,
    ignore=c("usr", "xaxp", "yaxp", "mai", "pin", "plt", "cex")) {
  default[ignore] <- NULL
  current <- par(no.readonly = TRUE)
  current[ignore] <- NULL
  isTRUE(all.equal(default, current))
}

# Function to compare current output with a previous .RData file
.compareValues <- function(oldFile=.imageFile, logFile=.logFile) {
  all.equal.function <- function(target, current, ...) {
    return(TRUE)  # skip check for functions (esp. nimble)
  }
  all.equal.jagsUI <- function(target, current, ...) {
    target$run.date <- NULL ; current$run.date <- NULL
    target$mcmc.info$elapsed.mins <- NULL ; current$mcmc.info$elapsed.mins <- NULL
    # Checking model text will show a mismatch even for tiny changes, eg in
    #   comments or white space
    target$model <- NULL ; current$model <- NULL
    all.equal.list(target, current, ...)
  }
  all.equal.cluster <- function(target, current, ...) {
    return(TRUE)  # skip check for clusters
  }

  if(file.exists(oldFile)) {
    outputNames <- ls(.GlobalEnv)
    oldResult <- new.env()
    load(oldFile, envir=oldResult)
    # check for missings in oldResult
    missing <- !outputNames %in% names(oldResult)
    if(sum(missing) > 0) {
      notGot <- paste(outputNames[missing], collapse=", ")
      cat("Old values not available for:", notGot, "\n",
            file = logFile, append = TRUE)
      outputNames <- outputNames[!missing]
    }
    same <- rep(TRUE, length(outputNames))
    for(i in seq_along(outputNames)) {
      this <- outputNames[i]
      new <- get(this, envir=.GlobalEnv)
      old <- get(this, envir=oldResult)
      # try first with default attribute check for relevant method
      tmp <- try(all.equal(new, old), silent=TRUE)
      if(!isTRUE(tmp)) {
        # try again without attribute check
        tmp <- try(all.equal(new, old,
            check.attributes=FALSE), silent=TRUE)
        if(inherits(tmp, "try-error")) {
          cat("NOTE: Can't compare results for:", this, "class:", class(new), "\n",
              file = logFile, append = TRUE)
          next
        }
      }
      same[i] <- isTRUE(tmp)
    }
    if(any(!same)) {
      bad <- paste(outputNames[!same], collapse=", ")
      cat("WARNING: Value differs from previous run for:", bad, "\n", file = logFile, append = TRUE)
    }
  } else {
    cat("NOTE: Old results not available\n", file = logFile, append = TRUE)
  }
}
# ..............................................................

# Get a listing of all .R files
files <- list.files(pattern = "[.][Rr]$", recursive = TRUE)
# Keep only file names not containing "#"
.ListOfFilesToCheck <- files[grep("#", files, invert=TRUE)]
head(.ListOfFilesToCheck)
length(.ListOfFilesToCheck)
# .ListOfFilesToCheck <- edit(.ListOfFilesToCheck)  # Edit the list if you wish
# length(.ListOfFilesToCheck)

# Attach rgdal before doing package count (does not unload properly) ## only if rgdal is needed!
# options("rgdal_show_exportToProj4_warnings"="none")
# library(rgdal)

# The old method of removing additional packages with package count no longer works as
# org:r-lib is inserted between  "Autoloads" and "package:base"
# Instead we check to see if the last item attached is the same:
( .lastPackage <- search()[2] )

# Set up log file
.logFile <- paste0("#check_", format(Sys.time(), "%Y-%m-%d_%H%M"), ".log")
cat("\nCheck log file for", getwd(), "\n", file = .logFile)
cat("started", format(Sys.time()), file = .logFile, append = TRUE)
cat("\n\n################################################################\n",
    file = .logFile, append = TRUE)

# Run the tests
.startTime <- Sys.time()
for(.i in seq_along(.ListOfFilesToCheck)) {
  cat("\n\n\n", .ListOfFilesToCheck[.i], "\n") ; flush.console()
  # clean up from previous run
  rm(list = ls())
  graphics.off()
  options(stringsAsFactors = FALSE)
  RNGkind("default", "default", "default")
  set.seed(42)  # for reproducible output when not set in the script
  while(search()[2] != .lastPackage)
    detach(pos=2)
  # create file names
  .fileStem <- substr(.ListOfFilesToCheck[.i], 1, nchar(.ListOfFilesToCheck[.i]) - 2)
  .pdfFile <- paste0(.fileStem, "#.pdf")
  .imageFile <- paste0(.fileStem, "#.Rdata")
  # run the script
  # cat("\n\n", .ListOfFilesToCheck[.i], "\n", file = .logFile, append = TRUE)
  cat("\n\n", .ListOfFilesToCheck[.i], "\nmodified", format(file.mtime(.ListOfFilesToCheck[.i])), "\n",
      file = .logFile, append = TRUE)
  cat(format(Sys.time()), "\n", file = .logFile, append = TRUE)
  pdf(.pdfFile)
  .defaultpars <- par(no.readonly = TRUE)

  .timing <- system.time(
    .returnValue <- try(source(.ListOfFilesToCheck[.i], chdir=TRUE)) )

  # Compare current output with a previous .RData file:
  .compareValues(oldFile=.imageFile, logFile=.logFile)
  # Check that par's have been restored:
  if(!.comparePars(.defaultpars)) {
    cat("NOTE: Plotting par's were not restored.\n", file = .logFile, append = TRUE)
    plot(0)
  }
  dev.off()
  sessionInfo <- sessionInfo()
  timeTaken <- .timing
  if(inherits(.returnValue, "try-error")) {
    cat("ERROR: ", .returnValue, file = .logFile, append = TRUE)
    .imageFile <- paste0(.fileStem, "_bad#.Rdata")
  }
  save.image(file = .imageFile)
  if(.timing[3] > 20)
    cat("Took", .hms(.timing[3]), file = .logFile, append = TRUE)
}
otime <- Sys.time() - .startTime  # overall time

# Check the log file for ERRORs, WARNINGs, NOTES
tmp <- readLines(.logFile)
nERROR <- sum(grepl("^ERROR", tmp))
nWARNING <- sum(grepl("^WARNING", tmp))
nNOTE <- sum(grepl("^NOTE", tmp))
issues <- nERROR > 0 || nWARNING > 0 || nNOTE > 0
if(issues)
  headsup <- paste("There were", nERROR, "ERRORs,", nWARNING, "WARNINGs,",
      nNOTE, "NOTEs.")


# Add overall time and sessionInfo to the end of the log file.
sink(.logFile, append=TRUE)
cat("\n\n#############################################\n")
cat("Completed", format(Sys.time()), file = .logFile, append = TRUE)
cat("\nOverall time taken", round(otime, 2), attr(otime, "units"), "\n")

if(issues)
  cat("\n", headsup, "\n", sep="", file = .logFile, append = TRUE)

cat("\nSessionInfo:\n\n")
sess <- sessionInfo()
sess$otherPkgs <- NULL
sess$loadedOnly <- NULL
print(sess)
cat("\nPackages used:\n")  # in alphabetical order
lns <- loadedNamespaces()
lns <- sort(lns[!(lns %in% sess$basePkgs)])
vers <- lapply(lns, packageVersion)
names(vers) <- lns
print(cbind(version=vers))
cat("\n")
sink(NULL)

# Display info in the Console (important if running multiple instances of R)
cat("\n\n##############################################\n")
cat("Checks completed for .R files in\n    ")
cat(getwd(), "\n")
cat("at", format(Sys.time()),
    " Overall time taken", round(otime, 2), attr(otime, "units"), "\n")

if(issues)
  cat("\n", headsup, "\n", sep="")
cat("For details see", dQuote(.logFile), "\n")
cat("##############################################\n\n")
