
# Source this script to check all the scripts in the current working directory
#   and its subdirectories.
# Files checked, errors, and times > 20 secs are listed in "#check_<date>.log".

# Objects with initial '.' are not removed by rm(list = ls()), so persist during
#   running of individual scripts.

# -------- preliminaries ----------------------------
# Function to convert secs to mins/hrs
.hms <- function(time) {
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

# Get a listing of all .R files
files <- list.files(pattern = "[.][Rr]$", recursive = TRUE)
# Keep only file names not containing "#"
.ListOfFilesToCheck <- files[grep("#", files, invert=TRUE)]
head(.ListOfFilesToCheck)
length(.ListOfFilesToCheck)
# .ListOfFilesToCheck <- edit(.ListOfFilesToCheck)  # Edit the list if you wish
# length(.ListOfFilesToCheck)

# Attach rgdal before doing package count (does not unload properly)
library(rgdal)
.oldPackageCount <- length(search())
.oldWD <- getwd()
# .ow <- options(warn=1)

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
  while(length(search()) > .oldPackageCount)
    detach(pos=2)
  # create file names
  .fileStem <- substr(.ListOfFilesToCheck[.i], 1, nchar(.ListOfFilesToCheck[.i]) - 2)
  .pdfFile <- paste0(.fileStem, "#.pdf")
  .imageFile <- paste0(.fileStem, "#.Rdata")
  # run the script
  cat("\n\n", .ListOfFilesToCheck[.i], "\n", file = .logFile, append = TRUE)
  cat(format(Sys.time()), "\n", file = .logFile, append = TRUE)
  pdf(.pdfFile)
  defaultpars <- par(no.readonly = TRUE)
  timing <- system.time(
    returnValue <- try(source(.ListOfFilesToCheck[.i], chdir=TRUE)) )
  # check that par's have been restored:
  if(!.comparePars(defaultpars)) {
    cat("Plotting par's were not restored.\n", file = .logFile, append = TRUE)
    plot(0)
  }
  dev.off()
  sessionInfo <- sessionInfo()
  save.image(file = .imageFile)
  if(inherits(returnValue, "try-error"))
    cat(returnValue, file = .logFile, append = TRUE)
  if(timing[3] > 20)
    cat("Took", .hms(timing[3]), file = .logFile, append = TRUE)
}
otime <- Sys.time() - .startTime  # overall time

# Add overall time and sessionInfo to the end of the log file.
sink(.logFile, append=TRUE)
cat("\n\n#############################################\n")
cat("Completed", format(Sys.time()), file = .logFile, append = TRUE)
cat("\nOverall time taken", round(otime, 2), attr(otime, "units"), "\n")
cat("\nSessionInfo:\n\n")
print(sessionInfo())
sink(NULL)

# Display info in the Console (important if running multiple instances of R)
cat("\n\n##############################################\n")
cat("Checks completed for .R files in\n")
cat(getwd(), "\n")
cat("at", format(Sys.time()), "\n")
cat("For details see", dQuote(.logFile), "\n")
cat("##############################################\n\n")
