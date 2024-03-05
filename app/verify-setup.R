# Custom added
# Verify the app setup

library(readr)
library(log4r)
library(RcppTOML)


# Settings from env vars
data_dir <- ""
log_level <- "DEBUG"
log_file <- "/rlogs/app.log"


tryCatch({
  print("Reading environment variables from .Renviron")

  readRenviron(".Renviron")

  data_dir <- Sys.getenv("DATA_DIR")
  log_level <- Sys.getenv("R_LOGLEVEL")
  log_file <- Sys.getenv("R_LOGFILE")

}, error=function(e) {
  print(paste("Error while trying to get env vars ", e))
}, warning=function(e) {
  print(paste("Warning while trying to get env vars ", e))
})

# Logging
print(paste("Running app with log level = ", log_level, ". Logging to ", log_file))

logger <- log4r::logger(log_level, appenders = file_appender(log_file))

log4r::info(logger, paste("START app.R. Logging enabled with log level =", log_level))

# Version
# Get the version nr from the project toml file
ver <- "unset"
project_filename <- "project.toml"

tryCatch({
  fd <- read_file(project_filename)
  toml <- parseTOML(fd, verbose = FALSE, fromFile=FALSE, includize=FALSE, escape=TRUE)
  ver <- toml["app"]$app$version
  log4r::info(logger, paste("Running app version =", ver))
}, error=function(e) {
  log4r::warn(logger, e)
}, warning=function(e) {
  log4r::warn(logger, e)
})


# Verify can create a directory and write to it
dir_appwork <- "appwork"

tryCatch({

  print(paste("print: Creating directory ", dir_appwork))
  log4r::debug(logger, paste("Creating directory ", dir_appwork))

  if (!dir.exists(dir_appwork)){
    dir.create(dir_appwork, showWarnings = TRUE, recursive = TRUE)
  } else {
    print("Directory appwork already exists.")
    log4r::debug(logger, paste("Directory appwork already exists: ", dir_appwork))
  }
}, error=function(e) {
    log4r::warn(logger, e)
}, warning=function(e) {
    log4r::warn(logger, e)
})

