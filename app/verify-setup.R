# Custom added
# Verify the app setup

library(readr)
library(log4r)
library(RcppTOML)


# Settings from env vars with default initial values
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


work_filename <- "df.csv"

df <- data.frame(var1=c(1, 2, 3, 4, 5),
                 var2=c(6, 7, 8, 9, 0))

tryCatch({

  work_path <- paste(dir_appwork, work_filename, sep = "/")
  print(paste("print: Creating file ", work_path))
  log4r::debug(logger, paste("Creating file ", work_path))
  fwrite(df, file=work_path)

}, error=function(e) {
  log4r::warn(logger, e)
}, warning=function(e) {
  log4r::warn(logger, e)
})


# Access data volume and read data file if configured
log4r::info(logger, paste("Data directory configuration set to =", data_dir))

if (data_dir == "") {
  print("data_dir var is not set")
  log4r::debug(logger, "data_dir var is not set")
} else {
  print("data_dir var is NOT empty string")
  log4r::debug(logger, "data_dir var is NOT empty string")

  # Verify that the data directory exists
  if (!dir.exists(data_dir)){
    print(paste("WARNING. Directory data_dir does NOT exist.", data_dir))
    log4r::warn(logger, paste("Directory data_dir does NOT exist.", data_dir))
    stop("ERROR. Incorrect app setup. A data directory is configured but was not found.")
  }

  # Attempt to read a file from the data directory
  data_file <- paste(data_dir, "existing.csv", sep = "/")
  if (file.exists(data_file)) {
    log4r::info(logger, paste("Data file was found. Exists at path =", data_file))

    # Open the data file
    tryCatch({

      df <- read_csv(data_file)
      df_summary <- summary(df)
      log4r::info(logger, paste("Data file summary =", df_summary))

    }, error=function(e) {
      log4r::warn(logger, e)
    }, warning=function(e) {
      log4r::warn(logger, e)
    })

  } else {
    log4r::warn(logger, paste("Data file was NOT found. Should exists at path =", data_file))
  }

}
