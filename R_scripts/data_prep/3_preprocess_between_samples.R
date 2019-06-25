# Match spectral peaks between MS images
#
# Author: Paolo Inglese, Imperial College London


options(stringsAsFactors = FALSE)

# SETUP ----

tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})

pacman::p_load("MALDIquant", "here")

# SOURCE ----
source(here("R_scripts", "functions", "misc", "set_dir.R"))
source(here("R_scripts", "functions", "misc", "load.R"))
source(here("R_scripts", "functions", "preprocessing", "match_peaks_between_samples.R"))

# START ----
DATA_DIR <- .setDir(here("DATA", "RData"), createIfNotExist = TRUE)

EXT <- "imzML"

TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)

# DIRS ----

samples_dirs <- c()
for (i in 1:NUM_TISSUES)
{
  tmp <- list.dirs(here(DATA_DIR, TISSUES[i]),
    full.names = FALSE,
    recursive = FALSE
  )
  samples_dirs <- c(samples_dirs, tmp)

  rm(tmp)
}
rm(i)

num_samples <- length(samples_dirs)
print(num_samples)

# MATCH PEAKS BETWEEN SAMPLES ----

cat("Matching the peaks between MS images...\n")

# Load the average spectra to be used as reference

cat("Loading the average spectra...\n")

ref_spectra <- .vector(mode = "list", length = num_samples, names = c(1:num_samples))
for (i in 1:num_samples)
{
  cat("MS image:", i, "\n")

  data_env <- .load(here(
    samples_dirs[i],
    "avg_spectrum_within_SPUTNIK.RData"
  )) # Load avg_spectrum
  ref_spectra[[i]] <- data_env$avg_spectrum

  rm(data_env)
  gc()
}
rm(i)

cat("Matching peaks...\n")

cmz <- matchPeaksBetweenSamples(
  samplesPath = samples_dirs,
  refPeaksList = ref_spectra,
  deiso = TRUE, # De-isotope
  freqThreshold = 1.0, # Only peaks present in all the MS images
  tolerance = 20, # Tolerance in PPM
  inFile = "X_matched_within_SPUTNIK.RData",
  outFile = "X_matched_between_SPUTNIK.RData",
  verbose = TRUE
)

save(cmz, file = here(DATA_DIR, "cmz_between_samples.RData"))

pacman::p_load("MALDIquant", "here")

cat("Done.\n")
