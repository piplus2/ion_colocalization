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

data_env <- .load(here("DATA", "RData", "split_idx_cv_test.RData"))
idx_samples <- .vector(mode = "list", length = 2, names = SETS)
idx_samples$cv <- data_env$idx_cv
idx_samples$test <- data_env$idx_test
rm(data_env)

cat("Loading the average spectra...\n")

ref_spectra <- .vector(mode = "list", length = length(idx_samples$cv),
                       names = c(1:length(idx_samples$cv)))
for (i in 1:length(idx_samples$cv))
{
  cat("MS image:", i, "\n")

  data_env <- .load(here(
    samples_dirs[idx_samples$cv[i]],
    "avg_spectrum_within_SPUTNIK.RData"
  )) # Load avg_spectrum
  ref_spectra[[i]] <- data_env$avg_spectrum

  rm(data_env)
  gc()
}
rm(i)

cat("Matching peaks...\n")

cmz <- matchPeaksBetweenSamples(
  samplesPath = samples_dirs[idx_samples$cv],
  refPeaksList = ref_spectra,
  deiso = TRUE, # De-isotope
  freqThreshold = 1.0, # Only peaks present in all the MS images
  tolerance = 20, # Tolerance in PPM
  inFile = "X_matched_within_SPUTNIK.RData",
  outFile = "X_matched_between_SPUTNIK.RData",
  verbose = TRUE
)

save(cmz, file = here(DATA_DIR, "cmz_between_samples.RData"))

# Match the test samples
matchPeaksWithCMZ(samples_dirs[idx_samples$test],
                  cmz,
                  tolPPM = 20,
                  inputFilename = 'X_matched_between_SPUTNIK.RData',
                  outFilename = 'X_matched_with_cmz.RData',
                  verbose = TRUE)

pacman::p_load("MALDIquant", "here")

cat("Done.\n")
