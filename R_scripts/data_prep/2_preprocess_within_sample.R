# Pre-process individual MS images
#
# Author: Paolo Inglese, Imperial College London

options(stringsAsFactors = FALSE)

# PACKAGES ----

tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})

pacman::p_load("SPUTNIK", "MALDIquant", "magrittr", "here")

# SETUP ----

FUNC_DIR <- here("R_scripts", "functions")
DATA_DIR <- here("DATA", "RAW_recal")
OUT_DIR <- here("DATA", "RData")

TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)

# SOURCE ----

source(here(FUNC_DIR, "misc", "load.R"))
source(here(FUNC_DIR, "preprocessing", "match_peaks_within_sample.R"))

# SAMPLES ----

samples_dirs <- c()
for (i in 1:length(TISSUES))
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

samples_names <- gsub(DATA_DIR, "", samples_dirs)
sapply(samples_names, function(z) {
  dir.create(here(OUT_DIR, z),
    recursive = TRUE
  )
})

# MATCH PEAKS WITHIN SAMPLES ----

for (i in 1:num_samples)
{
  cat("MS image", i, "/", num_samples, "\n")

  # Match the peaks within the sample

  cat("Matching peaks within the MS image...\n")

  peaks_matched_within <- matchPeaksWithinSample(
    imzMLPath = here(samples_dirs[i], "raw_recal.imzML"),
    minSignalPixelsFrac = 0.005, # Keep only if present in 0.5% of total pixels
    tolPPM = 10,
    deiso = T,
    verbose = T
  )

  # Calculate the average spectrum - it will be used as reference for matching
  # between samples

  cat("Calculate the average spectrum...\n")

  avg_spectrum <- createMassPeaks(
    mass = peaks_matched_within$mz,
    intensity = apply(peaks_matched_within$X, 2,
      mean,
      na.rm = T
    )
  )

  cat("Saving the matched peaks and the average spectrum...\n")

  save(avg_spectrum, file = here(
    OUT_DIR, samples_names[i],
    "avg_spectrum_within.RData"
  ))
  save(peaks_matched_within, file = here(
    OUT_DIR, samples_names[i],
    "X_matched_within.RData"
  ))

  rm(peaks_matched_within)
  gc()
}
rm(i)

# FILTER PEAKS ----

for (i in 1:num_samples)
{
  # Load sample
  cat(sprintf("%d/%d: %s\n", i, num_samples, samples_names[i]))

  data_env <- .load(here(
    OUT_DIR, samples_names[i],
    "X_matched_within.RData"
  )) ## Load X, mz, shape
  # Set the NA to 0
  data_env$peaks_matched_within$X[is.na(data_env$peaks_matched_within$X)] <- 0

  cat(
    "Dim X:", nrow(data_env$peaks_matched_within$X), "x",
    ncol(data_env$peaks_matched_within$X), "\n"
  )

  # Remove constant variables
  idx_const <- apply(data_env$peaks_matched_within$X, 2, var) == 0
  data_env$peaks_matched_within$X <- data_env$peaks_matched_within$X[, !idx_const]
  data_env$peaks_matched_within$mz <- data_env$peaks_matched_within$mz[!idx_const]

  # SPUTNIK
  cat("Creating MSI dataset...\n")
  msX <- msiDataset(
    data_env$peaks_matched_within$X,
    data_env$peaks_matched_within$mz,
    data_env$peaks_matched_within$shape[1],
    data_env$peaks_matched_within$shape[2]
  )
  cat("Pre-processing...\n")
  msX <- normIntensity(msX, "TIC")
  msX <- varTransform(msX, "log")

  rm(data_env)
  gc()

  # Extract the ROI
  ref_roi <- refAndROIimages(msiData = msX, roiMethod = "kmeans2")

  # Save ROI
  roi_mat <- ref_roi$ROI@values
  save(roi_mat, file = here(
    OUT_DIR, samples_names[i],
    "ROI_bin_before.RData"
  ))

  # Apply SPUTNIK

  # Global filter
  gpf <- globalPeaksFilter(msX, ref_roi$ROI)
  msX <- applyPeaksFilter(msX, gpf)

  # Pixel count filter
  pcf <- countPixelsFilter(msX, roiImage = ref_roi$ROI, aggressive = 1)
  msX <- applyPeaksFilter(msX, pcf)

  # Extract variables
  X <- getIntensityMat(msX)
  mz <- getMZ(msX)
  sz <- getShapeMSI(msX)
  cat(sprintf("Dim X: %d x %d\n", nrow(X), ncol(X)))

  # Average spectrum

  avg_spectrum <- createMassPeaks(
    mass = mz,
    intensity = apply(X, 2, mean, na.rm = T)
  )

  # Save the matched peaks and the average spectrum

  save(
    list = c("X", "sz", "mz"),
    file = paste0(
      OUT_DIR, samples_names[i],
      "X_matched_within_SPUTNIK.RData"
    )
  )
  save(avg_spectrum,
    file = paste0(
      OUT_DIR, samples_names[i],
      "avg_spectrum_within_SPUTNIK.RData"
    )
  )

  # Extract final ROI

  ref_roi <- refAndROIimages(msX, roiMethod = "kmeans2")
  roi_mat <- ref_roi$ROI@values
  save(roi_mat, file = paste0(
    OUT_DIR, samples_names[i],
    "ROI_bin_final.RData"
  ))

  rm(gpf, pcf, X, mz, sz)
  gc()
}
rm(i)

pacman::p_unload("SPUTNIK", "MALDIquant", "magrittr", "here")

cat("Done.\n")
