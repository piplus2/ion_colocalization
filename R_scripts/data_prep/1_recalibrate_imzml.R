# Re-calibrate the MS raw data using Palmitic acid (255.2330 m/z [M-H]-) as
# reference.
# A LOESS model is fitted on the distance between the theoretical and the
# observed reference m/z values and its predictions used to shift the entire
# m/z vector.
# Reference peak is matched within a search window of +/- 5ppm.
#
# Author: Paolo Inglese, Imperial College London

# PACKAGES ----

tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})

pacman::p_load("MALDIquant", "MALDIquantForeign", "here")

# SETUP ----

REF_MZ <- 255.2330 # Reference m/z for the re-calibration
DATA_DIR <- "DATA/RAW"
RECAL_DATA_DIR <- "DATA/RAW_recal"
if (!dir.exists(RECAL_DATA_DIR)) {
  cat("Creating dir for recalibrated raw data...\n")
  dir.create(RECAL_DATA_DIR, recursive = FALSE)
}
EXT <- "imzML" # raw file extension
TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)

# SOURCE ----

source(here("R_scripts", "functions", "misc", "vector.R"))
source(here("R_scripts", "functions", "preprocessing", "recalibrate.R"))

# DIRS AND FILES ----

# Read the data dirs and raw data paths

tissue_dirs <- lapply(TISSUES, function(z) {
  list.dirs(here(DATA_DIR, z), full.names = TRUE, recursive = FALSE)
})

stopifnot(all(unlist(lapply(tissue_dirs, length)) != 0)) # Check if empty

list_imzml_path <- .vector(mode = "list", length = length(TISSUES), names = TISSUES)
for (i in 1:NUM_TISSUES)
{
  list_imzml_path[[i]] <- array(NA, length(tissue_dirs[[i]]))

  for (j in 1:length(tissue_dirs[[i]]))
  {
    imzml_path <- list.files(
      path = tissue_dirs[[i]][j],
      pattern = paste0(EXT, "$"),
      full.names = TRUE, recursive = FALSE
    )

    stopifnot(length(imzml_path) == 1) # Check that only 1 raw file is present

    list_imzml_path[[i]][j] <- imzml_path

    rm(imzml_path)
  }

  rm(j)
}
rm(i)

num_samples <- unlist(lapply(list_imzml_path, length))
print(num_samples)

# START ----

for (i in 1:length(TISSUES))
{
  cat("Processing", TISSUES[i], "MS images...\n")

  for (j in 1:num_samples[i])
  {
    cat("MS image", j, "/", num_samples[i], "\n")

    peaks <- importImzMl(list_imzml_path[[i]][j], centroided = TRUE, verbose = FALSE)
    xy_coords <- matrix(NA, length(peaks), 2)
    for (k in 1:length(peaks))
    {
      xy_coords[k, ] <- c(peaks[[k]]@metaData$imaging$pos["x"],
                          peaks[[k]]@metaData$imaging$pos["y"])
    }

    image_shape <- c(length(unique(xy_coords[, 1])),
                     length(unique(xy_coords[, 2])))
    stopifnot(length(peaks) == prod(image_shape))

    print(list_imzml_path[[i]][j])

    sample_name <- strsplit(tissue_dirs[[i]][j], split = "\\/")[[1]]
    sample_name <- sample_name[length(sample_name)]

    save_dir <- here(RECAL_DATA_DIR, TISSUES[i], sample_name)
    if (!dir.exists(save_dir)) {
      cat("Creating the recalibrated MS image data directory...\n")
      dir.create(save_dir, recursive = TRUE)
    }

    recal_peaks <- recalibratePeaks(
      peaks,
      referenceMZ = REF_MZ,
      imageShape = image_shape,
      tolerancePPM = 10,
      genPlot = TRUE,
      saveRecal = TRUE,
      saveFilename = paste0(save_dir, "/raw_recal.imzML")
    )

    rm(peaks, image_shape, sample_name, save_dir, recal_peaks)
    gc()
  }
  rm(j)
}
rm(i)

pacman::p_unload("MALDIquant", "MALDIquantForeign", "here")

cat("Done.\n")
