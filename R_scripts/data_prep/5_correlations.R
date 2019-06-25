# Calculate and save the correlations using randomly sampled pixels
#
# Author: Paolo Inglese, Imperial College London

options(stringsAsFactors = FALSE)

# PACKAGES ----

tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})
pacman::p_load("parallel", "here", "psych")

# SETUP ----

NUM_CORES <- 5 # Num. of CPU cores used for the parallel computation -
# It must be smaller than total num. of cores
APPLY_SCALING <- T # UV-scaling before calculating the correlations

# SOURCE ----

source(here("R_scripts", "functions", "misc", "vector.R"))
source(here("R_scripts", "functions", "misc", "load.R"))

# START ----

RESULTS_DIR <- here("output", "cor_mats_rand_sampling")
if (!dir.exists(RESULTS_DIR)) {
  cat("Creating output directory...\n")
  dir.create(RESULTS_DIR, recursive = TRUE)
}
rm(RESULTS_DIR)

TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)

SETS <- c("cv", "test")
INPUT_FILENAME <- c(
  cv = "X_matched_between_SPUTNIK.RData",
  test = "X_matched_with_cmz.RData"
)

# DIRS/FILES

# Cross-validation / test MS images

data_env <- .load(here("DATA", "RData", "split_idx_cv_test.RData"))

idx_samples <- .vector(mode = "list", length = 2, names = SETS)
idx_samples$cv <- data_env$idx_cv
idx_samples$test <- data_env$idx_test

rm(data_env)
gc()

num_samples <- unlist(lapply(idx_samples, length))
print(num_samples)

# MS images directories and labels

list_samples_dirs <- c()
tissue_labels <- c()
for (i in 1:length(TISSUES))
{
  d_tmp <- list.dirs(here("DATA", "RData", TISSUES[i]),
    recursive = FALSE,
    full.names = TRUE
  )
  list_samples_dirs <- c(list_samples_dirs, d_tmp)
  tissue_labels <- c(tissue_labels, rep(TISSUES[i], length(d_tmp)))
}
rm(i)

# Split the sample directories in 'cv' and 'test' set

samples_dirs <- lapply(idx_samples, function(z) list_samples_dirs[z])

Y_CV <- factor(tissue_labels[idx_samples$cv])
Y_TEST <- factor(tissue_labels[idx_samples$test])

cat("Saving the tissue sections labels...\n")
save(list = c("Y_CV", "Y_TEST"), file = here("DATA", "RData", "labels.RData"))

# LOAD ----

cat("Loading the multi_data...\n")
load(here("DATA", "RData", "multi_data", ifelse(APPLY_SCALING, "_SCALED", ""), ".RData"))

# CORRELATIONS ----

# Determine the smallest ROI size. This is to define the maximum sample size
num_pixels_images <- c()
for (set in SETS)
{
  num_pixels_images <- c(
    num_pixels_images,
    sapply(1:num_samples[set], function(z)
      nrow(multi_data[[set]][[z]]))
  )
}
cat(paste0("Smallest ROI size: ", (min(num_pixels_images)), "\n"))

# Sample the pixels and save them for the correlation of the original and scaled
# intensities

cat("Calculate the correlation features. This may take a while...\n")

MAX_SAMPLE <- 1000 # It must be smaller than min(num_pixels_images)
NUM_PIXELS <- seq(100, MAX_SAMPLE, 100)
NUM_REPS <- 500 # Num. of repetitions

cat("Starting the parallel cluster...\n")
cl <- makeCluster(NUM_CORES)

for (n in 1:length(NUM_PIXELS))
{
  sample_size <- NUM_PIXELS[n]

  message(paste0("Sample size: ", sample_size, "\n"))

  # Correlations are split into two sets: 'cv' and 'test'. Each group contains
  # a list (one element per tissue section) of the replicated correlations
  # (num. replicates x num. correlations).

  cor_matrix <- .vector(mode = "list", length = 2, names = SETS)

  # CV and test sets
  for (set in SETS)
  {
    cor_matrix[[set]] <- vector(mode = "list", length = num_samples[set])

    for (i in 1:num_samples[set])
    {
      cat(paste0("set: ", set, ", MS image: ", i, "\n"))

      intensity_matrix <- multi_data[[set]][[i]]

      clusterExport(cl = cl, varlist = c("intensity_matrix", "sample_size"))

      # Calculate the repeated correlation matrices
      cor_matrix_tmp <- parSapply(cl = cl, X = 1:NUM_REPS, FUN = function(X) {
        require(psych)

        num_pixels <- nrow(intensity_matrix)

        sample_pixels <- sort(sample(num_pixels, sample_size))
        stopifnot(length(sample_pixels) == sample_size)

        cor_test <- corr.test(intensity_matrix[sample_pixels, ],
          use = "p",
          method = "spearman", adjust = "BH"
        )
        cor_test$r[is.na(cor_test$r)] <- 0

        rm(sample_pixels)

        # Calculate the significance of the correlations
        sym_p <- cor_test$p
        sym_p[lower.tri(sym_p, diag = T)] <- 0
        sym_p <- sym_p + t(sym_p)

        # Set the non-significant (p-value >= 0.05) to 0
        cor_test$r[sym_p >= 0.05] <- 0
        cor_test$r[upper.tri(cor_test$r, diag = FALSE)]
      })

      cor_matrix[[set]][[i]] <- t(cor_matrix_tmp)
    }
    rm(i)
  }
  rm(set)

  cat("Saving...\n")
  saveRDS(cor_matrix, file = here(
    "output",
    "cor_mats_rand_sampling",
    paste0(
      "cor_mat_", sample_size, "px_", NUM_REPS, "reps",
      ifelse(APPLY_SCALING, "_SCALED", ""), ".rds"
    )
  ))
}
rm(n)

stopCluster(cl = cl)
closeAllConnections()
gc(reset = TRUE)

pacman::p_unload("parallel", "here", "psych")

cat("Done.\n")
