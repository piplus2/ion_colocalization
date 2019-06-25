## Generate and save the multi_data and sampled pixels structures
#
# Author: Paolo Inglese, Imperial College London

options(stringsAsFactors = FALSE)

# PACKAGES ----
tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})

pacman::p_load("here")

# SETUP ----

APPLY_SCALING <- FALSE # Apply UV-scaling to the data (within image)

# SOURCE ----

source(here("R_scripts", "functions", "misc", "set_dir.R"))
source(here("R_scripts", "functions", "misc", "load.R"))
source(here("R_scripts", "functions", "misc", "vector.R"))

# START ----

DATA_DIR <- .setDir(here("DATA", "RData"), createIfNotExist = TRUE)
rm(DATA_DIR)

TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)

SETS <- c("cv", "test")

INPUT_FILENAME <- c(
  cv = "X_matched_between_SPUTNIK.RData",
  test = "X_matched_with_cmz.RData"
)

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
  rm(d_tmp)
}
rm(i)

# MULTI.DATA ----

cat("Loading the intensity matrices...\n")

# Split the sample directories in 'cv' and 'test' set

samples_dirs <- lapply(idx_samples, function(z) list_samples_dirs[z])

Y_CV <- tissue_labels[idx_samples$cv]
Y_TEST <- tissue_labels[idx_samples$test]

for (APPLY_SCALING in c(FALSE, TRUE))
{
  multi_data <- .vector(mode = "list", length = 2, names = SETS)

  for (set in SETS)
  {
    multi_data[[set]] <- .vector(
      mode = "list", length = num_samples[set],
      names = paste0("sample_", idx_samples[set])
    )

    for (j in 1:num_samples[set])
    {
      cat(paste0(
        "Set: ", set, ", MS image: ", j,
        ifelse(APPLY_SCALING, ", scaling...", ""), "\n"
      ))

      sample_dir <- samples_dirs[[set]][j]

      # Load the intensity matrix and the spatial dimensions

      data_env <- .load(here(sample_dir, INPUT_FILENAME[set]))
      data_env$X[is.na(data_env$X)] <- 0
      stopifnot(all(!is.na(data_env$X)))

      # Load the tissue mask

      roi_env <- .load(here(sample_dir, "ROI_bin_final.RData"))
      shape <- dim(roi_env$roi_mat)

      # Remove 282.2511 m/z manually (it corresponds to the column 3)
      cat("Removing column 3 = 282.2511 m/z\n")
      multi_data[[set]][[j]] <- data_env$X[c(roi_env$roi_mat) == 1, -3]

      # UV-scaling
      if (APPLY_SCALING) {
        multi_data[[set]][[j]] <- scale(multi_data[[set]][[j]])
      }
      if (set == "test") {
        multi_data[[set]][[j]][is.na(multi_data[[set]][[j]])] <- 0
      }
      stopifnot(!any(is.na(multi_data[[set]][[j]])))

      rm(data_env, roi_env)
      gc()
    }
    rm(j)
  }
  rm(set)

  cat("Saving...\n")
  save(multi_data, file = here(
    "DATA", "RData", "multi_data",
    ifelse(APPLY_SCALING, "_SCALED", ""), ".RData"
  ))
}

pacman::p_unload("here")

cat("Done.\n")
