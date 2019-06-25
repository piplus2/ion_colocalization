# Generate the correlation matrices after applying a random offset to the
# peak intensities.

options(stringsAsFactors = FALSE)

# PACKAGES ----

tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})
pacman::p_load("parallel", "here")

# SOURCE ----

source(here("R_scripts", "functions", "misc", "set_dir.R"))
source(here("R_scripts", "functions", "misc", "load.R"))
source(here("R_scripts", "functions", "misc", "vector.R"))

# SETUP ----

TISSUE_OF_INTEREST = "breast"  # Select one tissue type: "breast", "colorectal", "ovarian"
                               # Correlations are calculated only on one type
NUM_CORES <- 5 # Num. of CPU cores used for the parallel computation -
# It must be smaller than total num. of cores

APPLY_SCALING <- FALSE  # If TRUE the correlations are calculated from the
                        # UV-scaled peaks intensities

# START ----

# Set the parameters for the replicates
MAX_SAMPLE <- 1000  # It must be smaller than min(num_pixels_images)
NUM_PIXELS <- seq(100, MAX_SAMPLE, 100)
NUM_REPS <- 100  # Num. of repetitions

TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)

output_dir <- .setDir(here("output", "cor_mats_offset", TISSUE_OF_INTEREST),
                      createIfNotExist = TRUE)
rm(output_dir)

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
  d_tmp <- list.dirs(here("DATA", "RData",  TISSUES[i]),
    recursive = FALSE,
    full.names = TRUE
  )
  list_samples_dirs <- c(list_samples_dirs, d_tmp)
  tissue_labels <- c(tissue_labels, rep(TISSUES[i], length(d_tmp)))
}
rm(i)

# Split the sample directories in 'cv' and 'test' set

samples_dirs <- lapply(idx_samples, function(z) list_samples_dirs[z])

load(here("DATA", "RData", "labels.RData"))

# Load the multi_data structure
cat("Loading multi_data...\n")
multi_data_env <- .load(here("DATA",
                             "RData",
                             paste0("multi_data", ifelse(APPLY_SCALING, "_SCALED", ""), ".RData")))

# Select only one tissue type

cat(paste0("Analysing ", TISSUE_OF_INTEREST, " samples...\n"))

idx.samples.of.interest <- list(cv = sort(which(Y_CV == TISSUE_OF_INTEREST)))

# Generate the sampled pixels per tissue type ----

cat("Sampling pixels from images...\n")

# The same pixels will be used to calculate the correlation of the offset data
sample_pixels <- .vector(mode = "list", length = length(NUM_PIXELS), names = paste0(NUM_PIXELS, "px"))

for (n in 1:length(NUM_PIXELS))
{
  cat(paste0("num. pixels = ", NUM_PIXELS[n], "\n"))
  sample_size <- NUM_PIXELS[n]
  sample_pixels[[n]] <- vector(mode = "list", length = 2)
  names(sample_pixels[[n]]) <- SETS
  for (set in "cv")
  {
    sample_pixels[[n]][[set]] <- vector(mode = "list", length = length(idx.samples.of.interest[[set]]))
    names(sample_pixels[[n]][[set]]) <- paste0("sample", c(1:length(idx.samples.of.interest[[set]])))
    for (i in 1:length(idx.samples.of.interest[[set]]))
    {
      num_pixels <- nrow(multi_data_env$multi_data[[set]][[idx.samples.of.interest[[set]][i]]])

      sample_pixels[[n]][[set]][[i]] <- matrix(NA, NUM_REPS, sample_size)
      dimnames(sample_pixels[[n]][[set]][[i]]) <- list(paste0("rep", 1:NUM_REPS), NULL)
      cat(paste0("set: ", set, ", MS image: ", idx.samples.of.interest[[set]][i], "\n"))
      for (rx in 1:NUM_REPS)
      {
        sample_pixels[[n]][[set]][[i]][rx, ] <- sort(sample(num_pixels, sample_size))
      }
      rm(rx)
    }
    rm(i)
  }
  rm(set)
}
rm(n)

cat("Saving sampled pixels...\n")
save(sample_pixels, file = here(
  "output",
  "cor_mats_offset",
  TISSUE_OF_INTEREST,
  paste0("sample_pixels_", TISSUE_OF_INTEREST, ".RData")
))

# Original correlations ----

cat("Correlations from original data...\n")

cat(paste0("Starting the parallel cluster using ", NUM_CORES, " cores...\n"))
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
  for (set in "cv")
  {
    cor_matrix[[set]] <- vector(mode = "list", length = length(idx.samples.of.interest[[set]]))

    for (i in 1:length(idx.samples.of.interest[[set]]))
    {
      cat(paste0("set: ", set, ", MS image: ", idx.samples.of.interest[[set]][i], "\n"))

      intensity_matrix <- multi_data_env$multi_data[[set]][[idx.samples.of.interest[[set]][i]]]

      sample.pixels.mat <- sample_pixels[[n]][[set]][[i]]
      clusterExport(cl = cl, varlist = c("intensity_matrix", "sample_size", "sample.pixels.mat"))

      # Calculate the repeated correlation matrices
      cor_matrix_tmp <- parSapply(cl = cl, X = 1:NUM_REPS, FUN = function(X) {
        require(psych)

        sample_pixels <- sample.pixels.mat[X, ]
        stopifnot(length(sample_pixels) == sample_size)

        cor_test <- corr.test(intensity_matrix[sample_pixels, ],
          use = "p",
          method = "spearman", adjust = "BH"
        )
        cor_test$r[is.na(cor_test$r)] <- 0

        rm(sample_pixels)

        # Calculate the sigificance of the correlations
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
  save(cor_matrix, file = here(
    "output",
    "cor_mats_offset",
    TISSUE_OF_INTEREST,
    paste0("cor_mat_", TISSUE_OF_INTEREST, "_", sample_size, "px_", NUM_REPS, "reps",
           ifelse(APPLY_SCALING, "_SCALED", ""), ".RData")
    ))
}
rm(n)
cat("done.\n")

# Correlations offset data ----

cat("Correlations from data with offsets...\n")

# Add the offset to half of the samples
idx.samples.offset <- lapply(idx.samples.of.interest, function(z) {
  sample(z, round(length(z) / 2))
})
save(idx.samples.of.interest, file = here(
  "output",
  "cor_mats_offset",
  TISSUE_OF_INTEREST,
  paste0("idx_", TISSUE_OF_INTEREST, "_offset_samples.RData")
))

# Generate the sampled pixels ----
MU <- c(0.1, 0.5, seq(1, 5))   # offset distribution mean
SIG <- 0.1  # offset distribution standard deviation
num.vars <- ncol(multi_data_env$multi_data$cv[[1]])

# Apply the offset
for (m in 1:length(MU))
{
  cat(paste0("Offset distribution: mu=", MU[m], ", sigma=", SIG, "\n"))

  for (n in 1:length(NUM_PIXELS))
  {
    res.dir <- here(
      "output", "cor_mats_offset", TISSUE_OF_INTEREST,
      paste0("mean_", MU[m], "/sample_size_", NUM_PIXELS[n])
    )
    if (!dir.exists(res.dir)) {
      dir.create(res.dir, recursive = TRUE)
    }

    for (set in "cv")
    {
      for (i in 1:length(idx.samples.offset[[set]]))
      {
        curr.idx <- idx.samples.offset[[set]][i]
        curr.res.dir <- paste0(res.dir, "/", set, "/sample_", curr.idx)
        if (!dir.exists(curr.res.dir)) {
          dir.create(curr.res.dir, recursive = TRUE)
        }
      }
      rm(i)
    }
    rm(set)

    for (rx in 1:NUM_REPS)
    {
      cat(paste0("Sample size: ", NUM_PIXELS[n], ", replicate: ", rx, "\n"))

      for (set in "cv")
      {
        for (i in 1:length(idx.samples.of.interest[[set]]))
        {
          curr.idx <- idx.samples.of.interest[[set]][i]
          curr.res.dir <- paste0(res.dir, "/", set, "/sample_", curr.idx)

          if (curr.idx %in% idx.samples.offset[[set]]) {
            cat(paste0("set: ", set, ", MS image: ", curr.idx, "\n"))
            num.intensities <- num.vars * nrow(multi_data_env$multi_data[[set]][[curr.idx]])
            offset.mat <- matrix(
              rnorm(n = num.intensities, mean = MU[m], sd = SIG),
              nrow(multi_data_env$multi_data[[set]][[curr.idx]]), num.vars
            )
            multi_data_offset <- multi_data_env$multi_data[[set]][[curr.idx]] + offset.mat

            rm(offset.mat, num.intensities)
            gc()

            save(multi_data_offset, file = paste0(curr.res.dir, "/data_offset_rep_", rx, ".RData"))
            rm(multi_data_offset)
          }
          rm(curr.idx, curr.res.dir)
        }
        rm(i)
      }
      rm(set)
    }
    rm(rx)
  }
  rm(n)
}
rm(m)

gc()

for (m in 1:length(MU))
{
  cat(paste0("Offset distribution: mu=", MU[m], ", sigma=", SIG, "\n"))

  # CV and test sets
  cat(paste0("Starting the parallel cluster using ", NUM_CORES, "...\n"))
  cl <- makeCluster(NUM_CORES)

  for (n in 1:length(NUM_PIXELS))
  {
    sample_size <- NUM_PIXELS[n]
    # message(paste0('Sample size: ', sample_size, '\n'))

    res.dir <- here(
      "output", "cor_mats_offset", TISSUE_OF_INTEREST,
      paste0("/mean_", MU[m], "/sample_size_", NUM_PIXELS[n])
    )

    # Correlations are split into two sets: 'cv' and 'test'. Each group contains
    # a list (one element per tissue section) of the replicated correlations
    # (num. replicates x num. correlations).

    cor_matrix <- .vector(mode = "list", length = 1, n = "cv")

    for (set in "cv")
    {
      cor_matrix[[set]] <- vector(mode = "list", length = length(idx.samples.of.interest[[set]]))

      for (i in 1:length(idx.samples.of.interest[[set]]))
      {
        curr.idx <- idx.samples.of.interest[[set]][i]
        cat(paste0("offset: ", MU[m], ", sample size:, ", sample_size,
                   " set: ", set, ", MS image: ", curr.idx, "\n"))

        curr.res.dir <- paste0(res.dir, "/", set, "/sample_", curr.idx)

        intensity_matrix <- multi_data_env$multi_data[[set]][[curr.idx]]
        s <- sample_pixels[[n]][[set]][[i]]

        do.offset <- FALSE
        if (curr.idx %in% idx.samples.offset[[set]]) {
          do.offset <- TRUE
        }

        # Calculate the repeated correlation matrices

        clusterExport(cl = cl, varlist = c("intensity_matrix", "s", "do.offset", "curr.res.dir"))

        cor_matrix_tmp <- parSapply(cl = cl, X = 1:NUM_REPS, FUN = function(X) {
          require(psych)
          if (do.offset) {
            load(paste0(curr.res.dir, "/data_offset_rep_", X, ".RData"), envir = environment())
            cor_test <- corr.test(multi_data_offset[s[X, ], ],
              use = "p",
              method = "spearman", adjust = "BH"
            )
            rm(multi_data_offset)
          } else {
            cor_test <- corr.test(intensity_matrix[s[X, ], ],
              use = "p",
              method = "spearman", adjust = "BH"
            )
          }
          cor_test$r[is.na(cor_test$r)] <- 0

          rm(sample_pixels)

          # Calculate the sigificance of the correlations
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
    save(cor_matrix, file = paste0(
      save.dir, "/cor_mat_", TISSUE_OF_INTEREST, "_",
      sample_size, "px_", NUM_REPS, "reps",
      "_OFFSET_MU_", MU[m], ".RData"
    ))
  }
  rm(n)
  stopCluster(cl = cl)
  closeAllConnections()
  gc()
}
rm(m)
