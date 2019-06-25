# Classification of the colocalization features from the simulated offset

# PACKAGES ----
tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})
pacman::p_load("caret", "parallel", "here")


# FUNCTIONS ----

resultsClassif <- function(conf_mat) {
  num_classes <- nrow(conf_mat)
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  prec <- rec <- array(NA, num_classes, dimnames = list(colnames(conf_mat)))
  for (i in 1:num_classes)
  {
    tp <- conf_mat[i, i] / sum(conf_mat[, i])
    fp <- sum(conf_mat[i, -i]) / sum(conf_mat[, i])
    fn <- sum(conf_mat[-i, i]) / sum(conf_mat[, i])
    prec[i] <- tp / (tp + fp)
    rec[i] <- tp / (tp + fn)
  }
  return(list("acc" = accuracy, "prec" = prec, "rec" = rec))
}

source(here("R_scripts", "functions", "misc", "load.R"))
source(here("R_scripts", "functions", "misc", "set_dir.R"))
source(here("R_scripts", "functions", "misc", "vector.R"))

# SETUP ----

TISSUE_OF_INTEREST <- "breast" # Run the analysis on this class of samples
# It can be "breast", "colorectal", "ovarian"

NUM_CORES <- 5 # Num. of CPU cores used for the parallel computation -
# It must be smaller than total num. of cores
message(paste0("Using ", NUM_CORES, " cores."))

NUM_COR <- 2016 # num.ions * (num.ions - 1) / 2
MAX_SAMPLE <- 1000 # Largest random sample
NUM_PIXELS <- seq(100, MAX_SAMPLE, 100) # Size of random samples
NUM_REPS <- 100 # Num. of repetitions

# START ----

# DIRS/FILES ----

# Cross-validation / test samples MS images

SETS <- c("cv", "test")
INPUT_FILENAME <- c(
  cv = "X_matched_between_SPUTNIK.RData",
  test = "X_matched_with_cmz.RData"
)

output_dir <- .setDir(here("output", "classif", "cor_mats_offset", TISSUE_OF_INTEREST),
                      createIfNotExist = TRUE
)
rm(output_dir)

data_env <- .load(here("DATA", "RData", "split_idx_cv_test.RData"))

idx_samples <- .vector(mode = "list", length = 2, names = SETS)
idx_samples$cv <- data_env$idx_cv
idx_samples$test <- data_env$idx_test

rm(data_env)
gc()

# Create tissue labels (with/without offset)

# Load the idx of the samples with offset applied (var: idx.samples.of.interest)
load(here(
  "output",
  "cor_mats_offset",
  TISSUE_OF_INTEREST,
  paste0("idx_", TISSUE_OF_INTEREST, "_offset_samples.RData")
))

idx.samples.of.interest <- idx.samples.of.interest$cv
sort(idx.samples.of.interest)

num_samples <- length(idx.samples.of.interest)

# Samples with offset are always the same for the various simulations
# so we can identify them from the folder names of the first simulation
offset_samples <- list.dirs(here(
  "output",
  "cor_mats_offset",
  TISSUE_OF_INTEREST,
  "/mean_0.1/sample_size_100/cv"
), recursive = F, full.names = F)
offset_samples <- as.numeric(
  sapply(offset_samples, function(z) strsplit(z, split = "_")[[1]][2])
)

# Assign a label = 1 to samples with offset
Y <- array(0, length(idx.samples.of.interest))
Y[idx.samples.of.interest %in% offset_samples] <- 1
rm(offset_samples)

Y <- factor(Y)
print(table(Y))  # Check that the two groups are almost of the same size

# Classification ----

# Mean values of the offset distributions, SD is fixed to 0.1
MU <- c(0, 0.1, 0.5, seq(1, 5))

for (m in 1:length(MU))
{
  cat("Starting the parallel cluster...\n")
  cl <- makeCluster(NUM_CORES)

  for (n in 1:length(NUM_PIXELS))
  {
    sample_size <- NUM_PIXELS[n]

    message(paste0("Sample size: ", sample_size, "\n"))

    # Load the correlations

    # Load cor_matrix: list 2 elements ('cv', 'test'), each list element contains
    # a list of matrices (one per sample). Each matrix is NUM_REP x num. correlations

    # If MU = 0, use the original features
    if (MU[m] == 0) {
      data_env <- .load(here(
        "output",
        "cor_mats_sampling",
        "cor_mat_", sample_size, "px_500reps_SCALED.RData"
      ))
      # Select only the tissue-of-interest samples
      data_env$cor_matrix[["cv"]] <- data_env$cor_matrix[["cv"]][idx.samples.of.interest]
    } else {
      data_env <- .load(here(
        "output",
        "cor_mats_offset",
        TISSUE_OF_INTEREST,
        paste0("cor_mat_", TISSUE_OF_INTEREST, "_", sample_size, "px_", NUM_REPS,
               "reps_OFFSET_MU_", MU[m], ".RData")
      ))
    }

    conf_mat_cv <- .vector(mode = "list", length = NUM_REPS, names = paste0("nrep_", 1:NUM_REPS))

    for (ir in 1:NUM_REPS)
    {
      cat(paste0("Offset: ", MU[m], ", sample size: ", sample_size, ", rep: ", ir, "\n"))

      # Setup the correlation features matrices

      cor_feats_cv <- matrix(NA, num_samples, NUM_COR)

      for (i in 1:num_samples)
      {
        cor_feats_cv[i, ] <- data_env$cor_matrix[["cv"]][[i]][ir, ]
        stopifnot(!any(is.na(cor_feats_cv[i, ])))
      }
      rm(i)

      colnames(cor_feats_cv) <- make.names(1:ncol(cor_feats_cv))

      # Remove the cross-validation constant features from both the sets
      pre <- preProcess(cor_feats_cv, "zv")
      cor_feats_cv <- predict(pre, cor_feats_cv)
      rm(pre)

      # With a fixed sample size, we run PLS-DA with ncomp = 2:10

      start_time <- proc.time()

      clusterExport(cl = cl, varlist = c("Y", "cor_feats_cv"))

      conf_mat_cv[[ir]] <-
        parLapply(
          cl = cl, X = 2:10, # X = PLS-DA num. components
          fun = function(X) {
            require(caret)

            # LOOCV
            yPred <- array(NA, length(Y))

            for (cv in 1:length(Y))
            {
              xTrain <- cor_feats_cv[-cv, , drop = F]
              xTest <- cor_feats_cv[cv, , drop = F]
              yTrain <- Y[-cv]

              mdl <- plsda(xTrain, yTrain, ncomp = X)
              yPred[cv] <- as.character(predict(mdl, xTest))
            }

            return(confusionMatrix(
              factor(yPred, levels = levels(Y)),
              Y
            )$table)
          }
        )

      names(conf_mat_cv[[ir]]) <- paste0("ncomps_", 2:10)
      gc()

      diff_time <- proc.time() - start_time
      print(diff_time)
    }
    rm(ir)

    save(conf_mat_cv, file = here(
      "output", "classif", "cor_mats_offset", TISSUE_OF_INTEREST,
      paste0("conf_mat_sampling_MU_", MU[m], "_sample_", NUM_PIXELS[n], "pix.RData")
    ))
  }
  rm(n)

  stopCluster(cl)
  closeAllConnections()
  gc(reset = TRUE)
}
rm(m)
