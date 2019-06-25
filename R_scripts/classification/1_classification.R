# Classify the tissue sections using colocalization features
#
# Author: Paolo Inglese, Imperial College London

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

# Num. of CPU cores used for the parallel computation -
# It must be smaller than total num. of cores
NUM_CORES <- 5
message(paste0("Using ", NUM_CORES, " cores."))

NUM_COR <- 2016  # num.ions * (num.ions - 1) / 2
MAX_SAMPLE <- 1000  # Largest random sample
NUM_PIXELS <- seq(100, MAX_SAMPLE, 100)  # Size of random samples
NUM_REPS <- 500  # Num. of repetitions

APPLY_SCALING <- FALSE  # If TRUE correlations are UV-scaled within tissue section

# START ----

res_dir <- .setDir(here("output", "classif", "cor_mats_rand_sampling"),
                   createIfNotExist = TRUE)
rm(res_dir)

TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)

SETS <- c("cv", "test")
INPUT_FILENAME <- c(
  cv = "X_matched_between_SPUTNIK.RData",
  test = "X_matched_with_cmz.RData"
)

# DIRS/FILES ----

# Cross-validation / test samples MS images

data_env = .load(here("DATA", "RData", "split_idx_cv_test.RData"))

idx_samples <- .vector(mode = "list", length = 2, names = SETS)
idx_samples$cv <- data_env$idx_cv
idx_samples$test <- data_env$idx_test

rm(data_env)
gc()

num_samples <- unlist(lapply(idx_samples, length))

# CLASSIFICATION ----

# Load tissue labels (vars: Y_CV, Y_TEST))
load(here("DATA", "RData", "labels.RData"))

cat(paste0("Starting the parallel cluster using ", NUM_CORES, " cores...\n"))
cl <- makeCluster(NUM_CORES)

for (ir in 1:NUM_REPS)
{
  conf_mat_cv <- conf_mat_test <- .vector(
    mode = "list",
    length = length(NUM_PIXELS),
    names = paste0("npix_", NUM_PIXELS)
  )

  for (n in 1:length(NUM_PIXELS))
  {
    sample_size <- NUM_PIXELS[n]

    message(paste0("Rep: ", ir, ", sample size:", sample_size, "\n"))

    # Load the correlation matrices

    # Load cor_matrix: list 2 elements ('cv', 'test'), each list element contains
    # a list of matrices (one per sample).
    # Each matrix contains NUM_REP x NUM_COR elements
    .data_env <- .load(here(
      "output", "cor_mats_rand_sampling",
      paste0("cor_mat_", sample_size, "px_", NUM_REPS, "reps.RData")
    ))

    # Setup the correlation features matrices

    cor_feats_cv <- matrix(NA, num_samples["cv"], NUM_COR)
    cor_feats_test <- matrix(NA, num_samples["test"], NUM_COR)

    for (i in 1:num_samples["cv"])
    {
      if (APPLY_SCALING)
        data_env$cor_matrix[["cv"]][[i]][ir, ] <- scale(data_env$cor_matrix[["cv"]][[i]][ir, ])
      cor_feats_cv[i, ] <- data_env$cor_matrix[["cv"]][[i]][ir, ]
      stopifnot(!any(is.na(cor_feats_cv[i, ])))
    }
    rm(i)
    for (i in 1:num_samples["test"])
    {
      if (APPLY_SCALING)
        data_env$cor_matrix[["test"]][[i]][ir, ] <- scale(data_env$cor_matrix[["test"]][[i]][ir, ])
      cor_feats_test[i, ] <- data_env$cor_matrix[["test"]][[i]][ir, ]
      stopifnot(!any(is.na(cor_cor_feats_test[i, ])))
    }
    rm(i)

    colnames(cor_feats_cv) <- make.names(1:ncol(cor_feats_cv))
    colnames(cor_feats_test) <- make.names(1:ncol(cor_feats_test))

    rm(data_env)

    # Remove the cross-validation constant features from both the sets
    pre <- preProcess(cor_feats_cv, "zv")
    cor_feats_cv <- predict(pre, cor_feats_cv)
    cor_feats_test <- predict(pre, cor_feats_test)
    rm(pre)

    # With a fixed sample size, we run PLS-DA with ncomp = 2:10

    start_time <- proc.time()

    clusterExport(cl = cl, varlist = c("Y_CV", "cor_feats_cv"))

    conf_mat_cv[[n]] <-
      parLapply(
        cl = cl, X = 2:10, # X = PLS-DA num. components
        fun = function(X) {
          require(caret)

          # 10-fold cross-validation

          cv_folds <- createFolds(Y_CV, 10)
          y_pred <- array(NA, length(Y_CV))

          for (cv in 1:10)
          {
            xTrain <- cor_feats_cv[-cv_folds[[cv]], , drop = F]
            xTest <- cor_feats_cv[cv_folds[[cv]], , drop = F]
            yTrain <- Y_CV[-cv_folds[[cv]]]

            mdl <- plsda(xTrain, yTrain, ncomp = X)
            y_pred[cv_folds[[cv]]] <- as.character(predict(mdl, xTest))
          }

          # Calculate the confusion matrix over all predictions
          return(confusionMatrix(
            factor(y_pred, levels = levels(Y_CV)),
            Y_CV
          )$table)
        }
      )

    names(conf_mat_cv[[n]]) <- paste0("ncomps_", 2:10)
    gc()

    diff_time <- proc.time() - start_time
    print(diff_time)

    # Predictions on the test set using all cv features

    conf_mat_test[[n]] <- lapply(1:10, function(z) return(matrix(
        NA, NUM_TISSUES,
        NUM_TISSUES
      )))
    names(conf_mat_test[[n]]) <- paste0("ncomps_", 1:10)

    for (k in 2:10)
    {
      mdl <- plsda(cor_feats_cv, Y_CV, ncomp = k)
      y_pred <- as.character(predict(mdl, cor_feats_test))
      conf_mat_test[[n]][[k]] <-
        confusionMatrix(
          factor(y_pred, levels = levels(Y_CV)),
          factor(Y_TEST, levels = levels(Y_CV))
        )$table
    }
    rm(k)

    conf_mat_test[[n]] <- conf_mat_test[[n]][-1] # Remove the element corresponding to ncomp = 1
  }
  rm(n)

  save(conf_mat_cv, file = here(
    "output", "classif", "cor_mats_rand_sampling", "conf_mat_sampling_cv_", ir,
    ".RData"
  ))
  save(conf_mat_test, file = here(
    "output", "classif", "cor_mats_rand_sampling", "conf_mat_sampling_test_", ir,
    ".RData"
  ))
}
rm(ir)

stopCluster(cl = cl)
closeAllConnections()

# Calculate the macro-average of the confusion matrices

acc_cv <- acc_test <- array(NA, c(length(NUM_PIXELS), 9, NUM_REPS),
  dimnames = list(NUM_PIXELS, NULL, NULL)
)

for (ir in 1:NUM_REPS)
{
  cat("Repetition:", ir, "\n")

  data_cv_env <- .load(here(
    "output", "classif", "cor_mats_rand_sampling",
    paste0("conf_mat_sampling_cv_", ir, ".RData")
  ))

  data_test_env <- .load(here(
    "output", "classif", "cor_mats_rand_sampling",
    paste0("conf_mat_sampling_test_", ir, ".RData")
  ))

  for (n in 1:length(NUM_PIXELS))
  {
    for (k in 1:9) # k = ncomps - 1
    {
      res_cv <- resultsClassif(data_cv_env$conf_mat_cv[[n]][[k]])
      res_test <- resultsClassif(data_test_env$conf_mat_test[[n]][[k]])

      acc_cv[n, k, ir] <- res_cv$acc
      acc_test[n, k, ir] <- res_test$acc
    }
    rm(k)
  }
  rm(data_cv_env, data_test_env, n)
}
rm(ir)

acc_cv_mean <- acc_cv_sd <- matrix(NA, length(NUM_PIXELS), 9)
acc_test_mean <- acc_test_sd <- matrix(NA, length(NUM_PIXELS), 9)

for (n in 1:length(NUM_PIXELS))
{
  for (k in 1:9) # k = ncomps - 1
  {
    acc_cv_mean[n, k] <- mean(acc_cv[n, k, ])
    acc_cv_sd[n, k] <- sd(acc_cv[n, k, ])

    acc_test_mean[n, k] <- mean(acc_test[n, k, ])
    acc_test_sd[n, k] <- sd(acc_test[n, k, ])
  }
  rm(k)
}
rm(n)

write.csv(acc_cv_mean, file = here("output", "classif", "cor_mats_rand_sampling", "mean_acc_cv.csv"))
write.csv(acc_cv_sd, file = here("output", "classif", "cor_mats_rand_sampling", "sd_acc_cv.csv"))
write.csv(acc_test_mean, file = here("output", "classif", "cor_mats_rand_sampling", "mean_acc_test.csv"))
write.csv(acc_test_sd, file = here("output", "classif", "cor_mats_rand_sampling", "sd_acc_test.csv"))

gc()
