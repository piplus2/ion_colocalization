# Classification using the mean peaks intensities

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

# Num. of CPU cores used for the parallel computation
# It must be smaller than total num. of cores
NUM_CORES <- 5

message(paste0("Using ", NUM_CORES, " cores."))

NUM_REPS <- 500 # Num. of repetitions

# START ----

TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)
SETS <- c("cv", "test")
INPUT_FILENAME <- c(
  cv = "X_matched_between_SPUTNIK.RData",
  test = "X_matched_with_cmz.RData"
)

output_dir <- .setDir(here("output", "classif", "mean_spectrum"), createIfNotExist = T)
rm(output_dir)

# Cross-validation / test MS images

data_env <- .load(here("DATA", "RData", "split_idx_cv_test.RData"))

idx_samples <- .vector(mode = "list", length = 2, names = SETS)
idx_samples$cv <- data_env$idx_cv
idx_samples$test <- data_env$idx_test

rm(data_env)
gc()

num_samples <- unlist(lapply(idx_samples, length))
print(num_samples)

# MULTI DATA ----

# Load the multi_data and labels
multi_data_env <- .load(here("DATA", "RData", "multi_data.RData"))
load(here("DATA", "RData", "labels.RData"))

# Calculate the mean spectra
mean_spectrum <- .vector(mode = "list", length = length(SETS), names = SETS)
num_ions <- ncol(multi_data$cv[[1]])
for (set in SETS)
{
  mean_spectrum[[set]] <- matrix(NA, length(idx_samples[[set]]), num_ions,
    dimnames = list(NULL, make.names(1:num_ions))
  )
  for (i in 1:length(idx_samples[[set]]))
  {
    cat(paste0("Set: ", set, ", sample: ", i, "\n"))
    mean_spectrum[[set]][i, ] <- apply(multi_data_env$multi_data[[set]][[i]], 2, mean)
  }
  rm(i)
}
rm(set)

# CLASSIFICATION ----

conf_mat_cv <- .vector(mode = "list", length = NUM_REPS, names = paste0("rep_", 1:NUM_REPS))
cat("Starting the parallel cluster...\n")
cl <- makeCluster(NUM_CORES)

for (rx in 1:NUM_REPS)
{
  cat(paste0("Repetition:", rx, "/", NUM_REPS, "\n"))

  # Remove the cross-validation constant features from both the sets
  feats_cv <- mean_spectrum$cv
  pre <- preProcess(feats_cv, "zv")
  feats_cv <- predict(pre, feats_cv)
  rm(pre)

  # With a fixed sample size, we run PLS-DA with ncomp = 2:10

  start_time <- proc.time()

  clusterExport(cl = cl, varlist = c("Y_CV", "feats_cv"))

  conf_mat_cv[[rx]] <-
    parLapply(
      cl = cl, X = 2:10, # X = PLS-DA num. components
      fun = function(X) {
        require(caret)

        # 10-fold cross-validation

        cvFolds <- createFolds(Y_CV, 10)
        yPred <- array(NA, length(Y_CV))

        for (cv in 1:10)
        {
          xTrain <- feats_cv[-cvFolds[[cv]], , drop = F]
          xTest <- feats_cv[cvFolds[[cv]], , drop = F]
          yTrain <- Y_CV[-cvFolds[[cv]]]

          mdl <- plsda(xTrain, yTrain, ncomp = X)
          yPred[cvFolds[[cv]]] <- as.character(predict(mdl, xTest))
        }

        return(confusionMatrix(
          factor(yPred, levels = levels(Y_CV)),
          Y_CV
        )$table)
      }
    )

  names(conf_mat_cv[[rx]]) <- paste0("ncomps_", 2:10)
  gc()

  diff_time <- proc.time() - start_time
  print(diff_time)
}
rm(rx)

cat("Saving...\n")
save(conf_mat_cv, file = here("output",
                              "classif",
                              "mean_spectrum",
                              paste0("conf_mats_", NUM_REPS, "_reps.RData")))

closeAllConnections()

# Classification performance

acc.cv <- matrix(NA, 9, NUM_REPS,
  dimnames = list(
    paste0("comps_", c(2:10)),
    paste0("rep_", c(1:NUM_REPS))
  )
)
for (rx in 1:NUM_REPS)
{
  cat(paste0("Repetition: ", rx, "\n"))

  for (k in 1:9)  # k = ncomp - 1
  {
    res.cv <- resultsClassif(conf_mat_cv[[rx]][[k]])
    acc.cv[k, rx] <- res.cv$acc
  }
  rm(k)
}
rm(rx)

acc_cv_med <- apply(acc.cv, 1, median)
acc_cv_q1 <- apply(acc.cv, 1, function(z) quantile(z, 0.25))
acc_cv_q2 <- apply(acc.cv, 1, function(z) quantile(z, 0.75))

write.csv(acc.cv, file = here("output",
                              "classif",
                              "mean_spectrum",
                              paste0("accuracies_cv_", NUM_REPS, "reps.csv")))

# Test set
conf_mat_test <- .vector(mode = "list", length = 10, names = paste0("comps_", 1:10))
acc.test <- array(NA, 10)

# Remove the constant features
feats_cv <- mean_spectrum$cv
pre <- preProcess(feats_cv, "zv")
feats_cv <- predict(pre, feats_cv)
feats_test <- predict(pre, mean_spectrum$test)
rm(pre)

for (k in 2:10)
{
  mdl <- plsda(feats_cv, Y_CV, ncomp = k)
  yPred <- as.character(predict(mdl, feats_test))
  conf_mat_test[[k]] <-
    confusionMatrix(
      factor(yPred, levels = levels(Y_CV)),
      factor(Y_TEST, levels = levels(Y_CV))
    )$table
  acc.test[k] <- resultsClassif(conf_mat_test[[k]])$acc
}
rm(k)

write.csv(acc.test, file = here("output",
                                "classif",
                                "mean_spectrum",
                                paste0("accuracies_test.csv")))
