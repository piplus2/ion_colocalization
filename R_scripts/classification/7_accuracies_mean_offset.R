## Calculate the median and the 1st, 3rd quartile of the accuracies of the
## mean peaks intensity from the offset data

# PACKAGES ----

tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})

pacman::p_load("here")

# FUNCTIONS ----

resultsClassif <- function(conf_mat)
{
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
  return(list('acc' = accuracy, 'prec' = prec, 'rec' = rec))
}

source(here("R_scripts", "functions", "misc", "load.R"))
source(here("R_scripts", "functions", "misc", "set_dir.R"))
source(here("R_scripts", "functions", "misc", "vector.R"))

## Calculate the performances
MAX_SAMPLE <- 1000
NUM_PIXELS <- seq(100, MAX_SAMPLE, 100)
NUM_REPS <- 100 # Num. of repetitions
MU = c(0, 0.1, 0.5, seq(1, 5)) # Mean values of the offset distributions

# START ----

SETS <- c('cv', 'test')
INPUT_FILENAME <- c(cv = 'X_matched_between_SPUTNIK.RData',
                    test = 'X_matched_with_cmz.RData')


for (TISSUE_OF_INTEREST in c("breast", "colorectal", "ovarian"))
{
  message(paste0("Analyzing ", TISSUE_OF_INTEREST, " samples."))

  output_dir = .setDir(here("output", "classif", "mean_spectrum_offset", TISSUE_OF_INTEREST),,
                       createIfNotExist = TRUE)

  # DIRS/FILES ----

  # Cross-validation / test samples MS images
  data_env <- .load(here("DATA", "RData", "split_idx_cv_test.RData"))

  idx_samples <- .vector(mode = 'list', length = 2, names = SETS)
  idx_samples$cv <- data_env$idx_cv
  idx_samples$test <- data_env$idx_test

  rm(data_env)
  gc()

  # Load the idx of the samples with offset applied
  load(here("output",
            "cor_mats_offset",
            TISSUE_OF_INTEREST,
            paste0('idx_', TISSUE_OF_INTEREST, '_offset_samples.RData')))
  idx.samples.of.interest = idx.samples.of.interest$cv
  sort(idx.samples.of.interest)

  num_samples = length(idx.samples.of.interest)

  for (m in 1:length(MU))
  {
    acc_cv <- acc_test <- array(NA, c(length(NUM_PIXELS), 9, NUM_REPS),
                                dimnames = list(NUM_PIXELS, NULL, NULL))

    for (ir in 1:NUM_REPS)
    {
      cat('Repetition:', ir, '\n')

      data_cv_env <- new.env()
      data_cv_env <- .load(paste0(output_dir,
                                  '/conf_mat_sampling_MU_', MU[m], '_rep_', ir,
                                  '.RData'))

      for (n in 1:length(NUM_PIXELS))
      {
        for (k in 1:9) # k = ncomps - 1
        {
          res_cv <- resultsClassif(data_cv_env$conf_mat_cv[[n]][[k]])
          acc_cv[n, k, ir] <- res_cv$acc
        }
        rm(k)
      }

      rm(data_cv_env)
    }
    rm(ir)

    acc_cv_median <- acc_cv_q1 <- acc_cv_q2 <- matrix(NA, length(NUM_PIXELS), 9)

    for (n in 1:length(NUM_PIXELS))
    {
      for (nc in 1:9) # nc = ncomps - 1
      {
        acc_cv_median[n, nc] <- median(acc_cv[n, nc, ])
        acc_cv_q1[n, nc] <- quantile(acc_cv[n, nc, ], 0.25)
        acc_cv_q2[n, nc] = quantile(acc_cv[n, nc, ], 0.75)
      }
    }

    write.csv(acc_cv_median, file = paste0(output_dir, '/median_acc_MU_', MU[m], '.csv'))
    write.csv(acc_cv_q1, file = paste0(output_dir, '/q1_acc_MU_', MU[m], '.csv'))
    write.csv(acc_cv_q2, file = paste0(output_dir, '/q2_acc_MU_', MU[m], '.csv'))
  }
  rm(m)
}
rm(TISSUE_OF_INTEREST)

pacman::p_unload("here")
