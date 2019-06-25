# PACKAGES ----
tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})

pacman::p_load("here")

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

# START ----

NUM_REPS <- 500 # Num. of repetitions
NUM_PIXELS <- seq(100, 1000, 100)
MAX_PLSDA_COMPS <- 10

TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)
SETS <- c("cv", "test")
INPUT_FILENAME <- c(
  cv = "X_matched_between_SPUTNIK.RData",
  test = "X_matched_with_cmz.RData"
)

RESULTS_DIR <- .setDir(here("output", "cor_mats_rand_sampling"),
  createIfNotExist = FALSE
)

# ACCURACIES ----

# Calculate accuracies for CV set

acc.cv <- array(NA, c(MAX_PLSDA_COMPS - 1, length(NUM_PIXELS), NUM_REPS))

for (rx in 1:NUM_REPS)
{
  cat(paste0("Repetition: ", rx, "/", NUM_REPS, "\n"))
  conf.mat.env <- .load(paste0(RESULTS_DIR, "/conf_mat_sampling_cv", rx, ".RData"))

  for (n in 1:length(NUM_PIXELS))
  {
    for (k in 2:10) # k = ncomps
    {
      acc.cv[k - 1, n, rx] <- resultsClassif(conf.mat.env$conf_mat_cv[[n]][[k - 1]])$acc
    }
    rm(k)
  }
  rm(n)
}

acc.tab <- c()

for (k in 2:10)
{
  acc.tmp <- c()
  for (n in 1:length(NUM_PIXELS))
  {
    acc.med <- mean(acc.cv[k - 1, n, ])
    acc.q1 <- quantile(acc.cv[k - 1, n, ], probs = 0)
    acc.q2 <- quantile(acc.cv[k - 1, n, ], probs = 1)

    acc.tmp <- c(acc.tmp, c(acc.med, acc.q1, acc.q2))
  }
  rm(n)
  acc.tab <- rbind(acc.tab, acc.tmp)
}
rm(k)

acc.tab <- data.frame(acc.tab)
write.csv(acc.tab, file = paste0(RESULTS_DIR, "/acc_cv.csv"))


acc.test <- array(NA, c(MAX_PLSDA_COMPS - 1, length(NUM_PIXELS), NUM_REPS))

for (rx in 1:NUM_REPS)
{
  cat(paste0("Repetition: ", rx, "/", NUM_REPS, "\n"))
  conf.mat.env <- .load(paste0(RESULTS_DIR, "/conf_mat_sampling_test", rx, ".RData"))

  for (n in 1:length(NUM_PIXELS))
  {
    for (k in 2:10)
    {
      acc.test[k - 1, n, rx] <- resultsClassif(conf.mat.env$conf_mat_test[[n]][[k - 1]])$acc
    }
  }
}

acc.tab <- c()

for (k in 2:10)
{
  acc.tmp <- c()
  for (n in 1:length(NUM_PIXELS))
  {
    acc.med <- mean(acc.test[k - 1, n, ])
    acc.q1 <- quantile(acc.test[k - 1, n, ], probs = 0)
    acc.q2 <- quantile(acc.test[k - 1, n, ], probs = 1)

    acc.tmp <- c(acc.tmp, c(acc.med, acc.q1, acc.q2))
  }
  acc.tab <- rbind(acc.tab, acc.tmp)
}

acc.tab <- data.frame(acc.tab)
write.csv(acc.tab, file = paste0(RESULTS_DIR, "/acc_test.csv"))
