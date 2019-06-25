# Classification using the mean spectrum and offsets


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

MAX_SAMPLE <- 1000
NUM_PIXELS <- seq(100, MAX_SAMPLE, 100)
NUM_REPS <- 100 # Num. of repetitions

APPLY_SCALING <- FALSE # If TRUE it uses the UV-scaled peaks intensities

for (TISSUE_OF_INTEREST in c("breast", "colorectal", "ovarian"))
{
  cat(paste0("Analyzing ", TISSUE_OF_INTEREST, " samples\n"))

  # DIRS/FILES ----

  output_dir <- .setDir(here("offset", "classif", "mean_spectrum_offset", TISSUE_OF_INTEREST),
    createIfNotExist = TRUE
  )
  SETS <- c("cv", "test")
  INPUT_FILENAME <- c(
    cv = "X_matched_between_SPUTNIK.RData",
    test = "X_matched_with_cmz.RData"
  )

  # Cross-validation / test samples MS images

  data_env <- .load("./DATA/RData/split_idx_cv_test.RData")

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

  offset_samples <- list.dirs(here(
    "output",
    "cor_mats_offset",
    TISSUE_OF_INTEREST,
    "mean_0.1/sample_size_100/cv"
  ),
  recursive = FALSE,
  full.names = FALSE
  )

  offset_samples <- as.numeric(
    sapply(offset_samples, function(z) strsplit(z, split = "_")[[1]][2])
  )

  Y <- array(0, length(idx.samples.of.interest))
  Y[idx.samples.of.interest %in% offset_samples] <- 1
  rm(offset_samples)

  Y <- factor(Y)
  print(table(Y))

  # Load the multi_data
  load(here(
    "DATA",
    "RData",
    paste0("/multi_data_", ifelse(APPLY_SCALING, "SCALED", ""), ".RData")
  ))

  ## Load the sampled pixels
  load(here(
    "output",
    "cor_mats_offset",
    TISSUE_OF_INTEREST,
    paste0("sample_pixels_", TISSUE_OF_INTEREST, ".RData")
  ))

  # Classification
  MU <- c(0, 0.1, 0.5, seq(1, 5))
  NUM_VARS <- ncol(multi_data$cv[[1]])

  for (m in 1:length(MU))
  {
    cat("Starting the parallel cluster...\n")
    cl <- makeCluster(NUM_CORES)

    for (ir in 1:NUM_REPS)
    {
      conf_mat_cv <- .vector(
        mode = "list", length = length(NUM_PIXELS),
        names = paste0("npix_", NUM_PIXELS)
      )

      for (n in 1:length(NUM_PIXELS))
      {
        sample_size <- NUM_PIXELS[n]

        message(paste0("Rep: ", ir, ", sample size:", sample_size, "\n"))

        feats_cv <- matrix(NA, num_samples, NUM_VARS)

        for (i in 1:num_samples)
        {
          curr.idx <- idx.samples.of.interest[i]
          s <- sample_pixels[[n]][["cv"]][[i]][ir, ]
          if (MU[m] == 0) {
            feats_cv[i, ] <- apply(multi_data[["cv"]][[curr.idx]][s, ], 2, mean)
          } else {
            if (Y[i] == 1) {
              res.dir <- here(
                "output",
                "cor_mats_offset",
                TISSUE_OF_INTEREST,
                paste0("mean_", MU[m], "/sample_size_", NUM_PIXELS[n])
              )
              curr.res.dir <- paste0(res.dir, "/", "cv", "/sample_", curr.idx)
              load(paste0(curr.res.dir, "/data_offset_rep_", ir, ".RData"))
              # Calculate the mean intensities of the offset data
              feats_cv[i, ] <- apply(multi_data_offset[s, ], 2, mean)
              rm(multi_data_offset)
            } else {
              # Calculate the mean intensities of the original data
              feats_cv[i, ] <- apply(multi_data[["cv"]][[curr.idx]][s, ], 2, mean)
            }
          }
          stopifnot(!any(is.na(feats_cv[i, ])))
        }
        rm(i)

        colnames(feats_cv) <- make.names(1:ncol(feats_cv))

        # Remove the cross-validation constant features
        pre <- preProcess(feats_cv, "zv")
        feats_cv <- predict(pre, feats_cv)
        rm(pre)

        # With a fixed sample size, we run PLS-DA with ncomp = 2:10

        start_time <- proc.time()

        clusterExport(cl = cl, varlist = c("Y", "feats_cv"))

        conf_mat_cv[[n]] <-
          parLapply(
            cl = cl, X = 2:10, # X = PLS-DA num. components
            fun = function(X) {
              require(caret)

              # LOOCV
              yPred <- array(NA, length(Y))

              for (cv in 1:length(Y))
              {
                xTrain <- feats_cv[-cv, , drop = F]
                xTest <- feats_cv[cv, , drop = F]
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

        names(conf_mat_cv[[n]]) <- paste0("ncomps_", 2:10)
        gc()

        diff_time <- proc.time() - start_time
        print(diff_time)
      }

      save(conf_mat_cv, file = paste0(
        output_dir, "/conf_mat_sampling_MU_", MU[m], "_rep_", ir,
        ".RData"
      ))
    }
    rm(ir)

    stopCluster(cl)
    closeAllConnections()

    # Classification performance

    acc_cv <- acc_test <- array(NA, c(length(NUM_PIXELS), 9, NUM_REPS),
      dimnames = list(NUM_PIXELS, NULL, NULL)
    )

    for (ir in 1:NUM_REPS)
    {
      cat("Repetition:", ir, "\n")

      data_cv_env <- new.env()
      load(paste0(
        output_dir, "/conf_mat_sampling_MU_", MU[m], "_rep_", ir,
        ".RData"
      ),
      envir = data_cv_env
      )

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

    acc_cv_mean <- acc_cv_q1 <- acc_cv_q2 <- matrix(NA, length(NUM_PIXELS), 9)

    for (n in 1:length(NUM_PIXELS))
    {
      for (k in 1:9) # k = ncomps - 1
      {
        acc_cv_mean[n, k] <- median(acc_cv[n, k, ])
        acc_cv_q1[n, k] <- quantile(acc_cv[n, k, ], 0.25)
        acc_cv_q2[n, k] <- quantile(acc_cv[n, k, ], 0.75)
      }
      rm(k)
    }
    rm(n)

    write.csv(acc_cv_mean, file = paste0(output_dir, "/median_acc_MU_", MU[m], ".csv"))
    write.csv(acc_cv_q1, file = paste0(output_dir, "/q1_acc_MU_", MU[m], ".csv"))
    write.csv(acc_cv_q2, file = paste0(output_dir, "/q2_acc_MU_", MU[m], ".csv"))
  }
  rm(m)
}
rm(TISSUE_OF_INTEREST)
