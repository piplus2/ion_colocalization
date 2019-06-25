# Generate the adjacency matrices for the graphs of the significantly different
# correlations in the 3 contrasts.
#
# Author: Paolo Inglese, Imperial College London

# PACKAGES ----

tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})
pacman::p_load("here", "readr")

# SETUP ----

NUM_REPS <- 500 # Number repetitions
SAMPLE_SIZE <- 900 # Optimal sample size used for calculating the correlations

# SOURCE ----

source(here("R_scripts", "functions", "misc", "load.R"))

# START ----

SETS <- c('cv', 'test')
TISSUES <- c('breast', 'colorectal', 'ovarian')
NUM_TISSUES <- length(TISSUES)

# LOAD ----

table_molecules <- read_csv(here("DATA", "table_molecules.csv"))

cmz_env <- .load(here("DATA", "RData", "cmz_between_samples.RData"))
cmz_env$cmz <- cmz_env$cmz[-3] # Remove 282.2511 m/z manually
num_ions <- length(cmz_env$cmz)
num_cor <- num_ions * (num_ions - 1) / 2 # Number of correlation features

sets_env <- .load(here("DATA", "RData", "split_idx_cv_test.RData"))
num_samples <- c(length(sets_env$idx_cv), length(sets_env$idx_test))
names(num_samples) <- SETS

lbl_env <- .load(here("DATA", "RData", "labels.RData"))
Y_CV <- factor(lbl_env$Y_CV)
Y_TEST <- factor(lbl_env$Y_TEST)

rm(sets_env, lbl_env, cmz_env)
gc()

# Load the correlation matrices using the optimal sample size and calculate the
# average correlation within each sample (mean across the replicates as well)
# and their standard error.

cat('Loading the correlation features...\n')

data_env <- .load(here("output",
                       "cor_mats_rand_sampling",
                       paste0("cor_mat_", SAMPLE_SIZE, 'px_', NUM_REPS, 'reps.RData')))

cor_reps <- array(NA, c(num_samples['cv'], num_cor, NUM_REPS))
for (i in 1:num_samples['cv'])
{
  cat('Sample', i, '\n')
  cor_reps[i, , ] <- t(data_env$cor_matrix[['cv']][[i]])
}
rm(i)

rm(data_env)
gc()

# ADJACENCIES ----

# Calculate the mean of the correlations (mean of tissue sections and repetitions)
# within the 3 tissue types

# First calculate the mean correlation within the tissue in each repetition
cat("Mean correlations within repetition...\n")

cor_class_mean <- array(NA, c(NUM_TISSUES, num_cor, NUM_REPS))
for (i in 1:NUM_TISSUES)
{
  for (j in 1:num_cor)
  {
    for (k in 1:NUM_REPS)
    {
      cor_class_mean[i, j, k] <- mean(cor_reps[Y_CV == TISSUES[i], j, k])
    }
    rm(k)
  }
  rm(j)
}
rm(i)

# Then the grand mean across the repetitions
cat("Mean across repetitions...\n")

cor_grand_mean <- matrix(NA, NUM_TISSUES, num_cor)
cor_sd <- matrix(NA, NUM_TISSUES, num_cor)

for (i in 1:NUM_TISSUES)
{
  for (j in 1:num_cor)
  {
    cor_grand_mean[i, j] <- mean(cor_class_mean[i, j, ])
    cor_sd[i, j] <- sd(cor_class_mean[i, j, ])
  }
  rm(j)
}
rm(i)

save(list = c('cor_grand_mean', 'cor_sd'),
     file = here("output", "univariate_test", "cor_grand_mean_sd.RData"))

# Select the significant correlations in more than 95% of replicates

data_env <- .load(here("output", "univariate_test", "signif_cor_95perc_reps.RData"))

contr <- list(c('breast', 'colorectal'),
              c('breast', 'ovarian'),
              c('colorectal', 'ovarian'))

for (i in 1:length(contr))
{
  count <- 0
  sel_adj_1 <- sel_adj_2 <- matrix(0, num_ions, num_ions,
                                   dimnames = list(table_molecules$Symbol,
                                                   table_molecules$Symbol))
  for (j in 1:num_ions)
  {
    for (k in 1:num_ions)
    {
      if (k > j)
      {
        count <- count + 1
        if (any(data_env$sel_signif_cor[[i]] == count))
        {
          sel_adj_1[j, k] <- cor_grand_mean[TISSUES == contr[[i]][1], count]
          sel_adj_2[j, k] <- cor_grand_mean[TISSUES == contr[[i]][2], count]
        } else
        {
          sel_adj_1[j, k] <- NA
          sel_adj_2[j, k] <- NA
        }
      }
    }
    rm(k)
  }
  rm(j)
  sel_adj_1 <- sel_adj_1 + t(sel_adj_1)
  sel_adj_2 <- sel_adj_2 + t(sel_adj_2)

  disc_edges_1 <- apply(sel_adj_1, 2, function(z) all(is.na(z))) # Fully disconnected nodes
  sel_adj_1 <- sel_adj_1[!disc_edges_1, !disc_edges_1]

  disc_edges_2 <- apply(sel_adj_2, 2, function(z) all(is.na(z))) # Fully disconnected nodes
  sel_adj_2 <- sel_adj_2[!disc_edges_2, !disc_edges_2]

  rm(disc_edges_1, disc_edges_2)

  sel_adj_1[is.na(sel_adj_1)] <- 0
  sel_adj_2[is.na(sel_adj_2)] <- 0

  stopifnot(sum(sel_adj_1 != 0, na.rm = T) / 2 == length(data_env$sel_signif_cor[[i]]))
  stopifnot(sum(sel_adj_2 != 0, na.rm = T) / 2 == length(data_env$sel_signif_cor[[i]]))

  # Save the adjacencies of the two contrasts
  write_csv(as.data.frame(sel_adj_1),
            path = here("output", "univariate_test",
                        paste0("adj_1_signif_cor_contr_", contr[[i]][1], "_", contr[[i]][2], ".csv")))
  write_csv(as.data.frame(sel_adj_2),
            path = here("output", "univariate_test",
                        paste0("adj_2_signif_cor_contr_", contr[[i]][1], "_", contr[[i]][2], ".csv")))

  rm(sel_adj_1, sel_adj_2)
}
rm(i)
