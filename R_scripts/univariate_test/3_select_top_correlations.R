# Select the top-N correlations
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
TOP_N <- 50 # Number of correlations to select

# SOURCE ----

source(here("R_scripts", "functions", "misc", "load.R"))

# START ----

SETS <- c("cv", "test")
TISSUES <- c("breast", "colorectal", "ovarian")
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

# Load the significant correlations in at least 95% of the repetitions

sel_env <- .load(here("output", "univariate_test", "signif_cor_95perc_reps.RData"))
p_env <- .load(here("output", "univariate_test", "p_values_univ_kw_test.RData"))

num_signif_vars <- unlist(lapply(sel_env$sel_signif_cor, length))

contr <- list(
  c("breast", "colorectal"),
  c("breast", "ovarian"),
  c("colorectal", "ovarian")
)
num_contr <- length(contr)

contr_strings <- unlist(lapply(contr, function(z) return(paste0(z[1], "_", z[2]))))

# SELECT CORRELATIONS ----

# Calculate the average p-value of the selected correlations in the 500
# repetitions.
# p_env$p_adj is a matrix (number of repetitions x number of contrasts x number
# of correlations).

significance <- -log10(p_env$p_adj)
mean_signif <- sd_signif <- vector(mode = "list", length = NUM_TISSUES)

for (i in 1:num_contr)
{
  mean_signif[[i]] <- sd_signif[[i]] <- array(NA, num_signif_vars[i])
  for (j in 1:num_signif_vars[i])
  {
    mean_signif[[i]][j] <- mean(significance[, i, sel_env$sel_signif_cor[[i]][j]])
    sd_signif[[i]][j] <- sd(significance[, i, sel_env$sel_signif_cor[[i]][j]])
  }
}

# Select the N correlations with the largest mean significance across the
# repetitions.

s <- matrix(NA, num_contr, TOP_N)
for (i in 1:num_contr)
{
  s[i, ] <- sort(mean_signif[[i]], decreasing = TRUE, index.return = TRUE)$ix[1:TOP_N]
}
rm(i)

top_idx <- matrix(NA, num_contr, TOP_N)
for (i in 1:num_contr)
{
  top_idx[i, ] <- sel_env$sel_signif_cor[[i]][s[i, ]]
  mean_signif[[i]] <- mean_signif[[i]][s[i, ]]
  sd_signif[[i]] <- sd_signif[[i]][s[i, ]]
}
rm(i, s)

save(top_idx, file = here("output", "univariate_test", paste0("top_", TOP_N, "_cor_indices.RData")))

# Print the molecules corresponding to the top N correlations

# Grand mean of the top N correlations

cor_env <- .load(here("output", "univariate_test", "cor_grand_mean_sd.RData"))

top_cor_mean <- lapply(1:num_contr, function(z) {
  idx <- which(TISSUES %in% contr[[z]])
  tb <- rbind(
    cor_env$cor_grand_mean[idx[1], top_idx[z, ]],
    cor_env$cor_grand_mean[idx[2], top_idx[z, ]]
  )
  dimnames(tb) <- list(TISSUES[idx], paste0("top_", 1:TOP_N, "_cor"))
  return(tb)
})
names(top_cor_mean) <- contr_strings

save(top_cor_mean, file = here(
  "output", "univariate_test", paste0(
    "top_", TOP_N, "_signif_cor_",
    "values.RData"
  )
))

# First generate a matrix with all the molecule symbols
mol_names_mat <- matrix(NA, num_ions, num_ions)
for (i in 1:num_ions)
{
  for (j in 1:num_ions)
  {
    if (j > i) {
      mol_names_mat[i, j] <- paste0(
        table_molecules$Symbol[i], ",",
        table_molecules$Symbol[j]
      )
    }
  }
  rm(j)
}
rm(i)

# Vector of all the used pairs of molecules
mol_names_vec <- mol_names_mat[upper.tri(mol_names_mat, diag = FALSE)]
rm(mol_names_mat)

# Generate the final tables

for (i in 1:length(contr))
{
  fold_change <- top_cor_mean[[i]][1, ] / top_cor_mean[[i]][2, ]

  list_mol_pairs <- c()
  for (j in 1:ncol(top_idx))
  {
    list_mol_pairs <- rbind(
      list_mol_pairs,
      strsplit(mol_names_vec[top_idx[i, j]],
        split = ","
      )[[1]]
    )
  }
  rm(j)

  table_top_cor <- data.frame(
    mol1 = list_mol_pairs[, 1],
    mol2 = list_mol_pairs[, 2],
    fold_change = fold_change,
    mean_signif = mean_signif[[i]],
    sd_signif = sd_signif[[i]]
  )

  rm(list_mol_pairs, fold_change)

  write.csv(table_top_cor,
    file = here(
      "output", "univariate_test", paste0(
        "table_top_", TOP_N, "_cor_",
        contr[[i]][1], "_", contr[[i]][2], ".csv"
      )
    )
  )
}
rm(i)
