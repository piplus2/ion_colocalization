# Determine the significant correlations for each contrast using the univariate
# Kruskal-Wallis test
#
# Author: Paolo Inglese, Imperial College London

# PACKAGES ----

tryCatch({
  library("pacman")
}, error = function(e) {
  install.packages("pacman")
})

pacman::p_load("here", "ggplot2", "RColorBrewer", "reshape", "tikzDevice")

# SETUP ----

if (!dir.exists(here("output", "univariate_test")))
  dir.create(here("output", "univariate_test"))
if (!dir.exists(here("plots", "univariate_test")))
  dir.create(here("plots", "univariate_test"))

SETS <- c("cv", "test")
TISSUES <- c("breast", "colorectal", "ovarian")
NUM_TISSUES <- length(TISSUES)

NUM_REPS <- 500 # Number repetitions
SAMPLE_SIZE <- 900 # Optimal sample size used for calculating the correlations

# SOURCE ----

source(here("R_scripts", "functions", "misc", "load.R"))
source(here("R_scripts", "functions", "misc", "vector.R"))

# LOAD ----

# Load the common m/z vector

cmz_env <- .load(here("DATA", "RData", "cmz_between_samples.RData"))
cmz_env$cmz <- cmz_env$cmz[-3] # Remove 282.2511 m/z manually
num_ions <- length(cmz_env$cmz)

num_cor <- num_ions * (num_ions - 1) / 2 # Total number of correlation features

sets_env <- .load(here("DATA", "RData", "split_idx_cv_test.RData"))
num_samples <- c(length(sets_env$idx_cv), length(sets_env$idx_test))
names(num_samples) <- SETS

lbl_env <- .load(here("DATA", "RData", "labels.RData"))
Y_CV <- factor(lbl_env$Y_CV)
Y_TEST <- factor(lbl_env$Y_TEST)

rm(sets_env, cmz_env, lbl_env)
gc()

# UNIVARIATE ----

# Univariate test: Kruskal-Wallis test is applied to all the correlations for
# all the binary contrasts. The p-values are corrected within each repetition
# using Benjamini-Hochberg method.

# Load the correlation features for each repetition, using only the optimal
# sample size

cat("Loading the correlation features...\n")

data_env <- .load(here(
  "output", "cor_mats_rand_sampling",
  paste0("cor_mat_", SAMPLE_SIZE, "px_", NUM_REPS, "reps.RData")
))

cor_reps <- array(NA, c(num_samples["cv"], num_cor, NUM_REPS))
for (i in 1:num_samples["cv"])
{
  cat("Sample", i, "\n")
  cor_reps[i, , ] <- t(data_env$cor_matrix[["cv"]][[i]])
}
rm(i)

rm(data_env)
gc()

# Contrasts
contr <- list(
  c("breast", "colorectal"),
  c("breast", "ovarian"),
  c("colorectal", "ovarian")
)
num_contr <- length(contr)

p_raw <- p_adj <- array(NA, c(num_contr, num_cor, NUM_REPS))

for (r in 1:NUM_REPS)
{
  cat("Rep:", r, "\n")

  for (i in 1:num_contr)
  {
    cat(".")
    # Select the tissue sections of the current contrast
    ix <- which(Y_CV %in% contr[[i]])
    # Kruskal-Wallis test
    p_raw[i, , r] <- sapply(1:num_cor, function(z)
      kruskal.test(x = cor_reps[ix, z, r], g = factor(Y_CV[ix]))$p.value)
    rm(ix)
  }
  cat("\n")
  rm(i)

  p_adj[, , r] <- p.adjust(p_raw[, , r], method = "BH") # Num. elements = 3 x 500

  # Print the number of significant correlations
  for (i in 1:num_contr)
    print(sum(p_adj[i, , r] < 0.05))
  rm(i)
}
rm(r)

save(
  list = c("p_raw", "p_adj"),
  file = here("output", "univariate", "p_values_univ_kw_test.rda")
)

# Plot histograms of number of significant correlations across the replicates

# Prepare a matrix with counts of significant correlations in the 3 constrasts
# across the 500 repetitions
num_signif_vars <- matrix(NA, length(contr), NUM_REPS)
for (i in 1:num_contr)
{
  for (j in 1:NUM_REPS)
  {
    num_signif_vars[i, j] <- sum(p_adj[i, , j] < 0.05)
  }
  rm(j)
}
rm(i)

contr_strings <- c(
  "breast, colorectal",
  "breast, ovarian",
  "colorectal, ovarian"
)

df_num_signif <- melt(num_signif_vars)

for (i in 1:length(contr))
  df_num_signif$X1[df_num_signif$X1 == i] <- contr_strings[i]

df_num_signif$X1 <- factor(df_num_signif$X1)

options(
  tikzLatexPackages =
    c(getOption("tikzLatexPackages"), "\\usepackage{fixltx2e}")
)
tikz(here("plots", "univariate_test", "num_signif_cor.tex"),
  width = 8, height = 6, standAlone = TRUE
)

p <- ggplot(df_num_signif, aes(x = value, group = X1, fill = X1)) +
  geom_histogram(alpha = 0.5, bins = 50, color = "black") +
  scale_fill_brewer(name = "Contrast", palette = "Set1") +
  theme_bw(base_size = 14) +
  xlab("Num. significant correlations") +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
print(p)
dev.off()

system(paste0("pdflatex ", here("plots", "univariate_test", "num_signif_cor.tex")))

rm(p)

# Select the correlations that were significant in more than 95% of the
# repetitions.

# First count how many times the correlation was significant in the repetitions
count_signif_mat <- matrix(0, length(contr), num_cor)
for (i in 1:length(contr))
{
  for (j in 1:NUM_REPS)
  {
    count_signif_mat[i, p_adj[i, , j] < 0.05] <-
      count_signif_mat[i, p_adj[i, , j] < 0.05] + 1
  }
  rm(j)
}
rm(i)
saveRDS(count_signif_mat,
  file = here("output", "univariate_test", "count_signif_cor_reps.rds")
)
# Select those significant in more than 95% of times
sel_signif_cor <- .vector(mode = "list", length = length(contr), names = contr)
for (i in 1:length(contr))
{
  sel_signif_cor[[i]] <- which(count_signif_mat[i, ] >= 0.95 * NUM_REPS)
}
rm(i)
saveRDS(sel_signif_cor,
  file = here("output", "univariate_test", "signif_cor_95perc_reps.rds")
)

cat("Done.\n")

pacman::p_unload("here", "ggplot2", "RColorBrewer", "reshape", "tikzDevice")
