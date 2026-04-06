#!/usr/bin/env Rscript

required_packages <- c("brms", "ape", "posterior")
missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  message(
    "Error: Missing required R packages: ",
    paste(missing_packages, collapse = ", ")
  )
  message("Install with:")
  message("  install.packages(c('brms', 'ape', 'posterior'))")
  quit(status = 1)
}

for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
}

set.seed(42) # Same seed for reproducibility

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  message("Error: Exactly 5 arguments required")
  message(
    "Usage: Rscript create_brm.R <tree.nwk> <igs_summary.tsv> <summary.txt> <model.rds> <results_row.tsv>" # nolint
  )
  quit(status = 1)
}

tree_path <- args[1]
igs_path <- args[2]
summary_name <- args[3]
out_name <- args[4]
results_name <- args[5]
factor_levels <- c("same", "opposite")

if (!file.exists(tree_path)) {
  stop("File not found: ", tree_path)
}
if (!file.exists(igs_path)) {
  stop("File not found: ", igs_path)
}

message("Loading tree from: ", tree_path)
message("Loading igs data from: ", igs_path)

tree <- read.tree(tree_path)
tree <- multi2di(tree)

igr <- read.csv(igs_path, sep = "\t")
igr$taxon_tree <- factor(igr$taxon_tree)
igr$polarity_bin <- factor(igr$polarity_bin, levels = factor_levels)

# subsample
max_rows <- 500000
if (nrow(igr) > max_rows) {
  total_n <- nrow(igr)
  message(sprintf("Subsampling %d -> %d rows", total_n, max_rows))
  igr <- do.call(
    rbind,
    lapply(
      split(igr, igr$polarity_bin),
      function(grp) {
        n <- max(1, round(max_rows * nrow(grp) / total_n))
        grp[sample(nrow(grp), min(n, nrow(grp))), ]
      }
    )
  )
  igr$taxon_tree <- droplevels(igr$taxon_tree)
}

tree <- keep.tip(tree, levels(igr$taxon_tree))
tree <- multi2di(tree)
A <- vcv(tree, corr = FALSE) # nolint

n_cores <- min(4, max(1, parallel::detectCores(logical = FALSE) - 1))
message("Running brm with ", n_cores, " cores")

fit <- brm(
  log10_length ~ polarity_bin + (1 | gr(taxon_tree, cov = A)),
  data = igr,
  data2 = list(A = A),
  family = gaussian(),
  cores = n_cores,
  chains = 4,
  save_pars = save_pars(group = FALSE),
)

message("Saving brm model to: ", out_name)
saveRDS(fit, file = out_name)

ps <- posterior_summary(fit)
b <- ps["b_polarity_binopposite", "Estimate"]
lo <- ps["b_polarity_binopposite", "Q2.5"]
hi <- ps["b_polarity_binopposite", "Q97.5"]

message("Saving summary to: ", summary_name)
sink(summary_name)
cat("polarity_bin levels (first = baseline):\n")
print(levels(igr$polarity_bin))
cat("baseline polarity_bin:\n")
print(levels(igr$polarity_bin)[1])
cat("\n--- brms summary ---\n")
print(summary(fit))
sprintf("Delta PGLMM = %.3f (%.3f - %.3f)", b, lo, hi)
sprintf("Fold  PGLMM = %.3f (%.3f - %.3f)", 10^b, 10^lo, 10^hi)
sink()

# ---------- Results ----------
nd <- data.frame(polarity_bin = factor(factor_levels), levels = factor_levels)
mu <- fitted(fit, newdata = nd, re_formula = NA, summary = TRUE)[, c(
  "Estimate",
  "Q2.5",
  "Q97.5"
)]

mu_bp <- data.frame(
  polarity_bin = nd$polarity_bin,
  mean_bp = 10^mu[, "Estimate"],
  lo_bp = 10^mu[, "Q2.5"],
  hi_bp = 10^mu[, "Q97.5"]
)

# Raw medians from observed data (bp scale)
igr$length_bp <- 10^igr$log10_length
med_same <- median(igr$length_bp[igr$polarity_bin == "same"], na.rm = TRUE)
med_opp <- median(igr$length_bp[igr$polarity_bin == "opposite"], na.rm = TRUE)

# Fixed effect
fx <- fixef(fit)
if (!("polarity_binopposite" %in% rownames(fx))) {
  stop(
    "Expected coefficient 'polarity_binopposite' not found. Check polarity_bin levels/coding." # nolint
  )
}

beta_est <- fx["polarity_binopposite", "Estimate"]
beta_lo <- fx["polarity_binopposite", "Q2.5"]
beta_hi <- fx["polarity_binopposite", "Q97.5"]

# Fold-change (opposite / same) on bp scale
fold_est <- 10^beta_est
fold_lo <- 10^beta_lo
fold_hi <- 10^beta_hi

dr <- as_draws_df(fit)
p_beta_gt0 <- mean(dr$b_polarity_binopposite > 0)

sd_phylo <- as.numeric(VarCorr(fit)$taxon_tree$sd[1])
sigma_res <- as.numeric(VarCorr(fit)$residual__$sd[1])

get_mu <- function(pol) {
  mu_bp[mu_bp$polarity_bin == pol, ]
}

res <- data.frame(
  N_regions = nrow(igr),
  N_taxa = nlevels(igr$taxon_tree),

  median_same_bp = med_same,
  median_opposite_bp = med_opp,
  delta_median_bp = med_opp - med_same,

  mean_same_bp = get_mu("same")$mean_bp,
  lo_same_bp = get_mu("same")$lo_bp,
  hi_same_bp = get_mu("same")$hi_bp,

  mean_opposite_bp = get_mu("opposite")$mean_bp,
  lo_opposite_bp = get_mu("opposite")$lo_bp,
  hi_opposite_bp = get_mu("opposite")$hi_bp,

  beta_log10_opposite_vs_same = beta_est,
  beta_lo = beta_lo,
  beta_hi = beta_hi,

  fold_opposite_over_same = fold_est,
  fold_lo = fold_lo,
  fold_hi = fold_hi,

  P_beta_gt_0 = p_beta_gt0,

  sd_phylo_log10 = sd_phylo,
  sigma_log10 = sigma_res
)

message("Saving results row to: ", results_name)
write.table(
  res,
  file = results_name,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Done!")
