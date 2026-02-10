#!/usr/bin/env Rscript

required_packages <- c("brms", "ape")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  message("Error: Missing required R packages: ", paste(missing_packages, collapse = ", "))
  message("Install with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
  quit(status = 1)
}

library(brms)
library(ape)
library(posterior)

set.seed(42)  # Same seed for reproducibility

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  message("Error: Exactly 4 arguments required")
  message("Usage: Rscript create_brm.R <tree.nwk> <igs_summary.tsv> <summary.txt> <moldel.rds> <results_row.tsv>")
  quit(status = 1)
}

tree_path <- args[1]
igs_path <- args[2]
summary_name <- args[3]
out_name <- args[4]
results_name <- args[5]

if (!file.exists(tree_path)) stop("File not found: ", tree_path)
if (!file.exists(igs_path)) stop("File not found: ", igs_path)

message("Loading tree from: ", tree_path)
message("Loading igs data from: ", igs_path)

tree <- read.tree(tree_path)
tree <- multi2di(tree)
A <- vcv(tree, corr = FALSE)

igr <- read.csv(igs_path, sep = "\t")
igr$taxon_tree <- factor(igr$taxon_tree)
igr$polarity_bin <- factor(igr$polarity_bin, levels = c("same", "opposite"))

n_cores <- parallel::detectCores(logical = FALSE)
message("Running brm with ", n_cores, " cores")

fit <- brm(
  log10_length ~ polarity_bin + (1 | gr(taxon_tree, cov = A)),
  data = igr,
  data2 = list(A = A),
  family = gaussian(),
  cores = n_cores,
  chains = 10
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

# ---------- Results-row (for your PNAS Table 1) ----------
# Population-level predicted means per polarity (re_formula = NA), then back-transform to bp
nd <- data.frame(polarity_bin = factor(c("same", "opposite"), levels = c("same", "opposite")))
mu <- fitted(fit, newdata = nd, re_formula = NA, summary = TRUE)[, c("Estimate", "Q2.5", "Q97.5")]

mu_bp <- data.frame(
  polarity_bin = nd$polarity_bin,
  mean_bp = 10^mu[, "Estimate"],
  lo_bp   = 10^mu[, "Q2.5"],
  hi_bp   = 10^mu[, "Q97.5"]
)

# Raw medians from observed data (bp scale)
igr$length_bp <- 10^igr$log10_length
med_same <- median(igr$length_bp[igr$polarity_bin == "same"], na.rm = TRUE)
med_opp  <- median(igr$length_bp[igr$polarity_bin == "opposite"], na.rm = TRUE)

# Fixed effect (now should be polarity_binopposite, because baseline = same)
fx <- fixef(fit)
if (!("polarity_binopposite" %in% rownames(fx))) {
  stop("Expected coefficient 'polarity_binopposite' not found. Check polarity_bin levels/coding.")
}

beta_est <- fx["polarity_binopposite", "Estimate"]
beta_lo  <- fx["polarity_binopposite", "Q2.5"]
beta_hi  <- fx["polarity_binopposite", "Q97.5"]

# Fold-change (opposite / same) on bp scale
fold_est <- 10^beta_est
fold_lo  <- 10^beta_lo
fold_hi  <- 10^beta_hi

# Posterior probability that opposite > same (i.e., beta > 0)
dr <- as_draws_df(fit)
p_beta_gt0 <- mean(dr$b_polarity_binopposite > 0)

# Pull phylogenetic SD and residual sigma (optional, but useful as context)
sd_phylo <- as.numeric(VarCorr(fit)$taxon_tree$sd[1])
sigma_res <- as.numeric(VarCorr(fit)$residual__$sd[1])

res <- data.frame(
  N_regions = nrow(igr),
  N_taxa    = nlevels(igr$taxon_tree),

  median_same_bp     = med_same,
  median_opposite_bp = med_opp,
  delta_median_bp    = med_opp - med_same,

  mean_same_bp = mu_bp$mean_bp[mu_bp$polarity_bin == "same"],
  lo_same_bp   = mu_bp$lo_bp[mu_bp$polarity_bin == "same"],
  hi_same_bp   = mu_bp$hi_bp[mu_bp$polarity_bin == "same"],

  mean_opposite_bp = mu_bp$mean_bp[mu_bp$polarity_bin == "opposite"],
  lo_opposite_bp   = mu_bp$lo_bp[mu_bp$polarity_bin == "opposite"],
  hi_opposite_bp   = mu_bp$hi_bp[mu_bp$polarity_bin == "opposite"],

  beta_log10_opposite_vs_same = beta_est,
  beta_lo = beta_lo,
  beta_hi = beta_hi,

  fold_opposite_over_same = fold_est,
  fold_lo = fold_lo,
  fold_hi = fold_hi,

  P_beta_gt_0 = p_beta_gt0,

  sd_phylo_log10 = sd_phylo,
  sigma_log10    = sigma_res
)

message("Saving results row to: ", results_name)
write.table(res, file = results_name, sep = "\t", quote = FALSE, row.names = FALSE)

message("Done!")
