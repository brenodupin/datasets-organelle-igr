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
    "Usage: Rscript create_brm.R <tree.nwk> <igs_summary.tsv> <summary.txt> <moldel.rds> <results_row.tsv>" # nolint
  )
  quit(status = 1)
}

tree_path <- args[1]
igs_path <- args[2]
summary_name <- args[3]
out_name <- args[4]
results_name <- args[5]
factor_levels <- c("same", "convergent", "divergent")

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
  message(sprintf("Subsampling %d -> %d rows", nrow(igr), max_rows))
  igr <- do.call(rbind, lapply(
    split(igr, igr$polarity_bin),
    function(grp) {
      n <- max(1, round(max_rows * nrow(grp) / nrow(igr)))
      grp[sample(nrow(grp), min(n, nrow(grp))), ]
    }
  ))
  igr$taxon_tree <- droplevels(igr$taxon_tree)
}

tree <- keep.tip(tree, levels(igr$taxon_tree))
tree <- multi2di(tree)
A <- vcv(tree, corr = FALSE) # nolint

n_cores <- parallel::detectCores(logical = FALSE)
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
coefs <- c("b_polarity_binconvergent", "b_polarity_bindivergent")
for (coef in coefs) {
  b <- ps[coef, "Estimate"]
  lo <- ps[coef, "Q2.5"]
  hi <- ps[coef, "Q97.5"]
  message(sprintf("%s: %.3f (%.3f - %.3f)", coef, b, lo, hi))
}

message("Saving summary to: ", summary_name)
sink(summary_name)
cat("polarity_bin levels (first = baseline):\n")
print(levels(igr$polarity_bin))
cat("baseline polarity_bin:\n")
print(levels(igr$polarity_bin)[1])
cat("\n--- brms summary ---\n")
print(summary(fit))
for (coef in coefs) {
  b <- ps[coef, "Estimate"]
  lo <- ps[coef, "Q2.5"]
  hi <- ps[coef, "Q97.5"]
  cat(sprintf("Delta PGLMM %-30s = %.3f (%.3f - %.3f)\n", coef, b, lo, hi))
  cat(sprintf(
    "Fold  PGLMM %-30s = %.3f (%.3f - %.3f)\n",
    coef,
    10^b,
    10^lo,
    10^hi
  ))
}
sink()

# ---------- Fitted means (population-level, back-transformed) ----------
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

# ---------- Raw medians ----------
igr$length_bp <- 10^igr$log10_length
med_same <- median(igr$length_bp[igr$polarity_bin == "same"], na.rm = TRUE)
med_conv <- median(
  igr$length_bp[igr$polarity_bin == "convergent"],
  na.rm = TRUE
)
med_div <- median(igr$length_bp[igr$polarity_bin == "divergent"], na.rm = TRUE)

# ---------- Fixed effects ----------
fx <- fixef(fit)
for (coef in c("polarity_binconvergent", "polarity_bindivergent")) {
  if (!(coef %in% rownames(fx))) {
    stop(
      "Expected coefficient '",
      coef,
      "' not found. Check polarity_bin levels/coding."
    )
  }
}

beta_conv_est <- fx["polarity_binconvergent", "Estimate"]
beta_conv_lo <- fx["polarity_binconvergent", "Q2.5"]
beta_conv_hi <- fx["polarity_binconvergent", "Q97.5"]

beta_div_est <- fx["polarity_bindivergent", "Estimate"]
beta_div_lo <- fx["polarity_bindivergent", "Q2.5"]
beta_div_hi <- fx["polarity_bindivergent", "Q97.5"]

# ---------- Posterior probabilities ----------
dr <- as_draws_df(fit)
p_conv_gt0 <- mean(dr$b_polarity_binconvergent > 0)
p_div_gt0 <- mean(dr$b_polarity_bindivergent > 0)

# ---------- Phylogenetic / residual SD ----------
sd_phylo <- as.numeric(VarCorr(fit)$taxon_tree$sd[1])
sigma_res <- as.numeric(VarCorr(fit)$residual__$sd[1])

# ---------- Results row ----------
get_mu <- function(pol) {
  mu_bp[mu_bp$polarity_bin == pol, ]
}

res <- data.frame(
  N_regions = nrow(igr),
  N_taxa = nlevels(igr$taxon_tree),

  median_same_bp = med_same,
  median_convergent_bp = med_conv,
  median_divergent_bp = med_div,

  mean_same_bp = get_mu("same")$mean_bp,
  lo_same_bp = get_mu("same")$lo_bp,
  hi_same_bp = get_mu("same")$hi_bp,

  mean_convergent_bp = get_mu("convergent")$mean_bp,
  lo_convergent_bp = get_mu("convergent")$lo_bp,
  hi_convergent_bp = get_mu("convergent")$hi_bp,

  mean_divergent_bp = get_mu("divergent")$mean_bp,
  lo_divergent_bp = get_mu("divergent")$lo_bp,
  hi_divergent_bp = get_mu("divergent")$hi_bp,

  beta_log10_convergent_vs_same = beta_conv_est,
  beta_convergent_lo = beta_conv_lo,
  beta_convergent_hi = beta_conv_hi,
  fold_convergent_over_same = 10^beta_conv_est,
  fold_convergent_lo = 10^beta_conv_lo,
  fold_convergent_hi = 10^beta_conv_hi,
  P_beta_convergent_gt_0 = p_conv_gt0,

  beta_log10_divergent_vs_same = beta_div_est,
  beta_divergent_lo = beta_div_lo,
  beta_divergent_hi = beta_div_hi,
  fold_divergent_over_same = 10^beta_div_est,
  fold_divergent_lo = 10^beta_div_lo,
  fold_divergent_hi = 10^beta_div_hi,
  P_beta_divergent_gt_0 = p_div_gt0,

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
