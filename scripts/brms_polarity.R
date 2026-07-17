#!/usr/bin/env Rscript

options(warn = 1)
required_packages <- c("brms", "ape", "posterior")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

for (pkg in required_packages) {
  library(pkg, character.only = TRUE, quietly = TRUE)
}

log_time <- function(msg) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

log_duration <- function(label, expr) {
  t0 <- Sys.time()
  result <- force(expr)
  message(sprintf(
    "  -> %s finished in %.1f min",
    label,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
  result
}

script_start <- Sys.time()
log_time("Script started (polarity model)")

seed <- 42
set.seed(seed)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    "Exactly 5 arguments required\n",
    "Usage: Rscript brms_type_length.R <tree.nwk> <igs_summary.tsv> ",
    "<summary.txt> <models.rds> <results.tsv>"
  )
}

tree_path <- args[1]
igs_path <- args[2]
results_txt <- args[3]
out_name <- args[4]
results_tsv <- args[5]

factor_levels <- c("same", "conv", "div")

log_time(paste("Loading tree from: ", tree_path))
log_time(paste("Loading IGR from : ", igs_path))

tree <- read.tree(tree_path)
tree <- multi2di(tree)

igr <- read.csv(igs_path, sep = "\t", stringsAsFactors = FALSE)
igr$ncbi_taxid <- factor(igr$ncbi_taxid)
igr$polarity_3bin <- factor(igr$polarity_3bin, levels = factor_levels)

# subsample
max_rows <- 150000
if (nrow(igr) > max_rows) {
  total_n <- nrow(igr)
  log_time(sprintf(
    "Subsampling %d -> %d rows (stratified by polarity_3bin)",
    total_n,
    max_rows
  ))
  igr <- do.call(
    rbind,
    lapply(
      split(igr, igr$polarity_3bin),
      function(grp) {
        n <- max(1, round(max_rows * nrow(grp) / total_n))
        grp[sample(nrow(grp), min(n, nrow(grp))), ]
      }
    )
  )
  igr$ncbi_taxid <- droplevels(igr$ncbi_taxid)
}

tree <- keep.tip(tree, levels(igr$ncbi_taxid))
tree <- multi2di(tree)
A <- log_duration(
  "Building phylogenetic covariance matrix",
  vcv(tree, corr = TRUE)
) # nolint

message("N tips in tree: ", length(tree$tip.label))
message(
  "Size of A (covariance matrix): ",
  format(object.size(A), units = "auto")
)

total_cores <- parallel::detectCores(logical = FALSE) # 32
n_chains <- 4
iter_n <- 6000
warmup_n <- 2000
options(mc.cores = n_chains)
message("Cores available: ", total_cores)
message("Chains requested: ", n_chains)
message("Iterations per chain: ", iter_n)
message("Warmup iterations per chain: ", warmup_n)

# Draw subsampling for bayes_R2
n_r2_draws <- 2000L

cache_pol <- sub("\\.rds$", "_cache_polarity.rds", out_name)

log_time("Fitting model 1/1: polarity")
m_pol <- log_duration(
  "model 1/1 (polarity)",
  brm(
    log10_igr_len ~ polarity_3bin + (1 | gr(ncbi_taxid, cov = A)),
    data = igr,
    data2 = list(A = A),
    family = gaussian(),
    chains = n_chains,
    cores = n_chains,
    iter = iter_n,
    warmup = warmup_n,
    file = cache_pol,
    file_refit = "on_change",
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
)

ps <- posterior_summary(m_pol)
coefs <- c("b_polarity_3binconv", "b_polarity_3bindiv")
for (coef in coefs) {
  b <- ps[coef, "Estimate"]
  lo <- ps[coef, "Q2.5"]
  hi <- ps[coef, "Q97.5"]
  message(sprintf("%s: %.3f (%.3f - %.3f)", coef, b, lo, hi))
}

log_time("Computing Bayesian R^2 (marginal + conditional)")
set.seed(seed)
r2_marg <- log_duration(
  "polarity marginal R^2",
  brms::bayes_R2(m_pol, re_formula = NA, ndraws = n_r2_draws, summary = TRUE)
)
set.seed(seed)
r2_cond <- log_duration(
  "polarity conditional R^2",
  brms::bayes_R2(m_pol, re_formula = NULL, ndraws = n_r2_draws, summary = TRUE)
)

r2_marg_est <- unname(r2_marg[1, "Estimate"])
r2_marg_lo <- unname(r2_marg[1, "Q2.5"])
r2_marg_hi <- unname(r2_marg[1, "Q97.5"])
r2_cond_est <- unname(r2_cond[1, "Estimate"])
r2_cond_lo <- unname(r2_cond[1, "Q2.5"])
r2_cond_hi <- unname(r2_cond[1, "Q97.5"])
r2_phylo_contribution <- r2_cond_est - r2_marg_est

message(sprintf(
  "Marginal R^2    = %.4f (%.4f, %.4f)",
  r2_marg_est,
  r2_marg_lo,
  r2_marg_hi
))
message(sprintf(
  "Conditional R^2 = %.4f (%.4f, %.4f)",
  r2_cond_est,
  r2_cond_lo,
  r2_cond_hi
))
message(sprintf(
  "Phylo contribution (conditional - marginal) = %.4f",
  r2_phylo_contribution
))

message("Saving summary to: ", results_txt)
write_summary <- function(path) {
  con <- file(path, "w")
  sink(con)
  on.exit(
    {
      sink()
      close(con)
    },
    add = TRUE
  )

  cat("polarity_3bin levels (first = baseline):\n")
  print(levels(igr$polarity_3bin))
  cat("baseline polarity_3bin:\n")
  print(levels(igr$polarity_3bin)[1])
  cat("\n--- brms summary ---\n")
  print(summary(m_pol))
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

  cat("\n--- Bayesian R^2 (Gelman, Goodrich, Gabry & Vehtari 2019) ---\n")
  cat(sprintf(
    "\n  marginal R^2    = %.4f (%.4f, %.4f)\n",
    r2_marg_est,
    r2_marg_lo,
    r2_marg_hi
  ))
  cat(sprintf(
    "  conditional R^2 = %.4f (%.4f, %.4f)\n",
    r2_cond_est,
    r2_cond_lo,
    r2_cond_hi
  ))
  cat(sprintf(
    "  phylo contribution (conditional - marginal) = %.4f\n",
    r2_phylo_contribution
  ))
}
write_summary(results_txt)

# Results TSV
nd <- data.frame(polarity_3bin = factor(factor_levels, levels = factor_levels))
mu <- fitted(m_pol, newdata = nd, re_formula = NA, summary = TRUE)[, c(
  "Estimate",
  "Q2.5",
  "Q97.5"
)]

mu_bp <- data.frame(
  polarity_3bin = nd$polarity_3bin,
  mean_bp = 10^mu[, "Estimate"],
  lo_bp = 10^mu[, "Q2.5"],
  hi_bp = 10^mu[, "Q97.5"]
)

# Raw medians from observed data (bp scale)
igr$length_bp <- 10^igr$log10_igr_len
med_same <- median(igr$length_bp[igr$polarity_3bin == "same"], na.rm = TRUE)
med_conv <- median(
  igr$length_bp[igr$polarity_3bin == "conv"],
  na.rm = TRUE
)
med_div <- median(igr$length_bp[igr$polarity_3bin == "div"], na.rm = TRUE)

# Fixed effect
fx <- brms::fixef(m_pol)
for (coef in c("polarity_3binconv", "polarity_3bindiv")) {
  if (!(coef %in% rownames(fx))) {
    stop(
      "Expected coefficient '",
      coef,
      "' not found. Check polarity_3bin levels/coding."
    )
  }
}

beta_conv_est <- fx["polarity_3binconv", "Estimate"]
beta_conv_lo <- fx["polarity_3binconv", "Q2.5"]
beta_conv_hi <- fx["polarity_3binconv", "Q97.5"]

beta_div_est <- fx["polarity_3bindiv", "Estimate"]
beta_div_lo <- fx["polarity_3bindiv", "Q2.5"]
beta_div_hi <- fx["polarity_3bindiv", "Q97.5"]

dr <- posterior::as_draws_df(m_pol, variable = coefs)
p_conv_gt0 <- mean(dr$b_polarity_3binconv > 0)
p_div_gt0 <- mean(dr$b_polarity_3bindiv > 0)

sd_phylo <- unname(brms::VarCorr(m_pol)$ncbi_taxid$sd["Intercept", "Estimate"])
sigma_res <- unname(brms::VarCorr(m_pol)$residual__$sd[1, "Estimate"])

get_mu <- function(pol) {
  mu_bp[mu_bp$polarity_3bin == pol, ]
}

res <- data.frame(
  N_regions = nrow(igr),
  N_taxa = nlevels(igr$ncbi_taxid),

  median_same_bp = med_same,
  median_convergent_bp = med_conv,
  median_divergent_bp = med_div,

  mean_same_bp = get_mu("same")$mean_bp,
  lo_same_bp = get_mu("same")$lo_bp,
  hi_same_bp = get_mu("same")$hi_bp,

  mean_convergent_bp = get_mu("conv")$mean_bp,
  lo_convergent_bp = get_mu("conv")$lo_bp,
  hi_convergent_bp = get_mu("conv")$hi_bp,

  mean_divergent_bp = get_mu("div")$mean_bp,
  lo_divergent_bp = get_mu("div")$lo_bp,
  hi_divergent_bp = get_mu("div")$hi_bp,

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
  sigma_log10 = sigma_res,

  marginal_R2 = r2_marg_est,
  marginal_R2_lo = r2_marg_lo,
  marginal_R2_hi = r2_marg_hi,
  conditional_R2 = r2_cond_est,
  conditional_R2_lo = r2_cond_lo,
  conditional_R2_hi = r2_cond_hi,
  phylo_contribution = r2_phylo_contribution
)

message("Saving results row to: ", results_tsv)
write.table(
  res,
  file = results_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(sprintf(
  "Total runtime: %.1f min",
  as.numeric(difftime(Sys.time(), script_start, units = "mins"))
))
message("Done!")
