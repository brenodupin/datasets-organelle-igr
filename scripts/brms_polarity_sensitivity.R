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
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
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
    "Usage: Rscript brms_polarity.R <tree.nwk> <igs_summary.tsv> ",
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

# No subsampling: every IGR in the input file is used.
log_time(sprintf(
  "Using all %d IGRs from %d taxa (no subsampling)",
  nrow(igr),
  nlevels(igr$ncbi_taxid)
))

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

# ---------------------------------------------------------------------------
# Extra outputs for the revision (appended after the original columns so
# existing downstream code keeps working)
# ---------------------------------------------------------------------------
log_time("Computing extra outputs for the revision")

q3 <- function(x) {
  unname(c(median(x), quantile(x, 0.025), quantile(x, 0.975)))
}

# Posterior draws of the key parameters
dr_all <- posterior::as_draws_df(
  m_pol,
  variable = c(
    "b_Intercept", coefs, "sd_ncbi_taxid__Intercept", "sigma"
  )
)
b0 <- dr_all$b_Intercept
b_conv <- dr_all$b_polarity_3binconv
b_div <- dr_all$b_polarity_3bindiv
sd_tax <- dr_all$sd_ncbi_taxid__Intercept
sig <- dr_all$sigma

# ICC_taxon (R2.6, minor 6)
icc <- q3(sd_tax^2 / (sd_tax^2 + sig^2))

# Divergent vs convergent contrast (R1.3, R2.8)
d_dc <- b_div - b_conv
dc_beta <- q3(d_dc)
p_div_gt_conv <- mean(d_dc > 0)

# Model-implied absolute differences in bp for a typical taxon (R2.3)
same_bp_draws <- 10^b0
diff_conv_bp <- q3(10^(b0 + b_conv) - same_bp_draws)
diff_div_bp <- q3(10^(b0 + b_div) - same_bp_draws)
diff_div_conv_bp <- q3(10^(b0 + b_div) - 10^(b0 + b_conv))

# Raw per-class counts, medians and IQRs in bp (R2.3, R2.7)
raw_stats <- function(pol) {
  x <- igr$igr_len[igr$polarity_3bin == pol]
  c(
    n = length(x),
    median = median(x, na.rm = TRUE),
    q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
    q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
    n_len_le_1 = sum(x <= 1, na.rm = TRUE)
  )
}
rs_same <- raw_stats("same")
rs_conv <- raw_stats("conv")
rs_div <- raw_stats("div")

# Dataset composition (minor 6, R2.7)
genomes <- unique(igr[, c("AN", "genome_length")])

# Sampler diagnostics (minor 5; needed to report the new fits)
diag_vars <- c("b_Intercept", coefs, "sd_ncbi_taxid__Intercept", "sigma")
ess <- posterior::summarise_draws(
  posterior::as_draws_array(m_pol, variable = diag_vars),
  "ess_bulk", "ess_tail"
)
np <- brms::nuts_params(m_pol)
n_divergent <- sum(np$Value[np$Parameter == "divergent__"])
n_max_treedepth <- sum(np$Value[np$Parameter == "treedepth__"] >= 12)
max_rhat <- max(brms::rhat(m_pol), na.rm = TRUE)

extra <- data.frame(
  N_genomes = nrow(genomes),
  median_genome_length_bp = median(genomes$genome_length, na.rm = TRUE),
  median_igr_len_bp = median(igr$igr_len, na.rm = TRUE),

  N_same = rs_same[["n"]],
  N_convergent = rs_conv[["n"]],
  N_divergent = rs_div[["n"]],
  raw_q1_same_bp = rs_same[["q1"]],
  raw_q3_same_bp = rs_same[["q3"]],
  raw_q1_convergent_bp = rs_conv[["q1"]],
  raw_q3_convergent_bp = rs_conv[["q3"]],
  raw_q1_divergent_bp = rs_div[["q1"]],
  raw_q3_divergent_bp = rs_div[["q3"]],
  N_len_le_1bp_same = rs_same[["n_len_le_1"]],
  N_len_le_1bp_convergent = rs_conv[["n_len_le_1"]],
  N_len_le_1bp_divergent = rs_div[["n_len_le_1"]],

  diff_convergent_minus_same_bp = diff_conv_bp[1],
  diff_convergent_minus_same_lo = diff_conv_bp[2],
  diff_convergent_minus_same_hi = diff_conv_bp[3],
  diff_divergent_minus_same_bp = diff_div_bp[1],
  diff_divergent_minus_same_lo = diff_div_bp[2],
  diff_divergent_minus_same_hi = diff_div_bp[3],

  beta_log10_divergent_vs_convergent = dc_beta[1],
  beta_div_vs_conv_lo = dc_beta[2],
  beta_div_vs_conv_hi = dc_beta[3],
  fold_divergent_over_convergent = 10^dc_beta[1],
  fold_div_over_conv_lo = 10^dc_beta[2],
  fold_div_over_conv_hi = 10^dc_beta[3],
  P_divergent_gt_convergent = p_div_gt_conv,
  diff_divergent_minus_convergent_bp = diff_div_conv_bp[1],
  diff_divergent_minus_convergent_lo = diff_div_conv_bp[2],
  diff_divergent_minus_convergent_hi = diff_div_conv_bp[3],

  ICC_taxon = icc[1],
  ICC_taxon_lo = icc[2],
  ICC_taxon_hi = icc[3],

  max_rhat = max_rhat,
  min_ess_bulk = min(ess$ess_bulk),
  min_ess_tail = min(ess$ess_tail),
  n_divergent_transitions = n_divergent,
  n_max_treedepth_hits = n_max_treedepth
)
res <- cbind(res, extra)

# Posterior predictive check plot (minor 5); failure here never stops the run
ppc_path <- sub("\\.rds$", "_ppcheck.pdf", out_name)
tryCatch(
  {
    set.seed(seed)
    p <- brms::pp_check(m_pol, type = "dens_overlay", ndraws = 50)
    ggplot2::ggsave(ppc_path, p, width = 6, height = 4)
    message("Saved posterior predictive check to: ", ppc_path)
  },
  error = function(e) message("pp_check failed: ", conditionMessage(e))
)

# Append the new results and session info to the summary file
append_summary <- function(path) {
  con <- file(path, "a")
  sink(con)
  on.exit({
    sink()
    close(con)
  }, add = TRUE)
  cat("\n--- Divergent vs convergent ---\n")
  cat(sprintf(
    "Fold div/conv = %.3f (%.3f - %.3f); P(div > conv) = %.4f\n",
    10^dc_beta[1], 10^dc_beta[2], 10^dc_beta[3], p_div_gt_conv
  ))
  cat(sprintf(
    "ICC_taxon = %.3f (%.3f - %.3f)\n", icc[1], icc[2], icc[3]
  ))
  cat("\n--- Sampler diagnostics ---\n")
  cat(sprintf(
    "max Rhat = %.4f; min bulk ESS = %.0f; min tail ESS = %.0f\n",
    max_rhat, min(ess$ess_bulk), min(ess$ess_tail)
  ))
  cat(sprintf(
    "divergent transitions = %d; max-treedepth hits = %d\n",
    as.integer(n_divergent), as.integer(n_max_treedepth)
  ))
  cat("\n--- Session info ---\n")
  print(sessionInfo())
}
append_summary(results_txt)

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