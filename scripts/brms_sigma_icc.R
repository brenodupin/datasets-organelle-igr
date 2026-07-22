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
log_time("Script started (sigma + ICC companion)")

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
flank_col <- "log10_max_flanking_gene_len"

# get path from out_name
cache_flank <- sub("\\.rds$", "_cache_flank.rds", out_name)
cache_type <- sub("\\.rds$", "_cache_type.rds", out_name)
cache_both <- sub("\\.rds$", "_cache_both.rds", out_name)
cache_sigma <- sub("\\.rds$", "_cache_sigma.rds", out_name)

for (cf in c(cache_flank, cache_type, cache_both)) {
  if (!file.exists(cf)) {
    stop(
      "Expected cache not found: ",
      cf,
      "\n",
      "Run brms_type_length.R for this group first before running this script."
    )
  }
}

log_time(paste("Loading tree from:       ", tree_path))
log_time(paste("Loading IGR summary from:", igs_path))

tree <- read.tree(tree_path)
tree <- multi2di(tree)

igr <- read.csv(igs_path, sep = "\t", stringsAsFactors = FALSE)
igr$ncbi_taxid <- factor(igr$ncbi_taxid)
igr$polarity_3bin <- factor(igr$polarity_3bin, levels = factor_levels)

n_missing <- sum(is.na(igr$flanking_types))
if (n_missing > 0) {
  message(sprintf(
    "Warning: %d rows have a missing flanking_types value (dropping them).",
    n_missing
  ))
  igr <- igr[!is.na(igr$flanking_types), ]
}

ft_observed <- sort(unique(igr$flanking_types))
ft_levels <- c("gene-gene", setdiff(ft_observed, "gene-gene"))
igr$flanking_types <- factor(igr$flanking_types, levels = ft_levels)

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
  igr$flanking_types <- droplevels(igr$flanking_types)
}

mu_flank <- mean(igr[[flank_col]])
igr$flank <- igr[[flank_col]] - mu_flank

tree <- keep.tip(tree, levels(igr$ncbi_taxid))
tree <- multi2di(tree)
A <- log_duration(
  "Building phylogenetic covariance matrix",
  vcv(tree, corr = TRUE)
) # nolint

message("N tips in pruned tree: ", length(tree$tip.label))
message(
  "Size of A (covariance matrix): ",
  format(object.size(A), units = "auto")
)

total_cores <- parallel::detectCores(logical = FALSE)
n_chains <- 4
iter_n <- 6000
warmup_n <- 2000
options(mc.cores = n_chains)
message("Cores available: ", total_cores)
message("Chains requested: ", n_chains)
message("Iterations per chain: ", iter_n)
message("Warmup iterations per chain: ", warmup_n)

form_sigma <- bf(
  log10_igr_len ~ polarity_3bin +
    flank +
    flanking_types +
    (1 | gr(ncbi_taxid, cov = A)),
  sigma ~ polarity_3bin
)

log_time("Fitting distributional model: sigma ~ polarity_3bin")
m_sigma <- log_duration(
  "sigma model",
  brm(
    form_sigma,
    data = igr,
    data2 = list(A = A),
    family = gaussian(),
    chains = n_chains,
    cores = n_chains,
    iter = iter_n,
    warmup = warmup_n,
    file = cache_sigma,
    file_refit = "on_change",
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
)

log_time("Reloading cached fits: m_len, m_type, m_both")
m_len <- readRDS(cache_flank)
m_type <- readRDS(cache_type)
m_both <- readRDS(cache_both)

icc_of <- function(fit) {
  dr <- posterior::as_draws_df(fit)
  vp <- dr$sd_ncbi_taxid__Intercept^2
  ve <- dr$sigma^2
  icc <- vp / (vp + ve)
  c(
    median = median(icc),
    lo = as.numeric(quantile(icc, 0.025)),
    hi = as.numeric(quantile(icc, 0.975))
  )
}

icc_print <- function(fit, label) {
  icc <- icc_of(fit)
  cat(sprintf(
    "  [%-6s] ICC tax = %.4f (%.4f, %.4f)\n",
    label,
    icc["median"],
    icc["lo"],
    icc["hi"]
  ))
}

sigma_print <- function(m) {
  fx <- brms::fixef(m)
  rns <- grep("^sigma_polarity_3bin", rownames(fx), value = TRUE)
  if (length(rns) == 0) {
    cat("  (no sigma_polarity_3bin terms found)\n")
    return(invisible(NULL))
  }
  e <- fx[rns, c("Estimate", "Q2.5", "Q97.5"), drop = FALSE]
  rownames(e) <- sub("^sigma_polarity_3bin", "", rownames(e))
  cat("  SD ratio vs same-polarity baseline (exp of log-sigma coefficient):\n")
  print(round(exp(e), 3))
  cat(sprintf(
    "\n  Baseline residual SD (same-polarity): %.4f\n",
    exp(fx["sigma_Intercept", "Estimate"])
  ))
}

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

  cat("--- Proportion of residual variance among taxa ---\n")
  cat("--- (Intraclass correlation coefficient) ---\n")
  cat("ICC = var_phylo / (var_phylo + var_resid), from posterior draws\n")
  icc_print(m_len, "len")
  icc_print(m_type, "type")
  icc_print(m_both, "both")

  cat("\n--- Distributional model: sigma ~ polarity_3bin ---\n")
  sigma_print(m_sigma)

  cat("\n--- brms summary (distributional model) ---\n")
  print(summary(m_sigma))
}
write_summary(results_txt)

# Results TSV
dr_sigma <- posterior::as_draws_df(m_sigma)
fx_sigma <- brms::fixef(m_sigma)

sigma_ratio <- function(term) {
  e <- unname(fx_sigma[paste0("sigma_", term), c("Estimate", "Q2.5", "Q97.5")])
  c(ratio = exp(e[1]), lo = exp(e[2]), hi = exp(e[3]))
}
sr_conv <- sigma_ratio("polarity_3binconv")
sr_div <- sigma_ratio("polarity_3bindiv")

p_sigma_conv_gt0 <- mean(dr_sigma$b_sigma_polarity_3binconv > 0)
p_sigma_div_gt0 <- mean(dr_sigma$b_sigma_polarity_3bindiv > 0)

icc_len <- icc_of(m_len)
icc_type <- icc_of(m_type)
icc_both <- icc_of(m_both)

res <- data.frame(
  N_regions = nrow(igr),
  N_taxa = nlevels(igr$ncbi_taxid),

  icc_len_median = icc_len["median"],
  icc_len_lo = icc_len["lo"],
  icc_len_hi = icc_len["hi"],
  icc_type_median = icc_type["median"],
  icc_type_lo = icc_type["lo"],
  icc_type_hi = icc_type["hi"],
  icc_both_median = icc_both["median"],
  icc_both_lo = icc_both["lo"],
  icc_both_hi = icc_both["hi"],

  sigma_ratio_conv = sr_conv["ratio"],
  sigma_ratio_conv_lo = sr_conv["lo"],
  sigma_ratio_conv_hi = sr_conv["hi"],
  P_sigma_conv_gt_0 = p_sigma_conv_gt0,

  sigma_ratio_div = sr_div["ratio"],
  sigma_ratio_div_lo = sr_div["lo"],
  sigma_ratio_div_hi = sr_div["hi"],
  P_sigma_div_gt_0 = p_sigma_div_gt0,

  sigma_baseline_sd = exp(fx_sigma["sigma_Intercept", "Estimate"])
)

log_time(paste("Saving results to:", results_tsv))
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
