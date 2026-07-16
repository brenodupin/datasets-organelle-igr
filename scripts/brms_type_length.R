#!/usr/bin/env Rscript

options(warn = 1)
required_packages <- c("brms", "ape", "posterior", "loo", "callr")
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
log_time("Script started (brms type, length and both models)")

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

log_time(paste("Loading tree from:       ", tree_path))
log_time(paste("Loading IGR from :", igs_path))

tree <- read.tree(tree_path)
tree <- multi2di(tree)

igr <- read.csv(igs_path, sep = "\t", stringsAsFactors = FALSE)
igr$ncbi_taxid <- factor(igr$ncbi_taxid)
igr$polarity_3bin <- factor(igr$polarity_3bin, levels = factor_levels)

ft_observed <- sort(unique(igr$flanking_types))
ft_levels <- c("gene-gene", setdiff(ft_observed, "gene-gene"))
igr$flanking_types <- factor(igr$flanking_types, levels = ft_levels)

message("\nflanking_types distribution (pre-downsampling):")
print(table(igr$flanking_types))

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
  igr$flanking_types <- droplevels(igr$flanking_types)
}

# Center the chosen length summary on its post-downsampling mean
mu_flank <- mean(igr[[flank_col]])
igr$flank <- igr[[flank_col]] - mu_flank

message(sprintf("\nFlank column: %s", flank_col))
message(sprintf("  mu_flank = %.4f", mu_flank))

# Phylogenetic covariance matrix
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

total_cores <- parallel::detectCores(logical = FALSE) # 32
n_chains <- 4
iter_n <- 6000
warmup_n <- 2000
options(mc.cores = n_chains)
message("Cores available: ", total_cores)
message("Chains requested: ", n_chains)
message("Iterations per chain: ", iter_n)
message("Warmup iterations per chain: ", warmup_n)

k_folds <- 3

cache_flank <- sub("\\.rds$", "_cache_flank.rds", out_name)
cache_type <- sub("\\.rds$", "_cache_type.rds", out_name)
cache_both <- sub("\\.rds$", "_cache_both.rds", out_name)

##### Model 1: continuous gene-length predictor
form_len <- bf(
  log10_igr_len ~ polarity_3bin + flank + (1 | gr(ncbi_taxid, cov = A))
)

log_time(paste("Fitting model 1/3:", flank_col))
m_len <- log_duration(
  "model 1/3 (flank)",
  brm(
    form_len,
    data = igr,
    data2 = list(A = A),
    family = gaussian(),
    chains = n_chains,
    cores = n_chains,
    iter = iter_n,
    warmup = warmup_n,
    file = cache_flank,
    file_refit = "on_change",
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
)

##### Model 2: categorical gene-type predictor
form_type <- bf(
  log10_igr_len ~ polarity_3bin +
    flanking_types +
    (1 | gr(ncbi_taxid, cov = A))
)

log_time("Fitting model 2/3: flanking_types")
m_type <- log_duration(
  "model 2/3 (flanking_types)",
  brm(
    form_type,
    data = igr,
    data2 = list(A = A),
    family = gaussian(),
    chains = n_chains,
    cores = n_chains,
    iter = iter_n,
    warmup = warmup_n,
    file = cache_type,
    file_refit = "on_change",
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
)

##### Model 3: both predictors together
form_both <- bf(
  log10_igr_len ~ polarity_3bin +
    flank +
    flanking_types +
    (1 | gr(ncbi_taxid, cov = A))
)

log_time("Fitting model 3/3: both predictors")
m_both <- log_duration(
  "model 3/3 (both)",
  brm(
    form_both,
    data = igr,
    data2 = list(A = A),
    family = gaussian(),
    chains = n_chains,
    cores = n_chains,
    iter = iter_n,
    warmup = warmup_n,
    file = cache_both,
    file_refit = "on_change",
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
)

# We need to run kfold isolated, due to memory issues that
# 3 runs in a single R session can cause.
kfold_worker <- function(
  cache_file,
  K,
  seed,
  chains = 1,
  iter = iter_n,
  warmup = warmup_n
) {
  suppressPackageStartupMessages(library(brms))
  if (exists("mem.maxVSize")) {
    try(mem.maxVSize(vsize = Inf), silent = TRUE)
  }

  fit <- readRDS(cache_file)
  taxon <- fit$data$ncbi_taxid
  set.seed(seed)
  folds <- integer(length(taxon))
  for (idx in split(seq_along(taxon), taxon)) {
    folds[idx] <- (sample.int(length(idx)) %% K) + 1L
  }

  fit <- add_criterion(
    fit,
    "kfold",
    folds = folds,
    seed = seed,
    chains = chains,
    iter = iter,
    warmup = warmup
  )
  fit$criteria$kfold
}

run_kfold_isolated <- function(
  cache_file,
  K = k_folds,
  seed = 42,
  chains = 1,
  iter = iter_n,
  warmup = warmup_n
) {
  if (!file.exists(cache_file)) {
    stop("Cache not found: ", cache_file)
  }
  callr::r(
    kfold_worker,
    args = list(
      cache_file = cache_file,
      K = K,
      seed = seed,
      chains = chains,
      iter = iter,
      warmup = warmup
    ),
    show = TRUE
  )
}

# Free the models from the main process
rm(m_len, m_type, m_both)
gc()

log_time(paste0(
  "Computing K-fold (K=",
  k_folds,
  "): one isolated process per model"
))
kf_len <- log_duration(
  "kfold (flank)",
  run_kfold_isolated(cache_flank, K = k_folds, seed = seed)
)
kf_type <- log_duration(
  "kfold (flanking_types)",
  run_kfold_isolated(cache_type, K = k_folds, seed = seed)
)
kf_both <- log_duration(
  "kfold (both)",
  run_kfold_isolated(cache_both, K = k_folds, seed = seed)
)

# Reload the models for the results table + summary, attaching the CV criteria.
m_len <- readRDS(cache_flank)
m_len$criteria$kfold <- kf_len

m_type <- readRDS(cache_type)
m_type$criteria$kfold <- kf_type

m_both <- readRDS(cache_both)
m_both$criteria$kfold <- kf_both

loo_tab <- loo_compare(m_len, m_type, m_both, criterion = "kfold")

# Centralised kfold accessor — avoids repeating the string in every build_row
get_elpd <- function(fit) {
  est <- fit$criteria$kfold$estimates
  list(
    elpd = as.numeric(est["elpd_kfold", "Estimate"]),
    se = as.numeric(est["elpd_kfold", "SE"])
  )
}

fc_print <- function(m, label) {
  e <- brms::fixef(m)[
    c("polarity_3binconv", "polarity_3bindiv"),
    c("Estimate", "Q2.5", "Q97.5")
  ]
  cat(sprintf("\n  [%s]\n", label))
  print(round(10^e, 3))
}

flank_print <- function(m, label) {
  e <- brms::fixef(m)["flank", c("Estimate", "Q2.5", "Q97.5")]
  cat(sprintf(
    " [%-7s] beta = %.3f (%.3f, %.3f) fold per 10x flank = %.3f (%.3f, %.3f)\n",
    label,
    e[1],
    e[2],
    e[3],
    10^e[1],
    10^e[2],
    10^e[3]
  ))
}

type_contrasts_print <- function(m) {
  fx <- brms::fixef(m)
  rns <- grep("^flanking_types", rownames(fx), value = TRUE)
  if (length(rns) == 0) {
    cat("  (no flanking_types terms found)\n")
    return(invisible(NULL))
  }
  e <- fx[rns, c("Estimate", "Q2.5", "Q97.5"), drop = FALSE]
  rownames(e) <- sub("^flanking_types", "", rownames(e))
  cat("  Fold-change vs gene-gene (reference level):\n")
  print(round(10^e, 3))
}

n_r2_draws <- 2000L
calculate_r2 <- function(fit) {
  set.seed(seed)
  d_marg <- log_duration(
    "marginal R^2",
    brms::bayes_R2(fit, re_formula = NA, ndraws = n_r2_draws, summary = FALSE)
  )
  set.seed(seed)
  d_cond <- log_duration(
    "conditional R^2",
    brms::bayes_R2(fit, re_formula = NULL, ndraws = n_r2_draws, summary = FALSE)
  )

  d_marg <- as.numeric(d_marg[, 1])
  d_cond <- as.numeric(d_cond[, 1])
  d_diff <- d_cond - d_marg

  q <- function(x, p) unname(quantile(x, p))

  list(
    r2_marginal = mean(d_marg),
    r2_marginal_lo = q(d_marg, 0.025),
    r2_marginal_hi = q(d_marg, 0.975),
    r2_conditional = mean(d_cond),
    r2_conditional_lo = q(d_cond, 0.025),
    r2_conditional_hi = q(d_cond, 0.975),
    r2_phylo_contrib = mean(d_diff),
    r2_phylo_contrib_lo = q(d_diff, 0.025),
    r2_phylo_contrib_hi = q(d_diff, 0.975)
  )
}

r2_print <- function(r2) {
  cat(sprintf(
    "  Marginal R^2    = %.4f (%.4f, %.4f)\n",
    r2$r2_marginal,
    r2$r2_marginal_lo,
    r2$r2_marginal_hi
  ))
  cat(sprintf(
    "  Conditional R^2 = %.4f (%.4f, %.4f)\n",
    r2$r2_conditional,
    r2$r2_conditional_lo,
    r2$r2_conditional_hi
  ))
  cat(sprintf(
    "  Phylo contribution (conditional - marginal) = %.4f (%.4f, %.4f)\n",
    r2$r2_phylo_contrib,
    r2$r2_phylo_contrib_lo,
    r2$r2_phylo_contrib_hi
  ))
}

log_time("Computing Bayesian R^2 for all three models")
r2_len <- calculate_r2(m_len)
gc()
r2_type <- calculate_r2(m_type)
gc()
r2_both <- calculate_r2(m_both)
gc()

log_time(paste("Saving summary to:", results_txt))
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

  cat("\nFlank column: ", flank_col, "\n")
  cat(sprintf("  mu_flank = %.4f\n", mu_flank))

  cat("\nflanking_types distribution (post-downsampling):\n")
  print(table(igr$flanking_types))
  cat("  Reference level: gene-gene\n")

  cat(sprintf(
    "\n--- K-fold (K=%d) comparison (top = best by elpd) ---\n",
    k_folds
  ))
  print(loo_tab)
  cat("\n--- Polarity contrasts (all three models) ---\n")
  fc_print(m_len, flank_col)
  fc_print(m_type, "flanking_types")
  fc_print(m_both, "both")

  cat("\n--- Flank coefficient (gene-length and both models) ---\n")
  flank_print(m_len, "len")
  flank_print(m_both, "both")

  cat("\n--- Gene-type contrasts vs gene-gene ---\n")
  cat("\n  [flanking_types model]\n")
  type_contrasts_print(m_type)
  cat("\n  [both model — length-adjusted contrasts]\n")
  type_contrasts_print(m_both)

  cat("\n--- brms summary (gene-length model) ---\n")
  print(summary(m_len))
  cat("\n### Bayesian R^2 (gene-length model) ###\n")
  r2_print(r2_len)

  cat("\n--- brms summary (flanking_types model) ---\n")
  print(summary(m_type))
  cat("\n### Bayesian R^2 (flanking_types model) ###\n")
  r2_print(r2_type)

  cat("\n--- brms summary (both predictors model) ---\n")
  print(summary(m_both))
  cat("\n### Bayesian R^2 (both predictors model) ###\n")
  r2_print(r2_both)
}
write_summary(results_txt)

# Results TSV
build_row_len <- function(fit, flank_summary_name, mu_flank_val) {
  nd <- data.frame(
    polarity_3bin = factor(factor_levels, levels = factor_levels),
    flank = 0 # predict at centering point
  )
  mu <- fitted(fit, newdata = nd, re_formula = NA, summary = TRUE)[,
    c("Estimate", "Q2.5", "Q97.5")
  ]

  get_mu <- function(pol) {
    i <- which(nd$polarity_3bin == pol)
    list(
      mean = 10^mu[i, "Estimate"],
      lo = 10^mu[i, "Q2.5"],
      hi = 10^mu[i, "Q97.5"]
    )
  }

  fx <- brms::fixef(fit)
  beta_conv <- fx["polarity_3binconv", ]
  beta_div <- fx["polarity_3bindiv", ]
  beta_flank <- fx["flank", ]

  dr <- posterior::as_draws_df(fit)
  p_conv_gt0 <- mean(dr$b_polarity_3binconv > 0)
  p_div_gt0 <- mean(dr$b_polarity_3bindiv > 0)
  p_flank_gt0 <- mean(dr$b_flank > 0)

  sd_phylo <- as.numeric(brms::VarCorr(fit)$ncbi_taxid$sd[1])
  sigma_res <- as.numeric(brms::VarCorr(fit)$residual__$sd[1])
  ev <- get_elpd(fit)

  data.frame(
    predictor_type = "gene_length",
    flank_summary = flank_summary_name,
    mu_flank = mu_flank_val,
    N_regions = nrow(igr),
    N_taxa = nlevels(igr$ncbi_taxid),

    mean_same_bp = get_mu("same")$mean,
    lo_same_bp = get_mu("same")$lo,
    hi_same_bp = get_mu("same")$hi,
    mean_convergent_bp = get_mu("conv")$mean,
    lo_convergent_bp = get_mu("conv")$lo,
    hi_convergent_bp = get_mu("conv")$hi,
    mean_divergent_bp = get_mu("div")$mean,
    lo_divergent_bp = get_mu("div")$lo,
    hi_divergent_bp = get_mu("div")$hi,

    beta_log10_convergent_vs_same = beta_conv["Estimate"],
    beta_convergent_lo = beta_conv["Q2.5"],
    beta_convergent_hi = beta_conv["Q97.5"],
    fold_convergent_over_same = 10^beta_conv["Estimate"],
    fold_convergent_lo = 10^beta_conv["Q2.5"],
    fold_convergent_hi = 10^beta_conv["Q97.5"],
    P_beta_convergent_gt_0 = p_conv_gt0,

    beta_log10_divergent_vs_same = beta_div["Estimate"],
    beta_divergent_lo = beta_div["Q2.5"],
    beta_divergent_hi = beta_div["Q97.5"],
    fold_divergent_over_same = 10^beta_div["Estimate"],
    fold_divergent_lo = 10^beta_div["Q2.5"],
    fold_divergent_hi = 10^beta_div["Q97.5"],
    P_beta_divergent_gt_0 = p_div_gt0,

    beta_log10_flank = beta_flank["Estimate"],
    beta_flank_lo = beta_flank["Q2.5"],
    beta_flank_hi = beta_flank["Q97.5"],
    fold_per_10x_flank = 10^beta_flank["Estimate"],
    fold_per_10x_flank_lo = 10^beta_flank["Q2.5"],
    fold_per_10x_flank_hi = 10^beta_flank["Q97.5"],
    P_beta_flank_gt_0 = p_flank_gt0,

    kfold_elpd = ev$elpd,
    kfold_se = ev$se,
    sd_phylo_log10 = sd_phylo,
    sigma_log10 = sigma_res
  )
}

build_row_type <- function(fit) {
  nd <- data.frame(
    polarity_3bin = factor(factor_levels, levels = factor_levels),
    flanking_types = factor(
      "gene-gene",
      levels = levels(igr$flanking_types)
    )
  )
  mu <- fitted(fit, newdata = nd, re_formula = NA, summary = TRUE)[,
    c("Estimate", "Q2.5", "Q97.5")
  ]

  get_mu <- function(pol) {
    i <- which(nd$polarity_3bin == pol)
    list(
      mean = 10^mu[i, "Estimate"],
      lo = 10^mu[i, "Q2.5"],
      hi = 10^mu[i, "Q97.5"]
    )
  }

  fx <- brms::fixef(fit)
  beta_conv <- fx["polarity_3binconv", ]
  beta_div <- fx["polarity_3bindiv", ]

  dr <- posterior::as_draws_df(fit)
  p_conv_gt0 <- mean(dr$b_polarity_3binconv > 0)
  p_div_gt0 <- mean(dr$b_polarity_3bindiv > 0)

  sd_phylo <- as.numeric(brms::VarCorr(fit)$ncbi_taxid$sd[1])
  sigma_res <- as.numeric(brms::VarCorr(fit)$residual__$sd[1])
  ev <- get_elpd(fit)

  data.frame(
    predictor_type = "gene_type",
    flank_summary = "flanking_types",
    mu_flank = NA,
    N_regions = nrow(igr),
    N_taxa = nlevels(igr$ncbi_taxid),

    mean_same_bp = get_mu("same")$mean,
    lo_same_bp = get_mu("same")$lo,
    hi_same_bp = get_mu("same")$hi,
    mean_convergent_bp = get_mu("conv")$mean,
    lo_convergent_bp = get_mu("conv")$lo,
    hi_convergent_bp = get_mu("conv")$hi,
    mean_divergent_bp = get_mu("div")$mean,
    lo_divergent_bp = get_mu("div")$lo,
    hi_divergent_bp = get_mu("div")$hi,

    beta_log10_convergent_vs_same = beta_conv["Estimate"],
    beta_convergent_lo = beta_conv["Q2.5"],
    beta_convergent_hi = beta_conv["Q97.5"],
    fold_convergent_over_same = 10^beta_conv["Estimate"],
    fold_convergent_lo = 10^beta_conv["Q2.5"],
    fold_convergent_hi = 10^beta_conv["Q97.5"],
    P_beta_convergent_gt_0 = p_conv_gt0,

    beta_log10_divergent_vs_same = beta_div["Estimate"],
    beta_divergent_lo = beta_div["Q2.5"],
    beta_divergent_hi = beta_div["Q97.5"],
    fold_divergent_over_same = 10^beta_div["Estimate"],
    fold_divergent_lo = 10^beta_div["Q2.5"],
    fold_divergent_hi = 10^beta_div["Q97.5"],
    P_beta_divergent_gt_0 = p_div_gt0,

    # Not applicable for the categorical-only model
    beta_log10_flank = NA,
    beta_flank_lo = NA,
    beta_flank_hi = NA,
    fold_per_10x_flank = NA,
    fold_per_10x_flank_lo = NA,
    fold_per_10x_flank_hi = NA,
    P_beta_flank_gt_0 = NA,

    kfold_elpd = ev$elpd,
    kfold_se = ev$se,
    sd_phylo_log10 = sd_phylo,
    sigma_log10 = sigma_res
  )
}

# Predicted at gene-gene reference AND flank = 0 (centering point).
# Full flanking_types contrasts (all levels) are in the summary text file.
build_row_both <- function(fit, flank_summary_name, mu_flank_val) {
  nd <- data.frame(
    polarity_3bin = factor(factor_levels, levels = factor_levels),
    flank = 0,
    flanking_types = factor(
      "gene-gene",
      levels = levels(igr$flanking_types)
    )
  )
  mu <- fitted(fit, newdata = nd, re_formula = NA, summary = TRUE)[,
    c("Estimate", "Q2.5", "Q97.5")
  ]

  get_mu <- function(pol) {
    i <- which(nd$polarity_3bin == pol)
    list(
      mean = 10^mu[i, "Estimate"],
      lo = 10^mu[i, "Q2.5"],
      hi = 10^mu[i, "Q97.5"]
    )
  }

  fx <- brms::fixef(fit)
  beta_conv <- fx["polarity_3binconv", ]
  beta_div <- fx["polarity_3bindiv", ]
  beta_flank <- fx["flank", ]

  dr <- posterior::as_draws_df(fit)
  p_conv_gt0 <- mean(dr$b_polarity_3binconv > 0)
  p_div_gt0 <- mean(dr$b_polarity_3bindiv > 0)
  p_flank_gt0 <- mean(dr$b_flank > 0)

  sd_phylo <- as.numeric(brms::VarCorr(fit)$ncbi_taxid$sd[1])
  sigma_res <- as.numeric(brms::VarCorr(fit)$residual__$sd[1])
  ev <- get_elpd(fit)

  data.frame(
    predictor_type = "both",
    flank_summary = flank_summary_name,
    mu_flank = mu_flank_val,
    N_regions = nrow(igr),
    N_taxa = nlevels(igr$ncbi_taxid),

    mean_same_bp = get_mu("same")$mean,
    lo_same_bp = get_mu("same")$lo,
    hi_same_bp = get_mu("same")$hi,
    mean_convergent_bp = get_mu("conv")$mean,
    lo_convergent_bp = get_mu("conv")$lo,
    hi_convergent_bp = get_mu("conv")$hi,
    mean_divergent_bp = get_mu("div")$mean,
    lo_divergent_bp = get_mu("div")$lo,
    hi_divergent_bp = get_mu("div")$hi,

    beta_log10_convergent_vs_same = beta_conv["Estimate"],
    beta_convergent_lo = beta_conv["Q2.5"],
    beta_convergent_hi = beta_conv["Q97.5"],
    fold_convergent_over_same = 10^beta_conv["Estimate"],
    fold_convergent_lo = 10^beta_conv["Q2.5"],
    fold_convergent_hi = 10^beta_conv["Q97.5"],
    P_beta_convergent_gt_0 = p_conv_gt0,

    beta_log10_divergent_vs_same = beta_div["Estimate"],
    beta_divergent_lo = beta_div["Q2.5"],
    beta_divergent_hi = beta_div["Q97.5"],
    fold_divergent_over_same = 10^beta_div["Estimate"],
    fold_divergent_lo = 10^beta_div["Q2.5"],
    fold_divergent_hi = 10^beta_div["Q97.5"],
    P_beta_divergent_gt_0 = p_div_gt0,

    # Length-adjusted flank coefficient (controlling for type)
    beta_log10_flank = beta_flank["Estimate"],
    beta_flank_lo = beta_flank["Q2.5"],
    beta_flank_hi = beta_flank["Q97.5"],
    fold_per_10x_flank = 10^beta_flank["Estimate"],
    fold_per_10x_flank_lo = 10^beta_flank["Q2.5"],
    fold_per_10x_flank_hi = 10^beta_flank["Q97.5"],
    P_beta_flank_gt_0 = p_flank_gt0,

    kfold_elpd = ev$elpd,
    kfold_se = ev$se,
    sd_phylo_log10 = sd_phylo,
    sigma_log10 = sigma_res
  )
}

igr$length_bp <- 10^igr$log10_igr_len
med_same <- median(igr$length_bp[igr$polarity_3bin == "same"], na.rm = TRUE)
med_conv <- median(igr$length_bp[igr$polarity_3bin == "conv"], na.rm = TRUE)
med_div <- median(igr$length_bp[igr$polarity_3bin == "div"], na.rm = TRUE)

res <- rbind(
  build_row_len(m_len, flank_col, mu_flank),
  build_row_type(m_type),
  build_row_both(m_both, flank_col, mu_flank)
)
res$median_same_bp <- med_same
res$median_convergent_bp <- med_conv
res$median_divergent_bp <- med_div

res <- cbind(
  res,
  rbind(
    as.data.frame(r2_len),
    as.data.frame(r2_type),
    as.data.frame(r2_both)
  )
)
rownames(res) <- NULL

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
