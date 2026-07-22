#!/usr/bin/env Rscript

options(warn = 1)
required_packages <- c("brms", "posterior")
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

script_start <- Sys.time()
log_time("Script started (divergent - convergent posterior contrast)")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Exactly 5 arguments required\n",
    "Usage: Rscript brms_contrast.R <tree.nwk> <igr.tsv> ",
    "<summary.txt> <models.rds> <results.tsv>\n",
    "(tree and igr are accepted for runner compatibility but not used; ",
    "contrasts are read from the cached model fits.)"
  )
}

# args[1] (tree) and args[2] (igr) intentionally unused.
results_txt <- args[3]
out_name <- args[4]
results_tsv <- args[5]

cache_flank <- sub("\\.rds$", "_cache_flank.rds", out_name)
cache_type <- sub("\\.rds$", "_cache_type.rds", out_name)
cache_both <- sub("\\.rds$", "_cache_both.rds", out_name)
cache_sigma <- sub("\\.rds$", "_cache_sigma.rds", out_name)

get_b_draws <- function(fit) {}

summarise_contrast <- function(draws, link = c("log10", "log")) {
  link <- match.arg(link)
  f <- if (link == "log10") function(z) 10^z else exp
  est <- stats::median(draws)
  lo <- unname(stats::quantile(draws, 0.025))
  hi <- unname(stats::quantile(draws, 0.975))
  ess <- tryCatch(posterior::ess_bulk(draws), error = function(e) NA_real_)
  list(
    delta = est,
    delta_lo = lo,
    delta_hi = hi,
    ratio = f(est),
    ratio_lo = f(lo),
    ratio_hi = f(hi),
    p_div_gt_conv = mean(draws > 0),
    p_conv_gt_div = mean(draws < 0),
    ess_bulk = ess
  )
}

make_row <- function(model, parameter, scale, s, n) {
  data.frame(
    model = model,
    parameter = parameter,
    scale = scale,
    delta = s$delta,
    delta_lo = s$delta_lo,
    delta_hi = s$delta_hi,
    ratio_div_over_conv = s$ratio,
    ratio_lo = s$ratio_lo,
    ratio_hi = s$ratio_hi,
    P_div_gt_conv = s$p_div_gt_conv,
    P_conv_gt_div = s$p_conv_gt_div,
    ess_bulk_contrast = s$ess_bulk,
    N_regions = n,
    stringsAsFactors = FALSE
  )
}

conv_name <- "b_polarity_3binconv"
div_name <- "b_polarity_3bindiv"

mean_contrast_row <- function(cache, label) {
  if (!file.exists(cache)) {
    message("NOTE: cache not found, skipping: ", cache)
    return(NULL)
  }
  log_time(paste("Reading:", cache))
  fit <- readRDS(cache)
  dr <- posterior::as_draws_df(fit, variable = "^b_", regex = TRUE)
  if (!all(c(conv_name, div_name) %in% names(dr))) {
    stop(
      "Missing coefficients in ",
      cache,
      " (expected ",
      conv_name,
      " and ",
      div_name,
      "). Check polarity_3bin factor coding."
    )
  }
  contrast <- dr[[div_name]] - dr[[conv_name]] # log10 length scale
  s <- summarise_contrast(contrast, link = "log10")
  make_row(
    label,
    "beta_div - beta_conv (log10 length)",
    "log10",
    s,
    nrow(fit$data)
  )
}

rows <- list()
rows[["both"]] <- mean_contrast_row(cache_both, "both (adjusted mean; PRIMARY)")
rows[["len"]] <- mean_contrast_row(
  cache_flank,
  "length-adjusted mean (robustness)"
)
rows[["type"]] <- mean_contrast_row(
  cache_type,
  "type-adjusted mean (robustness)"
)

if (is.null(rows[["both"]])) {
  stop(
    "The primary 'both' model cache was not found: ",
    cache_both,
    "\nRun brms_type_length.R for this group first."
  )
}

sconv_name <- "b_sigma_polarity_3binconv"
sdiv_name <- "b_sigma_polarity_3bindiv"

log_time(paste("Reading:", cache_sigma))
m_sigma <- readRDS(cache_sigma)
dr_s <- posterior::as_draws_df(m_sigma, variable = "^b_", regex = TRUE)
contrast_s <- dr_s[[sdiv_name]] - dr_s[[sconv_name]] # log SD scale
ss <- summarise_contrast(contrast_s, link = "log")
rows[["sigma"]] <- make_row(
  "sigma (dispersion)",
  "sigma_div - sigma_conv (log SD)",
  "log",
  ss,
  nrow(m_sigma$data)
)

res <- do.call(
  rbind,
  rows[c("both", "len", "type", "sigma")[
    c("both", "len", "type", "sigma") %in% names(rows)
  ]]
)
rownames(res) <- NULL

log_time(paste("Saving human-readable summary to:", results_txt))
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

  cat("Posterior contrast: divergent vs convergent, within dataset\n")
  cat("Contrast formed from joint draws (accounts for posterior covariance).\n")
  cat("Positive delta => divergent effect exceeds convergent effect.\n\n")

  for (i in seq_len(nrow(res))) {
    r <- res[i, ]
    unit <- if (r$scale == "log10") {
      "length ratio (div/conv)"
    } else {
      "SD ratio (div/conv)"
    }
    cat(sprintf("[%s]\n", r$model))
    cat(sprintf("  %s\n", r$parameter))
    cat(sprintf(
      "  delta (%s scale) = %+.4f (%+.4f, %+.4f)\n",
      r$scale,
      r$delta,
      r$delta_lo,
      r$delta_hi
    ))
    cat(sprintf(
      "  %s = %.3f (%.3f, %.3f)\n",
      unit,
      r$ratio_div_over_conv,
      r$ratio_lo,
      r$ratio_hi
    ))
    cat(sprintf(
      "  P(divergent > convergent) = %.4f;  P(convergent > divergent) = %.4f\n",
      r$P_div_gt_conv,
      r$P_conv_gt_div
    ))
    cat(sprintf(
      "  bulk-ESS of contrast = %.0f;  N_regions = %d\n\n",
      r$ess_bulk_contrast,
      r$N_regions
    ))
  }
}
write_summary(results_txt)

log_time(paste("Saving machine-readable row(s) to:", results_tsv))
write.table(
  res,
  file = results_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(sprintf(
  "Total runtime: %.2f min",
  as.numeric(difftime(Sys.time(), script_start, units = "mins"))
))
message("Done!")
