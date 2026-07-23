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
log_time("Script started (divergent verification)")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: brms_divergence.R <polarity_dir> <polarity_out> <type_length_dir> <type_length_out>"
  )
}

polarity_path <- args[1]
polarity_out_name <- args[2]
type_length_path <- args[3]
type_length_out_name <- args[4]

check_one <- function(path) {
  log_time(paste("Checking", path))
  fit <- readRDS(path)
  np <- brms::nuts_params(fit)

  # read max_treedepth off the fit rather than assuming
  max_td <- tryCatch(
    {
      ctl <- attr(fit$fit@sim$samples[[1]], "args")$control
      if (is.null(ctl$max_treedepth)) 10 else ctl$max_treedepth
    },
    error = function(e) 10
  )

  div <- np$Value[np$Parameter == "divergent__"]
  td <- np$Value[np$Parameter == "treedepth__"]

  drw <- posterior::as_draws_array(fit)
  ndr <- posterior::ndraws(drw)

  pars <- setdiff(posterior::variables(drw), c("lp__", "lprior", "Intercept"))
  s <- posterior::summarise_draws(
    posterior::subset_draws(drw, variable = pars),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  data.frame(
    file = basename(path),
    n_nuts_rows = length(div),
    max_treedepth = max_td,
    divergences = sum(div),
    pct_div = round(100 * sum(div) / ndr, 3),
    treedepth_hits = sum(td >= max_td),
    max_td_reached = max(td),
    draws = ndr,
    max_rhat = round(max(s$rhat, na.rm = TRUE), 4),
    n_rhat_gt_101 = sum(s$rhat > 1.01, na.rm = TRUE),
    min_ess_bulk = round(min(s$ess_bulk, na.rm = TRUE)),
    min_ess_tail = round(min(s$ess_tail, na.rm = TRUE)),
    worst_par = s$variable[which.max(s$rhat)],
    n_params = nrow(s)
  )
}

polarity_models <- list.files(
  polarity_path,
  pattern = "cache.*\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

write.table(
  do.call(rbind, lapply(polarity_models, check_one)),
  file = polarity_out_name,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

type_length_models <- list.files(
  type_length_path,
  pattern = "cache.*\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

write.table(
  do.call(rbind, lapply(type_length_models, check_one)),
  file = type_length_out_name,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
