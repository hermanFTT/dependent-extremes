#!/usr/bin/env Rscript
# =============================================================================
# run_pitman.R -- run ONE grid point of the Pitman-closeness experiment.
#
# This is the unit of work Snakemake parallelises over: one job per true
# shape value on the grid.
#
# Usage:
#   Rscript R/run_pitman.R --index=7 --out=results/raw/pc_7.rds
#   Rscript R/run_pitman.R --index=1 --S=5000 --dep=markov --alpha=0.35 \
#           --t1="Weighted GPD" --t2="Conventional" --out=/tmp/pc_1.rds
#
# Every option below can be overridden as --name=value.
# =============================================================================

suppressWarnings(suppressMessages({
  ok <- requireNamespace("parallel", quietly = TRUE)
}))

# --- defaults (types here drive the coercion of command-line values) --------
defaults <- list(
  index        = 1L,           # which grid point (1 .. n_k)
  n_k          = 20L,          # number of grid points
  k_min        = -0.5,
  k_max        =  0.5,
  S            = 100L,         # replicates per grid point
  n            = 2000L,        # series length
  dep          = "markov",     # "ar1" | "markov"
  phi          = 0.95,         # AR(1) coefficient
  alpha        = 0.35,         # Markov dependence
  dep_model    = "log",
  sigma        = 1,
  p_thresh     = 0.95,
  threshold    = "empirical",  # "theoretical" | "empirical"
  run          = NA_integer_,  # NA -> auto_run_length() per replicate
  n_grid       = 50L,
  min_grid_pts = 100L,
  wip          = FALSE,
  t1           = "Weighted GPD",
  t2           = "Conventional",
  seed         = 4415L,
  ncores       = 1L,           # cores WITHIN this grid point
  gpdfit       = "weighted_gpdfit.R",
  core         = "R/pitman_core.R",
  out          = "results/raw/pc_1.rds"
)

# --- parse --key=value ------------------------------------------------------
parse_args <- function(defaults) {
  opts <- defaults
  args <- commandArgs(trailingOnly = TRUE)
  for (a in args) {
    if (!grepl("^--[^=]+=", a)) {
      warning("ignoring unparsable argument: ", a, call. = FALSE); next
    }
    key <- sub("^--([^=]+)=.*$", "\\1", a)
    val <- sub("^--[^=]+=", "", a)
    if (!key %in% names(defaults))
      stop("unknown option --", key, call. = FALSE)

    d <- defaults[[key]]
    opts[[key]] <-
      if (is.logical(d))                       as.logical(val)
      else if (is.integer(d))                  as.integer(val)
      else if (is.numeric(d))                  as.numeric(val)
      else                                     val
  }
  opts
}

opt <- parse_args(defaults)

source(opt$core)                       # pitman_one_k(), make_k_grid(), ...
source_gpd_functions(opt$gpdfit)       # estimators, without the scratch block

k_grid <- make_k_grid(opt$k_min, opt$k_max, opt$n_k)
if (opt$index < 1L || opt$index > length(k_grid))
  stop("--index=", opt$index, " outside 1..", length(k_grid))
k_true <- k_grid[opt$index]

cat(sprintf("[run_pitman] grid point %d/%d  k = %+.4f  S = %d  %s vs %s  (%d core%s)\n",
            opt$index, opt$n_k, k_true, opt$S, opt$t1, opt$t2,
            opt$ncores, if (opt$ncores > 1L) "s" else ""))

res <- pitman_one_k(
  k_true       = k_true,
  S            = opt$S,
  n            = opt$n,
  dep          = opt$dep,
  phi          = opt$phi,
  alpha        = opt$alpha,
  dep_model    = opt$dep_model,
  sigma        = opt$sigma,
  p_thresh     = opt$p_thresh,
  threshold    = opt$threshold,
  run          = if (is.na(opt$run)) NULL else opt$run,
  n_grid       = opt$n_grid,
  min_grid_pts = opt$min_grid_pts,
  wip          = opt$wip,
  t1           = opt$t1,
  t2           = opt$t2,
  # distinct but reproducible stream per grid point
  seed         = opt$seed + opt$index - 1L,
  ncores       = opt$ncores,
  verbose      = TRUE
)

res$index <- opt$index

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(res, opt$out)
cat("[run_pitman] wrote ", opt$out, "\n", sep = "")
