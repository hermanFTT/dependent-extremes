#!/usr/bin/env Rscript
# =============================================================================
# run_bias.R -- run ONE grid cell of the bias / sensitivity study.
#
# A cell is (dep_type, dep_par, k, n, r); this is the unit Snakemake
# parallelises over.
#
# Usage:
#   Rscript R/run_bias.R --dep_type=log --dep_par=0.35 --k=0.2 --n=4000 \
#           --r=auto --n_rep=1000 --out=results/bias_sensitivity/raw/cell.rds
#
# Any option below can be overridden as --name=value.  `--r=auto` means the
# run length is chosen per replicate by auto_run_length().
# =============================================================================

defaults <- list(
  # --- the cell ------------------------------------------------------------
  dep_type     = "log",        # "ar1" | "log" | "nlog"
  dep_par      = 0.35,         # phi (ar1) or alpha (log/nlog)
  k            = 0.2,          # true shape
  n            = 4000L,        # series length
  r            = "auto",       # run length: "auto" or an integer
  # --- Monte-Carlo ---------------------------------------------------------
  n_rep        = 1000L,        # S
  sigma        = 1,
  p_thresh     = 0.95,
  threshold    = "empirical",  # "empirical" | "theoretical"
  m_ret        = 1000,         # return period m for z_m
  # --- coverage ------------------------------------------------------------
  coverage     = TRUE,
  cred_mass    = 0.95,
  interval     = "hpd",        # "hpd" | "eti"
  n_draws      = 10000L,
  # --- fitting -------------------------------------------------------------
  n_grid       = 50L,
  min_grid_pts = 100L,
  wip          = FALSE,
  approaches   = "Conventional,Weighted GPD,Cluster maxima (POT)",
  # --- plumbing ------------------------------------------------------------
  seed         = 20250L,
  ncores       = 1L,
  gpdfit       = "weighted_gpdfit.R",
  core         = "R/bias_core.R",
  out          = "results/bias_sensitivity/raw/cell.rds"
)

parse_args <- function(defaults) {
  opts <- defaults
  for (a in commandArgs(trailingOnly = TRUE)) {
    if (!grepl("^--[^=]+=", a)) {
      warning("ignoring unparsable argument: ", a, call. = FALSE); next
    }
    key <- sub("^--([^=]+)=.*$", "\\1", a)
    val <- sub("^--[^=]+=", "", a)
    if (!key %in% names(defaults)) stop("unknown option --", key, call. = FALSE)
    d <- defaults[[key]]
    opts[[key]] <-
      if (is.logical(d))      as.logical(val)
      else if (is.integer(d)) as.integer(val)
      else if (is.numeric(d)) as.numeric(val)
      else                    val
  }
  opts
}

opt <- parse_args(defaults)

source(opt$core)
source_gpd_functions(opt$gpdfit)

# run length: "auto" -> NULL
r_val <- if (identical(tolower(opt$r), "auto")) NULL else as.integer(opt$r)

approaches <- trimws(strsplit(opt$approaches, ",", fixed = TRUE)[[1L]])

# --- per-cell seed --------------------------------------------------------
# Deterministic function of the master seed and the cell identity, so each
# cell is independent, reproducible, and unaffected by the grid's shape
# (adding a new n or k does not renumber anybody else's stream).
cell_id <- sprintf("%s|%.10g|%.10g|%d|%s",
                   opt$dep_type, opt$dep_par, opt$k, opt$n,
                   if (is.null(r_val)) "auto" else r_val)
chars     <- utf8ToInt(cell_id)
cell_seed <- opt$seed + as.integer(sum(chars * seq_along(chars)) %% 1000003L)

cat(sprintf("[run_bias] cell %s | S = %d | %d core%s | seed %d\n",
            cell_id, opt$n_rep, opt$ncores,
            if (opt$ncores > 1L) "s" else "", cell_seed))

res <- bias_one_cell(
  dep_type     = opt$dep_type,
  dep_par      = opt$dep_par,
  k            = opt$k,
  n            = opt$n,
  r            = r_val,
  S            = opt$n_rep,
  sigma        = opt$sigma,
  p_thresh     = opt$p_thresh,
  threshold    = opt$threshold,
  m            = opt$m_ret,
  approaches   = approaches,
  coverage     = opt$coverage,
  cred_mass    = opt$cred_mass,
  interval     = opt$interval,
  n_draws      = opt$n_draws,
  n_grid       = opt$n_grid,
  min_grid_pts = opt$min_grid_pts,
  wip          = opt$wip,
  seed         = cell_seed,
  ncores       = opt$ncores,
  verbose      = TRUE
)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(res, opt$out)
cat("[run_bias] wrote ", opt$out, "\n", sep = "")
