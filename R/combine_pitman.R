#!/usr/bin/env Rscript
# =============================================================================
# combine_pitman.R -- bind the per-grid-point RDS files into one tidy table.
#
# Usage:
#   Rscript R/combine_pitman.R --out_rds=results/pitman_results.rds \
#           --out_csv=results/pitman_results.csv  results/raw/pc_*.rds
#
# Output: a long data.frame, 2 rows per grid point (parameter = "k" / "sigma"),
# which is exactly the shape the plotting script wants.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

opt <- list(out_rds = "results/pitman_results.rds",
            out_csv = "results/pitman_results.csv")

is_opt <- grepl("^--[^=]+=", args)
for (a in args[is_opt]) {
  key <- sub("^--([^=]+)=.*$", "\\1", a)
  if (!key %in% names(opt)) stop("unknown option --", key)
  opt[[key]] <- sub("^--[^=]+=", "", a)
}

files <- args[!is_opt]
if (!length(files)) stop("no input .rds files given")
missing <- files[!file.exists(files)]
if (length(missing)) stop("missing input files: ", paste(missing, collapse = ", "))

res <- do.call(rbind, lapply(files, readRDS))

# stable, human-friendly ordering
res <- res[order(res$parameter, res$k_true), ]
rownames(res) <- NULL

dir.create(dirname(opt$out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(res, opt$out_rds)
utils::write.csv(res, opt$out_csv, row.names = FALSE)

cat(sprintf("[combine] %d files -> %d rows (%d grid points x %d parameters)\n",
            length(files), nrow(res),
            length(unique(res$k_true)), length(unique(res$parameter))))
cat("[combine] wrote ", opt$out_rds, " and ", opt$out_csv, "\n", sep = "")

if (any(res$n_ok < res$S))
  cat(sprintf("[combine] NOTE: %d row(s) had failed replicates (min n_ok = %d of S = %d)\n",
              sum(res$n_ok < res$S), min(res$n_ok), max(res$S)))
