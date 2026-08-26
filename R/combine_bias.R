#!/usr/bin/env Rscript
# =============================================================================
# combine_bias.R -- bind the per-cell RDS files into one tidy long table.
#
# Usage:
#   Rscript R/combine_bias.R --out_rds=results/bias_sensitivity/bias_results.rds \
#           --out_csv=results/bias_sensitivity/bias_results.csv  <cell files...>
#
# Output: one row per (dep_type, dep_par, k, n, r) x approach x parameter,
# which is what the plotting script filters and groups by.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

opt <- list(out_rds = "results/bias_sensitivity/bias_results.rds",
            out_csv = "results/bias_sensitivity/bias_results.csv")

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

res$approach  <- factor(res$approach,
                        levels = c("Conventional", "Weighted GPD",
                                   "Cluster maxima (POT)"))
res$parameter <- factor(res$parameter, levels = c("k", "sigma", "zm"))

res <- res[order(res$dep_type, res$dep_par, res$k_true, res$n,
                 res$r_label, res$parameter, res$approach), ]
rownames(res) <- NULL

dir.create(dirname(opt$out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(res, opt$out_rds)
utils::write.csv(res, opt$out_csv, row.names = FALSE)

cat(sprintf("[combine] %d cells -> %d rows\n", length(files), nrow(res)))
cat(sprintf("[combine] grid: dep_type=%s | dep_par=%d | k=%d | n=%d | r=%s\n",
            paste(unique(res$dep_type), collapse = ","),
            length(unique(res$dep_par)), length(unique(res$k_true)),
            length(unique(res$n)),
            paste(unique(res$r_label), collapse = ",")))
cat("[combine] wrote ", opt$out_rds, " and ", opt$out_csv, "\n", sep = "")

nbad <- sum(res$n_ok < res$S)
if (nbad)
  cat(sprintf("[combine] NOTE: %d row(s) lost replicates (min n_ok = %d of S = %d)\n",
              nbad, min(res$n_ok), max(res$S)))
