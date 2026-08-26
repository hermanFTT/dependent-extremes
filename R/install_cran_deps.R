#!/usr/bin/env Rscript
# =============================================================================
# install_cran_deps.R -- install R dependencies that conda does not provide.
#
#   Rscript R/install_cran_deps.R
#
# Idempotent: only installs what is actually missing, so it is safe to re-run
# (and safe to call from a Snakemake rule).
#
# POT is the only hard case -- it is not on conda-forge/bioconda, so it must
# come from CRAN.  The others are listed too so this script alone is enough to
# bootstrap a plain system R with no conda at all.
# =============================================================================

repos <- getOption("repos")
if (is.null(repos[["CRAN"]]) || repos[["CRAN"]] == "@CRAN@")
  repos <- c(CRAN = "https://cloud.r-project.org")

required <- c(
  POT         = "1.1.11",   # simmc(): Markov EV-copula simulation  [CRAN only]
  ggplot2     = "3.5.0",
  scales      = "1.3.0",
  matrixStats = "1.0.0",
  HDInterval  = "0.2.4"
)

have <- vapply(names(required), function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) "-" else v
}, character(1))

needed <- names(required)[
  vapply(names(required), function(p) {
    v <- tryCatch(utils::packageVersion(p), error = function(e) NULL)
    is.null(v) || v < package_version(required[[p]])
  }, logical(1))
]

cat("R ", as.character(getRversion()), " | library: ", .libPaths()[1L], "\n", sep = "")
for (p in names(required))
  cat(sprintf("  %-12s have %-8s need >= %s%s\n", p, have[[p]], required[[p]],
              if (p %in% needed) "   <-- installing" else ""))

if (!length(needed)) {
  cat("\nAll dependencies already satisfied.\n")
} else {
  cat("\nInstalling from ", repos[["CRAN"]], " ...\n", sep = "")
  utils::install.packages(needed, repos = repos)

  still <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still))
    stop("failed to install: ", paste(still, collapse = ", "),
         "\n(POT needs a working C/Fortran toolchain -- `conda install compilers`)")
  cat("\nDone.\n")
}
