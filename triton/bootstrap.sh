#!/bin/bash
# =============================================================================
# triton/bootstrap.sh -- one-time install of everything the pipeline needs.
#
# Run ONCE from anywhere in the repo, on the Triton LOGIN node:
#
#     bash triton/bootstrap.sh
#
# Installs:
#   1. a conda env with snakemake 9.4.1 + the Slurm executor plugin
#      (Snakemake >= 8 dropped --cluster in favour of executor plugins)
#   2. the R packages POT / matrixStats / HDInterval / ggplot2 / scales
#      into the private library defined by triton/env.sh
#
# Safe to re-run: both steps skip whatever is already present.
# =============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
echo "[bootstrap] repo root: $PWD"

# Slurm opens the --output file BEFORE running the job script, so logs/ must
# exist at submission time or the job dies with no diagnosable error. Git does
# not track empty directories, so a fresh clone will not have it.
mkdir -p logs

# shellcheck source=/dev/null
source triton/env.sh

# --- 1. snakemake env -------------------------------------------------------
if [ -d "$SNAKEMAKE_ENV" ]; then
    echo "[bootstrap] conda env already exists, skipping creation"
else
    echo "[bootstrap] creating conda env (takes a few minutes) ..."
    # `module load mamba` is used ONLY here, ONCE. It shares an Lmod family
    # with scicomp-r-env and unloads it as a side effect ("Lmod is
    # automatically replacing scicomp-r-env with mamba") -- harmless at this
    # point since we are not using R yet, but we must reload it below before
    # any further Rscript call.
    module load mamba
    # pulp is pinned: >= 2.9 breaks snakemake 9.4's scheduler ILP backend.
    mamba create -y -n dependent-extremes \
        -c conda-forge -c bioconda \
        python=3.12 \
        snakemake=9.4.1 \
        'pulp=2.8' \
        snakemake-executor-plugin-slurm
    # Restore R (this in turn unloads mamba again -- fine, the env is now on
    # disk and triton/env.sh reaches it via PATH, not via the module).
    module load scicomp-r-env
fi

export PATH="$SNAKEMAKE_ENV/bin:$PATH"
echo "[bootstrap] snakemake $(snakemake --version)"

# --- 2. R packages ----------------------------------------------------------
echo "[bootstrap] installing R packages into $R_LIBS_USER ..."
Rscript R/install_cran_deps.R

# --- 3. verify --------------------------------------------------------------
echo "[bootstrap] verifying that every package loads ..."
Rscript -e '
pkgs <- c("POT", "matrixStats", "HDInterval", "ggplot2", "scales")
for (p in pkgs)
  cat(sprintf("  %-12s %s\n", p,
      tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")))
bad <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(bad))
  stop("could not load: ", paste(bad, collapse = ", "),
       "\n(POT needs a working C/Fortran toolchain -- is `module load gcc` active?)")
cat("All R dependencies OK.\n")'

# --- 4. smoke test ----------------------------------------------------------
echo "[bootstrap] smoke-testing one tiny cell (this may take a minute) ..."
Rscript R/run_bias.R --dep_type=log --dep_par=0.5 --k=0.3 --n=2000 \
    --r=auto --n_rep=2 --coverage=FALSE \
    --out=/tmp/bootstrap_smoke_$USER.rds >/dev/null
rm -f "/tmp/bootstrap_smoke_$USER.rds"
echo "[bootstrap] smoke test passed."

cat <<'EOF'

[bootstrap] Done.

In every new shell, from the repo root:

    source triton/env.sh

Then, inside tmux:

    snakemake -s bias_sensitivity.smk    --profile profiles/triton-bias   -n
    snakemake -s bias_sensitivity.smk    --profile profiles/triton-bias

    snakemake -s pitman_closeness.smk    --profile profiles/triton-pitman -n
    snakemake -s pitman_closeness.smk    --profile profiles/triton-pitman

Calibrate the runtime first -- see triton/README.md.

EOF
