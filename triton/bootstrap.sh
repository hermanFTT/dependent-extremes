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

# --- 1. snakemake env + a full matching C/Fortran toolchain -----------------
# POT has C/Fortran sources. scicomp-r-env's own bundled compiler driver
# (x86_64-conda-linux-gnu-cc) is only what R itself was built with -- it is
# missing cc1, its actual backend, and fails with
# "cannot execute 'cc1': posix_spawnp: No such file or directory".
# Fix (already the plan in environment.yaml): conda-forge's `compilers`
# metapackage, a complete, self-contained gcc/gfortran. We install it into
# THIS conda env; since triton/env.sh always puts this env's bin/ ahead of
# scicomp-r-env's on PATH, its same-named compiler binary shadows the broken
# one -- R keeps calling "x86_64-conda-linux-gnu-cc" and now gets a working one.
if [ -d "$SNAKEMAKE_ENV" ]; then
    echo "[bootstrap] conda env already exists"
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
        snakemake-executor-plugin-slurm \
        compilers
    # Restore R (this in turn unloads mamba again -- fine, the env is now on
    # disk and triton/env.sh reaches it via PATH, not via the module).
    module load scicomp-r-env
fi

export PATH="$SNAKEMAKE_ENV/bin:$PATH"

# Env pre-existed (e.g. created before this script added `compilers` above)
# but lacks a working compiler: patch it in now. Idempotent -- a no-op once
# present.
if [ ! -x "$SNAKEMAKE_ENV/bin/x86_64-conda-linux-gnu-cc" ]; then
    echo "[bootstrap] env exists but has no working compiler -- installing 'compilers' ..."
    module load mamba
    mamba install -y -n dependent-extremes -c conda-forge compilers
    module load scicomp-r-env
    export PATH="$SNAKEMAKE_ENV/bin:$PATH"
fi

echo "[bootstrap] snakemake $(snakemake --version)"

# --- 2. force R to use that toolchain, via ~/.R/Makevars --------------------
# CONFIRMED NECESSARY: shadowing the compiler name earlier on PATH was not
# enough -- R's Makeconf appears to invoke the compiler by an absolute path
# baked in at build time, not a PATH lookup. A user Makevars file is always
# read by R and always overrides Makeconf's CC/CXX/FC, regardless of how
# Makeconf itself set them, so this works independent of that mechanism.
CC_BIN=$(command -v x86_64-conda-linux-gnu-cc  || command -v x86_64-conda-linux-gnu-gcc || true)
CXX_BIN=$(command -v x86_64-conda-linux-gnu-c++ || command -v x86_64-conda-linux-gnu-g++ || true)
FC_BIN=$(command -v x86_64-conda-linux-gnu-gfortran || true)
AR_BIN=$(command -v x86_64-conda-linux-gnu-ar || true)
RANLIB_BIN=$(command -v x86_64-conda-linux-gnu-ranlib || true)

if [ -z "$CC_BIN" ] || [ -z "$FC_BIN" ]; then
    echo "[bootstrap] ERROR: expected compiler binaries not found on PATH."
    echo "  looked for x86_64-conda-linux-gnu-{cc,gcc,c++,g++,gfortran,ar,ranlib} in:"
    echo "    $SNAKEMAKE_ENV/bin"
    echo "  actual contents matching 'gnu-':"
    ls "$SNAKEMAKE_ENV"/bin | grep -i gnu- || echo "    (none found)"
    exit 1
fi

echo "[bootstrap] pinning R's compiler via ~/.R/Makevars:"
echo "    CC = $CC_BIN"
echo "    FC = $FC_BIN"

mkdir -p "$HOME/.R"
if [ -f "$HOME/.R/Makevars" ] && ! grep -q "dependent-extremes bootstrap" "$HOME/.R/Makevars"; then
    cp "$HOME/.R/Makevars" "$HOME/.R/Makevars.bak.$(date +%s)"
    echo "[bootstrap] existing ~/.R/Makevars backed up (was not ours)"
fi

cat > "$HOME/.R/Makevars" <<EOF
# Written by dependent-extremes bootstrap ($(date))
# Overrides scicomp-r-env's own incomplete compiler driver with the
# conda-forge 'compilers' toolchain installed into $SNAKEMAKE_ENV.
CC     = $CC_BIN
CXX    = $CXX_BIN
FC     = $FC_BIN
F77    = $FC_BIN
AR     = $AR_BIN
RANLIB = $RANLIB_BIN
EOF

# --- 3. R packages ----------------------------------------------------------
echo "[bootstrap] installing R packages into $R_LIBS_USER ..."
Rscript R/install_cran_deps.R

# --- 4. verify --------------------------------------------------------------
echo "[bootstrap] verifying that every package loads ..."
Rscript -e '
pkgs <- c("POT", "matrixStats", "HDInterval", "ggplot2", "scales")
for (p in pkgs)
  cat(sprintf("  %-12s %s\n", p,
      tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")))
bad <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(bad))
  stop("could not load: ", paste(bad, collapse = ", "),
       "\n(check ~/.R/Makevars points at a working compiler -- see",
       "\ntriton/bootstrap.sh step 2)")
cat("All R dependencies OK.\n")'

# --- 5. smoke test ----------------------------------------------------------
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
