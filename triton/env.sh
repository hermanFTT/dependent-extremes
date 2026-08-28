# =============================================================================
# triton/env.sh -- environment for the dependent-extremes pipeline on Triton.
#
# SOURCE it (do not execute it) from the repo root:
#
#     source triton/env.sh
#
# Lmod modules are shell functions that set environment variables, so they only
# propagate to child processes -- sourcing is required, running is not enough.
# Source this in every new shell before launching snakemake, and from inside
# tmux before a long run.
# =============================================================================

# --- R ----------------------------------------------------------------------
# scicomp-r-env ships a large set of CRAN packages already compiled, so only a
# handful (POT, HDInterval, ...) have to be built into the private library.
# For a bare R instead, replace the next line with:  module load r
module load scicomp-r-env

# POT has C/Fortran sources and must be compiled from CRAN. scicomp-r-env's
# own bundled compiler is only what R itself was built with -- it is missing
# cc1, its actual backend -- so a real, separate toolchain module is needed.
module load gcc gmake binutils

# Private library for packages the module does not provide. Keyed by R
# major.minor so switching R versions can never mix incompatible builds.
_R_MM=$(Rscript -e 'cat(paste(R.version$major, sub("[.].*", "", R.version$minor), sep = "."))')
export R_LIBS_USER="$HOME/.local/lib/R/${_R_MM}"
mkdir -p "$R_LIBS_USER"

# The pipeline parallelises with parallel::mclapply and Snakemake gives each
# job exactly `threads` cores. Pin the BLAS to one thread so it does not
# separately oversubscribe the Slurm allocation.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# --- Snakemake --------------------------------------------------------------
# Conda envs live in WRKDIR: $HOME is only 10 GB and this env is ~1 GB.
#
# IMPORTANT: we do NOT `module load mamba` here. On Triton it shares an Lmod
# family with scicomp-r-env, so loading it silently UNLOADS R ("Lmod is
# automatically replacing scicomp-r-env with mamba"), which then breaks every
# `Rscript` call inside the snakemake rules. The mamba module is only ever
# needed once, to CREATE the env (see triton/bootstrap.sh) -- to USE an
# already-built env we just need its own bin/ on PATH, no module required.
export CONDA_ENVS_PATH="$WRKDIR/conda-envs"
export CONDA_PKGS_DIRS="$WRKDIR/conda-pkgs"
mkdir -p "$CONDA_ENVS_PATH" "$CONDA_PKGS_DIRS"

# Exported (not unset) -- triton/bootstrap.sh needs this to know where to
# create/find the env, including under `set -u`.
export SNAKEMAKE_ENV="$CONDA_ENVS_PATH/dependent-extremes"
if [ -d "$SNAKEMAKE_ENV/bin" ]; then
    export PATH="$SNAKEMAKE_ENV/bin:$PATH"
fi

echo "[env] R           : $(command -v Rscript) (${_R_MM})"
echo "[env] R_LIBS_USER : $R_LIBS_USER"
echo "[env] snakemake   : $(command -v snakemake 2>/dev/null || echo '(not installed -- run triton/bootstrap.sh)')"

unset _R_MM
