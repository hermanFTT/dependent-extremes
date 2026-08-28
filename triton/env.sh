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

# gcc: POT has C/Fortran sources and must be compiled from CRAN.
module load gcc

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
module load mamba
export CONDA_ENVS_PATH="$WRKDIR/conda-envs"
mkdir -p "$CONDA_ENVS_PATH"

# Activate only if it exists, so this file is safe to source before
# triton/bootstrap.sh has ever been run.
if [ -d "$CONDA_ENVS_PATH/dependent-extremes" ]; then
    source activate dependent-extremes
fi

echo "[env] R           : $(command -v Rscript) (${_R_MM})"
echo "[env] R_LIBS_USER : $R_LIBS_USER"
echo "[env] snakemake   : $(command -v snakemake 2>/dev/null || echo '(not installed -- run triton/bootstrap.sh)')"

unset _R_MM
