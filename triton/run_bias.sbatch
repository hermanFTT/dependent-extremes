#!/bin/bash
# =============================================================================
# triton/run_bias.sbatch -- run the WHOLE bias_sensitivity grid in ONE job.
#
#     sbatch triton/run_bias.sbatch
#
# Why one job rather than 76 Slurm jobs: a single cell takes ~4 min, so the
# entire 76-cell grid is only ~5 CPU-hours. Spread over 16 cores that is ~20
# min of wall time -- less than the queue wait for 76 separate submissions.
# Snakemake runs locally inside this job and parallelises across the cores.
#
# Use the Slurm executor (profiles/triton-bias) instead only when you scale
# n_rep or the grid up by an order of magnitude. See triton/README.md.
# =============================================================================

#SBATCH --job-name=bias_sens
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=6G
#SBATCH --output=logs/slurm-bias-%j.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# NB: logs/ must already exist at `sbatch` time -- Slurm opens the --output
# file above before this script runs. bootstrap.sh and logs/.gitkeep ensure
# that. This mkdir only covers the per-rule log subdirectories.
mkdir -p logs

# shellcheck source=/dev/null
source triton/env.sh

echo "[run_bias] $SLURM_CPUS_PER_TASK cores, started $(date)"

snakemake -s bias_sensitivity.smk \
    --cores "$SLURM_CPUS_PER_TASK" \
    --latency-wait 30 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds

echo "[run_bias] finished $(date)"
