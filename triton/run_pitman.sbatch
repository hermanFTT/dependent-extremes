#!/bin/bash
# =============================================================================
# triton/run_pitman.sbatch -- run the WHOLE pitman_closeness grid in ONE job.
#
#     sbatch triton/run_pitman.sbatch
#
# 20 grid points x ~3.5 min = ~1.2 CPU-hours; ~5 min wall on 16 cores.
# =============================================================================
#SBATCH --job-name=pitman
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --output=logs/slurm-pitman-%j.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# NB: logs/ must already exist at `sbatch` time -- Slurm opens the --output
# file above before this script runs. bootstrap.sh and logs/.gitkeep ensure
# that. This mkdir only covers the per-rule log subdirectories.
mkdir -p logs

# shellcheck source=/dev/null
source triton/env.sh

echo "[run_pitman] $SLURM_CPUS_PER_TASK cores, started $(date)"

snakemake -s pitman_closeness.smk \
    --cores "$SLURM_CPUS_PER_TASK" \
    --latency-wait 30 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds

echo "[run_pitman] finished $(date)"
