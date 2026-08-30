#!/bin/bash
#SBATCH --job-name=pitman_closeness
#SBATCH --time=04:00:00                   # UNVERIFIED starting guess -- benchmark first! This is
                                           # ONE job for the ENTIRE grid now. closeness_efficiency.md
                                           # measured ~2 min for n_k=20, S=1000, n=2000 on 16 cores;
                                           # the current config (n_rep=5000, n_obs=4000) is ~10x that.
                                           # If you scale n_k/n_rep/n_obs up further (per your plan
                                           # in config/pitman_closeness.yaml), raise this accordingly.
#SBATCH --cpus-per-task=20                # Real core budget for the WHOLE run (all grid points run
                                           # as local processes inside this one job/node). Must not
                                           # exceed a single node's core count -- check with
                                           # `sinfo -o "%c %m %P"`. Raise together with `n_k` if you
                                           # scale the grid up a lot.
#SBATCH --mem=20G                         # Sized generously for ~20 concurrent 1-core
                                           # `run_grid_point` processes (threads_per_job: 1 in
                                           # config/pitman_closeness.yaml). Tune after checking `seff`.
#SBATCH --output=logs/pitman_closeness/snakemake_slurm-%j.out
#SBATCH --error=logs/pitman_closeness/snakemake_slurm-%j.err

# =============================================================================
# pitman_closeness_slurm.sh -- run pitman_closeness.smk on Aalto Triton as a
# single Slurm job, using Snakemake's plain local `--cores` scheduling (no
# cluster executor plugin, no --profile).
#
# One-time setup (shared with bias_sensitivity_slurm.sh -- same
# environment.yaml env, only built once):
#
#   module load mamba
#   mamba env create --file environment.yaml --prefix env/
#
# Then submit with:
#
#   sbatch pitman_closeness_slurm.sh
#
# Every `run_grid_point` / combine / figures job runs as a plain local
# subprocess on THIS one node, packed by Snakemake's local scheduler to
# respect --cores. `run_grid_point`'s threads come from `threads_per_job`
# in config/pitman_closeness.yaml (currently 1), so with --cores 20 you get
# up to 20 grid points running concurrently, 1 core each. Override the
# per-point core count without editing any file:
#
#   snakemake -s pitman_closeness.smk --cores 20 --set-threads run_grid_point=4
#
# Useful before committing to a run:
#
#   # dry run: how many run_grid_point jobs will the current n_k produce?
#   snakemake -s pitman_closeness.smk --cores 20 -n
#
# This script is safely re-submittable: Snakemake skips any grid point whose
# output already exists, and --rerun-incomplete cleans up + reruns anything
# left partially written by a killed/timed-out job.
# =============================================================================

mkdir -p logs/pitman_closeness

module load mamba
source activate env/

module load triton/2024.1-gcc
module load r

# POT is not on conda-forge/bioconda. Idempotent -- cheap to re-run every
# submission, keeps the shared env self-healing.
Rscript R/install_cran_deps.R

snakemake \
  --snakefile pitman_closeness.smk \
  --cores 20 \
  --keep-going \
  --rerun-incomplete
