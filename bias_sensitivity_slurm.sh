#!/bin/bash
#SBATCH --job-name=bias_sensitivity
#SBATCH --time=1-00:00:00                 # UNVERIFIED starting guess -- benchmark first! This is
                                           # partition batch`). Benchmark a small subset first (see
                                           # below) before trusting this number.
#SBATCH --cpus-per-task=200                # Real core budget for the WHOLE run now (all cells run
                                           # as local processes inside this one job/node -- there is
                                           # no cluster-wide distribution anymore). Must not exceed a
                                           # single node's core count or the job will sit PENDING
                                           # forever -- check with `sinfo -o "%c %m %P"`, then set
                                           # this to whatever a `batch` node actually has (or less).
#SBATCH --mem=30G                         # Sized generously for ~20 concurrent 1-core `run_cell`
                                           # processes (threads default to 1 -- see
                                           # bias_sensitivity.smk). Tune down after checking `seff`
                                           # on a completed run; real usage per cell has been
                                           # observed well under 1 GB.
#SBATCH --output=logs/bias_sensitivity/snakemake_slurm-%j.out
#SBATCH --error=logs/bias_sensitivity/snakemake_slurm-%j.err

# =============================================================================
# bias_sensitivity_slurm.sh -- run bias_sensitivity.smk on Aalto Triton as a
# single Slurm job, using Snakemake's plain local `--cores` scheduling (no
# cluster executor plugin, no --profile).
#
# One-time setup (from the dependent-extremes/ project root):
#
#   module load mamba
#   mamba env create --file environment.yaml --prefix env/
#
# Then submit with:
#
#   sbatch bias_sensitivity_slurm.sh
#
# How this differs from a cluster-executor setup: every `run_cell` / combine
# / figures job runs as a plain local subprocess on THIS one node, packed by
# Snakemake's local scheduler to respect --cores. `run_cell`'s threads
# default to 1 (see bias_sensitivity.smk), so with --cores 20 you get up to
# 20 cells running concurrently, 1 core each -- as each finishes, the next
# pending cell takes its slot. Override the per-cell core count (fewer,
# heavier cells instead) WITHOUT editing any file:
#
#   snakemake -s bias_sensitivity.smk --cores 20 --set-threads run_cell=4
#
# Useful before committing to a long run:
#
#   # dry run: how many cells will this grid produce?
#   snakemake -s bias_sensitivity.smk --cores 20 -n
#
#   # benchmark ONE cheap cell for real timing/memory numbers, then scale up
#   # --cpus-per-task/--mem/--time above accordingly
#   snakemake -s bias_sensitivity.smk --cores 1 \
#     --config dep_par="[0.5]" k_grid="[0.2]" n_obs="[1000]" n_rep=1000 \
#     results/bias_sensitivity/raw/dp0.5_k0.2_n1000_rauto_m2000.rds
#
# This script is safely re-submittable: Snakemake skips any cell whose
# output already exists, and --rerun-incomplete cleans up + reruns anything
# left partially written by a killed/timed-out job.
# =============================================================================

mkdir -p logs/bias_sensitivity

module load mamba
source activate env/

module load triton/2024.1-gcc
module load r

# unlock, just in case the program didn't finished on previous run
#snakemake -s bias_sensitivity.smk --cores 1 --unlock

snakemake \
  --snakefile bias_sensitivity.smk \
  --until run_cell\
  --cores 200 \
  --set-threads run_cell=4 \
  --keep-going \
 --rerun-incomplete  
