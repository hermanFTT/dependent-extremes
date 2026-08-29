#!/bin/bash
#SBATCH --job-name=bias_sensitivity
#SBATCH --time=2-00:00:00                 # Driver only submits/monitors jobs, but with
                                           # thousands of `run_cell` cells and `jobs: 200`
                                           # in profiles/slurm/config.yaml, the whole
                                           # workflow can take a long time end-to-end.
                                           # Raise this if snakemake gets killed before the
                                           # grid finishes (check `sacct -j <jobid>`), or
                                           # shrink the grid in config/bias_sensitivity.yaml.
#SBATCH --cpus-per-task=2                 # Only for the snakemake driver process itself
#SBATCH --mem=2G                          # Only for the snakemake driver process itself
#SBATCH --output=logs/bias_sensitivity/snakemake_slurm-%j.out
#SBATCH --error=logs/bias_sensitivity/snakemake_slurm-%j.err

# =============================================================================
# bias_sensitivity_slurm.sh -- submit bias_sensitivity.smk to Aalto Triton.
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
# `--profile profiles/slurm/` (see that file) tells snakemake to submit each
# rule (run_cell / combine / figures) as its own Slurm job via the slurm
# executor plugin, instead of running the whole grid inside this one job.
#
#   snakemake -s bias_sensitivity.smk --profile profiles/slurm/ -n
#
# is a useful dry run first -- it reports exactly how many `run_cell` jobs
# the current config/bias_sensitivity.yaml grid will submit.
# =============================================================================

mkdir -p logs/bias_sensitivity

module load mamba
source activate env/

module load triton/2024.1-gcc
module load r

# POT (needed by weighted_gpdfit.R / the "markov" dependence simulator) is
# not on conda-forge/bioconda, so it isn't in environment.yaml. This script
# is idempotent -- it only installs what's missing -- so re-running it here
# on every submission is cheap and keeps the shared env self-healing.

snakemake \
  --snakefile bias_sensitivity.smk \
  --profile profiles/slurm/ \
  --cores 2 \
  --rerun-incomplete
