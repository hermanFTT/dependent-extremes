#!/bin/bash
#SBATCH --job-name=pitman_closeness
#SBATCH --time=04:00:00                   # Driver only submits/monitors jobs. pitman_closeness's
                                           # grid (n_k points) is far smaller/cheaper than
                                           # bias_sensitivity's, hence starting well below that
                                           # script's 2 days. Raise this if you scale n_k, n_rep or
                                           # n_obs up a lot in config/pitman_closeness.yaml and this
                                           # stops being enough (check with `sacct -j <jobid>`).
#SBATCH --cpus-per-task=2                 # Only for the snakemake driver process itself
#SBATCH --mem=2G                          # Only for the snakemake driver process itself
#SBATCH --output=logs/pitman_closeness/snakemake_slurm-%j.out
#SBATCH --error=logs/pitman_closeness/snakemake_slurm-%j.err

# =============================================================================
# pitman_closeness_slurm.sh -- submit pitman_closeness.smk to Aalto Triton.
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
# `--profile profiles/slurm_pitman/` (see that file) is a separate profile
# from bias_sensitivity's `profiles/slurm/`, sized for pitman_closeness's
# much smaller/cheaper grid (`run_grid_point`, one job per true-shape value).
#
#   snakemake -s pitman_closeness.smk --profile profiles/slurm_pitman/ -n
#
# is a useful dry run first -- it reports exactly how many `run_grid_point`
# jobs the current config/pitman_closeness.yaml grid (`n_k`) will submit,
# handy to check again after you scale the grid up.
# =============================================================================

mkdir -p logs/pitman_closeness

module load mamba
source activate env/

module load triton/2024.1-gcc
module load r


# NB on --cores: even for jobs submitted to Slurm, Snakemake clamps each
# job's threads to min(threads/set-threads, --cores) -- so this must stay
# >= the highest set-threads value used by any rule in
# profiles/slurm_pitman/config.yaml (currently 1 for every rule, hence 2 is
# fine here, but bump this together with set-threads if you ever raise
# run_grid_point above 2 -- see bias_sensitivity_slurm.sh, which hit exactly
# this with run_cell: 8 vs. an unraised --cores 2).
snakemake \
  --snakefile pitman_closeness.smk \
  --profile profiles/slurm_pitman/ \
  --cores 2 \
  --rerun-incomplete
