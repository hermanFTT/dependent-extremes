#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- plain-shell fallback for the pitman_closeness experiment.
#
# Use this if snakemake is not installed. It runs the same three stages with
# the same commands, parallelising the k-grid with xargs.
#
#   ./run_all.sh            # uses the defaults below
#   JOBS=8 NK=20 S=1000 ./run_all.sh
#
# Prefer the workflow when available (adds caching + resume):
#   snakemake -s pitman_closeness.smk --cores 16
# =============================================================================
set -euo pipefail

EXPERIMENT=pitman_closeness
RES_DIR="results/$EXPERIMENT"
LOG_DIR="logs/$EXPERIMENT"

JOBS=${JOBS:-15}                # grid points to run concurrently
NK=${NK:-20}
S=${S:-1000}
N=${N:-2000}
DEP=${DEP:-markov}
ALPHA=${ALPHA:-0.35}
PHI=${PHI:-0.95}
DEP_MODEL=${DEP_MODEL:-log}
THRESHOLD=${THRESHOLD:-empirical}
T1=${T1:-"Weighted GPD"}
T2=${T2:-"Conventional"}
SEED=${SEED:-4415}
NGRID=${NGRID:-50}
MINGRID=${MINGRID:-100}
DEVICE=${DEVICE:-pdf}           # figure format: pdf or png

mkdir -p "$RES_DIR/raw" "$RES_DIR/figures" "$LOG_DIR"

echo "[run_all] $NK grid points x S=$S, $JOBS concurrent jobs"

seq 1 "$NK" | xargs -P "$JOBS" -I{} sh -c \
  "Rscript R/run_pitman.R --index={} --n_k=$NK --S=$S --n=$N \
     --dep=$DEP --phi=$PHI --alpha=$ALPHA --dep_model=$DEP_MODEL \
     --threshold=$THRESHOLD --n_grid=$NGRID --min_grid_pts=$MINGRID \
     --t1='$T1' --t2='$T2' --seed=$SEED --ncores=1 \
     --out=$RES_DIR/raw/pc_{}.rds > $LOG_DIR/run_pitman_{}.log 2>&1"

echo "[run_all] combining"
Rscript R/combine_pitman.R \
  --out_rds="$RES_DIR/pitman_results.rds" \
  --out_csv="$RES_DIR/pitman_results.csv" \
  "$RES_DIR"/raw/pc_*.rds

echo "[run_all] plotting"
Rscript R/plot_pitman.R --results="$RES_DIR/pitman_results.rds" \
  --outdir="$RES_DIR/figures" --device="$DEVICE"

echo "[run_all] done -> $RES_DIR/figures"
