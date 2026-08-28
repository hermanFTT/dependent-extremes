# Running on Triton (Aalto HPC)

Everything below is run from the **repo root** — all script paths in the
workflows (`R/bias_core.R`, `weighted_gpdfit.R`) are relative to it.

## 1. Put the repo on Triton

`$HOME` is 10 GB and intended for code and configuration, not calculation
files. Work from `$WRKDIR` (large, fast, per-user Lustre — **not backed up**).

```bash
ssh triton.aalto.fi                      # requires Aalto VPN
cd $WRKDIR
git clone <your-remote> dependent-extremes
cd dependent-extremes
```

Or push from your laptop:

```bash
rsync -av --exclude .snakemake --exclude results --exclude .git \
  "~/Desktop/Git Repositories/dependent-extremes/" \
  triton.aalto.fi:'$WRKDIR/dependent-extremes/'
```

## 2. Bootstrap (once)

```bash
bash triton/bootstrap.sh
```

This creates a conda env with `snakemake=9.4.1`, `pulp=2.8` and
`snakemake-executor-plugin-slurm`, installs the five R packages the pipeline
needs into a private library, verifies they load, and smoke-tests one tiny
cell.

Why the plugin: Snakemake dropped `--cluster "sbatch ..."` in v8; distributed
execution now goes through executor plugins.

## 3. What to re-run, and when

| | Frequency | Why |
|---|---|---|
| `bash triton/bootstrap.sh` | **once, ever** | The conda env and the R library are written to disk and persist across logins. Re-run only after adding a package to `install_cran_deps.R`, or if Triton bumps R to a new major.minor (the library path is version-keyed). Idempotent. |
| `source triton/env.sh` | every new shell, **interactive work only** | Lmod modules and `PATH` / `R_LIBS_USER` live in the shell process and die on logout. |
| `sbatch triton/run_*.sbatch` | every run | Sources `env.sh` itself on the compute node — needs no manual setup. |

So a normal session is just:

```bash
ssh triton.aalto.fi
cd $WRKDIR/dependent-extremes
sbatch triton/run_bias.sbatch
```

You only need `source triton/env.sh` by hand for interactive things: a dry
run, a one-off `Rscript`, or the `--profile` route in tmux. Note it must be
**sourced, not executed** — Lmod modules are shell functions.

## 4. How big is this job, really?

Measured on one core at the current config (`n_rep: 5000`, `n_obs: 4000`,
`coverage: true`, `n_draws: 10000`, all three approaches):

| | per unit | units | total CPU | wall on 16 cores |
|---|---|---|---|---|
| `bias_sensitivity` | ~4 min, ~120 MB | 76 cells | ~5 CPU-h | **~20 min** |
| `pitman_closeness` | ~3.5 min, ~120 MB | 20 points | ~1.2 CPU-h | **~5 min** |

Timing is uniform across the grid — `dep_par` ∈ {0.1, 0.5, 1} × `k` ∈ {−0.3,
0.3} all landed within 1.2–1.3 s at `S=20`. Memory is flat in `n_rep`
(500 reps peaked at the same 115 MB as 20), because the loop accumulates
summaries, not draws.

**This is a small job.** Prefer the single-job route below; the Slurm executor
only pays off once you scale up by ~10x.

## 5. Recommended: one Slurm job, Snakemake parallelising inside it

```bash
sbatch triton/run_bias.sbatch        # ~20 min on 16 cores
sbatch triton/run_pitman.sbatch      # ~5 min
```

That's the whole thing. No executor plugin, no login-node babysitting, no
tmux, and it survives your laptop closing. Watch it with:

```bash
squeue -u $USER
tail -f logs/slurm-bias-<jobid>.out
```

Dry run first if you want to see the DAG:

```bash
source triton/env.sh
snakemake -s bias_sensitivity.smk --cores 16 -n
```

## 6. Alternative: one Slurm job per cell (for scaled-up runs)

Worth it when a single cell approaches an hour — e.g. `n_rep: 50000`, or a
much denser `dep_par`/`k_grid`/`n_obs` product. Then 76+ independent jobs
backfill across the cluster instead of queueing behind one 16-core allocation.

Keep Snakemake on the **login node** inside `tmux`: running it inside a Slurm
job leads to unpredictable behaviour, and the plugin detects this and warns.

```bash
tmux new -s bias
cd $WRKDIR/dependent-extremes
source triton/env.sh

snakemake -s bias_sensitivity.smk --profile profiles/triton-bias -n   # dry run
snakemake -s bias_sensitivity.smk --profile profiles/triton-bias      # go
# detach: Ctrl-b then d     reattach: tmux attach -t bias

snakemake -s pitman_closeness.smk --profile profiles/triton-pitman
```

If you scale `n_rep` up, scale `run_cell.runtime` in
`profiles/triton-bias/config.yaml` proportionally — it is currently 30 min,
sized for the ~4 min measured above.

## 7. Monitor

```bash
squeue -u $USER                                  # queued / running jobs
tail -f logs/bias_sensitivity/run_dp0.5_k0.3_n4000_rauto.log   # per-cell log
seff <jobid>                                     # actual mem/CPU used, post-hoc
```

`seff` after the first run tells you whether the `--mem` / `--time` in the
sbatch scripts can be tightened.

## Notes / gotchas

- **`clean` must run locally.** It has a shell command, so the profile would
  submit it to Slurm. Override:
  `snakemake -s bias_sensitivity.smk --executor local --cores 1 clean`
- **`detectCores()`** in `bias_core.R:176` / `pitman_core.R:148` reports the
  whole node, not the Slurm allocation. Harmless here — both `.smk` files
  always pass `--ncores={threads}` — but do not call those functions manually
  without `ncores`.
- **`OMP_NUM_THREADS=1`** is exported by `env.sh` so the BLAS does not
  oversubscribe on top of `parallel::mclapply`.
- **ggplot2 version drift.** `environment.yaml` pins 3.5.2; your laptop has
  4.0.3; `scicomp-r-env` may ship a third. Figures may differ from local
  output. Pin explicitly into `R_LIBS_USER` if that matters.
- **`fig_device: "pdf"`** — `plot_bias.R:98` falls back to `NULL` when
  `cairo_pdf` is unavailable, so headless nodes will not crash, just render
  differently.
- **Not needed on the cluster:** `patchwork`, `tidyr`, `latex2exp`, `glue`,
  `devtools`, `posterior`. Those belong to `sbc_exp.R` and
  `R/dependence_aware_unif_tests.R`, which neither workflow invokes. This also
  avoids `sbc_exp.R:13`'s `devtools::load_all("../posterior")`, whose path does
  not exist on Triton.
