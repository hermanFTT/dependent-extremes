# Task: Bias and Sensitivity Analysis

## Goal

Study the finite-sample bias, RMSE, and (optionally) coverage probability of
the shape parameter, scale parameter, and a return level, under different
degrees of temporal dependence, for three estimation approaches:

- Conventional POT (peaks-over-threshold)
- Weighted GPD
- Cluster maxima POT

## Dependence structures

Choose one dependence type per experiment run. Each has its own parameter
grid, and the grid endpoints correspond to different dependence regimes:

| Dependence type | Parameter | Grid | Regime at low end | Regime at high end |
|---|---|---|---|---|
| AR(1)-type | \u03c6 | (0.1, 0.98) | approx. independence | total dependence |
| Higher-order Markov ("log") | \u03b1 | (0.1, 0.98) | total dependence | approx. independence |
| Higher-order Markov ("nlog") | \u03b1 | (0.1, 30) | total independence | some dependence |

## Experimental design

**Fixed:**
- Scale: \u03c3 = 1
- Threshold: u = 95th percentile of the marginal distribution

**Varied:**
- Shape: k, grid over (-0.5, 0.5)
- Sample size: n = seq(1000, 20000, by = 1000)
- Run length r (used by the Weighted GPD and cluster-maxima POT approaches):
  either `NULL` (automatic selection) or a fixed value from
  {5, 7, 10, 13, 16, 20}

**Simulation:**
- Number of Monte Carlo replicates: S (e.g. S = 1000)
- Each replicate is a dependent series with GPD(\u03c3, k) margins, generated
  under the chosen dependence type and parameter value
- Return period: m (used to compute the return level z_m)

**Optional:**
- Coverage probability computation for k and z_m, with a user-specified
  nominal credible/confidence level, toggled on or off

## Estimation and metrics

For each parameter combination (dependence type, dependence parameter, k, n,
r, and any other varying setting), fit all three approaches to each of the S
replicates and compute, for each of k, \u03c3, and z_m:

- **Mean estimate**: average of the S fitted values
- **Bias**: mean estimate minus the true value,
  `bias = mean(estimate) - true_value`
- **RMSE**: `sqrt( mean( (estimate - true_value)^2 ) )`
- **Coverage probability** (if requested, for k and z_m only): the
  proportion of the S replicates whose interval at the specified nominal
  level contains the true value

## Implementation notes

- Write the simulation/estimation code as a single, optimized R script,
  separate from the plotting code, and drive it with a Snakemake workflow so
  that different parameter combinations run in parallel.
- The config file must allow each parameter to be given as either a single
  value or a set/grid, e.g. `n: [4000]` or `n: [1000, 4000]`, rather than
  forcing a full `seq()` range every time.
- Store results in a tidy, long-format table (one row per parameter
  combination × approach × target parameter) with columns for the
  dependence type, dependence parameter value, k, n, r, approach, target
  parameter (k / \u03c3 / z_m), mean estimate, bias, RMSE, and coverage
  probability (if computed). This makes all downstream plots simple filters
  and group-bys rather than custom reshaping.
  

## Plots (separate R script)

All plots overlay the three approaches on the same axes, with an option to
select a subset of approaches to compare. Unless noted, "observation size"
means the number of exceedances above the threshold, n·(1 \u2212 u).

| # | X-axis | Y-axis | Held fixed | Reference line |
|---|---|---|---|---|
| 1 | Dependence level (\u03c6 or \u03b1) | Bias (for chosen k, \u03c3, z_m) | Observation size | 0 |
| 2 | Observation size | Bias (for chosen k, \u03c3, z_m) | Dependence level | 0 |
| 3 | Dependence level | RMSE (for chosen k, \u03c3, z_m) | Observation size | \u2014 |
| 4 | Observation size | RMSE (for chosen k, \u03c3, z_m) | Dependence level | \u2014 |
| 5 | Dependence level | Coverage probability (for chosen k, z_m) | Observation size | nominal level |
| 6 | Observation size | Coverage probability (for chosen k, z_m) | Dependence level | nominal level |
| 7 | True k | Coverage probability | Observation size, dependence level | nominal level |

---

## Implementation

### Files

| File | Role |
|---|---|
| `bias_sensitivity.smk` | Snakemake workflow (one job per grid cell) |
| `config/bias_sensitivity.yaml` | all parameters; every varied one is a list |
| `R/bias_core.R` | simulation/estimation engine (functions only) |
| `R/run_bias.R` | CLI: runs **one** cell → `<cell>.rds` |
| `R/combine_bias.R` | binds all cells into one tidy long table |
| `R/plot_bias.R` | all figures, from the combined table |

Outputs are namespaced under `results/bias_sensitivity/` and
`logs/bias_sensitivity/`, alongside the other experiments.

### Pipeline

```
config/bias_sensitivity.yaml
          |
          v
      run_cell  x (dep_par x k x n x r)     <- parallel: one job per cell
     (R/run_bias.R)                            each runs S replicates
          |
          v   results/bias_sensitivity/raw/dp0.5_k0.2_n4000_rauto.rds ...
      combine
   (R/combine_bias.R)
          |
          v   results/bias_sensitivity/bias_results.{rds,csv}
      figures
    (R/plot_bias.R)
          |
          v   results/bias_sensitivity/figures/*.pdf
```

**One cell** = `(dep_type, dep_par, k, n, r)`.  Within a cell, `S` series are
simulated; **all three approaches are fitted to the same sample**, so the
comparison is paired.  A replicate is discarded for *every* approach if any
one of them fails, which keeps the pairing intact.

Cell file names are self-describing (`dp0.5_k0.2_n4000_rauto.rds`), so adding
a new `n` or `k` never renumbers existing cells and only the new ones are run.

The estimators come from `weighted_gpdfit.R`.  It is **not** `source()`d
directly — the file ends with a live scratch block (`sim <- sim_compare(...)`,
`x11()`) which would run a simulation and open a graphics device in every job;
`source_gpd_functions()` evaluates only its top-level function definitions.

### Targets and truths

| Target | Truth |
|---|---|
| `k` | the simulated shape |
| `sigma` | `sigma* = sigma + k*u`, the scale **above the threshold**, recorded per replicate (with an empirical threshold `u` varies between samples) |
| `zm` | `z_m = F^-1(1 - 1/m)`, from `m_ret` |

Coverage is computed for `k` and `z_m` only — `sigma*` has no posterior grid in
this implementation, so its `coverage` column is `NA` by construction.

### Output

`bias_results.csv` is long: one row per cell × approach × target.

```
approach, parameter, truth, mean_est, bias, rmse, mae, sd_est,
coverage, ci_width, n_ok, dep_type, dep_par, k_true, n, r, r_label,
n_exc, n_exc_nom, S, m, sigma, p_thresh, threshold, cred_mass,
interval, seed, secs
```

`n_exc_nom = n(1-p_thresh)` is the nominal observation size used on the plot
axes; `n_exc` is the realised mean number of exceedances.

**Reproducibility.** Each replicate draws from its own L'Ecuyer-CMRG substream,
so results are identical regardless of core count or scheduling.  Each cell
derives its own seed from `seed` plus a hash of the cell identity, so cells are
independent and unaffected by the shape of the rest of the grid.

## Running it

```bash
conda activate snakemake
snakemake -s bias_sensitivity.smk --cores 16
```

The default grid is 5 dep × 5 k × 5 n × 1 r = **125 cells**; at `n_rep: 1000`
this was measured at **33 min wall-clock on 16 cores** (~7 CPU-hours).  Cost is
dominated by the largest `n` (an `n = 20000` cell with `S = 1000` takes ~3 min
on one core, versus ~17 s at `n = 1000`).  Coverage roughly doubles the cost;
set `coverage: false` for a quick bias/RMSE-only pass.

Useful variants:

```bash
# dry run: how many cells will this grid produce?
snakemake -s bias_sensitivity.smk --cores 16 -n

# a small exploratory grid, without editing the config
snakemake -s bias_sensitivity.smk --cores 16 \
  --config dep_par="[0.1,0.5,0.98]" k_grid="[-0.2,0.2]" \
           n_obs="[1000,4000]" n_rep=100

# compare run lengths instead of using automatic selection
snakemake -s bias_sensitivity.smk --cores 16 --config run_len='["auto",5,10,20]'

# switch dependence structure (note nlog uses a (0.1, 30) grid)
snakemake -s bias_sensitivity.smk --cores 16 \
  --config dep_type=nlog dep_par="[0.1,1,5,15,30]"

# stop after the table; remove this experiment's outputs
snakemake -s bias_sensitivity.smk --cores 16 results/bias_sensitivity/bias_results.csv
snakemake -s bias_sensitivity.smk --cores 1 clean
```

Config keys must be ≥ 2 characters — snakemake rejects `--config S=…`, which is
why the keys are `n_rep`, `n_obs`, `m_ret` rather than `S`, `n`, `m`.

> **`--config` overrides are not remembered.**  They define the grid for *that
> invocation only*.  If you run a small grid with overrides and then issue any
> later command without them — even `--forcerun figures` — Snakemake targets the
> **full default grid** from the config file and will happily start computing
> all of its missing cells.  Either repeat the same `--config` flags on every
> command, or edit `config/bias_sensitivity.yaml` for anything you intend to
> run more than once.

### Re-making only the plots

```bash
# via the workflow
snakemake -s bias_sensitivity.smk --cores 1 --forcerun figures

# or directly - no snakemake, no conda
Rscript R/plot_bias.R \
  --results=results/bias_sensitivity/bias_results.rds \
  --outdir=results/bias_sensitivity/figures
```

The grid has more dimensions than a plot has axes.  Anything that is neither
on an axis nor in the legend is **faceted** (the "held fixed" column of the
table above: one panel per value) or **filtered** to a single value.  Filters
default to the middle value present in the data, and are always printed in the
subtitle so no figure silently hides a dimension.  Override them, and select a
subset of approaches, with:

```bash
Rscript R/plot_bias.R --results=... --outdir=... \
  --fix_k=0.2 --fix_dep=0.5 --fix_r=auto \
  --approaches="Weighted GPD,Conventional"
```

Figures are named `fig<N>_<metric>_vs_<x>_<target>.pdf`, matching the table
above (e.g. `fig5_cov_vs_dep_k.pdf`).  Only the applicable ones are written:
`sigma` has no coverage, so figures 5–7 are produced for `k` and `zm` only.
If more than one run length is simulated, an extra `fig8_rmse_by_runlength_k`
compares them.

> **Caveat on staleness.** As with the other experiments, changing a parameter
> re-runs the affected cells only if Snakemake's provenance metadata in
> `.snakemake/` is intact.  If results were copied in or moved, it reports
> *"Nothing to be done"* and silently keeps the old numbers — watch for the
> `N jobs have missing provenance/metadata` warning and use `--forceall`.
