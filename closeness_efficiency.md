# Task: Pitman Closeness and Relative efficiency for GPD Shape and scale Estimators

**Goal:** Compare two (scale) shape-parameter estimators for the Generalized Pareto
Distribution (scale fixed at 1) across a grid of true shape values
(-0.5,0.5), using the Pitman closeness criterion and relative efficiency.

**Estimators:** Default to weighted GPD vs. Conventional, but should also allow weighted GPD vs. Cluster maxima (POT).

**Criterions being estimated**, for true shape value $k$  and estimators T1, T2:

  * pitman closeness probabilities: $$ P_k( |T1 - k| < |T2 - k| ); \quad k\in (-0.5,0.5)$$
  i.e. the probability that T1 lands closer to the truth than T2 does, on the
        same sample. T1 is Pitman-closer than T2 at that k if this exceeds 0.5.

  * relative efficiency: $$e(T_1,T_2)=\frac{E[(T_2-k)^2]}{E[(T_1-k)^2]} $$
  
Note:  perform the same operations for the true scale parameter (above the threshold) as well (estimated for different values of k on the grid). 
  

### Steps ( function to compute closeness probabilities)

1. given a grid of true shape values from -0.5 to 0.5 (e.g , 20).
2. For each grid point:
   - Simulate S=100  dependent samples with GPD marginal ("ar1", "log",..)  with that true
     shape and scale 1. Fix and record the sample size and random seed.
   - Fit **both** estimators to the **same** simulated sample each time.
   - Record which estimator's estimate is closer (smaller absolute error)
     to the true shape value; split ties evenly.
   - Compute the win proportion as an average over the S simulated samples:

         P_hat(k) = (1/S) * sum_b  1{ |T1_b - k| < |T2_b - k| } 

     (with ties contributing 0.5). Its Monte Carlo standard error is

         SE = sqrt( P_hat * (1 - P_hat) / S)
	 - compute relative efficiency 
	          $$e(T_1,T_2)=\frac{E[(T_2-k)^2]}{E[(T_1-k)^2]} $$ 
			  where the expectations are obtained empirically. 
	 - perform the same operations for the scale parameter to judge closeness of the estimated scale to the  to the true scale above threshold.
	 - also record the rmse for k and sigma* ( for each grid point) . 
	 
	 
 ### Important notes: 
 
 * The code should be a separate fully optimized/efficient R script and bind to a snake make file to allow paralell processing over the grid of k values 
 
 * source the weighted_gpdfit.R file to use necessary functions for the task
 
 * store the results in a convenient way that will make it easier to make the following plots ( in a separate R script)
 
	 
 #### Plots
 
* For k :
	- Line plot : true shape value (x-axis) vs. win proportion (y-axis);  Add a horizontal reference line at 0.5; Optionally shade a standard-error band; Label which estimator is "winning" above 0.5.
	- line plot : true shape value (x-axis) vs. relative efficiency e (y-axis). Add a horizontal reference line at 1; Label which estimator is "winning" above 1.
	
	- line plot : true shape value (x-axis) vs. rmse ( for both estimators on the same plot)

* for true sigma* ( above threshold.) 
  - same plots as for k, only  k vary on the x-axis though, since true sigma* is fixed. 
  
###  Interpretation

* Pitman Closeness criterion
  - Above 0.5: first estimator closer to truth more often at that shape value.
  - Below 0.5: second estimator wins more often.
  - Consistently above 0.5 across the grid $\rightarrow$ first estimator is uniformly
	Pitman-closer; crossing the line $\rightarrow$ only a local/regional claim holds..
  - This is a sample-by-sample criterion, not average error \u2014 it can
   disagree with which estimator has the lower mean squared error.

* Relative efficiency

 - e(T_1,T_2)  being greater than one would indicate that T_1 is preferable. 


---

## Implementation

### Files

| File | Role |
|---|---|
| `pitman_closeness.smk` | Snakemake workflow (one job per grid point) |
| `config/pitman_closeness.yaml` | all parameters of the experiment |
| `R/pitman_core.R` | Monte-Carlo engine (functions only, runs nothing) |
| `R/run_pitman.R` | CLI: runs **one** grid point → `pc_<i>.rds` |
| `R/combine_pitman.R` | binds the per-`k` files into one tidy table |
| `R/plot_pitman.R` | all figures, from the combined table |
| `R/install_cran_deps.R` | installs `POT` (not available on conda) |
| `run_all.sh` | plain-shell fallback if Snakemake is unavailable |

Everything is namespaced by experiment (`results/pitman_closeness/`,
`logs/pitman_closeness/`) so further experiments can live alongside it.

### Pipeline

```
config/pitman_closeness.yaml
          |
          v
  run_grid_point   x n_k          <- parallel: one job per true shape value
   (R/run_pitman.R)                  each runs S replicates
          |
          v   results/pitman_closeness/raw/pc_1.rds ... pc_<n_k>.rds
      combine
   (R/combine_pitman.R)
          |
          v   results/pitman_closeness/pitman_results.{rds,csv}
      figures
   (R/plot_pitman.R)
          |
          v   results/pitman_closeness/figures/*.pdf
```

**One grid point** = simulate `S` dependent series with GPD margins at that
true `k`; set the threshold `u`; fit **both** estimators to the *same* sample;
record the errors of `k` and of `sigma* = sigma + k*u`.  Both parameters come
out of the same fits, so `sigma*` is free.

The estimator machinery is taken from `weighted_gpdfit.R`.  Note that it is
**not** `source()`d directly: the file ends with a live scratch block
(`sim <- sim_compare(...)`, `x11()`), which would launch a simulation and open
a graphics device in every job.  `source_gpd_functions()` therefore parses it
and evaluates only the top-level function definitions.

### Output

`pitman_results.csv` is long — 2 rows per grid point, `parameter ∈ {k, sigma}`:

```
parameter, k_true, truth, p_hat, se, rel_eff, rmse_t1, rmse_t2,
bias_t1, bias_t2, mean_t1, mean_t2, n_wins, n_ties, n_ok, S, t1, t2, ...
```

which is exactly the shape the plotting script consumes (`p_hat` → closeness
plots, `rel_eff` → efficiency plots, `rmse_t*` → RMSE plots).

**Reproducibility.** Each replicate draws from its own independent
L'Ecuyer-CMRG substream, so results are bit-identical regardless of the number
of cores or the order jobs complete in.  Grid point `i` uses `seed + i - 1`,
so every job is self-contained and resumable.

## Running it

```bash
conda activate snakemake                       # snakemake 9.4.1
snakemake -s pitman_closeness.smk --cores 16   # full pipeline
```

Roughly 2 minutes for `n_k = 20`, `S = 1000`, `n = 2000` on 16 cores.

Useful variants:

```bash
# dry run (show the DAG, execute nothing)
snakemake -s pitman_closeness.smk --cores 16 -n

# override parameters without editing the config
# (config keys must be >= 2 chars: snakemake rejects `S=`, hence `n_rep=`)
snakemake -s pitman_closeness.smk --cores 16 --config n_rep=5000 n_k=40

# stop after the combined table (no figures)
snakemake -s pitman_closeness.smk --cores 16 \
          results/pitman_closeness/pitman_results.csv

# remove this experiment's outputs only
snakemake -s pitman_closeness.smk --cores 1 clean
```

> **Caveat on staleness.**  Changing a parameter (`seed`, `n_rep`, `alpha`, …)
> normally triggers the right re-runs, because Snakemake records the parameters
> of each job in `.snakemake/`.  That tracking is lost, however, if the results
> were **not** produced by Snakemake — copied in, generated by `run_all.sh`, or
> moved between directories.  In that case it reports *"Nothing to be done"*
> and you silently keep the old numbers; the warning to watch for is
> `N jobs have missing provenance/metadata`.  Recover with `--forceall`
> (or `clean`).  Re-plotting alone is always safe.

### Re-making only the plots

The simulation is cached, so tweaking a figure never re-runs it:

```bash
# via the workflow (--forcerun needed: the figures are already up to date)
snakemake -s pitman_closeness.smk --cores 1 --forcerun figures

# or the plotting script alone - no snakemake, no conda
Rscript R/plot_pitman.R \
  --results=results/pitman_closeness/pitman_results.rds \
  --outdir=results/pitman_closeness/figures \
  --device=pdf
```

`plot_pitman.R` also accepts `--device=png`, `--width`, `--height`, `--dpi`.
While iterating on aesthetics, point `--outdir` somewhere scratch so Snakemake's
up-to-date tracking is not disturbed.

### Choosing the estimators

`t1` is the challenger, `t2` the baseline; both are picked from
`"Weighted GPD"`, `"Cluster maxima (POT)"`, `"Conventional"`.  Set them in
`config/pitman_closeness.yaml`, or on the command line:

```bash
snakemake -s pitman_closeness.smk --cores 16 \
          --config t1="Weighted GPD" t2="Cluster maxima (POT)"
```

Send the results elsewhere (or they will overwrite the previous pair's) by
editing `EXPERIMENT` in the `.smk`, or by running `run_all.sh` with different
`T1`/`T2` environment variables.

### Without Snakemake

```bash
Rscript R/install_cran_deps.R     # once: installs POT from CRAN
JOBS=15 ./run_all.sh              # same three stages, parallelised with xargs
```
