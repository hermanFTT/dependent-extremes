# =============================================================================
# pitman_closeness.smk
#
# Experiment: Pitman closeness & relative efficiency for GPD shape/scale
#             estimators under dependence.
#
# The k-grid is embarrassingly parallel: one job per true shape value.
#
# This workflow is NOT named `Snakefile`, so it must always be selected with
# `-s`.  That keeps several experiments side by side in one repository; each
# owns its own .smk file, its own config/<name>.yaml, and its own
# results/<name>/ + logs/<name>/ trees, so nothing can collide.
#
#   snakemake -s pitman_closeness.smk --cores 16            # run everything
#   snakemake -s pitman_closeness.smk --cores 16 -n         # dry run
#   snakemake -s pitman_closeness.smk --cores 1 --forcerun figures
#   snakemake -s pitman_closeness.smk --cores 16 --config S=5000 n_k=40
#   snakemake -s pitman_closeness.smk --dag | dot -Tpng > dag.png
#
# Parameters live in config/pitman_closeness.yaml.
# =============================================================================

EXPERIMENT = "pitman_closeness"

configfile: "config/" + EXPERIMENT + ".yaml"

NK      = int(config["n_k"])
INDICES = list(range(1, NK + 1))

RES_DIR = "results/" + EXPERIMENT
LOG_DIR = "logs/" + EXPERIMENT
RAW_DIR = RES_DIR + "/raw"
FIG_DIR = RES_DIR + "/figures"

FIGURES = expand(
    FIG_DIR + "/{fig}." + config["fig_device"],
    fig=["closeness_k",     "releff_k",     "rmse_k",
         "closeness_sigma", "releff_sigma", "rmse_sigma",
         "closeness_both",  "releff_both",  "rmse_both"],
)


# --- helper: turn config into the --key=value flags run_pitman.R expects -----
def sim_flags():
    c = config
    run = c.get("run", None)
    flags = [
        f"--n_k={c['n_k']}",
        f"--k_min={c['k_min']}",
        f"--k_max={c['k_max']}",
        f"--S={c['n_rep']}",
        f"--n={c['n_obs']}",
        f"--dep={c['dep']}",
        f"--phi={c['phi']}",
        f"--alpha={c['alpha']}",
        f"--dep_model={c['dep_model']}",
        f"--sigma={c['sigma']}",
        f"--p_thresh={c['p_thresh']}",
        f"--threshold={c['threshold']}",
        f"--n_grid={c['n_grid']}",
        f"--min_grid_pts={c['min_grid_pts']}",
        f"--wip={str(c['wip']).upper()}",
        f'--t1="{c["t1"]}"',
        f'--t2="{c["t2"]}"',
        f"--seed={c['seed']}",
    ]
    if run is not None:
        flags.append(f"--run={run}")
    return " ".join(flags)


rule all:
    input:
        FIGURES,
        RES_DIR + "/pitman_results.csv",


# --- one job per true shape value -------------------------------------------
rule run_grid_point:
    output:
        RAW_DIR + "/pc_{index}.rds",
    log:
        LOG_DIR + "/run_pitman_{index}.log",
    threads: int(config.get("threads_per_job", 1))
    params:
        flags=sim_flags(),
    shell:
        "mkdir -p " + LOG_DIR + " && "
        "Rscript R/run_pitman.R "
        "--index={wildcards.index} {params.flags} "
        "--ncores={threads} --out={output} > {log} 2>&1"


# --- bind the per-k results into one tidy table ------------------------------
rule combine:
    input:
        expand(RAW_DIR + "/pc_{index}.rds", index=INDICES),
    output:
        rds=RES_DIR + "/pitman_results.rds",
        csv=RES_DIR + "/pitman_results.csv",
    log:
        LOG_DIR + "/combine.log",
    shell:
        "mkdir -p " + LOG_DIR + " && "
        "Rscript R/combine_pitman.R "
        "--out_rds={output.rds} --out_csv={output.csv} {input} > {log} 2>&1"


# --- figures -----------------------------------------------------------------
rule figures:
    input:
        RES_DIR + "/pitman_results.rds",
    output:
        FIGURES,
    log:
        LOG_DIR + "/plot.log",
    params:
        device=config["fig_device"],
        dpi=config["fig_dpi"],
    shell:
        "mkdir -p " + LOG_DIR + " && "
        "Rscript R/plot_pitman.R --results={input} --outdir=" + FIG_DIR + " "
        "--device={params.device} --dpi={params.dpi} > {log} 2>&1"


# --- remove ONLY this experiment's outputs ----------------------------------
rule clean:
    shell:
        "rm -rf " + RES_DIR + " " + LOG_DIR
