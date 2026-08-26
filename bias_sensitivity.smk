# =============================================================================
# bias_sensitivity.smk
#
# Experiment: finite-sample bias, RMSE and coverage of the GPD shape, scale
#             and return level under temporal dependence, for three estimation
#             approaches.
#
# The unit of work is ONE CELL of the grid
#
#       dep_par  x  k_grid  x  n_obs  x  run_len
#
# and all cells are independent, so the whole product runs in parallel.
#
# This workflow is NOT named `Snakefile`; select it with `-s`:
#
#   snakemake -s bias_sensitivity.smk --cores 16            # run everything
#   snakemake -s bias_sensitivity.smk --cores 16 -n         # dry run
#   snakemake -s bias_sensitivity.smk --cores 1 --forcerun figures
#   snakemake -s bias_sensitivity.smk --cores 16 --config n_rep=200
#   snakemake -s bias_sensitivity.smk --cores 1 clean
#
# Parameters live in config/bias_sensitivity.yaml.
# =============================================================================

import itertools

EXPERIMENT = "bias_sensitivity"

configfile: "config/" + EXPERIMENT + ".yaml"

RES_DIR = "results/" + EXPERIMENT
LOG_DIR = "logs/" + EXPERIMENT
RAW_DIR = RES_DIR + "/raw"
FIG_DIR = RES_DIR + "/figures"


def as_list(x):
    """Accept either a single value or a list for any varied parameter."""
    return list(x) if isinstance(x, (list, tuple)) else [x]


def tok(x):
    """Filename-safe token for a parameter value (no spaces, keeps sign)."""
    s = str(x)
    return s.replace(" ", "").replace("/", "-")


DEP_TYPE = config["dep_type"]
DEP_PARS = as_list(config["dep_par"])
K_GRID   = as_list(config["k_grid"])
N_OBS    = as_list(config["n_obs"])
RUN_LEN  = as_list(config["run_len"])

# --- expand the cartesian product into one named cell per combination -------
# The cell name is self-describing, so results are readable on disk and adding
# a new value never renumbers the existing cells.
CELLS = {}
for dp, k, n, r in itertools.product(DEP_PARS, K_GRID, N_OBS, RUN_LEN):
    name = "dp{}_k{}_n{}_r{}".format(tok(dp), tok(k), tok(n), tok(r))
    CELLS[name] = dict(dep_par=dp, k=k, n=n, r=r)

CELL_NAMES = sorted(CELLS)


def cell_flags(wildcards):
    c = CELLS[wildcards.cell]
    return (
        f"--dep_type={DEP_TYPE} --dep_par={c['dep_par']} "
        f"--k={c['k']} --n={c['n']} --r={c['r']} "
        f"--n_rep={config['n_rep']} --sigma={config['sigma']} "
        f"--p_thresh={config['p_thresh']} --threshold={config['threshold']} "
        f"--m_ret={config['m_ret']} "
        f"--coverage={str(config['coverage']).upper()} "
        f"--cred_mass={config['cred_mass']} --interval={config['interval']} "
        f"--n_draws={config['n_draws']} "
        f"--n_grid={config['n_grid']} --min_grid_pts={config['min_grid_pts']} "
        f"--wip={str(config['wip']).upper()} "
        f'--approaches="{config["approaches"]}" '
        f"--seed={config['seed']}"
    )


wildcard_constraints:
    cell=r"[^/]+",


rule all:
    input:
        RES_DIR + "/bias_results.csv",
        RES_DIR + "/.figures.done",


# --- one job per grid cell --------------------------------------------------
rule run_cell:
    output:
        RAW_DIR + "/{cell}.rds",
    log:
        LOG_DIR + "/run_{cell}.log",
    threads: int(config.get("threads_per_job", 1))
    params:
        flags=cell_flags,
    shell:
        "mkdir -p " + LOG_DIR + " && "
        "Rscript R/run_bias.R {params.flags} "
        "--ncores={threads} --out={output} > {log} 2>&1"


# --- bind every cell into one tidy long table -------------------------------
rule combine:
    input:
        expand(RAW_DIR + "/{cell}.rds", cell=CELL_NAMES),
    output:
        rds=RES_DIR + "/bias_results.rds",
        csv=RES_DIR + "/bias_results.csv",
    log:
        LOG_DIR + "/combine.log",
    shell:
        "mkdir -p " + LOG_DIR + " && "
        "Rscript R/combine_bias.R "
        "--out_rds={output.rds} --out_csv={output.csv} {input} > {log} 2>&1"


# --- figures ----------------------------------------------------------------
# How many figures are produced depends on which dimensions actually vary, so
# the outputs are represented by a sentinel file rather than a fixed list.
rule figures:
    input:
        RES_DIR + "/bias_results.rds",
    output:
        touch(RES_DIR + "/.figures.done"),
    log:
        LOG_DIR + "/plot.log",
    params:
        device=config["fig_device"],
        dpi=config["fig_dpi"],
        approaches=config["approaches"],
        fix_k=config.get("fix_k", ""),
        fix_dep=config.get("fix_dep", ""),
        fix_r=config.get("fix_r", ""),
    shell:
        "mkdir -p " + LOG_DIR + " && "
        "Rscript R/plot_bias.R --results={input} --outdir=" + FIG_DIR + " "
        "--device={params.device} --dpi={params.dpi} "
        '--approaches="{params.approaches}" '
        '--fix_k="{params.fix_k}" --fix_dep="{params.fix_dep}" '
        '--fix_r="{params.fix_r}" > {log} 2>&1'


# --- remove ONLY this experiment's outputs ----------------------------------
rule clean:
    shell:
        "rm -rf " + RES_DIR + " " + LOG_DIR
