# Simulations (Supplementary Figs. S1–S6; Supplementary Tables S1, S2)

The simulation studies were run once, on a cluster, with the G-MoQ packages at the
Stage 1 release named in `config.yaml` (`package.stage1`). The code is in `code/`
(studies 1–6) and `lowsep_se_check/` (study 7). Their per-cell summaries are kept in
`inputs/` and imported into `results/08_simulations/summary/` by

```
Rscript analysis/08_simulations/import_simulation_summaries.R
```

which copies them unchanged, adds the mixture timing rows and recomputes the summary
sentence of study 7. The supplementary figures and tables are then drawn from `results/` by
`figures/figS01_*.R` … `figS06_*.R`, `tables/tableS1_*.R` and `tables/tableS2_*.R`.
Study 8 (the weight-design simulation) is in `analysis/02_mixtures/`.

## Studies

Scripts are under `code/validation/` unless a path is given.

| Study | Question | Script(s) | Summary file in `inputs/` | Paper |
|---|---|---|---|---|
| 1 | Does Stage 1 reproduce ordinary RIF–OLS, coefficient by coefficient? | `figures/run_stage1.py`, with `figures/stage1_methods.py`, `figures/export_stage1.py`, `benchmark.py`; job `figures/stage1.sbatch` | `simulation_summaries/TableS2_numerical_summary.csv` (S1 rows) | Supplementary Fig. S1; Table S2 |
| 2 | How fast is Stage 1, and how much memory does it use? | as study 1 (same run) | `simulation_summaries/S2_mean_timings.csv`; Table S2 (S2 rows) | Supplementary Fig. S2; Table S2 |
| 3 | Does Stage 2 recover simple simulated effects and their projection targets? | `figures/simulate_stage2.py`; job `figures/stage2.sbatch` | Table S2 (S3 rows) | Supplementary Fig. S3; Table S2 |
| 4 | Mixture recovery against direct mixture maximum likelihood | `figures/mixture_comparison.sh`, with `comparators.py` | Table S2 (S4 rows) | Supplementary Fig. S4a,b,d; Table S2 |
| 5 | Cost of G-MoQ against mixture maximum likelihood on matched panels | as study 4 (same run, panel arm) | `simulation_summaries/S5_mean_timings.csv`; Table S2 (S5 rows) | Supplementary Fig. S5a–c; Table S2 |
| 6 | Pooled against jackknife standard errors, common and rare variants | `figures/rare_variant_se.sh` | Table S2 (S6 rows) | Supplementary Fig. S6; Table S2 |
| 1–6 | Summary tables from the run outputs | `figures/plot_simulation_supplement.R` | all four files in `simulation_summaries/` | Supplementary Tables S1, S2 |
| 7 | Standard errors at the empirical mixture separations | `lowsep_se_check/run_lowsep_se_check.sh`, with `sim_lowsep_se.py`, `empirical_se_ratio.R`, `assemble_summary.py` | `lowsep_se_check/*.csv` | Supplementary Fig. S4c; Tables S1, S2 |

`simulation_summaries/TableS1_design_provenance.csv` gives, per study, the sample sizes,
numbers of variants and replicates, effect ranges and estimation target. The four files in
`simulation_summaries/` are byte-identical to those written by
`plot_simulation_supplement.R` on the cluster.

## The shared harness (`code/validation/`)

Every driver uses the same few files. None of them reimplements the estimator: Stage 1 runs
through the package's own `compute_rif`, `residualize`, pooled-covariance and
block-jackknife routines, and Stage 2 through `gmoq::param_gwas`.

| File | Role |
|---|---|
| `bootstrap.sh` | Builds the software environment once: compiles Stage 1 from a G-MoQ checkout at the pinned commit, installs Stage 2 into a private R library, writes `environment.json` and `activate.sh` |
| `run_examples.py` | Loads the packages from that environment; draws one data set and fits Stage 1 with both covariances (`_stage1_example`); runs Stage 2 (`_run_stage2`) |
| `truth.py` | The data-generating models (normal with location and log-variance effects; two-component normal mixture with membership and component-location effects) and the exact population projection target of each Stage 2 parameter. Independent of the package |
| `stage2.R` | Stage 2 for one simulated variant, with the pooled or the jackknife covariance |
| `comparators.py` | Direct maximum-likelihood fits; `fit_mixture_ml` is the comparator of studies 4 and 5 |
| `benchmark.py` | Timing harness of study 2: a fresh process per run, one core, GNU `time` |
| `config/examples.json` | Scenario definitions, the 45 quantile levels and the Stage 1 settings (density estimate, minimum minor allele count, 1,000 jackknife blocks) |
| `figures/plot_stage1.R`, `figures/plot_effect_recovery.R` | Working plots made at run time |

Each file starts with a header giving what it does, the study and figure it serves, its
inputs and outputs, and where and how long it runs.

## How each study was run

Every run used one core, with all numerical libraries limited to one thread. Before running,
set the user settings at the top of each job script, or export them: `GMOQ_SIM_ROOT`, the
directory that contains `validation/` (here `analysis/08_simulations/code`); `GMOQ_SIM_ENV`,
the environment made by `bootstrap.sh`; and `GMOQ_PACKAGE_ROOT`, the G-MoQ checkout. Commands
are run from `GMOQ_SIM_ROOT`; outputs go to new directories under its `.local/`.

```bash
# Once: build the environment (inside a compute-node job).
bash validation/bootstrap.sh --out <new environment directory>

# Studies 1 and 2: four workloads (N = 2,000 with 1 and 256 variants; N = 400,000 with 1 and
# 1,000 variants), three methods, three timed repeats each; data seeds 940001-940004.
# About 25 minutes.
sbatch validation/figures/stage1.sbatch 05 400000

# Study 3: four scenarios x 1,000 variants at N = 400,000, effects drawn from U(-0.1, 0.1),
# pooled and jackknife covariance; the seed of variant i in scenario s (both from 0) is
# 95,000,000 + 100,000 s + i. About 14.5 hours.
sbatch --time=16:00:00 validation/figures/stage2.sbatch 400000 1000

# Studies 4 and 5: 20 mixture data sets at each of N = 10,000, 50,000 and 400,000
# (seeds 96,000,000 + i, effect seeds 95,000,000 + i), then three matched panels of
# 100 null variants per N (seeds 97,000,000 + i). About 2.5 hours.
sbatch validation/figures/mixture_comparison.sh mixture-full-20260918-01 20 \
  10000,50000,400000 --panel-size 100 --panel-reps 3

# Study 6: N in {10,000, 50,000, 400,000} x MAF in {0.05, 0.01, 0.001} x {null, c = 0.1},
# 200 data sets per cell; the seed of replicate r in cell c (both from 0) is
# 98,100,000 + 10,000 c + r. About 2 hours.
sbatch validation/figures/rare_variant_se.sh rare-se-20260918-01 200

# Summary tables of studies 1-6 (reads the four run directories above).
Rscript --vanilla validation/figures/plot_simulation_supplement.R outputs/simulation_results
```

Study 7: `lowsep_se_check/run_lowsep_se_check.sh` runs everything (a Stage 1 rerun for JCV and
MCV VP1 keeping the jackknife payload; Stage 2 with both covariances; the simulation; the
combined table). `empirical_se_ratio.R` is the real-data arm and `sim_lowsep_se.py` the
simulation (two-component normal mixtures at the empirical separations, N = 10,000 and
400,000, MAF 0.30 and 0.05, null and an upper-component shift with |z| ≈ 5, 200 replicates in
each of 40 cells; the seed of replicate r in cell c, counting cells from 0, is
923,000,000 + 10,000 c + r). `assemble_summary.py` joins the two arms. The simulation uses
the harness in `code/validation/`.

The per-replicate outputs of the runs are not included because of their size; the summaries
in `inputs/` are what the paper uses.

## Running the code elsewhere

- Paths that are not user settings are the authors' cluster paths: the `#SBATCH --chdir` and
  `--output` lines; the module and conda lines of `bootstrap.sh`; the run directory names read
  by `plot_simulation_supplement.R`; and the paths in `config/examples.json` (JSON allows no
  comments; its `package_root` must match `GMOQ_PACKAGE_ROOT`).
- The Python drivers and plotting scripts check that they run inside a Slurm job on the
  run directory and write only below `.local/`.
  To run elsewhere, edit those checks (`guard()` in `simulate_stage2.py`; the first lines of
  `main()` in `run_stage1.py` and `stage1_methods.py`; `_compute_guard()` in
  `export_stage1.py`; the first lines of the three R plotting scripts). They do not affect the
  results.
- The harness also holds code for the eight small validation examples that preceded these
  studies (`run_examples.py` run as a script, `config/examples.json`) and a DGLM comparator
  (`comparators.fit_dglm`, whose `dglm.R` is not included). The paper uses neither.
