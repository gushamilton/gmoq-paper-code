# Analysis environment

The aggregate-only figure and table rebuild uses the checked-in `results/` and the
dependencies in the main README. Recomputing participant analyses or simulations also
requires the historical G-MoQ packages below. The current package repository has since
renamed both packages to `gmoq`; installing its default branch does not recreate the paper's
environment.

| Use | G-MoQ commit | Python module / R package |
|---|---|---|
| Stage 1 GWAS, sample loaders and simulations | `9c47e7ba56370e708810fff8c2c9ad3b723cbccf` | `fungwas_stage1` (distribution `fungwas-stage1`, 0.1.0) |
| Empirical Stage 2 and mixture weight-design simulation | `9fee67d1e510b81368d33b5b0ee223792bc6a3bd` | `fungwasStage2` 0.5.0 |
| Original `analysis/08_simulations/code/validation/` harness | `9c47e7ba56370e708810fff8c2c9ad3b723cbccf` (both stages) | `fungwasStage2` 0.4.0 |

The earlier empirical Stage 2 history is recorded separately in `config.yaml`; the original
simulation harness deliberately retains its own 0.4.0 pin. The later speed benchmark has
its own source and environment under `analysis/08_simulations/speed_benchmark/`.

Use separate source checkouts and R libraries. Run installation and analysis on a compute
node when using BluePebble. For example, after cloning the public G-MoQ repository:

```sh
git -C /path/to/GMoQ worktree add --detach /path/to/GMoQ-stage1 9c47e7ba56370e708810fff8c2c9ad3b723cbccf
git -C /path/to/GMoQ worktree add --detach /path/to/GMoQ-stage2 9fee67d1e510b81368d33b5b0ee223792bc6a3bd
# In a Python environment with the C++/Armadillo/OpenMP toolchain:
python -m pip install '/path/to/GMoQ-stage1/stage1[bgen]'
mkdir -p /path/to/Rlib-stage2-0.5 /path/to/Rlib-validation-0.4
R CMD INSTALL --library=/path/to/Rlib-stage2-0.5 /path/to/GMoQ-stage2/stage2
R CMD INSTALL --library=/path/to/Rlib-validation-0.4 /path/to/GMoQ-stage1/stage2
```

Stage 1 needs NumPy, SciPy and a compiled extension built with pybind11, Armadillo and
OpenMP; BGEN input also needs `bgen` (and `zstandard` for zstd-compressed BGENs). The Python
PRS and covariate helpers additionally use pandas. The original validation harness also
requires pytest and GNU `time`/`sha256sum`. Its cluster-specific `bootstrap.sh` builds the
extension in place and records source, executable and library hashes; update its user
settings and module/conda paths for your cluster. It refuses execution outside a compute
allocation and refuses a missing compiled extension.

The shared R analysis loader needs data.table and yaml. The mixture and HbA1c analyses
also use jsonlite, sn and mixsmsn; exact Davies/Imhof checks use CompQuadForm. Install the
Stage 2 DESCRIPTION dependencies, including Rcpp and RcppArmadillo, before `R CMD INSTALL`.
Other optional analyses declare their own dependencies in their scripts (for example
conquer for quantile-regression comparisons and Bioconductor annotation packages for gene
annotation). These are separate from the aggregate-only rebuild.

Set `package.stage1.root` and `package.stage1.python`, and the empirical `package.stage2.lib`
in `config.yaml`; update the private data paths before any participant analysis. R scripts
read `GMOQ_CONFIG` when set. Stage 2 callers use the configured R package and library,
check its version and resolve exported functions from that namespace. `GMOQ_STAGE2_LIB`
overrides the library for local checks. `GMOQ_STAGE1_PYTHON` selects the Stage 1 shell
driver's interpreter; Python helpers launched from R use `package.stage1.python`.

For an explicitly modern-package run, set `GMOQ_STAGE1_MODULE=gmoq` and
`GMOQ_STAGE2_PACKAGE=gmoq`, select the corresponding Python and R library, and update the
recorded package pin/version. This is a separate software environment. The original
validation harness enforces its historical pin and does not accept these overrides.
No CLI option translation is needed: the historical Stage 1 accepts the paper's `--snps`,
`--blocks`, `--covar-col`, `--loco-preds`, `--loco-mode` and `--cov-mode` options directly.
The `pooled` and `pooled-only` covariance modes retain their original meanings.

HbA1c producers export cohort labels under `sample_group`. The mixture analysis also
accepts the published aggregate `sample_counts.csv` schema, whose category column is
named `sample`. The result writer continues to reject identifier-like column names,
including `sample`, FID, IID and participant IDs.
