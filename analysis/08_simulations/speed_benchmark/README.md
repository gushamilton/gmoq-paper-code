# Later computational benchmarks

The cost results in Figure 2c, Supplementary Figures S2, S4d and S5, and the cost rows of Supplementary Table S2 use G-MoQ 0.8.1, source commit `166635d49527e27f8a714247d551c17093e02a55`. They are independent of the earlier recovery and coverage simulations. `results/08_simulations/speed_benchmark/` contains the aggregate measurements and package pin; `import_simulation_summaries.R` derives the paper summaries from those files without fitting models.

## Workloads and timing boundaries

All runs use one CPU with OpenMP and numerical-library threads limited to one. All inputs are simulated; none of these benchmarks times BGEN decoding or a multi-trait genome scan.

* **Stage 1 (S2):** N = 2,000 with 1 or 256 variants and N = 400,000 with 1 or 1,000 variants, data seeds 940001-940004, three fresh-process repeats per method/workload. Ordinary RIF-OLS uses one least-squares call per variant with 45 responses. All methods construct the same dense RIF before the fit timer. Fit time includes residualisation, plan construction, coefficients and standard errors. Whole-process time additionally includes imports, pin checks, input loading, RIF setup and output writing. Peak resident memory is the arithmetic mean of the per-process peaks in MiB. `stage1_panel_times.csv` is the source for all S2 timing rows.
* **Mixture panels (S4d, S5a-c):** N = 10,000, 50,000 and 400,000, three independent 100-variant null panels per N; phenotype seeds 97000000-97000002, genotype seeds plus 7919, binomial(2, 0.3) genotypes and a balanced normal mixture with means -2 and 2 and component SD 1. The baseline mixture is fitted to the observed phenotype. G-MoQ uses 45 quantile levels, model-density weights and pooled or 1,000-block jackknife covariance. The primary comparator fits the same conditional mixture by six-start ML with numerical-Hessian standard errors. Fit time sums fitting operations, including baseline, Stage 1, weights and Stage 2 for G-MoQ, and excludes interpreter startup and storage input/output. Whole-process time includes those costs. S4d uses medians of the three panel fit times; S5a-c and Table S2 use arithmetic means for fit time, whole-process wall time and peak RSS. All 100 variants in every arm completed; no variants are selected on fit quality when computing timing means. `mixture_panel_repeats.csv` and `mixture_panel_times.csv` retain the raw panel timings and means.
* **Marginal cost (Fig. 2c, S5d):** the separate component profiler uses the same null generator, sample sizes and seeds. Stage 1 seconds per additional variant is the 300-variant time minus the mean 100-variant time, divided by 200. Direct six-start ML seconds per variant is the arithmetic mean over 10, 10 and 5 fits at the three sample sizes. Stage 2 measures one calibrated call on 100,000 in-memory columns tiled from 100 simulated variants; dividing by 100,000 yields seconds per variant. The projection multiplies these marginal costs by 1.2 million and excludes per-trait setup and storage input/output. At N = 400,000, a separate call on 1.2 million tiled columns measured 6.181 seconds; this is distinct from the 7.02-second projection based on 100,000 columns. `component_costs.csv` records both denominators.

The additional phenotype-only warm-start ML arm is a sensitivity benchmark. All 900 paired fits completed; its predeclared output-equivalence check passed for 298/300 variants at N = 400,000 because two standard-error differences exceeded 1%. Its speed ratio at this N remains missing in `mixture_panel_speedups.csv`. These checks do not exclude the canonical six-start comparator or any timing observation. `mixture_comparator_checks.csv` and `mixture_quality_status.json` retain this result.

## Source and reproduction

The `scripts/` files are the exact later harness snapshot at source commit `36dd3ad68bdf0700a6da6473326931c29b506f7d`, including Stage 1 and paired-mixture runners and the aggregate calculations. `speed_profile_gmoq.py` and `stage2_profile.R` are the exact component profilers used by that run. These are source records, not a portable installation wrapper. Their scheduler wrappers preserve the original cluster paths; inspect and configure these paths before using them. Run fitting only on a compute node. Do not rerun completed simulations just to rebuild figures.

The harness imports the frozen generator and ML comparator from the paper worktree at `f6cb6f08dc281c1b3891d4a28bdce43cf6b54ead`. The corresponding historical sources are under `../code/validation/`. The package source, extension hash and R library are checked against the pin. The paired runner's `--paper-root` must identify a directory containing that `validation/` tree; `--pin-json` supplies the installed later package pin. Full reproduction also requires the package checkout and dependencies specified by the original wrappers. Saved synthetic input arrays, per-variant diagnostics and full component raw logs remain in the recorded cluster run; only aggregate public results are included here.

To regenerate the released summaries and artwork from those aggregate measurements, from the repository root:

```sh
Rscript analysis/08_simulations/import_simulation_summaries.R
Rscript figures/fig2_testosterone_recovery.R
Rscript figures/figS02_stage1_timing.R
Rscript figures/figS04_mixture_recovery.R
Rscript figures/figS05_mixture_timing.R
Rscript tables/tableS2_simulation_results.R
```

The pre-import `inputs/simulation_summaries/` files retain the historical simulation export. Importing replaces their S2/S5 cost columns with the later measurements and preserves the recovery/coverage values. This prevents old cost rows from silently returning during a rebuild.
