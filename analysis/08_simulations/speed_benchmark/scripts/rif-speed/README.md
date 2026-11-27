# Matched Stage 1 timing

This is the Supplementary Figure S2 speed benchmark, adapted from the frozen
`gmoq-paper-code/analysis/08_simulations/code/validation` runner to the current
`gmoq` package. It performs no recovery, coverage or power campaign. Coefficient
equality and the pooled versus ordinary OLS variance identity are benchmark
integrity checks.

The full benchmark retains the original four workloads:

| N | Variants | Data seed |
|---:|---:|---:|
| 2,000 | 1 | 940001 |
| 2,000 | 256 | 940002 |
| 400,000 | 1 | 940003 |
| 400,000 | 1,000 | 940004 |

Each workload has three fresh-process repeats of ordinary RIF-OLS, G-MoQ pooled
and G-MoQ jackknife. The ordinary baseline uses NumPy `lstsq` once per SNP with
45 responses. The pooled arm uses the production `make_pooled_bin_plan`,
`compute_pooled_binned` and `plan.sigma0()` APIs. The jackknife arm retains the
fast block plan/finalisation API and all 1,000 balanced blocks, including packed
covariance computation. Every arm returns the same coefficient and SE table.

The original timing boundary is preserved: all children load/check identical
inputs and compute the same dense RIF before fit timing. Residualisation, plan
construction, coefficients and SEs are inside fit timing. Whole-process GNU
time includes package import/pin checks, input reads, RIF setup and output
writing. Thus this is the matched published workload; pooled whole-process
time includes dense RIF materialisation that a minimal production pooled scan
does not need. Timings are not conditional quantile regression timings.

`benchmark.py` is copied unchanged from the existing benchmark harness. Its
interleave seed is 940100. Processes are pinned to one allowed CPU with all
numerical-library threads set to one. GNU time reports whole-process wall time
and peak RSS. Failed or timed-out arms remain on disk and prevent COMPLETE.
Output directories are exclusive and cannot be reused.

Run on a BP compute node under Slurm, using `stage1_python` from PIN.json:

```sh
"$GMOQ_STAGE1_PYTHON" /path/to/rif-speed/run_stage1.py \
  --pin-json /user/work/fh6520/gmoq-paper-v1/Rlib/gmoq-166635d4/PIN.json \
  --out /private/path/s2-smoke --preset smoke

"$GMOQ_STAGE1_PYTHON" /path/to/rif-speed/run_stage1.py \
  --pin-json /user/work/fh6520/gmoq-paper-v1/Rlib/gmoq-166635d4/PIN.json \
  --out /private/path/s2-n400000-m1000 --workload n400000_m1000
```

`--preset smoke` runs both original N=2,000 workloads (18 fresh processes).
`--preset full` runs all four (36 processes). `--workload` selects one original
size regardless of preset (9 processes), allowing the coordinator to schedule
different sizes in separate jobs. Each command retains exactly three repeats.
The default per-child timeout is 1,800 seconds. Every run verifies the full
source commit, a clean tracked checkout, exact extension SHA256 and interpreter
from PIN.json. No jobs are submitted by these scripts.

Outputs are `tables/S2_timings.csv` (each repeat), `S2_mean_timings.csv` (the
existing figure table schema), `S2_integrity_checks.csv`, and
`verification.json`, plus input hashes, command manifest, provenance, and the
fresh-process outputs/logs. Input arrays are simulated but remain on BP.
The full set occupies about 3.3 GB in input arrays; generation transiently
holds both integer and float64 genotype arrays, about 6.4 GB for the largest
cell. Original OLS timing was about 371 seconds per largest-cell repeat; use
the smoke run's measured RSS/elapsed and the coordinator's resource judgement
before scheduling full cells.
