# Paired mixture fitting timings

This runner imports the frozen null-mixture generator and direct maximum-likelihood comparator without modifying either. It uses the released, pinned G-MoQ production pooled kernel and public Stage 2. Four methods share the same subjects, phenotype and variants in every cell: pooled G-MoQ, G-MoQ with the original 1,000-block jackknife, the canonical six-start conditional mixture ML fit, and the phenotype-only warm-start ML sensitivity with canonical fallback.

There are nine production cells: N = 10,000, 50,000 and 400,000, each at frozen phenotype seeds 97,000,000, 97,000,001 and 97,000,002. Genotype seeds are the phenotype seed plus 7,919. Each genotype panel is drawn with shape `(N, 100)` exactly as in the frozen benchmark. Pilots take the first columns of that same panel. Each method starts in a fresh process with one OpenMP/BLAS thread, and fits the whole panel sequentially or as a production batch. The G-MoQ baseline mixture is fitted to the observed phenotype, and its fitted weights are included in the fitting cost. The warm ML baseline is also fitted once to the phenotype and included once. Neither direct arm receives generating effects, latent memberships, or fitted G-MoQ results.

`PIN.json` must contain `commit`, `stage1_root`, `stage1_extension_sha256` and `r_lib`; `stage1_python`, `rscript`, `r_package` and `package_version` are optional. By default the Python child uses `stage1_python` from the pin and the R child uses `Rscript` on PATH. The Stage 1 source commit, clean tracked state and compiled extension hash are checked before fitting. Generator, comparator, block-assignment, adapter and runner hashes are recorded. R is run with `--vanilla` and the pinned package library is supplied explicitly.

The analytic model-density weights reproduce the frozen comparison's convention. Stage 2 columns `gamma`, `beta_1` and `beta_2` map to membership logit slope, lower-component location slope `delta1` and upper-component location slope `delta2`. The membership weight is `-p1*(1-p1)*(F1-F2)/f`, matching the comparator's `logit(p1) + gamma*(g - 2*maf)`. Both direct ML arms estimate the same eight free parameters with unchanged bounds, optimiser settings, observed numerical Hessian, SEs and failure checks. The jackknife uses the frozen `balanced_block_ids` function and seed, the public fast block plan/kernel with panel width 16, and the frozen Stage 2 options (`n_blocks=1000`, pooled-tau OLS convention, saddlepoint Q calibration).

Run code only on a BP compute node. A small end-to-end smoke test is:

```bash
python paired_speed.py cell \
  --paper-root /user/work/fh6520/Gmoq_paper-speedprof-6bf1af5 \
  --pin-json /user/work/fh6520/gmoq-paper-v1/Rlib/gmoq-166635d4/PIN.json \
  --out /absolute/path/to/smoke --n 1000 --rep 1 --variants 2 --smoke \
  --validate-stage1
```

The smoke test's pooled dense-kernel identity check is outside its fitting timer, but its whole-process wall includes the check and is flagged ineligible as a production wall observation. Production runs omit that diagnostic. A first pilot has three independent N cells, one frozen seed and three variants each:

```bash
export PAIRED_RUNNER_DIR=/absolute/path/to/paired-speed
sbatch --array=0-2%1 --mem=<pilot-memory> --time=<pilot-time> \
  --output=/absolute/path/to/logs/pilot-%A_%a.out \
  /absolute/path/to/paired_speed.sbatch pilot /absolute/path/to/pilot
python paired_speed.py collate --out /absolute/path/to/pilot --reps 1 --variants 3
```

Review the pilot quality gate and measured elapsed/MaxRSS before choosing production resources. Set memory near 1.5 times the largest measured peak and a time limit based on the slowest cell's measured full-panel projection, with a documented margin. The production array is bounded explicitly; this example runs at most two cells concurrently:

```bash
sbatch --array=0-8%2 --mem=<measured-memory> --time=<measured-time> \
  --output=/absolute/path/to/logs/full-%A_%a.out \
  /absolute/path/to/paired_speed.sbatch full /absolute/path/to/full
python paired_speed.py collate --out /absolute/path/to/full --reps 3 --variants 100
```

The wrapper requires `PAIRED_RUNNER_DIR`, because Slurm copies the submission script to a spool directory. It defaults to the supplied frozen paper root and current package pin. `PAPER_ROOT`, `PIN_JSON`, `PYTHON_BIN`, `TIMEOUT_SECONDS` and `METHODS` can be set explicitly. `R_LIBS_USER` defaults to the existing shared CRAN library `/user/work/fh6520/r_packages`. The wrapper does not create an allocation or submit another job. Launch one paired cell directly with `cell --n <N> --rep <1|2|3> --variants 100`. `--methods gmoq_jackknife` can add that arm independently: completed other arms are discovered and retained when updating cell summaries. An interrupted method preserves its logs and variant checkpoints, and a retry starts that entire method in a new numbered attempt. Complete methods are reused without fitting again. The successful fresh attempt supplies its timing; incomplete attempts remain as failures, rather than silently stitching partial timing segments.

Each cell writes atomic input/provenance manifests, per-method logs and results, raw per-variant ML diagnostics, `timings.csv`, `estimates.csv`, `warm_quality.csv`, `quality.json`, and a completion record. All comparator diagnostics, fallback counts and unsuccessful fits are retained. An input/source/pin mismatch refuses to reuse the output. Concurrent jobs must have different cell IDs; do not launch two writers for the same cell.

`fit_seconds` excludes import/startup, data loading and checkpoint serialization. For direct ML it is the sum of per-variant fit plus Hessian timers, with the warm baseline added once. For G-MoQ it is the sum of baseline fitting, Stage 1 and covariance computation, weight construction, and the public Stage 2 call. `end_to_end_fit_seconds` additionally includes method-region Python/R dispatch and output/checkpoint IO. `process_wall_seconds` is the parent's complete fresh-process wall time, including interpreter startup, imports, source checks, input loading, fitting and output IO. The G-MoQ R subprocess wall and the individual Python/R components are retained separately in `result.json`. ML checkpoint logs flush after each SNP and are fsynced once at method completion; checkpoint work is excluded from `fit_seconds`. Linux GNU time records peak RSS, user and system CPU time.

Stage 2 reuse is measured separately after the full fitting arms, in its own R process, using the pooled Stage 1 summary and fitted weights. One warm-up fit precedes 20 timed calls. Raw timings, mean and median are retained, with exact estimate/SE reproducibility checks. Reuse excludes Stage 1, phenotype-only baseline fitting, JSON parsing, package loading and weight construction. It is not an end-to-end G-MoQ timing.

Warm speed comparisons are eligible only when every paired variant has accepted comparator fits, finite positive SEs, the expected start count, absolute log-likelihood difference per observation at most `1e-8`, maximum estimate difference at most `0.01` canonical SE, and maximum relative SE difference at most `0.01`. These thresholds are declared before fitting and are recorded with every cell. A failure is retained and stops the cell from being called complete; it does not cause an automatic relaxed threshold or another optimiser. This is a fitting-cost experiment on null panels, not a coverage, calibration or recovery study. The runner makes no figures or speed claims.
