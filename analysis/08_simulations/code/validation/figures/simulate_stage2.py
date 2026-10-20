#!/usr/bin/env python3
"""Stage 2 recovery of simple simulated effects.

What: for four scenarios (normal phenotype with a location effect, a dispersion effect, or
  both; two-component normal mixture with a shift of the lower component) and 1,000
  variants per scenario, draws one data set per variant with an effect drawn from
  U(-0.1, 0.1), runs Stage 1 and then Stage 2 with the pooled and the jackknife covariance,
  and records the estimates beside two targets: the generating effect and the exact
  population projection of the quantile effects onto the model's weight columns
  (truth.py).
Paper: Supplementary Fig. S3; Supplementary Table S2 (study 3)
Input: none (simulated data); config/examples.json (scenario definitions, quantile grid);
  the G-MoQ packages through the environment made by ../bootstrap.sh
Output: <out-dir>/results.csv (one row per variant x parameter x covariance), config.json,
  summary.json
Runs: on a cluster, one core, about 14.5 hours for N = 400,000, through stage2.sbatch:
  python validation/figures/simulate_stage2.py --n 400000 --snps 1000 \
    --environment-manifest $GMOQ_VALIDATION_ENV/environment.json --out-dir <new dir under .local/>
  Seed of variant i (from 0) in scenario s (from 0): 95,000,000 + 100,000 s + i; its effect
  is drawn with seed + 10,000,000.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import socket
import traceback
import time

import numpy as np

REPO = Path(__file__).resolve().parents[2]
PARAMS = {
    "normal": (("beta_mu", "a"), ("beta_sigma2", "c")),
    "mixture": (("gamma", "gamma"), ("beta_1", "delta1"), ("beta_2", "delta2")),
}
SCENARIOS = ("normal_location", "normal_dispersion", "normal_both", "mixture_left")
FIELDS = ("scenario", "snp_id", "seed", "n", "parameter", "covariance_source",
          "generating_effect", "projection_target", "estimate", "se", "ci_lower",
          "ci_upper", "p_value", "status", "reason")


def guard() -> None:
    """Formerly restricted runs to the authors' cluster; kept as a no-op."""


def case_with_effect(base, scenario, seed):
    rng = np.random.default_rng(seed + 10_000_000)
    c = dict(base)
    c["id"] = scenario
    if scenario == "normal_location":
        c.update(a=float(rng.uniform(-.10, .10)), c=0.0)
    elif scenario == "normal_dispersion":
        c.update(a=0.0, c=float(rng.uniform(-.10, .10)))
    elif scenario == "normal_both":
        c.update(a=float(rng.uniform(-.10, .10)), c=float(rng.uniform(-.10, .10)))
    else:
        c.update(gamma=0.0, delta1=float(rng.uniform(-.10, .10)), delta2=0.0)
    return c


def target_map(population):
    return dict(zip(population["theta_names"], map(float, population["theta"])))


def emit_rows(writer, scenario, snp_id, seed, n, params, case, fits=None,
              targets=None, reason=None):
    failed = reason is not None
    for source in ("pooled", "jackknife"):
        fit = None if fits is None else fits[source]
        for parameter, field in params:
            row = {"scenario": scenario, "snp_id": snp_id, "seed": seed, "n": n,
                   "parameter": parameter, "covariance_source": source,
                   "generating_effect": float(case[field]),
                   "projection_target": "" if targets is None else targets.get(field, ""),
                   "estimate": "", "se": "", "ci_lower": "", "ci_upper": "",
                   "p_value": "", "status": "failed" if failed else "success",
                   "reason": "" if reason is None else reason}
            if not failed:
                est = fit["estimate"]
                row.update(estimate=est["estimate"][parameter], se=est["se"][parameter],
                           ci_lower=est["ci_lower"][parameter], ci_upper=est["ci_upper"][parameter],
                           p_value=est["p_value"][parameter],
                           projection_target=targets[field])
            writer.writerow(row)


def run(args) -> int:
    import sys
    sys.path.insert(0, str(REPO / "validation"))
    import run_examples

    package_root = Path(run_examples.PINNED_PACKAGE_ROOT)
    environment = run_examples.load_environment_manifest(args.environment_manifest, package_root)
    core, package, sigma0 = run_examples.load_pinned_core(package_root, environment)
    truth = run_examples._load_truth_module(REPO)
    examples = json.loads((REPO / "validation/config/examples.json").read_text())
    bases = {c["id"]: c for c in examples["cases"]}
    taus = np.asarray(examples["taus"], dtype=np.float64)
    stage2 = REPO / examples["stage2_script"]
    source_commit = run_examples._git_head(REPO)

    out = args.out_dir.resolve()
    local = (REPO / ".local").resolve()
    if out == local or local not in out.parents:
        raise RuntimeError(f"--out-dir must be a new directory below {local}")
    if out.exists():
        raise RuntimeError(f"refusing to reuse output directory: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir()
    config = {
        "schema_version": 1, "source_commit": source_commit,
        "package_root": str(package_root), "package_commit": run_examples.PINNED_PACKAGE_COMMIT,
        "package_version": run_examples.PINNED_PACKAGE_VERSION,
        "stage1_distribution_version": run_examples.PINNED_STAGE1_VERSION,
        "package": package, "environment_manifest": environment["path"],
        "stage2_script": str(stage2), "scenarios": list(SCENARIOS),
        "case_definitions": {name: bases[name] for name in SCENARIOS},
        "snps": args.snps, "n": args.n, "seed_base": 95_000_000,
        "effect_seed_offset": 10_000_000, "covariance_sources": ["pooled", "jackknife"],
        "tau_cov_estimator": "ols",
    }
    run_examples._write_json(out / "config.json", config)

    started = time.perf_counter()
    attempted = failed = 0
    failed_by_scenario = {name: 0 for name in SCENARIOS}
    results = out / "results.csv"
    with results.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        handle.flush()
        for case_index, scenario in enumerate(SCENARIOS):
            base = bases[scenario]
            params = PARAMS[base["family"]]
            for i in range(args.snps):
                snp_id, seed = i + 1, 95_000_000 + case_index * 100_000 + i
                case_dir = out / scenario / f"snp-{snp_id:04d}"
                case_dir.mkdir(parents=True)
                case = case_with_effect(base, scenario, seed)
                population = None
                targets = None
                try:
                    population = truth.population_truth(case, taus)
                    targets = target_map(population)
                    payload = run_examples._stage1_example(
                        core, truth, case, args.n, seed, taus, examples,
                        sigma0, population=population)
                    payload["package"] = package
                    input_path = case_dir / "stage1_summary.json"
                    run_examples._write_json(input_path, payload)
                    fits = {}
                    for source in ("pooled", "jackknife"):
                        fits[source] = run_examples._run_stage2(
                            environment["rscript"], stage2, input_path,
                            case_dir / f"stage2-{source}.json", case_dir, source,
                            package_root, run_examples.PINNED_PACKAGE_COMMIT,
                            run_examples.PINNED_PACKAGE_VERSION,
                            Path(environment["path"]), environment["r_lib"], 600.0)
                    emit_rows(writer, scenario, snp_id, seed, args.n, params, case, fits, targets)
                except Exception as exc:
                    failed += 1
                    failed_by_scenario[scenario] += 1
                    reason = f"{type(exc).__name__}: {exc}".replace("\n", "; ")
                    run_examples._write_json(case_dir / "traceback.json", {
                        "scenario": scenario, "snp_id": snp_id, "seed": seed,
                        "error": reason, "traceback": traceback.format_exc(),
                    })
                    emit_rows(writer, scenario, snp_id, seed, args.n, params, case,
                              targets=targets, reason=reason)
                attempted += 1
                handle.flush()
                run_examples._write_json(out / "progress.json", {"attempted": attempted, "failed": failed, "scenario": scenario, "snp_id": snp_id, "elapsed_seconds": time.perf_counter()-started})
    summary = {"status": "ok" if failed == 0 else "failed", "attempted": attempted,
               "failed": failed, "failure_counts": failed_by_scenario, "elapsed_seconds": time.perf_counter()-started}
    run_examples._write_json(out / "summary.json", summary)
    return 0 if failed == 0 else 1


def main() -> int:
    guard()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snps", type=int, default=20)
    parser.add_argument("--n", type=int, default=50_000)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--environment-manifest", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.snps <= 1000 or args.n < 1000:
        raise ValueError("--snps must be 1..1000 and --n must be at least 1000")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
