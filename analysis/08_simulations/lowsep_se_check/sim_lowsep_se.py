#!/usr/bin/env python3
"""Standard errors at the empirical mixture separations: simulation arm.

What: for the three parameters of the two-component mixture model (gamma, membership;
beta_1, lower-component location; beta_2, upper-component location), at the separations
of the paper's empirical mixtures, compares pooled Stage 1 covariance standard errors
(Sigma_0 / d_resid propagated through OLS Stage 2) with 1,000-block delete-block jackknife
standard errors and with the Monte Carlo SD of the estimates.
Paper: Supplementary Information, study 7; Supplementary Fig. S4; Supplementary Tables S1, S2
Input: none (simulated data); the G-MoQ packages and the simulation harness
  in ../code/validation/ (run_examples.py, truth.py, stage2.R, config/examples.json);
  pass that directory as --validation-dir.
Output: <out-dir>/sim_summary.csv (one row per cell x parameter), sim_meta.json
Runs: on a cluster, 40 cells x 200 replicates, about 3 hours on 32 cores

Design (common variants only):
* Two-component normal mixture, sd1 = sd2 = 1, separation k within-component SDs,
  upper-component weight pi, marginal mean zero. The geometries match the empirical
  cases (testosterone unshrunk and shrunk by 0.2, 0.5 and 0.7; JCV VP1; MCV VP1), with
  the testosterone separations computed from the true sex labels.
* N = 10,000 for every geometry, and N = 400,000 as well for the four testosterone
  geometries. MAF 0.30 and 0.05.
* Two effect arms: null (gamma = delta1 = delta2 = 0), and a shift of the upper-component
  location sized for an expected |z| of about 5.
* 200 replicates per cell. Each replicate draws one data set, runs Stage 1 once (pooled
  Sigma_0 and the 1,000-block jackknife from the same RIF) and Stage 2 twice (pooled,
  jackknife).
* Coverage is of the fixed population projection target theta_pop.

Nothing here reimplements the estimator; Stage 1 and Stage 2 are the package code.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
import traceback
from multiprocessing import get_context
from pathlib import Path
from typing import Any

import numpy as np

TAU_GRID = np.linspace(0.02, 0.98, 45)

# (label, separation k, upper-component weight pi, N values)
GEOMETRIES = (
    ("testosterone_unshrunk", 5.89, 0.48, (10_000, 400_000)),
    ("testosterone_shrink_0.2", 4.75, 0.48, (10_000, 400_000)),
    ("testosterone_shrink_0.5", 2.99, 0.48, (10_000, 400_000)),
    ("testosterone_shrink_0.7", 1.81, 0.48, (10_000, 400_000)),
    ("jcv_vp1", 2.01, 0.43, (10_000,)),
    ("mcv_vp1", 3.07, 0.48, (10_000,)),
)
MAFS = (0.30, 0.05)
TARGET_Z = 5.0

_CTX: dict[str, Any] = {}


def build_cells(n_filter: set[int] | None) -> list[dict[str, Any]]:
    cells = []
    for label, sep, pi, ns in GEOMETRIES:
        for n in ns:
            if n_filter and n not in n_filter:
                continue
            for maf in MAFS:
                for effect in ("null", "upper_shift"):
                    delta2 = 0.0
                    if effect == "upper_shift":
                        # component-location SE ~ 1 / sqrt(N 2pq pi) at sd = 1
                        delta2 = TARGET_Z / math.sqrt(n * 2 * maf * (1 - maf) * pi)
                    case = {
                        "family": "mixture", "maf": maf, "a": 0.0, "c": 0.0,
                        "gamma": 0.0, "delta1": 0.0, "delta2": float(delta2),
                        "pi0": 1.0 - pi, "mu1": -pi * sep, "mu2": (1.0 - pi) * sep,
                        "sd1": 1.0, "sd2": 1.0,
                    }
                    cell = f"{label}_n{n}_maf{maf:.2f}_{effect}"
                    case["id"] = cell
                    cells.append({"cell": cell, "geometry": label, "separation": sep,
                                  "pi_upper": pi, "n": n, "maf": maf, "effect": effect,
                                  "delta2": float(delta2), "case": case})
    return cells


def _run_one(task: dict[str, Any]) -> dict[str, Any]:
    rx = _CTX["run_examples"]
    env = _CTX["environment"]
    spec = task["spec"]
    rep, seed = task["rep"], task["seed"]
    rep_dir = _CTX["out_root"] / "work" / spec["cell"] / f"rep-{rep:04d}"
    rep_dir.mkdir(parents=True, exist_ok=True)
    rec: dict[str, Any] = {"cell": spec["cell"], "rep": rep, "seed": seed}
    t0 = time.perf_counter()
    try:
        payload = rx._stage1_example(
            _CTX["core"], _CTX["truth"], dict(spec["case"]), spec["n"], seed,
            TAU_GRID, _CTX["config"], _CTX["sigma0"],
            population=_CTX["population"][spec["cell"]],
        )
        payload["package"] = _CTX["package_meta"]
        inp = rep_dir / "stage1_summary.json"
        rx._write_json(inp, payload)
        for source in ("pooled", "jackknife"):
            outp = rep_dir / f"stage2-{source}.json"
            res = rx._run_stage2(
                env["rscript"], _CTX["stage2_script"], inp, outp, rep_dir, source,
                _CTX["package_root"], _CTX["package_meta"]["package_commit"],
                _CTX["package_meta"]["package_version"], Path(env["path"]),
                env["r_lib"], 900.0,
            )
            est = res["estimate"]
            names = est["parameter_names"]
            rec[source] = {
                "names": names,
                "theta": [float(est["estimate"][k]) for k in names],
                "se": [float(est["se"][k]) for k in names],
                "lo": [float(est["ci_lower"][k]) for k in names],
                "hi": [float(est["ci_upper"][k]) for k in names],
                "df": [float(est["df"][k]) for k in names] if est.get("df") else None,
                "theta_pop": [float(res["target"]["theta_pop"][k]) for k in names],
            }
            for suffix in ("json", "stdout", "stderr"):
                p = rep_dir / f"stage2-{source}.{suffix}"
                if p.exists():
                    p.unlink()
        inp.unlink()
        rec["status"] = "ok"
    except Exception as exc:  # noqa: BLE001 - failures are retained and counted
        rec.update(status="failed", error=f"{type(exc).__name__}: {exc}",
                   traceback=traceback.format_exc())
        (rep_dir / "failure.json").write_text(json.dumps(rec, indent=2) + "\n")
    rec["seconds"] = time.perf_counter() - t0
    try:
        rep_dir.rmdir()
    except OSError:
        pass
    return rec


def aggregate(records: list[dict[str, Any]], cells: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for spec in cells:
        recs = [r for r in records if r["cell"] == spec["cell"]]
        ok = [r for r in recs if r["status"] == "ok"]
        if not ok:
            continue
        names = ok[0]["pooled"]["names"]
        R = len(ok)
        for j, name in enumerate(names):
            theta = np.array([r["pooled"]["theta"][j] for r in ok])
            target = float(np.mean([r["pooled"]["theta_pop"][j] for r in ok]))
            emp_sd = float(np.std(theta, ddof=1))
            se_p = np.array([r["pooled"]["se"][j] for r in ok])
            se_j = np.array([r["jackknife"]["se"][j] for r in ok])
            row = {
                "cell": spec["cell"], "geometry": spec["geometry"],
                "separation": spec["separation"], "pi_upper": spec["pi_upper"],
                "n": spec["n"], "maf": spec["maf"], "effect": spec["effect"],
                "delta2": spec["delta2"], "parameter": name, "replicates_ok": R,
                "replicates_failed": len(recs) - R, "target_theta_pop": target,
                "mean_estimate": float(np.mean(theta)),
                "bias": float(np.mean(theta) - target),
                "empirical_sd": emp_sd,
                "empirical_sd_rel_mcse": 1.0 / math.sqrt(2.0 * (R - 1)),
                "median_pooled_over_jackknife_se": float(np.median(se_p / se_j)),
                "q05_pooled_over_jackknife_se": float(np.quantile(se_p / se_j, 0.05)),
                "q95_pooled_over_jackknife_se": float(np.quantile(se_p / se_j, 0.95)),
            }
            for source, se in (("pooled", se_p), ("jackknife", se_j)):
                lo = np.array([r[source]["lo"][j] for r in ok])
                hi = np.array([r[source]["hi"][j] for r in ok])
                tp = np.array([r[source]["theta_pop"][j] for r in ok])
                cover = (lo <= tp) & (tp <= hi)
                cov = float(np.mean(cover))
                z = (theta - tp) / se
                row.update({
                    f"{source}_mean_se": float(np.mean(se)),
                    f"{source}_rms_se": float(np.sqrt(np.mean(se ** 2))),
                    f"{source}_rms_se_over_empirical_sd": float(np.sqrt(np.mean(se ** 2)) / emp_sd),
                    f"{source}_coverage95": cov,
                    f"{source}_coverage95_mcse": math.sqrt(max(cov * (1 - cov), 1e-12) / R),
                    f"{source}_sd_of_z": float(np.std(z, ddof=1)),
                })
            rows.append(row)
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--validation-dir", type=Path, required=True,
                    help="directory holding run_examples.py, truth.py, stage2.R and config/")
    ap.add_argument("--environment-manifest", type=Path, required=True,
                    help="environment.json naming the Python, Rscript, R library and compiled Stage 1")
    ap.add_argument("--replicates", type=int, default=200)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--blocks", type=int, default=1000)
    ap.add_argument("--n", type=int, action="append", default=None,
                    help="restrict to these N values (repeatable)")
    ap.add_argument("--cell-regex", default=None, help="restrict cells (smoke tests)")
    args = ap.parse_args(argv)

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    validation = args.validation_dir.resolve()
    sys.path.insert(0, str(validation))
    import run_examples  # noqa: PLC0415

    package_root = Path(run_examples.PINNED_PACKAGE_ROOT)
    environment = run_examples.load_environment_manifest(args.environment_manifest, package_root)
    core, provenance, sigma0 = run_examples.load_pinned_core(package_root, environment)
    truth = run_examples._load_truth_module(validation.parent)
    examples = json.loads((validation / "config" / "examples.json").read_text())
    config = dict(examples)
    config["n_blocks"] = int(args.blocks)

    cells = build_cells(set(args.n) if args.n else None)
    if args.cell_regex:
        import re  # noqa: PLC0415
        cells = [c for c in cells if re.search(args.cell_regex, c["cell"])]
    if not cells:
        raise RuntimeError("no cells selected")

    population = {}
    for spec in cells:
        population[spec["cell"]] = truth.population_truth(dict(spec["case"]), TAU_GRID)

    _CTX.update({
        "run_examples": run_examples, "core": core, "truth": truth, "sigma0": sigma0,
        "config": config, "out_root": out_dir, "environment": environment,
        "stage2_script": validation.parent / examples["stage2_script"], "package_root": package_root,
        "population": population,
        "package_meta": {
            "package_commit": run_examples.PINNED_PACKAGE_COMMIT,
            "package_version": run_examples.PINNED_PACKAGE_VERSION,
            "stage1_distribution_version": run_examples.PINNED_STAGE1_VERSION,
        },
    })

    tasks = []
    for ci, spec in enumerate(cells):
        for rep in range(1, args.replicates + 1):
            tasks.append({"spec": spec, "rep": rep, "seed": 923_000_000 + ci * 10_000 + rep})
    # Interleave large-N tasks so the pool stays balanced.
    tasks.sort(key=lambda t: (t["rep"], -t["spec"]["n"]))
    print(f"[setup] cells={len(cells)} tasks={len(tasks)} workers={args.workers} "
          f"blocks={args.blocks}", flush=True)

    started = time.perf_counter()
    records = []
    with (out_dir / "records.jsonl").open("w") as fh, get_context("fork").Pool(args.workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(_run_one, tasks, chunksize=1), start=1):
            records.append(rec)
            fh.write(json.dumps(rec) + "\n")
            if i % 100 == 0 or i == len(tasks):
                fh.flush()
                nf = sum(r["status"] != "ok" for r in records)
                print(f"[run] {i}/{len(tasks)} failed={nf} "
                      f"elapsed={time.perf_counter() - started:.0f}s", flush=True)

    rows = aggregate(records, cells)
    with (out_dir / "sim_summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    meta = {
        "replicates": args.replicates,
        "blocks": args.blocks, "workers": args.workers, "cells": len(cells),
        "failed": sum(r["status"] != "ok" for r in records),
        "elapsed_seconds": time.perf_counter() - started,
        "package_commit": run_examples.PINNED_PACKAGE_COMMIT,
        "package_version": run_examples.PINNED_PACKAGE_VERSION,
        "stage1_version": run_examples.PINNED_STAGE1_VERSION,
    }
    (out_dir / "sim_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"[done] {out_dir / 'sim_summary.csv'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
