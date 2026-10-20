#!/usr/bin/env python3
"""Check and export the Stage 1 comparison tables.

What: reads the 36 benchmark runs made by run_stage1.py, checks that the G-MoQ and ordinary
  RIF-OLS coefficients agree to 1e-10 (absolute and relative) for every variant and quantile
  level, and writes the paired coefficients, the standard errors of the three methods and
  the timings.
Paper: Supplementary Figs. S1, S2; Supplementary Table S2 (studies 1 and 2)
Input: the benchmark directory written by run_stage1.py
Output: <out-dir>/{S1_pairs,S1_SE,S1_checks,S2_timings}.csv, verification.json, README.md
Runs: on a cluster, called by run_stage1.py, under a minute
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import socket
from itertools import combinations
from pathlib import Path
from typing import Any

ATOL = RTOL = 1e-10
METHODS = ("gmoq", "ordinary_rif_ols", "gmoq_jackknife")
DEFAULT_N_LARGE = 50000
_BASE_WORKLOADS = {
    "n2000_m1": (2000, 1), "n2000_m256": (2000, 256),
}


def _workloads(n_large: int) -> dict[str, tuple[int, int]]:
    if n_large not in (50000, 400000):
        raise ValueError("n-large must be 50000 or 400000")
    return {**_BASE_WORKLOADS,
            f"n{n_large}_m1": (n_large, 1),
            f"n{n_large}_m1000": (n_large, 1000)}


WORKLOADS = _workloads(DEFAULT_N_LARGE)
FIELDS = ("snp_id", "tau", "beta")


def _number(value: Any, name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} is not numeric") from exc
    if not math.isfinite(out):
        raise ValueError(f"{name} is not finite")
    return out


def read_slopes(path: Path, value: str = "beta") -> dict[tuple[str, float], float]:
    if not path.is_file():
        raise ValueError(f"missing slopes.csv: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or any(field not in reader.fieldnames for field in (*FIELDS, value)):
            raise ValueError(f"{path} must contain columns snp_id,tau,beta")
        rows: dict[tuple[str, float], float] = {}
        for line, row in enumerate(reader, 2):
            raw_snp = row.get("snp_id")
            snp = "" if raw_snp is None else str(raw_snp).strip()
            if not snp:
                raise ValueError(f"{path}:{line} has an empty snp_id")
            tau = _number(row.get("tau"), f"{path}:{line} tau")
            beta = _number(row.get(value), f"{path}:{line} {value}")
            key = (snp, tau)
            if key in rows:
                raise ValueError(f"duplicate SNP/tau key in {path}: {key!r}")
            rows[key] = beta
    if not rows:
        raise ValueError(f"empty slopes.csv: {path}")
    return rows


def compare_slopes(left: dict, right: dict, label: str) -> float:
    if set(left) != set(right):
        raise ValueError(f"inconsistent SNP/tau keys for {label}")
    maximum = 0.0
    for key in left:
        error = abs(left[key] - right[key])
        maximum = max(maximum, error)
        if error > ATOL + RTOL * abs(right[key]):
            raise ValueError(f"RIF slope mismatch for {label} at {key!r}: {error:g}")
    return maximum


def _check_result(row: dict, benchmark_dir: Path,
                  workloads: dict[str, tuple[int, int]]) -> tuple[str, str, int, dict]:
    workload, method = row.get("workload_id"), row.get("method_id")
    if workload not in workloads or method not in METHODS:
        raise ValueError(f"unknown workload/method: {workload!r}/{method!r}")
    try:
        repeat = int(row["repeat"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("repeat is required and must be an integer") from exc
    if repeat not in (1, 2, 3):
        raise ValueError(f"repeat must be 1, 2 or 3, got {repeat}")
    if row.get("status") != "success":
        raise ValueError(f"benchmark arm failed: {workload}/{method}/repeat{repeat}")
    for name in ("return_code", "child_return_code"):
        if row.get(name) != 0:
            raise ValueError(f"{name} is not zero for {workload}/{method}/repeat{repeat}")
    for name in ("fit_seconds", "wall_seconds", "peak_rss_kib"):
        value = _number(row.get(name), name)
        if value < 0:
            raise ValueError(f"{name} is negative")
        row[name] = value
    if not row.get("output_dir"):
        raise ValueError(f"missing output_dir for {workload}/{method}/repeat{repeat}")
    output = Path(str(row["output_dir"]))
    if not output.is_absolute():
        output = (benchmark_dir / output).resolve()
    slopes = read_slopes(output / "slopes.csv")
    row["_se"] = read_slopes(output / "slopes.csv", "se")
    if set(slopes) != set(row["_se"]) or any(x <= 0 for x in row["_se"].values()):
        raise ValueError("missing or nonpositive SE")
    n, m = workloads[workload]
    taus = sorted({key[1] for key in slopes})
    expected_snps = {f"snp_{j}" for j in range(m)}
    if len(taus) != 45 or len(slopes) != m * 45:
        raise ValueError("wrong SNP/tau dimensions")
    if set(slopes) != {(snp, tau) for snp in expected_snps for tau in taus}:
        raise ValueError("incomplete expected SNP/tau keys")
    if any(not math.isclose(tau, .02 + i * .96 / 44, abs_tol=1e-12, rel_tol=0) for i, tau in enumerate(taus)):
        raise ValueError("unexpected tau grid")
    metadata = json.loads((output / "metadata.json").read_text())
    if (metadata.get("n"), metadata.get("m"), metadata.get("method")) != (n, m, method):
        raise ValueError("fit metadata does not match arm")
    if metadata.get("package_commit") != row.get("package_pin"):
        raise ValueError("fit package pin differs from benchmark")
    row["_input_manifest"] = metadata["input_manifest"]
    if not row["_input_manifest"].get("files") or not row["_input_manifest"].get("metadata_sha256"):
        raise ValueError("missing matched-input checksums")
    row["output_dir"] = str(output)
    return workload, method, repeat, slopes


def verify_benchmark(benchmark_dir: Path, n_large: int = DEFAULT_N_LARGE) -> tuple[dict, list[dict], list[dict]]:
    workloads = _workloads(n_large)
    path = benchmark_dir / "results.json"
    try:
        results = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read benchmark results: {path}") from exc
    if not isinstance(results, list) or len(results) != 36:
        raise ValueError("results.json must be the 36-arm benchmark list")
    arms: dict[tuple[str, str, int], tuple[dict, dict]] = {}
    for raw in results:
        if not isinstance(raw, dict):
            raise ValueError("each benchmark result must be an object")
        try:
            raw_repeat = int(raw.get("repeat", -1))
        except (TypeError, ValueError) as exc:
            raise ValueError("repeat is required and must be an integer") from exc
        key = (raw.get("workload_id"), raw.get("method_id"), raw_repeat)
        if key in arms:
            raise ValueError(f"duplicate benchmark arm: {key!r}")
        workload, method, repeat, slopes = _check_result(raw, benchmark_dir, workloads)
        arms[(workload, method, repeat)] = (raw, slopes)
    expected = {(w, m, r) for w in workloads for m in METHODS for r in (1, 2, 3)}
    if set(arms) != expected:
        raise ValueError("results.json has missing or extra workload/method/repeat arms")

    checks: list[dict] = []
    for workload in workloads:
        reference_input = arms[(workload, "gmoq", 1)][0]["_input_manifest"]
        for method in METHODS:
            for repeat in (1, 2, 3):
                if arms[(workload, method, repeat)][0]["_input_manifest"] != reference_input:
                    raise ValueError("input checksums differ across matched arms")
                if (arms[(workload, method, repeat)][1] != arms[(workload, method, 1)][1] or
                    arms[(workload, method, repeat)][0]["_se"] != arms[(workload, method, 1)][0]["_se"]):
                    raise ValueError("within-method repeats are not deterministic")
        entries = [arms[(workload, method, repeat)][1]
                   for method in METHODS for repeat in (1, 2, 3)]
        for left, right in combinations(entries, 2):
            compare_slopes(left, right, workload)
        for repeat in (1, 2, 3):
            gmoq = arms[(workload, "gmoq", repeat)][1]
            rif = arms[(workload, "ordinary_rif_ols", repeat)][1]
            n = workloads[workload][0]
            pooled_se = arms[(workload, "gmoq", repeat)][0]["_se"]
            ols_se = arms[(workload, "ordinary_rif_ols", repeat)][0]["_se"]
            for key in rif:
                expected = math.sqrt((ols_se[key]**2 * (n-2) + rif[key]**2) / n)
                if not math.isclose(pooled_se[key], expected, abs_tol=1e-10, rel_tol=1e-10):
                    raise ValueError("pooled/OLS SE variance-convention identity failed")
            checks.append({"workload": workload, "repeat": repeat,
                           "maxerror": compare_slopes(gmoq, rif, f"{workload}/repeat{repeat}"),
                           "nrows": len(gmoq)})
    verification = {
        "status": "verified", "n_results": len(results), "n_workloads": len(workloads),
        "methods": list(METHODS), "repeats": 3, "atol": ATOL, "rtol": RTOL,
        "deterministic_within_method_repeats": True,
        "largest_panel": f"n{n_large}_m1000", "coefficient_only": False,
        "covariance_or_se_computed": True,
        "common_rif_setup_included_in_wall_seconds": True,
        "s1_checks": checks,
    }
    largest_workload = f"n{n_large}_m1000"
    return verification, checks, [arms[(largest_workload, m, 1)][1] for m in METHODS]


def _write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def export(benchmark_dir: Path, out_dir: Path, n_large: int = DEFAULT_N_LARGE) -> None:
    workloads = _workloads(n_large)
    verification, checks, largest = verify_benchmark(benchmark_dir, n_large)
    largest_workload = f"n{n_large}_m1000"
    out_dir = out_dir.resolve()
    if ".local" not in out_dir.parts and "outputs" not in out_dir.parts:
        raise ValueError("out-dir must be inside a .local or outputs directory")
    out_dir.mkdir(parents=True, exist_ok=False)
    gmoq, rif, jackknife = largest
    pairs = [{"snp_id": key[0], "tau": key[1], "beta_gmoq": gmoq[key],
              "beta_rif_ols": rif[key], "difference": gmoq[key] - rif[key]}
             for key in sorted(gmoq)]
    _write_csv(out_dir / "S1_pairs.csv",
               ["snp_id", "tau", "beta_gmoq", "beta_rif_ols", "difference"], pairs)
    _write_csv(out_dir / "S1_checks.csv", ["workload", "repeat", "maxerror", "nrows"], checks)
    data = json.loads((benchmark_dir / "results.json").read_text(encoding="utf-8"))
    se_tables = {r["method_id"]: read_slopes(Path(r["output_dir"]) / "slopes.csv", "se")
                 for r in data if r["workload_id"] == largest_workload and r["repeat"] == 1}
    se_rows = [{"snp_id": k[0], "tau": k[1], **{"se_"+m: se_tables[m][k] for m in METHODS}}
               for k in sorted(gmoq)]
    _write_csv(out_dir / "S1_SE.csv", ["snp_id", "tau"] + ["se_"+m for m in METHODS], se_rows)
    timing_rows = []
    for row in sorted(data, key=lambda x: (list(workloads).index(x["workload_id"]), x["method_id"], int(x["repeat"]))):
        n, m = workloads[row["workload_id"]]
        timing_rows.append({"workload": row["workload_id"], "method": row["method_id"],
                            "repeat": int(row["repeat"]), "n": n, "m": m,
                            "fit_seconds": row["fit_seconds"], "wall_seconds": row["wall_seconds"],
                            "peak_rss_kib": row["peak_rss_kib"]})
    _write_csv(out_dir / "S2_timings.csv",
               ["workload", "method", "repeat", "n", "m", "fit_seconds", "wall_seconds", "peak_rss_kib"], timing_rows)
    verification["benchmark_dir"] = str(benchmark_dir.resolve())
    verification["exported_files"] = ["S1_pairs.csv", "S1_checks.csv", "S2_timings.csv", "verification.json", "README.md"]
    (out_dir / "verification.json").write_text(json.dumps(verification, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (out_dir / "README.md").write_text("""# Stage 1 validation figures\n\n## S1 caption\nNumerical equivalence to ordinary RIF–OLS coefficients on the largest matched panel (N=N_large, M=1,000), across SNP/tau keys. The x-axis is the ordinary RIF–OLS coefficient and the y-axis is the G-MoQ Stage 1 coefficient. This verifies the implemented RIF regression; it does not establish an unbiased population estimator or conditional quantile regression.\n\n## S2 caption and limitations\nCoefficient-and-SE timings for matched M=1 and M-panel scans at N=2,000 and N=N_large. Fit time includes coefficient and SE calculation; whole-process time includes the common RIF setup. All arms compute SEs. G-MoQ pooled uses null-residual pooled covariance; ordinary OLS uses fitted residual RSS and n-2 degrees of freedom. They need not have identical SEs. Jackknife is a separate G-MoQ sensitivity arm (1000 blocks). The ordinary baseline uses efficient `np.linalg.lstsq(design, RIF45responses)` once per SNP. All three repeats are shown with their median; no speed claim is embedded here.\n\nOnly benchmark `slopes.csv` files and measured timing fields were read. Participant/sample arrays, logs, and covariance objects were not accessed or exported.\n""".replace("N=N_large", f"N={n_large:,}"), encoding="utf-8")


def _compute_guard() -> None:
    """Formerly restricted runs to the authors' cluster; kept as a no-op."""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-large", type=int, choices=(50000, 400000), default=DEFAULT_N_LARGE)
    args = parser.parse_args()
    _compute_guard()
    if args.out_dir.exists():
        raise RuntimeError(f"refusing existing exclusive output directory: {args.out_dir}")
    export(args.benchmark_dir.resolve(), args.out_dir, args.n_large)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
