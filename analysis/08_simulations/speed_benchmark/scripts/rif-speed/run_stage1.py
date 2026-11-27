#!/usr/bin/env python3
"""Matched Stage 1 speed benchmark for Supplementary Figure S2.

The four original workloads, three methods, seeds and three fresh-process repeats
are frozen. Equality checks validate benchmark integrity only. No scientific
recovery, calibration or power simulations are performed here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import socket
import subprocess
import sys

from benchmark import run_benchmark

HERE = Path(__file__).resolve().parent
METHODS = ("gmoq", "ordinary_rif_ols", "gmoq_jackknife")
WORKLOADS = ((2000, 1, 940001), (2000, 256, 940002),
             (400000, 1, 940003), (400000, 1000, 940004))
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
               "RAYON_NUM_THREADS")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")


def write_csv(path, rows):
    with Path(path).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def export(results, specs, out):
    """Check matched inputs and deterministic coefficients before exporting S2."""
    import numpy as np
    expected = {(f"n{n}_m{m}", method, rep) for n, m, _ in specs
                for method in METHODS for rep in (1, 2, 3)}
    arms = {}
    for row in results:
        key = row["workload_id"], row["method_id"], row["repeat"]
        if key in arms or row["status"] != "success" or row["return_code"] != 0 or row["child_return_code"] != 0:
            raise ValueError(f"Duplicate or failed benchmark arm: {key}")
        path = Path(row["output_dir"])
        metadata = json.loads((path / "metadata.json").read_text())
        with (path / "slopes.csv").open(newline="") as f:
            slopes = list(csv.DictReader(f))
        keys = [(r["snp_id"], float(r["tau"])) for r in slopes]
        values = np.array([[float(r["beta"]), float(r["se"])] for r in slopes])
        n, m = metadata["n"], metadata["m"]
        grid = [(f"snp_{j}", float(t)) for j in range(m) for t in np.linspace(.02, .98, 45)]
        if keys != grid or not np.isfinite(values).all() or not (values[:, 1] > 0).all():
            raise ValueError(f"Invalid slopes/SEs/keys: {key}")
        if metadata["method"] != row["method_id"] or metadata["package_commit"] != row["package_pin"]:
            raise ValueError(f"Metadata differs: {key}")
        if row["workload_id"] != f"n{n}_m{m}":
            raise ValueError(f"Wrong workload dimensions: {key}")
        for name in ("fit_seconds", "wall_seconds", "peak_rss_kib"):
            if row[name] is None or not np.isfinite(row[name]) or row[name] < 0:
                raise ValueError(f"Invalid {name}: {key}")
        arms[key] = (row, metadata, values)
    if set(arms) != expected:
        raise ValueError("Missing or extra workload/method/repeat arms")
    checks, timing, means = [], [], []
    for n, m, seed in specs:
        wid = f"n{n}_m{m}"
        ref = arms[wid, "ordinary_rif_ols", 1]
        for method in METHODS:
            first = arms[wid, method, 1]
            for rep in (1, 2, 3):
                row, meta, values = arms[wid, method, rep]
                if meta["input_manifest"] != ref[1]["input_manifest"]:
                    raise ValueError(f"Inputs differ: {wid}/{method}/{rep}")
                if not np.array_equal(values, first[2]):
                    raise ValueError(f"Within-method repeats differ: {wid}/{method}")
                diff = float(np.max(np.abs(values[:, 0] - ref[2][:, 0])))
                if not np.allclose(values[:, 0], ref[2][:, 0], atol=1e-10, rtol=1e-10):
                    raise ValueError(f"Coefficient integrity failure: {wid}/{method}, max difference {diff}")
                checks.append({"workload": wid, "method": method, "repeat": rep,
                               "n_coefficients": m * 45, "max_beta_difference_from_ols": diff})
                timing.append({"workload": wid, "method": method, "repeat": rep, "n": n, "m": m,
                               "fit_seconds": row["fit_seconds"], "wall_seconds": row["wall_seconds"],
                               "peak_rss_kib": row["peak_rss_kib"], "seed": seed,
                               "input_seconds": row["input_seconds"],
                               "shared_setup_seconds": row["shared_setup_seconds"],
                               "write_seconds": row["write_seconds"], "order": row["order"]})
            for metric, field, divisor in (("fit + SE seconds", "fit_seconds", 1),
                                            ("whole-process seconds", "wall_seconds", 1),
                                            ("peak RSS MiB", "peak_rss_kib", 1024)):
                means.append({"workload": wid, "method": method, "n": n, "m": m,
                              "metric": metric, "mean_value": float(np.mean([
                                  arms[wid, method, rep][0][field] / divisor for rep in (1, 2, 3)])), "count": 3})
        pooled = arms[wid, "gmoq", 1][2]
        conventional = ref[2]
        expected_se = np.sqrt((conventional[:, 1] ** 2 * (n - 2) + conventional[:, 0] ** 2) / n)
        if not np.allclose(pooled[:, 1], expected_se, atol=1e-10, rtol=1e-10):
            raise ValueError(f"Pooled/ordinary OLS SE convention identity failure: {wid}")
    tables = out / "tables"
    tables.mkdir()
    write_csv(tables / "S2_timings.csv", timing)
    write_csv(tables / "S2_mean_timings.csv", means)
    write_csv(tables / "S2_integrity_checks.csv", checks)
    write_json(tables / "verification.json", {
        "status": "verified", "scope": "speed benchmark integrity only",
        "n_results": len(results), "methods": list(METHODS), "repeats": 3,
        "atol": 1e-10, "rtol": 1e-10, "input_checksums_equal": True,
        "within_method_repeats_equal": True, "pooled_ols_se_convention_identity": True,
        "max_beta_difference_from_ols": max(r["max_beta_difference_from_ols"] for r in checks),
        "common_dense_rif_setup_in_whole_process": True,
        "fit_includes_residualisation_plan_coefficients_se": True,
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pin-json", "--pin", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--preset", choices=("smoke", "full"), default="full")
    ap.add_argument("--workload", choices=tuple(f"n{n}_m{m}" for n, m, _ in WORKLOADS),
                    help="Run just one original size; all three methods and three repeats are retained")
    ap.add_argument("--timeout-seconds", type=float, default=1800)
    args = ap.parse_args()
    if not socket.gethostname().split(".")[0].startswith("bp1-compute") or not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run on a BluePebble compute node under Slurm")
    # Set numerical-library limits before importing NumPy, including generation/export.
    for name in THREAD_VARS:
        os.environ[name] = "1"
    from stage1_methods import load_pin
    pin_path = args.pin_json.expanduser().resolve()
    pin = load_pin(pin_path)
    if Path(sys.executable).resolve() != Path(pin["stage1_python"]).resolve():
        raise ValueError("Run with stage1_python from PIN.json")
    specs = list(WORKLOADS if args.preset == "full" else WORKLOADS[:2])
    if args.workload:
        specs = [s for s in WORKLOADS if f"n{s[0]}_m{s[1]}" == args.workload]
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=False)
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    os.environ["VALIDATION_COMPUTE_HOST"] = "1"
    topo = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
    write_json(out / "provenance.json", {
        "pin_json": str(pin_path), "pin_sha256": sha256(pin_path), "pin": pin,
        "host": socket.gethostname(), "slurm_job": os.environ["SLURM_JOB_ID"], "cpu": cpu,
        "core_id": (topo / "core_id").read_text().strip(),
        "socket_id": (topo / "physical_package_id").read_text().strip(),
        "platform": platform.platform(), "source_sha256": {p.name: sha256(p) for p in HERE.glob("*.py")},
        "source": "gmoq-paper-code analysis/08_simulations/code/validation/figures run_stage1/stage1_methods",
        "workloads": [{"n": n, "m": m, "seed": seed} for n, m, seed in specs],
        "scope": "Supplementary Figure S2 speed only; equality is an integrity check",
    })
    methods_script = HERE / "stage1_methods.py"
    workloads = []
    for n, m, seed in specs:
        wid = f"n{n}_m{m}"
        data = out / "inputs" / wid
        subprocess.run([sys.executable, str(methods_script), "generate", "--n", str(n), "--m", str(m),
                        "--seed", str(seed), "--output-dir", str(data)], check=True)
        workloads.append({"workload_id": wid, "execution_mode": "isolated_single_snp" if m == 1 else "amortized_scan",
                          "placeholders": {"input_path": str(data)}, "timing_file": "timing.json"})
    methods = [{"method_id": method, "command": [
        sys.executable, str(methods_script), "fit", "--input-dir", "{input_path}",
        "--method", method, "--output-dir", "{output_dir}", "--threads", "{threads}",
        "--package-commit", "{package_pin}", "--pin-json", str(pin_path)]} for method in METHODS]
    manifest = {"schema_version": "validation-benchmark-manifest-1", "experiment_id": out.name,
                "package_pin": pin["commit"], "output_root": str(out / "benchmark"),
                "compute_host_only": True, "physical_cores": 1, "physical_cpu_ids": [cpu], "threads": 1,
                "timeout_seconds": args.timeout_seconds, "interleave_seed": 940100,
                "warmup_count": 0, "repeats": 3, "workloads": workloads, "methods": methods}
    write_json(out / "command_manifest.json", manifest)
    results = run_benchmark(manifest)
    export(results, specs, out)
    (out / "COMPLETE").write_text("S2 speed benchmark and integrity checks complete.\n")
    print(json.dumps({"status": "complete", "out": str(out), "arms": len(results)}), flush=True)


if __name__ == "__main__":
    main()
