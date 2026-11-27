"""Stage 1 workloads and fitting methods for the matched S2 speed benchmark.

What: "generate" simulates one workload (genotypes binomial(2, 0.3); phenotype 0.2 x centred
  genotype of the first variant plus N(0, 1) noise) and writes it to a private directory.
  "fit" fits one workload with one method and records its coefficients, standard errors
  and fit time: G-MoQ Stage 1 with the pooled covariance ("gmoq"), G-MoQ Stage 1 with the
  1,000-block jackknife ("gmoq_jackknife"), or ordinary RIF-OLS, one numpy lstsq regression
  per variant on the same RIFs ("ordinary_rif_ols"). The pooled arm uses the current
  production bin-prefix kernel. All arms retain the original common dense-RIF setup.
Paper: Supplementary Fig. S2; equality is checked only for benchmark integrity.
Input: none (generate); a generated workload directory (fit)
Output: genotype/phenotype arrays (kept on the cluster); per-method slopes, SEs and timing.json
Runs: on a BP compute node, called by run_stage1.py through benchmark.py.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import os
import socket
from pathlib import Path
import sys
import subprocess
import time
from typing import Any

import numpy as np

TAUS = np.linspace(0.02, 0.98, 45, dtype=np.float64)
MAX_SNPS = 1000


def balanced_block_ids(n_samples: int, n_blocks: int, seed: int) -> np.ndarray:
    """The frozen benchmark's balanced block assignment, unchanged."""
    if not 1 <= n_blocks <= n_samples:
        raise ValueError("require 1 <= n_blocks <= n_samples")
    counts = np.full(n_blocks, n_samples // n_blocks, dtype=np.int64)
    counts[:n_samples % n_blocks] += 1
    labels = np.repeat(np.arange(n_blocks, dtype=np.int32), counts)
    return np.random.RandomState(seed + 7919).permutation(labels).astype(np.int32, copy=False)


def load_pin(path: Path) -> dict:
    pin = json.loads(path.read_text())
    for key in ("commit", "stage1_root", "stage1_extension_sha256", "stage1_python"):
        if not isinstance(pin.get(key), str) or not pin[key]:
            raise ValueError(f"PIN.json missing {key}")
    if len(pin["commit"]) != 40 or len(pin["stage1_extension_sha256"]) != 64:
        raise ValueError("PIN.json must contain full commit and extension SHA256")
    return pin


def load_core(pin: dict):
    """Require the exact checkout, interpreter and compiled binary in PIN.json."""
    root = Path(pin["stage1_root"]).resolve()
    if Path(sys.executable).resolve() != Path(pin["stage1_python"]).resolve():
        raise ValueError("Run with stage1_python from PIN.json")
    head = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    dirty = subprocess.check_output(["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"], text=True).strip()
    if head != pin["commit"] or dirty:
        raise ValueError(f"Pinned checkout differs or has tracked changes: {root}")
    sys.path.insert(0, str(root / "stage1"))
    core = importlib.import_module("gmoq.core")
    if not Path(core.__file__).resolve().is_relative_to(root / "stage1") or not core.HAVE_CPP:
        raise ValueError("gmoq.core was not loaded from the pinned compiled checkout")
    extension = Path(core.cpp_ext.__file__).resolve()
    if not extension.is_relative_to(root / "stage1") or _sha256(extension) != pin["stage1_extension_sha256"]:
        raise ValueError("Stage 1 extension path or SHA256 differs from PIN.json")
    if core.cpp_ext.openmp_info()["max_threads"] != 1:
        raise ValueError("Pinned extension must run with OMP_NUM_THREADS=1")
    return core, {"root": str(root), "commit": head, "core": str(core.__file__),
                  "extension": str(extension), "extension_sha256": pin["stage1_extension_sha256"],
                  "python": sys.executable, "numpy": np.__version__}


def _array_pair(g: Any, rif: Any) -> tuple[np.ndarray, np.ndarray]:
    g = np.asarray(g, dtype=np.float64)
    rif = np.asarray(rif, dtype=np.float64)
    if g.ndim != 2 or rif.ndim != 2 or g.shape[0] != rif.shape[0]:
        raise ValueError(f"g/rif must be 2-D with equal rows, got {g.shape}, {rif.shape}")
    if rif.shape[1] != 45 or not np.all(np.isfinite(g)) or not np.all(np.isfinite(rif)):
        raise ValueError("g and rif must be finite and rif must have 45 columns")
    return g, rif


def generate(n: int, m: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Generate an independent Binomial(2, .3) panel and one simple phenotype."""
    n, m = int(n), int(m)
    if n < 2 or not 1 <= m <= MAX_SNPS:
        raise ValueError(f"require n >= 2 and 1 <= m <= {MAX_SNPS}")
    rng = np.random.default_rng(int(seed))
    g = rng.binomial(2, 0.3, size=(n, m)).astype(np.float64)
    y = 0.2 * (g[:, 0] - 0.6) + rng.normal(size=n)
    return g, np.asarray(y, dtype=np.float64)


def ordinary_rif_ols_with_se(g: np.ndarray, rif: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Fit all RIF responses and conventional homoskedastic OLS SEs."""
    g, rif = _array_pair(g, rif)
    if g.shape[0] <= 2:
        raise ValueError("ordinary OLS SEs require n > 2")
    ones = np.ones((g.shape[0], 1), dtype=np.float64)
    centered_rif = rif - rif.mean(axis=0)
    tss = np.sum(centered_rif ** 2, axis=0)
    slopes = np.empty((g.shape[1], 45), dtype=np.float64)
    ses = np.empty_like(slopes)
    for j in range(g.shape[1]):
        design = np.column_stack((ones, g[:, j]))
        slopes[j] = np.linalg.lstsq(design, rif, rcond=None)[0][1]
        centered_g = g[:, j] - g[:, j].mean()
        d = np.dot(centered_g, centered_g)
        rss = tss - slopes[j] ** 2 * d
        ses[j] = np.sqrt(np.maximum(rss / (g.shape[0] - 2) / d, 0.0))
    return slopes, ses


def ordinary_rif_ols(g: np.ndarray, rif: np.ndarray) -> np.ndarray:
    """Fit all 45 RIF responses jointly for each SNP using ordinary OLS."""
    return ordinary_rif_ols_with_se(g, rif)[0]


def _pooled_components(core: Any, g: np.ndarray, rif: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    g, rif = _array_pair(g, rif)
    ones = np.ones((g.shape[0], 1), dtype=np.float64)
    rif_resid, q_basis = core.residualize(rif, ones)
    betas, d_resid = core.compute_pooled_panel(
        g, rif_resid, q_basis, min_mac=5.0, min_var=1e-8, n_threads=1, panel_size=16)
    betas, d_resid = np.asarray(betas, dtype=np.float64), np.asarray(d_resid, dtype=np.float64)
    expected = (g.shape[1], 45)
    if betas.shape != expected or d_resid.shape != (g.shape[1],):
        raise ValueError(f"compiled pooled panel returned {betas.shape}, {d_resid.shape}, expected {expected}")
    return betas, d_resid, rif_resid


def fast_rif_ols(core: Any, g: np.ndarray, rif: np.ndarray) -> np.ndarray:
    """Use the pinned compiled pooled panel kernel for the same coefficient."""
    return _pooled_components(core, g, rif)[0]


def gmoq(core: Any, g: np.ndarray, rif: np.ndarray, y: np.ndarray,
         q_tau: np.ndarray, fhat: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Production pooled bin-prefix kernel; setup remains within fit+SE."""
    g, rif = _array_pair(g, rif)
    # All methods intentionally retain the original common dense-RIF setup.
    # The production pooled kernel itself needs y/q/f/Q, not the dense RIF.
    _, q_basis = core.residualize(rif, np.ones((len(y), 1)))
    plan_start = time.perf_counter()
    plan = core.make_pooled_bin_plan(y, TAUS, q_basis, q_tau=q_tau, fhat=fhat, n_threads=1)
    plan_build_seconds = time.perf_counter() - plan_start
    betas, d_resid = core.compute_pooled_binned(g, plan, min_mac=5.0, min_var=1e-8, n_threads=1)
    betas, d_resid = np.asarray(betas, dtype=np.float64), np.asarray(d_resid, dtype=np.float64)
    if betas.shape != (g.shape[1], 45) or d_resid.shape != (g.shape[1],):
        raise ValueError("Pooled bin-prefix kernel returned wrong dimensions")
    sigma0 = np.asarray(plan.sigma0(), dtype=np.float64)
    if sigma0.shape != (45, 45) or not np.all(np.isfinite(sigma0)):
        raise ValueError(f"pooled Sigma_0 returned {sigma0.shape}, expected (45, 45)")
    ses = np.sqrt(np.diag(sigma0)[None, :] / d_resid[:, None])
    return betas, ses, plan_build_seconds


def gmoq_jackknife(core: Any, g: np.ndarray, rif: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Pinned fast delete-block G-MoQ betas/SEs; discard the full covariance."""
    g, rif = _array_pair(g, rif)
    ones = np.ones((g.shape[0], 1), dtype=np.float64)
    rif_resid, q_basis = core.residualize(rif, ones)
    n_blocks = 1000
    block_ids = balanced_block_ids(g.shape[0], n_blocks, seed=940500)
    plan_start = time.perf_counter()
    plan = core.make_block_plan(rif_resid, q_basis, block_ids, n_blocks, n_threads=1)
    plan_build_seconds = time.perf_counter() - plan_start
    result = core.compute_block_finalize_panel_fast(
        g, rif_resid, q_basis, block_ids, n_blocks,
        min_mac=5.0, min_var=1e-8, n_threads=1, panel_width=16, plan=plan)
    if not isinstance(result, (tuple, list)) or len(result) != 3:
        raise ValueError("fast jackknife panel kernel must return betas, SEs, packed covariance")
    betas, ses, _covariance = (np.asarray(value, dtype=np.float64) for value in result)
    expected = (g.shape[1], 45)
    if betas.shape != expected or ses.shape != expected:
        raise ValueError(f"fast jackknife returned {betas.shape}, {ses.shape}, expected {expected}")
    return betas, ses, plan_build_seconds


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_write(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _private_exclusive(path: Path, allow_empty: bool = False) -> Path:
    path = path.expanduser().resolve()
    if str(path) in {"/", "/tmp", "/user", "/user/work"}:
        raise ValueError("Refusing unsafe output directory")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if not allow_empty or not path.is_dir() or any(path.iterdir()):
            raise ValueError(f"refusing nonempty or reused output: {path}")
    else:
        path.mkdir(exist_ok=False)
    return path


def _write_generated(n: int, m: int, seed: int, output_dir: Path) -> None:
    output_dir = _private_exclusive(output_dir)
    g, y = generate(n, m, seed)
    np.save(output_dir / "g.npy", g, allow_pickle=False)
    np.save(output_dir / "y.npy", y, allow_pickle=False)
    _json_write(output_dir / "metadata.json", {
        "schema_version": 1,
        "kind": "stage1_methods_input",
        "n": int(n), "m": int(m), "seed": int(seed),
        "generation": "independent Binomial(2, 0.3) SNPs; y=0.2*(g[:,0]-0.6)+Normal(0,1)",
        "scope": "shared-input",
        "files": {name: _sha256(output_dir / name) for name in ("g.npy", "y.npy")},
    })


def _read_input_manifest(input_dir: Path) -> tuple[dict[str, Any], str]:
    input_dir = input_dir.expanduser().resolve()
    manifest_path = input_dir / "metadata.json"
    try:
        metadata = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid input metadata: {manifest_path}") from exc
    if metadata.get("kind") != "stage1_methods_input" or metadata.get("schema_version") != 1:
        raise ValueError("input metadata is not a stage1_methods input manifest")
    n, m = int(metadata.get("n", 0)), int(metadata.get("m", 0))
    files = metadata.get("files")
    if n < 2 or not 1 <= m <= MAX_SNPS or not isinstance(files, dict):
        raise ValueError("input metadata has invalid dimensions or file checksums")
    for name in ("g.npy", "y.npy"):
        if not isinstance(files.get(name), str) or len(files[name]) != 64:
            raise ValueError(f"input metadata has no SHA256 for {name}")
    return metadata, _sha256(manifest_path)


def _load_input(input_dir: Path, metadata: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, dict[str, str]]:
    paths = {name: input_dir / name for name in ("g.npy", "y.npy")}
    g = np.asarray(np.load(paths["g.npy"], allow_pickle=False))
    y = np.asarray(np.load(paths["y.npy"], allow_pickle=False))
    if g.shape != (int(metadata["n"]), int(metadata["m"])) or y.shape != (int(metadata["n"]),):
        raise ValueError("input arrays do not match their manifest dimensions")
    if g.dtype != np.float64 or y.dtype != np.float64:
        raise ValueError("input arrays must be float64")
    checksums = {name: _sha256(path) for name, path in paths.items()}
    if checksums != {name: metadata["files"][name] for name in paths}:
        raise ValueError("input array checksum does not match metadata.json")
    return g, y, checksums


def _write_slopes(path: Path, slopes: np.ndarray, ses: np.ndarray | None = None) -> None:
    if ses is None:
        ses = np.full_like(slopes, np.nan)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("snp_id", "tau", "beta", "se"))
        for j in range(slopes.shape[0]):
            for t, tau in enumerate(TAUS):
                writer.writerow((f"snp_{j}", f"{tau:.17g}", f"{slopes[j, t]:.17g}", f"{ses[j, t]:.17g}"))


def _compute_fit(args: argparse.Namespace) -> None:
    if args.threads != 1:
        raise ValueError("Stage 1 timing arms require exactly one thread")
    pin = load_pin(args.pin_json)
    if args.package_commit != pin["commit"]:
        raise ValueError("--package-commit differs from PIN.json")
    core, package_provenance = load_core(pin)
    input_dir = args.input_dir.expanduser().resolve()
    input_metadata, input_metadata_sha = _read_input_manifest(input_dir)
    output_dir = _private_exclusive(args.output_dir, allow_empty=True)

    analysis_start = time.perf_counter()
    input_start = time.perf_counter()
    g, y, input_checksums = _load_input(input_dir, input_metadata)
    input_seconds = time.perf_counter() - input_start
    shared_start = time.perf_counter()
    rif, q_tau, fhat = core.compute_rif(
        y, TAUS, density_method="silverman_gaussian", density_floor=1e-8,
    )
    shared_setup_seconds = time.perf_counter() - shared_start
    fit_start = time.perf_counter()
    if args.method == "gmoq":
        slopes, ses, plan_build_seconds = gmoq(core, g, rif, y, q_tau, fhat)
    elif args.method == "ordinary_rif_ols":
        slopes, ses = ordinary_rif_ols_with_se(g, rif)
        plan_build_seconds = 0.0
    else:
        slopes, ses, plan_build_seconds = gmoq_jackknife(core, g, rif)
    fit_seconds = time.perf_counter() - fit_start
    if not np.all(np.isfinite(slopes)) or not np.all(np.isfinite(ses)):
        raise ValueError("Stage 1 fit returned non-finite betas or SEs")

    conventions = {
        "gmoq": "production pooled bin-prefix plan Sigma_0 / D",
        "ordinary_rif_ols": "conventional homoskedastic OLS: RSS/(n-2)/centered genotype sum of squares",
        "gmoq_jackknife": "package delete-one-block jackknife covariance; diagonal SEs from 1000 balanced blocks",
    }
    write_start = time.perf_counter()
    _write_slopes(output_dir / "slopes.csv", slopes, ses)
    _json_write(output_dir / "metadata.json", {
        "schema_version": 1, "kind": "stage1_methods_fit",
        "method": args.method, "scope": "beta+SE", "inference_scope": "beta+SE",
        "n": int(g.shape[0]), "m": int(g.shape[1]), "n_taus": 45, "tau_count": 45,
        "covariance_computed": args.method != "ordinary_rif_ols", "standard_errors_computed": True,
        "computed_se": True, "covariance_convention": conventions[args.method],
        "full_jackknife_covariance": (args.method == "gmoq_jackknife"),
        "full_jackknife_covariance_note": (
            "full package packed covariance was computed and discarded; only diagonal SEs were written"
            if args.method == "gmoq_jackknife" else None
        ),
        "plan_build_seconds": plan_build_seconds,
        "taus": {"start": 0.02, "stop": 0.98, "count": 45},
        "density_method": "silverman_gaussian", "density_floor": 1e-8,
        "threads": 1, "panel_size": 16,
        "package_commit": pin["commit"],
        "pin_json": str(args.pin_json.resolve()),
        "pin_sha256": _sha256(args.pin_json),
        "kernel": "compute_pooled_binned" if args.method == "gmoq" else args.method,
        "timing_boundary": "common dense RIF outside fit; residualisation, plan, beta and SE inside fit",
        "package_provenance": package_provenance,
        "input_manifest": {
            "metadata_sha256": input_metadata_sha,
            "files": input_checksums,
            "n": int(g.shape[0]), "m": int(g.shape[1]),
        },
    })
    write_seconds = time.perf_counter() - write_start
    total_analysis_seconds = time.perf_counter() - analysis_start
    _json_write(output_dir / "timing.json", {
        "input_seconds": input_seconds,
        "shared_setup_seconds": shared_setup_seconds,
        "fit_seconds": fit_seconds,
        "write_seconds": write_seconds,
        "total_analysis_seconds": total_analysis_seconds,
        "scope": "beta+SE",
        "inference_scope": "beta+SE",
        "method": args.method,
        "plan_build_seconds": plan_build_seconds,
    })


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    gen = sub.add_parser("generate")
    gen.add_argument("--n", type=int, required=True)
    gen.add_argument("--m", type=int, required=True)
    gen.add_argument("--seed", type=int, required=True)
    gen.add_argument("--output-dir", type=Path, required=True)
    fit = sub.add_parser("fit")
    fit.add_argument("--input-dir", type=Path, required=True)
    fit.add_argument("--method", choices=("gmoq", "ordinary_rif_ols", "gmoq_jackknife"), required=True)
    fit.add_argument("--output-dir", type=Path, required=True)
    fit.add_argument("--threads", type=int, default=1)
    fit.add_argument("--threads1", action="store_true", help=argparse.SUPPRESS)
    fit.add_argument("--package-commit", required=True)
    fit.add_argument("--pin-json", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    if not socket.gethostname().split(".")[0].startswith("bp1-compute") or not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Benchmark execution requires a BluePebble compute node under Slurm")
    args = _parser().parse_args(argv)
    if args.mode == "generate":
        _write_generated(args.n, args.m, args.seed, args.output_dir)
    else:
        _compute_fit(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

__all__ = [
    "TAUS", "generate", "ordinary_rif_ols", "ordinary_rif_ols_with_se",
    "fast_rif_ols", "gmoq", "gmoq_jackknife", "balanced_block_ids",
]
