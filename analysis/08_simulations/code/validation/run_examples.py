#!/usr/bin/env python3
"""Shared simulation helpers: Stage 1 and Stage 2 through the G-MoQ packages.

What: the functions every simulation driver uses. load_environment_manifest and
  load_pinned_core load the compiled Stage 1 package from the environment made by
  bootstrap.sh; _stage1_example draws one data set (truth.draw), computes the RIFs,
  residualises them and fits Stage 1 with the pooled and the 1,000-block jackknife
  covariance, all through the package's own routines; _run_stage2 runs stage2.R on the
  result. Run as a script, it also runs the eight small validation examples of
  config/examples.json.
Paper: Supplementary Figs. S1-S6; Supplementary Tables S1, S2 (shared by studies 1-7)
Input: none (simulated data); the G-MoQ packages
Output: none when imported; small JSON summaries when run as a script
Runs: on a cluster. The package location (PINNED_PACKAGE_ROOT) and the output paths of the
  example profiles are the authors' cluster paths; set GMOQ_PACKAGE_ROOT to a checkout of the
  G-MoQ package at the Stage 1 commit named in config.yaml to override the first.

Simulated participant-level arrays never leave process memory; only one-variant JSON
summaries are written.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
from typing import Any, Mapping, Sequence

import numpy as np

# User setting: a checkout of the G-MoQ package at the Stage 1 commit below. The default is
# the authors' cluster path.
PINNED_PACKAGE_ROOT = Path(os.environ.get("GMOQ_PACKAGE_ROOT", "/user/work/fh6520/G-MoQ-validation-9c47e7b"))
PINNED_PACKAGE_COMMIT = "9c47e7ba56370e708810fff8c2c9ad3b723cbccf"
# The combined repository labels Stage 2 as 0.4.0 while its separately
# installable Stage 1 distribution still reports 0.1.0.  The Git SHA is the
# primary pin; keep both package-version facts explicit rather than rejecting
# the actual pinned Stage 1 metadata.
PINNED_PACKAGE_VERSION = "0.4.0"
PINNED_STAGE1_VERSION = "0.1.0"
ABS_REL_TOL = 1e-10
MAX_HANDOFF_ARRAY_ELEMENTS = 4096
EXPECTED_CASE_IDS = [
    "normal_null",
    "normal_location",
    "normal_dispersion",
    "normal_both",
    "mixture_null",
    "mixture_membership",
    "mixture_left",
    "mixture_right",
]
CASE_FIELDS = (
    "family", "maf", "a", "c", "gamma", "delta1", "delta2",
    "pi0", "mu1", "mu2", "sd1", "sd2",
)
CASE_SPECS = {
    "normal_null": ("normal", .30, 0., 0., 0., 0., 0., .5, -2., 2., 1., 1.),
    "normal_location": ("normal", .30, .05, 0., 0., 0., 0., .5, -2., 2., 1., 1.),
    "normal_dispersion": ("normal", .30, 0., .05, 0., 0., 0., .5, -2., 2., 1., 1.),
    "normal_both": ("normal", .30, .05, .05, 0., 0., 0., .5, -2., 2., 1., 1.),
    "mixture_null": ("mixture", .30, 0., 0., 0., 0., 0., .5, -2., 2., 1., 1.),
    "mixture_membership": ("mixture", .30, 0., 0., .10, 0., 0., .5, -2., 2., 1., 1.),
    "mixture_left": ("mixture", .30, 0., 0., 0., .10, 0., .5, -2., 2., 1., 1.),
    "mixture_right": ("mixture", .30, 0., 0., 0., 0., .10, .5, -2., 2., 1., 1.),
}


def _json_safe(value: Any) -> Any:
    """Convert small NumPy-containing objects to strict JSON values."""
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        if len(value) > MAX_HANDOFF_ARRAY_ELEMENTS:
            raise ValueError(
                "JSON hand-off contains an unexpectedly large array; participant "
                "arrays must remain in process memory/cluster scratch"
            )
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        if value.size > MAX_HANDOFF_ARRAY_ELEMENTS:
            raise ValueError(
                "JSON hand-off contains an unexpectedly large NumPy array; participant "
                "arrays must remain in process memory/cluster scratch"
            )
        return _json_safe(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError(f"non-finite value cannot enter JSON: {value!r}")
        return value
    return value


def _write_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(_json_safe(payload), indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def row_major_upper(matrix: np.ndarray) -> np.ndarray:
    """Return the packed upper triangle in Stage 2's row-major order."""
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"expected a square matrix, got {matrix.shape}")
    rows, cols = np.triu_indices(matrix.shape[0])
    return matrix[rows, cols]


def _git_head(package_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(package_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            f"cannot verify pinned package source at {package_root}; git rev-parse failed"
        ) from exc
    head = result.stdout.strip()
    if not head:
        raise RuntimeError(f"git returned an empty HEAD for package root {package_root}")
    return head


def _git_tracked_status(package_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(package_root), "status", "--porcelain=v1",
             "--untracked-files=no"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            f"cannot verify tracked source state at {package_root}"
        ) from exc
    return result.stdout.strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise RuntimeError(f"cannot read manifest-pinned file: {path}") from exc
    return digest.hexdigest()


def _manifest_string(manifest: Mapping[str, Any], *keys: str) -> str:
    value: Any = manifest
    for key in keys:
        if not isinstance(value, Mapping) or key not in value:
            raise RuntimeError("environment manifest is missing " + ".".join(keys))
        value = value[key]
    if not isinstance(value, str) or not value:
        raise RuntimeError("environment manifest has an invalid " + ".".join(keys))
    return value


def load_environment_manifest(path: Path, package_root: Path) -> dict[str, Any]:
    """Load and verify the compute-node environment made by bootstrap.sh."""
    path = path.resolve()
    if not path.is_file():
        raise RuntimeError(f"environment manifest is missing: {path}")
    try:
        manifest = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"cannot read environment manifest: {path}") from exc
    if not isinstance(manifest, Mapping):
        raise RuntimeError("environment manifest must contain a JSON object")

    package = manifest.get("package")
    if not isinstance(package, Mapping):
        raise RuntimeError("environment manifest is missing package provenance")
    if package.get("root") != str(package_root.resolve()):
        raise RuntimeError("environment manifest package root does not match the config")
    if package.get("commit") != PINNED_PACKAGE_COMMIT or package.get("head_verified") != PINNED_PACKAGE_COMMIT:
        raise RuntimeError("environment manifest does not pin the required package commit")
    if package.get("tracked_clean_before_build") is not True or package.get("tracked_clean_after_build") is not True:
        raise RuntimeError("environment manifest does not certify clean tracked package source")

    python_executable = Path(_manifest_string(manifest, "python", "executable")).resolve()
    if Path(sys.executable).resolve() != python_executable:
        raise RuntimeError(
            f"running Python is {Path(sys.executable).resolve()}, but the environment "
            f"manifest requires {python_executable}"
        )
    if not python_executable.is_file() or not os.access(python_executable, os.X_OK):
        raise RuntimeError(f"manifest Python executable is not executable: {python_executable}")

    extension_path = Path(_manifest_string(manifest, "stage1_extension", "path")).resolve()
    python_extension_path = Path(_manifest_string(manifest, "python", "extension", "path")).resolve()
    if extension_path != python_extension_path:
        raise RuntimeError("environment manifest has inconsistent Stage 1 extension paths")
    extension_sha = _manifest_string(manifest, "stage1_extension", "sha256")
    if not extension_path.is_file() or _sha256(extension_path) != extension_sha:
        raise RuntimeError("manifest-pinned Stage 1 extension path or SHA256 does not match")

    rscript = Path(_manifest_string(manifest, "r", "executables", "Rscript")).resolve()
    if not rscript.is_file() or not os.access(rscript, os.X_OK):
        raise RuntimeError(f"manifest Rscript is not executable: {rscript}")
    rlib = Path(_manifest_string(manifest, "stage2", "isolated_library")).resolve()
    r_gmoq = Path(_manifest_string(manifest, "r", "gmoq", "isolated_library")).resolve()
    if rlib != r_gmoq:
        raise RuntimeError("environment manifest has inconsistent private R library paths")
    installed = Path(_manifest_string(manifest, "r", "gmoq", "installed_path")).resolve()
    try:
        installed.relative_to(rlib)
    except ValueError as exc:
        raise RuntimeError("gmoq is not installed inside the manifest private R library") from exc
    if not rlib.is_dir():
        raise RuntimeError(f"manifest private R library is missing: {rlib}")

    return {
        "path": str(path),
        "rscript": str(rscript),
        "r_lib": str(rlib),
        "extension_path": str(extension_path),
        "extension_sha256": extension_sha,
        "manifest": manifest,
    }


def load_pinned_core(package_root: Path, environment: Mapping[str, Any] | None = None):
    """Import Stage 1 from the pinned checkout and reject every fallback path."""
    package_root = package_root.resolve()
    if package_root != PINNED_PACKAGE_ROOT:
        raise RuntimeError(
            f"the compact contract pins package_root={PINNED_PACKAGE_ROOT}, got {package_root}"
        )
    head = _git_head(package_root)
    if head != PINNED_PACKAGE_COMMIT:
        raise RuntimeError(
            f"refusing unpinned Stage 1 source: expected {PINNED_PACKAGE_COMMIT}, got {head}"
        )
    tracked_status = _git_tracked_status(package_root)
    if tracked_status:
        raise RuntimeError(
            "refusing Stage 1 source with tracked modifications:\n" + tracked_status
        )

    stage1_root = package_root / "stage1"
    if not stage1_root.is_dir():
        raise RuntimeError(f"pinned Stage 1 source directory is missing: {stage1_root}")
    sys.path.insert(0, str(stage1_root))
    # Do not silently reuse an installed package imported before this check.
    for name in list(sys.modules):
        if name == "gmoq" or name.startswith("gmoq."):
            del sys.modules[name]
    try:
        package = importlib.import_module("gmoq")
        core = importlib.import_module("gmoq.core")
        from gmoq.pooled_cov import compute_sigma0
    except Exception as exc:
        raise RuntimeError("failed to import Stage 1 from the pinned checkout") from exc

    package_file = Path(package.__file__).resolve()
    core_file = Path(core.__file__).resolve()
    try:
        package_file.relative_to(stage1_root.resolve())
        core_file.relative_to(stage1_root.resolve())
    except ValueError as exc:
        raise RuntimeError(
            "Stage 1 import did not resolve inside the pinned package checkout"
        ) from exc
    if getattr(package, "__version__", PINNED_STAGE1_VERSION) != PINNED_STAGE1_VERSION:
        raise RuntimeError(
            f"unexpected gmoq version: {getattr(package, '__version__', None)!r}"
        )
    if not bool(getattr(core, "HAVE_CPP", False)):
        raise RuntimeError(
            "the pinned Stage 1 compiled extension is unavailable; examples refuse "
            "the NumPy fallback"
        )
    if bool(getattr(core, "ALLOW_NUMPY_FALLBACK", False)):
        raise RuntimeError("NumPy fallback is enabled; compact examples require compiled code")
    cpp_ext = getattr(core, "cpp_ext", None)
    required_symbols = (
        "stage1_pooled_panel_static",
        "stage1_finalize_block_scores_panel_fast",
        "BlockPlan",
    )
    missing = [name for name in required_symbols if not hasattr(cpp_ext, name)]
    if missing:
        raise RuntimeError(
            "the pinned compiled extension lacks the required fast-panel symbols: "
            + ", ".join(missing)
        )
    cpp_file = Path(cpp_ext.__file__).resolve()
    try:
        cpp_file.relative_to(stage1_root.resolve())
    except ValueError as exc:
        raise RuntimeError(
            f"compiled extension resolved outside the pinned checkout: {cpp_file}"
        ) from exc
    extension_sha256 = _sha256(cpp_file)
    if environment is not None:
        expected_extension = Path(str(environment["extension_path"])).resolve()
        expected_sha256 = str(environment["extension_sha256"])
        if cpp_file != expected_extension or extension_sha256 != expected_sha256:
            raise RuntimeError(
                "loaded Stage 1 extension does not match the environment manifest "
                f"(path={cpp_file}, sha256={extension_sha256})"
            )
    return core, {
        "package_root": str(package_root),
        "package_commit": head,
        "package_version": PINNED_PACKAGE_VERSION,
        "stage1_distribution_version": getattr(package, "__version__", PINNED_STAGE1_VERSION),
        "python_package": str(package_file),
        "core_module": str(core_file),
        "compiled_extension": str(cpp_file),
        "compiled_extension_sha256": extension_sha256,
        "compiled_extension_required": True,
        "fast_panel_symbols": list(required_symbols),
        "numpy_fallback_allowed": False,
    }, compute_sigma0


def _load_truth_module(repo_root: Path):
    truth_path = repo_root / "validation" / "truth.py"
    if not truth_path.is_file():
        raise RuntimeError(f"required truth module is missing: {truth_path}")
    spec = importlib.util.spec_from_file_location("validation_truth", truth_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load truth module: {truth_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    required = ("draw", "population_truth", "fit_baseline", "independent_weights")
    missing = [name for name in required if not hasattr(module, name)]
    if missing:
        raise RuntimeError(f"truth.py does not expose the design API: {missing}")
    return module


def _assert_close(name: str, observed: np.ndarray, expected: np.ndarray) -> float:
    observed = np.asarray(observed, dtype=np.float64)
    expected = np.asarray(expected, dtype=np.float64)
    if observed.shape != expected.shape:
        raise RuntimeError(f"{name} shape mismatch: {observed.shape} versus {expected.shape}")
    if not np.all(np.isfinite(observed)) or not np.all(np.isfinite(expected)):
        raise RuntimeError(f"{name} contains a non-finite value")
    diff = np.abs(observed - expected)
    limit = ABS_REL_TOL + ABS_REL_TOL * np.maximum(1.0, np.abs(expected))
    if np.any(diff > limit):
        i = np.unravel_index(int(np.argmax(diff)), diff.shape)
        raise RuntimeError(
            f"{name} failed the independent RIF-OLS check at {i}: "
            f"difference={diff[i]:.6g}, limit={limit[i]:.6g}"
        )
    return float(np.max(diff))


def _finite_vector(name: str, value: Any, length: int | None = None) -> np.ndarray:
    if isinstance(value, Mapping):
        value = list(value.values())
    array = np.asarray(value, dtype=np.float64)
    if array.ndim != 1 or (length is not None and array.size != length):
        raise RuntimeError(f"{name} must be a vector of length {length}, got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"{name} contains a non-finite value")
    return array


def _validate_baseline(name: str, baseline: Mapping[str, Any], family: str) -> None:
    required = ("mu", "sd") if family == "normal" else ("p1", "mu1", "sd1", "mu2", "sd2")
    missing = [key for key in required if key not in baseline]
    if missing:
        raise RuntimeError(f"{name} is missing baseline fields: {missing}")
    for key in required:
        value = float(baseline[key])
        if not np.isfinite(value):
            raise RuntimeError(f"{name}.{key} is non-finite")
    if family == "normal" and float(baseline["sd"]) <= 0:
        raise RuntimeError(f"{name}.sd must be positive")
    if family == "mixture":
        if not 0 < float(baseline["p1"]) < 1:
            raise RuntimeError(f"{name}.p1 must lie strictly between zero and one")
        if float(baseline["sd1"]) <= 0 or float(baseline["sd2"]) <= 0:
            raise RuntimeError(f"{name} component standard deviations must be positive")
        if not float(baseline["mu1"]) < float(baseline["mu2"]):
            raise RuntimeError(f"{name} component means must be ordered lower then upper")


def validate_config(config: Mapping[str, Any]) -> None:
    """Validate one of the three fixed profiles before creating any output."""
    if config.get("schema_version") != 1:
        raise ValueError("validation config schema_version must be 1")
    if config.get("profile", "examples_v1") not in (
        "examples_v1", "normal_pilot_v1", "normal_replication_v1"
    ):
        raise ValueError("validation config has an unsupported profile")
    profile = config.get("profile", "examples_v1")
    pilot = profile == "normal_pilot_v1"
    replication = profile == "normal_replication_v1"
    expected = {
        "run_id": (
            "normal-replication-20260917-01" if replication else
            "normal-pilot-20260917-01" if pilot else
            "validation-examples-20260917-01"
        ),
        "output_root": (
            "/user/work/fh6520/Gmoq_paper-validation-20260917/.local/normal-replication-20260917-01"
            if replication else
            "/user/work/fh6520/Gmoq_paper-validation-20260917/.local/normal-pilot-20260917-01"
            if pilot else
            "/user/work/fh6520/Gmoq_paper-validation-20260917/.local/validation-examples-20260917-01"
        ),
        "n": 10000 if (pilot or replication) else 2000,
        "replicates": 1000 if replication else 20 if pilot else 2,
        "seed": 90000000 if replication else 20261017 if pilot else 20260917,
    }
    for key, value in expected.items():
        if config.get(key) != value:
            raise ValueError(f"{profile} requires {key}={value!r}")
    if config.get("package_root") != str(PINNED_PACKAGE_ROOT):
        raise ValueError("validation config has the wrong pinned package_root")
    if config.get("package_commit") != PINNED_PACKAGE_COMMIT:
        raise ValueError("validation config has the wrong pinned package_commit")
    if config.get("package_version") != PINNED_PACKAGE_VERSION:
        raise ValueError("validation config has the wrong package_version")
    if config.get("stage1_distribution_version") != PINNED_STAGE1_VERSION:
        raise ValueError("validation config has the wrong Stage 1 distribution version")
    if config.get("maf") != 0.3 or config.get("n_blocks") != 1000:
        raise ValueError("profiles require MAF 0.30 and 1000 blocks")
    if config.get("physical_cores") != 1 or config.get("threads") != 1 or config.get("panel_size") != 16:
        raise ValueError("profiles require one physical core/thread and panel_size 16")
    if config.get("stage2_script") != "validation/stage2.R":
        raise ValueError("profiles require the pinned Stage 2 adapter")
    if config.get("density_method") != "silverman_gaussian" or config.get("density_floor") != 1e-8:
        raise ValueError("profiles require the pinned Silverman Gaussian RIF KDE")
    if config.get("min_mac") != 5.0 or config.get("min_var") != 1e-8:
        raise ValueError("profiles require the fixed Stage 1 numerical floors")
    taus = np.asarray(config.get("taus", []), dtype=np.float64)
    if taus.shape != (45,) or not np.allclose(taus, np.linspace(0.02, 0.98, 45), atol=2e-12, rtol=0):
        raise ValueError("profiles require 45 equally spaced taus in [0.02, 0.98]")
    expected_ids = (
        ["normal_null", "normal_location"]
        if (pilot or replication) else EXPECTED_CASE_IDS
    )
    cases = config.get("cases")
    if not isinstance(cases, list) or [case.get("id") for case in cases] != expected_ids:
        raise ValueError(f"{profile} cases must be exactly {expected_ids}")
    for case in cases:
        if set(case) != set(("id",) + CASE_FIELDS):
            raise ValueError(f"{profile} case {case.get('id')!r} has the wrong schema")
        actual = tuple(case.get(key) for key in CASE_FIELDS)
        if actual != CASE_SPECS[case["id"]]:
            raise ValueError(f"{profile} case {case['id']} has non-contract effect or baseline values")
    if config["n"] < config["n_blocks"]:
        raise ValueError("profiles require n >= n_blocks")


def balanced_block_ids(n_samples: int, n_blocks: int, seed: int) -> np.ndarray:
    """Assign every requested block a near-equal count, then permute labels."""
    n_samples = int(n_samples)
    n_blocks = int(n_blocks)
    if n_samples < 1 or n_blocks < 1 or n_blocks > n_samples:
        raise ValueError(
            f"balanced blocks require 1 <= n_blocks <= n_samples, got "
            f"{n_blocks} and {n_samples}"
        )
    counts = np.full(n_blocks, n_samples // n_blocks, dtype=np.int64)
    counts[: n_samples % n_blocks] += 1
    labels = np.repeat(np.arange(n_blocks, dtype=np.int32), counts)
    return np.random.RandomState(int(seed) + 7919).permutation(labels).astype(np.int32, copy=False)


def _stage1_example(core, truth: Any, case: Mapping[str, Any], n: int,
                     seed: int, taus: np.ndarray, config: Mapping[str, Any],
                     compute_sigma0: Any, population: Mapping[str, Any] | None = None) -> dict:
    drawn = truth.draw(dict(case), n, seed)
    if not isinstance(drawn, Mapping) or "g" not in drawn or "y" not in drawn:
        raise RuntimeError("truth.draw must return a mapping with g and y")
    g = np.asarray(drawn["g"], dtype=np.float64)
    y = _finite_vector("drawn phenotype", drawn["y"], n)
    if g.ndim == 1:
        g = g[:, None]
    if g.shape != (n, 1) or not np.all(np.isfinite(g)):
        raise RuntimeError(f"truth.draw must return one finite SNP, got {g.shape}")

    family = str(case["family"])
    if population is None:
        population = truth.population_truth(dict(case), taus)
    if not isinstance(population, Mapping):
        raise RuntimeError("truth.population_truth must return a mapping")
    required_truth = ("q_tau", "f_tau", "b_tau", "baseline", "theta", "generating", "diagnostics")
    missing_truth = [key for key in required_truth if key not in population]
    if missing_truth:
        raise RuntimeError(f"population_truth is missing required fields: {missing_truth}")
    q_pop = _finite_vector("population q_tau", population["q_tau"], len(taus))
    _finite_vector("population f_tau", population["f_tau"], len(taus))
    _finite_vector("population b_tau", population["b_tau"], len(taus))
    _finite_vector("population theta", population["theta"])
    baseline_pop = population["baseline"]
    if not isinstance(baseline_pop, Mapping):
        raise RuntimeError("population baseline must be a mapping")
    _validate_baseline("population baseline", baseline_pop, family)
    w_independent = truth.independent_weights(taus, q_pop, baseline_pop, family)
    w_independent = np.asarray(w_independent, dtype=np.float64)
    if w_independent.ndim != 2 or w_independent.shape[0] != len(taus):
        raise RuntimeError(f"independent_weights returned {w_independent.shape}")
    if not np.all(np.isfinite(w_independent)):
        raise RuntimeError("independent population weights are non-finite")
    weight_rank = int(np.linalg.matrix_rank(w_independent))
    weight_condition = float(np.linalg.cond(w_independent))
    if weight_rank < w_independent.shape[1] or not np.isfinite(weight_condition) or weight_condition > 1e12:
        raise RuntimeError(
            f"independent population W is rank deficient or ill-conditioned: "
            f"rank={weight_rank}/{w_independent.shape[1]}, cond={weight_condition:.6g}"
        )

    baseline_fit = truth.fit_baseline(y, family)
    if not isinstance(baseline_fit, Mapping):
        raise RuntimeError("truth.fit_baseline must return a mapping")
    _validate_baseline("same-sample fitted baseline", baseline_fit, family)
    rif, q_sample, fhat = core.compute_rif(
        y,
        taus,
        density_floor=float(config["density_floor"]),
        density_method=str(config["density_method"]),
    )
    X = np.ones((n, 1), dtype=np.float64)
    rif_resid, Q = core.residualize(rif, X)
    core.check_orthonormal_basis(Q, context="compact example")

    # This is deliberately a direct fixed-nuisance RIF--OLS calculation, not a
    # second package implementation.  It checks the wiring of the fast kernel.
    g_resid = g - Q @ (Q.T @ g)
    d_direct = float(np.sum(g_resid[:, 0] ** 2))
    beta_direct = np.sum(g_resid[:, [0]] * rif_resid, axis=0) / d_direct
    if not np.isfinite(d_direct) or d_direct <= 0:
        raise RuntimeError("independent RIF-OLS denominator is not positive")

    pooled_betas, pooled_d = core.compute_pooled_panel(
        g, rif_resid, Q,
        min_mac=float(config["min_mac"]),
        min_var=float(config["min_var"]),
        n_threads=int(config["threads"]),
        panel_size=int(config["panel_size"]),
    )
    pooled_betas = np.asarray(pooled_betas, dtype=np.float64)
    pooled_d = np.asarray(pooled_d, dtype=np.float64)
    if pooled_betas.shape != (1, len(taus)) or pooled_d.shape != (1,):
        raise RuntimeError(f"unexpected pooled panel output shapes: {pooled_betas.shape}, {pooled_d.shape}")
    pooled_beta = _finite_vector("pooled beta", pooled_betas[0], len(taus))
    d_resid = float(pooled_d[0])
    pooled_beta_diff = _assert_close("pooled fast panel beta", pooled_beta, beta_direct)
    d_diff = abs(d_resid - d_direct)
    if d_diff > ABS_REL_TOL + ABS_REL_TOL * max(1.0, abs(d_direct)):
        raise RuntimeError(f"pooled fast panel d_resid failed independent check: {d_diff}")
    sigma0 = np.asarray(compute_sigma0(rif_resid), dtype=np.float64)
    if sigma0.shape != (len(taus), len(taus)) or not np.all(np.isfinite(sigma0)):
        raise RuntimeError("pooled Sigma_0 is not finite with the expected shape")
    sigma_pooled = sigma0 / d_resid
    se_pooled = _finite_vector(
        "pooled tau SE", np.sqrt(np.maximum(np.diag(sigma_pooled), 0.0)), len(taus)
    )
    if np.any(se_pooled <= 0):
        raise RuntimeError("pooled tau covariance has a non-positive diagonal")

    # The paired jackknife is also the restructured/fast panel kernel.  It is
    # retained as a sensitivity source; no legacy convenience call is made.
    n_blocks = int(config["n_blocks"])
    if n_blocks <= len(taus):
        raise RuntimeError(
            f"jackknife blocks ({n_blocks}) must exceed the number of taus ({len(taus)})"
        )
    block_ids = balanced_block_ids(n, n_blocks, seed)
    n_blocks_nonempty = n_blocks
    plan = core.make_block_plan(
        rif_resid, Q, block_ids, n_blocks, n_threads=int(config["threads"])
    )
    jk_out = core.compute_block_finalize_panel_fast(
        g, rif_resid, Q, block_ids, n_blocks,
        min_mac=float(config["min_mac"]),
        min_var=float(config["min_var"]),
        n_threads=int(config["threads"]),
        panel_width=int(config["panel_size"]),
        plan=plan,
    )
    if not isinstance(jk_out, (tuple, list)) or len(jk_out) != 3:
        raise RuntimeError("fast jackknife panel kernel must return betas, SEs, packed covariance")
    jk_betas, jk_se, jk_upper = (np.asarray(value, dtype=np.float64) for value in jk_out)
    if jk_betas.shape != (1, len(taus)) or jk_se.shape != (1, len(taus)):
        raise RuntimeError(f"unexpected jackknife output shapes: {jk_betas.shape}, {jk_se.shape}")
    expected_upper = len(taus) * (len(taus) + 1) // 2
    if jk_upper.shape != (1, expected_upper):
        raise RuntimeError(f"unexpected packed jackknife covariance shape: {jk_upper.shape}")
    jk_beta = _finite_vector("jackknife beta", jk_betas[0], len(taus))
    jk_se_one = _finite_vector("jackknife tau SE", jk_se[0], len(taus))
    if np.any(jk_se_one <= 0):
        raise RuntimeError("jackknife tau covariance has a non-positive diagonal")
    jk_upper_one = _finite_vector("jackknife covariance", jk_upper[0], expected_upper)
    jk_beta_diff = _assert_close("jackknife fast panel beta", jk_beta, beta_direct)

    # Only small summaries and tau-level arrays are handed to R.  No g, y, Q or
    # RIF rows are serialised.
    truth_fields = {
        key: population[key]
        for key in required_truth
    }
    truth_fields["W_independent"] = w_independent
    return {
        "case": dict(case),
        "n": n,
        "seed": seed,
        "family": family,
        "taus": taus,
        "q_tau": q_sample,
        "fhat_tau": fhat,
        "baseline_fit": dict(baseline_fit),
        "population_truth": truth_fields,
        "beta_pooled": pooled_beta,
        "se_tau_pooled": se_pooled,
        "beta_jackknife": jk_beta,
        "se_tau_jackknife": jk_se_one,
        "sigma0": sigma0,
        "cov_upper_jackknife": jk_upper_one,
        "d_resid": d_resid,
        "n_blocks_requested": n_blocks,
        "n_blocks": n_blocks,
        "n_blocks_nonempty": n_blocks_nonempty,
        "density_method": str(config["density_method"]),
        "density_convention": (
            "Stage 1 pinned core.compute_rif: Silverman-bandwidth Gaussian KDE "
            "at empirical type-8 q_tau, with the configured positive density floor."
        ),
        "stage2_weight_convention": (
            "Stage 2 pinned make_weights_* uses analytic model density/CDF at q_tau; "
            "it does not substitute a KDE for the model density."
        ),
        "population_weight_diagnostics": {
            "rank": weight_rank,
            "n_columns": int(w_independent.shape[1]),
            "condition_number": weight_condition,
            "condition_number_limit": 1e12,
        },
        "independent_rif_ols": {
            "beta_max_abs_diff_pooled": pooled_beta_diff,
            "beta_max_abs_diff_jackknife": jk_beta_diff,
            "d_resid_abs_diff": d_diff,
            "absolute_tolerance": ABS_REL_TOL,
            "relative_tolerance": ABS_REL_TOL,
        },
        "package": None,
    }


def _run_stage2(rscript: str, stage2_script: Path, input_path: Path, output_path: Path,
                case_dir: Path, source: str, package_root: Path, package_commit: str,
                package_version: str, environment_manifest: Path, r_lib: str,
                timeout: float) -> dict:
    command = [
        rscript,
        "--vanilla",
        str(stage2_script),
        "--input", str(input_path),
        "--output", str(output_path),
        "--cov-source", source,
        "--package-root", str(package_root),
        "--package-commit", package_commit,
        "--package-version", package_version,
        "--environment-manifest", str(environment_manifest),
        "--r-lib", r_lib,
    ]
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
    })
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
    except subprocess.TimeoutExpired as exc:
        (case_dir / f"stage2-{source}.stdout").write_text(exc.stdout or "")
        (case_dir / f"stage2-{source}.stderr").write_text(exc.stderr or "")
        raise RuntimeError(f"Stage 2 {source} timed out after {timeout:g} seconds") from exc
    (case_dir / f"stage2-{source}.stdout").write_text(completed.stdout)
    (case_dir / f"stage2-{source}.stderr").write_text(completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Stage 2 {source} exited with status {completed.returncode}; see "
            f"stage2-{source}.stderr"
        )
    if not output_path.is_file():
        raise RuntimeError(f"Stage 2 {source} returned success without {output_path.name}")
    try:
        result = json.loads(output_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Stage 2 {source} output is not valid JSON") from exc
    if result.get("status") != "ok":
        raise RuntimeError(f"Stage 2 {source} did not report status=ok")
    return result


def run(config: Mapping[str, Any], repo_root: Path, output_root: Path,
        environment_manifest: Path, rscript: str | None = None,
        r_lib: str | None = None, timeout: float = 600.0,
        selected_cases: Sequence[str] | None = None) -> int:
    validate_config(config)
    environment = load_environment_manifest(environment_manifest, Path(config["package_root"]))
    manifest_rscript = Path(environment["rscript"])
    if rscript is not None and Path(rscript).resolve() != manifest_rscript:
        raise RuntimeError("--rscript must match the Rscript path in --environment-manifest")
    rscript = str(manifest_rscript)
    manifest_rlib = Path(environment["r_lib"])
    if r_lib is not None and Path(r_lib).resolve() != manifest_rlib:
        raise RuntimeError("--r-lib must match the private R library in --environment-manifest")
    r_lib = str(manifest_rlib)
    output_root = output_root.resolve()
    try:
        relative_output = output_root.relative_to(repo_root.resolve())
    except ValueError:
        relative_output = None
    if relative_output is not None and relative_output.parts and relative_output.parts[0] != ".local":
        raise RuntimeError(
            f"refusing output inside the Git worktree outside .local: {output_root}"
        )
    if output_root.exists():
        raise RuntimeError(
            f"refusing to reuse output directory {output_root}; every example run "
            "needs a new exclusive directory"
        )
    output_root.parent.mkdir(parents=True, exist_ok=True)
    output_root.mkdir()
    _write_json(output_root / "config.json", config)

    manifest: dict[str, Any] = {
        "status": "running",
        "run_id": config["run_id"],
        "package_commit": PINNED_PACKAGE_COMMIT,
        "package_version": PINNED_PACKAGE_VERSION,
        "stage1_distribution_version": PINNED_STAGE1_VERSION,
        "environment_manifest": environment["path"],
        "rscript": rscript,
        "r_library": r_lib,
        "output_root": str(output_root),
        "records": [],
    }
    try:
        core, package_provenance, compute_sigma0 = load_pinned_core(
            Path(config["package_root"]), environment
        )
        truth = _load_truth_module(repo_root)
        manifest["package"] = package_provenance
    except Exception as exc:
        manifest.update({"status": "failed_setup", "error": str(exc)})
        _write_json(output_root / "manifest.json", manifest)
        _write_json(output_root / "setup_failure.json", {
            "error": str(exc),
            "traceback": traceback.format_exc(),
        })
        return 1

    configured_ids = [case["id"] for case in config["cases"]]
    requested = set(configured_ids if selected_cases is None else selected_cases)
    unknown = requested.difference(configured_ids)
    if unknown:
        raise ValueError(f"unknown case selection: {sorted(unknown)}")
    cases = [case for case in config["cases"] if case["id"] in requested]
    taus = np.asarray(config["taus"], dtype=np.float64)
    stage2_script = Path(config["stage2_script"])
    if not stage2_script.is_absolute():
        stage2_script = repo_root / stage2_script
    if not stage2_script.is_file():
        raise RuntimeError(f"Stage 2 adapter is missing: {stage2_script}")
    try:
        population_cache = {
            case["id"]: truth.population_truth(dict(case), taus) for case in cases
        }
    except Exception as exc:
        manifest.update({"status": "failed_setup", "error": str(exc)})
        _write_json(output_root / "manifest.json", manifest)
        _write_json(output_root / "setup_failure.json", {
            "error": str(exc),
            "traceback": traceback.format_exc(),
        })
        return 1

    for case in cases:
        canonical_case_index = EXPECTED_CASE_IDS.index(case["id"])
        for replicate in range(1, int(config["replicates"]) + 1):
            seed = int(config["seed"]) + (canonical_case_index + 1) * 10000 + replicate
            case_dir = output_root / str(case["id"]) / f"rep-{replicate:03d}"
            case_dir.mkdir(parents=True)
            record: dict[str, Any] = {
                "case": case["id"],
                "family": case["family"],
                "replicate": replicate,
                "seed": seed,
                "path": str(case_dir),
                "status": "running",
            }
            try:
                payload = _stage1_example(
                    core, truth, case, int(config["n"]), seed, taus, config,
                    compute_sigma0, population=population_cache[case["id"]]
                )
                payload["package"] = package_provenance
                input_path = case_dir / "stage1_summary.json"
                _write_json(input_path, payload)
                record["stage1"] = payload["independent_rif_ols"]
                stage2_records = {}
                for source in ("pooled", "jackknife"):
                    output_path = case_dir / f"stage2-{source}.json"
                    stage2_records[source] = _run_stage2(
                        rscript, stage2_script, input_path, output_path, case_dir,
                        source, Path(config["package_root"]), PINNED_PACKAGE_COMMIT,
                        PINNED_PACKAGE_VERSION, Path(environment["path"]), r_lib, timeout,
                    )
                record["stage2"] = stage2_records
                record["status"] = "ok"
            except Exception as exc:
                record.update({
                    "status": "failed",
                    "error": str(exc),
                })
                _write_json(case_dir / "failure.json", {
                    "record": record,
                    "traceback": traceback.format_exc(),
                })
            manifest["records"].append(record)
            # Large bank: retain each record immediately, but avoid quadratic
            # rewriting of the full manifest. Numerical processing is unchanged.
            if config.get("profile") == "normal_replication_v1":
                _write_json(case_dir / "record.json", record)
                _write_json(output_root / "progress.json", {
                    "status": "running", "n_attempted": len(manifest["records"]),
                    "last_case": case["id"], "last_replicate": replicate,
                    "last_status": record["status"],
                })
                if len(manifest["records"]) % 25 == 0 or record["status"] != "ok":
                    _write_json(output_root / "manifest.json", manifest)
            else:
                _write_json(output_root / "manifest.json", manifest)

    failures = sum(record["status"] != "ok" for record in manifest["records"])
    manifest["status"] = "ok" if failures == 0 else "failed_examples"
    manifest["n_attempted"] = len(manifest["records"])
    manifest["n_failed"] = failures
    _write_json(output_root / "manifest.json", manifest)
    if config.get("profile") == "normal_replication_v1":
        _write_json(output_root / "progress.json", {
            "status": manifest["status"], "n_attempted": len(manifest["records"]),
            "n_failed": failures,
        })
    return 0 if failures == 0 else 1


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    default_config = Path(__file__).resolve().parent / "config" / "examples.json"
    parser.add_argument("--config", type=Path, default=default_config)
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="new exclusive scratch directory; defaults to config output_root")
    parser.add_argument("--environment-manifest", type=Path, required=True,
                        help="bootstrap.sh environment.json; its Python, Rscript and private Rlib are enforced")
    parser.add_argument("--rscript", default=None,
                        help="explicit override, which must match the manifest Rscript")
    parser.add_argument("--r-lib", default=None,
                        help="explicit override, which must match the manifest private R library")
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--case", dest="cases", action="append", default=None,
                        help="optional case restriction for debugging; default runs all eight")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    config = json.loads(args.config.read_text())
    repo_root = Path(__file__).resolve().parents[1]
    output_root = args.out_dir or Path(config["output_root"])
    return run(
        config,
        repo_root=repo_root,
        output_root=output_root,
        environment_manifest=args.environment_manifest,
        rscript=args.rscript,
        r_lib=args.r_lib,
        timeout=args.timeout,
        selected_cases=args.cases,
    )


if __name__ == "__main__":
    raise SystemExit(main())
