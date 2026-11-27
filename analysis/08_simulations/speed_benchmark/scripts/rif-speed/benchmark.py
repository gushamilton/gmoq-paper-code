"""Timing harness for the Stage 1 comparison.

What: runs every workload x method x repeat as a separate process wrapped in GNU time
  (wall time and peak memory), pinned to one core with all numerical-library threads set
  to 1, in a seeded interleaved order, and records each run's real exit status. It contains
  no estimator.
Paper: Supplementary Fig. S2; Supplementary Table S2 (study 2)
Input: a command manifest built by figures/run_stage1.py
Output: one directory per run with timing.json, GNU time output and logs
Runs: on a cluster compute node (refuses to run unless VALIDATION_COMPUTE_HOST=1)
"""

from __future__ import annotations

import json
import math
import os
import random
import re
import shutil
import signal
import string
import subprocess
import time
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

SCHEMA_VERSION = "validation-benchmark-manifest-1"
RESULT_SCHEMA_VERSION = "validation-benchmark-result-1"
REQUIRED_REPEATS = 3
THREAD_ENV_NAMES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RAYON_NUM_THREADS",
)
_PACKAGE_PIN = re.compile(r"^[0-9a-fA-F]{40}$")
_SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


class BenchmarkError(ValueError):
    """Raised when a benchmark manifest is unsafe or incomplete."""


def _require_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise BenchmarkError(f"{name} must be an object")
    return value


def _required_text(mapping: Mapping[str, Any], name: str) -> str:
    value = mapping.get(name)
    if value is None or not str(value):
        raise BenchmarkError(f"{name} is required and cannot be empty")
    return str(value)


def _safe_id(value: Any, field: str) -> str:
    result = str(value)
    if not result or not _SAFE_ID.fullmatch(result):
        raise BenchmarkError(
            f"{field} must contain only letters, numbers, '.', '_' or '-': {result!r}"
        )
    return result


def _positive_int(value: Any, field: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise BenchmarkError(f"{field} must be a positive integer") from exc
    if result <= 0:
        raise BenchmarkError(f"{field} must be a positive integer")
    return result


def _nonnegative_int(value: Any, field: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise BenchmarkError(f"{field} must be a non-negative integer") from exc
    if result < 0:
        raise BenchmarkError(f"{field} must be a non-negative integer")
    return result


def _finite_positive(value: Any, field: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise BenchmarkError(f"{field} must be a positive number") from exc
    if not math.isfinite(result) or result <= 0:
        raise BenchmarkError(f"{field} must be a positive number")
    return result


def _contains_placeholder(command: Sequence[str], name: str) -> bool:
    needle = "{" + name + "}"
    return any(needle in token for token in command)


def _validate_command(command: Any, method_name: str) -> list[str]:
    if not isinstance(command, list) or not command or not all(
        isinstance(token, str) and token for token in command
    ):
        raise BenchmarkError(
            f"{method_name}.command must be a non-empty list of non-empty strings; "
            "shell command strings are not accepted"
        )
    # Do not silently turn a manifest into a shell script.  An explicit Python,
    # R, or compiled executable may still be the first token, but shell -c is
    # not an acceptable generic benchmark command.
    if any(token in {"sh", "bash", "zsh", "fish", "-c", "--command"} for token in command):
        raise BenchmarkError(f"{method_name}.command must not invoke a shell")
    formatter = string.Formatter()
    for token in command:
        try:
            fields = [field for _, field, _, _ in formatter.parse(token) if field]
        except ValueError as exc:
            raise BenchmarkError(f"invalid command placeholder in {method_name}") from exc
        for field in fields:
            if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", field):
                raise BenchmarkError(f"unsupported command placeholder {{{field}}}")
    return list(command)


def _validate_env(env: Any, owner: str) -> dict[str, str]:
    if env is None:
        return {}
    if not isinstance(env, Mapping):
        raise BenchmarkError(f"{owner}.env must be an object")
    result: dict[str, str] = {}
    for key, value in env.items():
        key = str(key)
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
            raise BenchmarkError(f"invalid environment variable name {key!r}")
        result[key] = str(value)
    return result


def validate_manifest(manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Validate and copy an explicit command manifest.

    Required command placeholders are ``{output_dir}``, ``{threads}``, and
    ``{package_pin}``.  Requiring them prevents an apparently matched benchmark
    from quietly writing elsewhere, using uncontrolled threads, or running an
    unpinned package.  Method-specific configurations can add placeholders
    listed by a workload's ``placeholders`` object.
    """

    manifest = _require_mapping(manifest, "manifest")
    schema_version = str(manifest.get("schema_version", SCHEMA_VERSION))
    if schema_version != SCHEMA_VERSION:
        raise BenchmarkError(f"unsupported benchmark schema_version {schema_version!r}")
    experiment_id = _safe_id(_required_text(manifest, "experiment_id"), "experiment_id")
    package_pin = _required_text(manifest, "package_pin")
    if not _PACKAGE_PIN.fullmatch(package_pin):
        raise BenchmarkError("package_pin must be a 40-character source commit SHA")
    output_root = Path(_required_text(manifest, "output_root"))
    if not output_root.is_absolute():
        raise BenchmarkError("output_root must be an absolute exclusive output path")
    if str(output_root) in {"/", "/tmp", "/var", "/home", "/user", "."}:
        raise BenchmarkError(f"refusing unsafe output_root {output_root}")

    repeats = _positive_int(manifest.get("repeats", REQUIRED_REPEATS), "repeats")
    if repeats != REQUIRED_REPEATS:
        raise BenchmarkError("the validation benchmark contract requires exactly three repeats")
    physical_cores = _positive_int(manifest.get("physical_cores"), "physical_cores")
    threads = _positive_int(manifest.get("threads", physical_cores), "threads")
    if threads > physical_cores:
        raise BenchmarkError("threads cannot exceed declared physical_cores")
    timeout_seconds = _finite_positive(
        manifest.get("timeout_seconds"), "timeout_seconds"
    )
    try:
        interleave_seed = int(manifest.get("interleave_seed"))
    except (TypeError, ValueError) as exc:
        raise BenchmarkError("interleave_seed is required and must be an integer") from exc
    warmup_count = _nonnegative_int(manifest.get("warmup_count", 0), "warmup_count")
    cpu_ids_raw = manifest.get("physical_cpu_ids")
    cpu_ids: list[int] | None
    if cpu_ids_raw is None:
        raise BenchmarkError(
            "physical_cpu_ids is required so the benchmark cannot silently "
            "oversubscribe or run on an unmatched CPU set"
        )
    else:
        if not isinstance(cpu_ids_raw, list) or not cpu_ids_raw:
            raise BenchmarkError("physical_cpu_ids must be a non-empty list when supplied")
        try:
            cpu_ids = [int(item) for item in cpu_ids_raw]
        except (TypeError, ValueError) as exc:
            raise BenchmarkError("physical_cpu_ids must contain integers") from exc
        if len(cpu_ids) != physical_cores or len(set(cpu_ids)) != len(cpu_ids):
            raise BenchmarkError(
                "physical_cpu_ids must contain one distinct ID per physical core"
            )
        if min(cpu_ids) < 0:
            raise BenchmarkError("physical_cpu_ids cannot contain negative IDs")

    workloads_raw = manifest.get("workloads")
    methods_raw = manifest.get("methods")
    if not isinstance(workloads_raw, list) or not workloads_raw:
        raise BenchmarkError("workloads must be a non-empty list")
    if not isinstance(methods_raw, list) or not methods_raw:
        raise BenchmarkError("methods must be a non-empty list")

    workloads: list[dict[str, Any]] = []
    workload_ids: set[str] = set()
    allowed_modes = {"isolated_single_snp", "amortized_scan"}
    for index, raw in enumerate(workloads_raw):
        workload = _require_mapping(raw, f"workloads[{index}]")
        workload_id = _safe_id(
            _required_text(workload, "workload_id"), f"workloads[{index}].workload_id"
        )
        if workload_id in workload_ids:
            raise BenchmarkError(f"duplicate workload_id {workload_id!r}")
        workload_ids.add(workload_id)
        mode = str(workload.get("execution_mode", ""))
        if mode not in allowed_modes:
            raise BenchmarkError(
                f"workloads[{index}].execution_mode must be one of {sorted(allowed_modes)}"
            )
        placeholders = workload.get("placeholders", {})
        if not isinstance(placeholders, Mapping):
            raise BenchmarkError(f"workloads[{index}].placeholders must be an object")
        workload_copy = dict(workload)
        workload_copy["workload_id"] = workload_id
        workload_copy["label"] = str(workload.get("label", workload_id))
        workload_copy["execution_mode"] = mode
        workload_copy["placeholders"] = {str(k): str(v) for k, v in placeholders.items()}
        covariance_glob = workload.get("covariance_glob")
        if covariance_glob is not None:
            covariance_glob = str(covariance_glob)
            if Path(covariance_glob).is_absolute() or ".." in Path(covariance_glob).parts:
                raise BenchmarkError("covariance_glob must be relative to output_dir")
            workload_copy["covariance_glob"] = covariance_glob
        timing_file = workload.get("timing_file")
        if timing_file is not None:
            timing_file = str(timing_file)
            if Path(timing_file).is_absolute() or ".." in Path(timing_file).parts:
                raise BenchmarkError("timing_file must be relative to output_dir")
            workload_copy["timing_file"] = timing_file
        workloads.append(workload_copy)

    methods: list[dict[str, Any]] = []
    method_ids: set[str] = set()
    base_allowed = {
        "experiment_id",
        "package_pin",
        "output_root",
        "output_dir",
        "workload_id",
        "workload_label",
        "execution_mode",
        "method_id",
        "method_label",
        "repeat",
        "order",
        "threads",
        "physical_cores",
        "input_path",
    }
    for index, raw in enumerate(methods_raw):
        method = _require_mapping(raw, f"methods[{index}]")
        method_id = _safe_id(
            _required_text(method, "method_id"), f"methods[{index}].method_id"
        )
        if method_id in method_ids:
            raise BenchmarkError(f"duplicate method_id {method_id!r}")
        method_ids.add(method_id)
        command = _validate_command(method.get("command"), f"methods[{index}]")
        for placeholder in ("output_dir", "threads", "package_pin"):
            if not _contains_placeholder(command, placeholder):
                raise BenchmarkError(
                    f"methods[{index}].command must contain {{{placeholder}}}"
                )
        method_copy = dict(method)
        method_copy["method_id"] = method_id
        method_copy["label"] = str(method.get("label", method_id))
        method_copy["command"] = command
        method_copy["env"] = _validate_env(method.get("env"), f"methods[{index}]")
        method_placeholders = method.get("placeholders", {})
        if not isinstance(method_placeholders, Mapping):
            raise BenchmarkError(f"methods[{index}].placeholders must be an object")
        method_copy["placeholders"] = {
            str(key): str(value) for key, value in method_placeholders.items()
        }
        cwd = method.get("cwd")
        if cwd is not None:
            cwd_path = Path(str(cwd))
            if not cwd_path.is_absolute():
                raise BenchmarkError("method cwd must be absolute")
            method_copy["cwd"] = str(cwd_path)
        # Reject fields that look like supplied performance estimates.  The
        # manifest is a command specification, never a speed table.
        forbidden = {"expected_seconds", "runtime_seconds", "speed", "peak_rss_kib"}
        present = forbidden.intersection(method)
        if present:
            raise BenchmarkError(
                "method manifest must not contain invented performance values: "
                + ", ".join(sorted(present))
            )
        methods.append(method_copy)

    global_env = _validate_env(manifest.get("env"), "manifest")
    global_placeholders_raw = manifest.get("placeholders", {})
    if not isinstance(global_placeholders_raw, Mapping):
        raise BenchmarkError("manifest.placeholders must be an object")
    global_placeholders = {
        str(key): str(value) for key, value in global_placeholders_raw.items()
    }
    compute_host_only = manifest.get("compute_host_only", True)
    if compute_host_only is not True:
        raise BenchmarkError("compute_host_only must remain true")
    time_path = manifest.get("gnu_time", "/usr/bin/time")
    if not str(time_path):
        raise BenchmarkError("gnu_time cannot be empty")

    return {
        "schema_version": SCHEMA_VERSION,
        "experiment_id": experiment_id,
        "package_pin": package_pin.lower(),
        "output_root": str(output_root),
        "repeats": repeats,
        "physical_cores": physical_cores,
        "threads": threads,
        "timeout_seconds": timeout_seconds,
        "interleave_seed": interleave_seed,
        "warmup_count": warmup_count,
        "physical_cpu_ids": cpu_ids,
        "compute_host_only": True,
        "gnu_time": str(time_path),
        "env": global_env,
        "placeholders": global_placeholders,
        "workloads": workloads,
        "methods": methods,
    }


def load_manifest(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.is_file():
        raise BenchmarkError(f"manifest does not exist: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise BenchmarkError(f"invalid JSON manifest: {path}") from exc
    return validate_manifest(_require_mapping(value, "manifest"))


def _schedule_items(manifest: Mapping[str, Any]) -> list[dict[str, Any]]:
    arms = [
        (workload, method)
        for workload in manifest["workloads"]
        for method in manifest["methods"]
    ]
    rng = random.Random(int(manifest["interleave_seed"]))
    schedule: list[dict[str, Any]] = []
    order = 0
    for warmup_index in range(int(manifest.get("warmup_count", 0))):
        indices = list(range(len(arms)))
        rng.shuffle(indices)
        for index in indices:
            workload, method = arms[index]
            schedule.append(
                {
                    "order": order,
                    "repeat": 0,
                    "warmup_index": warmup_index + 1,
                    "is_warmup": True,
                    "workload_id": workload["workload_id"],
                    "method_id": method["method_id"],
                }
            )
            order += 1
    for repeat in range(1, int(manifest["repeats"]) + 1):
        indices = list(range(len(arms)))
        rng.shuffle(indices)
        for index in indices:
            workload, method = arms[index]
            schedule.append(
                {
                    "order": order,
                    "repeat": repeat,
                    "warmup_index": 0,
                    "is_warmup": False,
                    "workload_id": workload["workload_id"],
                    "method_id": method["method_id"],
                }
            )
            order += 1
    return schedule


def build_schedule(manifest: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return the deterministic interleaved warm-up and three-repeat schedule."""

    checked = validate_manifest(manifest)
    return _schedule_items(checked)


class _PlaceholderMap(dict[str, str]):
    def __missing__(self, key: str) -> str:
        raise BenchmarkError(f"command uses unknown placeholder {{{key}}}")


def _expand_command(
    command: Sequence[str], values: Mapping[str, Any]
) -> list[str]:
    mapping = _PlaceholderMap({key: str(value) for key, value in values.items()})
    expanded: list[str] = []
    for token in command:
        try:
            expanded.append(token.format_map(mapping))
        except (KeyError, ValueError) as exc:
            raise BenchmarkError(f"cannot expand command token {token!r}") from exc
    if any(not token for token in expanded):
        raise BenchmarkError("command expansion produced an empty argument")
    return expanded


def _compute_host_allowed() -> bool:
    # The marker is set by the cluster job wrapper.  Requiring it is preferable
    # to guessing from a hostname, which can silently run a costly benchmark on
    # a login node or a laptop.
    return os.environ.get("VALIDATION_COMPUTE_HOST") == "1"


def _assert_compute_host() -> None:
    if not _compute_host_allowed():
        raise BenchmarkError(
            "benchmark execution is compute-host only; set "
            "VALIDATION_COMPUTE_HOST=1 from the cluster job wrapper"
        )


def _find_gnu_time(path_value: str) -> str:
    candidates = [path_value]
    if path_value == "/usr/bin/time":
        candidates.extend(["gtime", "time"])
    for candidate in candidates:
        resolved = candidate if os.path.isabs(candidate) else shutil.which(candidate)
        if not resolved or not os.access(resolved, os.X_OK):
            continue
        try:
            probe = subprocess.run(
                [resolved, "--version"],
                check=False,
                capture_output=True,
                text=True,
                timeout=5,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        version = (probe.stdout + probe.stderr).lower()
        if "gnu time" in version:
            return resolved
    raise BenchmarkError("GNU time is required for peak RSS measurement")


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _tree_bytes(path: Path) -> int | None:
    if not path.exists():
        return None
    total = 0
    for child in path.rglob("*"):
        if child.is_file():
            try:
                total += child.stat().st_size
            except OSError:
                pass
    return total


def _parse_time_file(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line:
            continue
        key, raw = line.split("=", 1)
        raw = raw.strip()
        try:
            if key in {"max_rss_kib", "exit_status"}:
                values[key] = int(float(raw))
            else:
                values[key] = float(raw)
        except ValueError:
            values[key] = raw
    return values


def _read_breakdown(output_dir: Path, workload: Mapping[str, Any]) -> dict[str, Any]:
    timing_file = workload.get("timing_file")
    result: dict[str, Any] = {
        "breakdown_source": None,
        "input_seconds": None,
        "shared_setup_seconds": None,
        "nuisance_seconds": None,
        "fit_seconds": None,
        "write_seconds": None,
    }
    if timing_file is None:
        return result
    path = output_dir / str(timing_file)
    if not path.is_file():
        result["breakdown_source"] = "declared_but_missing"
        return result
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        result["breakdown_source"] = "declared_but_invalid"
        return result
    if not isinstance(value, Mapping):
        result["breakdown_source"] = "declared_but_invalid"
        return result
    result["breakdown_source"] = str(path.name)
    for key in (
        "input_seconds",
        "shared_setup_seconds",
        "nuisance_seconds",
        "fit_seconds",
        "write_seconds",
    ):
        if key in value:
            try:
                numeric = float(value[key])
            except (TypeError, ValueError):
                continue
            if math.isfinite(numeric) and numeric >= 0:
                result[key] = numeric
    return result


def _covariance_bytes(output_dir: Path, workload: Mapping[str, Any]) -> int | None:
    pattern = workload.get("covariance_glob")
    if pattern is None:
        return None
    total = 0
    found = False
    for path in output_dir.glob(str(pattern)):
        if path.is_file():
            found = True
            total += path.stat().st_size
    return total if found else 0


def _kill_process_group(process: subprocess.Popen[Any]) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        try:
            process.terminate()
        except ProcessLookupError:
            return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            try:
                process.kill()
            except ProcessLookupError:
                pass
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            pass


def _run_one(
    manifest: Mapping[str, Any],
    item: Mapping[str, Any],
    workload: Mapping[str, Any],
    method: Mapping[str, Any],
    *,
    gnu_time: str,
    run_dir: Path,
) -> dict[str, Any]:
    run_dir.mkdir(parents=True, exist_ok=False)
    output_dir = run_dir / "output"
    output_dir.mkdir()
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    time_path = run_dir / "gnu_time.txt"
    values: dict[str, Any] = {
        "experiment_id": manifest["experiment_id"],
        "package_pin": manifest["package_pin"],
        "output_root": manifest["output_root"],
        "output_dir": str(output_dir),
        "workload_id": workload["workload_id"],
        "workload_label": workload.get("label", workload["workload_id"]),
        "execution_mode": workload["execution_mode"],
        "method_id": method["method_id"],
        "method_label": method.get("label", method["method_id"]),
        "repeat": item["repeat"],
        "order": item["order"],
        "threads": manifest["threads"],
        "physical_cores": manifest["physical_cores"],
        "physical_cpu_ids": manifest["physical_cpu_ids"],
        "input_path": workload.get("input_path", ""),
    }
    values.update(manifest.get("placeholders", {}))
    values.update(workload.get("placeholders", {}))
    values.update(method.get("placeholders", {}))
    command = _expand_command(method["command"], values)
    time_format = (
        "elapsed_seconds=%e\\nuser_seconds=%U\\nsys_seconds=%S\\n"
        "max_rss_kib=%M\\nexit_status=%x\\n"
    )
    wrapped = [gnu_time, "-f", time_format, "-o", str(time_path), *command]
    env = os.environ.copy()
    env.update({key: str(value) for key, value in manifest.get("env", {}).items()})
    env.update({key: str(value) for key, value in method.get("env", {}).items()})
    env["VALIDATION_PACKAGE_PIN"] = str(manifest["package_pin"])
    for name in THREAD_ENV_NAMES:
        env[name] = str(manifest["threads"])
    env["VALIDATION_PHYSICAL_CORES"] = str(manifest["physical_cores"])
    # Keep the child command free of a shell and start a process group so a
    # timeout covers descendants as well as the direct executable.
    preexec = None
    cpu_ids = manifest.get("physical_cpu_ids")
    if cpu_ids is not None:
        def set_affinity() -> None:  # pragma: no cover - executed only on the cluster.
            os.sched_setaffinity(0, set(int(value) for value in cpu_ids))
        preexec = set_affinity
    cwd = method.get("cwd")
    if cwd is None:
        cwd = str(output_dir)
    start = time.monotonic()
    timed_out = False
    launch_error: str | None = None
    return_code: int | None = None
    try:
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
            "w", encoding="utf-8"
        ) as stderr:
            process = subprocess.Popen(
                wrapped,
                cwd=cwd,
                env=env,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
                preexec_fn=preexec,
                text=False,
            )
            try:
                process.communicate(timeout=float(manifest["timeout_seconds"]))
            except subprocess.TimeoutExpired:
                timed_out = True
                _kill_process_group(process)
            return_code = process.returncode
    except (OSError, ValueError) as exc:
        launch_error = f"{type(exc).__name__}: {exc}"
    wall_measured = time.monotonic() - start
    timed = _parse_time_file(time_path)
    if timed_out:
        status = "timeout"
    elif launch_error is not None:
        status = "failed"
    elif return_code == 0:
        status = "success"
    else:
        status = "failed"
    result: dict[str, Any] = {
        "schema_version": RESULT_SCHEMA_VERSION,
        **values,
        "status": status,
        "timed_out": timed_out,
        "reason": (
            "whole_process_timeout"
            if timed_out
            else ("launch_error" if launch_error is not None else None)
        ),
        "launch_error": launch_error,
        # return_code is the GNU-time wrapper status; child_return_code is
        # GNU time's %x status for the actual benchmark command.
        "return_code": return_code,
        "child_return_code": timed.get("exit_status"),
        "return_status": (
            "timeout"
            if timed_out
            else ("launch_error" if launch_error is not None else str(return_code))
        ),
        "child_return_status": (
            "timeout"
            if timed_out
            else (
                "launch_error"
                if launch_error is not None
                else str(timed.get("exit_status"))
            )
        ),
        "wall_seconds": timed.get("elapsed_seconds", wall_measured),
        "wall_seconds_measured_by_wrapper": timed.get("elapsed_seconds"),
        "wall_seconds_measured_by_orchestrator": wall_measured,
        "process_user_seconds": timed.get("user_seconds"),
        "process_system_seconds": timed.get("sys_seconds"),
        "process_cpu_seconds": (
            timed.get("user_seconds", 0.0) + timed.get("sys_seconds", 0.0)
            if isinstance(timed.get("user_seconds"), (int, float))
            and isinstance(timed.get("sys_seconds"), (int, float))
            else None
        ),
        "peak_rss_kib": timed.get("max_rss_kib"),
        "gnu_time_exit_status": timed.get("exit_status"),
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
        "gnu_time_path": gnu_time,
        # These are measured after completion, never estimated from a method
        # label.  A missing covariance pattern stays null, not zero.
        "output_bytes": _tree_bytes(output_dir),
        "covariance_object_bytes": _covariance_bytes(output_dir, workload),
        "total_includes_input": True,
        "total_includes_shared_setup": True,
        "total_includes_nuisance_fit": True,
        "is_warmup": item["is_warmup"],
        "warmup_index": item["warmup_index"],
    }
    result.update(_read_breakdown(output_dir, workload))
    return result


def run_benchmark(
    manifest: Mapping[str, Any] | str | Path,
    *,
    dry_run: bool = False,
) -> list[dict[str, Any]]:
    """Execute a validated manifest, or return its schedule in dry-run mode.

    The output root must not already exist.  This is an exclusive-output guard:
    a rerun gets a new experiment ID/path instead of overwriting a failure.
    ``dry_run`` performs validation and schedule construction only and never
    starts a subprocess, GNU time, allocation, or shell.
    """

    if isinstance(manifest, (str, Path)):
        checked = load_manifest(manifest)
    else:
        checked = validate_manifest(manifest)
    schedule = _schedule_items(checked)
    if dry_run:
        return schedule
    _assert_compute_host()
    output_root = Path(checked["output_root"])
    if output_root.exists():
        raise BenchmarkError(
            f"exclusive benchmark output already exists; refusing overwrite: {output_root}"
        )
    output_root.mkdir(parents=True, exist_ok=False)
    _write_json(output_root / "manifest.json", checked)
    gnu_time = _find_gnu_time(checked["gnu_time"])
    workload_by_id = {item["workload_id"]: item for item in checked["workloads"]}
    method_by_id = {item["method_id"]: item for item in checked["methods"]}
    results: list[dict[str, Any]] = []
    for item in schedule:
        workload = workload_by_id[item["workload_id"]]
        method = method_by_id[item["method_id"]]
        run_name = f"{int(item['order']):04d}_{workload['workload_id']}_{method['method_id']}"
        run_dir = output_root / "runs" / run_name
        try:
            result = _run_one(
                checked,
                item,
                workload,
                method,
                gnu_time=gnu_time,
                run_dir=run_dir,
            )
        except Exception as exc:  # preserve an orchestrator-side failure.
            result = {
                "schema_version": RESULT_SCHEMA_VERSION,
                "experiment_id": checked["experiment_id"],
                "package_pin": checked["package_pin"],
                "workload_id": workload["workload_id"],
                "method_id": method["method_id"],
                "repeat": item["repeat"],
                "order": item["order"],
                "is_warmup": item["is_warmup"],
                "status": "failed",
                "timed_out": False,
                "reason": "orchestrator_error",
                "launch_error": f"{type(exc).__name__}: {exc}",
                "return_code": None,
                "return_status": "orchestrator_error",
                "wall_seconds": None,
                "peak_rss_kib": None,
                "output_bytes": None,
                "covariance_object_bytes": None,
            }
        results.append(result)
        _write_json(output_root / "results.json", results)
    _write_json(output_root / "results.json", results)
    return results


# Explicit aliases for integration wrappers.
run_manifest = run_benchmark


__all__ = [
    "BenchmarkError",
    "SCHEMA_VERSION",
    "RESULT_SCHEMA_VERSION",
    "validate_manifest",
    "load_manifest",
    "build_schedule",
    "run_benchmark",
    "run_manifest",
]


if __name__ == "__main__":  # pragma: no cover - integration adapter only.
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    result = run_benchmark(args.manifest, dry_run=args.dry_run)
    print(json.dumps(result, indent=2, sort_keys=True))
