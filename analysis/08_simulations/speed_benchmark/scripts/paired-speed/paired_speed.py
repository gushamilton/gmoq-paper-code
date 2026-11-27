#!/usr/bin/env python3
"""Paired single-thread timing of pooled G-MoQ and conditional mixture ML.

The null generator, deterministic six-start ML comparator and phenotype-only
warm-start adapter are imported unchanged. Each method runs in a fresh process
on the same frozen 100-variant panel; a smaller pilot uses its first columns.
"""
from __future__ import annotations

import time
ENTRY_TIME = time.perf_counter()

import argparse
import csv
import hashlib
import importlib
import json
import math
import os
import platform
import signal
import subprocess
import sys
import traceback
from contextlib import nullcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
METHODS = ("gmoq_pooled", "gmoq_jackknife", "mixture_ml", "mixture_ml_warm")
NS = (10000, 50000, 400000)
PARAMS = ("gamma", "delta1", "delta2")
THREAD_ENV = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
NULL_CASE = dict(family="mixture", maf=.3, pi0=.5, mu1=-2., mu2=2.,
                 sd1=1., sd2=1., gamma=0., delta1=0., delta2=0.)
FIT_CASE = dict(family="mixture", maf=.3)


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if hasattr(value, "tolist"):
        return clean(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")
    os.replace(tmp, path)


def read_json(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def csv_write(path, rows, fields=None):
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else ["status"])
    tmp = Path(str(path) + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(clean(rows))
    os.replace(tmp, path)


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def source_identity(a):
    paper = Path(a.paper_root).resolve()
    files = {
        "runner": HERE / "paired_speed.py", "stage2": HERE / "paired_stage2.R",
        "warm_adapter": HERE / "warm_ml.py", "truth": paper / "validation/truth.py",
        "comparator": paper / "validation/comparators.py",
        "block_generator": paper / "validation/run_examples.py",
    }
    return {
        "paper_root": str(paper), "paper_commit": git(paper, "rev-parse", "HEAD"),
        "files": {key: {"path": str(path), "sha256": sha(path)} for key, path in files.items()},
        "pin_json": str(Path(a.pin_json).resolve()), "pin_sha256": sha(a.pin_json),
        "pin": read_json(a.pin_json),
    }


def validation_imports(paper_root):
    sys.path.insert(0, str(Path(paper_root).resolve()))
    from validation import comparators, truth
    return comparators, truth


def load_core(pin):
    root = Path(pin["stage1_root"]).resolve()
    if git(root, "rev-parse", "HEAD") != pin["commit"]:
        raise RuntimeError("Stage 1 source commit differs from PIN.json")
    if git(root, "status", "--porcelain", "--untracked-files=no"):
        raise RuntimeError("Stage 1 has tracked modifications")
    sys.path.insert(0, str(root / "stage1"))
    core = importlib.import_module("gmoq.core")
    extension = Path(core.cpp_ext.__file__).resolve()
    if sha(extension) != pin["stage1_extension_sha256"]:
        raise RuntimeError("Stage 1 extension differs from PIN.json")
    if int(core.cpp_ext.openmp_info()["max_threads"]) != 1:
        raise RuntimeError("Stage 1 did not initialise with one OpenMP thread")
    for name in ("compute_rif_nuisance", "make_pooled_bin_plan", "compute_pooled_binned"):
        if not hasattr(core, name):
            raise RuntimeError(f"pinned production core lacks {name}")
    return core, {"root": str(root), "extension": str(extension), "extension_sha256": sha(extension)}


def input_panel(a, folder, identity):
    import numpy as np
    _, truth = validation_imports(a.paper_root)
    metadata = dict(schema_version=1, n=a.n, variants=a.variants, generator_variants=100,
                    rep=a.rep, phenotype_seed=97000000 + a.rep - 1,
                    genotype_seed=97000000 + a.rep - 1 + 7919,
                    generating_case=NULL_CASE, source=identity)
    folder.mkdir(parents=True, exist_ok=True)
    manifest = folder / "manifest.json"
    if manifest.exists():
        old = read_json(manifest)
        if {k: v for k, v in old.items() if k != "arrays"} != metadata:
            raise RuntimeError("existing cell has different inputs, source hashes or pin")
        for name, info in old["arrays"].items():
            if sha(folder / name) != info["sha256"]:
                raise RuntimeError(f"existing input hash mismatch: {name}")
        return old
    if any(folder.iterdir()):
        raise RuntimeError("input directory is incomplete; inspect it before removing/retrying")
    seed = metadata["phenotype_seed"]
    y = truth.draw(NULL_CASE, a.n, seed)["y"]
    # Draw all 100 columns exactly as the frozen benchmark does, even for pilots.
    # Drawing (N, pilot_M) would change the variants and would not be a matched pilot.
    G = np.random.default_rng(seed + 7919).binomial(2, .3, size=(a.n, 100)).astype(np.int8)
    np.save(folder / "y.npy", np.asarray(y, dtype=float), allow_pickle=False)
    np.save(folder / "G.npy", G[:, :a.variants], allow_pickle=False)
    del G, y
    metadata["arrays"] = {name: {"sha256": sha(folder / name)} for name in ("y.npy", "G.npy")}
    write_json(manifest, metadata)
    return metadata


def normalise_parameter(name):
    return {"beta_1": "delta1", "beta_2": "delta2", "delta_1": "delta1", "delta_2": "delta2"}.get(name, name)


def fit_valid(fit):
    return fit.get("status") == "ok" and all(
        math.isfinite(float(fit.get("estimates", {}).get(p, float("nan")))) and
        math.isfinite(float(fit.get("standard_errors", {}).get(p, float("nan")))) and
        float(fit["standard_errors"][p]) > 0 for p in PARAMS)


def child(a):
    """Exactly one method, one panel, one fresh single-thread process."""
    for name in THREAD_ENV:
        if os.environ.get(name) != "1":
            raise RuntimeError(f"{name} must be set to 1 before Python starts")
    import numpy as np
    comparators, truth = validation_imports(a.paper_root)
    pin = read_json(a.pin_json)
    work = Path(a.workdir)
    manifest = read_json(Path(a.input_dir) / "manifest.json")
    core = None
    package = None
    if a.method.startswith("gmoq_"):
        core, package = load_core(pin)
    elif a.method == "mixture_ml_warm":
        from warm_ml import warm_start
    import_ready = time.perf_counter()
    t = time.perf_counter()
    G = np.load(Path(a.input_dir) / "G.npy", allow_pickle=False)
    y = np.load(Path(a.input_dir) / "y.npy", allow_pickle=False)
    n, M = G.shape
    data_load = time.perf_counter() - t
    if y.shape != (n,) or n != manifest["n"] or M != manifest["variants"]:
        raise RuntimeError("input dimensions differ from the frozen manifest")
    started = time.perf_counter()
    components = {}
    fits = []
    if a.method in ("mixture_ml", "mixture_ml_warm"):
        if a.method.endswith("_warm"):
            context = warm_start(y, FIT_CASE)
        else:
            context = nullcontext(None)
        with (work / "variant_fits.jsonl").open("w") as log, context as warm:
            if warm is not None:
                components["baseline_fit_seconds"] = warm.baseline_wall_seconds
            for j in range(M):
                t = time.perf_counter()
                fit = warm(G[:, j]) if warm is not None else comparators.fit_mixture_ml(G[:, j], y, FIT_CASE)
                fit["snp_index"] = j
                fit["measured_fit_seconds"] = time.perf_counter() - t
                fit["valid_estimate_se"] = fit_valid(fit)
                # Preserve the comparator's complete Hessian, objective, nuisance,
                # status and fallback diagnostics; no label or DGP is passed in.
                fits.append(fit)
                log.write(json.dumps(clean(fit), allow_nan=False) + "\n")
                log.flush()
            os.fsync(log.fileno())
            warm_counts = dict(warm.counts) if warm is not None else None
            baseline_error = warm.baseline_error if warm is not None else None
        end_to_end_fit_seconds = time.perf_counter() - started
        fit_seconds = sum(fit["measured_fit_seconds"] for fit in fits) + components.get("baseline_fit_seconds", 0.)
        successful = sum(fit["valid_estimate_se"] for fit in fits)
        rows = [dict(snp_index=fit["snp_index"], parameter=p,
                     estimate=fit.get("estimates", {}).get(p),
                     se=fit.get("standard_errors", {}).get(p),
                     status=fit["status"], reason=fit.get("reason")) for fit in fits for p in PARAMS]
        extra = dict(warm_counts=warm_counts, baseline_error=baseline_error,
                     total_objective=sum(float(f["diagnostics"]["objective"]) for f in fits)
                         if successful == M else None)
    else:
        taus = np.linspace(.02, .98, 45)
        t = time.perf_counter()
        base = truth.fit_baseline(y, "mixture")
        components["baseline_fit_seconds"] = time.perf_counter() - t
        if a.method == "gmoq_pooled":
            t = time.perf_counter()
            q, fhat = core.compute_rif_nuisance(y, taus, density_method="silverman_gaussian", density_floor=1e-8)
            components["rif_nuisance_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            Q, _ = np.linalg.qr(np.ones((n, 1)))
            components["intercept_qr_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            plan = core.make_pooled_bin_plan(y, taus, Q, q_tau=q, fhat=fhat, n_threads=1)
            components["pooled_plan_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            beta, d = core.compute_pooled_binned(G, plan, min_mac=5., min_var=1e-8, n_threads=1)
            beta, d = np.asarray(beta, float), np.asarray(d, float)
            components["pooled_stage1_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            sigma = np.asarray(plan.sigma0(), float)
            se = np.sqrt(np.maximum(np.diag(sigma)[None, :] / d[:, None], 0.))
            components["pooled_covariance_seconds"] = time.perf_counter() - t
            cov = dict(source="pooled", cov_vec=sigma[np.triu_indices(45)],
                       offsets=np.zeros(M, dtype=int), cov_scales=1. / d)
        else:
            from validation.run_examples import balanced_block_ids
            t = time.perf_counter()
            rif, q, fhat = core.compute_rif(y, taus, density_method="silverman_gaussian", density_floor=1e-8)
            components["compute_rif_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            residual, Q = core.residualize(rif, np.ones((n, 1)))
            components["residualize_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            blocks = 1000
            bid = balanced_block_ids(n, blocks, manifest["phenotype_seed"])
            plan = core.make_block_plan(residual, Q, bid, blocks, n_threads=1)
            components["jackknife_plan_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            beta, se, upper = core.compute_block_finalize_panel_fast(G, residual, Q, bid, blocks,
                min_mac=5., min_var=1e-8, n_threads=1, panel_width=16, plan=plan)
            beta, se = np.asarray(beta, float), np.asarray(se, float)
            # n_eff is the genotype residual sum of squares, as in the frozen
            # panel's pooled preliminary fit; estimate/SE/covariance are jackknife.
            gsum = np.sum(G, axis=0, dtype=np.float64)
            d = np.sum(G*G, axis=0, dtype=np.float64) - gsum*gsum/n
            components["jackknife_stage1_seconds"] = time.perf_counter() - t
            L = len(taus) * (len(taus)+1)//2
            cov = dict(source="jackknife", cov_vec=np.asarray(upper, float).reshape(-1),
                       offsets=np.arange(M, dtype=int)*L*4, cov_scales=None)
        if not (np.all(np.isfinite(beta)) and np.all(np.isfinite(se)) and np.all(se > 0) and np.all(d > 0)):
            raise RuntimeError("pooled Stage 1 estimate/SE/denominator check failed")
        t = time.perf_counter()
        payload = dict(n=n, taus=taus, q_tau=q, baseline=base, beta=beta.T,
                       se_tau=se.T, n_eff=d, **cov)
        write_json(work / "stage2_input.json", payload)
        components["stage2_input_write_seconds"] = time.perf_counter() - t
        t = time.perf_counter()
        rscript = pin.get("rscript", "Rscript")
        with (work / "stage2.stdout.log").open("w") as out, (work / "stage2.stderr.log").open("w") as err:
            subprocess.run([rscript, "--vanilla", str(HERE / "paired_stage2.R"), a.pin_json,
                            str(work / "stage2_input.json"), str(work / "stage2_output.json"), "primary", "1"],
                           check=True, stdout=out, stderr=err)
        components["stage2_process_wall_seconds"] = time.perf_counter() - t
        stage2 = read_json(work / "stage2_output.json")
        # Primary fitting time ends before optional dense-kernel validation.
        end_to_end_fit_seconds = time.perf_counter() - started
        fit_seconds = (sum(value for key, value in components.items()
                           if key not in ("stage2_input_write_seconds", "stage2_process_wall_seconds"))
                       + stage2["timings"]["make_weights_seconds"]
                       + stage2["timings"]["first_param_gwas_seconds"])
        rows = []
        for one in stage2["rows"]:
            for name, estimate, error in zip(stage2["parameter_names"], one["estimates"], one["ses"]):
                rows.append(dict(snp_index=one["snp_index"], parameter=normalise_parameter(name),
                                 estimate=estimate, se=error, status="ok", reason=None))
        successful = M
        check = None
        if a.validate_stage1 and a.method == "gmoq_pooled":
            t = time.perf_counter()
            # A one-SNP identity check outside the primary fit timer. Allocation
            # is N x 45 only for this pilot diagnostic, never the production arm.
            rif, _, _ = core.compute_rif(y, taus, density_method="silverman_gaussian", density_floor=1e-8)
            residual, qr = core.residualize(rif, np.ones((n, 1)))
            dense_b, dense_d = core.compute_pooled_panel(G[:, :1], residual, qr,
                min_mac=5., min_var=1e-8, n_threads=1, panel_size=16)
            check = dict(max_abs_beta_diff=float(np.max(np.abs(np.asarray(dense_b) - beta[:1]))),
                         max_abs_d_diff=float(np.max(np.abs(np.asarray(dense_d) - d[:1]))),
                         validation_seconds=time.perf_counter() - t)
            check["passed"] = bool(np.allclose(dense_b, beta[:1], rtol=1e-9, atol=1e-11)
                                   and np.allclose(dense_d, d[:1], rtol=1e-11, atol=1e-9))
            if not check["passed"]:
                raise RuntimeError(f"production pooled kernel differs from dense reference: {check}")
        extra = dict(package=package, stage2_timings=stage2["timings"], stage2_version=stage2["package_version"],
                     stage1_identity_check=check, baseline=base)
    result = dict(schema_version=1, method=a.method, n=n, M=M, status="complete" if successful == M else "failed",
                  attempted_snp_fits=M, successful_snp_fits=successful, failures=M-successful,
                  fit_seconds=fit_seconds, startup_import_seconds=import_ready-ENTRY_TIME,
                  end_to_end_fit_seconds=end_to_end_fit_seconds,
                  data_load_seconds=data_load, child_elapsed_seconds=time.perf_counter()-ENTRY_TIME,
                  components=components, rows=rows, host=platform.node(), python=sys.version,
                  numpy=np.__version__, threads={name: os.environ.get(name) for name in THREAD_ENV}, **extra)
    write_json(work / "result.json", result)
    return 0 if result["status"] == "complete" else 1


def parse_gnu_time(path):
    output = {}
    if not path.exists():
        return output
    for line in path.read_text(errors="replace").splitlines():
        if "Maximum resident set size" in line:
            output["peak_rss_kib"] = int(line.rsplit(":", 1)[1].strip())
        elif "User time (seconds)" in line:
            output["user_seconds"] = float(line.rsplit(":", 1)[1].strip())
        elif "System time (seconds)" in line:
            output["system_seconds"] = float(line.rsplit(":", 1)[1].strip())
    return output


def launch(a, cell, method, manifest):
    arm = cell / method
    arm.mkdir(exist_ok=True)
    # Completed arms are immutable and reusable. Partial attempts are retained;
    # a retry restarts the whole arm, so it has an honest fresh-process total.
    attempts = sorted(arm.glob("attempt[0-9][0-9][0-9]"))
    if attempts:
        last = attempts[-1]
        if (last / "process.json").exists() and (last / "result.json").exists():
            process, result = read_json(last / "process.json"), read_json(last / "result.json")
            if process["return_code"] == 0 and result["status"] == "complete":
                return last, result, process
    work = arm / f"attempt{len(attempts)+1:03d}"
    work.mkdir()
    pin = manifest["source"]["pin"]
    python = pin.get("stage1_python", sys.executable)
    cmd = [python, str(HERE / "paired_speed.py"), "_child", "--paper-root", a.paper_root,
           "--pin-json", a.pin_json, "--method", method, "--input-dir", str(cell / "inputs"), "--workdir", str(work)]
    if a.validate_stage1:
        cmd.append("--validate-stage1")
    env = dict(os.environ)
    env.update({name: "1" for name in THREAD_ENV})
    timebin = Path("/usr/bin/time")
    use_time = platform.system() == "Linux" and timebin.exists()
    if use_time:
        cmd = [str(timebin), "-v", "-o", str(work / "time-v.txt"), *cmd]
    started = time.perf_counter()
    write_json(work / "launch.json", dict(command=cmd, host=platform.node(), slurm_job=os.getenv("SLURM_JOB_ID"),
                                         input_sha256={k: v["sha256"] for k, v in manifest["arrays"].items()}))
    with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
        proc = subprocess.Popen(cmd, env=env, stdout=out, stderr=err, start_new_session=True)
        timeout = False
        try:
            proc.wait(timeout=a.timeout_seconds)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            timeout = True
    process = dict(return_code=proc.returncode, timeout=timeout, process_wall_seconds=time.perf_counter()-started,
                   **parse_gnu_time(work / "time-v.txt"))
    write_json(work / "process.json", process)
    result = read_json(work / "result.json") if (work / "result.json").exists() else dict(
        method=method, n=a.n, M=a.variants, status="timeout" if timeout else "process_failure",
        attempted_snp_fits=a.variants, successful_snp_fits=0, failures=a.variants, rows=[])
    return work, result, process


def quality(a, cell, arms):
    canonical = arms.get("mixture_ml")
    warm = arms.get("mixture_ml_warm")
    rows = []
    complete = all(result["status"] == "complete" and process["return_code"] == 0
                   for _, result, process in arms.values()) and set(arms) == set(METHODS)
    if canonical and warm and canonical[1]["status"] == warm[1]["status"] == "complete":
        def fits(work):
            return [json.loads(line) for line in (work / "variant_fits.jsonl").read_text().splitlines()]
        cfits, wfits = fits(canonical[0]), fits(warm[0])
        if len(cfits) != a.variants or len(wfits) != a.variants:
            complete = False
        for c, w in zip(cfits, wfits):
            ce, we = c["estimates"], w["estimates"]
            cs, ws = c["standard_errors"], w["standard_errors"]
            cd, wd = c["diagnostics"], w["diagnostics"]
            objective_delta = float(wd["objective"]) - float(cd["objective"])
            est_delta = max(abs(we[p]-ce[p]) / cs[p] for p in PARAMS)
            se_delta = max(abs(ws[p]/cs[p]-1.) for p in PARAMS)
            nstarts_ok = cd.get("n_starts") == 6 and wd.get("n_starts") == (6 if wd.get("warm_fallback") else 1)
            passed = (fit_valid(c) and fit_valid(w) and c["snp_index"] == w["snp_index"] and nstarts_ok
                      and abs(objective_delta)/a.n <= a.loglik_per_obs_tolerance
                      and est_delta <= a.estimate_se_tolerance and se_delta <= a.se_relative_tolerance)
            rows.append(dict(snp_index=c["snp_index"], canonical_status=c["status"], warm_status=w["status"],
                             canonical_loglik=-cd["objective"], warm_loglik=-wd["objective"],
                             warm_minus_canonical_nll=objective_delta, abs_nll_delta_per_obs=abs(objective_delta)/a.n,
                             max_estimate_difference_in_canonical_se=est_delta, max_relative_se_difference=se_delta,
                             warm_fallback=wd.get("warm_fallback"), canonical_n_starts=cd.get("n_starts"),
                             warm_n_starts=wd.get("n_starts"), passed=passed))
    gate = dict(status="passed" if complete and len(rows) == a.variants and all(r["passed"] for r in rows) else "failed",
                n=a.n, M=a.variants, rep=a.rep, compared_variants=len(rows),
                passed_variants=sum(r["passed"] for r in rows), all_method_arms_complete=complete,
                thresholds=dict(abs_loglik_delta_per_obs=a.loglik_per_obs_tolerance,
                                max_estimate_delta_in_canonical_se=a.estimate_se_tolerance,
                                max_relative_se_delta=a.se_relative_tolerance),
                note="Warm-start sensitivity is timing-eligible only after the same fits, SEs and likelihoods pass these checks. Null panels assess cost and finite-SE fitting, not calibration.")
    csv_write(cell / "warm_quality.csv", rows)
    write_json(cell / "quality.json", gate)
    return gate


def reuse_stage2(a, cell, arms):
    work, result, _ = arms["gmoq_pooled"]
    if result["status"] != "complete":
        return None
    target = cell / "stage2_reuse.json"
    if target.exists():
        return read_json(target)
    pin = read_json(a.pin_json)
    cmd = [pin.get("rscript", "Rscript"), "--vanilla", str(HERE / "paired_stage2.R"), a.pin_json,
           str(work / "stage2_input.json"), str(target), "reuse", str(a.reuse_reps)]
    env = dict(os.environ, **{name: "1" for name in THREAD_ENV})
    started = time.perf_counter()
    with (cell / "stage2_reuse.stdout.log").open("w") as out, (cell / "stage2_reuse.stderr.log").open("w") as err:
        subprocess.run(cmd, env=env, check=True, stdout=out, stderr=err, timeout=a.timeout_seconds)
    reuse = read_json(target)
    reuse["process_wall_seconds"] = time.perf_counter()-started
    write_json(target, reuse)
    return reuse


def run_cell(a):
    if a.n not in NS and not a.smoke:
        raise ValueError("noncanonical N requires --smoke")
    if a.n < 100 or a.rep not in (1, 2, 3) or not 1 <= a.variants <= 100:
        raise ValueError("require N >= 100, rep in 1..3 and variants in 1..100")
    cell = Path(a.out).resolve() / f"n{a.n}_rep{a.rep}"
    cell.mkdir(parents=True, exist_ok=True)
    identity = source_identity(a)
    manifest = input_panel(a, cell / "inputs", identity)
    methods = a.methods.split(",")
    if len(set(methods)) != len(methods) or any(m not in METHODS for m in methods):
        raise ValueError("invalid or repeated method")
    # Rotate method order across cells and persist the actual order.
    rotation = ((NS.index(a.n) if a.n in NS else 0) + a.rep - 1) % len(methods)
    order = methods[rotation:] + methods[:rotation]
    write_json(cell / "cell.json", dict(n=a.n, rep=a.rep, M=a.variants, method_order=order,
                                        host=platform.node(), slurm_job=os.getenv("SLURM_JOB_ID")))
    arms = {}
    for method in METHODS:
        for work in sorted((cell / method).glob("attempt[0-9][0-9][0-9]"), reverse=True):
            if (work / "process.json").exists() and (work / "result.json").exists():
                process, result = read_json(work / "process.json"), read_json(work / "result.json")
                if process["return_code"] == 0 and result["status"] == "complete":
                    arms[method] = (work, result, process)
                    break
    for method in order:
        print(f"{time.strftime('%H:%M:%S')} N={a.n} rep={a.rep} M={a.variants} {method}", flush=True)
        arms[method] = launch(a, cell, method, manifest)
        write_json(cell / "progress.json", dict(status="running", completed_methods=list(arms),
                                                last_method=method, last_status=arms[method][1]["status"]))
    gate = quality(a, cell, arms)
    rows, timings = [], []
    for method, (work, result, process) in arms.items():
        startup = result.get("startup_import_seconds")
        # Parent process wall includes interpreter startup, imports, pin checks,
        # data loading, fitting, output/checkpoint IO and the R child if present.
        # Optional pilot-only validation is recorded separately, so its wall is
        # deliberately not eligible as a production full-process observation.
        timings.append(dict(n=a.n, rep=a.rep, M=a.variants, method=method, status=result["status"],
                            fit_seconds=result.get("fit_seconds"), process_wall_seconds=process["process_wall_seconds"],
                            end_to_end_fit_seconds=result.get("end_to_end_fit_seconds"),
                            startup_import_seconds=startup, data_load_seconds=result.get("data_load_seconds"),
                            peak_rss_kib=process.get("peak_rss_kib"), return_code=process["return_code"],
                            successful_snp_fits=result["successful_snp_fits"], attempted_snp_fits=a.variants,
                            quality_status=gate["status"], primary_wall_eligible=not a.validate_stage1,
                            attempt=work.name, log_dir=str(work)))
        for row in result.get("rows", []):
            rows.append(dict(n=a.n, rep=a.rep, M=a.variants, method=method, **row))
    csv_write(cell / "timings.csv", timings)
    csv_write(cell / "estimates.csv", rows)
    reuse = reuse_stage2(a, cell, arms) if "gmoq_pooled" in arms else None
    done = dict(status="complete" if gate["status"] == "passed" else "quality_failure", n=a.n, M=a.variants,
                rep=a.rep, quality=gate["status"], stage2_reuse=bool(reuse), methods=list(arms))
    write_json(cell / "COMPLETE.json", done)
    write_json(cell / "progress.json", done)
    print(json.dumps(done), flush=True)
    return 0 if done["status"] == "complete" else 1


def collate(a):
    root = Path(a.out).resolve()
    expected = [(n, rep) for n in NS for rep in range(1, a.reps+1)]
    all_timings, all_estimates, all_quality, reused, missing, bad = [], [], [], [], [], []
    for n, rep in expected:
        cell = root / f"n{n}_rep{rep}"
        if not (cell / "COMPLETE.json").exists():
            missing.append(str(cell))
            continue
        done = read_json(cell / "COMPLETE.json")
        if done["status"] != "complete" or done["M"] != a.variants:
            bad.append(dict(cell=str(cell), completion=done))
        for filename, target in (("timings.csv", all_timings), ("estimates.csv", all_estimates),
                                 ("warm_quality.csv", all_quality)):
            with (cell / filename).open(newline="") as handle:
                for row in csv.DictReader(handle):
                    if int(row.get("n", n)) != n or int(row.get("rep", rep)) != rep:
                        raise ValueError(f"wrong cell identifiers in {cell / filename}")
                    target.append({**row, "n": n, "rep": rep})
        if (cell / "stage2_reuse.json").exists():
            row = read_json(cell / "stage2_reuse.json")
            if int(row["n"]) != n:
                raise ValueError(f"wrong sample size in {cell / 'stage2_reuse.json'}")
            reused.append({**row, "n": n, "rep": rep})
    csv_write(root / "timings.csv", all_timings)
    csv_write(root / "estimates.csv", all_estimates)
    csv_write(root / "warm_quality.csv", all_quality)
    write_json(root / "stage2_reuse.json", reused)
    report = dict(status="complete" if not missing and not bad else "incomplete", expected_cells=len(expected),
                  available_cells=len(expected)-len(missing), variants_per_cell=a.variants, missing=missing, invalid=bad,
                  timing_rows=len(all_timings), estimate_rows=len(all_estimates), quality_rows=len(all_quality))
    write_json(root / "COMPLETE.json", report)
    print(json.dumps(report))
    return 0 if report["status"] == "complete" else 1


def parser():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    cell = sub.add_parser("cell", help="run or resume one paired N x seed cell")
    for name in ("paper-root", "pin-json", "out"):
        cell.add_argument("--"+name, required=True)
    cell.add_argument("--n", type=int, required=True)
    cell.add_argument("--rep", type=int, required=True)
    cell.add_argument("--variants", type=int, default=100)
    cell.add_argument("--methods", default=",".join(METHODS))
    cell.add_argument("--timeout-seconds", type=float, default=86400)
    cell.add_argument("--reuse-reps", type=int, default=20)
    cell.add_argument("--smoke", action="store_true", help="allow small noncanonical N for an end-to-end check")
    cell.add_argument("--validate-stage1", action="store_true", help="pilot-only dense identity check outside fit timer")
    cell.add_argument("--loglik-per-obs-tolerance", type=float, default=1e-8)
    cell.add_argument("--estimate-se-tolerance", type=float, default=.01)
    cell.add_argument("--se-relative-tolerance", type=float, default=.01)
    proc = sub.add_parser("_child", help=argparse.SUPPRESS)
    for name in ("paper-root", "pin-json", "method", "input-dir", "workdir"):
        proc.add_argument("--"+name, required=True)
    proc.add_argument("--validate-stage1", action="store_true")
    combine = sub.add_parser("collate", help="check and combine all expected cells")
    combine.add_argument("--out", required=True)
    combine.add_argument("--reps", type=int, choices=(1, 3), default=3)
    combine.add_argument("--variants", type=int, default=100)
    return ap


def main():
    a = parser().parse_args()
    if a.command == "_child":
        try:
            return child(a)
        except Exception as exc:
            write_json(Path(a.workdir) / "failure.json", dict(status="exception", reason=f"{type(exc).__name__}: {exc}"))
            traceback.print_exc()
            return 1
    return run_cell(a) if a.command == "cell" else collate(a)


if __name__ == "__main__":
    raise SystemExit(main())
