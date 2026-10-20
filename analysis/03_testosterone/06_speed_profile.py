#!/usr/bin/env python3
# Cost of G-MoQ and of direct mixture maximum likelihood, part 1: Stage 1 and direct ML
#
# What: times the parts of the two-component mixture benchmark separately, on one core.
#   It regenerates the benchmark's simulated panels (a null two-Normal mixture phenotype and
#   100 binomial(2, 0.3) genotypes, at N = 10,000, 50,000 and 400,000, three repeats, same seeds
#   and generators as the benchmark), and times each G-MoQ component: the baseline mixture fit,
#   the recentred influence function (quantiles and kernel density), residualisation, Stage 1
#   OLS for 100 and for 300 variants (their difference gives the marginal cost per variant),
#   the pooled covariance and the benchmark's own Stage 2 step. It also times direct mixture
#   maximum likelihood, one variant at a time, on the first 10, 10 and 5 variants. The Stage 2
#   input of each panel is kept for part 2 (06_speed_stage2.R). Simulated data only.
# Paper: Fig. 2c; Supplementary Fig. S5d
# Input: the mixture benchmark code (simulation truth, comparators, pinned Stage 1 loader) and
#   its validated software environment
# Output: <out>/python_timings.csv (n, rep, component, M, seconds, note) and
#   <out>/stage2_input_n<N>_rep<r>.json
# Runs: on a cluster node, single-threaded (06_speed.sh), in about 1 hour

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

TAUS = np.linspace(.02, .98, 45)
NULL_CASE = {"family": "mixture", "maf": .3, "pi0": .5, "mu1": -2., "mu2": 2., "sd1": 1., "sd2": 1.,
             "gamma": 0., "delta1": 0., "delta2": 0.}


def json_safe(x):
    if isinstance(x, dict): return {str(k): json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [json_safe(v) for v in x]
    if isinstance(x, np.ndarray): return json_safe(x.tolist())
    if isinstance(x, np.integer): return int(x)
    if isinstance(x, np.floating): return float(x)
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--benchmark-root", required=True)     # checkout holding validation/
    ap.add_argument("--environment", required=True)        # validated environment directory
    ap.add_argument("--rscript", required=True)
    ap.add_argument("--n", default="10000,50000,400000")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--ml-snps", default="10,10,5")
    ap.add_argument("--stage1-m-extra", type=int, default=300)
    a = ap.parse_args()

    bench = Path(a.benchmark_root); env_dir = Path(a.environment)
    sys.path.insert(0, str(bench))
    from validation import truth
    from validation.comparators import fit_mixture_ml
    from validation.run_examples import (PINNED_PACKAGE_COMMIT, PINNED_PACKAGE_ROOT,
                                         load_environment_manifest, load_pinned_core)
    stage2_script = bench / "validation/figures/mixture_comparison.sh"

    def panel(n, seed, M):
        """The benchmark's panel: null mixture phenotype y and binomial(2, 0.3) genotypes G."""
        y = truth.draw(NULL_CASE, n, seed)["y"]
        G = np.random.default_rng(seed + 7919).binomial(2, .3, size=(n, M)).astype(np.int8)
        return np.asarray(G, float), np.asarray(y, float)

    out = Path(a.out); out.mkdir(parents=True, exist_ok=False)
    ns = [int(v) for v in a.n.split(",")]; mls = [int(v) for v in a.ml_snps.split(",")]
    rows = []

    def rec(n, rep, component, seconds, m=None, note=""):
        rows.append({"n": n, "rep": rep, "component": component, "M": m, "seconds": seconds, "note": note})
        print(f"n={n} rep={rep} {component} M={m} {seconds:.4f}s", flush=True)

    # 1. Fixed process costs of the benchmark's Stage 2 child (shell, R start-up), once.
    for label, cmd in [
        ("child_bash_activate_only", ["bash", "-c", f"source {env_dir / 'activate.sh'} >/dev/null 2>&1; true"]),
        ("rscript_startup_vanilla", [a.rscript, "--vanilla", "-e", "invisible(NULL)"]),
        ("rscript_startup_plus_library", [a.rscript, "--vanilla", "-e",
                                          f".libPaths(c('{env_dir / 'Rlib'}', .libPaths())); suppressMessages(library(gmoq)); invisible(jsonlite::fromJSON('[1]'))"]),
    ]:
        for rep in range(3):
            t = time.perf_counter(); subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            rec(0, rep + 1, label, time.perf_counter() - t)

    # 2. Each G-MoQ component on each panel; direct ML on a few variants of the first panel.
    for n, n_ml in zip(ns, mls):
        for i in range(a.reps):
            seed = 97000000 + i; M = 100
            G, y = panel(n, seed, M)
            t_all = time.perf_counter()
            t = time.perf_counter(); env = load_environment_manifest(env_dir / "environment.json", Path(PINNED_PACKAGE_ROOT)); core, prov, sigma0fn = load_pinned_core(Path(PINNED_PACKAGE_ROOT), env); rec(n, i + 1, "g_env_core_load", time.perf_counter() - t)
            t = time.perf_counter(); base = truth.fit_baseline(y, "mixture"); rec(n, i + 1, "g_baseline_mixture_fit", time.perf_counter() - t)
            t = time.perf_counter(); rif, q, _ = core.compute_rif(y, TAUS, density_method="silverman_gaussian", density_floor=1e-8); rec(n, i + 1, "g_compute_rif_total", time.perf_counter() - t)
            t = time.perf_counter(); Qrif, Q = core.residualize(rif, np.ones((len(y), 1))); rec(n, i + 1, "g_residualize", time.perf_counter() - t)
            t = time.perf_counter(); beta, d = core.compute_pooled_panel(G, Qrif, Q, min_mac=5., min_var=1e-8, n_threads=1, panel_size=16); beta = np.asarray(beta, float); d = np.asarray(d, float); rec(n, i + 1, "g_stage1_ols_panel", time.perf_counter() - t, M)
            t = time.perf_counter(); sigma = np.asarray(sigma0fn(Qrif), float); se = np.sqrt(np.maximum(np.diag(sigma)[None, :] / d[:, None], 0.)); rec(n, i + 1, "g_sigma0", time.perf_counter() - t)
            cov = {"cov_vec": sigma[np.triu_indices(45)], "offsets": np.zeros(M, dtype=int), "cov_scales": 1. / d}
            t = time.perf_counter(); inp = out / f"stage2_input_n{n}_rep{i + 1}.json"
            inp.write_text(json.dumps(json_safe({"n": len(y), "taus": TAUS, "q_tau": q, "baseline": base, "beta": beta.T,
                                                 "se_tau": se.T, "n_eff": d, **cov}), indent=2, allow_nan=False) + "\n")
            rec(n, i + 1, "g_write_json", time.perf_counter() - t)
            outp = out / f"stage2_output_n{n}_rep{i + 1}.json"
            cmd = ["bash", str(stage2_script), "--panel-stage2", "--input", str(inp), "--output", str(outp), "--source", "pooled",
                   "--package-root", str(PINNED_PACKAGE_ROOT), "--package-commit", PINNED_PACKAGE_COMMIT, "--r-lib", env["r_lib"]]
            t = time.perf_counter(); subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1100); rec(n, i + 1, "g_stage2_child_process", time.perf_counter() - t, M)
            rec(n, i + 1, "g_total_reconstructed", time.perf_counter() - t_all, M)
            if i == 0:
                # Quantiles alone and the kernel density alone (the two halves of compute_rif).
                t = time.perf_counter(); qq = np.quantile(y, TAUS, method='median_unbiased'); rec(n, 1, "g_quantiles_only", time.perf_counter() - t)
                t = time.perf_counter(); h = 1.06 * np.std(y) * len(y) ** (-0.2)
                fh = np.array([np.mean(np.exp(-0.5 * ((y - v) / h) ** 2)) / (h * np.sqrt(2 * np.pi)) for v in qq]); rec(n, 1, "g_kde_only", time.perf_counter() - t)
                # A larger Stage 1 panel: the difference from M = 100 is the cost per variant.
                Mx = a.stage1_m_extra; Gx, _ = panel(n, seed, Mx)
                t = time.perf_counter(); core.compute_pooled_panel(Gx, Qrif, Q, min_mac=5., min_var=1e-8, n_threads=1, panel_size=16); rec(n, 1, "g_stage1_ols_panel", time.perf_counter() - t, Mx)
                del Gx
                # Direct mixture maximum likelihood, the benchmark's call, one variant at a time.
                for j in range(n_ml):
                    t = time.perf_counter(); fit = fit_mixture_ml(G[:, j], y, {"family": "mixture", "maf": .3})
                    rec(n, 1, "ml_fit_one_snp", time.perf_counter() - t, 1, str(fit["status"]))
            del G, y, rif, Qrif
    with (out / "python_timings.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["n", "rep", "component", "M", "seconds", "note"]); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
