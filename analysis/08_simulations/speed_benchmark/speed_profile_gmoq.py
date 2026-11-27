#!/usr/bin/env python
"""aux_speed_profile, part 1 (Python), for a released gmoq version: per-SNP costs for Fig 2C / Fig S5d.

Same panels as speed_profile.py (the frozen benchmark's null two-Normal mixture y and
binomial(2, .3) genotypes, same seeds), same single-thread settings, same components,
written with the same component names so collate.R is unchanged in kind:

  g_env_core_load          import of the pinned gmoq Stage 1 core (pin checked first)
  g_baseline_mixture_fit   two-Normal baseline mixture ML fit to y (validation.truth)
  g_compute_rif_total      core.compute_rif (quantiles + KDE + RIF)
  g_residualize            core.residualize on the intercept
  g_stage1_plan            the genotype-independent Stage 1 plan (one-off per phenotype)
  g_stage1_ols_panel       Stage 1 pooled OLS for M = 100 and M = 300 SNPs (production
                           pooled kernel, core.compute_pooled_binned); the per-SNP cost is
                           the difference / 200, as before
  g_sigma0                 pooled Sigma_0 and SEs
  g_write_json             Stage 2 input JSON (the payload stage2_profile.R reads)
  ml_fit_one_snp           direct mixture ML per SNP (validation.comparators.fit_mixture_ml)

The frozen benchmark's Stage 2 child process (mixture_comparison.sh) is not run: Stage 2
proper is timed in a warm R session by stage2_profile.R against the same gmoq pin.
Simulated data only; aggregate timings to <out>/python_timings.csv.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib, json, os, subprocess, sys, time
from pathlib import Path
import numpy as np

REPO = Path(os.environ.get("PAPER_REPO", Path(__file__).resolve().parents[3]))
sys.path.insert(0, str(REPO))
from validation import truth  # noqa: E402
from validation.comparators import fit_mixture_ml  # noqa: E402

TAUS = np.linspace(.02, .98, 45)
NULL_CASE = {"family": "mixture", "maf": .3, "pi0": .5, "mu1": -2., "mu2": 2., "sd1": 1., "sd2": 1.,
             "gamma": 0., "delta1": 0., "delta2": 0.}


def safe(x):
    if isinstance(x, dict): return {str(k): safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [safe(v) for v in x]
    if isinstance(x, np.ndarray): return safe(x.tolist())
    if isinstance(x, np.integer): return int(x)
    if isinstance(x, np.floating): return float(x)
    return x


def panel(n, seed, M):
    """Identical to speed_profile.py / the frozen benchmark: null mixture y and binomial(2, .3) G."""
    y = truth.draw(NULL_CASE, n, seed)["y"]
    G = np.random.default_rng(seed + 7919).binomial(2, .3, size=(n, M)).astype(np.int8)
    return np.asarray(G, float), np.asarray(y, float)


def load_core(root: Path, commit: str, so_sha256: str):
    """Import gmoq Stage 1 from the pinned checkout, refusing any other commit or build."""
    head = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], check=True,
                          capture_output=True, text=True).stdout.strip()
    if head != commit:
        raise SystemExit(f"gmoq checkout {root} is at {head}, not the pinned {commit}")
    dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
                           check=True, capture_output=True, text=True).stdout.strip()
    if dirty:
        raise SystemExit(f"gmoq checkout {root} has tracked changes")
    sys.path.insert(0, str(root / "stage1"))
    core = importlib.import_module("gmoq.core")
    so = Path(core.cpp_ext.__file__)
    got = hashlib.sha256(so.read_bytes()).hexdigest()
    if got != so_sha256:
        raise SystemExit(f"Stage 1 extension {so} sha256 {got} != pinned {so_sha256}")
    return core


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--pin-json", required=True, help="PIN.json of the gmoq install (commit, stage1_root, extension sha256)")
    ap.add_argument("--n", default="10000,50000,400000")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--ml-snps", default="10,10,5")
    ap.add_argument("--stage1-m-extra", type=int, default=300)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=False)
    pin = json.loads(Path(a.pin_json).read_text())
    ns = [int(v) for v in a.n.split(",")]; mls = [int(v) for v in a.ml_snps.split(",")]
    rows = []

    def rec(n, rep, component, seconds, m=None, note=""):
        rows.append({"n": n, "rep": rep, "component": component, "M": m, "seconds": seconds, "note": note})
        print(f"[{time.strftime('%H:%M:%S')}] n={n} rep={rep} {component} M={m} {seconds:.4f}s", flush=True)

    t = time.perf_counter()
    core = load_core(Path(pin["stage1_root"]), pin["commit"], pin["stage1_extension_sha256"])
    t_load = time.perf_counter() - t
    if core.cpp_ext.openmp_info()["max_threads"] != 1:
        raise SystemExit("expected a single OpenMP thread (OMP_NUM_THREADS=1)")

    for n, mlk in zip(ns, mls):
        for i in range(a.reps):
            seed = 97000000 + i; M = 100
            G, y = panel(n, seed, M)
            t_all = time.perf_counter()
            rec(n, i + 1, "g_env_core_load", t_load if (n == ns[0] and i == 0) else 0.0, note="once per process")
            t = time.perf_counter(); base = truth.fit_baseline(y, "mixture"); rec(n, i + 1, "g_baseline_mixture_fit", time.perf_counter() - t)
            t = time.perf_counter(); rif, q, fhat = core.compute_rif(y, TAUS, density_method="silverman_gaussian", density_floor=1e-8); rec(n, i + 1, "g_compute_rif_total", time.perf_counter() - t)
            t = time.perf_counter(); Qrif, Q = core.residualize(rif, np.ones((len(y), 1))); rec(n, i + 1, "g_residualize", time.perf_counter() - t)
            t = time.perf_counter(); plan = core.make_pooled_bin_plan(y, TAUS, Q, q_tau=q, fhat=fhat, n_threads=1); rec(n, i + 1, "g_stage1_plan", time.perf_counter() - t)
            t = time.perf_counter(); beta, d = core.compute_pooled_binned(G, plan, min_mac=5., min_var=1e-8, n_threads=1)
            beta = np.asarray(beta, float); d = np.asarray(d, float); rec(n, i + 1, "g_stage1_ols_panel", time.perf_counter() - t, M)
            t = time.perf_counter(); sigma = np.asarray(plan.sigma0(), float); se = np.sqrt(np.maximum(np.diag(sigma)[None, :] / d[:, None], 0.)); rec(n, i + 1, "g_sigma0", time.perf_counter() - t)
            cov = {"cov_vec": sigma[np.triu_indices(45)], "offsets": np.zeros(M, dtype=int), "cov_scales": 1. / d}
            t = time.perf_counter(); inp = out / f"stage2_input_n{n}_rep{i + 1}.json"
            inp.write_text(json.dumps(safe({"n": len(y), "taus": TAUS, "q_tau": q, "baseline": base, "beta": beta.T, "se_tau": se.T, "n_eff": d, **cov}), indent=2, allow_nan=False) + "\n")
            rec(n, i + 1, "g_write_json", time.perf_counter() - t)
            rec(n, i + 1, "g_total_reconstructed", time.perf_counter() - t_all, M, "Python side only; no Stage 2 child")
            if i == 0:
                # Identity check against the dense pooled kernel on the same panel (not a time).
                bd, dd = core.compute_pooled_panel(G, Qrif, Q, min_mac=5., min_var=1e-8, n_threads=1, panel_size=16)
                rel = float(np.nanmax(np.abs(np.asarray(bd) - beta)) / np.nanstd(np.asarray(bd)))
                rec(n, 1, "check_max_rel_diff_binned_vs_dense", rel, M, "not a time: |beta diff| / sd(beta)")
                Mx = a.stage1_m_extra; Gx, _ = panel(n, seed, Mx)
                t = time.perf_counter(); core.compute_pooled_binned(Gx, plan, min_mac=5., min_var=1e-8, n_threads=1); rec(n, 1, "g_stage1_ols_panel", time.perf_counter() - t, Mx)
                t = time.perf_counter(); core.compute_pooled_panel(Gx, Qrif, Q, min_mac=5., min_var=1e-8, n_threads=1, panel_size=16); rec(n, 1, "g_stage1_ols_panel_dense", time.perf_counter() - t, Mx, "diagnostic")
                t = time.perf_counter(); core.compute_pooled_panel(G, Qrif, Q, min_mac=5., min_var=1e-8, n_threads=1, panel_size=16); rec(n, 1, "g_stage1_ols_panel_dense", time.perf_counter() - t, M, "diagnostic")
                del Gx
                for j in range(mlk):
                    t = time.perf_counter(); fit = fit_mixture_ml(G[:, j], y, {"family": "mixture", "maf": .3})
                    rec(n, 1, "ml_fit_one_snp", time.perf_counter() - t, 1, str(fit["status"]))
            del G, y, rif, Qrif, plan
    with (out / "python_timings.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["n", "rep", "component", "M", "seconds", "note"]); w.writeheader(); w.writerows(rows)
    (out / "python_env.json").write_text(json.dumps({"host": os.uname().nodename, "slurm_job": os.environ.get("SLURM_JOB_ID"),
        "omp": os.environ.get("OMP_NUM_THREADS"), "openblas": os.environ.get("OPENBLAS_NUM_THREADS"), "python": sys.version,
        "numpy": np.__version__, "package": "gmoq", "package_commit": pin["commit"],
        "stage1_extension_sha256": pin["stage1_extension_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
