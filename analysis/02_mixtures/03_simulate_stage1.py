#!/usr/bin/env python3
# Weight-design simulation, part 1: simulate phenotypes and run Stage 1
#
# What: draws synthetic phenotypes from two-component skew-normal mixtures (and one shifted
#   lognormal, BMI-like control) with known per-allele variant effects on membership and on each
#   component's location, and runs G-MoQ Stage 1 on them with the gmoq-stage1 package
#   (recentred influence function at 45 quantile levels, Silverman Gaussian kernel density,
#   pooled covariance). Part 2 (03_stage2_designs.R) fits the mixture and compares the model
#   and hybrid weight designs on the same Stage 1 output. Synthetic data only.
#   gamma is the per-allele effect on logit P(lower component), with centred genotype; the
#   paper reports upper-component membership, -gamma.
# Paper: Supplementary Note 2 (Supplementary Note Table 1); Supplementary Fig. S7b
# Input: scenario definitions (03_weight_design_scenarios.json)
# Output: per replicate, in <out>/<scenario>/rep_NNNN/: y.f64 (phenotype), stage1.tsv
#   (d_resid and beta_tau_* per variant), sigma0.tsv (pooled covariance), nuis.json (quantile
#   levels, baseline quantiles, kernel density); consumed and deleted by part 2
# Runs: on a cluster (03_weight_design_simulation.sh); about 10 s per replicate at N = 100,000

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np


def load_stage1(package_root):
    """Import the pinned gmoq-stage1 package from its checkout."""
    sys.path.insert(0, str(Path(package_root) / "stage1"))
    from gmoq import core
    from gmoq.pooled_cov import compute_sigma0
    import gmoq
    core.require_cpp_extension("weight-design simulation")
    return core, compute_sigma0, str(Path(gmoq.__file__).resolve())


def rsn(rng, n, xi, omega, alpha):
    """Skew-normal draws by the additive representation."""
    delta = alpha / np.sqrt(1.0 + alpha * alpha)
    u0 = np.abs(rng.standard_normal(n))
    u1 = rng.standard_normal(n)
    return xi + omega * (delta * u0 + np.sqrt(1.0 - delta * delta) * u1)


def draw_component(rng, parts, n):
    """A component may itself be a small mixture of skew-normal parts (to misfit the gap)."""
    w = np.array([p["w"] for p in parts], dtype=float)
    w = w / w.sum()
    which = rng.choice(len(parts), size=n, p=w)
    out = np.empty(n)
    for j, p in enumerate(parts):
        m = which == j
        out[m] = rsn(rng, int(m.sum()), p["xi"], p["omega"], p["alpha"])
    return out


def simulate(scenario, variants, n, maf, rng):
    """Genotypes G (n x K) and phenotype y with the scenario's variant effects."""
    G = rng.binomial(2, maf, size=(n, len(variants))).astype(np.float64)
    Gc = G - 2.0 * maf
    if scenario["family"] == "mixture":
        # Membership of the lower component: logit P(lower) = logit(p_lower) + gamma' g_centred.
        eta = np.log(scenario["p_lower"] / (1 - scenario["p_lower"])) + \
            Gc @ np.array([v.get("gamma", 0.0) for v in variants])
        lower = rng.random(n) < 1.0 / (1.0 + np.exp(-eta))
        y = np.empty(n)
        n_lower = int(lower.sum())
        y[lower] = draw_component(rng, scenario["components"]["lower"], n_lower)
        y[~lower] = draw_component(rng, scenario["components"]["upper"], n - n_lower)
        # Location effects shift every value of the component the participant belongs to.
        shift_lower = G @ np.array([v.get("beta_left", 0.0) for v in variants])
        shift_upper = G @ np.array([v.get("beta_right", 0.0) for v in variants])
        y = y + np.where(lower, shift_lower, shift_upper)
        lower_share = float(lower.mean())
    elif scenario["family"] == "lognormal":
        ln = scenario["lognormal"]
        mult = G @ np.array([v.get("multiplicative", 0.0) for v in variants])
        add = G @ np.array([v.get("additive", 0.0) for v in variants])
        y = ln["shift"] + np.exp(ln["meanlog"] + ln["sdlog"] * rng.standard_normal(n) + mult) + add
        lower_share = None
    else:
        raise SystemExit(f"unknown family {scenario['family']}")
    return y, G, lower_share


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--scenario", required=True)
    ap.add_argument("--rep-start", type=int, required=True)
    ap.add_argument("--rep-end", type=int, required=True)   # exclusive
    ap.add_argument("--out", required=True)
    ap.add_argument("--package-root", required=True)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    scenario = cfg["scenarios"][a.scenario]
    variants = cfg["snp_sets"][scenario["snp_set"]]
    n = int(scenario.get("n", cfg["n"]))
    maf = float(cfg["maf"])
    taus = np.linspace(cfg["tau_from"], cfg["tau_to"], cfg["n_tau"])
    core, compute_sigma0, package_file = load_stage1(a.package_root)
    outdir = Path(a.out) / a.scenario
    outdir.mkdir(parents=True, exist_ok=True)

    for rep in range(a.rep_start, a.rep_end):
        t0 = time.time()
        # One seed per scenario and replicate, so replicates do not depend on how they are batched.
        seed = int(cfg["seed"]) + 1_000_003 * int(scenario["seed_index"]) + rep
        rng = np.random.default_rng(seed)
        y, G, lower_share = simulate(scenario, variants, n, maf, rng)

        # Stage 1: RIF at each quantile level, residualised on an intercept, then per-variant OLS.
        rif, q_tau, fhat = core.compute_rif(y, taus)
        X = np.ones((n, 1))
        rif_resid, Q = core.residualize(rif, X)
        core.check_orthonormal_basis(Q, context="weight-design simulation")
        betas, d = core.compute_pooled_panel(G, rif_resid, Q, min_mac=5.0, min_var=1e-8,
                                             n_threads=a.threads, panel_size=16)
        betas = np.asarray(betas, dtype=np.float64)
        d = np.asarray(d, dtype=np.float64)
        # Check against a direct RIF-OLS regression.
        g_resid = G - Q @ (Q.T @ G)
        d_direct = np.sum(g_resid ** 2, axis=0)
        b_direct = (g_resid.T @ rif_resid) / d_direct[:, None]
        err = float(np.max(np.abs(betas - b_direct) / (1.0 + np.abs(b_direct))))
        if not np.isfinite(err) or err > 1e-10:
            raise SystemExit(f"Stage 1 disagrees with direct RIF-OLS: {err}")
        sigma0 = np.asarray(compute_sigma0(rif_resid), dtype=np.float64)

        rd = outdir / f"rep_{rep:04d}"
        rd.mkdir(exist_ok=True)
        y.astype("<f8").tofile(rd / "y.f64")
        with open(rd / "stage1.tsv", "w") as fh:
            fh.write("snp\td_resid\t" + "\t".join(f"beta_tau_{t:.3f}" for t in taus) + "\n")
            for j, v in enumerate(variants):
                fh.write(v["id"] + "\t" + repr(float(d[j])) + "\t" +
                         "\t".join(repr(float(b)) for b in betas[j]) + "\n")
        np.savetxt(rd / "sigma0.tsv", sigma0, delimiter="\t", fmt="%.17g")
        (rd / "nuis.json").write_text(json.dumps({
            "taus": taus.tolist(), "q_tau": np.asarray(q_tau).tolist(), "fhat": np.asarray(fhat).tolist(),
            "bandwidth": float(1.06 * np.std(y) * n ** (-0.2)), "n": n, "mean": float(np.mean(y)),
            "sd": float(np.std(y, ddof=1)), "seed": seed, "lower_share": lower_share,
            "direct_check_max_rel_err": err, "stage1_seconds": time.time() - t0,
            "package_file": package_file}))
        (rd / "READY").write_text("ok\n")
        print(f"{a.scenario} rep {rep} n={n} {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
