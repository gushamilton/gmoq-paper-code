#!/usr/bin/env python
"""Quantiles, kernel density and mean of one trait on exactly the Stage 1 sample.

What: rebuilds a trait's Stage 1 analysis sample with the Stage 1 package's own loaders,
in the package's order (keep list, valid phenotype and complete covariates; presence in
the LOCO predictions for LOCO runs; match to the genotype sample file), then evaluates the
package's quantile rule (numpy median_unbiased, R type 8) and kernel density (Silverman
Gaussian, h = 1.06 sd N^-0.2) at the 45 quantile levels. These are the values the RIF of
the Stage 1 GWAS used. Called once per trait by baseline_quantiles.R.
Paper: Figs. 1-5 (baseline quantiles and densities used by Stage 2)
Input: phenotype, covariate, keep and LOCO files of the Stage 1 run; a BGEN file and the
  sample file (for the genotype sample match only)
Output: one JSON of aggregates (counts, 45 quantiles, 45 densities, mean, SD); no
  participant-level values or identifiers
Runs: on a cluster with UK Biobank data, in about a minute per trait
"""
import argparse
import json
import sys

import numpy as np

from gmoq import cli
from gmoq.core import compute_rif


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--pheno-col", required=True)
    ap.add_argument("--covar", required=True)
    ap.add_argument("--covar-cols", required=True, help="space-separated, as in the Stage 1 covariance metadata")
    ap.add_argument("--keep", required=True)
    ap.add_argument("--loco", default="none")
    ap.add_argument("--bgen", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--density-method", default="silverman_gaussian")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    taus = np.linspace(0.02, 0.98, 45)

    # 1. Keep list, valid phenotype and complete covariates (the package loader)
    keep_ids = cli.read_keep_ids(a.keep)
    ids, aliases, ys, X, cols = cli.load_analysis_inputs(
        a.pheno, [a.pheno_col], a.covar, keep_file=a.keep, covar_cols=a.covar_cols.split())
    y = np.asarray(ys[a.pheno_col], dtype=float)
    n_keep, n_valid = len(keep_ids), len(ids)

    # 2. LOCO runs: drop participants absent from the LOCO prediction header
    n_absent_loco = 0
    if a.loco != "none":
        with open(a.loco) as f:
            header = f.readline().strip().split()
        lookup = {s: i for i, s in enumerate(header[1:])}
        ok = np.array([cli.resolve_sample_column(lookup, al) is not None for al in aliases], dtype=bool)
        n_absent_loco = int((~ok).sum())
        ids = [x for x, k in zip(ids, ok) if k]
        y = y[ok]
    n_after_loco = len(ids)

    # 3. Match to the genotype sample file
    from bgen import BgenReader
    bf = BgenReader(a.bgen, sample_path=a.sample, delay_parsing=True)
    samples = list(bf.samples)
    if hasattr(bf, "close"):
        bf.close()
    valid, _ = cli.match_samples(samples, ids)
    y = y[valid]
    n = int(y.size)

    # 4. Quantiles and density exactly as the RIF used them; type 7 quantiles for comparison
    _, q, fhat = compute_rif(y, taus, density_method=a.density_method)
    out = {
        "trait": a.trait, "pheno_col": a.pheno_col, "covariates": cols,
        "n_keep": n_keep, "n_valid_pheno": n_valid, "n_absent_loco": n_absent_loco,
        "n_after_loco": n_after_loco, "n": n,
        "mean": float(y.mean()), "sd": float(y.std(ddof=1)), "sd_pop": float(y.std()),
        "bandwidth": float(1.06 * np.std(y) * n ** (-0.2)),
        "n_distinct": int(np.unique(y).size),
        "taus": taus.tolist(), "q_tau": q.tolist(), "fhat": fhat.tolist(),
        "q_tau_type7": np.quantile(y, taus, method="linear").tolist(),
        "quantile_method": "numpy median_unbiased (R type 8), as gmoq.core.compute_rif",
        "density_method": a.density_method,
    }
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"{a.trait}: n={n} (keep {n_keep} -> valid {n_valid} -> LOCO {n_after_loco} -> genotyped {n})")


if __name__ == "__main__":
    sys.exit(main())
