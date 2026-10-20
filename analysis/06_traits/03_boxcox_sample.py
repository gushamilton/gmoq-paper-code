#!/usr/bin/env python
"""Box-Cox summaries of one raw trait in exactly its Stage 1 analysis sample.

What: rebuilds the Stage 1 analysis sample of one trait with the Stage 1 package's own loaders
  (keep list, non-missing phenotype, complete covariates, presence in the LOCO prediction file
  when LOCO was used, presence in the genotype sample file), so the participants are exactly
  those of the quantile GWAS. A log-analysed trait (phenotype column log_*) is moved back to
  the raw trait, y = exp(log y). It writes one JSON of aggregates only (no participant values
  or identifiers):
    n, number of non-positive values, moment skewness of y and of log y;
    193 quantiles of y at p = 0.02, ..., 0.98 (R type 7), for the normalising power;
    the Box-Cox profile log-likelihood on lambda = -5, -4.99, ..., 5,
        ll(lambda) = -n/2 log(var(bc(y, lambda))) + (lambda - 1) sum log y,
    which is undefined (and not written) when any y <= 0.
Paper: Supplementary Table S5 (normalising and maximum-likelihood Box-Cox powers)
Input: phenotype, covariate, keep and LOCO files of the trait's Stage 1 run; the UK Biobank
  BGEN sample file
Output: one JSON per trait in the private work directory (read by 03_normalising_powers.R)
Runs: on a cluster with UK Biobank data, with the pinned Stage 1 package (gmoq)
"""
import argparse
import json
import sys

import numpy as np

from gmoq import cli


def bc(y, lam):
    return np.log(y) if abs(lam) < 1e-12 else (np.power(y, lam) - 1.0) / lam


def skew(v):
    c = v - v.mean()
    return float(np.mean(c ** 3) / np.mean(c ** 2) ** 1.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--pheno-col", required=True)
    ap.add_argument("--covar", required=True)
    ap.add_argument("--covar-cols", required=True)
    ap.add_argument("--keep", required=True)
    ap.add_argument("--loco", default="none")
    ap.add_argument("--bgen", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # 1. The Stage 1 analysis sample: keep list, valid phenotype, complete covariates.
    ids, aliases, ys, X, cols = cli.load_analysis_inputs(
        a.pheno, [a.pheno_col], a.covar, keep_file=a.keep, covar_cols=a.covar_cols.split())
    y = np.asarray(ys[a.pheno_col], dtype=float)
    # Participants without a LOCO prediction were dropped from a LOCO run.
    if a.loco != "none":
        with open(a.loco) as f:
            header = f.readline().strip().split()
        lookup = {s: i for i, s in enumerate(header[1:])}
        ok = np.array([cli.resolve_sample_column(lookup, al) is not None for al in aliases], dtype=bool)
        ids = [x for x, k in zip(ids, ok) if k]
        y = y[ok]
    # Participants absent from the genotype sample file were dropped by Stage 1.
    from bgen import BgenReader
    bf = BgenReader(a.bgen, sample_path=a.sample, delay_parsing=True)
    samples = list(bf.samples)
    if hasattr(bf, "close"):
        bf.close()
    valid, _ = cli.match_samples(samples, ids)
    y = y[valid]

    # 2. Raw trait and its aggregates.
    analysed_log = a.pheno_col.startswith("log_")
    yr = np.exp(y) if analysed_log else y
    n = int(yr.size)
    pp = np.linspace(0.02, 0.98, 193)
    qy = np.quantile(yr, pp, method="linear")          # R type 7
    grid = np.round(np.arange(-500, 501) / 100.0, 2)
    out = {"trait": a.trait, "pheno_col": a.pheno_col, "analysed_scale": "log" if analysed_log else "raw",
           "n": n, "n_nonpositive": int(np.sum(yr <= 0)), "skew_raw": skew(yr),
           "p_grid": pp.tolist(), "q_raw_type7": qy.tolist(), "lambda_grid": grid.tolist()}
    if np.all(yr > 0):
        out["skew_log"] = skew(np.log(yr))
        sum_log_y = float(np.sum(np.log(yr)))
        out["boxcox_loglik"] = [float(-n / 2.0 * np.log(np.var(bc(yr, lam))) + (lam - 1.0) * sum_log_y)
                                for lam in grid]
    with open(a.out, "w") as fh:
        json.dump(out, fh)
    print(f"{a.trait}: n={n} nonpositive={out['n_nonpositive']}")


if __name__ == "__main__":
    sys.exit(main())
