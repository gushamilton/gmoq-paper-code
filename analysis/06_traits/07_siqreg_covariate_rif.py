#!/usr/bin/env python
"""Quantile-effect profiles of non-genetic covariates and of the BMI polygenic score on raw BMI.

What: computes, with the Stage 1 package's own RIF code (no copy of it here), the unconditional
  quantile-effect profile b_S(tau) of each predictor S on RAW BMI, in the Stage 1 BMI sample:
    S in {age, age^2, sex, PC1-PC40}   the 43 non-genetic covariates SIQReg uses for BMI,
    S in {PRS01, ..., PRS22}           the BMI polygenic score leaving out chromosome c.
  Two designs per predictor:
    marginal  RIF_tau ~ 1 + S                           (SIQReg fits one covariate at a time)
    adjusted  RIF_tau ~ 1 + S + [age, sex, PC1-PC10] without S   (the Stage 1 covariates)
  By Frisch-Waugh, b_S = (S_r' R) / d_S with S_r the predictor residualised on the other
  columns and d_S = ||S_r||^2; its pooled covariance is Sigma0_S / d_S, with Sigma0_S from the
  package's compute_sigma0 on the RIF residuals. Because the RIF is on raw BMI, b_S is on the
  raw measurement scale, like the raw-scale lead profiles whose flat power 02_flat_powers.R finds.
Paper: Discussion (flat Box-Cox power of BMI from covariates and from the polygenic score);
  Methods "Change of scale"
Input: <work>/prep/rif_input.csv.gz (BMI, age, sex, PC1-PC40, PRS01-PRS22; no identifiers),
  written by 07_siqreg_flat_power.R
Output: <work>/rif/{profiles,sigma0,nuisance}.csv (aggregates; kept in the private work
  directory until 07_siqreg_flat_power.R exports the derived powers)
Runs: on a cluster with UK Biobank data, with the pinned Stage 1 package (gmoq)
Usage: 07_siqreg_covariate_rif.py --run <work> [--prep <file>] [--tau-from 0.02 --tau-to 0.98 --n-tau 45]
"""
import os
import sys

import numpy as np
import pandas as pd

from gmoq.cli import build_rif_residuals
from gmoq.core import compute_rif
from gmoq.pooled_cov import compute_sigma0


def main():
    argv = dict(zip(sys.argv[1:-1:2], sys.argv[2::2]))
    base = argv["--run"]
    n_tau = int(argv.get("--n-tau", 45))
    taus = np.linspace(float(argv.get("--tau-from", 0.02)), float(argv.get("--tau-to", 0.98)), n_tau)
    src = argv.get("--prep", os.path.join(base, "prep", "rif_input.csv.gz"))
    out = os.path.join(base, "rif")
    os.makedirs(out, exist_ok=True)

    D = pd.read_csv(src)
    n = len(D)
    STAGE1 = ["age", "sex"] + ["PC%d" % i for i in range(1, 11)]       # the Stage 1 covariates
    COVS = ["age", "age2", "sex"] + ["PC%d" % i for i in range(1, 41)]  # SIQReg's 43, in their order
    PRS = ["PRS%02d" % c for c in range(1, 23)]
    PRED = COVS + PRS
    D["age2"] = D["age"].values.astype(np.float64) ** 2
    if D[PRED + ["bmi"]].isna().any().any():
        sys.exit("missing values in the analysis sample")
    y = D["bmi"].values.astype(np.float64)
    if np.any(y <= 0):
        sys.exit("non-positive BMI: Box-Cox powers need a positive trait")

    # 1. RIF of raw BMI at every quantile level (Stage 1 default kernel density).
    rif, q_tau, fhat = compute_rif(y, taus)

    # 2. Adjustment sets: none (marginal), or the Stage 1 covariates without the predictor itself.
    #    age stays in the set for age2, whose adjusted profile is then the quadratic component.
    groups = {}
    for p in PRED:
        groups.setdefault(tuple(c for c in STAGE1 if c != p), []).append(p)

    prof, sig = [], []

    def run_block(mode, cov_names, preds):
        X = np.column_stack([np.ones(n)] + [D[c].values.astype(np.float64) for c in cov_names]) \
            if cov_names else np.ones((n, 1))
        # RIF residualised on X by the Stage 1 code path; Q is an orthonormal basis of X.
        rif_list, Q, q_list = build_rif_residuals({"bmi": y}, ["bmi"], X, taus)
        R = rif_list[0]
        assert np.max(np.abs(np.asarray(q_list[0]) - q_tau)) <= 1e-10
        for p in preds:
            S = D[p].values.astype(np.float64)
            Sr = S - Q @ (Q.T @ S)                     # predictor residualised on X
            dS = float(Sr @ Sr)
            b = (R.T @ Sr) / dS                        # one coefficient per quantile level
            S0 = compute_sigma0(R - np.outer(Sr, b))   # pooled covariance of the residual RIFs
            se = np.sqrt(np.diag(S0) / dS)
            for t in range(n_tau):
                prof.append(dict(predictor=p, mode=mode, n=n, tau=taus[t], b=b[t], se=se[t],
                                 d=dS, sd_predictor=float(np.std(S, ddof=1)),
                                 sd_resid_predictor=float(np.sqrt(dS / n))))
                sig.append(dict(predictor=p, mode=mode, tau=taus[t],
                                **{"s_%02d" % j: S0[t, j] for j in range(n_tau)}))

    run_block("marginal", (), PRED)
    for cov_names, preds in groups.items():
        run_block("adjusted", cov_names, preds)

    pd.DataFrame(prof).to_csv(os.path.join(out, "profiles.csv"), index=False)
    pd.DataFrame(sig).to_csv(os.path.join(out, "sigma0.csv"), index=False)
    pd.DataFrame([dict(tau=taus[t], q_raw=q_tau[t], fhat=fhat[t], n=n) for t in range(n_tau)]) \
        .to_csv(os.path.join(out, "nuisance.csv"), index=False)


if __name__ == "__main__":
    main()
