#!/usr/bin/env python
# Quantile profile of the BMI polygenic score itself
#
# What: runs the Stage 1 quantile regression (RIF regression at the 45 quantile levels) with the
#   leave-one-chromosome-out polygenic score in place of a variant, for each left-out
#   chromosome, on log BMI and BMI, in the full sample and in each random half. This is the
#   score's own quantile profile b_S(u): the reference against which the lead variants' common
#   tilt is measured. It uses the G-MoQ Stage 1 package itself (RIF, covariate
#   residualisation and pooled covariance), with covariates age, sex and ten genetic PCs and no
#   leave-one-chromosome-out prediction:
#     b_S(u)  = coefficient of the score in RIF_u ~ [1, covariates] + score   (Frisch-Waugh)
#     se_S(u) = sqrt(diag(Sigma0_S) / d_S), Sigma0_S from the RIF residualised on both,
#     d_S     = residual sum of squares of the score after the covariates.
# Paper: Fig. 6d; Methods (BMI x PRS interaction)
# Input: the analysis sample without identifiers written by 07_prs_scores.R (private area)
# Output: private aggregate rif/prs_profiles.csv, read by 09_prs_interaction.R
# Runs: on a cluster with UK Biobank data, in about 20 minutes:
#   python 08_prs_profile.py <private_root>/bmi_prs

import os
import sys

import numpy as np
import pandas as pd

from gmoq.cli import build_rif_residuals
from gmoq.pooled_cov import compute_sigma0

BASE = sys.argv[1]
OUT = os.path.join(BASE, "rif")
os.makedirs(OUT, exist_ok=True)
TAUS = np.linspace(0.02, 0.98, 45)
COV = ["age", "sex"] + ["PC%d" % i for i in range(1, 11)]
PRS = ["PRS%02d" % c for c in range(1, 23)]

data = pd.read_csv(os.path.join(BASE, "prep", "rif_input.csv.gz"))
profiles = []
for trait in ["log_bmi", "bmi"]:
    for sample in ["full", "h1", "h2"]:
        rows = np.ones(len(data), bool) if sample == "full" else (data["half"].values == int(sample[1]))
        y = data.loc[rows, trait].values.astype(np.float64)
        X = np.column_stack([np.ones(rows.sum()), data.loc[rows, COV].values.astype(np.float64)])
        # 1. RIF at each quantile, residualised on the covariates (Q: orthonormal basis of X)
        rif_list, Q, q_list = build_rif_residuals({trait: y}, [trait], X, TAUS)
        R = rif_list[0]
        n = R.shape[0]
        # 2. Each leave-one-chromosome-out score in turn
        for c in range(1, 23):
            S = data.loc[rows, PRS[c - 1]].values.astype(np.float64)
            S_resid = S - Q @ (Q.T @ S)
            d_S = float(S_resid @ S_resid)
            b = (R.T @ S_resid) / d_S
            se = np.sqrt(np.diag(compute_sigma0(R - np.outer(S_resid, b))) / d_S)
            for t in range(len(TAUS)):
                profiles.append(dict(scale=trait, sample=sample, left_out_chr=c, n=n, tau=TAUS[t],
                                     b_prs=b[t], se_prs=se[t], d_prs=d_S))
pd.DataFrame(profiles).to_csv(os.path.join(OUT, "prs_profiles.csv"), index=False)
