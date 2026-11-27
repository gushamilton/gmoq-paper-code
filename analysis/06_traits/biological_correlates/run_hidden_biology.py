#!/usr/bin/env python3
"""
ALT and LDL biological correlates:
  1. ALT (316 leads): link scale-invariant PC2 / heterogeneity / skew to UKB 81499
     recorded liver disease (ICD-10 K760/K758, First Occurrence K76/K74) and ALT x BMI interaction.
  2. LDL (230 leads): compare 45-quantile RIF-OLS profiles, whitened PC spectrum,
     parallel-analysis dimensionality, and heterogeneity tau_h between all LDL
     participants and participants reporting no lipid medication, plus statin-use GWAS.

The original numerical operations are preserved. The subgroup parallel-analysis
null and quantile convention differ from the primary analysis (README.md).
All participant-level intermediates stay in --root; only aggregate SNP-level and
trait-level CSVs are written to --root/aggregate/.
"""
import argparse
import csv
import glob
import math
import os
import socket
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

def configure_paths():
    """Restricted input locations are supplied by the authorised researcher."""
    global UKB_DIR, DATA_MAIN, DATA_FO, LINKER, COVAR_FILE, BGEN_SAMPLE
    UKB_DIR = Path(os.environ["GMOQ_UKB_DIR"])
    DATA_MAIN = UKB_DIR / "data.51913.csv"
    DATA_FO = UKB_DIR / "data.669894.csv"
    LINKER = UKB_DIR / "linker.81499.csv"
    COVAR_FILE = Path(os.environ["GMOQ_COVARIATES"])
    BGEN_SAMPLE = Path(os.environ["GMOQ_BGEN_SAMPLE"])
    global KEEP_ALT, KEEP_LDL
    KEEP_ALT = Path(os.environ["GMOQ_KEEP_ALT"])
    KEEP_LDL = Path(os.environ["GMOQ_KEEP_LDL"])

def check_compute_node():
    h = socket.gethostname()
    if not h.startswith("bp1-compute"):
        raise RuntimeError(f"Must run on a BP compute node, got {h}")

def read_keep_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                ids.add(parts[1])
    return ids

def run_prep(root: Path):
    check_compute_node()
    for d in ["keep", "pheno", "snps", "aggregate"]:
        (root / d).mkdir(parents=True, exist_ok=True)

    # 1. Read BGEN sample IDs, covariates, and European QC keep lists
    bgen_ids = set()
    with open(BGEN_SAMPLE) as f:
        next(f); next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                bgen_ids.add(parts[1])

    cov = pd.read_csv(COVAR_FILE, sep=r"\s+", dtype={"FID": str, "IID": str})
    cov = cov.dropna()
    cov_ids = set(cov["IID"])
    keep_alt_ids = read_keep_ids(KEEP_ALT) & bgen_ids & cov_ids
    keep_ldl_ids = read_keep_ids(KEEP_LDL) & bgen_ids & cov_ids
    keep_union = keep_alt_ids | keep_ldl_ids
    print(f"[prep] bgen={len(bgen_ids)}, cov_complete={len(cov_ids)}, keep_alt={len(keep_alt_ids)}, keep_ldl={len(keep_ldl_ids)}, union={len(keep_union)}", flush=True)

    # 2. Load 81499 linker (app eid -> IEU IID), keeping only union IDs
    app_to_ieu = {}
    with open(LINKER, newline="") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) >= 2 and row[0] in keep_union:
                app_to_ieu[row[1]] = row[0]
    print(f"[prep] linked union rows: {len(app_to_ieu)}", flush=True)

    # 3. Stream First Occurrences (data.669894.csv: 3.4 GB) for K70/K74/K75/K76
    fo_fields = {"131658-0.0": "fo_k70", "131666-0.0": "fo_k74", "131668-0.0": "fo_k75", "131670-0.0": "fo_k76"}
    fo_data = {}
    csv.field_size_limit(10_000_000)
    with open(DATA_FO, newline="") as f:
        r = csv.reader(f)
        hdr = next(r)
        eid_j = hdr.index("eid")
        col_j = {name: hdr.index(code) for code, name in fo_fields.items()}
        for row in r:
            ieu = app_to_ieu.get(row[eid_j])
            if ieu is None:
                continue
            fo_data[ieu] = {name: int(row[j] != "") for name, j in col_j.items()}
    print(f"[prep] scanned First Occurrences: matched {len(fo_data)} rows", flush=True)

    # 4. Stream main phenotype file (data.51913.csv) once
    med_codes = [f"6177-0.{i}" for i in range(3)] + [f"6153-0.{i}" for i in range(4)]
    records = {}
    with open(DATA_MAIN, newline="") as f:
        r = csv.reader(f)
        hdr = next(r)
        eid_j = hdr.index("eid")
        alt_j = hdr.index("30620-0.0")
        ldl_j = hdr.index("30780-0.0")
        bmi_j = hdr.index("21001-0.0")
        med_js = [j for j, h in enumerate(hdr) if h in med_codes]
        icd_js = [j for j, h in enumerate(hdr) if h.startswith("41270-")]
        print(f"[prep] main CSV columns: med={len(med_js)}, icd10={len(icd_js)}", flush=True)
        for n_i, row in enumerate(r, 1):
            ieu = app_to_ieu.get(row[eid_j])
            if ieu is None:
                continue
            alt_v = row[alt_j]
            ldl_v = row[ldl_j]
            bmi_v = row[bmi_j]
            statin = 1 if any(row[j] == "1" for j in med_js) else 0
            k760 = 0; k758 = 0; k74 = 0; k70 = 0
            for j in icd_js:
                code = row[j]
                if not code:
                    continue
                if code == "K760":
                    k760 = 1
                elif code == "K758":
                    k758 = 1
                elif code.startswith("K74"):
                    k74 = 1
                elif code.startswith("K70"):
                    k70 = 1
            fo = fo_data.get(ieu, {"fo_k70": 0, "fo_k74": 0, "fo_k75": 0, "fo_k76": 0})
            masld_k760 = 1 if (k760 or k758) else 0
            cirrhosis = 1 if (k74 or fo["fo_k74"]) else 0
            fo_k76 = fo["fo_k76"]
            alc = 1 if (k70 or fo["fo_k70"]) else 0
            masld_any = 1 if (masld_k760 or fo_k76 or cirrhosis) else 0
            masld_nonalc = 1 if (masld_any and not alc) else 0
            records[ieu] = (
                float(alt_v) if alt_v != "" else np.nan,
                float(ldl_v) if ldl_v != "" else np.nan,
                float(bmi_v) if bmi_v != "" else np.nan,
                statin, masld_k760, fo_k76, cirrhosis, masld_any, masld_nonalc, alc
            )
            if n_i % 100000 == 0:
                print(f"  [prep] scanned {n_i} rows of data.51913.csv", flush=True)

    df = pd.DataFrame.from_dict(
        records, orient="index",
        columns=["alt", "ldl", "bmi", "statin", "masld_k760", "fo_k76", "cirrhosis", "masld_any", "masld_nonalc", "alc"]
    ).reset_index().rename(columns={"index": "IID"})
    df["FID"] = df["IID"]

    # Build ALT cohort (exact non-missing alt > 0 in keep_alt_ids)
    df_alt = df[df["IID"].isin(keep_alt_ids) & df["alt"].notna() & (df["alt"] > 0)].copy()
    df_alt["log_alt"] = np.log(df_alt["alt"])
    # Build LDL cohort (exact non-missing ldl in keep_ldl_ids)
    df_ldl = df[df["IID"].isin(keep_ldl_ids) & df["ldl"].notna()].copy()

    # Save plink2 keep and pheno files (binary coded 1=ctrl, 2=case for plink2)
    df_alt[["FID", "IID"]].to_csv(root / "keep/keep_alt.txt", sep=" ", index=False, header=False)
    df_ldl[["FID", "IID"]].to_csv(root / "keep/keep_ldl.txt", sep=" ", index=False, header=False)
    pd.concat([df_alt[["FID", "IID"]], df_ldl[["FID", "IID"]]]).drop_duplicates().to_csv(
        root / "keep/keep_union.txt", sep=" ", index=False, header=False
    )

    p_alt = pd.DataFrame({
        "FID": df_alt["FID"], "IID": df_alt["IID"],
        "log_alt": df_alt["log_alt"], "bmi": df_alt["bmi"],
        "masld_k760": df_alt["masld_k760"] + 1,
        "fo_k76": df_alt["fo_k76"] + 1,
        "cirrhosis": df_alt["cirrhosis"] + 1,
        "masld_any": df_alt["masld_any"] + 1,
        "masld_nonalc": df_alt["masld_nonalc"] + 1,
    })
    p_alt.to_csv(root / "pheno/pheno_alt.txt", sep="\t", index=False, na_rep="NA")

    p_ldl = pd.DataFrame({
        "FID": df_ldl["FID"], "IID": df_ldl["IID"],
        "ldl": df_ldl["ldl"],
        "statin": df_ldl["statin"] + 1,
    })
    p_ldl.to_csv(root / "pheno/pheno_ldl.txt", sep="\t", index=False, na_rep="NA")

    # Save aggregate sample counts
    counts = [
        ("alt", "n_total", len(df_alt)),
        ("alt", "n_bmi_nonmissing", int(df_alt["bmi"].notna().sum())),
        ("alt", "cases_masld_k760_k758", int(df_alt["masld_k760"].sum())),
        ("alt", "cases_fo_k76", int(df_alt["fo_k76"].sum())),
        ("alt", "cases_cirrhosis_k74", int(df_alt["cirrhosis"].sum())),
        ("alt", "cases_masld_any", int(df_alt["masld_any"].sum())),
        ("alt", "cases_masld_nonalc", int(df_alt["masld_nonalc"].sum())),
        ("alt", "cases_alc_k70", int(df_alt["alc"].sum())),
        ("ldl", "n_total", len(df_ldl)),
        ("ldl", "n_statin_users", int(df_ldl["statin"].sum())),
        ("ldl", "n_untreated", int((df_ldl["statin"] == 0).sum())),
        ("ldl", "mean_ldl_all", float(df_ldl["ldl"].mean())),
        ("ldl", "mean_ldl_untreated", float(df_ldl.loc[df_ldl["statin"] == 0, "ldl"].mean())),
        ("ldl", "mean_ldl_treated", float(df_ldl.loc[df_ldl["statin"] == 1, "ldl"].mean())),
    ]
    pd.DataFrame(counts, columns=["cohort", "metric", "value"]).to_csv(root / "aggregate/sample_counts.csv", index=False)
    print(pd.DataFrame(counts, columns=["cohort", "metric", "value"]).to_string(index=False), flush=True)

    # Write per-chromosome SNP lists for bgenix
    alt_leads = pd.read_csv(root / "inputs/alt_leads.tsv", sep="\t")
    ldl_leads = pd.read_csv(root / "inputs/ldl_leads.tsv", sep="\t")
    for tag, ldf in [("alt", alt_leads), ("ldl", ldl_leads)]:
        for c, grp in ldf.groupby("chr"):
            (root / "snps" / f"{tag}_chr{int(c):02d}.txt").write_text("\n".join(grp["snp_id"].astype(str)) + "\n")
    print("[prep] done!", flush=True)

# ---- RIF-OLS, moments, whitened PCA and historical subgroup parallel analysis ----
def compute_rif_matrix(y, taus):
    n = len(y)
    qs = np.quantile(y, taus)
    sd = np.std(y, ddof=1)
    bw = 1.06 * sd * (n ** (-0.2))
    rif = np.empty((n, len(taus)), dtype=np.float64)
    fhat = np.empty(len(taus), dtype=np.float64)
    inv_sqrt2pi = 1.0 / math.sqrt(2.0 * math.pi)
    for j, (t, qt) in enumerate(zip(taus, qs)):
        z = (y - qt) / bw
        fj = np.mean(np.exp(-0.5 * z * z)) * inv_sqrt2pi / bw
        fhat[j] = fj
        rif[:, j] = qt + (t - (y <= qt).astype(np.float64)) / fj
    return rif, qs, fhat

def residualise_matrix(Y, X):
    # X includes intercept; QR projection
    Q, _ = np.linalg.qr(X, mode="reduced")
    return Y - Q @ (Q.T @ Y)

def load_dosages(geno_dir: Path, prefix: str, target_iids: np.ndarray, leads_df: pd.DataFrame):
    """Load plink2 --export A .raw files and align to target_iids and leads_df effect_allele."""
    iid_to_pos = {iid: idx for idx, iid in enumerate(target_iids)}
    n = len(target_iids)
    G = np.full((n, len(leads_df)), np.nan, dtype=np.float64)
    snp_to_col = {s: i for i, s in enumerate(leads_df["snp_id"].astype(str))}
    ea_map = dict(zip(leads_df["snp_id"].astype(str), leads_df["effect_allele"].astype(str)))
    oa_map = dict(zip(leads_df["snp_id"].astype(str), leads_df["other_allele"].astype(str)))

    found = 0
    for raw_path in sorted(geno_dir.glob(f"{prefix}_chr*.raw")):
        df_r = pd.read_csv(raw_path, sep=r"\s+", dtype={"FID": str, "IID": str})
        idx_keep = [i for i, iid in enumerate(df_r["IID"]) if iid in iid_to_pos]
        pos_dst = np.array([iid_to_pos[df_r["IID"].iat[i]] for i in idx_keep], dtype=int)
        for col in df_r.columns[6:]:
            # Format: rsID_COUNTEDALLELE (or rsID_ALLELE(/OTHER))
            parts = col.rsplit("_", 1)
            snp = parts[0]
            counted = parts[1].split("(")[0]
            if snp not in snp_to_col:
                continue
            vals = df_r[col].to_numpy(dtype=np.float64)[idx_keep]
            # Mean impute rare missing dosages if any
            if np.isnan(vals).any():
                vals = np.where(np.isnan(vals), np.nanmean(vals), vals)
            ea = ea_map[snp]; oa = oa_map[snp]
            if counted == ea:
                aligned = vals
            elif counted == oa:
                aligned = 2.0 - vals
            else:
                raise ValueError(f"Allele mismatch for {snp}: counted={counted}, ea={ea}, oa={oa}")
            G[pos_dst, snp_to_col[snp]] = aligned
            found += 1
    print(f"[dosages:{prefix}] loaded {found}/{len(leads_df)} leads across {n} samples", flush=True)
    return G

def rif_ols_leads(y, X_cov, G, taus):
    n, k = G.shape
    rif, qs, fhat = compute_rif_matrix(y, taus)
    # Residualise RIF and G on X_cov (with intercept)
    Q, _ = np.linalg.qr(X_cov, mode="reduced")
    Psi = rif - Q @ (Q.T @ rif)
    Sigma0 = (Psi.T @ Psi) / n
    G_res = G - Q @ (Q.T @ G)
    d_resid = np.sum(G_res * G_res, axis=0)
    B = (G_res.T @ Psi) / d_resid[:, None]  # k x 45
    return B, Sigma0, d_resid, qs, fhat

def trap_weights(u):
    w = np.empty_like(u)
    w[0] = 0.5 * (u[1] - u[0])
    w[-1] = 0.5 * (u[-1] - u[-2])
    w[1:-1] = 0.5 * (u[2:] - u[:-2])
    return w

def compute_moments_and_pca(B, Sigma0, d_resid, u, q, seed=20260923, n_pa=300):
    k, nt = B.shape
    cw = trap_weights(u)
    w_mu = cw
    w_l2 = cw * (2.0 * u - 1.0)
    w_l3 = cw * (6.0 * u * u - 6.0 * u + 1.0)

    # Shift + spread frame residual for skew_beyond_frame (matching P11)
    q_mu = float(np.sum(cw * q) / np.sum(cw))
    X_fr = np.column_stack([np.ones(nt), q - q_mu])
    P_res = np.eye(nt) - X_fr @ np.linalg.pinv(X_fr)
    w_sk = P_res.T @ w_l3

    def func_est(w):
        est = B @ w
        var0 = float(w @ Sigma0 @ w)
        se = np.sqrt(var0 / d_resid)
        return est, se, est / se

    d_mu, se_mu, z_mu = func_est(w_mu)
    d_l2, se_l2, z_l2 = func_est(w_l2)
    d_sk, se_sk, z_sk = func_est(w_sk)

    # Unit-SE orthonormal shift/spread basis for Cochran Q heterogeneity tau_h (matching P11)
    # Use c1 = w_mu / sqrt(w_mu' Sigma0 w_mu), c2 = Gram-Schmidt of w_l2 orthogonal to c1 under Sigma0
    s11 = float(w_mu @ Sigma0 @ w_mu)
    c1 = w_mu / math.sqrt(s11)
    s12 = float(c1 @ Sigma0 @ w_l2)
    w2_orth = w_l2 - s12 * c1
    s22 = float(w2_orth @ Sigma0 @ w2_orth)
    c2 = w2_orth / math.sqrt(s22)

    L = B @ c1; S = B @ c2
    se_c = 1.0 / np.sqrt(d_resid)  # unit-SE basis: Var(L_i) = Var(S_i) = 1/d_i, Cov(L_i, S_i) = 0!
    # Iteratively reweighted common slope r and heterogeneity tau_h (Q - (k - 1))
    r = float(np.sum(d_resid * L * S) / np.sum(d_resid * L * L))
    for _ in range(15):
        v_i = (1.0 + r * r) / d_resid
        w_i = 1.0 / v_i
        r_new = float(np.sum(w_i * L * S) / np.sum(w_i * L * L))
        if abs(r_new - r) < 1e-10:
            break
        r = r_new
    res = S - r * L
    v_i = (1.0 + r * r) / d_resid
    Q_het = float(np.sum(res * res / v_i))
    p_het = float(stats.chi2.sf(Q_het, df=k - 1))
    # Excess tau^2: E[sum (res_i^2 / v_i)] = (k - 1) + tau^2 * sum(L_i^2 / v_i)
    denom_tau = float(np.sum((L * L) / v_i))
    tau_h = math.sqrt(max(0.0, (Q_het - (k - 1)) / denom_tau)) if denom_tau > 0 else 0.0
    het_resid_z = res / np.sqrt(v_i)

    # Whitened unit-length PCA and historical subgroup parallel analysis
    # This fitted null differs from the primary procedure: rank reduction follows
    # unit normalisation here. See README.md before interpreting component counts.
    eigvals, eigvecs = np.linalg.eigh(Sigma0)
    W_ih = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T
    Y = (B @ W_ih) * np.sqrt(d_resid)[:, None]
    norms = np.linalg.norm(Y, axis=1)
    U = Y / norms[:, None]
    u_svd, s_svd, vt_svd = np.linalg.svd(U, full_matrices=False)
    pct = 100.0 * (s_svd ** 2) / np.sum(s_svd ** 2)

    # Parallel analysis for PC1..PC6
    rng = np.random.default_rng(seed)
    pa_rows = []
    n_comp = 0
    for c_test in range(1, 7):
        if c_test == 1:
            signal = np.zeros_like(Y)
        else:
            # Rank-(c_test - 1) reconstruction in whitened Y space scaled by norms
            U_rec = (u_svd[:, :c_test - 1] * s_svd[:c_test - 1]) @ vt_svd[:c_test - 1, :]
            signal = U_rec * norms[:, None]
        sim_pcts = np.empty(n_pa)
        for rep in range(n_pa):
            Y_sim = signal + rng.standard_normal(Y.shape)
            U_sim = Y_sim / np.linalg.norm(Y_sim, axis=1, keepdims=True)
            s_sim = np.linalg.svd(U_sim, compute_uv=False)
            sim_pcts[rep] = 100.0 * (s_sim[c_test - 1] ** 2) / np.sum(s_sim ** 2)
        n_med = float(np.median(sim_pcts))
        n_p95 = float(np.quantile(sim_pcts, 0.95))
        obs = float(pct[c_test - 1])
        pa_rows.append({
            "pc": c_test, "observed_pct": obs,
            "null_median_pct": n_med, "null_p95_pct": n_p95,
            "ratio_to_null_median": obs / n_med,
            "above_null_p95": bool(obs > n_p95),
        })
        if c_test == n_comp + 1 and obs > n_p95:
            n_comp = c_test

    # PC scores (unit-length coordinates = U @ v_j)
    pc1 = U @ vt_svd[0]
    pc2 = U @ vt_svd[1]
    pc3 = U @ vt_svd[2]
    if np.median(pc1) < 0:
        pc1 = -pc1

    per_lead = pd.DataFrame({
        "d_resid": d_resid,
        "delta_mu": d_mu, "delta_mu_se": se_mu, "delta_mu_z": z_mu,
        "delta_lambda2": d_l2, "delta_lambda2_se": se_l2, "delta_lambda2_z": z_l2,
        "skew_beyond_frame": d_sk, "skew_beyond_frame_se": se_sk, "skew_beyond_frame_z": z_sk,
        "het_resid_z": het_resid_z,
        "pc1_score": pc1, "pc2_score": pc2, "pc3_score": pc3,
    })
    summary = {
        "n_leads": k, "n_components": n_comp, "tau_h": tau_h, "Q_het": Q_het, "p_het": p_het,
        "pc1_pct": float(pct[0]), "pc2_pct": float(pct[1]), "pc3_pct": float(pct[2]),
        "pc2_ratio": pa_rows[1]["ratio_to_null_median"],
        "pc2_null_p95_ratio": pa_rows[1]["null_p95_pct"] / pa_rows[1]["null_median_pct"],
        "n_skew_bonf": int(np.sum(np.abs(z_sk) > stats.norm.isf(0.025 / k))),
        "n_het_bonf": int(np.sum(np.abs(het_resid_z) > stats.norm.isf(0.025 / k))),
    }
    return per_lead, pd.DataFrame(pa_rows), summary

def read_plink_glm(glm_dir: Path, pattern: str, leads_df: pd.DataFrame, col_prefix: str):
    dfs = []
    for p in sorted(glm_dir.glob(pattern)):
        df = pd.read_csv(p, sep="\t")
        if "TEST" in df.columns:
            df = df[df["TEST"] == "ADD"]
        dfs.append(df)
    if not dfs:
        raise FileNotFoundError(f"No GLM files matching {pattern} in {glm_dir}")
    g = pd.concat(dfs, ignore_index=True)
    # Harmonize to leads_df effect_allele
    ea_map = dict(zip(leads_df["snp_id"].astype(str), leads_df["effect_allele"].astype(str)))
    oa_map = dict(zip(leads_df["snp_id"].astype(str), leads_df["other_allele"].astype(str)))
    g = g[g["ID"].isin(ea_map)].copy()
    # In plink2 logistic, OR or BETA; detect column
    if "LOG(OR)_SE" in g.columns:
        beta = np.log(g["OR"].astype(float))
        se = g["LOG(OR)_SE"].astype(float)
    else:
        beta = g["BETA"].astype(float)
        se = g["SE"].astype(float)
    a1 = g["A1"].astype(str)
    signs = np.array([1.0 if a == ea_map[s] else (-1.0 if a == oa_map[s] else np.nan) for s, a in zip(g["ID"], a1)])
    beta = beta * signs
    z = beta / se
    out = pd.DataFrame({
        "snp_id": g["ID"].astype(str),
        f"{col_prefix}_beta": beta,
        f"{col_prefix}_se": se,
        f"{col_prefix}_z": z,
        f"{col_prefix}_p": g["P"].astype(float),
    })
    return out.drop_duplicates("snp_id")

def run_analyze(root: Path):
    check_compute_node()
    taus = np.linspace(0.02, 0.98, 45)
    cov = pd.read_csv(COVAR_FILE, sep=r"\s+", dtype={"FID": str, "IID": str}).dropna()
    cov_cols = ["age", "sex"] + [f"PC{i}" for i in range(1, 11)]

    # ===================== 1. ALT: MASLD + BMI interaction =====================
    alt_leads = pd.read_csv(root / "inputs/alt_leads.tsv", sep="\t")
    alt_df = alt_leads.copy()
    for pheno_name in ["masld_k760", "fo_k76", "cirrhosis", "masld_any", "masld_nonalc"]:
        gdf = read_plink_glm(root / "glm", f"alt_chr*.{pheno_name}.glm.logistic.hybrid", alt_leads, pheno_name)
        alt_df = alt_df.merge(gdf, on="snp_id", how="left")

    # ALT x BMI interaction with HC0 robust SE
    p_alt = pd.read_csv(root / "pheno/pheno_alt.txt", sep="\t", dtype={"FID": str, "IID": str})
    m_alt = p_alt.merge(cov[["IID"] + cov_cols], on="IID", how="inner").dropna(subset=["log_alt", "bmi"])
    G_alt = load_dosages(root / "geno", "alt", m_alt["IID"].to_numpy(), alt_leads)
    y_alt = m_alt["log_alt"].to_numpy(dtype=np.float64)
    bmi_std = (m_alt["bmi"].to_numpy(dtype=np.float64) - m_alt["bmi"].mean()) / m_alt["bmi"].std()
    X_base = np.column_stack([np.ones(len(m_alt)), m_alt[cov_cols].to_numpy(dtype=np.float64), bmi_std])
    Q_b, _ = np.linalg.qr(X_base, mode="reduced")
    y_r = y_alt - Q_b @ (Q_b.T @ y_alt)

    int_b = np.empty(len(alt_leads)); int_se = np.empty(len(alt_leads)); int_z = np.empty(len(alt_leads))
    for j in range(len(alt_leads)):
        g_j = G_alt[:, j]
        gx_j = g_j * bmi_std
        Z = np.column_stack([g_j, gx_j])
        Z_r = Z - Q_b @ (Q_b.T @ Z)
        ztz_inv = np.linalg.inv(Z_r.T @ Z_r)
        b_hat = ztz_inv @ (Z_r.T @ y_r)
        e = y_r - Z_r @ b_hat
        # HC0 sandwich
        meat = (Z_r * e[:, None]).T @ (Z_r * e[:, None])
        V = ztz_inv @ meat @ ztz_inv
        int_b[j] = b_hat[1]
        int_se[j] = math.sqrt(V[1, 1])
        int_z[j] = int_b[j] / int_se[j]
    alt_df["bmi_int_beta"] = int_b
    alt_df["bmi_int_se"] = int_se
    alt_df["bmi_int_z"] = int_z
    alt_df["bmi_int_p"] = 2.0 * stats.norm.sf(np.abs(int_z))
    alt_df.to_csv(root / "aggregate/alt_masld_per_lead.csv", index=False)

    # Summary correlations for ALT
    alt_sum_rows = []
    for xcol in ["pc2_score", "het_resid_z", "skew_beyond_frame_z"]:
        for ycol in ["masld_k760_z", "fo_k76_z", "cirrhosis_z", "masld_any_z", "masld_nonalc_z", "bmi_int_z"]:
            sub = alt_df[[xcol, ycol]].dropna()
            r_p, p_p = stats.pearsonr(sub[xcol], sub[ycol])
            r_s, p_s = stats.spearmanr(sub[xcol], sub[ycol])
            alt_sum_rows.append({"x": xcol, "y": ycol, "n": len(sub), "pearson_r": r_p, "pearson_p": p_p, "spearman_rho": r_s, "spearman_p": p_s})
    pd.DataFrame(alt_sum_rows).to_csv(root / "aggregate/alt_masld_summary.csv", index=False)
    print("\n=== ALT MASLD SUMMARY ===", flush=True)
    print(pd.DataFrame(alt_sum_rows).to_string(index=False), flush=True)

    # ===================== 2. LDL: All vs Untreated (statin-free) =====================
    ldl_leads = pd.read_csv(root / "inputs/ldl_leads.tsv", sep="\t")
    statin_glm = read_plink_glm(root / "glm", "ldl_chr*.statin.glm.logistic.hybrid", ldl_leads, "statin")
    p_ldl = pd.read_csv(root / "pheno/pheno_ldl.txt", sep="\t", dtype={"FID": str, "IID": str})
    m_ldl = p_ldl.merge(cov[["IID"] + cov_cols], on="IID", how="inner").dropna(subset=["ldl", "statin"])
    G_ldl = load_dosages(root / "geno", "ldl", m_ldl["IID"].to_numpy(), ldl_leads)

    arms = {
        "ldl_all":       np.ones(len(m_ldl), dtype=bool),
        "ldl_untreated": (m_ldl["statin"].to_numpy() == 1),  # 1 = control (no statin)
    }
    pca_tables = []
    arm_summaries = []
    profile_rows = []
    per_lead_merged = ldl_leads[["snp_id", "chr", "bp", "effect_allele", "other_allele",
                                 "delta_mu_z", "het_resid_z", "skew_beyond_frame_z", "pc1_score", "pc2_score"]].rename(
        columns={"delta_mu_z": "p11_delta_mu_z", "het_resid_z": "p11_het_resid_z",
                 "skew_beyond_frame_z": "p11_skew_z", "pc1_score": "p11_pc1", "pc2_score": "p11_pc2"}
    ).merge(statin_glm, on="snp_id", how="left")

    for arm_name, mask in arms.items():
        y_a = m_ldl.loc[mask, "ldl"].to_numpy(dtype=np.float64)
        X_a = np.column_stack([np.ones(int(mask.sum())), m_ldl.loc[mask, cov_cols].to_numpy(dtype=np.float64)])
        G_a = G_ldl[mask, :]
        B_a, S0_a, d_a, q_a, f_a = rif_ols_leads(y_a, X_a, G_a, taus)
        pl_a, pa_a, sum_a = compute_moments_and_pca(B_a, S0_a, d_a, taus, q_a)
        sum_a["arm"] = arm_name
        sum_a["n_samples"] = int(mask.sum())
        arm_summaries.append(sum_a)
        pa_a["arm"] = arm_name
        pca_tables.append(pa_a)

        # Merge per-lead stats
        pl_a.insert(0, "snp_id", ldl_leads["snp_id"].astype(str).to_numpy())
        rename_map = {c: f"{arm_name}_{c}" for c in pl_a.columns if c != "snp_id"}
        per_lead_merged = per_lead_merged.merge(pl_a.rename(columns=rename_map), on="snp_id", how="left")

        # Profiles table
        se_diag = np.sqrt(np.diag(S0_a))
        for i, sid in enumerate(ldl_leads["snp_id"].astype(str)):
            se_i = se_diag / math.sqrt(d_a[i])
            for t_idx, tv in enumerate(taus):
                profile_rows.append({
                    "arm": arm_name, "snp_id": sid, "tau": round(float(tv), 4),
                    "beta": float(B_a[i, t_idx]), "se": float(se_i[t_idx]),
                })

    # Verify all-LDL reproduction against P11
    r_repro_mu = stats.pearsonr(per_lead_merged["p11_delta_mu_z"], per_lead_merged["ldl_all_delta_mu_z"])[0]
    r_repro_het = stats.pearsonr(per_lead_merged["p11_het_resid_z"], per_lead_merged["ldl_all_het_resid_z"])[0]
    r_statin_all = stats.pearsonr(per_lead_merged["ldl_all_het_resid_z"], per_lead_merged["statin_z"])
    r_statin_unt = stats.pearsonr(per_lead_merged["ldl_untreated_het_resid_z"], per_lead_merged["statin_z"])
    print(f"\n=== LDL REPRODUCTION & STATIN CHECK ===", flush=True)
    print(f"Reproduction vs P11: r(delta_mu_z)={r_repro_mu:.5f}, r(het_resid_z)={r_repro_het:.5f}", flush=True)
    print(f"cor(het_resid_z, statin_z): ALL={r_statin_all[0]:.4f} (p={r_statin_all[1]:.2e}), UNTREATED={r_statin_unt[0]:.4f} (p={r_statin_unt[1]:.2e})", flush=True)
    print(pd.DataFrame(arm_summaries).to_string(index=False), flush=True)

    per_lead_merged.to_csv(root / "aggregate/ldl_statin_per_lead.csv", index=False)
    pd.concat(pca_tables, ignore_index=True).to_csv(root / "aggregate/ldl_statin_pca_parallel.csv", index=False)
    pd.DataFrame(arm_summaries).to_csv(root / "aggregate/ldl_statin_summary.csv", index=False)
    pd.DataFrame(profile_rows).to_csv(root / "aggregate/ldl_statin_profiles.csv", index=False)
    print("[analyze] all aggregate tables written cleanly!", flush=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["prep", "analyze"])
    ap.add_argument("--root", required=True)
    args = ap.parse_args()
    configure_paths()
    root = Path(args.root)
    if args.mode == "prep":
        run_prep(root)
    else:
        run_analyze(root)

if __name__ == "__main__":
    main()
