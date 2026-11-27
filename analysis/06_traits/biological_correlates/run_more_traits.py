#!/usr/bin/env python3
"""
Biological correlates for four traits:
  Modifier, disease and dosage comparisons for four UK Biobank traits:
    1. FEV1/FVC (421 leads): all vs never-smokers vs ever-smokers; G x smoking; COPD (J44).
    2. log urate (396 leads): all vs female vs male; G x sex; gout (M10).
    3. log bilirubin (183 leads): dominance deviation (recessive UGT1A1); all vs excluding UGT1A1 homozygotes.
    4. log CRP (389 leads): all vs excluding acute CRP > 10 mg/L; G x BMI.

The original numerical operations are preserved. The subgroup parallel-analysis
null and quantile convention differ from the primary analysis (README.md).
All participant-level intermediates stay in --root; only aggregate CSVs go to --root/aggregate/.
"""
import argparse
import csv
import math
import os
import socket
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
    global KEEP_DIR
    KEEP_DIR = Path(os.environ["GMOQ_KEEP_DIR"])

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

    bgen_ids = set()
    with open(BGEN_SAMPLE) as f:
        next(f); next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                bgen_ids.add(parts[1])

    cov = pd.read_csv(COVAR_FILE, sep=r"\s+", dtype={"FID": str, "IID": str}).dropna()
    base_ids = bgen_ids & set(cov["IID"])

    keeps = {
        "fev":   read_keep_ids(KEEP_DIR / "keep_fev1_fvc_ratio.txt") & base_ids,
        "urate": read_keep_ids(KEEP_DIR / "keep_urate.txt") & base_ids,
        "bil":   read_keep_ids(KEEP_DIR / "keep_bilirubin_total.txt") & base_ids,
        "crp":   read_keep_ids(KEEP_DIR / "keep_crp.txt") & base_ids,
    }
    union_ids = set().union(*keeps.values())
    print(f"[prep] union={len(union_ids)}, " + ", ".join(f"{k}={len(v)}" for k, v in keeps.items()), flush=True)

    app_to_ieu = {}
    with open(LINKER, newline="") as f:
        r = csv.reader(f); next(r)
        for row in r:
            if len(row) >= 2 and row[0] in union_ids:
                app_to_ieu[row[1]] = row[0]
    print(f"[prep] linked union={len(app_to_ieu)}", flush=True)

    csv.field_size_limit(10_000_000)
    records = {}
    with open(DATA_MAIN, newline="") as f:
        r = csv.reader(f); hdr = next(r)
        eid_j = hdr.index("eid")
        j_fev1 = hdr.index("20150-0.0")
        j_fvc  = hdr.index("20151-0.0")
        j_smok = hdr.index("20116-0.0")
        j_py   = hdr.index("20161-0.0")
        j_bmi  = hdr.index("21001-0.0")
        j_ura  = hdr.index("30880-0.0")
        j_bil  = hdr.index("30840-0.0")
        j_crp  = hdr.index("30710-0.0")
        icd_js = [j for j, h in enumerate(hdr) if h.startswith("41270-")]

        def pf(s):
            return float(s) if s != "" else np.nan

        for n_i, row in enumerate(r, 1):
            ieu = app_to_ieu.get(row[eid_j])
            if ieu is None:
                continue
            copd = 0; gout = 0
            for j in icd_js:
                c = row[j]
                if not c:
                    continue
                if c.startswith("J44"):
                    copd = 1
                elif c.startswith("M10"):
                    gout = 1
            fev1 = pf(row[j_fev1]); fvc = pf(row[j_fvc])
            ratio = fev1 / fvc if (not np.isnan(fev1) and not np.isnan(fvc) and fvc > 0) else np.nan
            smok = pf(row[j_smok])  # 0=Never, 1=Previous, 2=Current, -3=Prefer not to answer
            ever_smok = 0.0 if smok == 0.0 else (1.0 if smok in (1.0, 2.0) else np.nan)
            records[ieu] = (ratio, ever_smok, pf(row[j_py]), pf(row[j_bmi]), pf(row[j_ura]), pf(row[j_bil]), pf(row[j_crp]), copd, gout)
            if n_i % 100000 == 0:
                print(f"  [prep] scanned {n_i} rows", flush=True)

    df = pd.DataFrame.from_dict(
        records, orient="index",
        columns=["fev_fvc", "ever_smok", "pack_yrs", "bmi", "urate", "bil", "crp", "copd", "gout"]
    ).reset_index().rename(columns={"index": "IID"})
    df["FID"] = df["IID"]
    for k, kset in keeps.items():
        df[f"in_{k}"] = df["IID"].isin(kset).astype(int)

    df[["FID", "IID"]].to_csv(root / "keep/keep_union.txt", sep=" ", index=False, header=False)
    df.to_csv(root / "pheno/pheno_all.txt", sep="\t", index=False, na_rep="NA")

    # Separate binary pheno files for plink2 COPD & Gout GLMs
    df_fev = df[(df["in_fev"] == 1) & df["fev_fvc"].notna()].copy()
    df_fev[["FID", "IID"]].to_csv(root / "keep/keep_fev.txt", sep=" ", index=False, header=False)
    pd.DataFrame({"FID": df_fev["FID"], "IID": df_fev["IID"], "copd": df_fev["copd"] + 1}).to_csv(
        root / "pheno/pheno_copd.txt", sep="\t", index=False
    )
    df_ura = df[(df["in_urate"] == 1) & df["urate"].notna() & (df["urate"] > 0)].copy()
    df_ura[["FID", "IID"]].to_csv(root / "keep/keep_urate.txt", sep=" ", index=False, header=False)
    pd.DataFrame({"FID": df_ura["FID"], "IID": df_ura["IID"], "gout": df_ura["gout"] + 1}).to_csv(
        root / "pheno/pheno_gout.txt", sep="\t", index=False
    )

    # Write per-chromosome union SNP list (plus per-trait lists)
    all_snps = []
    for tag in ["fev", "urate", "bil", "crp"]:
        ldf = pd.read_csv(root / "inputs" / f"{tag}_leads.tsv", sep="\t")
        all_snps.append(ldf[["snp_id", "chr"]])
    uni_snps = pd.concat(all_snps, ignore_index=True).drop_duplicates("snp_id")
    for c, grp in uni_snps.groupby("chr"):
        (root / "snps" / f"union_chr{int(c):02d}.txt").write_text("\n".join(grp["snp_id"].astype(str)) + "\n")

    counts = [
        ("fev_fvc", "n_all", len(df_fev)),
        ("fev_fvc", "n_never_smoker", int((df_fev["ever_smok"] == 0).sum())),
        ("fev_fvc", "n_ever_smoker", int((df_fev["ever_smok"] == 1).sum())),
        ("fev_fvc", "cases_copd_j44", int(df_fev["copd"].sum())),
        ("urate", "n_all", len(df_ura)),
        ("urate", "cases_gout_m10", int(df_ura["gout"].sum())),
        ("bilirubin", "n_all", int(((df["in_bil"] == 1) & df["bil"].notna() & (df["bil"] > 0)).sum())),
        ("crp", "n_all", int(((df["in_crp"] == 1) & df["crp"].notna() & (df["crp"] > 0)).sum())),
        ("crp", "n_crp_le_10", int(((df["in_crp"] == 1) & df["crp"].notna() & (df["crp"] > 0) & (df["crp"] <= 10.0)).sum())),
        ("crp", "n_crp_gt_10", int(((df["in_crp"] == 1) & df["crp"].notna() & (df["crp"] > 10.0)).sum())),
    ]
    pd.DataFrame(counts, columns=["trait", "metric", "value"]).to_csv(root / "aggregate/sample_counts.csv", index=False)
    print(pd.DataFrame(counts, columns=["trait", "metric", "value"]).to_string(index=False), flush=True)
    print("[prep] done!", flush=True)

# ---- Reuse exact RIF-OLS + P11 moments + whitened PCA ----
def compute_rif_matrix(y, taus):
    n = len(y)
    qs = np.quantile(y, taus)
    bw = 1.06 * np.std(y, ddof=1) * (n ** (-0.2))
    rif = np.empty((n, len(taus)), dtype=np.float64)
    inv_sqrt2pi = 1.0 / math.sqrt(2.0 * math.pi)
    for j, (t, qt) in enumerate(zip(taus, qs)):
        z = (y - qt) / bw
        fj = np.mean(np.exp(-0.5 * z * z)) * inv_sqrt2pi / bw
        rif[:, j] = qt + (t - (y <= qt).astype(np.float64)) / fj
    return rif, qs

def rif_ols_leads(y, X_cov, G, taus):
    n = len(y)
    rif, qs = compute_rif_matrix(y, taus)
    Q, _ = np.linalg.qr(X_cov, mode="reduced")
    Psi = rif - Q @ (Q.T @ rif)
    Sigma0 = (Psi.T @ Psi) / n
    G_res = G - Q @ (Q.T @ G)
    d_resid = np.sum(G_res * G_res, axis=0)
    B = (G_res.T @ Psi) / d_resid[:, None]
    return B, Sigma0, d_resid, qs

def trap_weights(u):
    w = np.empty_like(u)
    w[0] = 0.5 * (u[1] - u[0])
    w[-1] = 0.5 * (u[-1] - u[-2])
    w[1:-1] = 0.5 * (u[2:] - u[:-2])
    return w

def moments_and_pca(B, Sigma0, d_resid, u, q, seed=20260923, n_pa=300):
    k, nt = B.shape
    cw = trap_weights(u)
    w_mu = cw; w_l2 = cw * (2.0 * u - 1.0); w_l3 = cw * (6.0 * u * u - 6.0 * u + 1.0)
    q_mu = float(np.sum(cw * q) / np.sum(cw))
    X_fr = np.column_stack([np.ones(nt), q - q_mu])
    w_sk = (np.eye(nt) - X_fr @ np.linalg.pinv(X_fr)).T @ w_l3

    def fe(w):
        est = B @ w; se = np.sqrt(float(w @ Sigma0 @ w) / d_resid); return est, se, est / se

    d_mu, se_mu, z_mu = fe(w_mu)
    d_l2, se_l2, z_l2 = fe(w_l2)
    d_sk, se_sk, z_sk = fe(w_sk)

    c1 = w_mu / math.sqrt(float(w_mu @ Sigma0 @ w_mu))
    w2o = w_l2 - float(c1 @ Sigma0 @ w_l2) * c1
    c2 = w2o / math.sqrt(float(w2o @ Sigma0 @ w2o))
    L = B @ c1; S = B @ c2
    r = float(np.sum(d_resid * L * S) / np.sum(d_resid * L * L))
    for _ in range(15):
        wi = d_resid / (1.0 + r * r)
        r_new = float(np.sum(wi * L * S) / np.sum(wi * L * L))
        if abs(r_new - r) < 1e-10: break
        r = r_new
    res = S - r * L; vi = (1.0 + r * r) / d_resid
    Q_het = float(np.sum(res * res / vi))
    p_het = float(stats.chi2.sf(Q_het, df=k - 1))
    tau_h = math.sqrt(max(0.0, (Q_het - (k - 1)) / float(np.sum((L * L) / vi))))
    het_z = res / np.sqrt(vi)

    eigvals, eigvecs = np.linalg.eigh(Sigma0)
    W_ih = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T
    Y = (B @ W_ih) * np.sqrt(d_resid)[:, None]
    norms = np.linalg.norm(Y, axis=1)
    U = Y / norms[:, None]
    u_svd, s_svd, vt_svd = np.linalg.svd(U, full_matrices=False)
    pct = 100.0 * (s_svd ** 2) / np.sum(s_svd ** 2)

    rng = np.random.default_rng(seed)
    n_comp = 0; pc2_ratio = np.nan; pc2_thr = np.nan
    for c_test in range(1, 5):
        signal = np.zeros_like(Y) if c_test == 1 else ((u_svd[:, :c_test - 1] * s_svd[:c_test - 1]) @ vt_svd[:c_test - 1, :]) * norms[:, None]
        sim_p = np.empty(n_pa)
        for rep in range(n_pa):
            Ys = signal + rng.standard_normal(Y.shape)
            Us = Ys / np.linalg.norm(Ys, axis=1, keepdims=True)
            ss = np.linalg.svd(Us, compute_uv=False)
            sim_p[rep] = 100.0 * (ss[c_test - 1] ** 2) / np.sum(ss ** 2)
        obs = float(pct[c_test - 1]); n_med = float(np.median(sim_p)); n_p95 = float(np.quantile(sim_p, 0.95))
        if c_test == 2:
            pc2_ratio = obs / n_med; pc2_thr = n_p95 / n_med
        if c_test == n_comp + 1 and obs > n_p95:
            n_comp = c_test

    pc1 = U @ vt_svd[0]; pc2 = U @ vt_svd[1]
    if np.median(pc1) < 0: pc1 = -pc1
    pl = pd.DataFrame({"delta_mu_z": z_mu, "delta_lambda2_z": z_l2, "skew_z": z_sk, "het_z": het_z, "pc1": pc1, "pc2": pc2})
    sm = {"n_leads": k, "n_components": n_comp, "tau_h": tau_h, "p_het": p_het,
          "pc1_pct": float(pct[0]), "pc2_pct": float(pct[1]), "pc3_pct": float(pct[2]),
          "pc2_ratio": pc2_ratio, "pc2_thr": pc2_thr,
          "n_skew_bonf": int(np.sum(np.abs(z_sk) > stats.norm.isf(0.025 / k))),
          "n_het_bonf": int(np.sum(np.abs(het_z) > stats.norm.isf(0.025 / k)))}
    return pl, sm

def load_union_dosages(geno_dir: Path, target_iids: np.ndarray, leads_df: pd.DataFrame):
    iid_to_pos = {iid: idx for idx, iid in enumerate(target_iids)}
    n = len(target_iids); k = len(leads_df)
    G = np.full((n, k), np.nan, dtype=np.float64)
    snp_to_col = {s: i for i, s in enumerate(leads_df["snp_id"].astype(str))}
    ea_map = dict(zip(leads_df["snp_id"].astype(str), leads_df["effect_allele"].astype(str)))
    oa_map = dict(zip(leads_df["snp_id"].astype(str), leads_df["other_allele"].astype(str)))
    for raw_path in sorted(geno_dir.glob("union_chr*.raw")):
        df_r = pd.read_csv(raw_path, sep=r"\s+", dtype={"FID": str, "IID": str})
        # Check if any columns needed
        wanted_cols = []
        for col in df_r.columns[6:]:
            snp = col.rsplit("_", 1)[0]
            if snp in snp_to_col:
                wanted_cols.append((col, snp, col.rsplit("_", 1)[1].split("(")[0]))
        if not wanted_cols:
            continue
        idx_keep = [i for i, iid in enumerate(df_r["IID"]) if iid in iid_to_pos]
        pos_dst = np.array([iid_to_pos[df_r["IID"].iat[i]] for i in idx_keep], dtype=int)
        for col, snp, counted in wanted_cols:
            vals = df_r[col].to_numpy(dtype=np.float64)[idx_keep]
            if np.isnan(vals).any():
                vals = np.where(np.isnan(vals), np.nanmean(vals), vals)
            ea = ea_map[snp]; oa = oa_map[snp]
            G[pos_dst, snp_to_col[snp]] = vals if counted == ea else (2.0 - vals)
    return G

def hc0_interaction(y, X_base, G, mod):
    Q, _ = np.linalg.qr(X_base, mode="reduced")
    yr = y - Q @ (Q.T @ y)
    k = G.shape[1]
    ib = np.empty(k); ise = np.empty(k); iz = np.empty(k)
    for j in range(k):
        gj = G[:, j]
        Z = np.column_stack([gj, gj * mod])
        Zr = Z - Q @ (Q.T @ Z)
        zi = np.linalg.inv(Zr.T @ Zr)
        bh = zi @ (Zr.T @ yr)
        e = yr - Zr @ bh
        V = zi @ ((Zr * e[:, None]).T @ (Zr * e[:, None])) @ zi
        ib[j] = bh[1]; ise[j] = math.sqrt(V[1, 1]); iz[j] = ib[j] / ise[j]
    return ib, ise, iz

def hc0_dominance(y, X_base, G):
    """Test dominance/recessiveness deviation: y ~ g + d_dom + covars, where d_dom = (g > 1.5) or (g - 1)^2."""
    Q, _ = np.linalg.qr(X_base, mode="reduced")
    yr = y - Q @ (Q.T @ y)
    k = G.shape[1]
    db = np.empty(k); dse = np.empty(k); dz = np.empty(k)
    for j in range(k):
        gj = G[:, j]
        # Quadratic deviation from additivity (positive = recessive / super-additive in homozygotes!)
        dj = (gj - 1.0) ** 2
        Z = np.column_stack([gj, dj])
        Zr = Z - Q @ (Q.T @ Z)
        zi = np.linalg.inv(Zr.T @ Zr)
        bh = zi @ (Zr.T @ yr)
        e = yr - Zr @ bh
        V = zi @ ((Zr * e[:, None]).T @ (Zr * e[:, None])) @ zi
        db[j] = bh[1]; dse[j] = math.sqrt(V[1, 1]); dz[j] = db[j] / dse[j]
    return db, dse, dz

def read_plink_glm(glm_dir: Path, pattern: str, leads_df: pd.DataFrame, prefix: str):
    dfs = []
    for p in sorted(glm_dir.glob(pattern)):
        df = pd.read_csv(p, sep="\t")
        if "TEST" in df.columns: df = df[df["TEST"] == "ADD"]
        dfs.append(df)
    g = pd.concat(dfs, ignore_index=True)
    ea_map = dict(zip(leads_df["snp_id"].astype(str), leads_df["effect_allele"].astype(str)))
    oa_map = dict(zip(leads_df["snp_id"].astype(str), leads_df["other_allele"].astype(str)))
    g = g[g["ID"].isin(ea_map)].copy()
    beta = np.log(g["OR"].astype(float)) if "LOG(OR)_SE" in g.columns else g["BETA"].astype(float)
    se   = g["LOG(OR)_SE"].astype(float) if "LOG(OR)_SE" in g.columns else g["SE"].astype(float)
    signs = np.array([1.0 if a == ea_map[s] else -1.0 for s, a in zip(g["ID"], g["A1"].astype(str))])
    return pd.DataFrame({"snp_id": g["ID"].astype(str), f"{prefix}_z": (beta * signs) / se, f"{prefix}_p": g["P"].astype(float)}).drop_duplicates("snp_id")

def run_analyze(root: Path):
    check_compute_node()
    taus = np.linspace(0.02, 0.98, 45)
    cov = pd.read_csv(COVAR_FILE, sep=r"\s+", dtype={"FID": str, "IID": str}).dropna()
    cov_cols = ["age", "sex"] + [f"PC{i}" for i in range(1, 11)]
    cov_nosex = ["age"] + [f"PC{i}" for i in range(1, 11)]
    ph = pd.read_csv(root / "pheno/pheno_all.txt", sep="\t", dtype={"FID": str, "IID": str}).merge(cov[["IID"] + cov_cols], on="IID", how="inner")

    all_arm_sums = []

    # ---- 1. FEV1/FVC: all vs never-smoker vs ever-smoker + G x smoking + COPD ----
    fev_leads = pd.read_csv(root / "inputs/fev_leads.tsv", sep="\t")
    sub = ph[(ph["in_fev"] == 1) & ph["fev_fvc"].notna()].reset_index(drop=True)
    G = load_union_dosages(root / "geno", sub["IID"].to_numpy(), fev_leads)
    copd_glm = read_plink_glm(root / "glm", "fev_chr*.copd.glm.logistic.hybrid", fev_leads, "copd")
    fev_out = fev_leads.merge(copd_glm, on="snp_id", how="left")

    # G x ever_smoker
    m_sm = sub["ever_smok"].notna().to_numpy()
    X_sm = np.column_stack([np.ones(int(m_sm.sum())), sub.loc[m_sm, cov_cols].to_numpy(dtype=np.float64), sub.loc[m_sm, "ever_smok"].to_numpy(dtype=np.float64)])
    _, _, iz_sm = hc0_interaction(sub.loc[m_sm, "fev_fvc"].to_numpy(dtype=np.float64), X_sm, G[m_sm, :], sub.loc[m_sm, "ever_smok"].to_numpy(dtype=np.float64))
    fev_out["smok_int_z"] = iz_sm

    for aname, mask in [("fev_all", np.ones(len(sub), dtype=bool)),
                        ("fev_never_smoker", (sub["ever_smok"] == 0).to_numpy()),
                        ("fev_ever_smoker",  (sub["ever_smok"] == 1).to_numpy())]:
        Xa = np.column_stack([np.ones(int(mask.sum())), sub.loc[mask, cov_cols].to_numpy(dtype=np.float64)])
        Ba, S0a, da, qa = rif_ols_leads(sub.loc[mask, "fev_fvc"].to_numpy(dtype=np.float64), Xa, G[mask, :], taus)
        pl, sm = moments_and_pca(Ba, S0a, da, taus, qa)
        sm["arm"] = aname; sm["n_samples"] = int(mask.sum()); all_arm_sums.append(sm)
        for c in pl.columns: fev_out[f"{aname}_{c}"] = pl[c].to_numpy()
    fev_out.to_csv(root / "aggregate/fev_smoking_per_lead.csv", index=False)

    # ---- 2. log Urate: all vs female vs male + G x sex + Gout ----
    ura_leads = pd.read_csv(root / "inputs/urate_leads.tsv", sep="\t")
    sub = ph[(ph["in_urate"] == 1) & ph["urate"].notna() & (ph["urate"] > 0)].reset_index(drop=True)
    y_u = np.log(sub["urate"].to_numpy(dtype=np.float64))
    G = load_union_dosages(root / "geno", sub["IID"].to_numpy(), ura_leads)
    gout_glm = read_plink_glm(root / "glm", "ura_chr*.gout.glm.logistic.hybrid", ura_leads, "gout")
    ura_out = ura_leads.merge(gout_glm, on="snp_id", how="left")
    X_u = np.column_stack([np.ones(len(sub)), sub[cov_cols].to_numpy(dtype=np.float64)])
    _, _, iz_sex = hc0_interaction(y_u, X_u, G, sub["sex"].to_numpy(dtype=np.float64))
    ura_out["sex_int_z"] = iz_sex

    for aname, mask, ccols in [("urate_all", np.ones(len(sub), dtype=bool), cov_cols),
                               ("urate_female", (sub["sex"] == 0).to_numpy(), cov_nosex),
                               ("urate_male",   (sub["sex"] == 1).to_numpy(), cov_nosex)]:
        Xa = np.column_stack([np.ones(int(mask.sum())), sub.loc[mask, ccols].to_numpy(dtype=np.float64)])
        Ba, S0a, da, qa = rif_ols_leads(y_u[mask], Xa, G[mask, :], taus)
        pl, sm = moments_and_pca(Ba, S0a, da, taus, qa)
        sm["arm"] = aname; sm["n_samples"] = int(mask.sum()); all_arm_sums.append(sm)
        for c in pl.columns: ura_out[f"{aname}_{c}"] = pl[c].to_numpy()
    ura_out.to_csv(root / "aggregate/urate_sex_per_lead.csv", index=False)

    # ---- 3. log Bilirubin: dominance deviation + excluding UGT1A1 homozygotes ----
    bil_leads = pd.read_csv(root / "inputs/bil_leads.tsv", sep="\t")
    sub = ph[(ph["in_bil"] == 1) & ph["bil"].notna() & (ph["bil"] > 0)].reset_index(drop=True)
    y_b = np.log(sub["bil"].to_numpy(dtype=np.float64))
    G = load_union_dosages(root / "geno", sub["IID"].to_numpy(), bil_leads)
    X_b = np.column_stack([np.ones(len(sub)), sub[cov_cols].to_numpy(dtype=np.float64)])
    _, _, dz_dom = hc0_dominance(y_b, X_b, G)
    bil_out = bil_leads.copy(); bil_out["dom_dev_z"] = dz_dom

    # Identify UGT1A1 lead rs11683974
    ugt_idx = int(np.where(bil_leads["snp_id"].astype(str) == "rs11683974")[0][0])
    no_ugt_hom = (G[:, ugt_idx] < 1.5)
    for aname, mask in [("bil_all", np.ones(len(sub), dtype=bool)),
                        ("bil_no_ugt1a1_hom", no_ugt_hom)]:
        Xa = np.column_stack([np.ones(int(mask.sum())), sub.loc[mask, cov_cols].to_numpy(dtype=np.float64)])
        Ba, S0a, da, qa = rif_ols_leads(y_b[mask], Xa, G[mask, :], taus)
        pl, sm = moments_and_pca(Ba, S0a, da, taus, qa)
        sm["arm"] = aname; sm["n_samples"] = int(mask.sum()); all_arm_sums.append(sm)
        for c in pl.columns: bil_out[f"{aname}_{c}"] = pl[c].to_numpy()
    bil_out.to_csv(root / "aggregate/bilirubin_dom_per_lead.csv", index=False)

    # ---- 4. log CRP: all vs CRP <= 10 mg/L + G x BMI ----
    crp_leads = pd.read_csv(root / "inputs/crp_leads.tsv", sep="\t")
    sub = ph[(ph["in_crp"] == 1) & ph["crp"].notna() & (ph["crp"] > 0)].reset_index(drop=True)
    y_c = np.log(sub["crp"].to_numpy(dtype=np.float64))
    G = load_union_dosages(root / "geno", sub["IID"].to_numpy(), crp_leads)
    crp_out = crp_leads.copy()
    m_bmi = sub["bmi"].notna().to_numpy()
    bmi_s = (sub.loc[m_bmi, "bmi"].to_numpy(dtype=np.float64) - sub.loc[m_bmi, "bmi"].mean()) / sub.loc[m_bmi, "bmi"].std()
    X_cb = np.column_stack([np.ones(int(m_bmi.sum())), sub.loc[m_bmi, cov_cols].to_numpy(dtype=np.float64), bmi_s])
    _, _, iz_bmi = hc0_interaction(y_c[m_bmi], X_cb, G[m_bmi, :], bmi_s)
    crp_out["bmi_int_z"] = iz_bmi

    for aname, mask in [("crp_all", np.ones(len(sub), dtype=bool)),
                        ("crp_le_10", (sub["crp"] <= 10.0).to_numpy())]:
        Xa = np.column_stack([np.ones(int(mask.sum())), sub.loc[mask, cov_cols].to_numpy(dtype=np.float64)])
        Ba, S0a, da, qa = rif_ols_leads(y_c[mask], Xa, G[mask, :], taus)
        pl, sm = moments_and_pca(Ba, S0a, da, taus, qa)
        sm["arm"] = aname; sm["n_samples"] = int(mask.sum()); all_arm_sums.append(sm)
        for c in pl.columns: crp_out[f"{aname}_{c}"] = pl[c].to_numpy()
    crp_out.to_csv(root / "aggregate/crp_acute_per_lead.csv", index=False)

    sum_df = pd.DataFrame(all_arm_sums)
    sum_df.to_csv(root / "aggregate/more_traits_arm_summary.csv", index=False)
    print("\n=== MORE TRAITS ARM SUMMARY ===", flush=True)
    print(sum_df.to_string(index=False), flush=True)

    # Correlation summary across all 4 traits
    corr_rows = []
    for tname, df_t, pairs in [
        ("fev_fvc", fev_out, [("pc2_score", "smok_int_z"), ("het_resid_z", "smok_int_z"), ("pc2_score", "copd_z"), ("het_resid_z", "copd_z")]),
        ("urate",   ura_out, [("pc2_score", "sex_int_z"), ("het_resid_z", "sex_int_z"), ("pc2_score", "gout_z"), ("het_resid_z", "gout_z")]),
        ("bilirubin", bil_out, [("pc2_score", "dom_dev_z"), ("het_resid_z", "dom_dev_z"), ("skew_beyond_frame_z", "dom_dev_z")]),
        ("crp",     crp_out, [("pc2_score", "bmi_int_z"), ("het_resid_z", "bmi_int_z")]),
    ]:
        for x, y in pairs:
            s = df_t[[x, y]].dropna()
            rp, pp = stats.pearsonr(s[x], s[y])
            rs, ps = stats.spearmanr(s[x], s[y])
            corr_rows.append({"trait": tname, "x": x, "y": y, "n": len(s), "pearson_r": rp, "pearson_p": pp, "spearman_rho": rs, "spearman_p": ps})
    cdf = pd.DataFrame(corr_rows)
    cdf.to_csv(root / "aggregate/more_traits_correlations.csv", index=False)
    print("\n=== MORE TRAITS CORRELATIONS ===", flush=True)
    print(cdf.to_string(index=False), flush=True)
    print("[analyze] ALL DONE!", flush=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["prep", "analyze"])
    ap.add_argument("--root", required=True)
    args = ap.parse_args()
    configure_paths()
    if args.mode == "prep":
        run_prep(Path(args.root))
    else:
        run_analyze(Path(args.root))

if __name__ == "__main__":
    main()
