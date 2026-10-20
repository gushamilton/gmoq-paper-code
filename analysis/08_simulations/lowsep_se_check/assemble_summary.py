#!/usr/bin/env python3
"""Standard errors at the empirical mixture separations: combined table.

What: joins the empirical arm (empirical_se_ratio.csv, separations.csv from
empirical_se_ratio.R) and the simulation arm (sim_summary.csv from sim_lowsep_se.py) into
one row per empirical case x Stage 2 parameter x simulated (MAF, effect) cell: fitted
separation, the pooled/jackknife SE ratio on the real Stage 2 output, and the simulated
pooled/jackknife ratio, SE / empirical SD and 95% coverage. Each empirical case is paired
with the simulated geometry of the same separation (testosterone at N = 400,000,
antibodies at N = 10,000).
Paper: Supplementary Information, study 7; Supplementary Table S2
Input: the aggregate CSVs above
Output: lowsep_se_summary.csv
Runs: anywhere, in under a second
"""
import argparse
import csv
from pathlib import Path

CASES = [  # empirical case, preferred empirical arm, sim geometry, sim N
    ("testosterone", "corrected_weights", "testosterone_unshrunk", 400000),
    ("testosterone_shrink_0.2", "corrected_weights", "testosterone_shrink_0.2", 400000),
    ("testosterone_shrink_0.5", "corrected_weights", "testosterone_shrink_0.5", 400000),
    ("testosterone_shrink_0.7", "corrected_weights", "testosterone_shrink_0.7", 400000),
    ("jcv_vp1_log_all", "scan_weights", "jcv_vp1", 10000),
    ("mcv_vp1_log_all", "scan_weights", "mcv_vp1", 10000),
]
SIM_PARAM = {"gamma": "gamma", "beta_left": "beta_1", "beta_right": "beta_2"}
OLD_ALIAS = {"jcv_vp1_log_all": "jcv_log_all", "mcv_vp1_log_all": "mcv_log_all"}


def read(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--empirical-dir", type=Path, required=True)
    ap.add_argument("--sim-summary", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    emp = read(a.empirical_dir / "empirical_se_ratio.csv")
    sep = read(a.empirical_dir / "separations.csv")
    sim = read(a.sim_summary)

    def e(case, arm, par, stratum):
        for r in emp:
            if (r["case"], r["arm"], r["parameter"], r["z_stratum"]) == (case, arm, par, stratum):
                return r
        return {}

    rows = []
    for case, arm, geom, n in CASES:
        s = next(r for r in sep if r["case"] == case and r["weight_source"] != "old_stage1_rif_density")
        old_case = OLD_ALIAS.get(case, case)
        for par, spar in SIM_PARAM.items():
            e0, e5 = e(case, arm, par, "|z_pooled|>0"), e(case, arm, par, "|z_pooled|>5")
            g0 = e(old_case, "old_weights_genome", par, "|z_pooled|>0")
            for c in sim:
                if c["geometry"] != geom or int(c["n"]) != n or c["parameter"] != spar:
                    continue
                rows.append({
                    "case": case, "parameter": par,
                    "separation_fitted": round(float(s["separation"]), 3),
                    "pi_upper_fitted": round(float(s["pi_upper"]), 3),
                    "n_samples": s["n_samples"],
                    "emp_arm": arm, "emp_n_snps": e0.get("n_snps"),
                    "emp_pooled_over_jk_median": e0.get("median_pooled_over_jk"),
                    "emp_pooled_over_jk_q05": e0.get("q05"),
                    "emp_pooled_over_jk_q95": e0.get("q95"),
                    "emp_z5_n_snps": e5.get("n_snps", 0),
                    "emp_z5_pooled_over_jk_median": e5.get("median_pooled_over_jk"),
                    "emp_genome_oldw_case": old_case,
                    "emp_genome_oldw_pooled_over_jk_median": g0.get("median_pooled_over_jk"),
                    "sim_separation": c["separation"], "sim_pi_upper": c["pi_upper"],
                    "sim_n": c["n"], "sim_maf": c["maf"], "sim_effect": c["effect"],
                    "sim_replicates": c["replicates_ok"], "sim_failed": c["replicates_failed"],
                    "sim_pooled_over_jk_median": c["median_pooled_over_jackknife_se"],
                    "sim_pooled_rmsse_over_empsd": c["pooled_rms_se_over_empirical_sd"],
                    "sim_jk_rmsse_over_empsd": c["jackknife_rms_se_over_empirical_sd"],
                    "sim_empsd_rel_mcse": c["empirical_sd_rel_mcse"],
                    "sim_pooled_coverage95": c["pooled_coverage95"],
                    "sim_jk_coverage95": c["jackknife_coverage95"],
                    "sim_coverage_mcse": c["pooled_coverage95_mcse"],
                    "sim_bias": c["bias"], "sim_empirical_sd": c["empirical_sd"],
                })
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {a.out}")


if __name__ == "__main__":
    main()
