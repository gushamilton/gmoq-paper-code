#!/usr/bin/env python3
"""Prepare aggregate lead metadata for the two original biological-correlate producers.

Run from the repository root. --root is the corresponding private analysis output
root on the cluster; the files written here contain variant summaries only.
"""
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", choices=["hidden", "more"], required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--moments", type=Path, default=Path("results/06_traits/moments/lead_moments.csv"))
    args = parser.parse_args()
    groups = {"hidden": {"alt": ("log_alt", 316), "ldl": ("ldl", 230)},
              "more": {"fev": ("fev1_fvc_ratio", 421), "urate": ("log_urate", 396),
                       "bil": ("log_bilirubin_total", 183), "crp": ("log_crp", 389)}}[args.analysis]
    columns = ["snp_id", "chr", "bp", "effect_allele", "other_allele", "delta_mu", "delta_mu_z",
               "delta_lambda2", "delta_lambda2_z", "skew_beyond_frame_z", "het_resid_z",
               "pc1_score", "pc2_score", "pc3_score"]
    with args.moments.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    out = args.root / "inputs"
    out.mkdir(parents=True, exist_ok=True)
    for tag, (arm, count) in groups.items():
        subset = [r for r in rows if r["arm"] == arm]
        if len(subset) != count or len({r["snp_id"] for r in subset}) != count:
            raise ValueError(f"Unexpected lead set for {arm}")
        with (out / f"{tag}_leads.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(subset)


if __name__ == "__main__":
    main()
