"""Summarize completed one-core timing runs without selecting successful variants."""

import csv
import json
import statistics
import sys
from pathlib import Path


def read_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def mean(rows, field):
    return statistics.mean(float(row[field]) for row in rows)


def main(root, out):
    out.mkdir(parents=True, exist_ok=True)
    sizes = [10000, 50000, 400000]
    methods = ["gmoq_pooled", "gmoq_jackknife", "mixture_ml", "mixture_ml_warm"]
    mixture_rows, checks, comparisons = [], [], []
    for n in sizes:
        raw, numerical = [], []
        for rep in (1, 2, 3):
            cell = root / "mixture-full" / f"n{n}_rep{rep}"
            rows = read_csv(cell / "timings.csv")
            assert sorted(row["method"] for row in rows) == sorted(methods)
            assert all(int(row["M"]) == 100 and int(row["rep"]) == rep and int(row["n"]) == n for row in rows)
            raw.extend(rows)
            q = read_csv(cell / "warm_quality.csv")
            numerical.extend(q)
            checks.append(dict(n=n, rep=rep, compared=len(q),
                passed=sum(row["passed"] == "True" for row in q),
                fallbacks=sum(row["warm_fallback"] == "True" for row in q),
                max_estimate_difference_in_se=max((float(row["max_estimate_difference_in_canonical_se"]) for row in q), default=None),
                max_relative_se_difference=max((float(row["max_relative_se_difference"]) for row in q), default=None)))
        group = {m: [row for row in raw if row["method"] == m] for m in methods}
        for m, rows in group.items():
            assert len(rows) == 3
            mixture_rows.append(dict(n=n, variants_per_panel=100, panels=3, method=m,
                successful_fits=sum(int(row["successful_snp_fits"]) for row in rows),
                mean_fit_seconds=mean(rows, "fit_seconds"),
                mean_wall_seconds=mean(rows, "process_wall_seconds"),
                min_wall_seconds=min(float(row["process_wall_seconds"]) for row in rows),
                max_wall_seconds=max(float(row["process_wall_seconds"]) for row in rows),
                mean_peak_rss_mib=mean(rows,"peak_rss_kib")/1024))
        for m in ("mixture_ml", "mixture_ml_warm"):
            full = all(int(row["successful_snp_fits"]) == 100 for row in group[m]+group["gmoq_pooled"])
            equivalent = len(numerical) == 300 and all(row["passed"] == "True" for row in numerical)
            eligible = full and (m != "mixture_ml_warm" or equivalent)
            comparisons.append(dict(n=n, reference=m, panels=3, variants_per_panel=100, eligible=eligible,
                reference_over_pooled_fit=mean(group[m],"fit_seconds")/mean(group["gmoq_pooled"],"fit_seconds") if eligible else None,
                reference_over_pooled_wall=mean(group[m],"process_wall_seconds")/mean(group["gmoq_pooled"],"process_wall_seconds") if eligible else None))
    write_csv(out/"mixture_panel_times.csv",mixture_rows)
    write_csv(out/"mixture_panel_speedups.csv",comparisons)
    write_csv(out/"mixture_comparator_checks.csv",checks)
    source = [root/"rif-smoke-20445185",root/"rif-n400000_m1-20445189",root/"rif-n400000_m1000-20445190"]
    rif = []
    for folder in source:
        assert (folder/"COMPLETE").is_file()
        assert json.loads((folder/"tables/verification.json").read_text())["status"]=="verified"
        rif.extend(read_csv(folder/"tables/S2_mean_timings.csv"))
    assert len(rif)==36
    write_csv(out/"stage1_panel_times.csv",rif)


if __name__ == "__main__":
    main(*(Path(arg) for arg in sys.argv[1:]))
