"""Summarize the original per-variant timing calculation on the current package."""

import csv
import statistics
import sys
from pathlib import Path


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def seconds(rows, component, size=None):
    return [
        float(row["seconds"])
        for row in rows
        if row["component"] == component
        and (size is None or float(row["M"]) == size)
    ]


def main(run, dest):
    py = read_table(run / "python_timings.csv")
    r = read_table(run / "stage2_timings_gmoq166635d4.csv")
    out = []
    for n in (10000, 50000, 400000):
        p = [row for row in py if int(row["n"]) == n]
        s = [row for row in r if int(float(row["n"])) == n]
        times = seconds(p, "ml_fit_one_snp")
        direct = statistics.mean(times)
        fit_100 = statistics.mean(seconds(p, "g_stage1_ols_panel", 100))
        fit_300 = statistics.mean(seconds(p, "g_stage1_ols_panel", 300))
        stage1 = (fit_300 - fit_100) / 200
        stage2 = seconds(s, "r_param_gwas_bench_args", 100000)[0] / 100000
        assert min(direct, stage1, stage2) > 0
        out.append(dict(
            n=n, direct_six_start_snps=len(times), direct_seconds_per_variant=direct,
            stage1_seconds_per_extra_variant=stage1,
            stage2_seconds_per_variant_given_summaries=stage2,
            direct_over_marginal_both_stages=direct / (stage1 + stage2),
            direct_over_stage2_given_summaries=direct / stage2,
            measured_stage2_seconds_1200000=(seconds(s, "r_param_gwas_bench_args", 1200000) or [""])[0],
        ))
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out[0])
        writer.writeheader()
        writer.writerows(out)


if __name__ == "__main__":
    main(*(Path(arg) for arg in sys.argv[1:]))
