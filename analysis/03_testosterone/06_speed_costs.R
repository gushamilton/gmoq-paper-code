#!/usr/bin/env Rscript
# Cost of G-MoQ and of direct mixture maximum likelihood, part 3: cost per variant
#
# What: turns the component timings of parts 1 and 2 into single-core costs per additional
#   variant at each sample size: Stage 1 (the difference between Stage 1 panels of 300 and 100
#   variants, per variant), Stage 2 given the Stage 1 output (param_gwas on 100,000 variants,
#   per variant) and direct mixture maximum likelihood (mean over the variants timed); and the
#   implied cost of 1.2 million variants. One-off costs per trait (baseline mixture fit,
#   quantiles and kernel density, residualisation, pooled covariance) are kept separate and
#   added only to the genome-wide Stage 1 total.
# Paper: Fig. 2c; Supplementary Fig. S5d (1.2 million variants at N = 400,000: direct ML
#   3,354 CPU-hours, G-MoQ Stages 1 and 2 5.3 CPU-hours, Stage 2 alone 3 CPU-minutes)
# Input: with a run directory of parts 1 and 2 as argument, its python_timings.csv and
#   stage2_timings_*.csv; without one, the same timings already in
#   results/03_testosterone/speed/{python_component_timings,stage2_timings}.csv
# Output: results/03_testosterone/speed/per_snp_costs.csv (and, from a run directory, the two
#   timing tables above)
# Runs: locally from results/ in about 2 seconds

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
out_dir <- "03_testosterone/speed"
run <- commandArgs(trailingOnly = TRUE)[1]

## 1. Read the timings ---------------------------------------------------------------------
if (!is.na(run)) {
  py <- fread(file.path(run, "python_timings.csv"))
  s2 <- rbindlist(lapply(sort(list.files(run, pattern = "^stage2_timings_.*\\.csv$", full.names = TRUE)), fread), fill = TRUE)
  write_result(py, out_dir, "python_component_timings")
  write_result(s2, out_dir, "stage2_timings")
} else {
  py <- as.data.table(read_result(out_dir, "python_component_timings"))
  s2 <- as.data.table(read_result(out_dir, "stage2_timings"))
}

## 2. One-off Stage 1 costs per trait (mean over repeats, 100-variant panel) ---------------
one_off <- py[component %in% c("g_env_core_load", "g_baseline_mixture_fit", "g_compute_rif_total",
                               "g_residualize", "g_sigma0"),
              .(mean_seconds = mean(seconds)), by = .(n, component)][, .(stage1_one_off = sum(mean_seconds)), by = n]

## 3. Marginal cost per variant -----------------------------------------------------------
# Stage 1: (time for 300 variants - time for 100 variants) / 200.
st1 <- dcast(py[component == "g_stage1_ols_panel", .(s = mean(seconds)), by = .(n, M)], n ~ M, value.var = "s")
m_sizes <- setdiff(names(st1), "n")
st1[, stage1_per_snp := (get(m_sizes[2]) - get(m_sizes[1])) / (as.numeric(m_sizes[2]) - as.numeric(m_sizes[1]))]
# Direct mixture maximum likelihood, one variant per fit.
ml <- py[component == "ml_fit_one_snp", .(ml_per_snp = mean(seconds), ml_snps_timed = .N, ml_ok = sum(note == "ok")), by = n]
# Stage 2 on 100,000 variants, with the package version the benchmark used.
s2m <- s2[label == "pin9c47e7b" & M == 100000 &
            component %in% c("r_param_gwas_bench_args", "r_param_gwas_chisq_naive", "r_bare_projection_plus_se"),
          .(per_snp = seconds / M), by = .(n, component)]
s2w <- dcast(s2m, n ~ component, value.var = "per_snp")
setnames(s2w, c("r_param_gwas_bench_args", "r_param_gwas_chisq_naive", "r_bare_projection_plus_se"),
         c("stage2_per_snp_param_gwas", "stage2_per_snp_no_q_calibration", "stage2_per_snp_bare_projection"))

## 4. Genome-wide projections (1.2 million variants, one core) ----------------------------
costs <- Reduce(function(x, y) merge(x, y, by = "n"), list(st1[, .(n, stage1_per_snp)], ml, s2w, one_off))
n_variants <- 1.2e6
costs[, `:=`(ml_over_stage2_per_snp = ml_per_snp / stage2_per_snp_param_gwas,
             ml_over_stage1_plus_stage2_per_snp = ml_per_snp / (stage1_per_snp + stage2_per_snp_param_gwas),
             genome_1p2M_ml_cpu_hours = n_variants * ml_per_snp / 3600,
             genome_1p2M_stage1_cpu_hours = (n_variants * stage1_per_snp + stage1_one_off) / 3600,
             genome_1p2M_stage2_cpu_seconds = n_variants * stage2_per_snp_param_gwas)]
write_result(costs, out_dir, "per_snp_costs")
