#!/usr/bin/env Rscript
# Import the simulation summaries
#
# What: the simulation studies (Stage 1 equivalence and speed, Stage 2 effect recovery,
#   mixture recovery against full mixture maximum likelihood, mixture speed, pooled versus
#   jackknife standard errors, and the standard-error check at the empirical mixture
#   separations) were run once with the G-MoQ packages; their per-cell summaries are kept
#   in inputs/. This script copies them unchanged into results/, adds the mixture timing
#   rows, and derives the one-sentence summary of the low-separation standard-error check
#   quoted in the Supplementary Information, so that the sentence and the table agree.
# Paper: Supplementary Figs. S1-S6; Supplementary Tables S1, S2
# Input: analysis/08_simulations/inputs/ (aggregate per-cell summaries);
#   results/03_testosterone/recovery/timings.csv (timings of one panel of 100 variants)
# Output: results/08_simulations/summary/{table_s1_design_provenance,
#   table_s2_numerical_summary, s2_mean_timings, s5_mean_timings, mixture_timings,
#   lowsep_empirical_se_ratio, lowsep_lowsep_se_summary, lowsep_reproducibility_checks,
#   lowsep_separations, lowsep_sim_summary, lowsep_citation}.csv
# Runs: locally from results/ in about 2 seconds

.repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(.repo, "R", "load.R"))

inputs <- repo_path("analysis", "08_simulations", "inputs")
out <- "08_simulations/summary"

# Counts in these tables are simulation replicates and timing repeats, not participants, so
# the small-count check of write_result() does not apply.
save_summary <- function(x, name) write_result(x, out, name, count_cols = character(0))

## 1. Simulation studies 1-6: copy the per-cell summaries unchanged
# Table S1 (design of each study), Table S2 (bias, RMSE, coverage, timing per cell) and the
# mean timings plotted in Supplementary Figs. S2 and S5.
copy_summary <- function(file, name) {
  x <- fread(file.path(inputs, "simulation_summaries", file))
  save_summary(x, name)
}
copy_summary("TableS1_design_provenance.csv", "table_s1_design_provenance")
copy_summary("TableS2_numerical_summary.csv", "table_s2_numerical_summary")
copy_summary("S2_mean_timings.csv", "s2_mean_timings")
copy_summary("S5_mean_timings.csv", "s5_mean_timings")

## 2. Mixture timings: G-MoQ against full mixture maximum likelihood
# Median time to fit the estimates and their standard errors for one panel of 100 variants,
# over three panels per sample size (G-MoQ with the pooled or the jackknife covariance, and
# mixture maximum likelihood refitted at every variant).
timings <- as.data.table(read_result("03_testosterone/recovery", "timings"))
mixture_timings <- timings[method %in% c("gmoq_pooled", "mixture_ml", "gmoq_jackknife"),
  .(study = "mixture_timing", n, M, method, n_panels, median_fit_seconds, min_fit_seconds,
    max_fit_seconds, median_wall_seconds, ratio_to_gmoq_pooled)]
setorder(mixture_timings, n, method)
save_summary(mixture_timings, "mixture_timings")

## 3. Study 7: standard errors at the empirical mixture separations
# Pooled against 1,000-block jackknife standard errors of the mixture parameters, on the
# real data (empirical_se_ratio) and in simulation (sim_summary), combined per case and
# parameter in lowsep_se_summary (see lowsep_se_check/).
for (file in c("empirical_se_ratio", "lowsep_se_summary", "reproducibility_checks",
               "separations", "sim_summary")) {
  x <- fread(file.path(inputs, "lowsep_se_check", paste0(file, ".csv")))
  save_summary(x, paste0("lowsep_", file))
}

## 4. The one-sentence summary of study 7, recomputed from the combined table
lowsep <- fread(file.path(inputs, "lowsep_se_check", "lowsep_se_summary.csv"))
empirical <- unique(lowsep[, .(case, parameter, separation_fitted, emp_pooled_over_jk_median)])
null_cells <- lowsep[sim_effect == "null"]
citation <- data.table(
  separation_min = min(empirical$separation_fitted),
  separation_max = max(empirical$separation_fitted),
  emp_median_ratio_min = min(empirical$emp_pooled_over_jk_median, na.rm = TRUE),
  emp_median_ratio_max = max(empirical$emp_pooled_over_jk_median, na.rm = TRUE),
  sim_null_cells = nrow(null_cells),
  sim_null_pooled_coverage_min = min(null_cells$sim_pooled_coverage95),
  sim_null_pooled_coverage_max = max(null_cells$sim_pooled_coverage95),
  sim_null_jk_coverage_min = min(null_cells$sim_jk_coverage95),
  sim_null_jk_coverage_max = max(null_cells$sim_jk_coverage95))
citation[, citation := sprintf(paste0(
  "At the mixture separations of our empirical cases (%.1f-%.1f within-component SD), pooled and ",
  "1,000-block jackknife SEs for the membership and component-location parameters agreed to within 1%% (median ratio ",
  "%.3f-%.3f genome-wide; simulated 95%% coverage %.3f-%.3f under the null for both), the only departure being the ",
  "upper-component location at low MAF with a large true effect at N = 10,000, which is independent of separation."),
  separation_min, separation_max, emp_median_ratio_min, emp_median_ratio_max,
  min(sim_null_pooled_coverage_min, sim_null_jk_coverage_min),
  max(sim_null_pooled_coverage_max, sim_null_jk_coverage_max))]
save_summary(citation, "lowsep_citation")
