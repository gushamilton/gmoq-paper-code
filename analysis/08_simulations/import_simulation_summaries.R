#!/usr/bin/env Rscript
# Import the simulation summaries
#
# What: the simulation studies (Stage 1 equivalence and speed, Stage 2 effect recovery,
#   mixture recovery against full mixture maximum likelihood, mixture speed, pooled versus
#   jackknife standard errors, and the standard-error check at the empirical mixture
#   separations) were run once with the G-MoQ packages; their per-cell summaries are kept
#   in inputs/. This script imports them, replaces cost rows with the later timing
#   benchmarks, and derives the one-sentence summary of the low-separation standard-error check
#   quoted in the Supplementary Information, so that the sentence and the table agree.
# Paper: Supplementary Figs. S1-S6; Supplementary Tables S1, S2
# Input: analysis/08_simulations/inputs/ (aggregate per-cell summaries);
#   results/08_simulations/speed_benchmark/ (later single-thread timing benchmarks)
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
# Table S1 records each study's own source. Table S2 retains historical recovery
# and coverage rows, with the later costs used in Supplementary Figs. S2 and S5.
design <- fread(file.path(inputs, "simulation_summaries", "TableS1_design_provenance.csv"))
timing_pin <- jsonlite::fromJSON(repo_path("results", "08_simulations", "speed_benchmark", "software_pin.json"))$commit
design[study %in% c("S2", "S5"), `:=`(package_pin = timing_pin, source = "results/08_simulations/speed_benchmark")]
design[study == "S2", `:=`(N = "2000/400000", M_MAF = "M=1/256/1000")]
design[study == "S5", `:=`(DGP = "Independent balanced null-mixture panels", M_MAF = "MAF .30; M=100 panel",
                          effect_range = "null; shared panel genotypes", tau45 = "yes")]
save_summary(design, "table_s1_design_provenance")
# Recovery/coverage rows retain the historical simulations. Cost rows use the later
# benchmark throughout; do not infer timing from the effect-recovery runs.
summary_rows <- fread(file.path(inputs, "simulation_summaries", "TableS2_numerical_summary.csv"))
stage1 <- as.data.table(read_result("08_simulations/speed_benchmark", "stage1_panel_times"))
panel <- as.data.table(read_result("08_simulations/speed_benchmark", "mixture_panel_times"))[
  method %in% c("gmoq_pooled", "gmoq_jackknife", "mixture_ml")]
stopifnot(nrow(stage1) == 36L, all(stage1$count == 3L), nrow(panel) == 9L,
          all(panel$variants_per_panel == 100L), all(panel$panels == 3L), all(panel$successful_fits == 300L))
mixture_means <- rbind(
  panel[, .(n, M = variants_per_panel, method, metric = "fit + SE seconds", mean_value = mean_fit_seconds, count = panels)],
  panel[, .(n, M = variants_per_panel, method, metric = "whole-process seconds", mean_value = mean_wall_seconds, count = panels)],
  panel[, .(n, M = variants_per_panel, method, metric = "peak RSS MiB", mean_value = mean_peak_rss_mib, count = panels)])
save_summary(stage1, "s2_mean_timings")
save_summary(mixture_means, "s5_mean_timings")
metric_columns <- c("fit + SE seconds" = "timing_fit_seconds", "whole-process seconds" = "timing_wall_seconds",
                    "peak RSS MiB" = "peak_rss_mib")
for (i in which(summary_rows$study %in% c("S2", "S5"))) {
  row <- summary_rows[i]
  t <- if (row$study == "S2") stage1[workload == row$workload & method == row$method] else
    mixture_means[n == row$n & method == row$method]
  stopifnot(nrow(t) == 3L, identical(sort(t$metric), sort(names(metric_columns))))
  for (metric_name in names(metric_columns)) set(summary_rows, i, metric_columns[[metric_name]], t[metric == metric_name, mean_value])
  set(summary_rows, i, "source", "results/08_simulations/speed_benchmark")
}
save_summary(summary_rows, "table_s2_numerical_summary")

## 2. Mixture timings: G-MoQ against full mixture maximum likelihood
# Median time to fit the estimates and their standard errors for one panel of 100 variants,
# over three panels per sample size (G-MoQ with the pooled or the jackknife covariance, and
# mixture maximum likelihood refitted at every variant).
timings <- as.data.table(read_result("08_simulations/speed_benchmark", "mixture_panel_repeats"))
stopifnot(nrow(timings) == 36L, all(timings$status == "complete"), all(timings$M == 100L),
          all(timings$successful_snp_fits == timings$attempted_snp_fits),
          all(tolower(as.character(timings$primary_wall_eligible)) == "true"))
mixture_timings <- timings[method %in% c("gmoq_pooled", "mixture_ml", "gmoq_jackknife"),
  .(study = "mixture_timing", n_panels = .N, median_fit_seconds = median(fit_seconds),
    min_fit_seconds = min(fit_seconds), max_fit_seconds = max(fit_seconds),
    median_wall_seconds = median(process_wall_seconds)), by = .(n, M, method)]
mixture_timings[, ratio_to_gmoq_pooled := median_fit_seconds / median_fit_seconds[method == "gmoq_pooled"], by = n]
setorder(mixture_timings, n, method)
save_summary(mixture_timings, "mixture_timings")

## 3. Study 7: standard errors at the empirical mixture separations
# Pooled against 1,000-block jackknife standard errors of the mixture parameters, on the
# real data (empirical_se_ratio) and in simulation (sim_summary), combined per case and
# parameter in lowsep_se_summary (see lowsep_se_check/).
for (file in c("empirical_se_ratio", "lowsep_se_summary", "reproducibility_checks",
               "separations", "sim_summary")) {
  x <- fread(file.path(inputs, "lowsep_se_check", paste0(file, ".csv")))
  name <- paste0("lowsep_", file)
  path <- file.path(Sys.getenv("GMOQ_OUT", repo_path("results")), out, paste0(name, ".csv"))
  # Keep the frozen unrelated export when importing only changes binary-rounding
  # at the last printed digit; changed inputs still regenerate their summary.
  unchanged <- file.exists(path) && isTRUE(all.equal(as.data.frame(x), utils::read.csv(path, check.names = FALSE),
                                                    check.attributes = FALSE, tolerance = 1e-12))
  if (!unchanged) save_summary(x, name)
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
