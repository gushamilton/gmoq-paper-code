#!/usr/bin/env Rscript
# Supplementary Table S2: numerical results of the simulation studies
#
# What: the full numerical summary behind Supplementary Figs S1-S7, one row per study,
#   scenario, sample size, allele frequency, parameter, method and covariance (design in
#   Table S1). Two kinds of row are added to the stored summary: the median time per
#   100-variant mixture panel for G-MoQ against direct mixture maximum likelihood, and for
#   study 7 the pooled/jackknife SE ratio in the empirical fits with the range of null coverage
#   in the matching simulated cells.
# Paper: Supplementary Table S2; Supplementary Note 1
# Input: results/08_simulations/summary/{table_s2_numerical_summary,mixture_timings,
#   lowsep_lowsep_se_summary}.csv
# Output: paper_assets/tableS2_simulation_results.csv, _caption.md
# Runs: locally from results/ in about 1 second

suppressPackageStartupMessages({library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load ----------------------------------------------------------------------------------
summary_rows <- as.data.table(read_result("08_simulations/summary", "table_s2_numerical_summary"))
timings <- as.data.table(read_result("08_simulations/summary", "mixture_timings"))
lowsep <- as.data.table(read_result("08_simulations/summary", "lowsep_lowsep_se_summary"))

## 2. Mixture timing rows: G-MoQ pooled against direct mixture ML -------------------------
timing_row <- function(n_val) {
  t <- timings[n == n_val & M == 100]
  data.table(
    study = "mixture timing",
    quantity = paste0("median fit seconds, G-MoQ pooled vs direct mixture ML (N = ",
                      format(n_val, big.mark = ",", scientific = FALSE), ", 100-SNP panels)"),
    n = n_val,
    timing_fit_seconds = paste0(round(t[method == "gmoq_pooled", median_fit_seconds][1], 2), " vs ",
                                round(t[method == "mixture_ml", median_fit_seconds][1], 1)))
}
timing_rows <- rbind(timing_row(10000), timing_row(400000))

## 3. Study 7 rows: pooled against jackknife SE in the empirical fits ---------------------
# One row per case and parameter. Empirical arms: testosterone with model-design weights on
# chromosomes 19-22, and JCV/MCV VP1 with the 44-antigen scan weights on chromosomes 6 and 22.
arm_label <- c(corrected_weights = "testosterone model-design weights, chromosomes 19-22",
               scan_weights = "antibody-scan weights, chromosomes 6 and 22")
stopifnot(all(lowsep$emp_arm %in% names(arm_label)))
lowsep_rows <- lowsep[, .(study = "low-separation SE check",
  quantity = sprintf("pooled / jackknife SE ratio in the empirical fit, median over %s variants (%s); separation %.2f SD",
                     format(emp_n_snps[1], big.mark = ","), arm_label[[emp_arm[1]]], separation_fitted[1]),
  scenario = case[1], n = n_samples[1], parameter = parameter[1], covariance_source = "pooled vs jackknife",
  se_agreement = sprintf("%.4f", emp_pooled_over_jk_median[1]),
  coverage = sprintf("null sim: pooled %.3f-%.3f, jackknife %.3f-%.3f",
                     min(sim_pooled_coverage95[sim_effect == "null"]), max(sim_pooled_coverage95[sim_effect == "null"]),
                     min(sim_jk_coverage95[sim_effect == "null"]), max(sim_jk_coverage95[sim_effect == "null"]))),
  by = .(grp_case = case, grp_param = parameter)][, !c("grp_case", "grp_param")]
stopifnot(nrow(lowsep_rows) == uniqueN(lowsep[, .(case, parameter)]))

## 4. Combine and save ----------------------------------------------------------------------
summary_rows[, coverage := as.character(coverage)]; summary_rows[, se_agreement := as.character(se_agreement)]
out <- rbind(summary_rows, timing_rows, lowsep_rows, fill = TRUE)
out[, source := NULL]                                  # file paths of the simulation outputs
save_table(
  as.data.frame(out), "tableS2_simulation_results",
  caption = paste(
    "Table S2. Numerical results of the simulation studies. One row per study, scenario, sample size, minor allele",
    "frequency, parameter, method and covariance (design in Supplementary Table S1; Supplementary Note 1). Bias, RMSE,",
    "empirical s.d. (of the estimates over replicates) and mean SE are relative to the study's estimand; for study 3 the",
    "empirical s.d. is taken across variants with different true effects and is not an error s.d. Coverage is of nominal",
    "95% confidence intervals, with its Monte Carlo standard error. Attempted, succeeded and failed count replicates",
    "(variants in study 3, coefficients in study 1). SE agreement is the ratio of G-MoQ to ordinary RIF–OLS standard",
    "errors (study 1) or, for study 7, the median ratio of pooled to jackknife standard errors in the empirical fit,",
    "taken over the variants of the stated chromosomes (testosterone with model-design weights; JC and Merkel cell",
    "polyomavirus VP1 with the weights of the 44-antigen scan); the coverage column of study 7 gives the range of null",
    "coverage for each covariance over that case's simulated cells (200 replicates each). Timing columns give the time to",
    "fit with standard errors and the whole-process wall time (s) and peak resident memory (MiB), as means of the stated",
    "repeats; mixture timing rows give the median over three 100-variant panels for G-MoQ with pooled covariance against",
    "direct mixture maximum likelihood.",
    "MAF, minor allele frequency; MiB, mebibyte; ML, maximum likelihood; RIF–OLS, recentred influence function ordinary",
    "least squares; RMSE, root mean square error; s.d., standard deviation; SE, standard error."),
  results_used = "08_simulations/summary")
