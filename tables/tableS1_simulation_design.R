#!/usr/bin/env Rscript
# Supplementary Table S1: design of the simulation studies
#
# What: one row per simulation study, described in the ADEMP framing (aim, data-generating
#   mechanism, estimand, methods compared, performance measures; Morris, White & Crowther
#   2019). Studies 1-6 are Supplementary Figs S1-S6, study 7 the pooled against jackknife
#   check at the empirical mixture separations (Fig. S4c) and study 8 the mixture
#   weight-design simulation (Fig. S7b). The design text is written here; sample sizes,
#   replicate counts and grids are read from the simulation results and checked against it.
# Paper: Supplementary Table S1; Supplementary Note 1
# Input: results/08_simulations/summary/{table_s1_design_provenance,table_s2_numerical_summary,
#   lowsep_sim_summary}.csv, results/02_mixtures/weight_design/sim/{operating_characteristics,
#   fit_summary}.csv
# Output: paper_assets/tableS1_simulation_design.csv, _caption.md
# Runs: locally from results/ in about 1 second

suppressPackageStartupMessages({library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load ----------------------------------------------------------------------------------
design <- as.data.table(read_result("08_simulations/summary", "table_s1_design_provenance"))
sim <- as.data.table(read_result("08_simulations/summary", "table_s2_numerical_summary"))
lowsep <- as.data.table(read_result("08_simulations/summary", "lowsep_sim_summary"))
oc <- as.data.table(read_result("02_mixtures/weight_design/sim", "operating_characteristics"))
fs <- as.data.table(read_result("02_mixtures/weight_design/sim", "fit_summary"))

## 2. Numbers read from the results, and checks against the fixed text ----------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
nlist <- function(x) paste(fmt(sort(unique(as.numeric(x)))), collapse = ", ")
attempted <- function(st) { a <- unique(sim[study == st & !is.na(attempted), attempted]); stopifnot(length(a) == 1); a }
stopifnot(identical(design$study, paste0("S", 1:6)),
          attempted("S1") == 45000, attempted("S3") == 1000, attempted("S4") == 20, attempted("S6") == 200,
          setequal(unique(sim[study == "S6", n]), c(1e4, 5e4, 4e5)),
          setequal(unique(sim[study == "S6", maf]), c(0.05, 0.01, 0.001)),
          setequal(unique(sim[study == "S4", n]), c(1e4, 5e4, 4e5)))
stopifnot(all(lowsep$replicates_failed == 0), uniqueN(lowsep$replicates_ok) == 1, uniqueN(lowsep$cell) == 40)
mix_scenarios <- unique(oc[role == "gamma" & param == "gamma", scenario])
n_null <- uniqueN(oc[role == "null", snp])
stopifnot(length(mix_scenarios) == 9, uniqueN(oc$n_reps) == 1)

## 3. Table ---------------------------------------------------------------------------------
out <- data.table(
  Study = paste("Study", 1:8),
  `Supplementary figure` = c("S1", "S2", "S3", "S4a,b", "S4d, S5", "S6", "S4c", "S7b"),
  Aim = c(
    "Check that Stage 1 reproduces separately fitted RIF–OLS regressions",
    "Compare the computational cost of Stage 1 with separately fitted RIF–OLS regressions",
    "Check that Stage 2 recovers location, dispersion and mixture-component effects without bias and with calibrated standard errors",
    "Compare Stage 2 mixture estimates with direct mixture maximum likelihood",
    "Compare the cost of the two mixture methods",
    "Compare pooled and 1,000-block jackknife standard errors for rare variants",
    "Compare pooled and jackknife standard errors at the separations of the empirical mixtures",
    "Compare the model and hybrid weight designs for the mixture membership parameter"),
  `Data-generating mechanism` = c(
    sprintf("Quantitative phenotype; N = %s; %s variants", fmt(4e5), fmt(1000)),
    sprintf("Matched workloads: N = %s with 1 or 256 variants, N = %s with 1 or %s variants", fmt(2000), fmt(4e5), fmt(1000)),
    sprintf("Normal phenotype with an effect on location, dispersion or both; two-component normal mixture with an effect on the lower-component location; N = %s; minor allele frequency 0.3; per-allele effects drawn from U(−0.1, 0.1)", fmt(4e5)),
    sprintf("Two-component mixture with effects on the lower-component location drawn from U(−0.1, 0.1); N = %s; minor allele frequency 0.3", nlist(sim[study == "S4", n])),
    sprintf("The panels of study 4; N = %s", nlist(sim[study == "S4", n])),
    sprintf("Normal phenotype with no effect or a dispersion effect (c = 0.1); N = %s; minor allele frequency %s",
            nlist(sim[study == "S6", n]), paste(sort(unique(sim[study == "S6", maf]), decreasing = TRUE), collapse = ", ")),
    sprintf("Two-component normal mixtures (component s.d. 1) with the separations (%.2f–%.2f s.d.) and upper-component proportions (%.2f–%.2f) of the testosterone and antibody mixtures; N = %s (%s for the testosterone geometries only); minor allele frequency %s; no effect or an upper-component shift sized for |z| ≈ 5",
            min(lowsep$separation), max(lowsep$separation), min(lowsep$pi_upper), max(lowsep$pi_upper), fmt(1e4), fmt(4e5),
            paste(sort(unique(lowsep$maf), decreasing = TRUE), collapse = " or ")),
    sprintf("Two-component skew-normal mixtures in %d scenarios, from well separated (with and without misfit of the density in the gap between components) to heavily overlapping; N = %s; minor allele frequency 0.3; in each replicate, variants with effects on each component location, on membership, on all three and on both locations, and %d null variants",
            length(mix_scenarios), nlist(fs[scenario %in% mix_scenarios, n]), n_null)),
  Estimand = c(
    "Coefficients and standard errors of ordinary RIF–OLS at each of 45 quantile levels",
    "None (cost only)",
    "Population projection target of the generating effect onto the model's weights",
    "Each method's own known target",
    "None (cost only)",
    "Population projection targets of location and dispersion effects",
    "Population projection targets of membership and component locations",
    "Generating effects: lower-component location 0.055, upper-component location 0.04, membership 0.08"),
  `Methods compared` = c(
    "Stage 1 (pooled covariance; 1,000-block jackknife) against ordinary RIF–OLS",
    "Ordinary RIF–OLS; Stage 1 with pooled or jackknife covariance",
    "Stage 2 by ordinary least squares with pooled or jackknife covariance",
    "Stage 1 and mixture Stage 2 (pooled or jackknife covariance) against mixture maximum likelihood refitted at every variant",
    "As study 4",
    "Pooled against jackknife covariance",
    "Pooled against jackknife covariance; also compared in the empirical fits",
    "Model weights (every column divided by the fitted mixture density) against hybrid weights (membership column divided by the Stage 1 kernel density)"),
  `Performance measures` = c(
    "Coefficient difference (mean and root mean square); ratio of standard errors",
    "Time to fit with standard errors; whole-process time; peak memory (means of 3 repeats, 1 thread)",
    "Bias over RMSE; mean SE over RMSE; 95% coverage; Monte Carlo intervals",
    "Bias; mean SE over RMSE; Monte Carlo intervals",
    "Time to fit; whole-process time; peak memory (3 panels per N); implied genome-wide cost",
    "Bias and mean SE over the empirical s.d.; 95% coverage; Monte Carlo intervals",
    "Pooled over jackknife SE; SE over empirical s.d.; 95% coverage with Monte Carlo SE",
    "Mean estimate over truth; leakage into the lower-component location; 95% coverage; s.d. of z for null variants"),
  Replicates = c(
    sprintf("1 dataset (%s coefficients)", fmt(attempted("S1"))),
    "3 timing repeats per workload",
    sprintf("%s variants per scenario", fmt(attempted("S3"))),
    sprintf("%d panels of 100 variants per N", attempted("S4")),
    "3 panels per N",
    sprintf("%d per cell", attempted("S6")),
    sprintf("%d per cell (%d cells)", lowsep$replicates_ok[1], uniqueN(lowsep$cell)),
    sprintf("%d per scenario", oc$n_reps[1])))

## 4. Save ----------------------------------------------------------------------------------
save_table(
  as.data.frame(out), "tableS1_simulation_design",
  caption = paste(
    "Table S1. Design of the simulation studies. One row per study, described by its aim, data-generating mechanism,",
    "estimand, the methods compared and the performance measures (Supplementary Note 1). Every Stage 1 analysis used the",
    "45 quantile levels 0.02–0.98 of the empirical analyses and, for the jackknife, 1,000 blocks with individuals",
    "assigned at random; studies 1–6 used 1 CPU thread for timings. Numerical results are in Supplementary Table S2.",
    "N, simulated sample size; RIF–OLS, recentred influence function ordinary least squares; RMSE, root mean square",
    "error; s.d., standard deviation; SE, standard error; U(a, b), uniform distribution on [a, b]."),
  results_used = c("08_simulations/summary", "02_mixtures/weight_design/sim"))
