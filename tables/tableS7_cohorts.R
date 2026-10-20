#!/usr/bin/env Rscript
# Supplementary Table S7: cohorts and variant counts for every trait with a Stage 1 quantile GWAS
#
# What: one row per trait and scale: the participant counts at each step from the analysis
#   sample to the number analysed, the covariates, the number of variants analysed, the
#   number of lead variants, and where in the paper the trait is reported.
# Paper: Supplementary Table S7
# Input: results/01_data/sample_sizes/cohort_funnel.csv,
#   results/01_data/stage1_gwas_index/stage1_traits.csv,
#   results/01_data/lead_variants/lead_sets.csv, results/06_traits/survey/table1.csv
# Output: paper_assets/tableS7_cohorts.csv, _caption.md
# Runs: locally from results/ in about 1 second

suppressPackageStartupMessages({library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load ----------------------------------------------------------------------------------
cohort <- as.data.table(read_result("01_data/sample_sizes", "cohort_funnel"))
traits <- as.data.table(read_result("01_data/stage1_gwas_index", "stage1_traits"))
leads <- as.data.table(read_result("01_data/lead_variants", "lead_sets"))[primary == TRUE, .(trait, n_leads)]
# The Fig. 5 traits: the survey traits less Lp(a), which is censored at both assay limits.
excluded <- unlist(read_config()$estimation$het_pc$exclude_traits)
fig5_traits <- setdiff(as.data.table(read_result("06_traits/survey", "table1"))[in_fig6 == TRUE, trait], excluded)
stopifnot(length(fig5_traits) == 21)

## 2. Build ---------------------------------------------------------------------------------
out <- cohort[, .(trait, group, N_keep = as.integer(n_keep), valid_phenotype = as.integer(n_valid_pheno),
                  after_LOCO = as.integer(n_after_loco), analysed = as.integer(n_analysed))]
out[traits, on = "trait", `:=`(n_covar = i.n_covariates, covar_names = i.covariates, n_variants = i.n_variants)]
out[, covariates := paste0(n_covar, ": ", covar_names)]          # e.g. "13: sex, age, ..."
out[leads, on = "trait", n_leads := i.n_leads]
out[, paper_use := fcase(
  trait %in% excluded, "Excluded (censored assay)",
  grepl("^testosterone", trait), "Figs 1, 2",
  trait %in% c("jcv_vp1_log_all", "mcv_vp1_log_all"), "Fig. 3",
  trait == "log_bmi", "Figs 4–6",
  trait == "bmi", "Supplementary Figs S12, S14; Supplementary Table S6",
  trait %in% fig5_traits, "Fig. 5",
  default = "Not reported")]
stopifnot(sum(out$paper_use == "Excluded (censored assay)") == 1)
out <- out[, .(trait, group, N_keep, valid_phenotype, after_LOCO, analysed, covariates,
               n_variants = as.integer(n_variants), n_leads_primary = as.integer(n_leads), paper_use)]
setorder(out, group, trait)

## 3. Save ----------------------------------------------------------------------------------
save_table(
  as.data.frame(out), "tableS7_cohorts",
  caption = paste(
    "Table S7. Cohorts and variant counts for all traits with a Stage 1 quantile GWAS. One row per trait and scale",
    "(log-transformed traits are marked ‘(log)’). N (sample): participants in the analysis sample with the trait;",
    "N phenotyped: with a valid phenotype; N after LOCO: with leave-one-chromosome-out predictions, where used; N analysed:",
    "the final number analysed. Covariates: the covariates, with their number in brackets; weight is the height-adjusted",
    "model (sex, age, four height terms and genetic principal components 1–10). Variants: HapMap3 variants analysed",
    "(minor allele frequency above 0.05 additionally required for the antibodies). Leads: lead variants (REGENIE where",
    "available, otherwise the Stage 1 mean effect; Methods); antibody leads are those of a uniform REGENIE clump (JCV VP1",
    "3, MCV VP1 4), and the testosterone analyses use sex-stratified lead sets (Methods). Use: where the trait is",
    "reported; lipoprotein(a) was analysed but is excluded from the paper because its measurements are censored at both",
    "assay limits, and traits marked ‘Not reported’ were analysed but are not used. JCV, JC polyomavirus; LOCO,",
    "leave-one-chromosome-out; MCV, Merkel cell polyomavirus; PC, genetic principal component; VP1, major capsid protein",
    "VP1."),
  results_used = c("01_data/sample_sizes", "01_data/stage1_gwas_index", "01_data/lead_variants", "06_traits/survey"))
