#!/usr/bin/env Rscript
# Supplementary Table S3: testosterone mixture components against sex-stratified reference GWAS
#
# What: for each testosterone phenotype (unshrunk and three shrinkage settings) and each
#   Stage 2 weight design, compares the female, male and membership components with their
#   references: genomic inflation and LD score regression intercepts, correlation and slope
#   with the sex-stratified REGENIE effects, LDSC genetic correlation and internal Mendelian
#   randomisation in both directions. Membership (gamma) is compared only with a GWAS of
#   chromosomal sex in the same participants; its other cells are blank.
# Paper: Supplementary Table S3
# Input: results/03_testosterone/recovery/{table_s4_inputs,lambda_gc}.csv
# Output: paper_assets/tableS3_testosterone_components.csv, _caption.md
# Runs: locally from results/ in about 1 second

suppressPackageStartupMessages({library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

cfg <- read_config()

## 1. Load and join -------------------------------------------------------------------------
d <- as.data.table(read_result("03_testosterone/recovery", "table_s4_inputs"))
lambda <- as.data.table(read_result("03_testosterone/recovery", "lambda_gc"))
lambda <- lambda[reference %in% c("combined", "sex"), .(weight_design, arm, param, lambda_gc_regenie_all,
                                                        ldsc_gmoq_intercept, ldsc_gmoq_intercept_se, ldsc_gmoq_ratio,
                                                        ldsc_ref_intercept, ldsc_ref_intercept_se)]
d <- merge(d, lambda, by = c("weight_design", "arm", "param"), all.x = TRUE)
# Genetic correlation of membership with the GWAS of chromosomal sex in all 462,869 genotyped
# participants (the sensitivity reference quoted in the text), next to the Stage 1-sample one.
rg_full <- as.data.table(read_result("03_testosterone/sex_gwas", "ldsc_rg_gamma_vs_sex"))[sex_gwas == "old_n462869",
  .(weight_design, arm, param, rg_full = rg, rg_full_se = rg_se)]
d <- merge(d, rg_full, by = c("weight_design", "arm", "param"), all.x = TRUE)
# Reported design per arm (config.yaml): hybrid for unshrunk and F = 0.2, model otherwise.
d[, reported := weight_design == mixture_design(cfg, arm)]

## 2. Format --------------------------------------------------------------------------------
arm_labels <- c(testosterone = "unshrunk", testosterone_shrink_0_2 = "shrink 0.2",
                testosterone_shrink_0_5 = "shrink 0.5", testosterone_shrink_0_7 = "shrink 0.7")
comp_labels <- c(beta_left = "female (beta_left)", beta_right = "male (beta_right)", gamma = "sex/mixture (gamma)")
rd <- function(x, digits = 3) ifelse(is.na(x), NA_character_, formatC(round(x, digits), format = "f", digits = digits))
pair <- function(b, se, digits = 3) ifelse(is.na(b) | is.na(se), NA_character_, paste0(rd(b, digits), " (", rd(se, digits), ")"))
pair_n <- function(b, se, n, digits = 3) ifelse(is.na(b) | is.na(se), NA_character_,
                                                paste0(rd(b, digits), " (", rd(se, digits), "), n=", n))
out <- d[, .(
  reported_design = ifelse(reported, "reported (main text)", "alternative (S7)"),
  weight_design = weight_design,
  arm = factor(arm_labels[arm], levels = unname(arm_labels)),
  component = factor(comp_labels[param], levels = unname(comp_labels)),
  n_gws = n_gws,
  lambda_GC = rd(lambda_gc),
  lambda_GC_regenie = rd(lambda_gc_regenie_all),
  ldsc_intercept_gmoq = pair(ldsc_gmoq_intercept, ldsc_gmoq_intercept_se),
  ldsc_intercept_regenie = pair(ldsc_ref_intercept, ldsc_ref_intercept_se),
  r_all_vs_reference = rd(r_all_combined),
  slope_all_vs_reference = rd(slope_all_combined),
  r_gws_vs_reference = rd(r_gws_combined),
  slope_gws_vs_reference = rd(slope_gws_combined),
  r_all_paper_matched = rd(r_all_paper_matched),
  r_all_sex_gwas = rd(r_all_sex),
  ldsc_rg = pair(ldsc_rg, ldsc_rg_se),
  ldsc_rg_sex_gwas_full_cohort = pair(rg_full, rg_full_se),
  mr_regenie_to_gmoq = pair_n(mr_regenie_to_gmoq_beta, mr_regenie_to_gmoq_se, mr_regenie_to_gmoq_nsnp),
  mr_gmoq_to_regenie = pair_n(mr_gmoq_to_regenie_beta, mr_gmoq_to_regenie_se, mr_gmoq_to_regenie_nsnp)
)]
# Membership is weakly identified and its sign is a convention: it is compared only with the
# chromosomal-sex GWAS (r_all_sex_gwas, ldsc_rg); every column referring to the sex-stratified
# references, or to a slope, is blank for it. Conversely r_all_sex_gwas is blank for the
# female and male components.
is_gamma <- out$component == "sex/mixture (gamma)"
blank_for_gamma <- c("n_gws", "lambda_GC", "lambda_GC_regenie", "r_all_vs_reference", "slope_all_vs_reference",
                     "r_gws_vs_reference", "slope_gws_vs_reference", "r_all_paper_matched",
                     "mr_regenie_to_gmoq", "mr_gmoq_to_regenie")
for (cc in blank_for_gamma) out[[cc]][is_gamma] <- NA
out$r_all_sex_gwas[!is_gamma] <- NA
out[, reported_design := factor(reported_design, levels = c("reported (main text)", "alternative (S7)"))]
setorder(out, reported_design, arm, component)

## 3. Save ----------------------------------------------------------------------------------
save_table(
  as.data.frame(out), "tableS3_testosterone_components",
  caption = paste(
    "Table S3. Testosterone mixture components against sex-stratified reference GWAS, by phenotype and mixture weight",
    "design. One row per phenotype, weight design and mixture parameter; n = 399,710 (199,572 female, 200,138 male).",
    "Stage 2 by ordinary least squares with the pooled Stage 1 covariance. Weights: ‘Model’, every column divided by",
    "the fitted mixture density; ‘Hybrid’, only the membership (γ) column divided by the Stage 1 kernel density",
    "(Supplementary Note 2). Design: ‘Reported’, the design used in Figs 1 and 2 and Supplementary Fig. S8 (hybrid",
    "for unshrunk and F = 0.2, model for F = 0.5 and 0.7); ‘Alternative’, the other design. Phenotype: unshrunk, or male values",
    "shifted towards the female mean by F = 0.2, 0.5 or 0.7 of the difference. Components: female (lower component",
    "location), male (upper component location) and membership (γ). GWS variants: two-sided P < 5 × 10⁻⁸. λGC:",
    "median genomic inflation of the component and of the reference. r and slope: Pearson correlation and",
    "least-squares slope against the reference over all variants and over reference GWS variants; LOCO-matched",
    "reference: fitted in the same participants with the same covariates and leave-one-chromosome-out predictions.",
    "IVW: inverse-variance weighted Mendelian randomization slope (s.e., number of instruments) in each direction.",
    sprintf("Membership is compared only with a GWAS of chromosomal sex in the same %s participants (sign aligned to",
            format(cfg$p08$sex_reference_n, big.mark = ",")),
    "male = 1) and, in the last LDSC column, with a GWAS of chromosomal sex in all 462,869 genotyped participants;",
    "its other cells are blank. GWAS, genome-wide association study; GWS, genome-wide significant; IVW,",
    "inverse-variance weighted; LDSC, LD score regression; LOCO, leave-one-chromosome-out; rg, genetic correlation;",
    "s.e., standard error."),
  results_used = c("03_testosterone/recovery", "03_testosterone/sex_gwas"))
