#!/usr/bin/env Rscript
# Main Table 1: biological correlates of the primary component results.
# Reads aggregate results only. Primary PC2 ratios come from the same one-shape
# null as Figure 6b. Modifier correlations use unnormalised whitened projections;
# HbA1c correlations use unit-length whitened profiles. Subgroup component counts
# from the historical auxiliary procedure are not used for dimension comparisons.
suppressPackageStartupMessages(library(data.table))
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

mom <- as.data.table(read_result("06_traits/moments", "lead_moments"))
pa <- as.data.table(read_result("06_traits/moments", "pca_parallel"))[component_tested == 2]
traits <- c("log_hba1c", "hba1c_nodiab", "log_alt", "log_urate", "fev1_fvc_ratio", "log_crp", "ldl", "log_bilirubin_total")
labels <- c("Log HbA1c", "HbA1c without diabetes (raw)", "Log ALT", "Log urate", "FEV1/FVC", "Log CRP", "LDL", "Log bilirubin")
lead_n <- mom[arm %in% traits, .(n_leads = .N), by = .(trait = arm)]
ratio <- pa[match(traits, arm), .(trait = arm, dimension = n_components,
                                pc2_ratio = observed_pct / null_median_pct)]
stopifnot(nrow(ratio) == 8L, !anyNA(ratio), !anyDuplicated(ratio$trait))
base <- merge(ratio, lead_n, by = "trait", sort = FALSE)[match(traits, trait)]
base[, trait_label := labels]

# One row per correlation, with its actual lead denominator and score definition.
cor_rows <- list()
add_cor <- function(trait, endpoint, x, y, score) {
  ok <- is.finite(x) & is.finite(y)
  test <- cor.test(x[ok], y[ok], method = "pearson", alternative = "two.sided")
  cor_rows[[length(cor_rows) + 1L]] <<- data.table(trait, endpoint, score,
    n_leads = sum(ok), pearson_r = unname(test$estimate), nominal_p = test$p.value)
}
for (s in list(
  c("log_alt", "alt_masld_per_lead", "bmi_int_z", "SNP x BMI interaction"),
  c("log_alt", "alt_masld_per_lead", "masld_any_z", "Recorded liver disease"),
  c("log_urate", "urate_sex_per_lead", "sex_int_z", "SNP x sex interaction"),
  c("fev1_fvc_ratio", "fev_smoking_per_lead", "smok_int_z", "SNP x ever-smoking interaction"),
  c("log_crp", "crp_acute_per_lead", "bmi_int_z", "SNP x BMI interaction"),
  c("log_bilirubin_total", "bilirubin_dom_per_lead", "dom_dev_z", "Quadratic dosage term adjusted for additive dosage"))) {
  x <- merge(mom[arm == s[1], .(snp_id, pc2 = whitened_norm * pc2_score)],
             as.data.table(read_result("06_traits/hidden_biology", s[2])), by = "snp_id")
  stopifnot(nrow(x) == base[trait == s[1], n_leads])
  add_cor(s[1], s[4], x$pc2, x[[s[3]]], "unnormalised whitened PC2 projection")
}
# Align the diabetes association to the effect allele of each profile set.
diabetes <- as.data.table(read_result("07_hba1c/skew_diabetes", "per_lead"))[
  , .(snp_id, diab_z, diab_effect_allele = effect_allele, diab_other_allele = other_allele)]
for (a in c("log_hba1c", "hba1c_nodiab")) {
  x <- merge(mom[arm == a], diabetes, by = "snp_id")
  same <- x$effect_allele == x$diab_effect_allele & x$other_allele == x$diab_other_allele
  swap <- x$effect_allele == x$diab_other_allele & x$other_allele == x$diab_effect_allele
  stopifnot(all(same | swap))
  z <- x$diab_z * ifelse(same, 1, -1)
  pc2 <- x$pc2_score
  if (cor(pc2, z) < 0) pc2 <- -pc2  # PC sign convention in Figure 5d.
  add_cor(a, "Diabetes", pc2, z, "unit-length whitened PC2 direction")
  if (a == "hba1c_nodiab")
    add_cor(a, "Diabetes", x$het_resid_z, z, "heterogeneity residual z")
}
correlations <- rbindlist(cor_rows)

# The two asymmetry denominators use the independently selected log-scale panels.
skew <- mom[arm %in% c("log_hba1c", "hba1c_nodiab_log"),
            .(n_leads = .N, n_asymmetry = sum(2 * pnorm(-abs(skew_beyond_frame_z)) < 0.05 / .N)),
            by = .(trait = arm)]
ldl <- as.data.table(read_result("06_traits/hidden_biology", "ldl_statin_summary"))
stopifnot(identical(skew$n_asymmetry[match(c("log_hba1c", "hba1c_nodiab_log"), skew$trait)], c(69L, 4L)),
          all(ldl$n_leads == 230), identical(ldl$n_skew_bonf, c(6L, 0L)))

# Numeric source rows retain full precision and identify the analysis behind each value.
source <- rbindlist(list(
  base[, .(trait, quantity = "primary_pc2_ratio", value = pc2_ratio, n_leads)],
  base[, .(trait, quantity = "primary_dimension", value = as.numeric(dimension), n_leads)],
  correlations[, .(trait, quantity = paste(score, "vs", endpoint, "Pearson r"), value = pearson_r, n_leads)],
  correlations[, .(trait, quantity = paste(score, "vs", endpoint, "nominal two-sided P"), value = nominal_p, n_leads)],
  skew[, .(trait, quantity = "log-scale asymmetry hits", value = as.numeric(n_asymmetry), n_leads)],
  ldl[, .(trait = arm, quantity = "asymmetry hits", value = as.numeric(n_skew_bonf), n_leads)]
), fill = TRUE)

r_of <- function(a, e, score_type = "unnormalised whitened PC2 projection") {
  correlations[trait == a & endpoint == e & score == score_type, pearson_r]
}
base[, proposed_interpretation := c(
  "Diabetes-enriched mixture", "Glycaemic and red-cell turnover pathways",
  "Recorded liver disease and SNP x BMI interaction", "SNP x sex interaction despite sex adjustment",
  "SNP x ever-smoking interaction", "SNP x BMI interaction",
  "Profile differences by reported lipid medication", "Non-additive genotype-dosage effects")]
base[, key_loci := c("TCF7L2, HLA-DOB, CDKAL1", "TCF7L2, GCK; HK1, HFE, ANK1",
  "PNPLA3, TM6SF2, HSD17B13, GPAM, MBOAT7", "SLC2A9, ABCG2", "CHRNA3/CHRNA5 (rs12914385)",
  "LEPR, IL1RN, CRP", "APOE, PCSK9, LDLR; APOB, ABCG8", "UGT1A1, SLCO1B1")]
base[, quantitative_evidence := c(
  sprintf("PC2 vs diabetes r = %.2f (472 leads); log-scale asymmetry: 69/472 with diabetes included and 4/533 after exclusion",
          r_of("log_hba1c", "Diabetes", "unit-length whitened PC2 direction")),
  sprintf("PC2 vs diabetes r = %.2f across 266 shared leads (heterogeneity residual r = %.2f)",
          r_of("hba1c_nodiab", "Diabetes", "unit-length whitened PC2 direction"),
          r_of("hba1c_nodiab", "Diabetes", "heterogeneity residual z")),
  sprintf("PC2 vs SNP x BMI r = %.2f; vs recorded liver disease r = %.2f", r_of("log_alt", "SNP x BMI interaction"), r_of("log_alt", "Recorded liver disease")),
  sprintf("PC2 vs SNP x sex r = %.2f", r_of("log_urate", "SNP x sex interaction")),
  sprintf("PC2 vs SNP x ever-smoking r = %.2f; rs12914385 PC2 = %.2f", r_of("fev1_fvc_ratio", "SNP x ever-smoking interaction"),
          mom[arm == "fev1_fvc_ratio" & snp_id == "rs12914385", whitened_norm * pc2_score]),
  sprintf("PC2 vs SNP x BMI r = %.2f across 389 leads", r_of("log_crp", "SNP x BMI interaction")),
  "Asymmetry hits: 6/230 in all participants and 0/230 in participants without reported cholesterol-lowering medication",
  sprintf("PC2 vs quadratic dosage term (g - 1)^2, adjusted for additive dosage, r = %.2f", r_of("log_bilirubin_total", "Quadratic dosage term adjusted for additive dosage")))]
base[, pc2_ratio := round(pc2_ratio, 2)]
cap <- paste("Table 1. Biological correlates of scale-invariant dimensionality in UK Biobank quantitative traits.",
  "Dimensions and PC2 ratios use the primary one-shape null, as in Figure 6b.",
  "Correlations are Pearson r across selected lead variants in the same cohort; two-sided P values in the source data are nominal.",
  "HbA1c correlations use unit-length whitened PC2 directions; other correlations use unnormalised whitened PC2 projections.",
  "The no-diabetes HbA1c component ratio uses raw HbA1c; the separate asymmetry count uses log HbA1c.",
  "The restricted-panel asymmetry counts do not imply a change in primary component dimensionality.")
save_table(as.data.frame(base[, .(trait = trait_label, n_leads, dimension, pc2_ratio,
                                 proposed_interpretation, key_loci, quantitative_evidence)]),
           "table1_biological_correlates", cap,
           results_used = c("06_traits/moments", "06_traits/hidden_biology", "07_hba1c/skew_diabetes"))
fwrite(source, file.path(repo_path("paper_assets"), "table1_biological_correlates_source.csv"))
