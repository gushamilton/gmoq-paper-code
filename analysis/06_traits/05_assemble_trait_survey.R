#!/usr/bin/env Rscript
# Trait survey: one summary row per trait, and the numbers quoted in the text
#
# What: Collects, for each trait, the per-trait results of 01_moments_and_principal_components.R
#   and the Box-Cox powers of 03_normalising_powers.R into one table: sample size, leads,
#   components above noise, the common ratio r of spread to location and the between-lead
#   heterogeneity tau_h (95% bootstrap CI, P, detection floor), spread and skew-beyond counts and
#   lack of fit of the shift + spread frame. No new statistics. A trait "fails the frame" when
#   more than 20% of its leads have lack-of-fit P < 0.05; an "exception" has more than two
#   components, fails the frame, or has skew beyond the frame at more than 5% of its leads.
#   The HbA1c rows with and without diabetes are set side by side.
# Paper: Fig. 5a,b; Supplementary Table S5; Supplementary Table S7 (trait list); Results (trait
#   survey and HbA1c counts)
# Input: results/06_traits/moments/*.csv, results/06_traits/scale/{boxcox_lambda,
#   scale_comparison}.csv, results/01_data/lead_variants/lead_sets.csv
# Output: results/06_traits/survey/{table1,fig6_components,fig6_exceptions,survey_statements}.csv
# Runs: locally from results/ in a few seconds

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

moments <- function(name) as.data.table(read_result("06_traits/moments", name))
arms <- moments("arms"); het <- moments("heterogeneity"); floors <- moments("detection_floor")
skew <- moments("skew_counts"); lof <- moments("lof_summary"); spectrum <- moments("pca_spectrum")
powers <- as.data.table(read_result("06_traits/scale", "boxcox_lambda"))
scales <- as.data.table(read_result("06_traits/scale", "scale_comparison"))
lead_sets <- as.data.table(read_result("01_data/lead_variants", "lead_sets"))

## 1. One row per trait
# Heterogeneity is reported on the unit-SE shift/spread basis with excess Q - (n - 1).
het_main <- het[basis == "unit_se" & component == "c2_on_c1" & excess == "n_minus_1"]
het_exact <- het[basis == "exact" & excess == "n_minus_1"]
T1 <- merge(arms[group %in% c("bmi", "weight", "survey", "survey_extra"),
                 .(arm, trait, group, scale, lead_set, n_sample, n_leads, n_components)],
            het_main[, .(arm, common_tilt_r = r, common_tilt_r_lo95 = r_lo95, common_tilt_r_hi95 = r_hi95,
                         tau, tau_lo95, tau_hi95, p_het = p, Q_het = Q, df_het = df, n_resid_z_gt_3)], by = "arm")
T1 <- merge(T1, floors[, .(arm, tau_floor_80, observed_above_floor)], by = "arm", all.x = TRUE)
T1 <- merge(T1, skew[, .(arm, n_spread_lambda2_bonf, n_skew_raw_bonf, n_skew_beyond_bonf, n_skew_beyond_nominal,
                         expected_nominal, median_spread_lambda2_z)], by = "arm")
T1 <- merge(T1, lof[frame == "shift_spread", .(arm, lof_frac_p_lt_0.05 = frac_p_lt_0.05,
                                               lof_n_bonf_rejected = n_bonf_rejected)], by = "arm")
T1 <- merge(T1, het_exact[, .(arm, tau_over_r_exact = tau_over_r)], by = "arm")
T1 <- merge(T1, spectrum[pc <= 2, .(pc1_pct = pct[1], pc2_pct = pct[2]), by = arm], by = "arm")
# Box-Cox powers belong to the raw trait: log BMI takes those of BMI.
T1[, raw_trait := fifelse(trait == "log_bmi", "bmi", trait)]
T1 <- merge(T1, powers[, .(raw_trait = trait, lambda_boxcox_main = lambda_main, lambda_boxcox_ml = lambda_ml)],
            by = "raw_trait", all.x = TRUE)
T1[, lead_rule := lead_sets$rule[match(lead_set, lead_sets$lead_set)]]
T1[, frame_failure := lof_frac_p_lt_0.05 > 0.20]
T1[, skew_beyond_share := n_skew_beyond_bonf / n_leads]
T1[, at_most_two_components := n_components <= 2]
T1[, exception := n_components > 2 | frame_failure | skew_beyond_share > 0.05]
T1[, in_fig6 := group %in% c("bmi", "weight", "survey") & arm != "bmi"]   # BMI appears as log BMI
setorder(T1, -in_fig6, group, trait)

## 2. Per-trait components and exceptions (Fig. 5)
F6 <- T1[in_fig6 == TRUE, .(trait = arm, scale, n_leads, n_components, n_skew_beyond_bonf, skew_beyond_share,
                            lof_frac_p_lt_0.05, frame_failure, exception)]
hba1c <- T1[arm %in% c("log_hba1c", "hba1c_nodiab", "hba1c_nodiab_log"),
            .(item = paste0("hba1c:", arm), trait, scale, n_leads, n_components, n_skew_beyond_bonf, lof_frac_p_lt_0.05,
              note = fifelse(arm == "log_hba1c", "with diabetes (log)", fifelse(arm == "hba1c_nodiab_log",
                "diabetes excluded, converted to log (like-for-like)", "diabetes excluded, analysed raw")))]
lpa <- T1[arm == "log_lipoprotein_a", .(item = "lpa:frame_failure", trait, scale, n_leads, n_components,
          n_skew_beyond_bonf, lof_frac_p_lt_0.05, note = "shift + spread frame rejected (lack of fit)")]
other <- T1[in_fig6 == TRUE & exception == TRUE & !arm %in% c("log_hba1c", "log_lipoprotein_a"),
            .(item = paste0("other:", arm), trait, scale, n_leads, n_components, n_skew_beyond_bonf, lof_frac_p_lt_0.05,
              note = paste0(fifelse(n_components > 2, "components > 2; ", ""),
                            fifelse(frame_failure, "frame lack of fit > 20%; ", ""),
                            fifelse(skew_beyond_share > 0.05, "skew beyond > 5% of leads", "")))]
exceptions <- rbind(hba1c, lpa, other, fill = TRUE)

## 3. Numbers quoted in the text
f6 <- T1[in_fig6 == TRUE]
statements <- data.table(
  quantity = c("n_traits_fig6", "n_traits_at_most_two_components", "n_traits_one_component",
               "n_traits_two_components", "n_traits_more_than_two", "n_traits_frame_failure",
               "hba1c_skew_beyond_with_diabetes_log", "hba1c_skew_beyond_without_diabetes_log",
               "hba1c_skew_beyond_without_diabetes_raw", "log_bmi_tau", "log_bmi_tau_floor",
               "bmi_spread_leads_raw", "bmi_spread_leads_log", "bmi_spread_leads_boxcox"),
  value = c(nrow(f6), sum(f6$n_components <= 2), sum(f6$n_components == 1), sum(f6$n_components == 2),
            sum(f6$n_components > 2), sum(f6$frame_failure),
            T1[arm == "log_hba1c"]$n_skew_beyond_bonf, T1[arm == "hba1c_nodiab_log"]$n_skew_beyond_bonf,
            T1[arm == "hba1c_nodiab"]$n_skew_beyond_bonf, T1[arm == "log_bmi"]$tau, T1[arm == "log_bmi"]$tau_floor_80,
            scales[trait == "bmi" & scale == "raw"]$n_spread_lambda2_bonf,
            scales[trait == "bmi" & scale == "log"]$n_spread_lambda2_bonf,
            scales[trait == "bmi" & scale == "boxcox"]$n_spread_lambda2_bonf))
print(statements)

## 4. Write
out <- "06_traits/survey"
write_result(T1, out, "table1")
write_result(F6, out, "fig6_components")
write_result(exceptions, out, "fig6_exceptions")
write_result(statements, out, "survey_statements")
