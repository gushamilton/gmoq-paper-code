#!/usr/bin/env Rscript
# Supplementary Table S4: lead variants for JCV and MCV VP1 antibody levels
#
# What: one row per lead variant and antigen (plus the Fig. 3 examples that are not leads),
#   with the Stage 2 mixture decomposition under the reported weight design (model for JCV,
#   hybrid for MCV): membership of the upper component, lower and upper locations, joint test,
#   mixture lack of fit; the same under the other design; the location/spread coordinates with their
#   detection bounds, and the lack of fit of the shift + spread frame.
# Paper: Supplementary Table S4; Supplementary Note 3
# Input: results/04_antibodies/decomposition/{lead_decomposition,detection_bounds_per_lead,
#   lack_of_fit_per_lead}.csv
# Output: paper_assets/tableS4_antibody_leads.csv, _caption.md
# Runs: locally from results/ in about 1 second

suppressPackageStartupMessages({library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load and join -------------------------------------------------------------------------
# Reported weight design per antibody (config.yaml): 'model' for JCV VP1, 'hybrid' for MCV VP1.
all_leads <- as.data.table(read_result("04_antibodies/decomposition", "lead_decomposition"))
all_leads[, reported := weight_design == mixture_design(read_config(), trait)]
detection <- as.data.table(read_result("04_antibodies/decomposition", "detection_bounds_per_lead"))
lack_of_fit <- as.data.table(read_result("04_antibodies/decomposition", "lack_of_fit_per_lead"))
alternative <- all_leads[reported == FALSE, .(trait, snp_id, alt_design = weight_design, alt_gamma = gamma,
                                               alt_gamma_se = gamma_se, alt_label = label, alt_hotelling = hotelling_mlog10p,
                                               alt_Q_mixture = Q_mixture)]
leads <- all_leads[reported == TRUE][is_p04_lead == TRUE | is_example == TRUE]
leads <- merge(leads, alternative, by = c("trait", "snp_id"), all.x = TRUE)
m <- merge(leads, detection[, .(trait, snp_id, chr, bp, L, S, z_location, z_spread, bound_nominal, bound_bonferroni)],
           by = c("trait", "snp_id", "chr", "bp"), all.x = TRUE)
m <- merge(m, lack_of_fit[, .(trait, snp_id, chr, bp, Q_ols_location_scale, p_davies_location_scale, davies_floored_location_scale)],
           by = c("trait", "snp_id", "chr", "bp"), all.x = TRUE)

## 2. Format --------------------------------------------------------------------------------
rd <- function(x, digits = 3) ifelse(is.na(x), NA_character_, formatC(round(x, digits), format = "f", digits = digits))
pair <- function(b, se, digits = 3) ifelse(is.na(b) | is.na(se), NA_character_, paste0(rd(b, digits), " (", rd(se, digits), ")"))
pct <- function(x, digits = 2) ifelse(is.na(x), NA_character_, paste0(rd(x * 100, digits), "%"))
zpair <- function(v, z, digits = 3) ifelse(is.na(v) | is.na(z), NA_character_, paste0(rd(v, digits), " (z=", rd(z, 2), ")"))
trait_labels <- c(mcv_vp1_log_all = "MCV VP1 (log)", jcv_vp1_log_all = "JCV VP1 (log)")
out <- m[, .(
  trait = trait_labels[trait],
  snp_id, chr, bp,
  MHC = mhc,
  effect_allele, other_allele,
  # The stored gamma is the log odds of the lower component; the table reports membership of
  # the upper component, so the sign is flipped (as in Fig. 3).
  gamma_upper = pair(-gamma, gamma_se),
  beta_left = pair(beta_left, beta_left_se),
  beta_right = pair(beta_right, beta_right_se),
  hotelling_mlog10p = rd(hotelling_mlog10p),
  mean_effect = pair(mean_effect, mean_se),
  regenie_beta = pair(regenie_beta, regenie_se),
  mean_over_regenie = rd(mean_over_regenie),
  label = ifelse(!is_p04_lead, paste0(label, " (example, not a lead)"), label),
  mixture_lack_of_fit_Q = rd(Q_mixture, 1),
  alternative_design_gamma_upper = pair(-alt_gamma, alt_gamma_se),
  alternative_design_label = alt_label,
  alternative_design_hotelling_mlog10p = rd(alt_hotelling),
  alternative_design_mixture_Q = rd(alt_Q_mixture, 1),
  location_L = zpair(L, z_location),
  spread_S = zpair(S, z_spread),
  bound_nominal_pct = pct(bound_nominal),
  bound_bonferroni_pct = pct(bound_bonferroni),
  lack_of_fit_Q = rd(Q_ols_location_scale, 1),
  lack_of_fit_p = ifelse(is.na(p_davies_location_scale), NA_character_,
                         ifelse(davies_floored_location_scale %in% TRUE, "<1e-12",
                                formatC(p_davies_location_scale, format = "e", digits = 2)))
)]
setorder(out, trait, chr, bp)

## 3. Save ----------------------------------------------------------------------------------
save_table(
  as.data.frame(out), "tableS4_antibody_leads",
  caption = paste(
    "Table S4. Lead variants for JC polyomavirus and Merkel cell polyomavirus VP1 antibody levels. One row per lead",
    "variant and antigen; n = 9,145; log mean fluorescence intensity. Leads follow the rule used for every trait",
    "(REGENIE two-sided P < 5 × 10⁻⁸, clumped at r² < 0.01 within 10 Mb; Methods); Fig. 3 examples that are not leads",
    "are flagged ‘example, not a lead’. Stage 2 by ordinary least squares with the pooled Stage 1 covariance and",
    "mixture weights of the model design for JCV and the hybrid design for MCV (Supplementary Note 2); the alternative",
    "columns repeat γ, label, joint test and mixture Q under the other design. γ, lower location, upper location: per-allele effects (s.e.) on membership of the upper",
    "component (log odds) and on the lower- and upper-component locations. Joint test: −log₁₀ P of the Hotelling",
    "test of the three parameters. Mean / REGENIE: ratio of the mixture mean effect to the REGENIE effect. Label:",
    "‘level’ if |z| < 1.96 for membership and for the component difference; ‘membership’ if |z| ≥ 5.45 for",
    "membership; otherwise ‘mixed’. Location L and spread S: exact location and scale coordinates (z). Detection",
    "bounds: smallest departure from the common spread-to-location line detectable at nominal and Bonferroni-corrected",
    "thresholds (%). Shift + spread lack of fit: residual of the profile from the location + scale frame (43 degrees of",
    "freedom), with exact P (Davies' method), printed ‘<1e-12’ below the computable floor. Mixture Q: residual from the",
    "fitted mixture (42 degrees of freedom). Membership is weakly identified: its scale rests on the kernel",
    "density near the gap, but its sign and z hardly depend on the design. Leads in the major histocompatibility complex (MHC,",
    "chromosome 6) are in strong linkage disequilibrium and are not independent signals. JCV, JC polyomavirus; MCV,",
    "Merkel cell polyomavirus; s.e., standard error; VP1, major capsid protein VP1."),
  results_used = "04_antibodies/decomposition")
