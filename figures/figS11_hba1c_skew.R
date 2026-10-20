#!/usr/bin/env Rscript
# Supplementary Figure S11: asymmetry of HbA1c lead profiles against their association with diabetes
#
# What: for each of the 472 log-HbA1c lead variants, the asymmetry left in its quantile-effect
#   profile after a least-squares fit of shift and spread (the residual third L-moment, as a
#   z statistic) is plotted against the variant's association with diabetes. Leads with
#   asymmetry beyond the Bonferroni threshold are mostly diabetes-associated. The six leads of
#   Fig. 5c are circled.
#   Single panel: asymmetry z against diabetes z, one point per lead.
# Paper: Supplementary Fig. S11
# Input: results/07_hba1c/skew_diabetes/{examples,per_lead,tests}.csv
# Output: paper_assets/figS11_hba1c_skew.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load ----------------------------------------------------------------------------------
examples <- as.data.table(read_result("07_hba1c/skew_diabetes", "examples"))
per_lead <- as.data.table(read_result("07_hba1c/skew_diabetes", "per_lead"))
tests <- as.data.table(read_result("07_hba1c/skew_diabetes", "tests"))
tv <- function(k) tests[quantity == k, value]
z_bonf <- qnorm(1 - 0.025 / nrow(per_lead))       # two-sided Bonferroni threshold over the leads

## 2. Plot ----------------------------------------------------------------------------------
per_lead[, grp := fifelse(skew_hit, "hit", "other")]
fig <- ggplot(per_lead, aes(diab_z, skew_beyond_z)) +
  geom_hline(yintercept = c(-z_bonf, z_bonf), linetype = "22", linewidth = 0.3, colour = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_point(aes(colour = grp), size = 0.9, alpha = 0.8) +
  geom_point(data = per_lead[snp_id %in% examples$snp_id], shape = 21, size = 2.2, colour = "black", fill = NA, stroke = 0.4) +
  ggrepel::geom_text_repel(data = merge(per_lead, examples[, .(snp_id, role)], by = "snp_id"),
                           aes(label = snp_id), size = 6 / .pt, colour = "grey20", min.segment.length = 0,
                           seed = 1, max.overlaps = 20) +
  scale_colour_manual(values = c(hit = paper_palette[["spread"]], other = "grey60"),
                      labels = c(hit = "Asymmetry", other = "Other"), breaks = c("hit", "other"), name = NULL) +
  labs(x = "Diabetes association (z)", y = "Asymmetry beyond shift and spread (z)") +
  theme_paper() + legend_inside(0.02, 0.98, c(0, 1))

## 3. Save ----------------------------------------------------------------------------------
src <- per_lead[, .(snp_id, diab_z, skew_beyond_z, skew_hit, fig5c_example = snp_id %in% examples$snp_id)]
n_case <- tv("diab_cases_log_hba1c_sample"); n_ctrl <- tv("diab_controls_log_hba1c_sample")
fmtn <- function(x) format(x, big.mark = ",")
cap <- paste(
  "Figure S11. Asymmetry beyond shift and spread in HbA1c lead variants against their association with diabetes.",
  sprintf("Each point is one of %d lead variants for log HbA1c (HbA1c-increasing allele; n = %s participants).", nrow(per_lead),
          fmtn(n_case + n_ctrl)),
  "The vertical axis is the z statistic for asymmetry beyond shift and spread, the third L-moment of the residual after",
  "a least-squares fit of shift and spread to the lead's quantile-effect profile (Methods); the horizontal axis is the z",
  sprintf("statistic for association with diabetes (logistic regression; %s cases, %s controls).", fmtn(n_case), fmtn(n_ctrl)),
  sprintf("Dashed lines mark the two-sided Bonferroni threshold for asymmetry (P = 0.05/%d); Asymmetry, leads beyond it;", nrow(per_lead)),
  "circled and labelled points are the six leads shown in Fig. 5c.",
  sprintf("Of the %d leads with asymmetry beyond shift and spread, %d (%.0f%%) were associated with diabetes at two-sided",
          tv("n_skew_hits"), tv("n_diab_p_lt_bonf472_hit"), 100 * tv("frac_diab_p_lt_bonf472_hit")),
  sprintf("P < 0.05/%d, against %d (%.0f%%) of the other %d; Spearman correlation of the two z statistics, ρ = %.2f.",
          nrow(per_lead), tv("n_diab_p_lt_bonf472_nonhit"), 100 * tv("frac_diab_p_lt_bonf472_nonhit"), nrow(per_lead) - tv("n_skew_hits"),
          tv("spearman_skew_z_vs_diab_z")),
  sprintf("In a separate analysis of the %s participants without diabetes, %d of these %d leads retained the asymmetry.",
          fmtn(n_ctrl), tv("n_nd_skew_hit_same_snp_hit"), tv("n_skew_hits")),
  "Both tests are two-sided and Bonferroni-corrected for the number of leads. HbA1c, glycated haemoglobin.")
save_figure(fig, "figS11_hba1c_skew", src, cap, width = 120, height = 100,
            results_used = "07_hba1c/skew_diabetes")
