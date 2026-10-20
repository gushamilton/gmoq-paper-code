#!/usr/bin/env Rscript
# Supplementary Figure S8: internal Mendelian randomisation between G-MoQ sex components and REGENIE
#
# What: a check that the female and male components recovered from the sex-blind analysis are
#   on the same scale as the sex-stratified reference. For each testosterone arm, an
#   inverse-variance-weighted slope relates the G-MoQ component effect to the reference effect
#   in the same sex, in both directions (instruments from REGENIE leads, or from G-MoQ leads);
#   a slope of 1 is exact agreement.
#   Single panel: IVW slope with 95% CI by arm, sex (rows) and direction (columns).
# Paper: Supplementary Fig. S8
# Input: results/03_testosterone/recovery/internal_mr.csv
# Output: paper_assets/figS08_internal_mr.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load (reported weight design per arm: hybrid unshrunk, model for the shrunk arms) ------
mr <- pick_design(read_result("03_testosterone/recovery", "internal_mr"), "arm")
arm_labels <- c(testosterone = "Unshrunk", testosterone_shrink_0_2 = "Shrink 0.2",
                testosterone_shrink_0_5 = "Shrink 0.5", testosterone_shrink_0_7 = "Shrink 0.7")
mr[, arm_label := factor(arm_labels[arm], levels = rev(unname(arm_labels)))]
mr[, component_label := factor(ifelse(component == "female", "Female", "Male"), levels = c("Male", "Female"))]
stopifnot(all(mr$direction %in% c("REGENIE -> G-MoQ", "G-MoQ -> REGENIE")))
# Plotmath strip labels so the arrow renders in the PDF as well as the PNG.
dir_labels <- c(`REGENIE -> G-MoQ` = '"REGENIE" %->% "G-MoQ"', `G-MoQ -> REGENIE` = '"G-MoQ" %->% "REGENIE"')
mr[, direction_f := factor(dir_labels[direction], levels = unname(dir_labels))]

## 2. Plot: slopes with 95% CI; numbers are the instrument counts ---------------------------
fig <- ggplot(mr, aes(x = beta, y = arm_label, colour = component_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.22, linewidth = 0.45, orientation = "y") +
  geom_point(size = 1.9) +
  geom_text(aes(label = nsnp, x = ci_high), hjust = -0.25, size = 6 / .pt,
            colour = "grey30", show.legend = FALSE) +
  facet_grid(component_label ~ direction_f, labeller = labeller(direction_f = label_parsed)) +
  scale_colour_manual(values = c(Female = paper_palette[["female"]], Male = paper_palette[["male"]]), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.2))) +
  labs(x = "IVW slope (95% CI)", y = NULL) +
  theme_paper() +
  theme(panel.spacing.x = unit(1.1, "lines"))

## 3. Save ----------------------------------------------------------------------------------
src <- mr[, .(arm, weight_design, component, direction, reference, nsnp, beta, se, ci_low, ci_high)]
cap <- paste(
  "Figure S8. Internal Mendelian randomisation between G-MoQ sex components and sex-stratified REGENIE.",
  "Each point is the inverse-variance-weighted slope (fixed effect, through the origin) relating the G-MoQ Stage 2",
  "female or male component effect on testosterone (rows) to the sex-stratified REGENIE effect in the same sex,",
  "for the unshrunk analysis and three shrinkage settings (Methods); bars, 95% confidence intervals; numbers,",
  "instrument variants per estimate. Columns give the direction: REGENIE → G-MoQ uses independent REGENIE lead variants",
  "as instruments and G-MoQ → REGENIE uses independent genome-wide significant (P < 5 × 10⁻⁸) Stage 2",
  "lead variants (Methods). Dashed line, slope = 1, exact agreement with the sex-stratified reference.",
  "n = 399,710 participants (199,572 female, 200,138 male); the reference uses the same phenotype, covariates and",
  "leave-one-chromosome-out predictions as Stage 1. Comparisons with the reference under both mixture weight designs",
  "are given in Supplementary Table S3; the designs are compared in Supplementary Fig. S7.",
  "CI, confidence interval; IVW, inverse-variance weighted.")
save_figure(fig, "figS08_internal_mr", source = src, caption = cap, width = 136, height = 100,
            results_used = "03_testosterone/recovery")
