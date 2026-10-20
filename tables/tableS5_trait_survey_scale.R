#!/usr/bin/env Rscript
# Supplementary Table S5: trait survey and change of scale
#
# What: one row per trait of Fig. 5 (19 survey traits, log BMI and height-adjusted weight):
#   lead and principal-component counts and skew counts on the analysed scale; the Box-Cox
#   power that best normalises the raw trait, its maximum-likelihood estimate, and the power
#   closest to the implied natural scale (pooled and per strong lead); and the counts of
#   spread leads, components and skew-beyond leads when every lead is converted exactly to the
#   raw, log and Box-Cox scales.
# Paper: Supplementary Table S5
# Input: results/06_traits/survey/table1.csv,
#   results/06_traits/scale/{boxcox_lambda,flat_lambda_summary,scale_comparison}.csv
# Output: paper_assets/tableS5_trait_survey_scale.csv, _caption.md
# Runs: locally from results/ in about 1 second

suppressPackageStartupMessages({library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Traits of Fig. 5 ----------------------------------------------------------------------
# Lp(a) is excluded from the paper: censored at both assay limits (config.yaml lists it).
excluded <- unlist(read_config()$estimation$het_pc$exclude_traits)
survey <- as.data.table(read_result("06_traits/survey", "table1"))[in_fig6 == TRUE]
survey <- survey[!arm %in% excluded]
stopifnot(nrow(survey) == 21, sum(survey$n_leads) == 10598,
          sum(survey$n_components == 1) == 5, sum(survey$n_components == 2) == 15, sum(survey$n_components == 3) == 1)

## 2. Scale results -------------------------------------------------------------------------
boxcox <- as.data.table(read_result("06_traits/scale", "boxcox_lambda"))
flat <- as.data.table(read_result("06_traits/scale", "flat_lambda_summary"))
by_scale <- as.data.table(read_result("06_traits/scale", "scale_comparison"))
by_scale_wide <- dcast(by_scale, trait ~ scale, value.var = c("n_spread_lambda2_bonf", "n_components", "n_skew_beyond_bonf"))

## 3. Join (Box-Cox results refer to the raw trait) -------------------------------------------
out <- merge(survey[, .(trait = arm, raw_trait, analysed_scale = scale, n_leads, n_components, n_skew_raw_bonf, n_skew_beyond_bonf)],
             boxcox[, .(raw_trait = trait, lambda_normalising = lambda_main, lambda_at_grid_edge = lambda_main_at_grid_edge,
                        lambda_ml, lambda_ml_lo95, lambda_ml_hi95, lambda_flat_pooled)], by = "raw_trait", all.x = TRUE)
out <- merge(out, flat[, .(raw_trait = trait, n_strong_leads = n_strong, lambda_flat_lead_median_strong = round(lambda_v_median_strong, 2),
                           lambda_flat_lead_iqr_strong = sprintf("%.2f to %.2f", lambda_v_q25_strong, lambda_v_q75_strong))],
             by = "raw_trait", all.x = TRUE)
out <- merge(out, by_scale_wide[, c("trait", grep("^n_", names(by_scale_wide), value = TRUE)), with = FALSE],
             by.x = "raw_trait", by.y = "trait", all.x = TRUE)
setcolorder(out, c("trait", "raw_trait", "analysed_scale"))
setorder(out, trait)

## 4. Save ----------------------------------------------------------------------------------
cap <- paste(
  "Table S5. Trait survey and scale. One row per trait: the 21 traits of Fig. 5 (19 survey traits, log BMI and",
  "height-adjusted weight), on the analysed scale; Box–Cox powers refer to the raw trait (raw BMI for log BMI).",
  "Leads: lead variants. Components: principal components of the whitened lead profiles above the parallel-analysis",
  "null. Skew leads, raw: leads with a third L-moment effect; skew-beyond leads: leads with skew beyond the shift and",
  "spread frame (Bonferroni within trait, analysed scale). λ normalising: Box–Cox power that best normalises the trait",
  "(maximum normal probability-plot correlation over 193 quantile levels, searched from −5 to 5 in steps of 0.01;",
  "‘at search edge’ flags the limit of the search). λ ML: maximum-likelihood power with 95% profile interval (undefined",
  "when any value is ≤ 0). λ flat, pooled: Box–Cox power closest to the implied natural scale, under which the leads'",
  "common profile is a pure shift. Strong leads: leads with |mean z| > 10 on the raw scale; λ flat per lead: median",
  "and interquartile range over strong leads of the power at which a lead's converted profile has no L-moment spread.",
  "Spread leads, components and skew-beyond leads by scale: counts with every lead converted exactly to the Box–Cox,",
  "log and raw scales (Bonferroni within trait). Component counts do not depend on the scale, because a change of scale",
  "rotates all whitened profiles together; spread counts do (BMI 152 raw, 51 log, 11 Box–Cox). BMI, body mass index;",
  "IQR, interquartile range; ML, maximum likelihood.")
save_table(as.data.frame(out), "tableS5_trait_survey_scale", cap,
           results_used = c("06_traits/survey", "06_traits/scale"))
