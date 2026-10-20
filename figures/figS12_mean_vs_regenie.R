#!/usr/bin/env Rscript
# Supplementary Figure S12: agreement of the Stage 1 mean effect with REGENIE
#
# What: the plain average of a variant's 45 Stage 1 quantile effects (weight W = 1) estimates
#   its mean effect. This is compared with REGENIE run in the same participants with the same
#   covariates. Over all variants the slope of one on the other follows the ratio of their
#   standard errors (regression dilution of the noisier estimate); among genome-wide
#   significant variants the slope is close to 1.
#   a: Stage 1 mean against REGENIE per variant, for five traits.
#   b: per trait, all-variant and significant-variant slopes against the median SE ratio.
# Paper: Supplementary Fig. S12
# Input: results/01_data/stage1_mean_vs_regenie/{mean_vs_regenie_points,mean_vs_regenie_summary}.csv
# Output: paper_assets/figS12_mean_vs_regenie.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 6 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

trait_lab <- c(bmi = "BMI", log_bmi = "log BMI", weight = "Weight",
               jcv_vp1_log_all = "JCV VP1", mcv_vp1_log_all = "MCV VP1")
trait_lab_all <- c(trait_lab, testosterone_shrink_0_2 = "Testosterone shrink 0.2", testosterone_shrink_0_5 = "Testosterone shrink 0.5",
                   testosterone_shrink_0_7 = "Testosterone shrink 0.7")

## 1. Load ----------------------------------------------------------------------------------
# Points: 5,000 randomly sampled variants plus every variant with REGENIE P < 1e-5.
points <- as.data.table(read_result("01_data/stage1_mean_vs_regenie", "mean_vs_regenie_points"))[trait %in% names(trait_lab)]
summary_ <- as.data.table(read_result("01_data/stage1_mean_vs_regenie", "mean_vs_regenie_summary"))
points[, trait_f := factor(trait_lab[trait], levels = trait_lab)]
points[, grp := factor(fifelse(sampled, "random", "sig"), levels = c("random", "sig"))]
stopifnot(all(points[sampled == TRUE, .N, by = trait]$N == 5000))
# One row per trait (unshrunk testosterone is in the legend only): slopes and SE ratios over
# all variants and over REGENIE P < 5e-8.
slopes <- dcast(summary_[subset %in% c("all", "regenie_p_lt_5e-8") & trait != "testosterone"], trait ~ subset,
                value.var = c("slope", "median_se_ratio", "n"))
setnames(slopes, gsub("regenie_p_lt_5e-8", "sig", names(slopes)))
slopes[, trait_f := factor(trait_lab_all[trait], levels = trait_lab_all)]

## 2. Panel a: per-variant scatter ------------------------------------------------------
# Corner label: Pearson r over all tested variants.
r_lab <- summary_[subset == "all" & trait %in% names(trait_lab), .(trait_f = factor(trait_lab[trait], levels = trait_lab),
                                                                   lab = sprintf("r = %.3f", pearson_r))]
pA <- ggplot(points, aes(beta_regenie, mean_w1)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_point(aes(colour = grp), size = 0.35, alpha = 0.5) +
  geom_text(data = r_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.4, size = 6 / .pt,
            colour = "grey20") +
  facet_wrap(~trait_f, scales = "free", ncol = 3) +
  scale_colour_manual(values = c(random = "grey55", sig = paper_palette[["gmoq"]]),
                      labels = c(random = "Random", sig = expression(italic(P) < 10^-5)), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = "REGENIE per-allele effect", y = "Stage 1 mean per-allele effect", tag = "a") +
  theme_paper() +
  legend_inside(0.84, 0.25, c(0.5, 0.5))          # key in the empty sixth facet slot

## 3. Panel b: slope against SE ratio ---------------------------------------------------
pB <- ggplot(slopes, aes(median_se_ratio_all, slope_all)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_segment(aes(xend = median_se_ratio_all, yend = slope_sig), colour = "grey75", linewidth = 0.3) +
  geom_point(aes(shape = "all"), size = 1.8, colour = paper_palette[["spread"]]) +
  geom_point(aes(y = slope_sig, shape = "sig"), size = 1.8, colour = paper_palette[["gmoq"]]) +
  ggrepel::geom_text_repel(aes(label = trait_f), size = 6 / .pt, colour = "grey20", min.segment.length = 0,
                           max.overlaps = 30, seed = 1) +
  scale_shape_manual(values = c(all = 16, sig = 17),
                     labels = c(all = "All variants", sig = expression(italic(P) < 5 %*% 10^-8)), name = NULL) +
  guides(shape = guide_legend(override.aes = list(colour = c(paper_palette[["spread"]], paper_palette[["gmoq"]])))) +
  labs(x = "Median SE ratio (Stage 1 mean / REGENIE)", y = "Slope of Stage 1 mean on REGENIE effect",
       tag = "b") +
  theme_paper() + legend_inside(0.55, 0.03, c(0.5, 0))

## 4. Assemble and save ---------------------------------------------------------------------
fig <- pA / pB + plot_layout(heights = c(1.5, 1))
src <- rbind(points[, .(panel = "A", trait, snp_id, x = beta_regenie, y = mean_w1, group = as.character(grp))],
             slopes[, .(panel = "B", trait, snp_id = NA_character_, x = median_se_ratio_all, y = slope_all, group = "all")],
             slopes[, .(panel = "B", trait, snp_id = NA_character_, x = median_se_ratio_all, y = slope_sig, group = "p<5e-8")])
g <- function(t, col) slopes[trait == t][[col]]
testo <- function(col, sub) summary_[trait == "testosterone" & subset == sub][[col]]
fmt <- function(x) format(x, big.mark = ",")
n_sig <- points[sampled == FALSE, .N, by = trait][, setNames(N, trait)]
n_all <- summary_[subset == "all" & trait %in% names(trait_lab), setNames(n, trait)]
cap <- paste(
  "Figure S12. Agreement of the Stage 1 mean effect with the REGENIE effect.",
  "The Stage 1 mean effect (W = 1) is the unweighted average of the 45 quantile effects; REGENIE was run in the",
  "same participants with the same covariates (Methods).",
  "**a**, Each point is a variant: 5,000 randomly sampled variants (grey) and every variant with REGENIE",
  sprintf("P < 10⁻⁵ (blue; BMI %s, log BMI %s, weight %s, JCV VP1 %s, MCV VP1 %s).",
          fmt(n_sig[["bmi"]]), fmt(n_sig[["log_bmi"]]), fmt(n_sig[["weight"]]), fmt(n_sig[["jcv_vp1_log_all"]]),
          fmt(n_sig[["mcv_vp1_log_all"]])),
  sprintf("Dashed line, identity; r, Pearson correlation over all tested variants (%s for BMI and weight traits, %s for VP1).",
          fmt(n_all[["bmi"]]), fmt(n_all[["jcv_vp1_log_all"]])),
  "**b**, Each pair of points is one trait: the least-squares slope of the Stage 1 mean effect on the REGENIE effect",
  "over all variants (circles) and over variants with REGENIE P < 5 × 10⁻⁸ (triangles; grey segment joins the two),",
  "against the median ratio of the two standard errors; testosterone is shown at three shrinkage settings (Methods).",
  "Dashed line, identity.",
  sprintf("All-variant slopes (%.3f–%.3f) follow the SE ratio (%.3f–%.3f; JCV VP1 slope %.3f, SE ratio %.3f);",
          min(slopes$slope_all), max(slopes$slope_all), min(slopes$median_se_ratio_all), max(slopes$median_se_ratio_all),
          g("jcv_vp1_log_all", "slope_all"), g("jcv_vp1_log_all", "median_se_ratio_all")),
  sprintf("slopes at P < 5 × 10⁻⁸ are %.3f–%.3f.", min(slopes$slope_sig), max(slopes$slope_sig)),
  sprintf("Unshrunk testosterone (not shown) has an all-variant slope of %.2f with SE ratio %.2f, and a slope of %.3f at",
          testo("slope", "all"), testo("median_se_ratio", "all"), testo("slope", "regenie_p_lt_5e-8")),
  "P < 5 × 10⁻⁸.",
  "BMI, body mass index; JCV, JC polyomavirus; MCV, Merkel cell polyomavirus; SE, standard error; VP1, viral capsid protein 1.")
save_figure(fig, "figS12_mean_vs_regenie", src, cap, width = 183, height = 190,
            results_used = "01_data/stage1_mean_vs_regenie")
