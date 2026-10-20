#!/usr/bin/env Rscript
# Figure 5: heterogeneity of genetic effects between lead variants across quantitative traits
#
# What: for each trait, do the lead variants share one quantile-effect shape, or do they
#   differ? Panel a measures how much the leads' change in spread per unit change in
#   location varies around the trait's common ratio; panel b asks whether a second shape
#   (principal component) stands out from noise; panel c shows, for log HbA1c, leads matched
#   on mean effect whose shapes differ, the diabetes-associated ones being asymmetric.
#   a: between-lead heterogeneity tau per trait, on the analysed and on other scales.
#   b: share of shape variance in the second principal component relative to its null.
#   c: six log-HbA1c lead profiles with single-mechanism least-squares fits.
# Paper: Fig. 5
# Input: results/06_traits/survey/fig6_components.csv, results/06_traits/moments/
#   {heterogeneity,pca_parallel}.csv, results/06_traits/scale/scale_comparison.csv,
#   results/06_traits/pc_heterogeneity/pc2_share_ratio.csv,
#   results/07_hba1c/skew_diabetes/{examples,profiles}.csv,
#   results/07_hba1c/mixture/weight_matrices.csv
# Output: paper_assets/fig5_heterogeneity_across_traits.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Traits and labels -------------------------------------------------------------------
traits <- as.data.table(read_result("06_traits/survey", "fig6_components"))
# Lp(a) is censored at both reportable limits (10.2% below, 7.0% above), so its profile shape
# is an assay artefact rather than a genetic one; it is excluded.
traits <- traits[trait != "log_lipoprotein_a"]
n_traits <- nrow(traits)

nice <- function(x) {
  x <- sub("^log_", "log ", x); x <- gsub("_", " ", x)
  x <- sub("lipoprotein a", "Lp(a)", x); x <- sub("hba1c nodiab", "HbA1c (no diabetes)", x); x <- sub("hba1c", "HbA1c", x)
  x <- sub("fev1 fvc ratio", "FEV1/FVC", x); x <- sub("bilirubin total", "bilirubin", x); x <- sub("ldl", "LDL", x)
  x <- sub("hdl", "HDL", x); x <- sub("crp", "CRP", x); x <- sub("alt", "ALT", x); x <- sub("ggt", "GGT", x)
  x <- sub("igf1", "IGF-1", x); x <- sub("bmi", "BMI", x); x <- sub("mch", "MCH", x); x <- sub("monocyte pct", "monocyte %", x)
  x
}
comp_col <- c(`1` = paper_palette[["shift"]], `2` = "grey30", `3` = paper_palette[["spread"]])

## 2. Panel a: between-lead heterogeneity --------------------------------------------------
# tau is the standard deviation across leads of spread per unit shift (unit-SE basis) around
# the trait's common ratio, from Cochran's Q against chi-square(n - 1).
het <- as.data.table(read_result("06_traits/moments", "heterogeneity"))[
  basis == "unit_se" & component == "c2_on_c1" & excess == "n_minus_1" & arm %in% traits$trait]
stopifnot(nrow(het) == n_traits)
scales_cmp <- as.data.table(read_result("06_traits/scale", "scale_comparison"))
d <- merge(het[, .(trait = arm, tau, lo = tau_lo95, hi = tau_hi95, p)], traits[, .(trait, n_leads, n_components)],
           by = "trait")
d[, comp := factor(n_components, levels = 1:3)]
setorder(d, tau)
d[, label_f := factor(sprintf("%s (%d)", nice(trait), n_leads), levels = sprintf("%s (%d)", nice(trait), n_leads))]
by_scale <- scales_cmp[, .(trait = fifelse(trait == "bmi", "log_bmi", trait), scale, tau_scale = tau_unit_se)][trait %in% d$trait]
by_scale <- merge(by_scale, d[, .(trait, label_f)], by = "trait")
by_scale[, scale_f := factor(scale, levels = c("raw", "log", "boxcox"), labels = c("raw", "log", "Box-Cox"))]
pHet <- ggplot(d, aes(y = label_f)) +
  geom_point(data = by_scale, aes(x = tau_scale, shape = scale_f), colour = "grey60", size = 1.1, stroke = 0.35) +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = comp), width = 0, linewidth = 0.5, orientation = "y") +
  geom_point(aes(x = tau, colour = comp), size = 1.8) +
  scale_colour_manual(values = comp_col, name = "Components") +
  scale_shape_manual(values = c(raw = 1, log = 2, `Box-Cox` = 5), name = "Scale") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Between-lead heterogeneity (tau)", y = NULL, tag = "a") +
  theme_paper() + legend_inside(0.98, 0.02, c(1, 0)) +
  theme(legend.box = "horizontal", legend.title = element_text(size = 7),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

## 3. Panel b: the second shape relative to noise -----------------------------------------
# Share of whitened shape variance in PC2 over its parallel-analysis null median (a change of
# scale rotates all whitened profiles together, so this ratio does not depend on the scale).
pc2 <- as.data.table(read_result("06_traits/moments", "pca_parallel"))[component_tested == 2 & arm %in% d$trait]
pc2 <- merge(pc2[, .(trait = arm, ratio = observed_pct / null_median_pct, thr = null_p95_pct / null_median_pct)],
             d[, .(trait, label_f, comp)], by = "trait")
# 95% bootstrap intervals over leads, with the null held fixed.
ci <- as.data.table(read_result("06_traits/pc_heterogeneity", "pc2_share_ratio"))[, .(trait, ratio_boot = ratio, lo = lo95, hi = hi95)]
pc2 <- merge(pc2, ci, by = "trait")
stopifnot(nrow(pc2) == n_traits, max(abs(pc2$ratio - pc2$ratio_boot)) < 1e-8)
pPC <- ggplot(pc2, aes(y = label_f)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_point(aes(x = thr), shape = 124, size = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = comp), width = 0, linewidth = 0.5, orientation = "y") +
  geom_point(aes(x = ratio, colour = comp), size = 1.8) +
  scale_colour_manual(values = comp_col, guide = "none") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3), labels = c("1", "1.5", "2", "3")) +
  labs(x = "Second shape relative to noise", y = NULL, tag = "b") +
  theme_paper() + theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank())

## 4. Panel c: HbA1c leads with matched mean effects but different shapes ------------------
# Pairs fixed in advance: the three diabetes-associated leads with the largest skew beyond
# shift + spread (bottom row), each matched on mean-effect z to a lead with no skew beyond
# that frame (top row).
examples <- as.data.table(read_result("07_hba1c/skew_diabetes", "examples"))
hba1c_prof <- as.data.table(read_result("07_hba1c/skew_diabetes", "profiles"))[analysis == "all_participants_log"]
examples[, grp := fifelse(role == "skew_hit", "diabetes", "other")]
setorder(examples, role, rank)                 # matched leads first (top row), then the skew hits
stopifnot(nrow(examples) == 6, all(examples[role == "non_skew_matched", matched_to] == examples[role == "skew_hit", snp_id]))
ex <- merge(hba1c_prof[snp_id %in% examples$snp_id, .(snp_id, tau, beta, se = se_pooled)],
            examples[, .(snp_id, grp, nearest_gene)], by = "snp_id")
strip_lab <- examples[, sprintf('italic("%s")~"%s"', nearest_gene, snp_id)]
ex[, panel_f := factor(sprintf('italic("%s")~"%s"', nearest_gene, snp_id), levels = strip_lab)]
ex[, grp_f := factor(grp, levels = c("other", "diabetes"))]
# One reference fit per row, by least squares on the lead's own profile: a pure shift (a
# constant) for the top row; for the bottom row the HbA1c mixture's membership weight alone,
# i.e. a variant that changes only membership of the upper (diabetes-enriched) component.
W <- as.data.table(read_result("07_hba1c/mixture", "weight_matrices"))[mixture == "sn2_all" & weight_design == "model"]
W[, tk := round(tau, 3)]; ex[, tk := round(tau, 3)]
ex[, w_mem := W$gamma[match(tk, W$tk)]]
stopifnot(!anyNA(ex$w_mem))
ex[, fit := if (grp_f[1] == "other") rep(mean(beta), .N) else w_mem * sum(w_mem * beta) / sum(w_mem^2), by = snp_id]
ex_col <- c(other = paper_palette[["shift"]], diabetes = paper_palette[["spread"]])
pEx <- ggplot(ex, aes(tau, beta)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_line(aes(y = fit), colour = "grey55", linewidth = 0.9, linetype = "22") +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se, fill = grp_f), colour = NA, alpha = 0.2) +
  geom_line(aes(colour = grp_f), linewidth = 0.7, lineend = "round") +
  facet_wrap(~panel_f, nrow = 2, scales = "free_y", labeller = label_parsed) +
  scale_colour_manual(values = ex_col, guide = "none") + scale_fill_manual(values = ex_col, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1")) +
  labs(x = "Quantile level tau", y = "Per-allele effect on log HbA1c", tag = "c") +
  theme_paper()

## 5. Assemble and save ---------------------------------------------------------------------
fig <- ((pHet | pPC) + plot_layout(widths = c(1.4, 1))) / pEx + plot_layout(heights = c(1.6, 1.25))
src <- rbindlist(list(
  d[, .(panel = "A", trait, series = "tau_analysed_scale", x = tau, lo, hi, p)],
  by_scale[, .(panel = "A", trait, series = paste0("tau_", scale), x = tau_scale)],
  pc2[, .(panel = "B", trait, series = "pc2_over_null_median", x = ratio, lo, hi, threshold = thr)],
  ex[, .(panel = "C", trait = as.character(grp_f), series = snp_id, x = tau, y = beta,
         lo = beta - 1.96 * se, hi = beta + 1.96 * se, fit)]), fill = TRUE)
# Upper bound for the caption that is never below the largest P: mantissa rounded up.
p_max <- max(d$p)
p_e <- floor(log10(p_max)); p_m <- ceiling(p_max / 10^p_e * 10) / 10
stopifnot(p_m * 10^p_e >= p_max)
superscript <- function(n) paste(c("⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹")[as.integer(strsplit(as.character(n), "")[[1]]) + 1], collapse = "")
p_bound <- sprintf("%.1f × 10⁻%s", p_m, superscript(abs(p_e)))
cap <- paste(
  "Figure 5. Heterogeneity of genetic effects between lead variants across quantitative traits.",
  sprintf("%d UK Biobank trait analyses (20 traits); %s lead variants in total, with leads per trait in brackets.", n_traits,
          format(sum(d$n_leads), big.mark = ",")),
  "**a**, Between-lead heterogeneity: the standard deviation across leads of the change in spread per unit change in",
  "location, around the trait's common ratio, on the analysed scale with 95% bootstrap confidence intervals; grey",
  "symbols, the same estimate on the raw, log and Box-Cox scales. Colour, principal components above noise.",
  sprintf("Heterogeneity was tested by a one-sided Cochran Q test; every trait has P ≤ %s, not adjusted for multiple",
          p_bound),
  "comparisons. **b**, Share of shape variance explained by the second principal component relative to its",
  "expectation if every lead shared one shape (parallel-analysis null median; tick, null 95th percentile), with 95%",
  "bootstrap intervals over leads (Methods); 1 indicates no second shape, on any scale.",
  "**c**, Per-allele effects of six lead variants on log HbA1c, with 95% confidence intervals.",
  "Bottom, the three diabetes-associated leads with the largest asymmetry beyond shift and spread (diabetes,",
  "logistic regression in the same participants, two-sided P < 0.05/472); top, a lead with no asymmetry beyond the",
  "frame matched to each on mean effect. Dashed grey lines are least-squares fits of a single mechanism: a pure shift",
  "(top) and a change only in membership of the diabetes-enriched mixture component (bottom) (Methods).")
save_figure(fig, "fig5_heterogeneity_across_traits", src, cap, width = 183, height = 190,
            results_used = c("06_traits/survey", "06_traits/moments", "06_traits/scale", "06_traits/pc_heterogeneity",
                             "07_hba1c/skew_diabetes", "07_hba1c/mixture"))
