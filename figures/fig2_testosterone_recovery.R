#!/usr/bin/env Rscript
# Figure 2: recovery of sex-specific testosterone effects from a sex-blind analysis
#
# What: the mixture of Fig. 1 turns one sex-blind quantile GWAS into a female effect, a male
#   effect and a mixing (membership) effect for every variant. Here those are checked against
#   sex-stratified REGENIE GWAS in the same participants, and the cost is compared with
#   fitting the mixture directly by maximum likelihood at each variant.
#   a: G-MoQ female and male effects against the sex-stratified reference at reference leads.
#   b: LD score regression genetic correlation of each component with its reference.
#   c: computing time per additional variant against sample size.
#   d, e: two exemplar variants chosen by a rule fixed in advance (estimates and profiles).
# Paper: Fig. 2
# Input: results/03_testosterone/recovery/{concordance_leads,concordance,ldsc_rg,
#   exemplar_selected,exemplar_estimates,exemplar_profiles}.csv,
#   results/03_testosterone/exemplars/exemplar_sex_reference.csv,
#   results/03_testosterone/speed/per_snp_costs.csv
# Output: paper_assets/fig2_testosterone_recovery.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 3 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

cfg <- read_config()
arm_lab <- c(testosterone = "No shrinkage", testosterone_shrink_0_2 = "Shrink 0.2",
             testosterone_shrink_0_5 = "Shrink 0.5", testosterone_shrink_0_7 = "Shrink 0.7")
comp_col <- c(female = paper_palette[["female"]], male = paper_palette[["male"]],
              gamma = paper_palette[["mixture"]], total = "grey10")
unshrunk <- "testosterone"                         # the arm shown in panels d, e
unshrunk_design <- mixture_design(cfg, unshrunk)

# Testosterone results, keeping the reported weight design of each arm.
recovery <- function(name) {
  x <- as.data.table(read_result("03_testosterone/recovery", name))
  if ("arm" %in% names(x)) pick_design(x, "arm", cfg = cfg) else pick_design(x, fit = unshrunk, cfg = cfg)
}

## 1. Panel a: G-MoQ against the sex-stratified reference at reference leads ----------------
leads <- recovery("concordance_leads")[reference == "combined" & component %in% c("female", "male")]
leads[, arm_f := factor(arm_lab[arm], levels = arm_lab)]
# r and the least-squares slope are computed on the plotted lead variants themselves.
conc_leads <- leads[, .(pearson_beta = cor(beta_regenie, beta_gmoq),
                        slope_gmoq_on_regenie = unname(coef(lm(beta_gmoq ~ beta_regenie))[2]), n = .N),
                    by = .(arm, arm_f, component)]
comp_lab <- c("Female effect", "Male effect")
leads[, comp_f := factor(component, levels = c("female", "male"), labels = comp_lab)]
conc_leads[, comp_f := factor(component, levels = c("female", "male"), labels = comp_lab)]
# Arms on hybrid weights are outlined.
hybrid_arms <- names(arm_lab)[mixture_design(cfg, names(arm_lab)) == "hybrid"]
outline <- CJ(arm_f = factor(arm_lab[hybrid_arms], levels = arm_lab), comp_f = factor(comp_lab, levels = comp_lab))
pA <- ggplot(leads, aes(beta_regenie, beta_gmoq, colour = component)) +
  geom_rect(data = outline, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
            fill = NA, colour = "grey45", linetype = "22", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_point(size = 0.7, alpha = 0.65, stroke = 0) +
  geom_text(data = conc_leads, aes(x = -Inf, y = Inf, label = sprintf("r = %.3f\nslope = %.2f", pearson_beta,
            slope_gmoq_on_regenie)), hjust = -0.08, vjust = 1.12, size = 2.0, colour = "grey15", lineheight = 0.95) +
  facet_grid(comp_f ~ arm_f) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05)) +
  scale_colour_manual(values = comp_col, guide = "none") +
  labs(x = "Sex-stratified reference REGENIE beta", y = "G-MoQ Stage 2 beta", tag = "a") +
  theme_paper() + theme(strip.text.y.right = element_text(angle = -90, hjust = 0.5), panel.spacing = unit(4, "pt"))

## 2. Panel b: genetic correlation with the reference -------------------------------------
# Female and male components against the sex-stratified GWAS; membership (-gamma) against a
# GWAS of chromosomal sex in the same participants.
rg <- recovery("ldsc_rg")
arm_disp <- ifelse(names(arm_lab) %in% hybrid_arms, paste0(arm_lab, "\n(hybrid)"), arm_lab)
rg[, arm_f := factor(arm_disp[match(arm, names(arm_lab))], levels = arm_disp)]
rg[, comp_f := factor(component, levels = rev(c("female", "male", "sex")),
                      labels = rev(c("Female effect", "Male effect", "Membership (-gamma)\nvs sex GWAS")))]
pB <- ggplot(rg, aes(arm_f, comp_f, fill = rg)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", rg), colour = rg > 0.93), size = 2.5) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
  # rg can exceed 1 by sampling error, so the colour scale is capped at 1.
  scale_fill_gradient(low = "#f1e6d8", high = "#2f6690", limits = c(min(0.7, min(rg$rg)), 1),
                      oob = scales::squish, name = "rg") +
  labs(x = NULL, y = NULL, tag = "b") +
  theme_paper() + theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                        legend.position = "right", legend.key.height = unit(1.2, "lines"),
                        legend.key.width = unit(0.5, "lines"))

## 3. Panel c: cost per additional variant ------------------------------------------------
# Seconds per variant on one core; one-off per-trait costs (baseline mixture fit, Stage 1
# density) are excluded.
costs <- as.data.table(read_result("03_testosterone/speed", "per_snp_costs"))
methods <- c("Direct mixture ML", "G-MoQ Stage 1 + Stage 2", "Stage 2 alone (given Stage 1)")
method_col <- setNames(c(paper_palette[["mixture"]], paper_palette[["gmoq"]], "#7fb3d5"), methods)
timing <- rbind(costs[, .(n, method = methods[1], s = ml_per_snp)],
                costs[, .(n, method = methods[2], s = stage1_per_snp + stage2_per_snp_param_gwas)],
                costs[, .(n, method = methods[3], s = stage2_per_snp_param_gwas)])
timing[, method := factor(method, levels = methods)]
fmt_time <- function(s) ifelse(s >= 1, sprintf("%.0f s", s), ifelse(s >= 1e-3, sprintf("%.0f ms", s * 1e3), sprintf("%.2f ms", s * 1e3)))
at_400k <- timing[n == 4e5][order(method)]
pC <- ggplot(timing, aes(n, s, colour = method)) +
  geom_line(linewidth = 0.7, lineend = "round") + geom_point(size = 1.3) +
  geom_text(data = at_400k, aes(label = fmt_time(s)), hjust = -0.3, size = 1.95, show.legend = FALSE) +
  scale_x_log10(breaks = c(1e4, 5e4, 4e5), labels = c("10k", "50k", "400k"), expand = expansion(mult = c(0.05, 0.22))) +
  scale_y_log10(breaks = 10^(-4:1), labels = c("0.1 ms", "1 ms", "10 ms", "0.1 s", "1 s", "10 s"), limits = c(5e-5, 15)) +
  scale_colour_manual(values = method_col, guide = "none") +
  coord_cartesian(clip = "off") +   # end labels sit in the right-hand x expansion
  labs(x = "Sample size", y = "Seconds per SNP", tag = "c") +
  theme_paper() + theme(plot.caption = element_text(hjust = 0, size = 5.6, lineheight = 0.95, colour = "grey15"),
                        plot.caption.position = "plot")

## 4. Panels d, e: two exemplar variants (unshrunk arm) -----------------------------------
selected <- recovery("exemplar_selected")
estimates <- recovery("exemplar_estimates")[arm == unshrunk]
profiles <- recovery("exemplar_profiles")[arm == unshrunk]
sex_ref <- as.data.table(read_result("03_testosterone/exemplars", "exemplar_sex_reference"))
# Selection categories of the pre-written rule: a variant with female, male and mixing
# effects all |z| > 1.96, and a male-specific variant.
exemplars <- c(shared_concordant = "Mixed-component signal", male_specific = "Male-specific signal")
forest_col <- c(ref = "grey25", female = comp_col[["female"]], male = comp_col[["male"]], gamma = comp_col[["gamma"]])
profile_col <- c(total = "grey10", female = comp_col[["female"]], male = comp_col[["male"]], gamma = comp_col[["gamma"]])
profile_lab <- c(total = "total", female = "female", male = "male", gamma = "mixing")
ex_rows <- list(); src_ex <- list()
for (k in seq_along(exemplars)) {
  s <- selected[category == names(exemplars)[k]]
  e <- estimates[snp_id == s$snp_id]; pr <- profiles[snp_id == s$snp_id]
  sx <- sex_ref[snp_id == s$snp_id & weight_design == unshrunk_design]
  # Membership is shown as -gamma, signed to match the sex GWAS (male = 1).
  fd <- rbind(
    data.table(label = "Female: reference REGENIE", est = e$regenie_combined_female_beta, se = e$regenie_combined_female_se, ck = "ref", block = "f"),
    data.table(label = "Female: G-MoQ", est = e$beta_left, se = e$beta_left_se, ck = "female", block = "f"),
    data.table(label = "Male: reference REGENIE", est = e$regenie_combined_male_beta, se = e$regenie_combined_male_se, ck = "ref", block = "m"),
    data.table(label = "Male: G-MoQ", est = e$beta_right, se = e$beta_right_se, ck = "male", block = "m"),
    data.table(label = "Membership: sex GWAS", est = sx$sex_beta, se = sx$sex_se, ck = "ref", block = "g"),
    data.table(label = "Membership: G-MoQ (-gamma)", est = -e$gamma, se = e$gamma_se, ck = "gamma", block = "g"))
  obs <- pr[, .(tau, beta = beta_stage1, se = se_stage1)]
  fits <- rbind(pr[, .(tau, series = "total", value = fitted_total)], pr[, .(tau, series = "female", value = fitted_beta_left)],
                pr[, .(tau, series = "male", value = fitted_beta_right)], pr[, .(tau, series = "gamma", value = fitted_gamma)])
  # Few participants lie between the modes, so the Stage 1 SE near tau = 0.5 is huge; the
  # display range is set from the quantile levels away from that gap.
  away <- pr[abs(tau - 0.5) > 0.03]
  yl <- range(c(away$beta_stage1 - 1.96 * away$se_stage1, away$beta_stage1 + 1.96 * away$se_stage1, pr$fitted_total[abs(pr$tau - 0.5) > 0.03]))
  yl <- yl + c(-0.08, 0.08) * diff(yl)
  pf <- snp_forest(fd, forest_col, xlab = "Per-allele effect on log testosterone (95% CI)")
  pp <- snp_profile(obs, fits, profile_col, labels = profile_lab, ylim = yl,
                    legend_pos = if (k == 1) c(0.02, 0.98) else c(0.02, 0.03), legend_just = if (k == 1) c(0, 1) else c(0, 0)) +
    labs(y = "Stage 1 beta_tau (95% CI)")
  # The rsID is printed above the forest plot (the theme blanks titles, so it is re-enabled here).
  ex_rows[[k]] <- (pf + labs(tag = letters[3 + k]) + ggtitle(s$snp_id) +
                     theme(plot.title = element_text(size = 7, face = "plain", hjust = 0,
                                                     margin = margin(b = 2))) | pp) +
    plot_layout(widths = c(1, 1.45))
  src_ex[[k]] <- rbind(fd[, .(panel = LETTERS[3 + k], group = s$snp_id, series = label, x = NA_real_, y = est, se)],
                       obs[, .(panel = LETTERS[3 + k], group = s$snp_id, series = "stage1", x = tau, y = beta, se)],
                       fits[, .(panel = LETTERS[3 + k], group = s$snp_id, series = paste0("fit_", series), x = tau, y = value, se = NA_real_)])
}

## 5. Assemble and save ---------------------------------------------------------------------
fig <- pA / ((pB | pC) + plot_layout(widths = c(1.15, 1))) / ex_rows[[1]] / ex_rows[[2]] +
  plot_layout(heights = c(1.25, 0.8, 1, 1))
src <- rbind(
  leads[, .(panel = "A", group = paste(arm, component), series = snp_id, x = beta_regenie, y = beta_gmoq, se = se_gmoq)],
  rg[, .(panel = "B", group = arm, series = component, x = NA_real_, y = rg, se = rg_se)],
  timing[, .(panel = "C", group = as.character(method), series = "seconds_per_snp", x = n, y = s, se = NA_real_)],
  rbindlist(src_ex))
cap <- paste(
  "Figure 2. Recovery of sex-specific testosterone effects from a sex-blind analysis.",
  "n = 399,710; the sex-stratified reference is REGENIE in the same 199,572 female and 200,138 male participants.",
  "Stage 2 is the mixture model of Fig. 1, fitted by ordinary least squares with the pooled Stage 1 covariance.",
  "**a**, G-MoQ female and male estimates against the sex-stratified reference at reference lead variants",
  "(two-sided P < 5 x 10-8), by degree of shrinkage; r is the Pearson correlation and the slope is from ordinary least",
  "squares, both over the plotted leads; the dashed line is y = x; outlined panels use hybrid weights (Methods). **b**, LD score regression genetic correlation of each G-MoQ component with its",
  "reference: sex-stratified GWAS for the female and male components, and a GWAS of chromosomal sex in the same",
  "participants for membership. **c**, Computational cost per additional variant on one core against sample size, for",
  "direct mixture maximum likelihood, G-MoQ Stages 1 and 2, and Stage 2 alone given the Stage 1 output.",
  "**d**,**e**, Two exemplar loci selected by a rule fixed in advance: rs1260326, with effects on all three",
  "parameters, and rs7097461, with a male-specific effect. Left, per-allele estimates with 95% confidence intervals;",
  "right, Stage 1 quantile effects (grey, 95% confidence intervals) with the fitted profile (black) and each",
  "parameter's contribution. Few participants lie between the modes, so effects near tau = 0.5 are very imprecise and",
  "values beyond the axis are drawn at its edge as triangles. Membership is plotted as -gamma, signed to match the",
  "GWAS of sex.")
save_figure(fig, "fig2_testosterone_recovery", src, cap, width = 183, height = 235,
            results_used = c("03_testosterone/recovery", "03_testosterone/exemplars", "03_testosterone/speed"))
