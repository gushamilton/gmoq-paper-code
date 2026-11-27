#!/usr/bin/env Rscript
# Fig 7: Scale-invariant PC2 uncovers distinct biological mechanisms across 2D traits.
#   a  ALT (316 leads): scale-invariant PC2 score against SNP x BMI interaction z (r = 0.75),
#      coloured by recorded liver-disease association (Bonferroni P < 0.05/316); liver-disease genes labelled.
#   b  Urate (396 leads, sex-adjusted): scale-invariant PC2 score against SNP x sex interaction z
#      (r = 0.82); SLC2A9 and ABCG2 labelled.
#   c  FEV1/FVC (421 leads): scale-invariant PC2 score against SNP x ever-smoking interaction z
#      (r = -0.27), isolating the smoking-dependent CHRNA3/CHRNA5 locus (rs12914385).
#   d  LDL (230 leads): quantile-effect profiles at APOE (rs7254892) and PCSK9 (rs11591147) in
#      all participants (showing lower-tail statin distortion) vs participants reporting no lipid medication.
# Reads aggregate results/ only (P11_moments, aux_hidden_biology, aux_more_traits).
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))

d0 <- file.path(.repo, "results", "06_traits", "moments")
d1 <- file.path(.repo, "results", "06_traits", "hidden_biology")
d2 <- file.path(.repo, "results", "06_traits", "hidden_biology")
BL <- paper_palette[["shift"]]
CO <- paper_palette[["spread"]]

lm_dt <- fread(file.path(d0, "lead_moments.csv"),
               select = c("arm", "snp_id", "whitened_norm", "pc2_score"))
lm_dt[, pc2_z := whitened_norm * pc2_score]

# ---- a: ALT (MASLD + G x BMI) ----
alt <- merge(fread(file.path(d1, "alt_masld_per_lead.csv")),
             lm_dt[arm == "log_alt", .(snp_id, pc2_z)], by = "snp_id")
alt[, masld_sig := factor(fifelse(masld_any_p < 0.05 / nrow(alt), "Liver-disease associated", "Other"),
                          levels = c("Other", "Liver-disease associated"))]
r_alt <- cor(alt$pc2_z, alt$bmi_int_z, use = "complete.obs")

alt_labs <- data.table(
  snp_id = c("rs738408", "rs10401969", "rs6832102", "rs10787429", "rs641738", "rs2980875"),
  gene   = c("PNPLA3",   "TM6SF2",     "HSD17B13",  "GPAM",       "MBOAT7",   "TRIB1"),
  lx     = c(20.5,       12.2,         9.8,         11.0,         -2.8,       13.5),
  ly     = c(28.7,       17.5,         13.0,        8.2,          9.2,        3.0),
  hj     = c(1,          0,            0,           0,            1,          0)
)
alt_lab_dt <- merge(alt_labs, alt, by = "snp_id")

pA <- ggplot(alt, aes(bmi_int_z, pc2_z)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(data = alt[masld_sig == "Other"], aes(colour = masld_sig), size = 1.05, alpha = 0.55, stroke = 0) +
  geom_point(data = alt[masld_sig == "Liver-disease associated"], aes(colour = masld_sig), size = 1.75, alpha = 0.9, stroke = 0) +
  geom_point(data = alt_lab_dt, shape = 21, fill = NA, colour = "grey15", size = 2.4, stroke = 0.5) +
  geom_segment(data = alt_lab_dt, aes(x = bmi_int_z, y = pc2_z, xend = lx, yend = ly),
               colour = "grey45", linewidth = 0.25) +
  geom_text(data = alt_lab_dt, aes(x = lx, y = ly, label = gene, hjust = hj),
            fontface = "italic", colour = "grey15", size = 2.1) +
  scale_colour_manual(values = c("Other" = BL, "Liver-disease associated" = CO), name = NULL) +
  coord_cartesian(xlim = c(-8, 26), ylim = c(-8.5, 31)) +
  labs(x = "SNP \u00d7 BMI interaction z statistic", y = "ALT scale-invariant PC2 score") +
  theme_paper() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.04),
        legend.justification = c(1, 0), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.65, "lines")) +
  labs(tag = "a")

# ---- b: Urate (G x Sex despite sex adjustment) ----
ura <- merge(fread(file.path(d2, "urate_sex_per_lead.csv")),
             lm_dt[arm == "log_urate", .(snp_id, pc2_z)], by = "snp_id")
ura[, sex_sig := factor(fifelse(abs(sex_int_z) > qnorm(1 - 0.025 / nrow(ura)), "Sex-interacting", "Other"),
                        levels = c("Other", "Sex-interacting"))]
r_ura <- cor(ura$pc2_z, ura$sex_int_z, use = "complete.obs")

ura_labs <- data.table(
  snp_id = c("rs10805346", "rs7655090", "rs2728109"),
  gene   = c("SLC2A9",     "SLC2A9",    "ABCG2"),
  lx     = c(-33.5,        -14.5,       7.0),
  ly     = c(-28.0,        -16.0,       10.8),
  hj     = c(0,            0,           0)
)
ura_lab_dt <- merge(ura_labs, ura, by = "snp_id")

pB <- ggplot(ura, aes(sex_int_z, pc2_z)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(data = ura[sex_sig == "Other"], aes(colour = sex_sig), size = 1.05, alpha = 0.55, stroke = 0) +
  geom_point(data = ura[sex_sig == "Sex-interacting"], aes(colour = sex_sig), size = 1.65, alpha = 0.9, stroke = 0) +
  geom_point(data = ura_lab_dt, shape = 21, fill = NA, colour = "grey15", size = 2.4, stroke = 0.5) +
  geom_segment(data = ura_lab_dt, aes(x = sex_int_z, y = pc2_z, xend = lx, yend = ly),
               colour = "grey45", linewidth = 0.25) +
  geom_text(data = ura_lab_dt, aes(x = lx, y = ly, label = gene, hjust = hj),
            fontface = "italic", colour = "grey15", size = 2.1) +
  scale_colour_manual(values = c("Other" = BL, "Sex-interacting" = CO), name = NULL) +
  coord_cartesian(xlim = c(-40, 15), ylim = c(-33, 13.5)) +
  labs(x = "SNP \u00d7 sex interaction z statistic", y = "Urate scale-invariant PC2 score") +
  theme_paper() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.04),
        legend.justification = c(1, 0), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.65, "lines")) +
  labs(tag = "b")

# ---- c: FEV1/FVC (G x Smoking, CHRNA3/CHRNA5) ----
fev <- merge(fread(file.path(d2, "fev_smoking_per_lead.csv")),
             lm_dt[arm == "fev1_fvc_ratio", .(snp_id, pc2_z)], by = "snp_id")
fev[, sm_sig := factor(fifelse(abs(smok_int_z) > qnorm(1 - 0.025 / nrow(fev)), "Smoking-interacting", "Other"),
                       levels = c("Other", "Smoking-interacting"))]
r_fev <- cor(fev$pc2_z, fev$smok_int_z, use = "complete.obs")
fev_lab <- fev[snp_id == "rs12914385"]

pC <- ggplot(fev, aes(smok_int_z, pc2_z)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(data = fev[sm_sig == "Other"], aes(colour = sm_sig), size = 1.05, alpha = 0.55, stroke = 0) +
  geom_point(data = fev[sm_sig == "Smoking-interacting"], aes(colour = sm_sig), size = 1.8, alpha = 0.9, stroke = 0) +
  geom_point(data = fev_lab, shape = 21, fill = NA, colour = "grey15", size = 2.6, stroke = 0.55) +
  geom_segment(data = fev_lab, aes(x = smok_int_z, y = pc2_z, xend = 9.0, yend = -5.6),
               colour = "grey45", linewidth = 0.25) +
  annotate("text", x = 8.9, y = -5.4, label = "CHRNA3/CHRNA5", fontface = "italic", hjust = 1, size = 2.15, colour = "grey15") +
  scale_colour_manual(values = c("Other" = BL, "Smoking-interacting" = CO), name = NULL) +
  coord_cartesian(xlim = c(-4.2, 12.2), ylim = c(-8.2, 4.8)) +
  labs(x = "SNP \u00d7 ever-smoking interaction z statistic", y = "FEV1/FVC scale-invariant PC2 score") +
  theme_paper() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.65, "lines")) +
  labs(tag = "c")

# ---- d: LDL statin distortion at APOE (rs7254892) and PCSK9 (rs11591147) ----
pr <- fread(file.path(d1, "ldl_statin_profiles.csv"))[snp_id %in% c("rs7254892", "rs11591147")]
pr[, locus_f := factor(fifelse(snp_id == "rs7254892", 'italic("APOE")~"rs7254892"', 'italic("PCSK9")~"rs11591147"'),
                       levels = c('italic("APOE")~"rs7254892"', 'italic("PCSK9")~"rs11591147"'))]
pr[, arm_f := factor(fifelse(arm == "ldl_all", "All participants", "No reported lipid medication"),
                     levels = c("All participants", "No reported lipid medication"))]
pr[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]

pD <- ggplot(pr, aes(tau, beta, colour = arm_f, fill = arm_f)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(data = smooth_curve(pr, "tau", "beta", by = c("locus_f", "arm_f")),
            aes(x, y, colour = arm_f), inherit.aes = FALSE, linewidth = 0.85, lineend = "round") +
  facet_wrap(~locus_f, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = c("All participants" = CO, "No reported lipid medication" = BL), name = NULL) +
  scale_fill_manual(values = c("All participants" = CO, "No reported lipid medication" = BL), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  coord_cartesian(ylim = c(0.20, 0.60)) +
  labs(x = expression("Quantile level " * tau), y = "Per-allele effect on LDL (mmol/L)") +
  theme_paper() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.65, "lines")) +
  labs(tag = "d")

fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1))

src <- rbindlist(list(
  alt[, .(panel = "a", trait = "ALT", snp_id, x = bmi_int_z, y = pc2_z, group = as.character(masld_sig))],
  ura[, .(panel = "b", trait = "Urate", snp_id, x = sex_int_z, y = pc2_z, group = as.character(sex_sig))],
  fev[, .(panel = "c", trait = "FEV1/FVC", snp_id, x = smok_int_z, y = pc2_z, group = as.character(sm_sig))],
  pr[,  .(panel = "d", trait = "LDL", snp_id, series = as.character(arm_f), x = tau, y = beta, lo, hi)]
), fill = TRUE)

cap <- paste(
  "Figure 7. Scale-invariant PC2 scores uncover distinct biological mechanisms across two-dimensional traits.",
  "UK Biobank lead variants; scale-invariant PC2 scores are computed from Stage 1 quantile-effect profiles blind to",
  "clinical diagnoses, environmental modifiers and medication use (Methods).",
  "**a**, ALT (316 lead variants; n = 439,662): scale-invariant PC2 score against SNP x BMI interaction z (Pearson r = 0.745,",
  "nominal two-sided P = 3.8 x 10^-57), coloured by association with recorded liver disease (13,882 cases among 441,106 participants; logistic",
  "regression, two-sided Bonferroni P < 0.05/316). Labelled loci include PNPLA3, TM6SF2 and HSD17B13.",
  "**b**, Log urate (396 lead variants; n = 440,735; sex-adjusted in Stage 1): scale-invariant PC2 score against SNP x sex",
  "interaction z (Pearson r = 0.822, nominal two-sided P = 2.9 x 10^-98), coloured by Bonferroni SNP x sex significance (two-sided",
  "P < 0.05/396); female-biased SLC2A9 variants separate cleanly from ABCG2.",
  "**c**, FEV1/FVC (421 lead variants; n = 344,396): scale-invariant PC2 score against SNP x ever-smoking interaction z",
  "(Pearson r = -0.269, nominal two-sided P = 2.0 x 10^-8), coloured by Bonferroni SNP x ever-smoking significance (two-sided",
  "P < 0.05/421); the smoking-dependent CHRNA3/CHRNA5 locus (rs12914385; PC2 = -7.35, interaction z = 11.4) stands apart",
  "from all other variants.",
  "**d**, Per-allele quantile-effect profiles (pointwise 95% confidence intervals) on LDL at APOE (rs7254892) and PCSK9 (rs11591147)",
  "in all participants (n = 440,443) and those without reported cholesterol-lowering medication (n = 364,354, excluding 76,089 participants reporting cholesterol-lowering",
  "medication). Profiles at both loci change after exclusion. ALT, alanine aminotransferase; BMI, body mass index; FEV1/FVC, forced",
  "expiratory volume in 1 s to forced vital capacity ratio; LDL, low-density lipoprotein cholesterol."
)

save_figure(fig, "fig7_hidden_biology", src, cap, width = 183, height = 148,
            results_used = c("06_traits/moments", "06_traits/hidden_biology"))
cat("Saved fig7_hidden_biology cleanly!\n")
