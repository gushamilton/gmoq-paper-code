#!/usr/bin/env Rscript
# Fig 5: Scale-invariant dimensionality of genetic effects (1D, 2D and 3D archetypes).
#   a  1D trait (IGF-1, 644 leads): profiles across three transformations (log, untransformed,
#      Box-Cox lambda = 0.57). All leads rotate together onto a flat pure shift.
#   b  2D trait (GGT, 454 leads): profiles across three transformations (untransformed, log,
#      Box-Cox lambda = -0.67). Two high-precision leads (rs2650000 and rs12373325) diverge in
#      opposite directions on the trait-flattening transformation; no monotonic transformation
#      can flatten both.
#   c  Exact invariance of whitened principal-component spectra (PCs 1-6) across untransformed,
#      log and Box-Cox transformations for IGF-1 (1D), BMI (2D), GGT (2D), HbA1c (3D) and
#      HbA1c without diabetes (2D), with parallel-analysis null 95th percentiles.
#   d  Scale-invariant whitened PC1 and PC2 scores of the 472 HbA1c leads, separating
#      diabetes-associated loci (TCF7L2, HLA-DOB, CDKAL1) from baseline shift loci (FN3KRP).
# Reads aggregate results/ only (P03_baseline, P05_lead_profiles, P11_moments, P13_scale,
# aux_hba1c_skew_diabetes).
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))

tau <- tau_grid()
bq  <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))
lp  <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
lm  <- as.data.table(read_result("06_traits/moments", "lead_moments"))
pa  <- as.data.table(read_result("06_traits/moments", "pca_parallel"))
diab <- as.data.table(read_result("07_hba1c/skew_diabetes", "per_lead"))
exs  <- as.data.table(read_result("07_hba1c/skew_diabetes", "examples"))
bc  <- grep("^beta_tau_", names(lp), value = TRUE)
stopifnot(length(bc) == length(tau))

BL <- paper_palette[["shift"]]
CO <- paper_palette[["spread"]]
TE <- "#2a7f62"
GR <- "grey30"

# Helper: transform a trait's lead profiles and compute normalized curves + median/IQR + SEs
make_fan <- function(lead_set_name, trait_name, sigma_file, transforms, hi_snps = character(0), n_bg = 120) {
  sub <- lp[lead_set == lead_set_name]
  B0  <- as.matrix(sub[, ..bc])
  s0  <- fread(file.path(.repo, "results", "01_data", "lead_profiles", sigma_file))[chr == 1]
  S0  <- as.matrix(s0[, -c(1, 2), with = FALSE])
  se0_diag <- sqrt(diag(S0))
  q   <- bq[trait == trait_name][order(tau)]$q
  top_idx <- order(-abs(sub$delta_mu_z))[seq_len(min(n_bg, nrow(sub)))]

  bg_list <- list(); sum_list <- list(); hi_list <- list()
  for (nm in names(transforms)) {
    dvec <- transforms[[nm]](q)
    Bt <- sweep(B0, 2, dvec, "*")
    mu <- rowMeans(Bt)
    St <- Bt / mu
    # Background leads (top n_bg by precision)
    M_bg <- St[top_idx, , drop = FALSE]
    bg_list[[nm]] <- data.table(
      scale = nm,
      snp_id = rep(sub$snp_id[top_idx], each = length(tau)),
      tau = rep(tau, nrow(M_bg)),
      val = as.vector(t(M_bg))
    )
    # Pointwise median and IQR across all leads
    sum_list[[nm]] <- data.table(
      scale = nm,
      tau = tau,
      med = apply(St, 2, median),
      q25 = apply(St, 2, quantile, probs = 0.25),
      q75 = apply(St, 2, quantile, probs = 0.75)
    )
    # Highlighted high-precision leads with 95% CI
    if (length(hi_snps)) {
      for (sid in hi_snps) {
        idx <- which(sub$snp_id == sid)
        stopifnot(length(idx) == 1)
        b_i  <- Bt[idx, ]
        se_i <- dvec * se0_diag / sqrt(sub$d_resid[idx])
        mu_i <- mean(b_i)
        hi_list[[paste(nm, sid)]] <- data.table(
          scale = nm, snp_id = sid, tau = tau,
          val = b_i / mu_i,
          lo  = (b_i - 1.96 * se_i) / mu_i,
          hi  = (b_i + 1.96 * se_i) / mu_i
        )
      }
    }
  }
  list(
    bg  = rbindlist(bg_list),
    sum = rbindlist(sum_list),
    hi  = if (length(hi_list)) rbindlist(hi_list) else data.table()
  )
}

# ---- a: 1D trait (IGF-1) across log, untransformed, Box-Cox lambda = 0.57 ----
igf_tr <- list(
  '"Log"'                            = function(q) rep(1, length(q)),
  '"Untransformed"'                  = function(q) exp(q),
  '"Box-Cox,"~lambda == 0.57'        = function(q) exp(0.57 * q)
)
igf_f <- make_fan("log_igf1_A_deltamu", "log_igf1", "sigma0_log_igf1.csv", igf_tr, n_bg = 120)
igf_lvls <- names(igf_tr)
igf_f$bg[,  scale_f := factor(scale, levels = igf_lvls)]
igf_f$sum[, scale_f := factor(scale, levels = igf_lvls)]

pA <- ggplot() +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.35) +
  geom_line(data = igf_f$bg, aes(tau, val, group = snp_id), colour = BL, alpha = 0.08, linewidth = 0.22) +
  geom_ribbon(data = igf_f$sum, aes(tau, ymin = q25, ymax = q75), fill = BL, alpha = 0.22) +
  geom_line(data = smooth_curve(igf_f$sum, "tau", "med", by = "scale_f"),
            aes(x, y), colour = BL, linewidth = 0.9, lineend = "round") +
  facet_wrap(~scale_f, nrow = 1, labeller = label_parsed) +
  coord_cartesian(ylim = c(0.15, 2.05)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  labs(x = expression("Quantile level " * tau), y = "IGF-1 effect / mean") +
  theme_paper() + labs(tag = "a")

# ---- b: 2D trait (GGT) across untransformed, log, Box-Cox lambda = -0.67 ----
ggt_tr <- list(
  '"Untransformed"'                  = function(q) exp(q),
  '"Log"'                            = function(q) rep(1, length(q)),
  '"Box-Cox,"~lambda == -0.67'       = function(q) exp(-0.67 * q)
)
ggt_snps <- c("rs2650000", "rs12373325")
ggt_f <- make_fan("log_ggt_A_deltamu", "log_ggt", "sigma0_log_ggt.csv", ggt_tr, hi_snps = ggt_snps, n_bg = 120)
ggt_lvls <- names(ggt_tr)
ggt_f$bg[,  scale_f := factor(scale, levels = ggt_lvls)]
ggt_f$sum[, scale_f := factor(scale, levels = ggt_lvls)]
ggt_f$hi[,  scale_f := factor(scale, levels = ggt_lvls)]
ggt_f$hi[,  snp_f   := factor(snp_id, levels = ggt_snps)]
hi_col <- c(rs2650000 = BL, rs12373325 = CO)

pB <- ggplot() +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.35) +
  geom_line(data = ggt_f$bg, aes(tau, val, group = snp_id), colour = "grey55", alpha = 0.08, linewidth = 0.22) +
  geom_line(data = smooth_curve(ggt_f$sum, "tau", "med", by = "scale_f"),
            aes(x, y), colour = "grey25", linewidth = 0.65, linetype = "32") +
  geom_ribbon(data = ggt_f$hi, aes(tau, ymin = lo, ymax = hi, fill = snp_f), alpha = 0.20, colour = NA) +
  geom_line(data = smooth_curve(ggt_f$hi, "tau", "val", by = c("scale_f", "snp_f")),
            aes(x, y, colour = snp_f), linewidth = 0.9, lineend = "round") +
  facet_wrap(~scale_f, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = hi_col, name = NULL) +
  scale_fill_manual(values = hi_col, guide = "none") +
  coord_cartesian(ylim = c(0.0, 2.7)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  labs(x = expression("Quantile level " * tau), y = "GGT effect / mean") +
  theme_paper() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.96),
        legend.justification = c(0, 1), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.65, "lines")) +
  labs(tag = "b")

# ---- c: Scale-invariant whitened PC spectra across 3 transformations ----
calc_spec <- function(lead_set_name, trait_name, sigma_file, lam, label) {
  sub <- lp[lead_set == lead_set_name]
  B0  <- as.matrix(sub[, ..bc])
  s0  <- fread(file.path(.repo, "results", "01_data", "lead_profiles", sigma_file))[chr == 1]
  S0  <- as.matrix(s0[, -c(1, 2), with = FALSE])
  q   <- bq[trait == trait_name][order(tau)]$q
  is_log <- startsWith(trait_name, "log_")
  Ds <- if (is_log) {
    list(Untransformed = exp(q), Log = rep(1, 45), `Box-Cox` = exp(lam * q))
  } else {
    list(Untransformed = rep(1, 45), Log = 1 / q, `Box-Cox` = q^(lam - 1))
  }
  rbindlist(lapply(names(Ds), function(sc) {
    dvec <- Ds[[sc]]
    Bt <- sweep(B0, 2, dvec, "*")
    St <- (dvec %o% dvec) * S0
    ee <- eigen(St, symmetric = TRUE)
    Wih <- ee$vectors %*% diag(1 / sqrt(ee$values)) %*% t(ee$vectors)
    Y <- sweep(Bt %*% Wih, 1, sqrt(sub$d_resid), "*")
    U <- Y / sqrt(rowSums(Y^2))
    sv <- svd(U)
    pct <- 100 * sv$d^2 / sum(sv$d^2)
    data.table(trait = label, arm = trait_name, scale = sc, pc = 1:6, pct = pct[1:6])
  }))
}

spec_dt <- rbindlist(list(
  calc_spec("log_igf1_A_deltamu",     "log_igf1",     "sigma0_log_igf1.csv",     0.57,  "IGF-1 (1D)"),
  calc_spec("bmi_A_regenie",          "bmi",          "sigma0_bmi.csv",         -1.51,  "BMI (2D)"),
  calc_spec("log_ggt_A_deltamu",      "log_ggt",      "sigma0_log_ggt.csv",     -0.67,  "GGT (2D)"),
  calc_spec("log_hba1c_A_deltamu",    "log_hba1c",    "sigma0_log_hba1c.csv",   -1.41,  "HbA1c (3D)"),
  calc_spec("hba1c_nodiab_A_deltamu", "hba1c_nodiab", "sigma0_hba1c_nodiab.csv", 0.73,  "HbA1c, no diabetes (2D)")
))
# Verify invariance to floating-point precision
stopifnot(spec_dt[, max(pct) - min(pct), by = .(trait, pc)][, max(V1)] < 1e-10)

pa_sub <- pa[arm %in% c("log_igf1", "bmi", "log_ggt", "log_hba1c", "hba1c_nodiab") & component_tested <= 5,
             .(arm, pc = component_tested, null_med = null_median_pct, null_p95 = null_p95_pct)]
spec_r <- merge(spec_dt[pc >= 2 & pc <= 5], pa_sub, by = c("arm", "pc"))
spec_r[, ratio := pct / null_med]
spec_r[, thr   := null_p95 / null_med]

tr_lvls <- c("IGF-1 (1D)", "BMI (2D)", "GGT (2D)", "HbA1c, no diabetes (2D)", "HbA1c (3D)")
spec_r[, trait_f := factor(trait, levels = tr_lvls)]
spec_r[, scale_f := factor(scale, levels = c("Untransformed", "Log", "Box-Cox"))]
tr_col <- c("IGF-1 (1D)" = BL, "BMI (2D)" = "grey35", "GGT (2D)" = TE,
            "HbA1c, no diabetes (2D)" = "#d98b2b", "HbA1c (3D)" = CO)
dodge_map <- setNames(seq(-0.24, 0.24, length.out = length(tr_lvls)), tr_lvls)
spec_r[, x_pos := pc + dodge_map[as.character(trait_f)]]
thr_dt <- spec_r[scale == "Log"]

pC <- ggplot(spec_r, aes(x_pos, ratio, colour = trait_f)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_segment(data = thr_dt, aes(x = x_pos - 0.055, xend = x_pos + 0.055, y = thr, yend = thr),
               colour = "grey30", linewidth = 0.45) +
  geom_segment(data = thr_dt, aes(x = x_pos, xend = x_pos, y = 1, yend = ratio), linewidth = 0.45, alpha = 0.65) +
  geom_point(aes(shape = scale_f), size = 1.85, stroke = 0.55) +
  scale_colour_manual(values = tr_col, name = NULL) +
  scale_shape_manual(values = c(Untransformed = 16, Log = 2, `Box-Cox` = 5), name = NULL) +
  scale_x_continuous(breaks = 2:5, labels = paste0("PC", 2:5)) +
  scale_y_log10(breaks = c(1, 1.5, 2, 3, 4), labels = c("1", "1.5", "2", "3", "4")) +
  coord_cartesian(ylim = c(0.85, 4.35)) +
  labs(x = "Principal component", y = "Variance share / null median") +
  theme_paper() +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.62, "lines"),
        legend.box = "horizontal", legend.spacing.x = unit(4, "pt")) +
  labs(tag = "c")

# ---- d: Scale-invariant whitened PC1 vs PC2 in HbA1c separates diabetes loci ----
hb_d <- merge(
  lm[arm == "log_hba1c", .(snp_id, pc1 = pc1_score, pc2 = pc2_score, pc3 = pc3_score, z_mu = delta_mu_z)],
  diab[, .(snp_id, diab_z, diab_p, skew_hit, nearest_gene)],
  by = "snp_id"
)
if (median(hb_d$pc1) < 0) hb_d[, pc1 := -pc1]
if (cor(hb_d$pc2, hb_d$diab_z, use = "complete.obs") < 0) hb_d[, pc2 := -pc2]
hb_d[, diab_sig := factor(fifelse(diab_p < 0.05 / nrow(hb_d), "Diabetes-associated", "Other"),
                          levels = c("Other", "Diabetes-associated"))]
lab_snps <- data.table(
  snp_id = c("rs7903146", "rs2071479", "rs7756992", "rs1046917"),
  gene   = c("TCF7L2", "HLA-DOB", "CDKAL1", "FN3KRP"),
  lx     = c(0.81, 0.13, 0.91, 0.92),
  ly     = c(0.65, 0.72, 0.50, -0.33),
  hj     = c(0,    0,    0,    0.5)
)
lab_dt <- merge(lab_snps, hb_d, by = "snp_id")

pD <- ggplot(hb_d, aes(pc1, pc2)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70", linewidth = 0.3) +
  geom_point(data = hb_d[diab_sig == "Other"], aes(colour = diab_sig), size = 1.1, alpha = 0.55, stroke = 0) +
  geom_point(data = hb_d[diab_sig == "Diabetes-associated"], aes(colour = diab_sig), size = 1.5, alpha = 0.85, stroke = 0) +
  geom_point(data = lab_dt, shape = 21, fill = NA, colour = "grey15", size = 2.5, stroke = 0.55) +
  geom_segment(data = lab_dt[gene %in% c("TCF7L2", "CDKAL1", "FN3KRP")],
               aes(x = pc1, y = pc2, xend = lx, yend = ly),
               colour = "grey45", linewidth = 0.25) +
  geom_text(data = lab_dt, aes(x = lx, y = ly, label = gene, hjust = hj),
            fontface = "italic", colour = "grey15", size = 2.15,
            vjust = ifelse(lab_dt$gene == "FN3KRP", 1.15, 0.35)) +
  scale_colour_manual(values = c("Other" = BL, "Diabetes-associated" = CO), name = NULL) +
  coord_cartesian(xlim = c(0.03, 1.07), ylim = c(-0.52, 0.78)) +
  labs(x = "Whitened PC1 score", y = "Whitened PC2 score") +
  theme_paper() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.03),
        legend.justification = c(0, 0), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.65, "lines")) +
  labs(tag = "d")

fig <- pA / pB / ((pC | pD) + plot_layout(widths = c(1.15, 1))) +
  plot_layout(heights = c(1, 1.05, 1.25))

src <- rbindlist(list(
  igf_f$sum[, .(panel = "a", trait = "IGF-1", scale, x = tau, y = med, lo = q25, hi = q75)],
  ggt_f$hi[,  .(panel = "b", trait = "GGT",   scale, series = snp_id, x = tau, y = val, lo, hi)],
  spec_r[,    .(panel = "c", trait,           scale, x = pc, y = ratio, threshold = thr)],
  hb_d[,      .(panel = "d", trait = "HbA1c", series = snp_id, x = pc1, y = pc2, group = as.character(diab_sig))]
), fill = TRUE)

cap <- paste(
  "Figure 5. Scale-invariant dimensionality separates scale-reducible, heterogeneous and mixture traits.",
  "UK Biobank lead variants; quantile profiles are converted across transformations directly from Stage 1 output (Methods).",
  "**a**, One-dimensional trait: IGF-1 (644 lead variants; n = 438,938) on log, untransformed and Box-Cox (lambda = 0.57)",
  "scales, each profile divided by its mean effect. Thin lines show the 120 most precise leads; heavy line and band show",
  "the pointwise median and interquartile range across all 644 leads; dashed line, pure shift. All leads rotate together",
  "onto a flat shift under one transformation.",
  "**b**, Two-dimensional trait: GGT (454 lead variants; n = 441,089) on untransformed, log and Box-Cox (lambda = -0.67)",
  "scales. Grey lines show the 120 most precise leads and dashed grey line the pointwise median across all leads;",
  "coloured curves and bands show two high-precision leads (rs2650000 and rs12373325, 95% confidence intervals).",
  "Log transformation flattens rs2650000 while leaving rs12373325 steeply tilted; the trait-flattening Box-Cox",
  "transformation bends the two variants in opposite directions, showing that no monotonic transformation flattens both.",
  "**c**, Ratio of observed whitened principal-component variance share to its parallel-analysis null median for PC2 to PC5",
  "(ticks, null 95th percentile; dashed line, 1) across untransformed, log and Box-Cox transformations for five traits.",
  "Symbols for all three transformations coincide identically (maximum difference < 10^-10), demonstrating exact scale",
  "invariance and distinguishing one-dimensional (IGF-1), two-dimensional (BMI, GGT, HbA1c without diabetes) and",
  "three-dimensional (HbA1c) genetic architectures.",
  "**d**, Scale-invariant whitened PC1 and PC2 scores for the 472 log HbA1c leads (n = 441,144), coloured by association",
  "with diabetes (logistic regression in the same participants, two-sided Bonferroni P < 0.05/472; Pearson r = 0.81",
  "between PC2 and diabetes z statistic). Exemplar diabetes loci (TCF7L2, HLA-DOB, CDKAL1) separate cleanly from a",
  "mean-matched shift locus (FN3KRP) purely from quantile-profile geometry. GGT, gamma-glutamyltransferase; IGF-1,",
  "insulin-like growth factor 1; PC, principal component."
)

save_figure(fig, "fig5_scale_invariance", src, cap, width = 183, height = 165,
            results_used = c("01_data/baseline_quantiles", "01_data/lead_profiles",
                             "06_traits/moments", "06_traits/moments",
                             "07_hba1c/skew_diabetes"))
cat("Saved fig5_scale_invariance cleanly!\n")
