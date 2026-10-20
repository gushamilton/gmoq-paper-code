#!/usr/bin/env Rscript
# BMI variant effects grow with polygenic background: summaries across the 639 leads
#
# What: from the per-lead interaction estimates, (1) estimates the common relative interaction
#   kappa = delta / beta across leads (each lead's interaction delta against its main effect
#   beta, with the measurement-error-corrected heterogeneity model) on four scales of BMI;
#   (2) contrasts the FTO effect in the top and bottom score deciles; (3) runs the cross-fitted
#   link test: does a lead with an unusual interaction in one half of the sample also have an
#   unusual tilt, relative to the score's own quantile profile, in the other half? and
#   (4) compares the leads' common tilt with the tilt predicted if each allele's effect grew with
#   a modifier correlated with BMI: BMI itself, genetic background at the LDSC SNP heritability,
#   genetic background at twin heritability, or the measured score. Bootstrap intervals
#   resample leads (2,000 resamples).
# Paper: Fig. 6; Supplementary Fig. S14; Supplementary Table S6
# Input: per-lead regressions, FTO deciles, outcome summary and the score's quantile profile
#   (results/05_bmi/prs_interaction, from 09_prs_interaction.R); lead profiles and pooled
#   covariance, LOCO none (results/01_data/lead_profiles); Stage 1 quantiles
#   (results/01_data/baseline_quantiles); LDSC SNP heritability of log BMI
#   (results/05_bmi/genomewide_models/h2.csv)
# Output: results/05_bmi/prs_interaction/{kappa_by_scale,fto_prs_contrast,link_test,
#   fig5d_summary,fig5d_curves,fig5d_inputs,per_lead_prs}.csv
# Runs: locally from results/ in about a minute

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
source(file.path(repo, "analysis", "05_bmi", "bmi_common.R"))
cfg <- read_config()
P <- cfg$prs
u <- tau_grid(cfg); z <- qnorm(u); n_tau <- length(u)
N_BOOT <- as.integer(P$boot); SEED <- as.integer(cfg$estimation$boot_seed)
OUT <- "05_bmi/prs_interaction"

## 1. Inputs -----------------------------------------------------------------------------
R <- as.data.table(read_result(OUT, "per_lead_regressions"))
fto_deciles <- as.data.table(read_result(OUT, "fto_prs_deciles"))
outcomes <- as.data.table(read_result(OUT, "outcome_summary"))
prs_profile <- as.data.table(read_result(OUT, "prs_quantile_profile"))
baseline_q <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))
baseline_summary <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_summary"))
ldsc_h2 <- as.data.table(read_result("05_bmi/genomewide_models", "h2"))[trait == "regenie__log_bmi"]

# Lead profiles without a leave-one-chromosome-out adjustment (that prediction is itself a
# polygenic score), ordered by position and aligned to the BMI-increasing allele.
profiles <- read_lead_profiles()
log_set <- read_lead_set(profiles, P$lead_set, u, order_by_position = TRUE)
ref <- log_set$meta[, .(snp_id, chr, bp, effect_allele, other_allele, d_resid)]
stopifnot(identical(R$snp_id, ref$snp_id))
n_leads <- nrow(ref)
raw_set <- read_lead_set(profiles, "bmi_prs_none", u, order_by_position = TRUE)
leads <- list(log_bmi = align_to_reference(ref, log_set$B, ref),
              bmi = align_to_reference(raw_set$meta, raw_set$B, ref))
Sigma0 <- list(log_bmi = read_lead_sigma0("log_bmi_prs_none", u), bmi = read_lead_sigma0("bmi_prs_none", u))
prs_profile_matrix <- function(trait) {
  x <- prs_profile[scale == trait][order(left_out_chr, tau)]
  matrix(x$b_prs, nrow = 22, byrow = TRUE)
}

## 2. Common relative interaction kappa, by outcome scale ---------------------------------------
# Heterogeneity model with L = beta_i and S = delta_i: S_i = kappa_i L_i + e_i,
# kappa_i ~ (kappa, tau^2), using each lead's (co)variances of (beta, delta).
kappa_rows <- list(); kappa_fit <- list()
for (yc in c("raw", "log", "int", "natural")) for (se in c("h", "m")) {
  L <- R[[paste0(yc, "_b")]]; S <- R[[paste0(yc, "_g")]]
  vL <- R[[paste0(yc, "_vb_", se)]]; vS <- R[[paste0(yc, "_vg_", se)]]; cLS <- R[[paste0(yc, "_cbg_", se)]]
  ht <- het_test(L, S, vL, vS, cLS); bt <- het_boot(L, S, vL, vS, cLS, B = N_BOOT, seed = SEED)
  kappa_fit[[paste(yc, se)]] <- ht
  kappa_rows[[length(kappa_rows) + 1L]] <- data.table(
    scale = yc, se_type = c(h = "HC0", m = "model")[[se]], n_leads = ht$n,
    kappa = ht$r, kappa_se_boot = bt$r_se, kappa_lo = bt$r_ci[1], kappa_hi = bt$r_ci[2],
    Q = ht$Q, df = ht$df, p_het = ht$p, tau = ht$tau, tau_lo = bt$tau_ci[1], tau_hi = bt$tau_ci[2],
    tau_excess = ht$excess, mean_z_int = mean(S / sqrt(vS)), n_gamma_pos = sum(S > 0),
    median_hc0_over_model_se_gamma = median(sqrt(R[[paste0(yc, "_vg_h")]] / R[[paste0(yc, "_vg_m")]])),
    lambda_mean_b_prs = mean(R[[paste0(yc, "_b_prs")]]), var_y = outcomes[scale == yc]$variance)
}
kappa_by_scale <- rbindlist(kappa_rows)
kappa_by_scale[, kappa_over_kappa_raw := kappa / kappa_by_scale[scale == "raw" & se_type == "HC0", kappa]]

## 3. Lead profiles against the score's own quantile profile ------------------------------------
# Full sample: alpha (level) and rho (extra tilt) per lead and their heterogeneity fit; and the
# fit of each profile to the line [1, q - qbar] on this scale's own quantiles, whose common
# slope/level ratio r_obs is the observed common tilt used in step 6.
tilt <- list()
for (tr in c("log_bmi", "bmi")) {
  x <- leads[[tr]]
  tt <- prs_referenced_tilt(x$B, x$d, ref$chr, Sigma0[[tr]], prs_profile_matrix(tr), u)
  q <- baseline_q[trait == tr][order(tau)]$q; qbar <- mean(q)
  line <- project_profiles(x$B, cbind(1, q - qbar), lapply(as.character(ref$chr), function(k) Sigma0[[tr]][[k]]), x$d)
  tilt[[tr]] <- list(tt = tt, ht = het_test(tt$alpha, tt$rho, tt$vaa, tt$vrr, tt$var), line = line, q = q, qbar = qbar)
}

## 4. Cross-fitted link test ------------------------------------------------------------------
# Residuals from the common interaction (per half) and from the common PRS-referenced tilt (per
# half); each half's interaction residuals are correlated with the other half's tilt residuals,
# and the two cross-fits averaged. Full-sample correlations share estimation error and are
# shown for reference only.
link_rows <- list(); link_resid <- list()
for (tr in c("log_bmi", "bmi")) {
  yc <- c(log_bmi = "log", bmi = "raw")[[tr]]
  r_int <- lapply(c("h1", "h2"), function(hh) het_test(R[[paste0(hh, "_", yc, "_b")]], R[[paste0(hh, "_", yc, "_g")]],
    R[[paste0(hh, "_", yc, "_vb_h")]], R[[paste0(hh, "_", yc, "_vg_h")]], R[[paste0(hh, "_", yc, "_cbg_h")]])$resid_z)
  r_tilt <- lapply(c("h1", "h2"), function(hh) { col <- function(k) R[[paste("tilt", tr, hh, k, sep = "_")]]
    het_test(col("alpha"), col("rho"), col("vaa"), col("vrr"), col("var"))$resid_z })
  cross <- function(ix) mean(c(cor(r_int[[1]][ix], r_tilt[[2]][ix]), cor(r_int[[2]][ix], r_tilt[[1]][ix])))
  set.seed(SEED)
  boot <- replicate(N_BOOT, cross(sample.int(n_leads, n_leads, TRUE)))
  estimate <- cross(seq_len(n_leads)); se_boot <- sd(boot)
  r_int_full <- kappa_fit[[paste(yc, "h")]]$resid_z; r_tilt_full <- tilt[[tr]]$ht$resid_z
  link_resid[[tr]] <- list(r_int = r_int, r_tilt = r_tilt, r_int_full = r_int_full, r_tilt_full = r_tilt_full)
  link_rows[[tr]] <- data.table(scale = tr,
    test = c("PRIMARY cross-fitted: mean(cor(int_h1, tilt_h2), cor(int_h2, tilt_h1))",
             "cross-fit component int_h1 x tilt_h2", "cross-fit component int_h2 x tilt_h1",
             "full-sample int x PRS-referenced tilt (non-inferential: shared estimation error)"),
    r = c(estimate, cor(r_int[[1]], r_tilt[[2]]), cor(r_int[[2]], r_tilt[[1]]), cor(r_int_full, r_tilt_full)),
    se_boot = c(se_boot, NA, NA, NA), lo = c(quantile(boot, 0.025), NA, NA, NA), hi = c(quantile(boot, 0.975), NA, NA, NA),
    n_leads = n_leads)
  link_rows[[tr]][1, `:=`(z = r / se_boot, p = 2 * pnorm(-abs(r / se_boot)))]
}

## 5. FTO: top versus bottom decile of the score ---------------------------------------------------
fto <- R[snp_id == P$fto_snp]
fto_prs_contrast <- rbindlist(lapply(c("raw", "log"), function(yc) {
  top <- fto_deciles[scale == yc & prs_decile == 10]; bottom <- fto_deciles[scale == yc & prs_decile == 1]
  log_ratio <- log(top$beta / bottom$beta)
  se_log_ratio <- sqrt((top$se_hc0 / top$beta)^2 + (bottom$se_hc0 / bottom$beta)^2)
  b <- fto[[paste0(yc, "_b")]]; g <- fto[[paste0(yc, "_g")]]
  kappa <- kappa_by_scale[scale == yc & se_type == "HC0", kappa]
  diff_z <- (top$beta - bottom$beta) / sqrt(top$se_hc0^2 + bottom$se_hc0^2)
  data.table(snp_id = P$fto_snp, scale = yc, beta_top10 = top$beta, se_top10 = top$se_hc0, n_top10 = top$n,
             mean_prs_top10 = top$mean_prs, beta_bottom10 = bottom$beta, se_bottom10 = bottom$se_hc0, n_bottom10 = bottom$n,
             mean_prs_bottom10 = bottom$mean_prs, ratio_top_over_bottom = top$beta / bottom$beta,
             ratio_lo = exp(log_ratio - 1.96 * se_log_ratio), ratio_hi = exp(log_ratio + 1.96 * se_log_ratio),
             diff_z = diff_z, diff_p = 2 * pnorm(-abs(diff_z)), fto_beta_full = b, fto_gamma_full = g,
             predicted_ratio_own_interaction = (b + g * top$mean_prs) / (b + g * bottom$mean_prs),
             predicted_ratio_common_kappa = (1 + kappa * top$mean_prs) / (1 + kappa * bottom$mean_prs))
}))

## 6. Common tilt predicted by candidate modifiers (log BMI) ----------------------------------
# If each allele's effect grows with a modifier M by kappa per SD of the score, and the score
# has main effect lambda, an allele's quantile profile tilts as 1 + omega k (q - mu), where
# k = kappa / lambda and omega is the share of BMI variance the modifier carries: 1 for BMI
# itself, h2_SNP or twin h2 for genetic background, R2_S = lambda^2 / Var(log BMI) for the
# measured score. Projected on [1, q - qbar] the predicted slope/level ratio is
# r_H = omega k / (1 + omega k (qbar - mu)); F = r_H / r_obs is the fraction of the observed
# common tilt accounted for. Leads are resampled jointly for kappa and r_obs.
q <- tilt$log_bmi$q; qbar <- tilt$log_bmi$qbar
mu <- baseline_summary[trait == "log_bmi"]$mean
kappa <- kappa_fit[["log h"]]$r
lambda <- mean(R$log_b_prs); var_y <- outcomes[scale == "log"]$variance; R2_S <- lambda^2 / var_y
k <- kappa / lambda
line <- tilt$log_bmi$line
a_line <- line$theta[, 1]; s_line <- line$theta[, 2]; v_aa <- line$V[, 1, 1]; v_ss <- line$V[, 2, 2]; v_as <- line$V[, 1, 2]
r_obs <- het_test(a_line, s_line, v_aa, v_ss, v_as)$r
h2 <- ldsc_h2$h2; h2_se <- ldsc_h2$h2_se
hypotheses <- data.table(hypothesis = c("BMI (own phenotype) as modifier", "genetic background only (h2_SNP, LDSC)",
                                        "genetic background, twin h2 (sensitivity)", "measured PRS only (R2_S)"),
  hyp_key = c("bmi", "h2snp", "twin", "prs"),
  omega = c(1, h2, cfg$headline$fig5d_twin_h2, R2_S),
  omega_source = c("E[M | Y] = Y - mu", sprintf("LDSC h2_SNP %.4f (SE %.4f), REGENIE log BMI N 461,320, w_hm3, 1000G EUR", h2, h2_se),
                   "config headline.fig5d_twin_h2", "lambda^2 / Var(log BMI), in-sample"))
r_H <- function(omega, kk) omega * kk / (1 + omega * kk * (qbar - mu))
set.seed(SEED)
boot <- t(replicate(N_BOOT, {
  i <- sample.int(n_leads, n_leads, TRUE)
  c(kappa = het_test(R$log_b[i], R$log_g[i], R$log_vb_h[i], R$log_vg_h[i], R$log_cbg_h[i])$r,
    r_obs = het_test(a_line[i], s_line[i], v_aa[i], v_ss[i], v_as[i])$r)
}))
fig5d_summary <- rbindlist(lapply(seq_len(nrow(hypotheses)), function(j) {
  om <- hypotheses$omega[j]
  r_boot <- r_H(om, boot[, "kappa"] / lambda); F_boot <- r_boot / boot[, "r_obs"]
  data.table(hypotheses[j], k = k, lambda = lambda, kappa = kappa,
             r_H = r_H(om, k), r_H_lo = quantile(r_boot, 0.025), r_H_hi = quantile(r_boot, 0.975),
             r_obs = r_obs, r_obs_lo = quantile(boot[, "r_obs"], 0.025), r_obs_hi = quantile(boot[, "r_obs"], 0.975),
             F = r_H(om, k) / r_obs, F_lo = quantile(F_boot, 0.025), F_hi = quantile(F_boot, 0.975), n_boot = N_BOOT)
}))
band <- sapply(seq_len(n_tau), function(t) quantile(1 + boot[, "r_obs"] * (q[t] - qbar), c(0.025, 0.975)))
fig5d_curves <- data.table(tau = u, q = q, q_minus_mu = q - mu, observed = 1 + r_obs * (q - qbar),
                           observed_lo = band[1, ], observed_hi = band[2, ])
for (j in seq_len(nrow(hypotheses))) set(fig5d_curves, j = paste0("pred_", hypotheses$hyp_key[j]),
                                         value = 1 + hypotheses$omega[j] * k * (q - mu))
fig5d_inputs <- data.table(
  item = c("kappa_log_HC0", "kappa_log_se_boot", "lambda_mean_b_prs_log", "lambda_sd_over_leads", "var_log_bmi",
           "R2_S = lambda^2/Var(Y)", "k = kappa/lambda", "h2_SNP_log_bmi_LDSC", "h2_SNP_log_bmi_LDSC_se",
           "twin_h2_sensitivity", "mu_log_bmi_P03", "qbar_mean_q_tau", "n", "n_leads", "lead_profiles", "q_source"),
  value = c(kappa, kappa_by_scale[scale == "log" & se_type == "HC0", kappa_se_boot], lambda, sd(R$log_b_prs), var_y,
            R2_S, k, h2, h2_se, cfg$headline$fig5d_twin_h2, mu, qbar, outcomes[scale == "log"]$n, n_leads,
            "lead profiles log_bmi_prs_none (LOCO none), Sigma0 per chromosome",
            "Stage 1 log BMI quantiles (type 8, exact levels) on the 461,320 sample"))

## 7. Per-lead table -----------------------------------------------------------------------
z_none <- matrix(NA_real_, n_leads, 3)                     # unit-SE level/tilt/bend components
for (cc in unique(ref$chr)) {
  ii <- which(ref$chr == cc)
  z_none[ii, ] <- unit_se_components(leads$log_bmi$B[ii, , drop = FALSE], Sigma0$log_bmi[[as.character(cc)]], leads$log_bmi$d[ii], u)$Z
}
tl <- tilt$log_bmi$tt; lk <- link_resid$log_bmi
per_lead_prs <- data.table(ref[, .(snp_id, chr, bp, effect_allele, other_allele)],
  beta_log = R$log_b, se_beta_log_hc0 = sqrt(R$log_vb_h), gamma_log = R$log_g, se_gamma_log_hc0 = sqrt(R$log_vg_h),
  se_gamma_log_model = sqrt(R$log_vg_m), b_prs_log = R$log_b_prs,
  beta_raw = R$raw_b, gamma_raw = R$raw_g, se_gamma_raw_hc0 = sqrt(R$raw_vg_h),
  gamma_int = R$int_g, se_gamma_int_hc0 = sqrt(R$int_vg_h), gamma_natural = R$natural_g, se_gamma_natural_hc0 = sqrt(R$natural_vg_h),
  alpha_rel_prs_log = tl$alpha, rho_rel_prs_log = tl$rho, se_rho_rel_prs_log = sqrt(tl$vrr),
  resid_z_int_log_full = lk$r_int_full, resid_z_tilt_log_full = lk$r_tilt_full,
  resid_z_int_log_h1 = lk$r_int[[1]], resid_z_int_log_h2 = lk$r_int[[2]],
  resid_z_tilt_log_h1 = lk$r_tilt[[1]], resid_z_tilt_log_h2 = lk$r_tilt[[2]],
  line_a_log = a_line, line_s_log = s_line, line_se_s_log = sqrt(v_ss), line_s_over_a_log = s_line / a_line,
  z_c1_none = z_none[, 1], z_c2_none = z_none[, 2])

write_result(kappa_by_scale, OUT, "kappa_by_scale")
write_result(fto_prs_contrast, OUT, "fto_prs_contrast")
write_result(rbindlist(link_rows, fill = TRUE), OUT, "link_test")
write_result(fig5d_summary, OUT, "fig5d_summary")
write_result(fig5d_curves, OUT, "fig5d_curves")
write_result(fig5d_inputs, OUT, "fig5d_inputs")
write_result(per_lead_prs, OUT, "per_lead_prs")
