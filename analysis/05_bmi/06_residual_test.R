#!/usr/bin/env Rscript
# Is there anything in BMI lead-variant profiles beyond a shift and a spread?
#
# What: for each of the 639 BMI lead variants, removes the best-fitting shift + spread
#   (a constant and a linear term in the normal quantile z) from its Stage 1 quantile profile
#   and tests what is left, Q = d r' Sigma^-1 r, against its exact null distribution (a weighted
#   sum of chi-square(1) variables; Davies' method). Profiles are first converted to normal
#   scores, beta_R(u) = beta(u) fhat(u) / phi(z_u), which makes the test invariant to any
#   monotone rescaling of BMI; raw and log scales are reported as secondary. The test is
#   calibrated on 20,001 null variants and on 20,000 simulated pure shift + spread profiles.
#   A sensitivity analysis removes the top two principal components of the whitened profiles
#   instead, estimated leaving out the lead's chromosome.
# Paper: Methods (residual test of BMI profiles)
# Input: lead profiles and pooled covariance (results/01_data/lead_profiles); Stage 1
#   quantiles and density (results/01_data/baseline_quantiles); for the null calibration, the
#   genome-wide Stage 1 quantile GWAS of BMI (one file per chromosome) and the REGENIE BMI GWAS
#   (allele frequencies)
# Output: results/05_bmi/residual_test/{lead_residual,pc_sensitivity,top10_residual}.csv
#   (from results/), and {calibration,qq,summary,null_counts_by_chr}.csv (needs the genome-wide
#   Stage 1)
# Runs: leads and simulation locally from results/ in about a minute; the null calibration on
#   a cluster with UK Biobank data

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
source(file.path(repo, "analysis", "05_bmi", "bmi_common.R"))
cfg <- read_config()
u <- tau_grid(cfg); n_tau <- length(u); z <- qnorm(u); trap_w <- trap_weights(u)
SEED <- 20260923L; N_NULL <- 20000L; N_SIM <- 20000L
OUT <- "05_bmi/residual_test"
set.seed(SEED)

## 1. Inputs ---------------------------------------------------------------------------
profiles <- read_lead_profiles()
baseline_q <- read_result("01_data/baseline_quantiles", "baseline_quantiles")
leads_raw <- read_lead_set(profiles, "bmi_A_regenie", u)
leads_log <- read_lead_set(profiles, "log_bmi_A_regenie", u)
Sigma0_raw <- read_lead_sigma0("bmi", u); Sigma0_log <- read_lead_sigma0("log_bmi", u)
q_raw <- stage1_quantiles(baseline_q, "bmi", u); q_log <- stage1_quantiles(baseline_q, "log_bmi", u)
to_normal_raw <- q_raw$fhat / dnorm(z)                       # d(normal score) / d(BMI) at Q(u)
to_normal_log <- q_log$fhat / dnorm(z)                       # d(normal score) / d(log BMI)
mu_raw <- sum(trap_w * q_raw$q) / sum(trap_w); mu_log <- sum(trap_w * q_log$q) / sum(trap_w)

# Each scale: the conversion factor k(u) (beta -> k beta, Sigma0 -> D Sigma0 D), the
# shift + spread frame X, and the Stage 1 run whose profiles are used.
SCALES <- list(
  normal_score = list(k = to_normal_raw, X = cbind(1, z), Sigma0 = Sigma0_raw, src = "raw"),
  raw = list(k = rep(1, n_tau), X = cbind(1, q_raw$q - mu_raw), Sigma0 = Sigma0_raw, src = "raw"),
  log = list(k = rep(1, n_tau), X = cbind(1, q_log$q - mu_log), Sigma0 = Sigma0_log, src = "log"),
  normal_score_from_log = list(k = to_normal_log, X = cbind(1, z), Sigma0 = Sigma0_log, src = "log"))
# Null variants come from the raw-BMI Stage 1; for the log scale they are converted exactly.
SCALES_NULL <- list(normal_score = SCALES$normal_score, raw = SCALES$raw,
                    log = list(k = 1 / q_raw$q, X = cbind(1, q_log$q - mu_log), Sigma0 = Sigma0_raw, src = "raw"))
scaled <- function(S, k) change_scale(matrix(0, 1, n_tau), S, k)$Sigma0
spectrum_cache <- new.env()
null_weights <- function(scale_name, sc, chr) {             # chi-square(1) weights of Q's null
  id <- paste(scale_name, chr)
  if (is.null(spectrum_cache[[id]])) spectrum_cache[[id]] <- lof_spectrum(sc$X, scaled(sc$Sigma0[[chr]], sc$k))$lambda
  spectrum_cache[[id]]
}
# Q, Davies p and residual share (Q / total chi-square) of profiles B (rows), per chromosome.
residual_test <- function(B, d, chr, sc) {
  out <- data.table(Q = NA_real_, p = NA_real_, floored = NA, share = NA_real_, T = NA_real_)[rep(1, nrow(B))]
  for (k in unique(chr)) {
    i <- which(chr == k)
    BR <- sweep(B[i, , drop = FALSE], 2, sc$k, `*`); SR <- scaled(sc$Sigma0[[k]], sc$k)
    fit <- lof_davies(BR, SR, d[i], sc$X)
    total <- d[i] * rowSums((BR %*% solve((SR + t(SR)) / 2)) * BR)
    set(out, i, "Q", fit$Q); set(out, i, "p", fit$p_davies); set(out, i, "floored", fit$davies_floored)
    set(out, i, "T", total); set(out, i, "share", fit$Q / total)
  }
  out
}

## 2. Lead variants -------------------------------------------------------------------------
chr <- as.character(leads_raw$meta$chr); d_raw <- leads_raw$meta$d_resid
m <- match(leads_raw$meta$snp_id, leads_log$meta$snp_id); stopifnot(!anyNA(m))
B_log <- leads_log$B[m, ]; d_log <- leads_log$meta$d_resid[m]
lead_tests <- lapply(SCALES, function(sc)
  if (sc$src == "raw") residual_test(leads_raw$B, d_raw, chr, sc) else residual_test(B_log, d_log, chr, sc))
n <- nrow(leads_raw$B); bonferroni <- 0.05 / n
se_location <- sqrt(vapply(seq_len(n), function(i) functional_cov(trap_w, trap_w, Sigma0_raw[[chr[i]]], d_raw[i]), 0))
w_spread <- trap_w * (2 * u - 1)
se_spread <- sqrt(vapply(seq_len(n), function(i) functional_cov(w_spread, w_spread, Sigma0_raw[[chr[i]]], d_raw[i]), 0))
lead_residual <- data.table(snp_id = leads_raw$meta$snp_id, chr = leads_raw$meta$chr, bp = leads_raw$meta$bp,
  raising_allele = leads_raw$meta$effect_allele, d_resid = d_raw,
  location_kgm2 = as.numeric(leads_raw$B %*% trap_w), location_z = as.numeric(leads_raw$B %*% trap_w) / se_location,
  spread_z = as.numeric(leads_raw$B %*% w_spread) / se_spread)
for (sn in names(SCALES)) {
  x <- lead_tests[[sn]]
  set(lead_residual, j = paste0("Q_", sn), value = x$Q); set(lead_residual, j = paste0("p_", sn), value = x$p)
  set(lead_residual, j = paste0("share_", sn), value = x$share); set(lead_residual, j = paste0("bonf_", sn), value = x$p < bonferroni)
  set(lead_residual, j = paste0("floored_", sn), value = x$floored)
}
lead_residual[, total_chi2_normal_score := lead_tests$normal_score$T]

## 3. Sensitivity: residual after the top two principal components ---------------------------
# Whitened normal-score profiles, scaled to unit length; the top two PCs are estimated without
# the lead's own chromosome; the residual is referred to chi-square on 43 df.
Sigma_R_mean <- scaled(Sigma0_raw[["mean"]], to_normal_raw)
B_normal <- sweep(leads_raw$B, 2, to_normal_raw, `*`)
Y <- whiten(B_normal, Sigma_R_mean, d_raw)
top2 <- function(Ym) { U <- Ym / sqrt(rowSums(Ym^2)); eigen(crossprod(U), symmetric = TRUE)$vectors[, 1:2] }
Q_pc <- numeric(n)
for (k in unique(chr)) {
  ii <- which(chr == k); V <- top2(Y[-ii, , drop = FALSE])
  R <- Y[ii, , drop = FALSE] - Y[ii, , drop = FALSE] %*% V %*% t(V); Q_pc[ii] <- rowSums(R^2)
}
lead_residual[, `:=`(Q_PC = Q_pc, p_PC = pchisq(Q_pc, n_tau - 2L, lower.tail = FALSE))][, bonf_PC := p_PC < bonferroni]
# principal angles between the shift + spread frame and the top-2 PC plane (whitened)
Sm <- sym_pow(Sigma_R_mean, -0.5)
principal_angles <- sort(acos(pmin(1, svd(t(qr.Q(qr(Sm %*% cbind(1, z)))) %*% qr.Q(qr(top2(Y))))$d)) * 180 / pi)
size_adjusted <- function(a) residuals(lm(log(a) ~ log(lead_residual$total_chi2_normal_score)))
pc_sensitivity <- data.table(
  quantity = c("principal_angle_1_deg", "principal_angle_2_deg", "n_bonf_frame", "n_bonf_PC", "n_bonf_both",
               "spearman_resid_logQ_vs_logQPC_on_logT"),
  value = c(principal_angles, sum(lead_residual$bonf_normal_score), sum(lead_residual$bonf_PC),
            sum(lead_residual$bonf_normal_score & lead_residual$bonf_PC),
            cor(size_adjusted(lead_residual$Q_normal_score), size_adjusted(lead_residual$Q_PC), method = "spearman")))

## 4. Calibration by simulation ---------------------------------------------------------------
# 20,000 pure shift + spread normal-score profiles: (a, c) resampled from the leads' own fits
# to [1, z], plus noise N(0, Sigma_R / d) of the drawn lead.
X_frame <- cbind(1, z)
theta <- B_normal %*% t(solve(crossprod(X_frame), t(X_frame)))
draw <- sample.int(n, N_SIM, replace = TRUE)
B_sim <- theta[draw, 1] + outer(theta[draw, 2], z)
for (k in unique(chr[draw])) {
  jj <- which(chr[draw] == k)
  L <- chol((scaled(Sigma0_raw[[k]], to_normal_raw) + t(scaled(Sigma0_raw[[k]], to_normal_raw))) / 2)
  B_sim[jj, ] <- B_sim[jj, ] + (matrix(rnorm(length(jj) * n_tau), length(jj)) %*% L) / sqrt(d_raw[draw][jj])
}
sim_scale <- list(k = rep(1, n_tau), X = X_frame, Sigma0 = lapply(Sigma0_raw, scaled, k = to_normal_raw), src = "raw")
sim_tests <- residual_test(B_sim, d_raw[draw], chr[draw], sim_scale)

setorder(lead_residual, p_normal_score, -Q_normal_score)
top10 <- lead_residual[1:10, .(snp_id, chr, bp, location_z, spread_z, Q_normal_score, p_normal_score,
                              share_normal_score, p_raw, p_log, p_PC)]
write_result(lead_residual, OUT, "lead_residual")
write_result(pc_sensitivity, OUT, "pc_sensitivity")
write_result(top10, OUT, "top10_residual")

## 5. Calibration on null variants (genome-wide Stage 1) -------------------------------------
# Null variants: MAF >= 0.01 (REGENIE allele frequency), |location z| < 1, more than 1 Mb from
# any lead, complete profile; 6,667 sampled in each of three MAF bands (0.01-0.05, 0.05-0.2,
# 0.2-0.5). The same variants give the raw and log (converted) calibration.
if (!dir.exists(cfg$cohorts$bmi$root)) {
  message("Genome-wide Stage 1 not available: lead-level results written; null calibration skipped.")
  quit(save = "no")
}
stage1_dir <- function(cc) file.path(cfg$cohorts$bmi$root, sprintf("chr%02d", cc))
regenie_comparator <- regenie_files(cfg, "bmi", "comparator")
# the density used here must be the one Stage 1 used on every chromosome
fhat_dev <- max(vapply(1:22, function(cc) {
  j <- jsonlite::read_json(file.path(stage1_dir(cc), "rif_nuisance.json"), simplifyVector = TRUE)
  max(abs(j$fhat - q_raw$fhat) / q_raw$fhat) }, 0))
lead_bp <- split(leads_raw$meta$bp, as.character(leads_raw$meta$chr))
eligible <- list(); counts <- list()
for (cc in 1:22) {
  x <- fread(cmd = paste("zcat", shQuote(file.path(stage1_dir(cc), "scan.stage1.tsv.gz"))), showProgress = FALSE)
  x[, snp_id := as.character(snp_id)]
  beta_cols <- grep("^beta_tau_", names(x), value = TRUE)
  beta_cols <- beta_cols[match(round(u, 3), round(as.numeric(sub("^beta_tau_", "", beta_cols)), 3))]
  freq <- fread(regenie_comparator[chr == cc]$path, select = c("ID", "A1FREQ"), showProgress = FALSE)
  freq <- freq[!(ID %in% ID[duplicated(ID)])]
  af <- freq$A1FREQ[match(x$snp_id, freq$ID)]; maf <- pmin(af, 1 - af)
  Bx <- as.matrix(x[, beta_cols, with = FALSE]); storage.mode(Bx) <- "double"
  location_z <- as.numeric(Bx %*% trap_w) /
    sqrt(as.numeric(t(trap_w) %*% Sigma0_raw[[as.character(cc)]] %*% trap_w) / x$d_resid)
  lb <- sort(lead_bp[[as.character(cc)]])
  far <- if (length(lb)) {
    j <- findInterval(x$bp, lb)
    below <- ifelse(j >= 1, x$bp - lb[pmax(j, 1)], Inf); above <- ifelse(j < length(lb), lb[pmin(j + 1, length(lb))] - x$bp, Inf)
    pmin(below, above) > 1e6 } else rep(TRUE, nrow(x))
  maf_ok <- !is.na(maf) & maf >= 0.01
  ok <- maf_ok & far & abs(location_z) < 1 & rowSums(is.na(Bx)) == 0 & x$d_resid > 0
  counts[[cc]] <- data.table(chr = cc, n_scan = nrow(x), n_no_af = sum(is.na(maf)), n_maf_ok = sum(maf_ok),
                             n_far = sum(maf_ok & far), n_eligible = sum(ok))
  eligible[[cc]] <- list(meta = data.table(snp_id = x$snp_id[ok], chr = cc, bp = x$bp[ok], d_resid = x$d_resid[ok],
                                           maf = maf[ok], loc_z = location_z[ok]), B = Bx[ok, , drop = FALSE])
  rm(x, Bx); invisible(gc(verbose = FALSE))
}
null_counts <- rbindlist(counts)
EM <- rbindlist(lapply(eligible, `[[`, "meta")); EB <- do.call(rbind, lapply(eligible, `[[`, "B")); rm(eligible)
EM[, maf_bin := cut(maf, c(0.01, 0.05, 0.2, 0.5 + 1e-9), right = FALSE, labels = c("0.01-0.05", "0.05-0.2", "0.2-0.5"))]
per_bin <- ceiling(N_NULL / 3)
pick <- unlist(lapply(split(seq_len(nrow(EM)), EM$maf_bin), function(ii) if (length(ii) <= per_bin) ii else sample(ii, per_bin)))
null_meta <- EM[pick]; null_B <- EB[pick, , drop = FALSE]; rm(EB)
null_chr <- as.character(null_meta$chr)
null_tests <- lapply(SCALES_NULL, function(sc) residual_test(null_B, null_meta$d_resid, null_chr, sc))

# Davies against Imhof on 300 null variants (normal scores)
davies_vs_imhof <- max(vapply(1:300, function(i) {
  lam <- null_weights("normal_score", SCALES$normal_score, null_chr[i])
  abs(CompQuadForm::davies(null_tests$normal_score$Q[i], lam, acc = 1e-12, lim = 1e5)$Qq -
      CompQuadForm::imhof(null_tests$normal_score$Q[i], lam, epsabs = 1e-14, epsrel = 1e-14, limit = 1e4)$Qq) }, 0))

# E[Q], and its value under the |location z| < 1 selection, which truncates the location part
EQ_of <- function(scale_name, sc, chrs) median(vapply(unique(chrs), function(k) sum(null_weights(scale_name, sc, k)), 0))
EQ_normal <- EQ_of("normal_score", SCALES$normal_score, chr)
SR <- scaled(Sigma0_raw[["mean"]], to_normal_raw); M <- diag(n_tau) - X_frame %*% solve(crossprod(X_frame), t(X_frame))
a_R <- trap_w / to_normal_raw
g_sel <- as.numeric(t(a_R) %*% SR %*% M %*% solve(SR) %*% M %*% SR %*% a_R) / as.numeric(t(a_R) %*% SR %*% a_R)
EQ_selected <- EQ_normal - g_sel * (1 - (1 - 2 * dnorm(1) / (2 * pnorm(1) - 1)))

calibration_row <- function(Q, p, set, scale, EQ) data.table(set = set, scale = scale, n = length(Q),
  mean_Q_over_EQ = mean(Q) / EQ, frac_p_lt_0.05 = mean(p < 0.05), frac_p_lt_1e3 = mean(p < 1e-3),
  frac_p_lt_1e4 = mean(p < 1e-4), lambda_gc = lambda_gc(p))
calibration <- rbindlist(c(
  lapply(names(null_tests), function(sn) calibration_row(null_tests[[sn]]$Q, null_tests[[sn]]$p, "null_all", sn,
                                                         EQ_of(sn, SCALES_NULL[[sn]], null_chr))),
  lapply(levels(null_meta$maf_bin), function(b) { i <- null_meta$maf_bin == b
    calibration_row(null_tests$normal_score$Q[i], null_tests$normal_score$p[i], paste0("null_maf_", b), "normal_score", EQ_normal) }),
  list(calibration_row(sim_tests$Q, sim_tests$p, "sim_shift_spread", "normal_score", EQ_normal))))
calibration[, `:=`(E_Q_normal_score = EQ_normal, E_Q_under_locz_selection = EQ_selected)]

qq <- rbindlist(c(
  list(data.table(set = "leads", qq_quantiles(-log10(lead_residual$p_normal_score)))),
  list(data.table(set = "null_all", qq_quantiles(-log10(null_tests$normal_score$p)))),
  lapply(levels(null_meta$maf_bin), function(b) data.table(set = paste0("null_maf_", b),
    qq_quantiles(-log10(null_tests$normal_score$p[null_meta$maf_bin == b])))),
  list(data.table(set = "sim_shift_spread", qq_quantiles(-log10(sim_tests$p))))))

summary_by_scale <- rbindlist(lapply(names(SCALES), function(sn) data.table(scale = sn, n_lead = n, bonferroni = bonferroni,
  n_bonf = sum(lead_residual[[paste0("bonf_", sn)]]), n_nominal_0.05 = sum(lead_residual[[paste0("p_", sn)]] < 0.05),
  expected_nominal = 0.05 * n, median_share_leads = median(lead_residual[[paste0("share_", sn)]]),
  n_floored = sum(lead_residual[[paste0("floored_", sn)]]))))
summary_by_scale[scale == "normal_score", `:=`(median_share_null = median(null_tests$normal_score$share),
                                                median_share_sim = median(sim_tests$share))]
checks <- data.table(quantity = c("n_null_eligible", "n_null_sampled", "n_sim", "davies_imhof_max_abs_diff_300",
                                  "fhat_max_rel_dev_vs_rif_nuisance", "E_Q_normal_score", "E_Q_under_locz_selection", "n_bonf_PC"),
                     value = c(nrow(EM), nrow(null_meta), N_SIM, davies_vs_imhof, fhat_dev, EQ_normal, EQ_selected,
                               sum(lead_residual$bonf_PC)))

write_result(calibration, OUT, "calibration")
write_result(rbind(summary_by_scale, checks, fill = TRUE), OUT, "summary")
write_result(null_counts, OUT, "null_counts_by_chr")
write_result(qq, OUT, "qq")
