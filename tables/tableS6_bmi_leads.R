#!/usr/bin/env Rscript
# Supplementary Table S6: moments, polygenic-score interaction and model scores of the BMI leads
#
# What: one row per BMI lead variant (639; BMI-increasing allele). Columns: moment functionals
#   of the log-BMI Stage 1 profile with their standard errors; the per-lead regression on the
#   polygenic score (PRS) excluding the variant's chromosome (main effect beta, interaction
#   delta, and kappa_i = delta / beta); and the AIC-like and BIC-like scores of the seven
#   Stage 2 models on raw and log BMI, with the AIC-like winner, runner-up and gap.
#   (The PRS results store the interaction as gamma_log; the paper calls it delta.)
# Paper: Supplementary Table S6
# Input: results/06_traits/moments/lead_moments.csv,
#   results/05_bmi/prs_interaction/{per_lead_prs,kappa_by_scale}.csv,
#   results/05_bmi/model_scores/{lead_scores,lead_winners}.csv
# Output: paper_assets/tableS6_bmi_leads.csv, _caption.md
# Runs: locally from results/ in about 1 second

suppressPackageStartupMessages({library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load ----------------------------------------------------------------------------------
moments <- as.data.table(read_result("06_traits/moments", "lead_moments"))[arm == "log_bmi"]
prs <- as.data.table(read_result("05_bmi/prs_interaction", "per_lead_prs"))
kappa_log <- as.data.table(read_result("05_bmi/prs_interaction", "kappa_by_scale"))[scale == "log" & se_type == "HC0"]
scores <- as.data.table(read_result("05_bmi/model_scores", "lead_scores"))[set == "primary"]
winners <- as.data.table(read_result("05_bmi/model_scores", "lead_winners"))[set == "primary" & criterion == "AIC"]
stopifnot(nrow(moments) == 639, nrow(prs) == 639, setequal(moments$snp_id, prs$snp_id),
          setequal(unique(scores$snp_id), moments$snp_id), setequal(unique(scores$scale), c("raw", "log")))

## 2. Column blocks ---------------------------------------------------------------------------
s4 <- function(x) signif(x, 4)
mom <- moments[, .(snp_id, chr, bp, effect_allele, other_allele,
                   delta_mu = s4(delta_mu), delta_mu_se = s4(delta_mu_se),
                   delta_var = s4(delta_var), delta_var_se = s4(delta_var_se),
                   delta_lambda2 = s4(delta_lambda2), delta_lambda2_se = s4(delta_lambda2_se),
                   delta_lambda3 = s4(delta_lambda3), delta_lambda3_se = s4(delta_lambda3_se),
                   skew_beyond_frame_z = round(skew_beyond_frame_z, 3),
                   lof_p_davies = signif(lof_p_davies, 3), het_resid_z = round(het_resid_z, 3))]
prs_cols <- prs[, .(snp_id, prs_beta = s4(beta_log), prs_beta_se = s4(se_beta_log_hc0), prs_delta = s4(gamma_log),
                    prs_delta_se = s4(se_gamma_log_hc0), kappa_i = s4(gamma_log / beta_log))]
MODELS <- c("additive", "log_additive", "shifted", "both", "logscale_vqtl", "shifted_lognormal_full", "three_components")
stopifnot(setequal(unique(scores$model), MODELS))
# One column per model for a criterion (AIC or BIC) on one scale, e.g. AIC_additive_log.
score_block <- function(crit, scl) {
  w <- dcast(scores[scale == scl], snp_id ~ model, value.var = crit)
  setcolorder(w, c("snp_id", MODELS))
  setnames(w, MODELS, sprintf("%s_%s_%s", crit, MODELS, scl))
  w[, (names(w)[-1]) := lapply(.SD, round, 1), .SDcols = names(w)[-1]]
  w
}
winner_block <- function(scl) {
  w <- winners[scale == scl, .(snp_id, winner, runner = runner_up, gap = round(gap, 2))]
  setnames(w, c("winner", "runner", "gap"), paste0(c("winner_", "runner_", "gap_"), scl))
}

## 3. Join and save -------------------------------------------------------------------------
out <- Reduce(function(a, b) merge(a, b, by = "snp_id"),
              list(mom, prs_cols, score_block("AIC", "raw"), winner_block("raw"), score_block("AIC", "log"), winner_block("log"),
                   score_block("BIC", "raw"), score_block("BIC", "log")))
stopifnot(nrow(out) == 639)
setorder(out, chr, bp)
near_tie <- c(raw = mean(out$gap_raw < 2), log = mean(out$gap_log < 2))
cap <- paste(
  "Table S6. Moments, polygenic-score interaction and model scores of the 639 BMI lead variants. One row per lead",
  "variant (effect allele, the BMI-increasing allele; n = 461,320). Moments of the Stage 1 quantile-effect profile on",
  "log BMI, with pooled-covariance standard errors (Methods): mean effect, change in variance, and changes in L-moment",
  "spread and skew; skew beyond shift and spread (z); lack of fit of the shift and spread frame (exact one-sided P,",
  "Davies' method, not adjusted for multiple comparisons); and the residual from the common tilt of the 639 leads (z).",
  "Polygenic-score interaction: regression of log BMI on the variant, a polygenic risk score excluding its chromosome",
  "(per s.d.), their product, age, sex and ten genetic principal components, with heteroskedasticity-robust",
  "standard errors; β is the main effect, δ the interaction and κᵢ = δ/β (noisy for single leads; common κ",
  sprintf("%.3f, s.e. %.3f; Fig. 6). Model scores, on raw and log BMI: AIC-like and BIC-like scores (full-covariance",
          kappa_log$kappa, kappa_log$kappa_se_boot),
  "residual plus 2 or log n per parameter; Methods) of seven Stage 2 models: additive; log-additive (multiplicative",
  "about zero); shifted log-additive (about a fitted origin); additive + shifted log-additive; log-scale vQTL;",
  "additive + log-scale vQTL; and shift + spread + bend. Winner, runner-up and gap refer to the AIC-like score;",
  sprintf("the runner-up is within 2 units for %.0f%% (raw) and %.0f%% (log) of leads (Supplementary Fig. S13).",
          100 * near_tie[["raw"]], 100 * near_tie[["log"]]),
  "AIC, Akaike information criterion; BIC, Bayesian information criterion; BMI, body mass index; s.d., standard",
  "deviation; s.e., standard error; vQTL, variance quantitative trait locus.")
save_table(as.data.frame(out), "tableS6_bmi_leads", cap,
           results_used = c("06_traits/moments", "05_bmi/prs_interaction", "05_bmi/model_scores"))
