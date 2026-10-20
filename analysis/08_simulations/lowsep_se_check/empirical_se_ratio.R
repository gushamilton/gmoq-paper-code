#!/usr/bin/env Rscript
# Standard errors at the empirical mixture separations: real-data arm
#
# What: for each empirical mixture case and each Stage 2 mixture parameter (gamma,
#   beta_left, beta_right), the ratio of the pooled standard error (Sigma_0 / d_resid
#   through OLS) to the 1,000-block jackknife standard error, on the SAME Stage 1 betas and
#   the SAME weight matrix W, so the point estimates are identical and only the covariance
#   source differs. Also each case's mixture separation from the fitted component
#   parameters stored with the weights.
# Paper: Supplementary Information, study 7; Supplementary Table S2
# Input: Stage 2 mixture results computed with the pooled and with the jackknife covariance
#   (one file per chromosome), and the fitted mixture weights, for three arms:
#     old_weights_genome  earlier genome-wide pooled + jackknife result pairs (weights
#                         built on the Stage 1 kernel density), all autosomes;
#     corrected_weights   testosterone x4 on the model-density weights, chromosomes 19-22;
#     scan_weights        JCV/MCV VP1 on the 44-antigen scan weights, chromosomes 6 and 22,
#                         Stage 1 rerun with the jackknife payload.
# Output: <out-dir>/{empirical_se_ratio,separations,reproducibility_checks}.csv (aggregates)
# Runs: on a cluster, after parts A and B of run_lowsep_se_check.sh

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(flag, " is required")
  args[i + 1L]
}

## User-set inputs (passed by run_lowsep_se_check.sh)
run_dir <- getarg("--run-dir")          # parts A and B of this check (Stage 2 reruns)
out <- getarg("--out-dir")
# Earlier genome-wide runs with pooled and jackknife Stage 2, as sprintf templates
# (%1$s = trait), comma-separated. The second layout holds the chromosomes that the first
# run reused from an earlier run.
old_run_templates <- strsplit(getarg("--old-run-templates"), ",", fixed = TRUE)[[1]]
corrected_dir <- getarg("--corrected-dir")   # testosterone model-density weights and Stage 2
scan_dir <- getarg("--scan-dir")             # 44-antigen scan output (JCV/MCV VP1)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

PARAMS <- c("gamma", "beta_left", "beta_right")
KEY <- c("snp_id", "effect_allele", "other_allele")

## 1. Helpers
read_res <- function(f) {
  cols <- c(KEY, PARAMS, paste0(PARAMS, "_se"))
  fread(cmd = paste("gzip -dc", shQuote(f)), select = cols, showProgress = FALSE)
}

# Per variant and parameter: |difference in estimate| (0 by construction), pooled z, and
# the pooled / jackknife SE ratio.
pair_rows <- function(pooled_file, jackknife_file) {
  p <- read_res(pooled_file); j <- read_res(jackknife_file)
  m <- merge(p, j, by = KEY, suffixes = c(".p", ".j"))
  rbindlist(lapply(PARAMS, function(k) {
    data.table(parameter = k,
               beta_diff = abs(m[[paste0(k, ".p")]] - m[[paste0(k, ".j")]]),
               z = m[[paste0(k, ".p")]] / m[[paste0(k, "_se.p")]],
               ratio = m[[paste0(k, "_se.p")]] / m[[paste0(k, "_se.j")]])
  }))
}

# Distribution of the ratio among all variants and among those with |z_pooled| > 3 or > 5.
summarise_ratio <- function(d, case, arm, chrs) {
  d <- d[is.finite(ratio) & is.finite(z)]
  rbindlist(lapply(c(0, 3, 5), function(zc) {
    s <- d[abs(z) > zc]
    s[, .(case = case, arm = arm, chromosomes = chrs, z_stratum = paste0("|z_pooled|>", zc),
          n_snps = .N,
          median_pooled_over_jk = if (.N) median(ratio) else NA_real_,
          q05 = if (.N) quantile(ratio, 0.05) else NA_real_,
          q95 = if (.N) quantile(ratio, 0.95) else NA_real_,
          mean_pooled_over_jk = if (.N) mean(ratio) else NA_real_,
          frac_ratio_below_0.9 = if (.N) mean(ratio < 0.9) else NA_real_,
          max_abs_beta_diff = if (.N) max(beta_diff) else NA_real_),
      by = parameter]
  }))
}

# Separation S = (m_upper - m_lower) / sqrt((v_lower + v_upper) / 2), with skew-normal
# moments m = xi + omega d sqrt(2/pi), v = omega^2 (1 - 2 d^2/pi), d = alpha/sqrt(1+alpha^2).
sn_moments <- function(xi, omega, alpha) {
  d <- alpha / sqrt(1 + alpha^2)
  list(m = xi + omega * d * sqrt(2 / pi), v = omega^2 * (1 - 2 * d^2 / pi))
}
separation <- function(weights_file, case, source) {
  w <- readRDS(weights_file); cp <- as.data.table(w$components)
  lo <- cp[component == "lower"]; up <- cp[component == "upper"]
  a <- sn_moments(lo$xi, lo$omega, lo$alpha); b <- sn_moments(up$xi, up$omega, up$alpha)
  data.table(case = case, weight_source = source, weight_file_chr = w$provenance$chromosome,
             density_weight_mode = as.character(w$provenance$density_weight_mode),
             n_samples = w$provenance$n_samples, pi_upper = up$weight,
             mean_lower = a$m, mean_upper = b$m, sd_lower = sqrt(a$v), sd_upper = sqrt(b$v),
             separation = (b$m - a$m) / sqrt((a$v + b$v) / 2))
}

# Largest difference between a rerun and the stored result of the same pooled analysis.
rerun_check <- function(new_file, stored_file, case, chr, label) {
  a <- read_res(new_file); b <- read_res(stored_file); m <- merge(a, b, by = KEY)
  data.table(case = case, chr = chr, check = label, n = nrow(m),
    max_abs_beta_diff = max(sapply(PARAMS, function(k) max(abs(m[[paste0(k, ".x")]] - m[[paste0(k, ".y")]])))),
    max_rel_se_diff = max(sapply(PARAMS, function(k) max(abs(m[[paste0(k, "_se.x")]] / m[[paste0(k, "_se.y")]] - 1)))))
}

res <- list(); seps <- list(); checks <- list()

## 2. Earlier genome-wide pooled + jackknife pairs (Stage 1 kernel-density weights)
old_cases <- list(
  testosterone = c("testosterone", "skew"),
  testosterone_shrink_0.2 = c("testosterone_shrink_0_2", "skew"),
  testosterone_shrink_0.5 = c("testosterone_shrink_0_5", "skew"),
  testosterone_shrink_0.7 = c("testosterone_shrink_0_7", "skew"),
  jcv_log_all = c("jcv_log_all", "mixture"),
  mcv_log_all = c("mcv_log_all", "mixture"))
for (case in names(old_cases)) {
  dn <- old_cases[[case]][1]; model <- old_cases[[case]][2]
  rows <- list(); used <- integer()
  for (c in 1:22) {
    cc <- sprintf("chr%02d", c)
    cands <- file.path(sprintf(old_run_templates, dn), cc)
    for (d in cands) {
      pf <- file.path(d, "stage2", model, "pooled", "result.tsv.gz")
      jf <- file.path(d, "stage2", model, "jackknife", "result.tsv.gz")
      if (file.exists(pf) && file.exists(jf)) { rows[[cc]] <- pair_rows(pf, jf); used <- c(used, c); break }
    }
    if (c == 22) {
      wcands <- file.path(cands, "weights", paste0("weights_", model, ".rds"))
      wf <- wcands[file.exists(wcands)][1]
      if (!is.na(wf)) seps[[length(seps) + 1L]] <- separation(wf, case, "old_stage1_rif_density")
    }
  }
  if (length(rows)) res[[length(res) + 1L]] <- summarise_ratio(
    rbindlist(rows), case, "old_weights_genome", paste(range(used), collapse = "-"))
  message(case, ": earlier-weight chromosomes ", length(used))
}

## 3. Testosterone on the model-density weights (part B), chromosomes 19-22
for (t in c("testosterone", "testosterone_shrink_0.2", "testosterone_shrink_0.5", "testosterone_shrink_0.7")) {
  rows <- list()
  for (c in 19:22) {
    d <- file.path(run_dir, "testosterone_corrected", t, sprintf("chr%02d", c))
    pf <- file.path(d, "pooled", "result.tsv.gz"); jf <- file.path(d, "jackknife", "result.tsv.gz")
    if (!(file.exists(pf) && file.exists(jf))) { message("missing corrected ", t, " chr", c); next }
    rows[[length(rows) + 1L]] <- pair_rows(pf, jf)
    # The pooled rerun must reproduce the stored Stage 2 result on the same weights.
    tu <- gsub(".", "_", t, fixed = TRUE)
    ref <- file.path(corrected_dir, "stage2_corrected", tu, "full", tu,
                     sprintf("chr%02d", c), "stage2", "skew", "pooled", "result.tsv.gz")
    if (file.exists(ref))
      checks[[length(checks) + 1L]] <- rerun_check(pf, ref, t, c, "pooled_vs_stored_corrected_reprojection")
  }
  if (length(rows)) res[[length(res) + 1L]] <- summarise_ratio(rbindlist(rows), t, "corrected_weights", "19-22")
  wf <- file.path(corrected_dir, "corrected_weights", t, "chr22", "weights_skew.rds")
  seps[[length(seps) + 1L]] <- separation(wf, t, "corrected_parametric_model_density")
}

## 4. JCV/MCV VP1 on the 44-antigen scan weights (part A), chromosomes 6 and 22
for (t in c("jcv_vp1_log_all", "mcv_vp1_log_all")) {
  rows <- list()
  for (c in c(6, 22)) {
    cc <- sprintf("chr%02d", c)
    d <- file.path(run_dir, "antibody", t, cc, "stage2")
    pf <- file.path(d, "pooled", "result.tsv.gz"); jf <- file.path(d, "jackknife", "result.tsv.gz")
    if (!(file.exists(pf) && file.exists(jf))) { message("missing antibody ", t, " ", cc); next }
    rows[[length(rows) + 1L]] <- pair_rows(pf, jf)
    # The pooled arm of the Stage 1 rerun must reproduce the stored scan result.
    ref <- file.path(scan_dir, t, cc, "stage2", "result.tsv.gz")
    checks[[length(checks) + 1L]] <- rerun_check(pf, ref, t, c, "rerun_pooled_vs_stored_scan")
  }
  if (length(rows)) res[[length(res) + 1L]] <- summarise_ratio(rbindlist(rows), t, "scan_weights", "6,22")
  seps[[length(seps) + 1L]] <- separation(
    file.path(scan_dir, t, "chr22", "weights", "weights_mixture.rds"), t, "scan_parametric_model_density")
}

## 5. Write the aggregates
fwrite(rbindlist(res), file.path(out, "empirical_se_ratio.csv"))
fwrite(rbindlist(seps, fill = TRUE), file.path(out, "separations.csv"))
fwrite(rbindlist(checks), file.path(out, "reproducibility_checks.csv"))
message("wrote ", out)
