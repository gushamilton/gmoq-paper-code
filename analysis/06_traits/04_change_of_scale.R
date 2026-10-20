#!/usr/bin/env Rscript
# Change of scale of Stage 1 output: check against a separate GWAS, then every trait on three scales
#
# What: For any strictly increasing transformation T of the phenotype, Stage 1 output converts
#   exactly without rerunning the GWAS: beta_T(u) = T'(Q(u)) beta(u) and Sigma_T = D Sigma0 D
#   with D = diag(T'(Q(u))); only the kernel density inside the RIF is not exact. (1) We check
#   this by converting the raw-BMI Stage 1 profiles of the 639 BMI leads to the log scale and
#   comparing them with the separately run log-BMI Stage 1. (2) Every trait is then converted
#   to the raw, log and Box-Cox (normalising power) scales and given the same moment reading as
#   01_moments_and_principal_components.R (without bootstrap or detection floor): spread and
#   skew counts, lack of fit, heterogeneity tau_h and principal components. Component counts
#   cannot change with the scale (in the whitened metric a change of scale is a rotation shared
#   by all leads); spread counts and tau_h can.
# Paper: Fig. 5a (grey symbols); Supplementary Table S5; Methods "Change of scale"
# Input: aligned lead-variant profiles and per-chromosome pooled covariance
#   (results/01_data/lead_profiles), baseline quantiles (results/01_data/baseline_quantiles),
#   normalising Box-Cox powers (results/06_traits/scale/boxcox_lambda.csv)
# Output: results/06_traits/scale/{identity_summary,identity_per_lead,scale_comparison}.csv
# Runs: locally from results/ in about 2 minutes on one core

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
u <- tau_grid(cfg); n_tau <- length(u)
n_null_pca <- as.integer(cfg$estimation$pca_parallel_draws)
seed <- as.integer(cfg$estimation$boot_seed)

lead_profiles <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
baseline <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))
baseline_quantiles <- function(trait_name) {
  x <- baseline[baseline$trait == trait_name]
  as.numeric(x$q[match(round(u, 3), round(as.numeric(x$tau), 3))])
}
boxcox_powers <- as.data.table(read_result("06_traits/scale", "boxcox_lambda"))

## 1. Raw-BMI Stage 1 converted to log against the log-BMI Stage 1
bmi_raw <- lead_set_matrix(lead_profiles, "bmi_A_regenie", u)
bmi_log <- lead_set_matrix(lead_profiles, "log_bmi_A_regenie", u)
m <- match(bmi_raw$meta$snp_id, bmi_log$meta$snp_id)
stopifnot(!anyNA(m), all(bmi_raw$meta$effect_allele == bmi_log$meta$effect_allele[m]))
B_log <- bmi_log$B[m, ]; B_raw <- bmi_raw$B
q_raw <- baseline_quantiles("bmi"); q_log <- baseline_quantiles("log_bmi")
converted <- sweep(B_raw, 2, 1 / q_raw, `*`)           # T = log, T'(Q) = 1 / Q
rel_error <- abs(converted - B_log) / pmax(abs(B_log), .Machine$double.eps)
fto <- match(cfg$prs$fto_snp, bmi_raw$meta$snp_id)     # rs1421085
# Leads are aligned to the BMI-increasing allele; the correlation is also given on the Stage 1
# effect-allele coding, where the between-lead sign spread makes it larger (errors are the same).
s1 <- ifelse(as.logical(bmi_raw$meta$aligned_flip), -1, 1)
identity_summary <- data.table(
  quantity = c("n_leads", "n_cells", "mean_abs_rel_error_pct", "median_abs_rel_error_pct",
               "p95_abs_rel_error_pct", "max_abs_rel_error_pct", "overall_r", "overall_slope",
               "fto_rs1421085_mean_abs_rel_error_pct", "fto_rs1421085_max_abs_rel_error_pct",
               "quantile_map_max_abs_diff_log_q", "overall_r_stage1_allele_coding"),
  value = c(nrow(B_raw), length(B_raw), 100 * mean(rel_error), 100 * median(rel_error),
            100 * quantile(rel_error, 0.95), 100 * max(rel_error),
            cor(as.numeric(converted), as.numeric(B_log)),
            unname(coef(lm(as.numeric(B_log) ~ as.numeric(converted)))[2]),
            100 * mean(rel_error[fto, ]), 100 * max(rel_error[fto, ]),
            max(abs(log(q_raw) - q_log)),                 # quantiles map exactly: log Q_raw = Q_log
            cor(as.numeric(s1 * converted), as.numeric(s1 * B_log))))
identity_per_lead <- data.table(
  snp_id = bmi_raw$meta$snp_id, chr = bmi_raw$meta$chr,
  mean_abs_rel_error_pct = 100 * rowMeans(rel_error), max_abs_rel_error_pct = 100 * apply(rel_error, 1, max),
  profile_r = vapply(seq_len(nrow(B_raw)), function(i) cor(converted[i, ], B_log[i, ]), 0))
print(identity_summary)

## 2. Every trait on the raw, log and Box-Cox scales
survey <- cfg$cohorts$survey$traits
survey_traits <- names(survey)[vapply(survey, function(x) !identical(x$fig6, FALSE) &&
                                        !is.null(x[["n", exact = TRUE]]), logical(1))]
traits <- data.table(trait = c("bmi", "weight", survey_traits),
                     lead_set = c("bmi_A_regenie", "weight_A_regenie", paste0(survey_traits, "_A_deltamu")))

# Profiles and covariance on the raw measurement scale (a log-analysed trait is converted back
# with D = Q_raw), from which every Box-Cox scale is one more exact conversion.
raw_profiles <- function(trait_name, lead_set) {
  leads <- lead_set_matrix(lead_profiles, lead_set, u)
  q_analysed <- baseline_quantiles(trait_name)
  is_log <- grepl("^log_", trait_name)
  q_raw <- if (is_log) exp(q_analysed) else q_analysed
  B_raw <- if (is_log) sweep(leads$B, 2, q_raw, `*`) else leads$B
  sigma0 <- p05_sigma0(leads$source, u)
  sigma0_raw <- if (is_log) lapply(sigma0, function(S) change_scale(matrix(0, 1, n_tau), S, q_raw)$Sigma0) else sigma0
  list(meta = leads$meta, q_raw = q_raw, B_raw = B_raw, sigma0_raw = sigma0_raw)
}

on_scale <- function(trait_name, lead_set, scale) {
  rp <- raw_profiles(trait_name, lead_set)
  lambda <- switch(scale, raw = 1, log = 0, boxcox = boxcox_powers[trait == trait_name]$lambda_main)
  if (!is.finite(lambda) || (lambda != 1 && any(rp$q_raw <= 0))) return(NULL)   # Box-Cox needs Y > 0
  D <- rp$q_raw^(lambda - 1)                          # T'(Q_raw) for the Box-Cox power lambda
  B <- sweep(rp$B_raw, 2, D, `*`)
  sigma0 <- lapply(rp$sigma0_raw, function(S) change_scale(matrix(0, 1, n_tau), S, D)$Sigma0)
  q_T <- if (scale == "raw") rp$q_raw else boxcox(rp$q_raw, lambda)
  md <- moment_decomposition(paste(trait_name, scale, sep = "@"), rp$meta, B, sigma0, q_T, u,
                             full = FALSE, npa = n_null_pca, seed = seed)
  het <- md$het[basis == "unit_se" & component == "c2_on_c1" & excess == "n_minus_1"]
  het_exact <- md$het[basis == "exact" & excess == "n_minus_1"]
  data.table(trait = trait_name, scale = scale, lambda = lambda, n_leads = md$n, n_components = md$n_components,
             pc1_pct = md$spec$pct[1], pc2_pct = md$spec$pct[2], pc3_p = md$par$p_exceeds_null[3],
             md$skew[, .(n_spread_lambda2_bonf, n_spread_exact_bonf, n_skew_raw_bonf, n_skew_beyond_bonf,
                         median_spread_lambda2_z, cor_spread_location)],
             lof_frac_p_lt_0.05 = md$lof[frame == "shift_spread"]$frac_p_lt_0.05,
             lof_n_bonf_rejected = md$lof[frame == "shift_spread"]$n_bonf_rejected,
             tau_unit_se = het$tau, r_unit_se = het$r, p_het = het$p, tau_over_r_exact = het_exact$tau_over_r)
}
scale_comparison <- rbindlist(lapply(seq_len(nrow(traits)), function(i)
  rbindlist(lapply(c("raw", "log", "boxcox"), function(s) on_scale(traits$trait[i], traits$lead_set[i], s)))))
print(scale_comparison[trait == "bmi"])

## 3. Write
out <- "06_traits/scale"
write_result(identity_summary, out, "identity_summary")
write_result(identity_per_lead, out, "identity_per_lead")
write_result(scale_comparison, out, "scale_comparison")
