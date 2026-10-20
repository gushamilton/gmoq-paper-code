# Shared helpers for the BMI analyses
#
# What: small functions used by several scripts in this folder: reading the Stage 1
#   quantile profiles of the 639 BMI lead variants and their pooled covariance, and
#   building the weight matrices of the seven BMI models.
# Paper: Fig. 4; Supplementary Figs. S13, S14; Supplementary Table S6
# Input: results/01_data/lead_profiles/{lead_profiles,sigma0_<source>}.csv
# Output: none (sourced by the other scripts)
# Runs: sourced; no work of its own

# Lead-variant profiles of one lead set (e.g. "log_bmi_A_regenie"): the n x 45 matrix of
# Stage 1 quantile effects B (rows aligned to the BMI-increasing allele), the per-lead
# metadata (chromosome, position, alleles, d_resid = residual sum of squares of the dosage,
# which scales the pooled covariance Sigma0 / d), and the Stage 1 run it came from.
# Rows keep the order of the file unless `order_by_position` is TRUE.
read_lead_set <- function(profiles, set, u, order_by_position = FALSE) {
  x <- data.table::as.data.table(profiles)[lead_set == set]
  if (!nrow(x)) stop("lead set not found: ", set)
  if (order_by_position) x <- x[order(chr, bp)]
  beta_cols <- paste0("beta_tau_", sprintf("%.3f", u))
  B <- as.matrix(x[, beta_cols, with = FALSE]); storage.mode(B) <- "double"
  list(meta = x[, .(lead_set, snp_id, chr, bp, effect_allele, other_allele, d_resid, delta_mu_z)],
       B = B, source = unique(x$stage1_source))
}

read_lead_profiles <- function() {
  data.table::fread(file.path(repo_path("results"), "01_data", "lead_profiles", "lead_profiles.csv"),
                    showProgress = FALSE)
}

# Pooled Stage 1 covariance Sigma0 (45 x 45, per unit d) of one Stage 1 run, one matrix per
# chromosome ("1" ... "22") plus "mean", the average over chromosomes.
read_lead_sigma0 <- function(source, u) {
  s <- data.table::fread(file.path(repo_path("results"), "01_data", "lead_profiles",
                                   paste0("sigma0_", source, ".csv")), showProgress = FALSE)
  s_cols <- grep("^s_", names(s), value = TRUE)
  col_index <- match(round(u, 3), round(as.numeric(sub("^s_", "", s_cols)), 3))
  out <- lapply(split(s, s$chr), function(x) {
    row_index <- match(round(u, 3), round(as.numeric(x$tau), 3))
    M <- as.matrix(x[, s_cols, with = FALSE])[row_index, col_index, drop = FALSE]
    (M + t(M)) / 2
  })
  names(out) <- sub("^0+", "", names(out))
  out
}

# Weight matrices W (45 x k) of the seven BMI models, on the raw scale. Each column is the
# derivative of the baseline quantile function q(u) with respect to one model parameter,
# evaluated at the fitted shifted log-normal baseline q(u) = c + X(u), X(u) = exp(m + s z_u):
#   additive                 {1}            each allele adds a constant
#   log_additive             {Q(u)}         multiplicative about zero: for Y = exp(m + e), dq/dm = q,
#                                           exact for any error law, so it uses the Stage 1
#                                           quantiles Q(u) and needs no baseline fit
#   shifted                  {X}            multiplicative above the fitted origin c
#   both                     {1, X}         additive + multiplicative (= raw-scale mean + variance)
#   logscale_vqtl            {X, X z}       log-scale location and log-scale SD
#   shifted_lognormal_full   {1, X, X z}    all three baseline parameters
#   three_components         {1, z, z^2}    model-free level, tilt and bend
bmi_model_weights <- function(fit, q_stage1, u) {
  z <- stats::qnorm(u); nt <- length(u)
  X <- shifted_lognormal_weights(fit, u)[, "multiplicative"]
  list(additive = cbind(additive = rep(1, nt)),
       log_additive = cbind(multiplicative_c0 = q_stage1),
       shifted = cbind(multiplicative = X),
       both = cbind(additive = 1, multiplicative = X),
       logscale_vqtl = cbind(meanlog = X, sdlog = X * z),
       shifted_lognormal_full = cbind(additive = 1, meanlog = X, sdlog = X * z),
       three_components = cbind(level = 1, tilt = z, bend = z^2))
}

# Stage 1 quantiles Q(u) and kernel density fhat(u) of a trait on the 45-level grid.
stage1_quantiles <- function(baseline, trait, u) {
  x <- baseline[baseline$trait == trait, ]
  i <- match(round(u, 3), round(as.numeric(x$tau), 3))
  list(q = as.numeric(x$q[i]), fhat = as.numeric(x$fhat[i]))
}

# Lead profiles measured against the polygenic score's own quantile profile. For a lead on
# chromosome c, h(u) = b_S(u) / mean(b_S) is the profile of the score that excludes c, and the
# lead's profile is projected by OLS on [h, h z]: alpha is its level relative to the score and
# rho its extra tilt (a pure change of scale shared with the score predicts rho = 0). Returns
# alpha, rho and their pooled (co)variances per lead.
#   prs_profile: 22 x 45 matrix, row c = profile of the score that leaves out chromosome c.
prs_referenced_tilt <- function(B, d, chrs, Sigma0, prs_profile, u) {
  z <- stats::qnorm(u)
  alpha <- rho <- v_aa <- v_rr <- v_ar <- numeric(nrow(B))
  for (i in seq_len(nrow(B))) {
    h <- prs_profile[chrs[i], ] / mean(prs_profile[chrs[i], ])
    fit <- project_profiles(B[i, ], cbind(h, h * z), Sigma0[[as.character(chrs[i])]], d[i])
    alpha[i] <- fit$theta[1, 1]; rho[i] <- fit$theta[1, 2]
    v_aa[i] <- fit$V[1, 1, 1]; v_rr[i] <- fit$V[1, 2, 2]; v_ar[i] <- fit$V[1, 1, 2]
  }
  list(alpha = alpha, rho = rho, vaa = v_aa, vrr = v_rr, var = v_ar)
}

# Align a lead set's profiles to a reference allele (sign flip where the alleles are swapped).
align_to_reference <- function(info, B, reference) {
  m <- match(reference$snp_id, info$snp_id); stopifnot(!anyNA(m))
  info <- info[m]; B <- B[m, , drop = FALSE]
  same <- info$effect_allele == reference$effect_allele
  swap <- info$effect_allele == reference$other_allele & info$other_allele == reference$effect_allele
  stopifnot(all(same | swap))
  B[swap, ] <- -B[swap, ]
  list(info = reference, d = info$d_resid, B = B)
}
