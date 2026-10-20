#!/usr/bin/env Rscript
# Can a single BMI lead variant tell an additive effect from a multiplicative one?
#
# What: for pairs of one-parameter models of equal rank (additive vs multiplicative about
#   zero, additive vs multiplicative about a fitted origin, and the two multiplicative
#   models), compares each lead's lack-of-fit statistics Q_B - Q_A (the difference that both
#   the AIC-like and BIC-like scores rank by). The null distribution of that difference is a
#   generalised chi-square, obtained by simulation from the pooled Stage 1 covariance under
#   a pure model-A and a pure model-B profile with the lead's own mean effect. This gives,
#   per lead, how surprising the observed difference is under each model (used for the share
#   of leads at which a purely additive profile is rejected) and the power of a size-5% test.
# Paper: Fig. 4 (legend: additive rejected in favour of multiplicative at 89% of leads, the
#   reverse at one lead)
# Input: lead profiles and pooled covariance (results/01_data/lead_profiles); model weights
#   and scores (results/05_bmi/model_scores)
# Output: results/05_bmi/model_separability/{lead_pair_power,pair_power_summary}.csv
# Runs: locally from results/ in about 6 minutes

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
source(file.path(repo, "analysis", "05_bmi", "bmi_common.R"))
cfg <- read_config()
u <- tau_grid(cfg); n_tau <- length(u); trap_w <- trap_weights(u)
N_SIM <- 20000L; SEED <- 20260928L
OUT <- "05_bmi/model_separability"

profiles <- read_lead_profiles()
weight_table <- as.data.table(read_result("05_bmi/model_scores", "weight_matrices"))
lead_scores <- as.data.table(read_result("05_bmi/model_scores", "lead_scores"))
MODELS <- c("additive", "log_additive", "shifted", "both", "logscale_vqtl", "shifted_lognormal_full", "three_components")
PAIRS <- list(c("additive", "log_additive"), c("additive", "shifted"), c("log_additive", "shifted"))

# Weight matrix of each model on one scale, rebuilt from the model-scores table.
weights_for_scale <- function(sc) {
  x <- weight_table[trait == "bmi" & scale == sc]
  wcols <- c("additive", "multiplicative_c0", "multiplicative", "meanlog", "sdlog", "level", "tilt", "bend")
  setNames(lapply(MODELS, function(m) {
    y <- x[model == m][order(as.numeric(tau))]
    cols <- wcols[vapply(wcols, function(k) all(!is.na(y[[k]])), TRUE)]
    W <- as.matrix(y[, cols, with = FALSE]); storage.mode(W) <- "double"; W
  }), MODELS)
}
# Q_M = d beta' K_M beta with K_M = M' Sigma0^-1 M and M = I - H_M the OLS residual maker.
residual_form <- function(X, Sigma0_inv) {
  M <- diag(n_tau) - X %*% solve(crossprod(X), t(X)); t(M) %*% Sigma0_inv %*% M
}
quad_form <- function(B, K) rowSums((B %*% K) * B)

## 1. Shared noise draws -----------------------------------------------------------------
# One set of standard-normal draws, turned into N(0, Sigma0) noise per chromosome and shared
# across leads and model pairs (common random numbers), so that power curves are smooth.
set.seed(SEED)
noise <- matrix(rnorm(N_SIM * n_tau), N_SIM, n_tau)

## 2. Per lead and model pair ---------------------------------------------------------------
# With beta = s mu + e, e ~ N(0, Sigma0 / d), the statistic is
#   dQ(s) = Q_B - Q_A = s^2 a + 2 s b + C,
# where a = d mu' D mu, b = sqrt(d) e' D mu, C = e' D e and D = K_B - K_A. The truth under A is
# the pure model-A profile with the lead's observed mean effect (trapezoid functional of its
# profile), and likewise for B. The size-5% test of B against A uses the 95th percentile of dQ
# under the pure-B truth; its power is taken under the pure-A truth, and z80 is the mean-effect
# z at which that power reaches 80% (amplitude scaled by bisection, shape and d fixed).
# The observed dQ is placed in both null distributions (pA_obs, pB_obs).
lead_rows <- list(); q_check <- list()
for (sc in c("raw", "log")) {
  leads <- read_lead_set(profiles, if (sc == "raw") "bmi_A_regenie" else "log_bmi_A_regenie", u)
  Sigma0 <- read_lead_sigma0(leads$source, u)
  meta <- leads$meta; B <- leads$B; W <- weights_for_scale(sc)
  meta[, delta_mu := as.numeric(B %*% trap_w)]
  for (k in unique(as.character(meta$chr))) {
    i <- which(as.character(meta$chr) == k); d <- meta$d_resid[i]
    Sk <- Sigma0[[k]]; Si <- solve(Sk); Si <- (Si + t(Si)) / 2
    noise_k <- noise %*% chol(Sk)                                     # rows ~ N(0, Sigma0)
    Sk_inv_half <- { e <- eigen(Sk, symmetric = TRUE); e$vectors %*% diag(1 / sqrt(e$values)) %*% t(e$vectors) }
    K <- lapply(W, residual_form, Sigma0_inv = Si)
    Q <- sapply(MODELS, function(m) d * quad_form(B[i, , drop = FALSE], K[[m]]))
    if (length(i) == 1L) Q <- matrix(Q, 1, dimnames = list(NULL, MODELS))
    q_check[[length(q_check) + 1L]] <- data.table(scale = sc, snp_id = meta$snp_id[i], Q)
    for (pair in PAIRS) {
      A <- pair[1]; Bm <- pair[2]; wA <- as.numeric(W[[A]]); wB <- as.numeric(W[[Bm]])
      D <- K[[Bm]] - K[[A]]; noise_D <- noise_k %*% D; C <- rowSums(noise_D * noise_k)
      # angle between the two models in this chromosome's whitened metric (d cancels)
      yA <- Sk_inv_half %*% wA; yB <- Sk_inv_half %*% wB
      sin_angle <- sqrt(max(0, 1 - sum(yA * yB)^2 / (sum(yA^2) * sum(yB^2))))
      for (j in seq_along(i)) {
        ii <- i[j]; dm <- meta$delta_mu[ii]; abs_z <- abs(meta$delta_mu_z[ii])
        muA <- wA * dm / sum(trap_w * wA); muB <- wB * dm / sum(trap_w * wB)
        aA <- d[j] * sum(muA * (D %*% muA)); aB <- d[j] * sum(muB * (D %*% muB))
        bA <- sqrt(d[j]) * as.numeric(noise_D %*% muA); bB <- sqrt(d[j]) * as.numeric(noise_D %*% muB)
        dQ <- function(s, a, b) s^2 * a + 2 * s * b + C
        power <- function(s, decide) {            # decide "A": dQ large; decide "B": dQ small
          qa <- dQ(s, aA, bA); qb <- dQ(s, aB, bB)
          if (decide == "A") mean(qa > quantile(qb, 0.95, names = FALSE)) else mean(qb < quantile(qa, 0.05, names = FALSE))
        }
        scale_for_80 <- function(decide) {
          if (power(100, decide) < 0.8) return(Inf)
          lo <- log(1e-3); hi <- log(100)
          for (it in 1:30) { mid <- (lo + hi) / 2; if (power(exp(mid), decide) >= 0.8) hi <- mid else lo <- mid }
          exp(hi)
        }
        qa1 <- dQ(1, aA, bA); qb1 <- dQ(1, aB, bB); observed <- Q[j, Bm] - Q[j, A]
        zw_A <- sqrt(d[j]) * sqrt(sum((Sk_inv_half %*% muA)^2))     # whitened length of the pure-A truth
        zw_B <- sqrt(d[j]) * sqrt(sum((Sk_inv_half %*% muB)^2))
        lead_rows[[length(lead_rows) + 1L]] <- data.table(scale = sc, snp_id = meta$snp_id[ii], chr = meta$chr[ii],
          delta_mu_z = meta$delta_mu_z[ii], model_a = A, model_b = Bm,
          lambda_truthA = d[j] * sum(muA * (K[[Bm]] %*% muA)), lambda_truthB = d[j] * sum(muB * (K[[A]] %*% muB)),
          mean_dQ_truthA = mean(qa1), mean_dQ_truthB = mean(qb1), tr_DSigma = mean(C),
          p_correct_truthA = mean(qa1 > 0), p_correct_truthB = mean(qb1 < 0),
          p_nontie_correct_truthA = mean(qa1 > 2), p_nontie_correct_truthB = mean(qb1 < -2),
          power_test_truthA = power(1, "A"), power_test_truthB = power(1, "B"),
          z80_test_truthA = scale_for_80("A") * abs_z, z80_test_truthB = scale_for_80("B") * abs_z,
          whitened_norm_truthA = zw_A, whitened_norm_truthB = zw_B,
          angle_whitened_deg = asin(min(1, sin_angle)) * 180 / pi,
          # the targeted 1-df approximation 2.80 / sin(angle) of 01_model_profiles.R, for comparison
          z80_P10_whitened = 2.8016 / sin_angle, z80_P10_dmz_truthA = 2.8016 / sin_angle * abs_z / zw_A,
          z80_P10_dmz_truthB = 2.8016 / sin_angle * abs_z / zw_B,
          dQ_obs = observed, pA_obs = mean(qa1 <= observed), pB_obs = mean(qb1 >= observed))
      }
    }
  }
}
lead_pair_power <- rbindlist(lead_rows)

# The Q statistics recomputed here must equal the model scores.
q_here <- melt(rbindlist(q_check), id.vars = c("scale", "snp_id"), variable.name = "model", value.name = "Q_here")
q_here[, model := as.character(model)]
q_here <- merge(q_here, lead_scores[set == "primary", .(scale, snp_id, model, Q_full)], by = c("scale", "snp_id", "model"))
stopifnot(nrow(q_here) == length(MODELS) * nrow(lead_pair_power) / length(PAIRS),
          max(abs(q_here$Q_here - q_here$Q_full) / pmax(q_here$Q_full, 1e-12)) < 1e-6)

## 3. Summary over leads ---------------------------------------------------------------------
pair_power_summary <- lead_pair_power[, .(n_leads = .N,
  median_dmz = median(abs(delta_mu_z)),
  frac_power80_test_truthA = mean(power_test_truthA >= 0.8), frac_power80_test_truthB = mean(power_test_truthB >= 0.8),
  frac_pcorrect80_truthA = mean(p_correct_truthA >= 0.8), frac_pcorrect80_truthB = mean(p_correct_truthB >= 0.8),
  frac_pnontie80_truthA = mean(p_nontie_correct_truthA >= 0.8), frac_pnontie80_truthB = mean(p_nontie_correct_truthB >= 0.8),
  median_power_test_truthA = median(power_test_truthA), median_power_test_truthB = median(power_test_truthB),
  median_z80_test_truthA = median(z80_test_truthA), median_z80_test_truthB = median(z80_test_truthB),
  median_z80_P10_dmz_truthA = median(z80_P10_dmz_truthA), median_z80_P10_whitened = median(z80_P10_whitened),
  median_angle_whitened_deg = median(angle_whitened_deg),
  median_ratio_whitened_to_dmz = median(whitened_norm_truthA / abs(delta_mu_z)),
  median_tr_DSigma = median(tr_DSigma),
  frac_obs_rejects_B = mean(pB_obs < 0.05), frac_obs_rejects_A = mean(pA_obs < 0.05),
  frac_obs_rejects_neither = mean(pA_obs >= 0.05 & pB_obs >= 0.05), frac_obs_rejects_both = mean(pA_obs < 0.05 & pB_obs < 0.05),
  frac_obs_absdQ_ge2 = mean(abs(dQ_obs) >= 2), median_obs_dQ = median(dQ_obs),
  q10_obs_dQ = quantile(dQ_obs, 0.1), q90_obs_dQ = quantile(dQ_obs, 0.9)), by = .(scale, model_a, model_b)]

write_result(lead_pair_power, OUT, "lead_pair_power")
write_result(pair_power_summary, OUT, "pair_power_summary")
