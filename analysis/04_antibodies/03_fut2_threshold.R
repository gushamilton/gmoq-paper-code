#!/usr/bin/env Rscript
# Why a level shift looks like a serostatus effect after thresholding (FUT2)
#
# What: at the JCV FUT2-region lead rs681343 G-MoQ estimates a shift in the level of both
#   mixture components with little effect on membership, yet the thresholded REGENIE analyses
#   find a much stronger association with seropositivity than with IgG among seropositives.
#   This script asks whether a pure level shift plus the assay threshold (MFI > 250, c = ln 250)
#   explains that, for the example loci (FUT2 and two HLA variants) in JCV and MCV:
#   (1) where the threshold falls: its quantile level tau_c, the density there, agreement of
#       the assay flag with Y > c;
#   (2) first-order predictions from a quantile-effect profile beta(u):
#         dP(Y > c)/dg = f(c) beta(tau_c),
#         dmu_+/dg     = mean of beta(u) over u > tau_c - (mu_+ - c) f(c) beta(tau_c) / (1 - tau_c),
#       for the Stage 1 profile, the G-MoQ fit and its parts, a pure level shift and a pure
#       membership effect with the same mean effect, converted to expected Wald z statistics;
#   (3) a semi-synthetic check of (2) on the observed antibody levels with simulated genotypes;
#   (4) a simulation from the fitted mixture without linearisation: participants assigned to
#       components by their mixture posterior, genotype effects applied to membership
#       (logit P(upper) for the simulated dpi_upper; the code's gamma is logit P(lower)) and to
#       each component's level, under the G-MoQ estimates, a pure shift and pure membership
#       (300 replicates each), each analysed by the two thresholded regressions.
# Paper: Supplementary Note 3 (Supplementary Note Table 2)
# Input: UK Biobank antibody phenotypes, serostatus and covariates of the analysed cohort
#   (02_serostatus_gwas.R); Stage 1 kernel density (rif_nuisance.json); results of
#   01_decompose_antibody_associations.R and 02_serostatus_gwas.R
# Output: results/04_antibodies/fut2_threshold/{threshold_summary, predictions, semisynthetic,
#   semisynthetic_mixture, mixture_split}.csv
# Runs: on a cluster with UK Biobank data, in about 10 minutes

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
pb <- cfg$p09b
TRAITS <- unlist(cfg$headline$antibody_traits)
NSIM <- 300L
SEED <- 20260928L
c_thr <- log(pb$seropositivity_threshold)                  # ln 250 = 5.52 log MFI
work <- private_dir(cfg, "04_antibodies", "serostatus")
out_dir <- "04_antibodies/fut2_threshold"
prof <- as.data.table(read_result("04_antibodies/decomposition", "example_profiles"))
dec <- as.data.table(read_result("04_antibodies/decomposition", "lead_decomposition"))
sero <- as.data.table(read_result("04_antibodies/serostatus_gwas", "serostatus_leads"))
COMP <- as.data.table(read_result("04_antibodies/decomposition", "mixture_components"))

## 1. The analysed cohort and where the threshold falls ------------------------------------
ph <- fread(pb$phenotype_file); sr <- fread(file.path(work, "phenotypes_serostatus.txt")); cv <- fread(pb$covar_file)
kp <- fread(file.path(work, "keep_analysed.txt"), header = FALSE)
for (d in list(ph, sr, cv)) d[, IID := as.character(IID)]
ids <- as.character(kp[[2]])
covariates <- c("age", "sex", paste0("PC", 1:10))                 # as the REGENIE comparators
X <- as.matrix(cv[match(ids, IID), ..covariates])
if (anyNA(X)) stop("covariates missing in the analysed cohort")
resid_sd <- function(y, X) { f <- lm.fit(cbind(1, X), y); sqrt(sum(f$residuals^2) / (length(y) - ncol(X) - 1)) }

SUM <- list(); DIST <- list()
for (tr in TRAITS) {
  stem <- sub("_log_all$", "", tr)
  y <- ph[match(ids, IID)][[tr]]; s <- sr[match(ids, IID)][[paste0(stem, "_seropos")]]
  if (anyNA(y) || anyNA(s) || length(y) != cfg$cohorts[[tr]]$n) stop(tr, ": cohort mismatch")
  rn <- read_json(file.path(cfg$cohorts[[tr]]$root, "chr19", "rif_nuisance.json"))
  tau_c <- mean(y <= c_thr); p <- mean(s == 1)
  f_rif <- approx(rn$q_tau, rn$fhat, xout = c_thr)$y        # Stage 1 kernel density at c
  f_kde <- { d <- density(y, n = 4096); approx(d$x, d$y, xout = c_thr)$y }
  pos <- s == 1; above <- y > c_thr
  SUM[[tr]] <- data.table(
    trait = tr, n = length(y), n_seropositive_flag = sum(pos), n_y_gt_c = sum(above),
    agreement_flag_vs_threshold = mean(pos == above), c_log_mfi = c_thr,
    tau_c_threshold = tau_c, tau_c_flag = 1 - p, p_seropositive = p,
    f_c_rif = f_rif, f_c_kde = f_kde, f_mode_rif = max(rn$fhat),
    q_mode_rif = rn$q_tau[which.max(rn$fhat)], tau_mode_rif = tau_grid(cfg)[which.max(rn$fhat)],
    mean_y = mean(y), sd_y = sd(y), median_y = median(y),
    mean_y_seropos_flag = mean(y[pos]), mean_y_gt_c = mean(y[above]),
    sd_y_seropos_flag = sd(y[pos]), resid_sd_y_seropos_flag = resid_sd(y[pos], X[pos, ]),
    resid_sd_seropos_binary = resid_sd(as.numeric(pos), X),
    trunc_attenuation_pure_shift = 1 - (mean(y[above]) - c_thr) * f_rif / (1 - tau_c))
  DIST[[tr]] <- list(y = y, u = (frank(y, ties.method = "average") - 0.5) / length(y))
}
SUM <- rbindlist(SUM)

## 2. First-order predictions from a quantile-effect profile -------------------------------
profile_fun <- function(tau, b) function(u) approx(tau, b, xout = u, rule = 2)$y
predict_thresholded <- function(bf, S, q_all, q_pos, n_all, n_pos) {
  tc <- S$tau_c_threshold; fc <- S$f_c_rif; m <- S$mean_y_gt_c; p <- S$p_seropositive
  uu <- seq(tc, 1, length.out = 4001); w <- trap_weights(uu) / (1 - tc)
  within <- sum(w * bf(uu)); selection <- -(m - S$c_log_mfi) * fc * bf(tc) / (1 - tc)
  dP <- fc * bf(tc); dm <- within + selection
  data.table(beta_at_tau_c = bf(tc), mean_effect = sum(trap_weights(seq(0, 1, length.out = 4001)) * bf(seq(0, 1, length.out = 4001))),
             pred_dP = dP, pred_logodds = dP / (p * (1 - p)),
             pred_z_sero = dP / (S$resid_sd_seropos_binary / sqrt(n_all * 2 * q_all * (1 - q_all))),
             pred_igg_within = within, pred_igg_selection = selection, pred_igg_beta = dm,
             pred_z_igg = dm / (S$resid_sd_y_seropos_flag / sqrt(n_pos * 2 * q_pos * (1 - q_pos))))
}
z_from <- function(beta, mlog10p) sign(beta) * qnorm(-mlog10p * log(10) - log(2), lower.tail = FALSE, log.p = TRUE)

PRED <- list(); SIM_INPUT <- list()
keys <- unique(prof[, .(weight_design, trait, snp_id)])
for (k in seq_len(nrow(keys))) {
  key <- keys[k]
  P <- prof[key, on = .(weight_design, trait, snp_id)][order(tau)]
  D <- dec[key, on = .(weight_design, trait, snp_id)]
  R1 <- sero[trait == key$trait & snp_id == key$snp_id & analysis == "seropositivity"]
  R2 <- sero[trait == key$trait & snp_id == key$snp_id & analysis == "seropositive_igg"]
  if (nrow(D) != 1 || nrow(R1) != 1 || nrow(R2) != 1) stop("join failed: ", paste(key, collapse = " "))
  if (D$effect_allele != R1$effect_allele || D$effect_allele != R2$effect_allele) stop("allele mismatch ", key$snp_id)
  S <- SUM[trait == key$trait]
  gamma_shape <- P$fitted_gamma / D$gamma
  mean_of <- function(b) sum(trap_weights(seq(0, 1, length.out = 4001)) * profile_fun(P$tau, b)(seq(0, 1, length.out = 4001)))
  profiles <- list(
    stage1_raw = P$beta_stage1, gmoq_fitted_total = P$fitted_total,
    gmoq_part_gamma = P$fitted_gamma, gmoq_part_locations = P$fitted_beta_left + P$fitted_beta_right,
    alt_pure_level_shift = rep(D$mean_effect, nrow(P)),
    alt_pure_membership = gamma_shape * D$mean_effect / mean_of(gamma_shape))
  for (nm in names(profiles)) {
    pr <- predict_thresholded(profile_fun(P$tau, profiles[[nm]]), S, R1$a1freq, R2$a1freq, R1$n, R2$n)
    PRED[[length(PRED) + 1]] <- cbind(key, profile = nm, effect_allele = D$effect_allele, pr)
    if (key$weight_design == "model" && nm %in% c("gmoq_fitted_total", "alt_pure_level_shift", "alt_pure_membership"))
      SIM_INPUT[[length(SIM_INPUT) + 1]] <- list(key = key, profile = nm, tau = P$tau, b = profiles[[nm]], q = R1$a1freq)
  }
}
PRED <- rbindlist(PRED)
OBS <- sero[, .(trait, snp_id, analysis, beta, se, mlog10p, n, a1freq)][, z := z_from(beta, mlog10p)]
OBS <- dcast(OBS, trait + snp_id ~ analysis, value.var = c("beta", "se", "z", "n"))
PRED <- merge(PRED, OBS, by = c("trait", "snp_id"))
PRED[, `:=`(ratio_pred_obs_z_sero = pred_z_sero / z_seropositivity,
            ratio_pred_obs_z_igg = pred_z_igg / z_seropositive_igg,
            pred_z_igg_over_z_sero = pred_z_igg / pred_z_sero,
            obs_z_igg_over_z_sero = z_seropositive_igg / z_seropositivity)]
# The Stage 1 quantile effect at tau_c is the effect on the indicator 1{Y > q(tau_c)}, rescaled.
S1Z <- prof[, {
  tc <- SUM[trait == .BY$trait]$tau_c_threshold
  .(stage1_beta_tau_c = approx(tau, beta_stage1, tc)$y, stage1_z_tau_c = approx(tau, beta_stage1 / se_stage1, tc)$y)
}, by = .(weight_design, trait, snp_id)]
PRED <- merge(PRED, S1Z, by = c("weight_design", "trait", "snp_id"))

## 3. Semi-synthetic check: observed levels, simulated genotypes ---------------------------
# Y' = Y + beta(rank(Y)) g with g ~ Bin(2, q); then logistic regression of Y' > c on g, and
# OLS of Y' on g among Y' > c. The random-number stream continues into section 4.
set.seed(SEED)
SIM <- rbindlist(lapply(SIM_INPUT, function(si) {
  Dd <- DIST[[si$key$trait]]; bu <- approx(si$tau, si$b, xout = Dd$u, rule = 2)$y; n <- length(Dd$y)
  r <- t(vapply(seq_len(NSIM), function(i) {
    g <- rbinom(n, 2, si$q); y2 <- Dd$y + bu * g; s2 <- as.numeric(y2 > c_thr)
    Xg <- cbind(1, g); fl <- glm.fit(Xg, s2, family = binomial()); mu <- fl$fitted.values
    V <- solve(crossprod(Xg * sqrt(mu * (1 - mu))))
    pos <- s2 == 1; fo <- lm.fit(cbind(1, g[pos]), y2[pos])
    se <- sqrt(sum(fo$residuals^2) / (sum(pos) - 2) / sum((g[pos] - mean(g[pos]))^2))
    c(fl$coefficients[2], fl$coefficients[2] / sqrt(V[2, 2]), fo$coefficients[2], fo$coefficients[2] / se, mean(pos))
  }, numeric(5)))
  cbind(si$key, profile = si$profile, nsim = NSIM,
        sim_logodds_mean = mean(r[, 1]), sim_z_sero_mean = mean(r[, 2]), sim_z_sero_sd = sd(r[, 2]),
        sim_igg_beta_mean = mean(r[, 3]), sim_z_igg_mean = mean(r[, 4]), sim_z_igg_sd = sd(r[, 4]),
        sim_p_seropos_mean = mean(r[, 5]))
}))

## 4. Simulation from the fitted mixture, without linearisation ---------------------------
# Each participant's observed level is split between the components by the posterior
# w_i = P(upper | y_i); r_l and r_u are the resulting empirical component distributions, so at
# baseline the simulation reproduces the observed levels. For genotype g ~ Bin(2, q), the
# component is upper with probability pi0 + dpi (g - 2q), the level is drawn from that
# component and shifted by beta_k g. Scenarios with the same mean effect:
#   gmoq            beta_left, beta_right and dpi from the fitted membership part
#                   (mean of fitted_gamma = dpi (mu_u - mu_l));
#   pure_shift      beta_l = beta_u = mean effect, dpi = 0;
#   pure_membership beta = 0, dpi = mean effect / (mu_u - mu_l).
u01 <- seq(0, 1, length.out = 4001); tw01 <- trap_weights(u01)
MIXSIM <- list(); MIXPAR <- list()
for (tr in TRAITS) {
  y <- DIST[[tr]]$y; n <- length(y); cm <- COMP[trait == tr]
  lo <- cm[component == "lower"]; up <- cm[component == "upper"]
  dl <- lo$weight * dsn_density(y, lo$xi, lo$omega, lo$alpha); du <- up$weight * dsn_density(y, up$xi, up$omega, up$alpha)
  w <- du / (dl + du); pi0 <- mean(w); rl <- (1 - w) / sum(1 - w); ru <- w / sum(w)
  mul <- sum(rl * y); muu <- sum(ru * y)
  MIXPAR[[tr]] <- data.table(trait = tr, pi0_upper = pi0, mu_lower = mul, mu_upper = muu,
                             upper_share_above_c = sum(ru[y > c_thr]), lower_share_above_c = sum(rl[y > c_thr]))
  for (sn in unique(prof[trait == tr & weight_design == "model"]$snp_id)) {
    P <- prof[weight_design == "model" & trait == tr & snp_id == sn][order(tau)]
    D <- dec[weight_design == "model" & trait == tr & snp_id == sn]
    q <- sero[trait == tr & snp_id == sn & analysis == "seropositivity"]$a1freq
    gm <- sum(tw01 * approx(P$tau, P$fitted_gamma, xout = u01, rule = 2)$y)
    scenarios <- list(gmoq = c(D$beta_left, D$beta_right, gm / (muu - mul)),
                      pure_shift = c(D$mean_effect, D$mean_effect, 0),
                      pure_membership = c(0, 0, D$mean_effect / (muu - mul)))
    for (nm in names(scenarios)) {
      th <- scenarios[[nm]]
      r <- t(vapply(seq_len(NSIM), function(i) {
        g <- rbinom(n, 2, q); pg <- pmin(pmax(pi0 + th[3] * (g - 2 * q), 0), 1); k <- runif(n) < pg
        yy <- numeric(n); yy[k] <- y[sample.int(n, sum(k), TRUE, ru)]; yy[!k] <- y[sample.int(n, sum(!k), TRUE, rl)]
        y2 <- yy + ifelse(k, th[2], th[1]) * g; s2 <- as.numeric(y2 > c_thr)
        Xg <- cbind(1, g); fl <- glm.fit(Xg, s2, family = binomial()); mu <- fl$fitted.values
        V <- solve(crossprod(Xg * sqrt(mu * (1 - mu))))
        pos <- s2 == 1; fo <- lm.fit(cbind(1, g[pos]), y2[pos])
        se <- sqrt(sum(fo$residuals^2) / (sum(pos) - 2) / sum((g[pos] - mean(g[pos]))^2))
        fm <- lm.fit(Xg, y2)$coefficients[2]
        c(fl$coefficients[2], fl$coefficients[2] / sqrt(V[2, 2]), fo$coefficients[2], fo$coefficients[2] / se, fm)
      }, numeric(5)))
      MIXSIM[[length(MIXSIM) + 1]] <- data.table(
        trait = tr, snp_id = sn, scenario = nm, beta_lower = th[1], beta_upper = th[2], dpi_upper = th[3], nsim = NSIM,
        sim_mean_effect = mean(r[, 5]), sim_logodds_mean = mean(r[, 1]), sim_z_sero_mean = mean(r[, 2]),
        sim_z_sero_sd = sd(r[, 2]), sim_igg_beta_mean = mean(r[, 3]), sim_z_igg_mean = mean(r[, 4]),
        sim_z_igg_sd = sd(r[, 4]))
    }
  }
}
MIXSIM <- merge(rbindlist(MIXSIM), OBS[, .(trait, snp_id, beta_seropositivity, z_seropositivity,
                                          beta_seropositive_igg, z_seropositive_igg)], by = c("trait", "snp_id"))

## 5. Write -------------------------------------------------------------------------------
write_result(rbindlist(MIXPAR), out_dir, "mixture_split")
write_result(MIXSIM, out_dir, "semisynthetic_mixture")
write_result(SUM, out_dir, "threshold_summary")
write_result(PRED, out_dir, "predictions")
write_result(SIM, out_dir, "semisynthetic")
