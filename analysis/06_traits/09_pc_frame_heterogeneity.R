#!/usr/bin/env Rscript
# Between-lead heterogeneity measured in the principal-component frame, which no change of scale can alter
#
# What: The fixed-basis heterogeneity tau_h (spread against location, 01 and 04) depends on the
#   scale. It moves by 30-60% between raw and log for strongly skewed traits because its basis is
#   anchored on a "pure shift", and it rejects for every trait. Here the same model
#   (S_i = rho_i L_i + e_i, rho_i ~ (r, tau^2)) is fitted with L and S each lead's PC1 and PC2
#   scores in the whitened metric y = sqrt(d) Sigma0^-1/2 beta, which a change of scale only
#   rotates. The PC directions are the top two eigenvectors of the noise-corrected uncentred
#   second moment of the whitened profiles, estimated on one random half of the leads and used
#   to score the other half, then swapped (cross-fitting keeps a lead's own noise out of its
#   direction). One common ratio r per half; Q = Q_1 + Q_2 ~ chi2(n - 2). Because the two halves'
#   Q are slightly positively correlated, the reference distribution is also simulated per trait
#   (one shared shape, tau = 0) to give a calibrated P; the detection floor (80% power at 5%)
#   and the attenuation of tau by noisy directions (recovery at the observed tau) are simulated
#   with the same statistic. Each trait is run on the raw, log and Box-Cox scales, by exact
#   conversion of its profiles, to show the PC-frame result does not move, and the fixed-basis
#   tau_h is recomputed on the same converted profiles for comparison.
# Paper: scale-invariant check of the between-lead heterogeneity in Fig. 5a (not drawn in a figure)
# Input: aligned lead-variant profiles and per-chromosome pooled covariance
#   (results/01_data/lead_profiles), baseline quantiles (results/01_data/baseline_quantiles),
#   results/06_traits/moments/{arms,pca_parallel,detection_floor}.csv,
#   results/06_traits/scale/{boxcox_lambda,scale_comparison}.csv,
#   results/06_traits/survey/fig6_components.csv, results/06_traits/pc_heterogeneity/pc2_share_ratio.csv
# Output: results/06_traits/pc_heterogeneity/{het_pc,het_pc_floor,het_pc_power_curve,
#   het_pc_comparison,het_pc_invariance,het_pc_checks}.csv
# Runs: locally from results/ in about an hour on one core (trait x scale jobs are independent
#   and each sets its own seeds, so they can be run in parallel)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
u <- tau_grid(cfg); n_tau <- length(u)
hp <- cfg$estimation$het_pc
seed <- as.integer(hp$split_seed)
n_boot <- as.integer(cfg$estimation$het_boot)          # 2,000
n_split <- as.integer(hp$n_split)                      # further random splits (split sensitivity)
n_sim <- as.integer(hp$floor_nsim)                     # simulations for the null and the floor
n_recovery <- as.integer(hp$recovery_nsim)
n_folds <- as.integer(hp$folds)                        # 2: halves
floor_grid <- seq(0, 0.5, by = as.numeric(hp$floor_grid_step))

lead_profiles <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
baseline <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))
baseline_quantiles <- function(trait_name) {
  x <- baseline[baseline$trait == trait_name]
  as.numeric(x$q[match(round(u, 3), round(as.numeric(x$tau), 3))])
}
arms <- as.data.table(read_result("06_traits/moments", "arms"))
parallel <- as.data.table(read_result("06_traits/moments", "pca_parallel"))
fixed_floor <- as.data.table(read_result("06_traits/moments", "detection_floor"))
powers <- as.data.table(read_result("06_traits/scale", "boxcox_lambda"))
fixed_by_scale <- as.data.table(read_result("06_traits/scale", "scale_comparison"))
pc2 <- as.data.table(read_result("06_traits/pc_heterogeneity", "pc2_share_ratio"))
traits <- setdiff(as.data.table(read_result("06_traits/survey", "fig6_components"))$trait,
                  unlist(hp$exclude_traits))             # the 21 traits of Fig. 5
stopifnot(length(traits) == 21L)

## 1. Each trait's profiles on the raw, log and Box-Cox scales
# Exact conversion through the raw scale: beta_T = T'(Q) beta, Sigma_T = D Sigma0 D. Box-Cox uses
# the normalising power of the raw trait (log BMI takes that of BMI).
jobs <- list()
for (tr in traits) {
  leads <- lead_set_matrix(lead_profiles, arms[arm == tr]$lead_set, u)
  is_log <- grepl("^log_", tr)
  q_raw <- if (is_log) exp(baseline_quantiles(tr)) else baseline_quantiles(tr)
  B_raw <- if (is_log) sweep(leads$B, 2, q_raw, `*`) else leads$B
  sigma0 <- p05_sigma0(leads$source, u)
  sigma0_raw <- if (is_log) lapply(sigma0, function(S) change_scale(matrix(0, 1, n_tau), S, q_raw)$Sigma0) else sigma0
  lambda_bc <- powers[trait == ifelse(tr == "log_bmi", "bmi", tr)]$lambda_main
  for (sc in c("raw", "log", "boxcox")) {
    lambda <- switch(sc, raw = 1, log = 0, boxcox = lambda_bc)
    if (length(lambda) != 1L || !is.finite(lambda) || (lambda != 1 && any(q_raw <= 0))) next
    D <- q_raw^(lambda - 1)
    jobs[[length(jobs) + 1L]] <- list(trait = tr, scale = sc, lambda = lambda,
      analysed = sc == (if (is_log) "log" else "raw"), B = sweep(B_raw, 2, D, `*`),
      S0 = lapply(sigma0_raw, function(S) change_scale(matrix(0, 1, n_tau), S, D)$Sigma0),
      d = leads$meta$d_resid, chr = as.character(leads$meta$chr))
  }
}

## 2. One trait on one scale
one_job <- function(jb) {
  # PC-frame test: cross-fitted halves, bootstrap within halves with the scores held fixed
  # (re-estimating directions per resample, *_reest, attenuates tau), and n_split further splits.
  f <- het_test_pc(jb$B, jb$S0, jb$d, jb$chr, seed = seed, nboot = n_boot, n_split = n_split, K = n_folds)
  # Fixed-basis comparator on the same converted profiles: unit-SE shift/spread components.
  Z <- matrix(NA_real_, nrow(jb$B), 3)
  for (k in unique(jb$chr)) {
    i <- which(jb$chr == k)
    Z[i, ] <- unit_se_components(jb$B[i, , drop = FALSE], jb$S0[[k]], jb$d[i], u)$Z
  }
  h <- het_test(Z[, 1], Z[, 2], excess = "n_minus_1")
  hb <- het_boot(Z[, 1], Z[, 2], B = n_boot, seed = as.integer(cfg$estimation$boot_seed))
  row <- data.table(trait = jb$trait, scale = jb$scale, lambda = jb$lambda, analysed_scale = jb$analysed,
    n_leads = f$n, whitening_rank = f$rank,
    tau_pc = f$tau, tau_pc_lo95 = f$tau_ci[1], tau_pc_hi95 = f$tau_ci[2], tau_pc_boot_median = median(f$boot[, "tau"]),
    tau_pc_lo95_reest = f$tau_ci_reest[1], tau_pc_hi95_reest = f$tau_ci_reest[2],
    tau_pc_boot_median_reest = median(f$boot_reest[, "tau"]),
    Q_pc = f$Q, df_pc = f$df, p_pc = f$p, log10_p_pc = pchisq(f$Q, f$df, lower.tail = FALSE, log.p = TRUE) / log(10),
    r_fold1 = f$r_fold[1], r_fold2 = f$r_fold[2], tau_fold1 = f$tau_fold[1], tau_fold2 = f$tau_fold[2],
    converged = f$converged, n_resid_z_gt_3 = sum(abs(f$resid_z) > 3),
    split_tau_median = median(f$split[, "tau"]), split_tau_min = min(f$split[, "tau"]),
    split_tau_max = max(f$split[, "tau"]), split_p_median = median(f$split[, "p"]),
    split_p_max = max(f$split[, "p"]), split_frac_p_lt_0.05 = mean(f$split[, "p"] < 0.05),
    eig1 = f$full_dirs$values[1], eig2 = f$full_dirs$values[2], eig3 = f$full_dirs$values[3],
    tau_old = h$tau, tau_old_lo95 = hb$tau_ci[1], tau_old_hi95 = hb$tau_ci[2], tau_old_boot_median = median(hb$draws[, "tau"]),
    r_old = h$r, Q_old = h$Q, df_old = h$df, p_old = h$p,
    log10_p_old = pchisq(h$Q, h$df, lower.tail = FALSE, log.p = TRUE) / log(10), seconds = NA_real_)
  if (!jb$analysed) return(list(row = row))
  # The statistic does not depend on the scale, so the null, floor and recovery are simulated
  # once per trait, on the analysed scale.
  null <- het_pc_null(f, nsim = n_sim, seed = seed + 3L)
  fl <- het_pc_detection_floor(f, grid = floor_grid, nsim = n_sim, seed = seed, q_crit = null$q_crit)
  power_at <- function(t) { i <- which(abs(fl$curve$tau - t) < 1e-9); if (length(i)) fl$curve$power[i] else NA_real_ }
  rec <- if (f$tau > 0) het_pc_recovery(f, tau = f$tau, nsim = n_recovery, seed = seed) else
    list(mean = NA_real_, sd = NA_real_, ratio = NA_real_)
  floor_row <- data.table(trait = jb$trait, scale = jb$scale, n_leads = f$n, alpha = 0.05, target_power = 0.8,
    n_sim = n_sim, r0 = fl$r0, q_crit_calibrated = null$q_crit, q_crit_chi2 = qchisq(0.95, f$df),
    null_mean_Q = null$mean, null_sd_Q = null$sd, null_dispersion = null$dispersion, size_of_chi2_test = null$size_chi2,
    p_pc_empirical = null$p_empirical, p_pc_calibrated = null$p_calibrated, log10_p_pc_calibrated = null$log10_p_calibrated,
    tau_pc_floor_80 = fl$floor, tau_pc_observed = f$tau, observed_above_floor = is.finite(fl$floor) && f$tau >= fl$floor,
    power_at_0 = power_at(0), power_at_0.05 = power_at(0.05), power_at_0.10 = power_at(0.10),
    recovery_nsim = n_recovery, recovery_mean_at_observed = rec$mean, recovery_sd_at_observed = rec$sd,
    recovery_ratio = rec$ratio, tau_pc_attenuation_corrected = f$tau / rec$ratio)
  list(row = row, floor = floor_row, curve = data.table(trait = jb$trait, tau = fl$curve$tau, power = fl$curve$power))
}
results <- lapply(jobs, function(jb) {
  t0 <- proc.time()[["elapsed"]]
  r <- one_job(jb)
  r$row[, seconds := proc.time()[["elapsed"]] - t0]
  message(jb$trait, " (", jb$scale, "): tau_pc ", signif(r$row$tau_pc, 3), ", p ", signif(r$row$p_pc, 3))
  r
})
het <- rbindlist(lapply(results, `[[`, "row"))[order(match(trait, traits), match(scale, c("raw", "log", "boxcox")))]
floors <- rbindlist(lapply(results, `[[`, "floor"))
curves <- rbindlist(lapply(results, `[[`, "curve"))
n_by_trait <- het[analysed_scale == TRUE, .(trait, n_leads)]
curve_order <- n_by_trait[order(-n_leads, match(trait, traits))]$trait        # largest traits first
curves <- curves[order(match(trait, curve_order))]

## 3. PC-frame against fixed-basis results
floors <- merge(floors, fixed_floor[, .(trait = arm, tau_old_floor_80 = tau_floor_80)], by = "trait", all.x = TRUE)[order(match(trait, traits))]
comparison <- merge(het[, .(trait, scale, lambda, n_leads, tau_pc, tau_pc_lo95, tau_pc_hi95, p_pc, log10_p_pc,
                            tau_old, tau_old_lo95, tau_old_hi95, p_old, log10_p_old)],
                    parallel[k_assumed == 1L, .(trait = arm, pc2_parallel_p = p_exceeds_null, n_components)],
                    by = "trait", all.x = TRUE)
comparison <- merge(comparison, floors[, .(trait, p_pc_calibrated, log10_p_pc_calibrated, p_pc_empirical, null_dispersion,
                                           recovery_ratio, tau_pc_floor_80, tau_old_floor_80)], by = "trait", all.x = TRUE)
comparison[, `:=`(pc_above_floor = tau_pc >= tau_pc_floor_80, pc_reject_0.05 = p_pc_calibrated < 0.05,
                  old_reject_0.05 = p_old < 0.05, pc2_parallel_reject_0.05 = pc2_parallel_p < 0.05)]
comparison <- comparison[order(match(trait, traits), match(scale, c("raw", "log", "boxcox")))]

## 4. Checks: invariance of the PC-frame result across scales, reproduction of earlier runs
invariance <- het[, .(max_abs_dev_tau_pc = max(abs(tau_pc - tau_pc[analysed_scale])),
                      max_rel_dev_Q_pc = max(abs(Q_pc - Q_pc[analysed_scale]) / Q_pc[analysed_scale]),
                      max_abs_dev_tau_pc_ci = max(abs(c(tau_pc_lo95 - tau_pc_lo95[analysed_scale],
                                                        tau_pc_hi95 - tau_pc_hi95[analysed_scale]))),
                      range_tau_old = max(tau_old) - min(tau_old),
                      rel_range_tau_old = (max(tau_old) - min(tau_old)) / tau_old[analysed_scale]), by = trait]
# The scale comparison converts raw BMI; here log BMI is converted (its own Stage 1, so the
# kernel density differs slightly).
vs_scale <- merge(het[, .(trait, scale, tau_old)],
                  fixed_by_scale[, .(trait = ifelse(trait == "bmi", "log_bmi", trait), scale, tau_scale = tau_unit_se)],
                  by = c("trait", "scale"))
vs_scale[, abs_diff := abs(tau_old - tau_scale)]
checks <- data.table(
  check = c("max over traits/scales |tau_pc - tau_pc(analysed)|", "max relative |Q_pc - Q_pc(analysed)|",
            "max |tau_pc CI - CI(analysed)|", "max |tau_old - scale comparison tau_unit_se| (same Stage 1)",
            "max |tau_old - scale comparison tau_unit_se| (log_bmi vs bmi rows)",
            "max relative range of tau_old over scales", "all fold IRLS converged",
            "max |PC2 share - moments observed_pct|", "max |null median recomputed - moments| / moments",
            "PC2 share ratio inside its raw fixed-null percentile CI (number of traits)",
            "PC2 share ratio inside its duplicated-noise-null CI (number of traits)",
            "median bootstrap shift of the ratio (fixed null)"),
  value = c(max(invariance$max_abs_dev_tau_pc), max(invariance$max_rel_dev_Q_pc), max(invariance$max_abs_dev_tau_pc_ci),
            max(vs_scale[trait != "log_bmi"]$abs_diff), max(vs_scale[trait == "log_bmi"]$abs_diff),
            max(invariance$rel_range_tau_old), as.numeric(all(het$converged)),
            max(abs(pc2$observed_pc2_pct - pc2$observed_pc2_pct_p11)),
            max(abs(pc2$null_median_recomputed_pct - pc2$null_median_pct) / pc2$null_median_pct),
            sum(pc2$ratio >= pc2$lo95_fixed_null_raw & pc2$ratio <= pc2$hi95_fixed_null_raw),
            sum(pc2$ratio >= pc2$lo95_dup_null & pc2$ratio <= pc2$hi95_dup_null), median(pc2$boot_shift)))
print(checks)
print(comparison[, .(trait, scale, tau_pc, p_pc_calibrated, tau_old, p_old, pc2_parallel_p)])

## 5. Write
out <- "06_traits/pc_heterogeneity"
write_result(het, out, "het_pc")
write_result(floors, out, "het_pc_floor")
write_result(curves, out, "het_pc_power_curve")
write_result(comparison, out, "het_pc_comparison")
write_result(invariance, out, "het_pc_invariance")
write_result(checks, out, "het_pc_checks")
