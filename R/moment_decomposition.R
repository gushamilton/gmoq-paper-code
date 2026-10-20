# The moment and principal-component reading of one set of lead profiles.
#
# One function, used for every trait and every measurement scale, that turns the Stage 1
# quantile-effect profiles of a trait's lead variants into the summaries reported in the
# paper:
#   1. linear functionals (change in mean, variance, L-moments, Cornish-Fisher level, tilt
#      and bend) with exact standard errors from the per-chromosome pooled Sigma0;
#   2. lack of fit of each profile to the shift + spread frame [1, Q - mu_b] (primary) and to
#      the normal-score frame [1, z, z^2 - 1], with exact Davies p-values;
#   3. counts of leads with significant spread, raw skew and skew beyond shift + spread
#      (Bonferroni within the lead set);
#   4. between-locus heterogeneity of shape relative to location (het_test()) on several
#      bases, with bootstrap intervals and a detection floor;
#   5. principal components of the whitened unit-length profiles, the number retained by
#      parallel analysis, and how well each loading is described by shift, spread and bend
#      bases (R^2 and angles in the quantile and whitened metrics).
# Requires data.table (attached by R/load.R).

# Acute angle in degrees between v and its projection on the column span of X (Euclidean).
angle_to_span <- function(v, X) {
  p <- X %*% qr.solve(X, v); cs <- sqrt(sum(p^2) / sum(v^2)); acos(min(1, cs)) * 180 / pi
}

# Arguments:
#   label: name of the lead set or analysis arm (written to every output table);
#   meta: data.frame with lead_set, snp_id, chr, bp, effect_allele, other_allele, d_resid;
#   B: leads x quantile levels, profiles aligned to the trait-increasing allele;
#   S0: named list of pooled Sigma0 by chromosome, optionally with "mean" (the whitening frame);
#   q: baseline quantiles on the analysed scale; u: quantile levels;
#   full: also compute bootstrap intervals and the detection floor (slow);
#   nboot, npa, nfloor: bootstrap resamples, parallel-analysis draws, floor simulations;
#   seed: RNG seed for every random step.
# Value: list of data.tables: lead (one row per lead), lof, skew, het, floor, pcurve, spec,
#   par, load, bf, ang, sig; plus n_components, mu_b, n_chr, n.
moment_decomposition <- function(label, meta, B, S0, q, u, full = TRUE, nboot = 2000L, npa = 300L,
                                 seed = 20260923L, nfloor = 1000L) {
  nt <- length(u); z <- stats::qnorm(u); cw <- trap_weights(u)
  d <- meta$d_resid; chr <- as.character(meta$chr); n <- nrow(B)
  if (!all(unique(chr) %in% names(S0))) stop(label, ": Sigma0 missing for some chromosomes")
  S0bar <- if ("mean" %in% names(S0)) S0[["mean"]] else Reduce(`+`, S0[unique(chr)]) / length(unique(chr))
  chrs <- unique(chr)
  byc <- split(seq_len(n), chr)
  mu_b <- sum(cw * q) / sum(cw)                       # baseline mean on the grid

  ## 1. Functionals with exact standard errors (Sigma0 of each lead's chromosome)
  W <- functional_weights(u, q)
  fn <- c("delta_mu", "delta_var", "delta_lambda2", "delta_lambda3", "skew_beyond_frame",
          "cf_level", "cf_tilt", "cf_bend", "mean_grid")
  Wm <- do.call(cbind, W[fn])
  EST <- SE <- matrix(NA_real_, n, length(fn), dimnames = list(NULL, fn))
  # (location, shape) pairs whose measurement covariances feed the heterogeneity tests
  vcov_pairs <- list(exact = c("delta_mu", "delta_var"), lmoment = c("delta_mu", "delta_lambda2"),
                     cornish_fisher = c("cf_level", "cf_tilt"))
  VC <- lapply(vcov_pairs, function(p) matrix(NA_real_, n, 3, dimnames = list(NULL, c("vLL", "vSS", "vLS"))))
  for (k in names(byc)) {
    i <- byc[[k]]
    fe <- functional_estimate(B[i, , drop = FALSE], S0[[k]], d[i], Wm)
    EST[i, ] <- fe$estimate; SE[i, ] <- fe$se
    for (pn in names(vcov_pairs)) {
      p <- vcov_pairs[[pn]]
      VC[[pn]][i, "vLL"] <- functional_cov(W[[p[1]]], W[[p[1]]], S0[[k]], d[i])
      VC[[pn]][i, "vSS"] <- functional_cov(W[[p[2]]], W[[p[2]]], S0[[k]], d[i])
      VC[[pn]][i, "vLS"] <- functional_cov(W[[p[1]]], W[[p[2]]], S0[[k]], d[i])
    }
  }
  Zf <- EST / SE

  ## 2. Lack of fit: shift + spread frame (primary) and normal-score frame, exact Davies
  frames <- list(shift_spread = cbind(1, q - mu_b), cf3 = cbind(1, z, z^2 - 1))
  LOF <- list(); LOFS <- list()
  for (fr in names(frames)) {
    out <- data.frame(Q = rep(NA_real_, n), p_davies = NA_real_, floored = NA, p_satt = NA_real_)
    lam1 <- NULL
    for (k in names(byc)) {
      i <- byc[[k]]
      lf <- lof_davies(B[i, , drop = FALSE], S0[[k]], d[i], frames[[fr]])
      out$Q[i] <- lf$Q; out$p_davies[i] <- lf$p_davies; out$floored[i] <- lf$davies_floored
      out$p_satt[i] <- lf$p_satterthwaite
      lam <- attr(lf, "lambda"); lam1 <- rbind(lam1, c(EQ = sum(lam), df_eff = sum(lam)^2 / sum(lam^2)))
    }
    LOF[[fr]] <- out
    bonf <- 0.05 / n
    LOFS[[fr]] <- data.table(arm = label, frame = fr, n_leads = n, df_nominal = nt - ncol(frames[[fr]]),
      E_Q_median_chr = median(lam1[, "EQ"]), df_eff_median_chr = median(lam1[, "df_eff"]),
      median_Q = median(out$Q), frac_p_lt_0.05 = mean(out$p_davies < 0.05),
      n_bonf_rejected = sum(out$p_davies < bonf), frac_adequate_bonf = mean(out$p_davies >= bonf),
      n_davies_floored = sum(out$floored), frac_satt_p_lt_0.05 = mean(out$p_satt < 0.05))
  }

  ## 3. Counts: spread, raw skew and skew beyond shift + spread (two-sided Bonferroni)
  zb <- qnorm(1 - 0.025 / n)
  sk <- Zf[, "skew_beyond_frame"]
  SKW <- data.table(arm = label, n_leads = n, z_bonferroni = zb,
    n_spread_lambda2_bonf = sum(abs(Zf[, "delta_lambda2"]) > zb),
    n_spread_exact_bonf = sum(abs(Zf[, "delta_var"]) > zb),
    n_spread_lambda2_nominal = sum(abs(Zf[, "delta_lambda2"]) > 1.96),
    n_skew_raw_bonf = sum(abs(Zf[, "delta_lambda3"]) > zb),
    n_skew_beyond_bonf = sum(abs(sk) > zb), n_skew_beyond_pos = sum(sk > zb), n_skew_beyond_neg = sum(sk < -zb),
    n_skew_beyond_nominal = sum(abs(sk) > 1.96), expected_nominal = 0.05 * n,
    mean_chi2_skew_beyond = mean(sk^2),
    median_spread_lambda2_z = median(Zf[, "delta_lambda2"]), median_spread_exact_z = median(Zf[, "delta_var"]),
    cor_spread_location = cor(EST[, "delta_mu"], EST[, "delta_lambda2"]),
    median_delta_mu_z = median(Zf[, "delta_mu"]))

  ## 4. Heterogeneity of shape relative to location
  # Unit-variance components c1 (level), c2 (tilt given level), c3 (bend given both).
  ZC <- matrix(NA_real_, n, 3, dimnames = list(NULL, c("z_c1", "z_c2", "z_c3")))
  gsd <- NULL
  for (k in names(byc)) {
    i <- byc[[k]]
    uc <- unit_se_components(B[i, , drop = FALSE], S0[[k]], d[i], u)
    ZC[i, ] <- uc$Z; gsd <- c(gsd, uc$max_abs_var_minus_1, uc$max_abs_offdiag)
  }
  # One output row per basis: het_test() plus percentile bootstrap intervals for r, tau and
  # tau / |r| when full = TRUE.
  het_row <- function(basis, component, L, S, vLL = 1, vSS = 1, vLS = 0, excess = "n_minus_1", boot = TRUE) {
    h <- het_test(L, S, vLL, vSS, vLS, excess = excess)
    ci <- c(NA, NA); tr_ci <- c(NA, NA); r_ci <- c(NA, NA)
    if (boot && full) {
      bt <- het_boot(L, S, vLL, vSS, vLS, B = nboot, seed = seed, excess = excess)
      ci <- bt$tau_ci; r_ci <- bt$r_ci
      tr_ci <- unname(quantile(bt$draws[, "tau"] / abs(bt$draws[, "r"]), c(0.025, 0.975), na.rm = TRUE))
    }
    list(fit = h, row = data.table(arm = label, basis = basis, component = component, excess = excess,
      n_leads = h$n, r = h$r, r_lo95 = r_ci[1], r_hi95 = r_ci[2], Q = h$Q, df = h$df, p = h$p,
      tau = h$tau, tau_lo95 = ci[1], tau_hi95 = ci[2], tau_over_r = h$tau_over_r,
      tau_over_r_lo95 = tr_ci[1], tau_over_r_hi95 = tr_ci[2], converged = h$converged,
      n_resid_z_gt_1.96 = sum(abs(h$resid_z) > 1.96), n_resid_z_gt_3 = sum(abs(h$resid_z) > 3)))
  }
  h_main <- het_row("unit_se", "c2_on_c1", ZC[, 1], ZC[, 2])
  HET <- rbindlist(list(
    h_main$row,
    het_row("unit_se", "c2_on_c1", ZC[, 1], ZC[, 2], excess = "n", boot = TRUE)$row,
    het_row("unit_se", "c3_on_c1", ZC[, 1], ZC[, 3])$row,
    het_row("exact", "delta_var_on_delta_mu", EST[, "delta_mu"], EST[, "delta_var"],
            VC$exact[, "vLL"], VC$exact[, "vSS"], VC$exact[, "vLS"])$row,
    het_row("lmoment", "lambda2_on_delta_mu", EST[, "delta_mu"], EST[, "delta_lambda2"],
            VC$lmoment[, "vLL"], VC$lmoment[, "vSS"], VC$lmoment[, "vLS"])$row,
    het_row("cornish_fisher", "cf_tilt_on_cf_level", EST[, "cf_level"], EST[, "cf_tilt"],
            VC$cornish_fisher[, "vLL"], VC$cornish_fisher[, "vSS"], VC$cornish_fisher[, "vLS"])$row))
  HET[, gs_max_abs_var_minus_1 := max(gsd[c(TRUE, FALSE)])]

  # Detection floor of the primary test (80% power at alpha 0.05).
  FLOOR <- NULL; PCURVE <- NULL
  if (full) {
    fl <- het_detection_floor(ZC[, 1], h_main$fit$r, nsim = nfloor, seed = seed)
    pw_at <- function(t) { i <- which.min(abs(fl$curve$tau - t)); if (abs(fl$curve$tau[i] - t) < 1e-9) fl$curve$power[i] else NA_real_ }
    FLOOR <- data.table(arm = label, basis = "unit_se", component = "c2_on_c1", n_leads = n,
      alpha = 0.05, target_power = 0.8, n_sim = nfloor, q_crit = fl$qcrit, r_used = h_main$fit$r,
      tau_floor_80 = fl$floor, tau_observed = h_main$fit$tau,
      observed_above_floor = !is.na(fl$floor) && h_main$fit$tau >= fl$floor,
      power_at_0.05 = pw_at(0.05), power_at_0.10 = pw_at(0.10), power_at_0.15 = pw_at(0.15),
      power_at_0 = pw_at(0))
    PCURVE <- data.table(arm = label, tau = fl$curve$tau, power = fl$curve$power)
  }

  ## 5. Principal components of the whitened profiles (one frame per trait: mean Sigma0)
  wf <- whiten_frame(S0bar)
  Sih <- wf$to_white; Sh <- wf$to_tau
  Y <- (B %*% t(Sih)) * sqrt(d)
  pc <- pca_unit(Y)
  pa <- parallel_analysis(Y, k_max = 4L, draws = npa, seed = seed)
  # Sign loading k to agree with the kth Hermite direction: shift 1, tilt z, bend z^2 - 1,
  # cubic z^3 - 3z (in the whitened metric).
  H4 <- cbind(shift = 1, tilt = z, bend = z^2 - 1, fourth = z^3 - 3 * z)
  Hw <- Sih %*% H4
  V <- pc$vectors[, 1:4]
  for (k in 1:4) if (sum(V[, k] * Hw[, k]) < 0) V[, k] <- -V[, k]
  Vtau <- Sh %*% V                                     # loadings on the quantile scale
  U <- Y / pc$norms
  SC <- U %*% V
  SPEC <- data.table(arm = label, pc = 1:10, eigenvalue = pc$values[1:10], pct = pc$pct[1:10],
                     cum_pct = cumsum(pc$pct)[1:10], isotropic_pct = 100 / nt)
  PAR <- data.table(arm = label, as.data.table(pa$table), n_components = pa$n_components)
  # Bases and single directions that the loadings are compared with.
  bases <- list(shift_spread_z = cbind(1, z), shift_spread_exact = cbind(1, q - mu_b),
                shift_spread_bend_z = cbind(1, z, z^2), lmoment_3 = cbind(1, 2 * u - 1, 6 * u^2 - 6 * u + 1))
  dirs <- list(shift = rep(1, nt), tilt_z = z, bend_z2m1 = z^2 - 1, exact_spread = q - mu_b,
               lmoment_spread = 2 * u - 1, lmoment_skew = 6 * u^2 - 6 * u + 1)
  BF <- list(); AN <- list()
  for (k in 1:4) {
    for (bn in names(bases)) {
      X <- bases[[bn]]; r2 <- basis_fit_r2(Vtau[, k], X)
      BF[[length(BF) + 1L]] <- data.table(arm = label, pc = k, pct = pc$pct[k], basis = bn, basis_rank = ncol(X),
        r2_centred = r2[["r2_centred"]], r2_uncentred = r2[["r2_uncentred"]],
        angle_to_span_tau_deg = angle_to_span(Vtau[, k], X),
        angle_to_span_whitened_deg = angle_to_span(V[, k], Sih %*% X))
    }
    for (dn in names(dirs)) {
      AN[[length(AN) + 1L]] <- data.table(arm = label, pc = k, direction = dn,
        angle_whitened_deg = angle_deg(V[, k], Sih %*% dirs[[dn]]),
        angle_tau_deg = angle_deg(Vtau[, k], dirs[[dn]]))
    }
  }
  LOAD <- rbindlist(lapply(1:4, function(k) data.table(arm = label, pc = k, tau = u, z = z,
    loading_tau = Vtau[, k] / sqrt(sum(Vtau[, k]^2)), loading_whitened = as.numeric(wf$E %*% V[, k]))))
  # How far each chromosome's Sigma0 is from the trait mean used for whitening.
  SIGI <- data.table(arm = label, chr = chrs, sigma0_rank = vapply(chrs, function(k) whiten_frame(S0[[k]])$rank, 0L),
    whitening_rank = wf$rank,
    max_rel_dev_from_mean = vapply(chrs, function(k) max(abs(S0[[k]] - S0bar)) / max(abs(S0bar)), 0))

  ## 6. One row per lead
  LT <- data.table(arm = label, meta[, c("lead_set", "snp_id", "chr", "bp", "effect_allele", "other_allele", "d_resid")])
  for (j in fn) {
    set(LT, j = j, value = EST[, j]); set(LT, j = paste0(j, "_se"), value = SE[, j])
    set(LT, j = paste0(j, "_z"), value = Zf[, j])
  }
  LT[, `:=`(lof_Q = LOF$shift_spread$Q, lof_p_davies = LOF$shift_spread$p_davies,
            lof_davies_floored = LOF$shift_spread$floored, lof_p_satterthwaite = LOF$shift_spread$p_satt,
            lof_cf3_Q = LOF$cf3$Q, lof_cf3_p_davies = LOF$cf3$p_davies)]
  LT[, `:=`(z_c1 = ZC[, 1], z_c2 = ZC[, 2], z_c3 = ZC[, 3], het_resid_z = h_main$fit$resid_z,
            whitened_norm = pc$norms, pc1_score = SC[, 1], pc2_score = SC[, 2], pc3_score = SC[, 3])]
  list(lead = LT, lof = rbindlist(LOFS), skew = SKW, het = HET, floor = FLOOR, pcurve = PCURVE,
       spec = SPEC, par = PAR, load = LOAD, bf = rbindlist(BF), ang = rbindlist(AN), sig = SIGI,
       n_components = pa$n_components, mu_b = mu_b, n_chr = length(chrs), n = n)
}
