#!/usr/bin/env Rscript
# Decomposition of JCV and MCV antibody associations into membership and level
#
# What: log IgG antibody levels against the VP1 capsid protein of JC polyomavirus (JCV) and
#   Merkel cell polyomavirus (MCV) are fitted as a two-component skew-normal mixture, roughly
#   seronegative-like (lower) and seropositive-like (upper) participants, neither observed. For
#   each lead variant of the REGENIE GWAS this script splits the association into membership
#   (gamma = logit P(lower component); the paper reports -gamma, upper-component membership)
#   and the level of each component (beta_left, beta_right), and gives:
#   - a joint Hotelling test, the effect on the mean antibody level (the exact gradient of the
#     mixture mean, functional_grad()) and the share of it carried by membership, a test of
#     equal shifts of the two components, a descriptive label (level: |z_gamma| < 1.96 and
#     |z_left - right| < 1.96; membership: |z_gamma| >= 5.45; otherwise mixed), and the
#     mixture model's lack of fit (Davies);
#   - Stage 1 profiles with their fitted parts at the example loci (Fig. 3);
#   - genome-wide overlap of the joint and mean tests with REGENIE, Manhattan and QQ data;
#   - the fitted mixture densities and phenotype histograms (Fig. 3a);
#   - for each lead, how small a change of spread relative to location could be detected, and
#     the lack of fit of a shift + spread frame, with 2,000 random null variants as calibration.
#   Everything that depends on the Stage 2 weights is given under both weight designs.
# Paper: Fig. 3; Supplementary Fig. S9; Supplementary Table S4
# Input: per-variant Stage 2 results and weight matrices (02_mixtures); Stage 1 output of the
#   antibody scan; the antibody REGENIE GWAS (N = 9,145); lead variants
#   (results/01_data/lead_variants/leads.csv); UK Biobank antibody phenotype files for the
#   histograms
# Output: results/04_antibodies/decomposition/{lead_decomposition, example_profiles, overlap,
#   qq_lambda, qq_<design>, manhattan_points_<trait>_<design>, manhattan_bins_<trait>_<design>,
#   mixture_density, mixture_components, phenotype_histograms, detection_bounds_per_lead,
#   lack_of_fit_per_lead, lack_of_fit_summary}.csv
# Runs: on a cluster with UK Biobank data, in about 30 minutes

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
suppressPackageStartupMessages(library(gmoq))   # the paper used a development version of this package (Stage 2 commit 9fee67d)
cfg <- read_config()
p9 <- cfg$p09
TRAITS <- unlist(cfg$headline$antibody_traits)              # jcv_vp1_log_all, mcv_vp1_log_all
PRIMARY <- unique(mixture_design(cfg, TRAITS))              # model
DESIGNS <- c(PRIMARY, setdiff(WEIGHT_DESIGNS, PRIMARY))
EXAMPLES <- unlist(p9$example_loci)                         # FUT2 rs681343, HLA rs9271055, rs9276815
LEADS <- as.data.table(read_result("01_data/lead_variants", "leads"))
GWS <- -log10(5e-8)
out_dir <- "04_antibodies/decomposition"
u <- tau_grid(cfg); cw <- trap_weights(u)
mhc <- function(chr, bp) chr == 6 & bp >= 25e6 & bp <= 34e6
# Mean of the two-component skew-normal mixture, as a function of its parameters.
mixture_mean <- function(p) {
  d <- p$alpha / sqrt(1 + p$alpha^2); E <- p$xi + p$omega * d * sqrt(2 / pi)
  p$weight[1] * E[1] + p$weight[2] * E[2]
}

OUT <- list()
add <- function(name, x, design) {
  x <- as.data.table(x); x[, weight_design := design]; setcolorder(x, "weight_design")
  OUT[[name]] <<- c(OUT[[name]], list(x))
}
NULL_IDS <- list(); S_primary <- NULL

for (design in DESIGNS) {
  S <- list()
  for (tr in TRAITS) {
    n_an <- as.numeric(cfg$cohorts[[tr]]$n)

    ## 1. Genome-wide: Stage 2, mean-effect test and REGENIE, per chromosome --------------
    G <- list(); CH <- list()
    for (chr in 1:22) {
      g <- fread(stage2_file(cfg, design, tr, chr), showProgress = FALSE)
      wobj <- readRDS(weight_file(cfg, design, tr, chr)); W <- wobj$W
      Sigma0 <- gw_read_sigma0(file.path(wobj$stage1_dir, "scan.sigma0.tsv.gz"), wobj$taus)
      A <- solve(crossprod(W), t(W)); V_unit <- A %*% Sigma0 %*% t(A)
      # Effect on the mean antibody level: g' theta, with g the gradient of the mixture mean
      # with respect to (gamma, beta_left, beta_right); variance g' A Sigma0 A' g / d.
      grad <- gmoq::functional_grad(list(W = W, components = wobj$components), mixture_mean)[colnames(W)]
      theta <- as.matrix(g[, colnames(W), with = FALSE])
      g[, mean_eff := as.numeric(theta %*% grad)]
      g[, mean_se := sqrt(as.numeric(t(grad) %*% V_unit %*% grad) / d_resid)]
      g[, mlog10p_mean := mlog10_p_hotelling((mean_eff / mean_se)^2, 1, n_an)]
      rg <- fread(file.path(p9$regenie_dir, sprintf("chr%02d", chr), sprintf("antigen_qt_chr%02d_%s.regenie", chr, tr)),
                  select = c("ID", "ALLELE0", "ALLELE1", "N", "BETA", "SE", "LOG10P"), showProgress = FALSE)
      m <- match(g$snp_id, rg$ID)
      flip <- allele_sign(g$effect_allele, g$other_allele, rg$ALLELE1[m], rg$ALLELE0[m])
      g[, `:=`(regenie_beta = rg$BETA[m] * flip, regenie_se = rg$SE[m], mlog10p_regenie = rg$LOG10P[m],
               regenie_n = rg$N[m])]
      G[[chr]] <- g
      CH[[chr]] <- list(W = W, S0 = Sigma0, Vu = V_unit, grad = grad, components = wobj$components,
                        s1d = wobj$stage1_dir, rn = read_json(file.path(wobj$stage1_dir, "rif_nuisance.json")))
    }
    G <- rbindlist(G)
    q <- as.numeric(CH[[1]]$rn$q_tau)                    # baseline quantiles Q(u)

    ## 2. Leads, examples and a null sample -------------------------------------------------
    lead_ids <- LEADS[lead_set == paste0(tr, "_A_regenie"), snp_id]
    L <- G[snp_id %in% lead_ids]
    if (nrow(L) != length(unique(lead_ids))) stop(tr, ": lead variants missing from the G-MoQ scan")
    L[, `:=`(trait = tr, mhc = mhc(chr, bp))]
    # 2,000 variants with neither a REGENIE nor a joint signal (P > 0.1), drawn once under the
    # reported design, calibrate the lack-of-fit statistics.
    if (is.null(NULL_IDS[[tr]])) {
      set.seed(cfg$estimation$boot_seed)
      pool <- G[is.finite(mlog10p_regenie) & mlog10p_regenie < 1 & mlog10p_joint < 1]
      NULL_IDS[[tr]] <- pool$snp_id[sample.int(nrow(pool), min(p9$n_null, nrow(pool)))]
    }
    want <- unique(c(L$snp_id, EXAMPLES, NULL_IDS[[tr]]))
    B <- rbindlist(lapply(1:22, function(chr) {
      s1 <- read_stage1(CH[[chr]]$s1d)[snp_id %in% want]; s1[, chr := chr]; s1
    }))
    S[[tr]] <- list(G = G, L = L, B = B, bcols = grep("^beta_tau_", names(B), value = TRUE), CH = CH, q = q,
                    n = n_an, null_ids = NULL_IDS[[tr]])
  }

  ## 3. Decomposition of each lead and example --------------------------------------------
  dec <- list(); prof <- list()
  for (tr in TRAITS) {
    s <- S[[tr]]
    for (id in unique(c(s$L$snp_id, EXAMPLES))) {
      g <- s$G[snp_id == id]; b <- s$B[snp_id == id]
      if (!nrow(g) || !nrow(b)) next
      ch <- s$CH[[g$chr]]
      th <- as.numeric(g[, MIXTURE_PARAMS, with = FALSE]); V <- ch$Vu / g$d_resid
      wt <- wald_test(th, V, s$n)
      beta <- as.numeric(b[1, s$bcols, with = FALSE]); se_t <- sqrt(diag(ch$S0) / g$d_resid)
      fit <- as.numeric(ch$W %*% th)
      rms_z <- sqrt(mean(((beta - fit) / se_t)^2))
      lof <- lof_davies(beta, ch$S0, g$d_resid, ch$W)
      contrast <- c(0, 1, -1)                              # beta_left - beta_right
      z_diff <- sum(contrast * th) / sqrt(as.numeric(t(contrast) %*% V %*% contrast))
      z <- th / sqrt(diag(V)); z_mean <- g$mean_eff / g$mean_se
      label <- if (abs(z[1]) < 1.96 && abs(z_diff) < 1.96) "level" else if (abs(z[1]) >= 5.45) "membership" else "mixed"
      is_lead <- id %in% s$L$snp_id
      dec[[length(dec) + 1L]] <- data.table(trait = tr, snp_id = id, chr = g$chr, bp = g$bp,
        effect_allele = g$effect_allele, other_allele = g$other_allele, is_regenie_lead = is_lead, is_p04_lead = is_lead,
        is_example = id %in% EXAMPLES, mhc = mhc(g$chr, g$bp), d_resid = g$d_resid,
        gamma = th[1], gamma_se = sqrt(V[1, 1]), beta_left = th[2], beta_left_se = sqrt(V[2, 2]),
        beta_right = th[3], beta_right_se = sqrt(V[3, 3]), z_gamma = z[1], z_left = z[2], z_right = z[3],
        z_left_minus_right = z_diff, hotelling_T = wt$T, hotelling_k = 3L, hotelling_nu = s$n,
        hotelling_mlog10p = mlog10_p_hotelling(wt$T, 3, s$n),
        mean_effect = g$mean_eff, mean_se = g$mean_se, chi_mean = z_mean^2, frac_along_mean = z_mean^2 / wt$T,
        share_of_mean_from_gamma = ch$grad[["gamma"]] * th[1] / g$mean_eff,
        regenie_beta = g$regenie_beta, regenie_se = g$regenie_se, regenie_mlog10p = g$mlog10p_regenie,
        mean_over_regenie = g$mean_eff / g$regenie_beta, rms_z = rms_z,
        Q_mixture = lof$Q, p_mixture_davies = lof$p_davies, df_mixture = lof$df_nominal, label = label)
      if (id %in% EXAMPLES) {
        pr <- data.table(trait = tr, snp_id = id, tau = u, beta_stage1 = beta, se_stage1 = se_t, fitted_total = fit)
        for (j in seq_along(MIXTURE_PARAMS)) pr[[paste0("fitted_", MIXTURE_PARAMS[j])]] <- ch$W[, j] * th[j]
        prof[[length(prof) + 1L]] <- pr
      }
    }
  }
  add("lead_decomposition", rbindlist(dec), design)
  add("example_profiles", rbindlist(prof), design)

  ## 4. Genome-wide overlap with REGENIE, Manhattan and QQ data -----------------------------
  ov <- list(); qqs <- list(); lam <- list()
  for (tr in TRAITS) {
    G <- S[[tr]]$G[is.finite(mlog10p_regenie)]
    G[, mlog10p_any := pmax(mlog10p_joint, mlog10p_mean)]
    sets <- list(gmoq_joint = G$mlog10p_joint >= GWS, gmoq_mean = G$mlog10p_mean >= GWS,
                 gmoq_any = G$mlog10p_any >= GWS, regenie = G$mlog10p_regenie >= GWS)
    # Loci: 500 kb distance clumps; a locus is shared if a lead of the other set is within 500 kb.
    near <- function(x, y) if (!nrow(x) || !nrow(y)) rep(FALSE, nrow(x)) else
      vapply(seq_len(nrow(x)), function(i) any(y$chr == x$chr[i] & abs(y$bp - x$bp[i]) <= p9$lead_window_bp), logical(1))
    for (gs in c("gmoq_joint", "gmoq_mean", "gmoq_any")) {
      a <- sets[[gs]]; r <- sets$regenie
      ov[[length(ov) + 1L]] <- data.table(trait = tr, level = "snp", gmoq_test = gs, n_snps_tested = nrow(G),
        n_gmoq = sum(a), n_regenie = sum(r), shared = sum(a & r), gmoq_only = sum(a & !r), regenie_only = sum(!a & r))
      la <- as.data.table(distance_clump(G, paste0("mlog10p_", sub("gmoq_", "", gs)), 5e-8, p9$lead_window_bp))
      lr <- as.data.table(distance_clump(G, "mlog10p_regenie", 5e-8, p9$lead_window_bp))
      ov[[length(ov) + 1L]] <- data.table(trait = tr, level = "locus_500kb", gmoq_test = gs, n_snps_tested = nrow(G),
        n_gmoq = nrow(la), n_regenie = nrow(lr), shared = sum(near(la, lr)), gmoq_only = sum(!near(la, lr)),
        regenie_only = sum(!near(lr, la)))
    }
    mt <- manhattan_thin(G, c("mlog10p_joint", "mlog10p_mean", "mlog10p_regenie"),
                         keep_p = cfg$p07$manhattan_keep_p, bin_bp = cfg$p07$manhattan_bin_bp)
    write_result(as.data.table(mt$points)[, weight_design := design], out_dir, paste0("manhattan_points_", tr, "_", design))
    write_result(as.data.table(mt$bins)[, weight_design := design], out_dir, paste0("manhattan_bins_", tr, "_", design))
    for (t in c("joint", "mean", "regenie")) {
      x <- as.data.table(qq_quantiles(G[[paste0("mlog10p_", t)]])); x[, `:=`(trait = tr, test = t)]
      qqs[[length(qqs) + 1L]] <- x
    }
    lam[[tr]] <- data.table(trait = tr, lambda_joint = median(G$joint_T) / qchisq(0.5, 3),
      lambda_mean = median((G$mean_eff / G$mean_se)^2) / qchisq(0.5, 1),
      lambda_regenie = median((G$regenie_beta / G$regenie_se)^2, na.rm = TRUE) / qchisq(0.5, 1), n_snps = nrow(G))
  }
  add("overlap", rbindlist(ov), design)
  write_result(rbindlist(qqs)[, weight_design := design], out_dir, paste0("qq_", design))
  add("qq_lambda", rbindlist(lam), design)
  if (design == PRIMARY) S_primary <- S
}
for (nm in names(OUT)) write_result(rbindlist(OUT[[nm]], fill = TRUE), out_dir, nm)
S <- S_primary

## 5. Fitted mixture densities and phenotype histograms (Fig. 3a) -------------------------
dens <- list(); comps <- list(); hist_rows <- list()
for (tr in TRAITS) {
  s <- S[[tr]]; cmp <- as.data.frame(s$CH[[1]]$components)
  comps[[tr]] <- cbind(data.table(trait = tr), as.data.table(cmp))
  x <- seq(min(s$q) - 1.5, max(s$q) + 1.5, length.out = 400)
  lo <- cmp[cmp$component == "lower", ]; up <- cmp[cmp$component == "upper", ]
  dl <- lo$weight * dsn_density(x, lo$xi, lo$omega, lo$alpha); du <- up$weight * dsn_density(x, up$xi, up$omega, up$alpha)
  dens[[tr]] <- data.table(trait = tr, x = x, density_total = dl + du, density_lower = dl, density_upper = du)
  # Histogram of the analysed sample; inner edges between the 0.5th and 99.5th percentiles.
  inp <- stage1_inputs(cfg, tr)
  y <- stage1_analysed_sample(inp$pheno, inp$trait, inp$covar, inp$keep, inp$loco, inp$bgen_sample)$y
  if (length(y) != s$n) stop(tr, ": analysed sample has ", length(y), " participants, not ", s$n)
  edges <- as.numeric(stats::quantile(y, c(0.005, 0.995)))
  h <- privacy_histogram(y, edges[1], edges[2], p9$histogram_bins, cfg$privacy$min_cell_count)
  hist_rows[[tr]] <- data.table(trait = tr, h[, .(bin_lower, bin_upper, count)], n = length(y), density = h$density)
}
write_result(rbindlist(dens), out_dir, "mixture_density")
write_result(rbindlist(comps, fill = TRUE), out_dir, "mixture_components")
write_result(rbindlist(hist_rows), out_dir, "phenotype_histograms")

## 6. Shift and spread of each lead: detection bounds and lack of fit ---------------------
# Location L = change in the mean, int beta(u) du; spread S = change in the variance divided by
# twice the baseline s.d. Leads are aligned to L > 0. With the common ratio r of S to L across
# leads (het_test), a lead's spread is detectably different from r L when
# |S - r L| > z sqrt(v0), v0 = Var(S - r L); the bound is z sqrt(v0) / |L| with z = 1.96 or the
# Bonferroni value over 2n tests. Lack of fit of the shift + spread frame X = [1, Q(u) - mu]
# uses the OLS residual and the exact Davies reference; location-only and
# location + spread + third-L-moment frames, the Imhof and Satterthwaite p-values, the earlier
# whitened chi-square statistic and Cornish-Fisher level and tilt are given for context.
det <- list(); lof <- list(); lsum <- list()
for (tr in TRAITS) {
  s <- S[[tr]]; q <- s$q
  fw <- functional_weights(u, q)
  mu_b <- sum(cw * q) / sum(cw); sd_b <- sqrt(sum(cw * (q - mu_b)^2) / sum(cw))
  wL <- fw$delta_mu; wS <- fw$delta_var / (2 * sd_b)
  qs <- (q - mu_b) / sd_b; z <- qnorm(u)
  X3 <- cbind(1, z, z^2 - 1); G3 <- solve(crossprod(X3 * sqrt(cw)), t(X3 * cw))
  kC <- as.numeric(G3 %*% qs)[2]; wCFl <- G3[1, ]; wCFs <- G3[2, ] / kC
  frames <- list(location_only = matrix(1, length(u), 1), location_scale = cbind(1, q - mu_b),
                 location_scale_third = cbind(1, q - mu_b, 6 * u^2 - 6 * u + 1))
  per_variant <- function(ids, role) {
    rows <- list()
    for (id in ids) {
      b <- s$B[snp_id == id]; g <- s$G[snp_id == id]; if (!nrow(b) || !nrow(g)) next
      ch <- s$CH[[b$chr[1]]]; beta <- as.numeric(b[1, s$bcols, with = FALSE]); d <- b$d_resid[1]
      r <- data.table(trait = tr, snp_id = id, chr = b$chr[1], bp = b$bp[1], role = role,
                      mhc = mhc(b$chr[1], b$bp[1]), d_resid = d)
      Lx <- sum(wL * beta); Sx <- sum(wS * beta); fl <- if (Lx < 0) -1 else 1
      r[, `:=`(L = fl * Lx, S = fl * Sx, vLL = as.numeric(t(wL) %*% ch$S0 %*% wL) / d,
               vSS = as.numeric(t(wS) %*% ch$S0 %*% wS) / d, vLS = as.numeric(t(wL) %*% ch$S0 %*% wS) / d,
               L_cf = fl * sum(wCFl * beta), S_cf = fl * sum(wCFs * beta))]
      Sw <- sym_pow(ch$S0, -0.5); yw <- sqrt(d) * as.numeric(Sw %*% beta)
      for (nm in names(frames)) {
        X <- frames[[nm]]
        o <- lof_davies(beta, ch$S0, d, X, imhof_check = (nm == "location_scale"))
        r[[paste0("Q_ols_", nm)]] <- o$Q; r[[paste0("p_davies_", nm)]] <- o$p_davies
        r[[paste0("davies_floored_", nm)]] <- o$davies_floored
        r[[paste0("df_", nm)]] <- o$df_nominal
        if (nm == "location_scale") { r$p_imhof_location_scale <- o$p_imhof; r$p_satterthwaite_location_scale <- o$p_satterthwaite }
        Xw <- Sw %*% X; H <- Xw %*% solve(crossprod(Xw), t(Xw))
        r[[paste0("chisq_gls_whitened_", nm)]] <- sum((yw - H %*% yw)^2)
      }
      rows[[length(rows) + 1L]] <- r
    }
    rbindlist(rows)
  }
  LD <- per_variant(s$L$snp_id, "lead"); NU <- per_variant(s$null_ids, "null")
  h <- het_test(LD$L, LD$S, LD$vLL, LD$vSS, LD$vLS, excess = "n_minus_1")
  LD[, v0 := LD$vSS - 2 * h$r * LD$vLS + h$r^2 * LD$vLL]
  z_bonf <- qnorm(1 - 0.05 / (2 * 2 * nrow(LD)))
  LD[, `:=`(bound_nominal = 1.96 * sqrt(v0) / abs(L), bound_bonferroni = z_bonf * sqrt(v0) / abs(L),
            resid_z_from_line = h$resid_z, z_location = L / sqrt(vLL), z_spread = S / sqrt(vSS))]
  det[[tr]] <- LD[, .(trait, snp_id, chr, bp, mhc, L, S, vLL, vSS, vLS, z_location, z_spread, v0,
                      resid_z_from_line, bound_nominal, bound_bonferroni)]
  lof[[tr]] <- rbind(LD, NU, fill = TRUE)
  for (rl in c("lead", "null")) {
    dd <- lof[[tr]][role == rl]
    lsum[[length(lsum) + 1L]] <- data.table(trait = tr, role = rl, n = nrow(dd), df_location_scale = 43,
      median_Q_ols_location_scale = median(dd$Q_ols_location_scale),
      median_Q_ols_location_only = median(dd$Q_ols_location_only),
      median_Q_ols_location_scale_third = median(dd$Q_ols_location_scale_third),
      median_chisq_gls_location_scale = median(dd$chisq_gls_whitened_location_scale),
      frac_p_davies_lt_0.05 = mean(dd$p_davies_location_scale < 0.05),
      frac_p_davies_lt_bonferroni = mean(dd$p_davies_location_scale < 0.05 / max(1, nrow(dd))),
      median_p_davies = median(dd$p_davies_location_scale),
      max_abs_log10_davies_vs_imhof = max(abs(log10(dd$p_davies_location_scale) -
                                               log10(pmax(dd$p_imhof_location_scale, 1e-300))), na.rm = TRUE))
  }
}
write_result(rbindlist(det), out_dir, "detection_bounds_per_lead")
LOF <- rbindlist(lof, fill = TRUE)
write_result(LOF[role == "lead"], out_dir, "lack_of_fit_per_lead")
write_result(rbindlist(lsum), out_dir, "lack_of_fit_summary")
