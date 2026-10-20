#!/usr/bin/env Rscript
# Second principal component relative to noise, per trait, with a bootstrap interval
#
# What: For each of the 21 traits, the share of whitened shape variance carried by the second
#   principal component is divided by its median under the parallel-analysis null in which
#   every lead shares one shape (rank-1 signal plus unit noise; from
#   01_moments_and_principal_components.R). Because a change of scale rotates all whitened
#   profiles together, the ratio does not depend on the scale; 1 means no second shape. The
#   95% interval comes from 2,000 bootstrap resamples of leads with the null median held fixed.
#   Resampling with replacement duplicates leads, and duplicated rows carry duplicated noise,
#   which by itself inflates the PC2 share; so the percentile interval is re-centred on the
#   point estimate (shifted by the bootstrap median minus the estimate) and only its spread is
#   used. A sensitivity interval that redraws the null with the noise duplicated as the rows
#   are (50 null draws per resample) is reported beside it.
# Paper: Fig. 5b; Methods "Moment decomposition and principal components"
# Input: aligned lead-variant profiles and per-chromosome pooled covariance
#   (results/01_data/lead_profiles), baseline quantiles (results/01_data/baseline_quantiles),
#   results/06_traits/moments/{arms,pca_parallel}.csv, results/06_traits/survey/fig6_components.csv
# Output: results/06_traits/pc_heterogeneity/pc2_share_ratio.csv
# Runs: locally from results/ in about 2.5 hours on one core (traits are independent and can be
#   run in parallel)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
u <- tau_grid(cfg); n_tau <- length(u)
n_boot <- as.integer(cfg$estimation$het_boot)                                   # 2,000
seed <- as.integer(cfg$estimation$het_pc$split_seed)
n_null <- as.integer(cfg$estimation$pca_parallel_draws)                         # 300
n_null_per_resample <- as.integer(cfg$estimation$het_pc$pc2_share_null_draws_per_resample)   # 50

lead_profiles <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
baseline <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))
baseline_quantiles <- function(trait_name) {
  x <- baseline[baseline$trait == trait_name]
  as.numeric(x$q[match(round(u, 3), round(as.numeric(x$tau), 3))])
}
arms <- as.data.table(read_result("06_traits/moments", "arms"))
parallel <- as.data.table(read_result("06_traits/moments", "pca_parallel"))
fig5_traits <- as.data.table(read_result("06_traits/survey", "fig6_components"))$trait

# The 21 traits of Fig. 5: lipoprotein(a) is excluded (censored assay).
traits <- setdiff(fig5_traits, unlist(cfg$estimation$het_pc$exclude_traits))
stopifnot(length(traits) == 21L)

## 1. Whitened lead profiles on the analysed scale
# The profiles are taken through the raw scale and back (an exact identity for the analysed
# scale), in the same way as the scale comparison, then whitened with the trait's
# chromosome-mean Sigma0: y = sqrt(d) Sigma0^-1/2 beta, not scaled to unit length.
whitened_profiles <- function(trait_name) {
  leads <- lead_set_matrix(lead_profiles, arms[arm == trait_name]$lead_set, u)
  q_analysed <- baseline_quantiles(trait_name)
  is_log <- grepl("^log_", trait_name)
  q_raw <- if (is_log) exp(q_analysed) else q_analysed
  B_raw <- if (is_log) sweep(leads$B, 2, q_raw, `*`) else leads$B
  sigma0 <- p05_sigma0(leads$source, u)
  sigma0_raw <- if (is_log) lapply(sigma0, function(S) change_scale(matrix(0, 1, n_tau), S, q_raw)$Sigma0) else sigma0
  D <- q_raw^((if (is_log) 0 else 1) - 1)
  B <- sweep(B_raw, 2, D, `*`)
  sigma0 <- lapply(sigma0_raw, function(S) change_scale(matrix(0, 1, n_tau), S, D)$Sigma0)
  pc_cloud(B, sigma0, leads$meta$d_resid, as.character(leads$meta$chr))$Y
}

## 2. Ratio and bootstrap interval
pc2_share <- function(Y) pca_unit(Y)$pct[2]
# Median PC2 share under the one-shape null: the rank-1 reconstruction of Y plus unit isotropic
# noise, decomposed as the data are (parallel analysis with k = 1). With idx (a bootstrap
# resample) the noise is drawn once per original lead and indexed by idx, so duplicated leads
# carry duplicated noise, as their observed rows do.
null_median_pc2 <- function(Y, draws, idx = NULL) {
  n <- nrow(Y); r <- ncol(Y)                         # a resample has as many rows as the data
  sv <- svd(Y, nu = 1L, nv = 1L); Y1 <- sv$d[1] * sv$u %*% t(sv$v)
  100 * stats::median(replicate(draws, {
    E <- matrix(stats::rnorm(n * r), n)
    if (!is.null(idx)) E <- E[idx, , drop = FALSE]
    ev <- pca_unit(Y1 + E)$values
    ev[2] / sum(ev)
  }))
}
pc2_ratio <- function(trait_name) {
  Y <- whitened_profiles(trait_name); n <- nrow(Y)
  par2 <- parallel[arm == trait_name & component_tested == 2L]
  observed <- pc2_share(Y)
  null_median <- par2$null_median_pct
  set.seed(seed + 7L)
  boot <- t(vapply(seq_len(n_boot), function(b) {
    i <- sample.int(n, n, replace = TRUE)
    Yb <- Y[i, , drop = FALSE]
    share <- pc2_share(Yb)
    c(share = share, ratio_dup = share / null_median_pc2(Yb, n_null_per_resample, idx = i))
  }, numeric(2)))
  set.seed(seed + 9L)
  null_median_recomputed <- null_median_pc2(Y, n_null)          # check against the moments run
  q <- function(x, p) unname(stats::quantile(x, p))
  ratio_boot <- boot[, "share"] / null_median
  shift <- stats::median(ratio_boot) - observed / null_median   # the duplication artefact
  message(trait_name, ": ratio ", round(observed / null_median, 3))
  data.table(trait = trait_name, ratio = observed / null_median,
             lo95 = q(ratio_boot, 0.025) - shift, hi95 = q(ratio_boot, 0.975) - shift,
             null_median_pct = null_median, null_p95_pct = par2$null_p95_pct,
             ci_method = "fixed null; bootstrap percentile interval re-centred on the point estimate",
             boot_shift = shift, lo95_fixed_null_raw = q(ratio_boot, 0.025), hi95_fixed_null_raw = q(ratio_boot, 0.975),
             lo95_dup_null = q(boot[, "ratio_dup"], 0.025), hi95_dup_null = q(boot[, "ratio_dup"], 0.975),
             boot_dup_null_median = stats::median(boot[, "ratio_dup"]),
             observed_pc2_pct = observed, observed_pc2_pct_p11 = par2$observed_pct, n_leads = n, n_boot = n_boot,
             null_draws_per_resample = n_null_per_resample, null_median_recomputed_pct = null_median_recomputed)
}
ratios <- rbindlist(lapply(traits, pc2_ratio))

## 3. Write
write_result(ratios, "06_traits/pc_heterogeneity", "pc2_share_ratio")
