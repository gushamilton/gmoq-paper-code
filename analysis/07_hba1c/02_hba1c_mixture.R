#!/usr/bin/env Rscript
# A two-component mixture for log HbA1c, fitted blind to diabetes, and genetic effects on it
#
# What: log HbA1c has three supported principal components of genetic effect. As a candidate
#   explanation, a two-component skew-normal mixture is fitted to log HbA1c in the Stage 1 sample
#   without using diabetes status (mixsmsn, as for testosterone; best of three seeds). Its upper
#   component is then compared with recorded diabetes (prevalence, posterior membership by
#   status). With G-MoQ Stage 2, each of the 472 lead profiles is mapped to effects on three
#   mixture parameters: membership (gamma = logit of the lower component's weight; reported as
#   gamma_up = -gamma, positive = more mass in the upper component) and the locations of the
#   lower and upper components. Weights are built by gmoq with every column divided by
#   the mixture density ("model") and, as a sensitivity, with the membership column divided by
#   the Stage 1 kernel density ("hybrid"). theta = (W'W)^-1 W' beta with
#   Var(theta) = A Sigma0 A' / d. Membership z is then set against the diabetes z of
#   01_skew_and_diabetes.R (Spearman), and lower-component location against the no-diabetes
#   HbA1c effect.
# Paper: Fig. 5c (membership weight column, dashed fits); Results (HbA1c mixture: upper weight,
#   prevalence, Spearman correlations)
# Input: log HbA1c and diabetes status of the Stage 1 sample (private file written by
#   01_skew_and_diabetes.R); aggregate results/01_data/{lead_profiles,baseline_quantiles} and
#   results/07_hba1c/skew_diabetes/{per_lead,sample_counts}.csv
# Output: results/07_hba1c/mixture/{mixture_params,mixture_seeds,classification,weight_matrices,
#   per_lead,stage2_summary}.csv
# Runs: on a cluster with UK Biobank data (the three mixture fits take about an hour each; run in
#   parallel)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
suppressPackageStartupMessages(library(gmoq))

cfg <- read_config()
pheno_file <- file.path(cfg$bp$private_root, "07_hba1c", "skew_diabetes", "pheno", "pheno.txt")
seeds <- c(20260915L, 1L, 2L)
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
u <- tau_grid(cfg)
z_bonf <- qnorm(1 - 0.025 / 472)

## 1. Sample: log HbA1c and diabetes status (participant level, in memory only)
ph <- fread(pheno_file, select = c("diab", "log_hba1c_all"))
y <- ph$log_hba1c_all; diabetes <- ph$diab == 2L
rm(ph)
counts <- as.data.table(read_result("07_hba1c/skew_diabetes", "sample_counts"))
stopifnot(length(y) == cfg$cohorts$survey$traits$log_hba1c$n, !anyNA(y),
          sum(diabetes) == counts[sample == "log_hba1c" & stage == "diabetes_cases"]$n_participants)
baseline <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))[trait == "log_hba1c"]
baseline <- baseline[match(round(u, 3), round(as.numeric(baseline$tau), 3))]
q_tau <- as.numeric(baseline$q); f_kde <- as.numeric(baseline$fhat)   # Stage 1 quantiles and kernel density
# The sample reproduces the Stage 1 quantiles (type 8, as Stage 1 computes them).
stopifnot(max(abs(quantile(y, u, type = 8, names = FALSE) - q_tau)) <= 1e-5)

## 2. Two-component skew-normal mixture, blind to diabetes
# Skew-normal moments of each component, for reporting.
sn_moments <- function(xi, om, al) {
  de <- al / sqrt(1 + al^2)
  list(mean = xi + om * de * sqrt(2 / pi), sd = sqrt(om^2 * (1 - 2 * de^2 / pi)),
       skew = (4 - pi) / 2 * (de * sqrt(2 / pi))^3 / (1 - 2 * de^2 / pi)^1.5)
}
fit_mixture <- function(seed) {
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(mixsmsn::smsn.mix(y = y, g = 2, family = "Skew.normal", nu = 1, get.init = TRUE,
                                  iter.max = 500, error = 1e-4, calc.im = FALSE, criteria = FALSE),
                error = function(e) e)
  if (inherits(f, "error")) return(list(ok = FALSE, msg = conditionMessage(f), seed = seed))
  comp <- data.table(weight = as.numeric(f$pii), xi = as.numeric(f$mu), omega = sqrt(as.numeric(f$sigma2)),
                     alpha = as.numeric(f$shape))
  mo <- sn_moments(comp$xi, comp$omega, comp$alpha)
  comp[, `:=`(mean = mo$mean, sd = mo$sd, skew = mo$skew)]
  comp <- comp[order(mean)]                               # lower component first
  dens <- sapply(1:2, function(k) comp$weight[k] * sn::dsn(y, xi = comp$xi[k], omega = comp$omega[k], alpha = comp$alpha[k]))
  list(ok = TRUE, comp = comp, loglik = sum(log(rowSums(dens))), iter = if (!is.null(f$iter)) f$iter else NA_integer_,
       seed = seed, seconds = proc.time()[["elapsed"]] - t0)
}
runs <- parallel::mclapply(seeds, fit_mixture, mc.cores = max(1L, min(n_cores, length(seeds))), mc.preschedule = FALSE)
mixture_seeds <- data.table(fit = "sn2_all", g = 2L, family = "Skew.normal", sample = "all", seed = seeds,
  ok = vapply(runs, function(r) isTRUE(r$ok), logical(1)),
  loglik = vapply(runs, function(r) if (isTRUE(r$ok)) r$loglik else NA_real_, numeric(1)),
  iter = vapply(runs, function(r) if (isTRUE(r$ok)) as.numeric(r$iter) else NA_real_, numeric(1)),
  seconds = vapply(runs, function(r) if (isTRUE(r$ok)) r$seconds else NA_real_, numeric(1)),
  msg = vapply(runs, function(r) if (isTRUE(r$ok)) "" else r$msg, character(1)))
best <- runs[[which(mixture_seeds$ok)[which.max(mixture_seeds$loglik[mixture_seeds$ok])]]]
comp <- best$comp
n_par <- 4L * 2L - 1L
mixture_params <- copy(comp)[, `:=`(fit = "sn2_all", family = "Skew.normal", sample = "all", n = length(y),
  component = c("lower", "upper"), mean_mmol_mol = exp(mean), loglik = best$loglik, n_par = n_par,
  bic = -2 * best$loglik + n_par * log(length(y)), seed = best$seed)]
setcolorder(mixture_params, c("fit", "family", "sample", "n", "component"))
print(mixture_params)

## 3. Upper component against recorded diabetes
component_density <- function(x, k) comp$weight[k] * sn::dsn(x, xi = comp$xi[k], omega = comp$omega[k], alpha = comp$alpha[k])
posterior_upper <- function(x) { d1 <- component_density(x, 1); d2 <- component_density(x, 2); d2 / (d1 + d2) }
auc <- function(s, cl) { r <- frank(s); n1 <- as.numeric(sum(cl)); n0 <- as.numeric(sum(!cl)); (sum(r[cl]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
overlap <- function() {          # integral of min(f_lower, f_upper), component densities unweighted
  x <- seq(min(y), max(y), length.out = 20001)
  sum(pmin(sn::dsn(x, comp$xi[1], comp$omega[1], comp$alpha[1]), sn::dsn(x, comp$xi[2], comp$omega[2], comp$alpha[2]))) * diff(x[1:2])
}
# Where the model posterior crosses 0.5, on a fixed grid (a property of the model, not the data).
crossings <- function() {
  xg <- seq(log(10), log(250), length.out = 20001); p <- posterior_upper(xg) - 0.5
  i <- which(diff(sign(p)) != 0); (xg[i] + xg[i + 1]) / 2
}
p <- posterior_upper(y); upper <- p >= 0.5
classification <- list()
add <- function(q, v, note = "") classification[[length(classification) + 1L]] <<- data.table(fit = "sn2_all", quantity = q, value = as.numeric(v), note = note)
add("upper_weight", comp$weight[2]); add("overlap_coef_lower_upper", overlap(), "integral min(f_lower, f_upper), unweighted")
add("diabetes_prevalence", mean(diabetes))
add("auc_posterior_upper_vs_diabetes", auc(p, diabetes))
add("auc_log_hba1c_vs_diabetes", auc(y, diabetes), "reference: raw log HbA1c")
add("sensitivity_at_0.5", mean(upper[diabetes])); add("specificity_at_0.5", mean(!upper[!diabetes]))
add("frac_assigned_upper", mean(upper)); add("n_assigned_upper", sum(upper))
add("ppv_frac_cases_among_assigned_upper", mean(diabetes[upper]))
add("frac_cases_in_upper_by_weight", sum(diabetes & upper) / sum(diabetes))
add("mean_posterior_cases", mean(p[diabetes])); add("mean_posterior_controls", mean(p[!diabetes]))
add("sum_posterior_cases_over_sum_posterior", sum(p[diabetes]) / sum(p), "expected fraction of upper-component mass that is diabetes")
xc <- crossings()
for (j in seq_along(xc)) add(paste0("posterior_0.5_crossing_", j, "_hba1c_mmol_mol"), exp(xc[j]), "model-based (grid), not data")
add("posterior_monotone_in_y", as.numeric(!is.unsorted(p[order(y)])))
classification <- rbindlist(classification)
rm(p, upper); invisible(gc())

## 4. Stage 2 weights for the mixture parameters
# Callbacks for gmoq: mixture CDF and density, and the derivatives of the CDF with
# respect to gamma = logit(p_lower) and the two component locations (finite differences).
params <- list(p1 = comp$weight[1], xi1 = comp$xi[1], omega1 = comp$omega[1], alpha1 = comp$alpha[1],
               xi2 = comp$xi[2], omega2 = comp$omega[2], alpha2 = comp$alpha[2])
mix_cdf <- function(x, z) z$p1 * sn::psn(x, xi = z$xi1, omega = z$omega1, alpha = z$alpha1) +
  (1 - z$p1) * sn::psn(x, xi = z$xi2, omega = z$omega2, alpha = z$alpha2)
mix_pdf <- function(x, z) z$p1 * sn::dsn(x, xi = z$xi1, omega = z$omega1, alpha = z$alpha1) +
  (1 - z$p1) * sn::dsn(x, xi = z$xi2, omega = z$omega2, alpha = z$alpha2)
grad_p1 <- make_fd_grad("p1", mix_cdf)
grads <- list(gamma = function(x, z) grad_p1(x, z) * z$p1 * (1 - z$p1),   # chain rule for logit
              beta_left = make_fd_grad("xi1", mix_cdf), beta_right = make_fd_grad("xi2", mix_cdf))
tol <- cfg$estimation$weight_gates$beta_left_plus_right_tol

## 5. Effects of the 472 lead variants on the mixture parameters
lead_profiles <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
leads <- lead_set_matrix(lead_profiles, "log_hba1c_A_deltamu", u)   # HbA1c-raising allele
sigma0 <- p05_sigma0(leads$source, u)
M <- as.data.table(leads$meta); B <- leads$B
stopifnot(nrow(B) == 472L)
skew_leads <- as.data.table(read_result("07_hba1c/skew_diabetes", "per_lead"))
H <- skew_leads[match(M$snp_id, snp_id)]
stopifnot(all(H$effect_allele == M$effect_allele))
per_lead <- data.table(snp_id = M$snp_id, chr = M$chr, bp = M$bp, effect_allele = M$effect_allele,
                       other_allele = M$other_allele, d_resid = M$d_resid, skew_hit = H$skew_hit,
                       skew_beyond_z = H$skew_beyond_z, delta_mu_z = H$delta_mu_z, diab_z = H$diab_z,
                       nd_delta_mu_z = H$nd_delta_mu_z, lin_nd_log_z = H$lin_nd_log_z, nearest_gene = H$nearest_gene)
weight_matrices <- list(); summary_rows <- list()
for (design in c("model", "hybrid")) {
  W <- make_weights_generic(taus = u, q_tau = q_tau, dist_cdf = mix_cdf, dist_pdf = mix_pdf,
    params = params, grad_funcs = grads, density = design,
    stage1_density = if (design == "model") NULL else
      list(taus = u, q_tau = q_tau, fhat = f_kde, source = "Stage 1 kernel density at the Stage 1 quantiles"),
    stage1_columns = "gamma", unit_sum = list(c("beta_left", "beta_right")), identity_tol = tol)
  stopifnot(identical(weight_density_mode(W), design))
  W <- unclass(W); attributes(W) <- attributes(W)[c("dim", "dimnames")]
  # A shift of both locations together is a pure shift: the two location columns sum to 1.
  stopifnot(max(abs(rowSums(W[, c("beta_left", "beta_right")]) - 1)) <= tol, qr(W)$rank == 3L)
  weight_matrices[[design]] <- data.table(mixture = "sn2_all", weight_design = design, tau = round(u, 3), W)
  A <- solve(crossprod(W), t(W))                            # theta = A beta (OLS)
  theta <- B %*% t(A); se <- matrix(NA_real_, nrow(B), 3)
  for (k in unique(as.character(M$chr))) {
    i <- which(as.character(M$chr) == k)
    se[i, ] <- sqrt(outer(1 / M$d_resid[i], diag(A %*% sigma0[[k]] %*% t(A))))
  }
  z <- theta / se
  pre <- paste0("sn2_all_", design, "_")
  per_lead[[paste0(pre, "gamma_up")]] <- -theta[, 1]; per_lead[[paste0(pre, "gamma_up_se")]] <- se[, 1]
  per_lead[[paste0(pre, "gamma_up_z")]] <- -z[, 1]
  per_lead[[paste0(pre, "beta_left")]] <- theta[, 2]; per_lead[[paste0(pre, "beta_left_se")]] <- se[, 2]
  per_lead[[paste0(pre, "beta_left_z")]] <- z[, 2]
  per_lead[[paste0(pre, "beta_right")]] <- theta[, 3]; per_lead[[paste0(pre, "beta_right_se")]] <- se[, 3]
  per_lead[[paste0(pre, "beta_right_z")]] <- z[, 3]

  ## 6. Membership against diabetes; lower location against HbA1c without diabetes
  gz <- -z[, 1]; lz <- z[, 2]; rz <- z[, 3]; h <- per_lead$skew_hit
  sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
  lf <- lm(gz ~ per_lead$delta_mu_z + per_lead$diab_z)
  r <- list(
    pearson_gamma_up_z_vs_diab_z = cor(gz, per_lead$diab_z, use = "complete.obs"),
    spearman_gamma_up_z_vs_diab_z = sp(gz, per_lead$diab_z),
    spearman_gamma_up_z_vs_diab_z_hits = sp(gz[h], per_lead$diab_z[h]),
    spearman_gamma_up_z_vs_diab_z_others = sp(gz[!h], per_lead$diab_z[!h]),
    lm_gamma_up_z_on_diab_z_adj_mean_z_coef = coef(lf)[3],
    lm_gamma_up_z_on_diab_z_adj_mean_z_p = summary(lf)$coefficients[3, 4],
    spearman_gamma_up_z_vs_skew_beyond_z = sp(gz, per_lead$skew_beyond_z),
    spearman_gamma_up_z_vs_delta_mu_z = sp(gz, per_lead$delta_mu_z),
    median_gamma_up_z_hits = median(gz[h]), median_gamma_up_z_others = median(gz[!h]),
    n_gamma_up_bonf472_pos_hits = sum(gz[h] > z_bonf), n_gamma_up_bonf472_pos_others = sum(gz[!h] > z_bonf),
    n_gamma_up_bonf472_neg_all = sum(gz < -z_bonf), n_gamma_bonf472_all = sum(abs(gz) > z_bonf),
    pearson_beta_left_z_vs_nd_delta_mu_z = cor(lz, per_lead$nd_delta_mu_z, use = "complete.obs"),
    spearman_beta_left_z_vs_nd_delta_mu_z = sp(lz, per_lead$nd_delta_mu_z),
    pearson_beta_left_z_vs_lin_nd_log_z = cor(lz, per_lead$lin_nd_log_z, use = "complete.obs"),
    pearson_beta_right_z_vs_diab_z = cor(rz, per_lead$diab_z, use = "complete.obs"),
    pearson_beta_right_z_vs_lin_nd_log_z = cor(rz, per_lead$lin_nd_log_z, use = "complete.obs"),
    median_abs_beta_left_z = median(abs(lz)), median_abs_beta_right_z = median(abs(rz)),
    median_abs_gamma_z = median(abs(gz)))
  summary_rows[[design]] <- data.table(mixture = "sn2_all", weight_design = design, quantity = names(r),
                                       value = as.numeric(unlist(r)))
}
stage2_summary <- rbindlist(summary_rows)
print(dcast(stage2_summary, quantity ~ weight_design, value.var = "value"))

## 7. Write
out <- "07_hba1c/mixture"
write_result(mixture_params, out, "mixture_params")
write_result(mixture_seeds, out, "mixture_seeds")
write_result(classification, out, "classification")
write_result(rbindlist(weight_matrices), out, "weight_matrices")
write_result(per_lead, out, "per_lead")
write_result(stage2_summary, out, "stage2_summary")
