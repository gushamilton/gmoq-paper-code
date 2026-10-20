#!/usr/bin/env Rscript
# Moments, lack of fit, between-lead heterogeneity and principal components of lead profiles
#
# What: For every trait, each lead variant's Stage 1 quantile-effect profile b(u) is summarised
#   without a parametric model. (1) Linear functionals w'b with exact pooled-covariance SEs:
#   location (delta_mu), exact spread (delta_var), L-moment spread and skew, and "skew beyond
#   shift + spread" (the third L-moment of the residual after fitting the frame [1, Q(u) - mu]).
#   (2) Lack of fit of that frame, with exact Davies p-values. (3) Between-lead heterogeneity of
#   spread relative to location, tau_h, with a 2,000-resample bootstrap CI and a simulated
#   detection floor. (4) Uncentred PCA of the whitened, unit-length profiles, with the number of
#   components chosen by sequential parallel analysis (300 null draws). All four are computed
#   by moment_decomposition() in R/.
# Paper: Fig. 5a,b; Supplementary Fig. S10; Supplementary Tables S5 and S6; Methods "Moment
#   decomposition and principal components"
# Input: aligned lead-variant profiles and per-chromosome pooled covariance
#   (results/01_data/lead_profiles), baseline quantiles Q(u) on the Stage 1 sample
#   (results/01_data/baseline_quantiles)
# Output: results/06_traits/moments/{arms,lead_moments,lof_summary,skew_counts,heterogeneity,
#   detection_floor,pca_spectrum,pca_parallel,pca_loadings}.csv
# Runs: locally from results/ in about 2 minutes on one core

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
u <- tau_grid(cfg)                                   # 45 levels, 0.02 ... 0.98
n_tau <- length(u)
n_boot <- as.integer(cfg$estimation$het_boot)        # 2,000 bootstrap resamples over leads
n_null_pca <- as.integer(cfg$estimation$pca_parallel_draws)   # 300 parallel-analysis draws
seed <- as.integer(cfg$estimation$boot_seed)
n_floor_sims <- 1000L                                # simulations per tau for the detection floor

lead_profiles <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
baseline <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))
baseline_n <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_summary"))

# Baseline quantiles Q(u) of a trait on the analysed scale, on the 45-level grid.
baseline_quantiles <- function(trait_name) {
  x <- baseline[baseline$trait == trait_name]
  as.numeric(x$q[match(round(u, 3), round(as.numeric(x$tau), 3))])
}

## 1. The analyses ("arms")
# Raw BMI, log BMI (the headline scale: the same 639 leads, BMI-increasing allele) and
# height-adjusted weight use REGENIE-clumped leads. The survey traits use leads clumped on the
# Stage 1 location effect; traits marked fig6: false in config.yaml are not analysed (mass
# points and lack of fit), nor are traits without a completed Stage 1 (n: null). The extra arm
# hba1c_nodiab_log is HbA1c without diabetes, converted exactly to the log scale so that it can
# be set beside log HbA1c (with diabetes).
survey <- cfg$cohorts$survey$traits
survey_traits <- names(survey)[vapply(survey, function(x) !identical(x$fig6, FALSE) &&
                                        !is.null(x[["n", exact = TRUE]]), logical(1))]
arms <- rbind(
  data.table(arm = c("bmi", "log_bmi", "weight"), trait = c("bmi", "log_bmi", "weight"),
             lead_set = c("bmi_A_regenie", "log_bmi_A_regenie", "weight_A_regenie"),
             group = c("bmi", "bmi", "weight"), transform = "none", scale = c("raw", "log", "raw")),
  data.table(arm = survey_traits, trait = survey_traits, lead_set = paste0(survey_traits, "_A_deltamu"),
             group = "survey", transform = "none",
             scale = ifelse(grepl("^log_", survey_traits), "log", "raw")),
  data.table(arm = "hba1c_nodiab_log", trait = "hba1c_nodiab", lead_set = "hba1c_nodiab_A_deltamu",
             group = "survey_extra", transform = "log", scale = "log (converted from raw)"))
arms[, n_sample := baseline_n$n[match(trait, baseline_n$trait)]]
arms[, source := "lead_profiles"]

## 2. The moment reading of each arm
analyse_arm <- function(a) {
  leads <- lead_set_matrix(lead_profiles, a$lead_set, u)
  B <- leads$B
  sigma0 <- p05_sigma0(leads$source, u)              # per chromosome, plus "mean" (the PCA frame)
  q <- baseline_quantiles(a$trait)
  if (a$transform == "log") {
    # Exact change of scale to log: beta_log(u) = beta(u) / Q(u), Sigma_log = D Sigma0 D with
    # D = diag(1 / Q(u)); the quantiles map as log Q(u).
    D <- 1 / q
    B <- change_scale(B, sigma0[[1]], D)$B
    sigma0 <- lapply(sigma0, function(S) change_scale(matrix(0, 1, n_tau), S, D)$Sigma0)
    q <- log(q)
  }
  md <- moment_decomposition(a$arm, leads$meta, B, sigma0, q, u, full = TRUE, nboot = n_boot,
                             npa = n_null_pca, seed = seed, nfloor = n_floor_sims)
  md$arm <- data.table(as.data.table(a), n_requested = nrow(leads$meta), n_leads = md$n,
                       n_chr = md$n_chr, mu_b = md$mu_b, n_components = md$n_components)
  message(a$arm, ": ", md$n, " leads, ", md$n_components, " components")
  md
}
# Each arm sets its own seeds inside the bootstrap, floor and parallel-analysis functions, so
# arms can be run in any order (or in parallel) with identical results.
results <- lapply(seq_len(nrow(arms)), function(i) analyse_arm(as.list(arms[i])))
collect <- function(part) rbindlist(lapply(results, `[[`, part), fill = TRUE)

## 3. Write
out <- "06_traits/moments"
write_result(collect("arm"), out, "arms")
write_result(collect("lead"), out, "lead_moments")          # one row per lead: functionals, LoF, PC scores
write_result(collect("lof"), out, "lof_summary")            # lack of fit of shift + spread (and Cornish-Fisher) frame
write_result(collect("skew"), out, "skew_counts")           # spread / skew / skew-beyond lead counts, Bonferroni within trait
write_result(collect("het"), out, "heterogeneity")          # tau_h by basis; unit_se, c2_on_c1, excess n_minus_1 is reported
write_result(collect("floor"), out, "detection_floor")      # tau_h at 80% power, same statistic
write_result(collect("spec"), out, "pca_spectrum")
write_result(collect("par"), out, "pca_parallel")
write_result(collect("load"), out, "pca_loadings")
