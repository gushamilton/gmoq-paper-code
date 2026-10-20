#!/usr/bin/env Rscript
# Baseline quantiles and density of each trait on the Stage 1 sample
#
# What: Stage 2 needs, for every trait, the population quantile function Q(u) and density
#   f(Q(u)) at the 45 quantile levels: the moment and model weights are built from them, and
#   the Stage 1 RIF divided by the same kernel density. This script recomputes them on
#   exactly the participants of each Stage 1 GWAS, with the Stage 1 package's own quantile
#   rule (R type 8) and Silverman Gaussian kernel density, so they are the values the RIF
#   used. It checks the rebuilt N against the Stage 1 N and, where the Stage 1 run saved its
#   quantiles and density (rif_nuisance.json), that they agree on every chromosome.
#   mu_trap = sum_t c_t Q_t / sum_t c_t (trapezoid weights over 0.02-0.98) is reported
#   because the location functional is centred on it.
# Paper: inputs to Stage 2 for every trait (Figs. 1-5; Supplementary Figs. S10-S14)
# Input: the Stage 1 sample of each trait (phenotype, covariate, keep and LOCO files named
#   in results/01_data/sample_sizes/cohort_funnel.csv); one BGEN file and the sample file
# Output: results/01_data/baseline_quantiles/{baseline_quantiles,baseline_summary}.csv
# Runs: on a cluster with UK Biobank data, about 10 minutes (traits run 4 at a time)

.repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(.repo, "R", "load.R"))
cfg <- read_config()

python <- cfg$package$stage1$python                  # Python with the gmoq package
bgen <- gsub("{chr}", "22", cfg$ukb$bgen, fixed = TRUE) # any chromosome: only its sample list is read
sample_file <- cfg$ukb$sample
helper <- repo_path("analysis", "01_data", "baseline_sample.py")
work <- file.path(tempdir(), "baseline")             # per-trait JSON of aggregates
dir.create(work, showWarnings = FALSE)

u <- tau_grid(cfg)
cw <- trap_weights(u)

## 1. Traits and their Stage 1 inputs
# Traits whose Stage 1 run did not complete (no covariance metadata) are skipped.
cohorts <- as.data.table(read_result("01_data/sample_sizes", "cohort_funnel"))
cohorts <- cohorts[!is.na(n_sigma0_meta) & !is.na(covariates) & nzchar(covariates)]

## 2. Rebuild each Stage 1 sample and evaluate Q(u) and f(Q(u)) (Python, package loaders)
run_trait <- function(i) {
  r <- cohorts[i]
  out <- file.path(work, paste0(r$trait, ".json"))
  loco <- if (is.na(r$loco_file_chr1) || !nzchar(r$loco_file_chr1)) "none" else r$loco_file_chr1
  system2(python, c(helper, "--trait", r$trait, "--pheno", r$pheno_file, "--pheno-col", r$pheno_col,
                    "--covar", r$covar_file, "--covar-cols", shQuote(r$covariates), "--keep", r$keep_file,
                    "--loco", loco, "--bgen", bgen, "--sample", sample_file,
                    "--density-method", r$density_method, "--out", out))
  out
}
json_files <- unlist(parallel::mclapply(seq_len(nrow(cohorts)), run_trait, mc.cores = 4, mc.preschedule = FALSE))

## 3. Tables, with agreement checks
quantiles <- list(); summary_rows <- list()
for (k in seq_len(nrow(cohorts))) {
  r <- cohorts[k]
  j <- jsonlite::read_json(json_files[k], simplifyVector = TRUE)
  q <- j$q_tau; f <- j$fhat

  # The quantiles, density and mean that each chromosome's Stage 1 run saved, if any.
  saved <- na.omit(stage1_chr_paths(cfg, r$trait)$rif_nuisance)
  dq <- df <- dm <- NA_real_; n_ok <- NA
  if (length(saved)) {
    L <- lapply(saved, jsonlite::read_json, simplifyVector = TRUE)
    dq <- max(vapply(L, function(z) max(abs(z$q_tau - q)), 0))
    df <- max(vapply(L, function(z) max(abs(z$fhat - f) / f), 0))
    dm <- max(vapply(L, function(z) abs(z$mean - j$mean), 0))
    n_ok <- all(vapply(L, function(z) z$n == j$n, TRUE))
  }
  # The earlier trait-survey quantiles used R type 7; compared for information only.
  survey_q <- file.path(dirname(cfg$cohorts$survey$root), "qtau", paste0(r$trait, "_qtau.csv"))
  d7 <- if (file.exists(survey_q)) max(abs(fread(survey_q)$q - j$q_tau_type7)) else NA_real_

  quantiles[[k]] <- data.table(trait = r$trait, tau = round(u, 3), q = q, fhat = f, q_type7 = j$q_tau_type7)
  summary_rows[[k]] <- data.table(
    trait = r$trait, group = r$group, n_keep = j$n_keep, n_valid_pheno = j$n_valid_pheno,
    n_absent_loco = j$n_absent_loco, n_after_loco = j$n_after_loco, n = j$n,
    n_stage1 = r$n_analysed, n_matches_stage1 = j$n == r$n_analysed,
    mean = j$mean, sd = j$sd, sd_pop = j$sd_pop, bandwidth = j$bandwidth,
    mu_trap = sum(cw * q) / sum(cw), q_050 = q[which.min(abs(u - 0.5))],
    n_distinct = j$n_distinct, quantile_method = j$quantile_method, density_method = j$density_method,
    n_rif_nuisance_files = length(saved), rif_nuisance_n_ok = n_ok,
    max_abs_q_diff_rif_nuisance = dq, max_rel_fhat_diff_rif_nuisance = df,
    abs_mean_diff_rif_nuisance = dm, max_abs_q7_diff_old_survey_qtau = d7)
}
baseline_quantiles <- rbindlist(quantiles)
baseline_summary <- rbindlist(summary_rows)
print(baseline_summary[, .(trait, n, n_stage1, n_matches_stage1, mean = signif(mean, 5),
                           dq = signif(max_abs_q_diff_rif_nuisance, 3),
                           df = signif(max_rel_fhat_diff_rif_nuisance, 3))])

write_result(baseline_quantiles, "01_data/baseline_quantiles", "baseline_quantiles")
write_result(baseline_summary, "01_data/baseline_quantiles", "baseline_summary")

# The rebuilt sample must be the Stage 1 sample, and the saved quantiles must agree.
stopifnot(all(baseline_summary$n_matches_stage1),
          all(is.na(baseline_summary$max_abs_q_diff_rif_nuisance) | baseline_summary$max_abs_q_diff_rif_nuisance <= 1e-9),
          all(is.na(baseline_summary$max_rel_fhat_diff_rif_nuisance) | baseline_summary$max_rel_fhat_diff_rif_nuisance <= 1e-9))
