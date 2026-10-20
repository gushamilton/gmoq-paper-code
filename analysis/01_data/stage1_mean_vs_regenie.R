#!/usr/bin/env Rscript
# Stage 1 mean effect against REGENIE
#
# What: a check that Stage 1 is an ordinary GWAS underneath. Projecting a quantile profile
#   on a constant (W = 1) gives the mean effect, which here is the plain average of the 45
#   quantile effects, m = (1/45) sum_t beta_t, with pooled SE sqrt(1' Sigma0 1 / 45^2 /
#   d_resid) (OLS on a constant W reduces to this exactly, so no Stage 2 fit is needed).
#   Every variant is matched to the trait's additive REGENIE GWAS by chromosome, position
#   and rsID with alleles identical or swapped, and per trait the script reports the slope,
#   intercept and correlation of m on the REGENIE beta and the ratio of their SEs, for all
#   variants and for REGENIE P < 5e-8 and P < 1e-5. The location effect delta_mu =
#   sum_t c_t beta_t (trapezoid weights, covering 0.96 of the distribution) is reported
#   alongside. Testosterone is compared with sex-combined REGENIE without LOCO, the only
#   sex-combined REGENIE run.
# Paper: Supplementary Fig. S12
# Input: Stage 1 quantile GWAS and the matching additive REGENIE GWAS of BMI, log BMI,
#   weight, testosterone x4, JCV VP1 and MCV VP1 (one file per chromosome)
# Output: results/01_data/stage1_mean_vs_regenie/{mean_vs_regenie_summary,
#   mean_vs_regenie_points,match_counts}.csv
# Runs: on a cluster with the GWAS files, about 30 minutes (traits run 4 at a time)

.repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(.repo, "R", "load.R"))
cfg <- read_config()
u <- tau_grid(cfg)
cw <- trap_weights(u)                    # location functional delta_mu
w1 <- rep(1 / length(u), length(u))      # W = 1 mean: the grid average
out <- "01_data/stage1_mean_vs_regenie"

TRAITS <- c("bmi", "log_bmi", "weight", "testosterone", "testosterone_shrink_0_2",
            "testosterone_shrink_0_5", "testosterone_shrink_0_7", "jcv_vp1_log_all", "mcv_vp1_log_all")

## 1. Per trait and chromosome: Stage 1 mean and delta_mu, matched to REGENIE
match_trait <- function(tr) {
  stage1 <- stage1_chr_paths(cfg, tr)[stage1_exists == TRUE][order(chr)]
  regenie <- regenie_files(cfg, tr, "comparator")[file.exists(path)]
  genome_wide_file <- nrow(regenie) == 1L && identical(regenie$layout, "combined")
  rows <- list(); counts <- list()
  for (j in seq_len(nrow(stage1))) {
    cc <- stage1$chr[j]
    st <- stage1_taus(stage1$stage1[j])
    x <- fread(cmd = sprintf("zcat %s", shQuote(stage1$stage1[j])),
               select = c("snp_id", "chr", "bp", "effect_allele", "other_allele", "d_resid", st$beta_cols),
               showProgress = FALSE)
    Sigma0 <- read_sigma0(stage1$sigma0[j], st$taus)
    B <- as.matrix(x[, st$beta_cols, with = FALSE])
    g <- data.table(snp_id = x$snp_id, chr = as.integer(x$chr), bp = as.integer(x$bp),
                    ea = x$effect_allele, oa = x$other_allele, d_resid = x$d_resid,
                    mean_w1 = as.numeric(B %*% w1),
                    mean_w1_se = sqrt(as.numeric(t(w1) %*% Sigma0 %*% w1) / x$d_resid),
                    delta_mu = as.numeric(B %*% cw),
                    delta_mu_se = sqrt(as.numeric(t(cw) %*% Sigma0 %*% cw) / x$d_resid))
    rg <- if (genome_wide_file) read_regenie(regenie$path, chr = cc) else read_regenie(regenie[chr == cc]$path)
    m <- merge(g, rg, by.x = c("chr", "bp", "snp_id"), by.y = c("CHROM", "GENPOS", "ID"))
    m[, sgn := allele_sign(ea, oa, ALLELE1, ALLELE0)]
    counts[[j]] <- data.table(trait = tr, chr = cc, n_stage1 = nrow(g), n_regenie = nrow(rg),
                              n_position_match = nrow(m), n_allele_ok = sum(!is.na(m$sgn)),
                              n_swapped = sum(m$sgn == -1, na.rm = TRUE))
    m <- m[!is.na(sgn)]
    m[, `:=`(beta_regenie = sgn * BETA, se_regenie = SE, log10p_regenie = LOG10P, n_regenie = N)]
    rows[[j]] <- m[, .(trait = tr, snp_id, chr, bp, effect_allele = ea, other_allele = oa, d_resid,
                       mean_w1, mean_w1_se, delta_mu, delta_mu_se, beta_regenie, se_regenie,
                       log10p_regenie, n_regenie)]
  }
  list(matched = rbindlist(rows), counts = rbindlist(counts))
}
matched <- parallel::mclapply(TRAITS, match_trait, mc.cores = 4, mc.preschedule = FALSE)
for (x in matched) if (inherits(x, "try-error")) stop(x)

## 2. Agreement summaries
summarise_agreement <- function(d, tr, subset) {
  if (nrow(d) < 10) return(NULL)
  fit <- lm(mean_w1 ~ beta_regenie, data = d)
  ratio <- d$mean_w1_se / d$se_regenie
  data.table(trait = tr, subset = subset, n = nrow(d),
             slope = unname(coef(fit)[2]), slope_se = unname(summary(fit)$coefficients[2, 2]),
             intercept = unname(coef(fit)[1]), pearson_r = cor(d$mean_w1, d$beta_regenie),
             slope_origin = sum(d$mean_w1 * d$beta_regenie) / sum(d$beta_regenie^2),
             slope_reverse = unname(coef(lm(beta_regenie ~ mean_w1, data = d))[2]),
             delta_mu_slope = unname(coef(lm(delta_mu ~ beta_regenie, data = d))[2]),
             delta_mu_pearson_r = cor(d$delta_mu, d$beta_regenie),
             median_se_ratio = median(ratio), se_ratio_p05 = unname(quantile(ratio, 0.05)),
             se_ratio_p95 = unname(quantile(ratio, 0.95)),
             median_n_regenie = median(d$n_regenie),
             n_sig_regenie = sum(d$log10p_regenie >= -log10(5e-8)),
             n_sig_gmoq_w1 = sum(2 * pnorm(-abs(d$mean_w1 / d$mean_w1_se)) < 5e-8))
}

## 3. Points for the figure: every variant with REGENIE P < 1e-5, plus a fixed random 5,000
## of the others per trait (the seed is set once, before the first trait)
set.seed(cfg$estimation$boot_seed)
summaries <- list(); points <- list(); counts <- list()
for (i in seq_along(TRAITS)) {
  M <- matched[[i]]$matched; tr <- TRAITS[i]; counts[[tr]] <- matched[[i]]$counts
  summaries[[length(summaries) + 1L]] <- summarise_agreement(M, tr, "all")
  summaries[[length(summaries) + 1L]] <- summarise_agreement(M[log10p_regenie >= -log10(5e-8)], tr, "regenie_p_lt_5e-8")
  summaries[[length(summaries) + 1L]] <- summarise_agreement(M[log10p_regenie >= 5], tr, "regenie_p_lt_1e-5")
  significant <- M[log10p_regenie >= 5]
  rest <- M[log10p_regenie < 5]
  points[[tr]] <- rbind(significant[, sampled := FALSE], rest[sample(.N, min(.N, 5000L))][, sampled := TRUE])
}
summary_table <- rbindlist(summaries)
point_table <- rbindlist(points)
print(summary_table[, .(trait, subset, n, slope = round(slope, 4), r = round(pearson_r, 4),
                        se_ratio = round(median_se_ratio, 4))])

write_result(summary_table, out, "mean_vs_regenie_summary")
num6 <- c("mean_w1", "mean_w1_se", "delta_mu", "delta_mu_se", "beta_regenie", "se_regenie", "log10p_regenie")
point_table[, (num6) := lapply(.SD, signif, 6), .SDcols = num6]
write_result(point_table[, .(trait, snp_id, chr, bp, effect_allele, other_allele, mean_w1, mean_w1_se,
                             delta_mu, delta_mu_se, beta_regenie, se_regenie, log10p_regenie, sampled)],
             out, "mean_vs_regenie_points")
write_result(rbindlist(counts), out, "match_counts")
