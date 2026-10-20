#!/usr/bin/env Rscript
# Pooled Stage 1 standard errors against the exact per-individual HC0 sandwich (real data)
#
# What: Stage 1 reports the standard error of each quantile effect from a pooled covariance,
#   se_pooled = sqrt(Sigma0[t,t] / D), Sigma0 = Psi'Psi / n, D = sum of squared residualised
#   genotype, which assumes the genotype does not change the spread of the influence functions.
#   This script compares it with the exact heteroscedasticity-robust (HC0) sandwich,
#   var_HC0 = sum_i g_i^2 (psi_i,t - g_i beta_t)^2 / D^2, on real data: BMI in UK Biobank
#   (the Stage 1 sample, n = 461,320), chromosome 22, 500 randomly chosen HapMap3 variants,
#   all 45 quantile levels. The recentred influence functions are rebuilt from the quantiles
#   and kernel densities Stage 1 itself used (rif_nuisance.json), with the same covariates and
#   leave-one-chromosome-out prediction, and the rebuilt estimates and pooled standard errors
#   are checked against the Stage 1 output before the HC0 comparison is made.
# Paper: Methods (Stage 1: "Against an exact per-individual HC0 sandwich the pooled standard
#   errors had a median ratio of ...")
# Input: Stage 1 output for BMI chromosome 22; phenotype and covariate files; REGENIE LOCO
#   prediction for BMI; UK Biobank imputed BGEN for chromosome 22 (via plink2 --export A)
# Output: <out_dir>/se_ratio_summary.csv, se_ratio_by_tau.csv (aggregate; copied to
#   results/01_data/stage1_se_check/). Result: median ratio 1.000, 5-95% range 0.995-1.005.
# Runs: on a cluster with UK Biobank data, one core, about 6 GB, under 30 minutes
# Usage: Rscript stage1_se_check.R <out_dir>

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1]] else stop("usage: stage1_se_check.R <out_dir>")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## User settings: the BMI Stage 1 output for chromosome 22 and its inputs (authors' paths).
stage1_dir <- "/path/to/stage1/bmi/chr22"                      # scan.stage1.tsv.gz, rif_nuisance.json
pheno_file <- "/path/to/phenotypes.txt"                        # FID IID bmi ...
covar_file <- "/path/to/covariates.txt"                        # FID IID age sex PC1-PC10
loco_file  <- "/path/to/regenie_step1/bmi_step1_1.loco"        # REGENIE LOCO predictions
keep_file  <- "/path/to/bmi_stage1_sample.txt"                 # FID IID of the Stage 1 sample
bgen <- "/path/to/ukb_imp_chr22.bgen"
bgen_sample <- "/path/to/ukb_imp.sample"
n_variants <- 500L
set.seed(20261020)

## 1. Sample, phenotype, covariates and the chromosome-22 LOCO prediction --------------------
keep <- fread(keep_file, header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
ph <- fread(pheno_file, select = c("FID", "IID", "bmi"), colClasses = list(character = c("FID", "IID")))
cv <- fread(covar_file, colClasses = list(character = c("FID", "IID")))
loco <- fread(loco_file)                                  # row per chromosome, column per FID_IID
loco22 <- data.table(id = names(loco)[-1], loco = as.numeric(unlist(loco[loco[[1]] == 22, -1])))
d <- merge(keep, ph, by = c("FID", "IID"))
d <- merge(d, cv, by = c("FID", "IID"))
d[, id := paste(FID, IID, sep = "_")]
d <- merge(d, loco22, by = "id")
d <- d[complete.cases(d)]
cat("analysis sample:", nrow(d), "\n")

## 2. Recentred influence functions, exactly as Stage 1 built them ---------------------------
nuis <- fromJSON(file.path(stage1_dir, "rif_nuisance.json"))
stopifnot(nrow(d) == nuis$n)                               # same sample as Stage 1
taus <- as.numeric(nuis$taus); q <- as.numeric(nuis$q_tau); fhat <- as.numeric(nuis$fhat)
y <- d$bmi
Psi <- sapply(seq_along(taus), function(t) q[t] + (taus[t] - (y <= q[t])) / fhat[t])
X <- cbind(1, as.matrix(d[, c("age", "sex", paste0("PC", 1:10), "loco"), with = FALSE]))
qrX <- qr(X)
Psi_r <- qr.resid(qrX, Psi)                                # covariate-residualised RIF
Sigma0_diag <- colSums(Psi_r^2) / nrow(Psi_r)

## 3. Dosages for 500 random chromosome-22 variants of the Stage 1 scan ---------------------
s1 <- fread(cmd = paste("zcat", shQuote(file.path(stage1_dir, "scan.stage1.tsv.gz"))))
pick <- s1[sample(.N, n_variants)]
writeLines(pick$snp_id, file.path(out_dir, "variants.txt"))
plink_out <- file.path(out_dir, "dosages")
system2("plink2", c("--memory", "3000", "--bgen", bgen, "ref-first", "--sample", bgen_sample,
                    "--extract", file.path(out_dir, "variants.txt"), "--keep", keep_file,
                    "--rm-dup", "force-first", "--export", "A", "--out", plink_out))
G <- fread(paste0(plink_out, ".raw"), colClasses = list(character = c("FID", "IID")))
G[, id := paste(FID, IID, sep = "_")]
G <- G[match(d$id, id)]
stopifnot(identical(G$id, d$id))
gcols <- setdiff(names(G), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "id"))

## 4. Per variant: estimate, pooled and HC0 standard errors -----------------------------------
res <- rbindlist(lapply(gcols, function(gc) {
  g <- G[[gc]]; g[is.na(g)] <- mean(g, na.rm = TRUE)
  g_r <- qr.resid(qrX, g)
  D <- sum(g_r^2)
  beta <- colSums(g_r * Psi_r) / D
  se_pooled <- sqrt(Sigma0_diag / D)
  E <- Psi_r - outer(g_r, beta)
  se_hc0 <- sqrt(colSums((g_r^2) * E^2)) / D
  snp <- sub("_[ACGT]+$", "", gc)
  data.table(snp_id = snp, tau = taus, beta = beta, se_pooled = se_pooled, se_hc0 = se_hc0)
}))

## 5. Check the rebuild against Stage 1, then summarise the ratio -----------------------------
s1_long <- melt(pick, id.vars = "snp_id", measure.vars = patterns(beta = "^beta_tau_", se = "^se_tau_"),
                variable.name = "t")
s1_long[, tau := taus[as.integer(t)]]
chk <- merge(res, s1_long[, .(snp_id, tau, beta_s1 = beta, se_s1 = se)], by = c("snp_id", "tau"))
check <- chk[, .(n_variant_tau = .N,
                 cor_abs_beta_vs_stage1 = cor(abs(beta), abs(beta_s1)),
                 median_se_pooled_over_stage1 = median(se_pooled / se_s1))]
print(check)
ratio <- res[, se_pooled / se_hc0]
summary_tab <- data.table(trait = "bmi", chromosome = 22L, n = nrow(d), n_variants = length(gcols),
                          n_variant_tau = length(ratio), median_ratio = median(ratio),
                          p05 = quantile(ratio, 0.05), p95 = quantile(ratio, 0.95),
                          p01 = quantile(ratio, 0.01), p99 = quantile(ratio, 0.99),
                          check[, .(cor_abs_beta_vs_stage1, median_se_pooled_over_stage1)])
by_tau <- res[, .(median_ratio = median(se_pooled / se_hc0), p05 = quantile(se_pooled / se_hc0, 0.05),
                  p95 = quantile(se_pooled / se_hc0, 0.95)), by = tau]
fwrite(summary_tab, file.path(out_dir, "se_ratio_summary.csv"))
fwrite(by_tau, file.path(out_dir, "se_ratio_by_tau.csv"))
print(summary_tab)
unlink(paste0(plink_out, c(".raw", ".log")))
