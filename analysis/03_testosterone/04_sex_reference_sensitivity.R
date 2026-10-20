#!/usr/bin/env Rscript
# Membership against two GWAS of sex: the Stage 1 sample and the full cohort
#
# What: the testosterone membership parameter (-gamma, membership of the upper, male
#   component; gamma itself is logit P(lower component)) is compared with two GWAS of
#   chromosomal sex: the one on exactly the 399,710 Stage 1 participants (02_sex_gwas.R, the
#   reference used everywhere else) and the earlier one on all 462,869 European-ancestry
#   participants, which overlaps the Stage 1 sample by 86%. For each arm and weight design it
#   gives the per-variant concordance with both, and, for the design reported for each arm, the
#   genetic correlation with both side by side.
# Paper: Supplementary Table S3
# Input: per-variant Stage 2 results (02_mixtures/02_genomewide_stage2.R); munged -gamma
#   summary statistics from 03_compare_with_sex_stratified_gwas.R; both sex GWAS (raw and
#   munged; paths in config.yaml, section p08)
# Output: results/03_testosterone/sex_gwas/{concordance,ldsc_rg_side_by_side}.csv
# Runs: on a cluster with UK Biobank-derived summary statistics, in about 30 minutes

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
p8 <- cfg$p08
REPORTED <- setNames(mixture_design(cfg, p8$arms), p8$arms)   # hybrid for unshrunk, model otherwise
out_dir <- "03_testosterone/sex_gwas"
SEX_GWAS <- list(
  old_n462869 = list(raw = p8$sex_reference_full_cohort, munged = p8$sex_reference_full_cohort_munged),
  new_n399710 = list(raw = p8$sex_reference, munged = p8$sex_reference_munged))

## 1. Genetic correlation of -gamma with each sex GWAS (reported design per arm) ----------
rg_rows <- list()
for (a in p8$arms) {
  design <- REPORTED[[a]]
  ldsc_dir <- private_dir(cfg, "03_testosterone", "ldsc", design)
  gamma_munged <- file.path(ldsc_dir, sprintf("gmoq_%s_gamma.sumstats.gz", a))   # sign -1 applied
  for (ref in names(SEX_GWAS)) {
    log <- ldsc_rg(cfg, gamma_munged, SEX_GWAS[[ref]]$munged,
                   file.path(private_dir(cfg, "03_testosterone", "sex_gwas", "rg"), sprintf("%s__%s__%s", design, a, ref)))
    v <- ldsc_rg_summary(log)
    rg_rows[[length(rg_rows) + 1L]] <- data.table(weight_design = design, arm = a, sex_gwas = ref,
      rg = as.numeric(v$rg), rg_se = as.numeric(v$se), gcov_int = as.numeric(v$gcov_int))
  }
}
side_by_side <- dcast(rbindlist(rg_rows), weight_design + arm ~ sex_gwas, value.var = c("rg", "rg_se", "gcov_int"))
side_by_side[, rg_diff_new_minus_old := rg_new_n399710 - rg_old_n462869]
write_result(side_by_side, out_dir, "ldsc_rg_side_by_side")

## 2. Per-variant concordance of -gamma with each sex GWAS (both designs) -----------------
read_ldsc_format <- function(path) {
  x <- fread(path, select = c("SNP", "A1", "A2", "N", "P", "BETA", "SE"), showProgress = FALSE)
  data.table(snp_id = x$SNP, a1 = x$A1, a2 = x$A2, n = x$N, beta = x$BETA, se = x$SE,
             mlog10p = -log10(pmax(x$P, 1e-320)))
}
REF <- lapply(SEX_GWAS, function(s) read_ldsc_format(s$raw))
rows <- list()
for (design in WEIGHT_DESIGNS) for (a in p8$arms) {
  g <- read_stage2(cfg, design, a, select = c("snp_id", "chr", "bp", "effect_allele", "other_allele", "gamma", "gamma_se"))
  for (ref in names(REF)) {
    m <- harmonise_with_reference(g, "gamma", p8$gamma_sign_vs_sex, REF[[ref]])
    for (set in c("all", "reference_p_lt_5e-8")) {
      mm <- if (set == "all") m else m[mlog10p >= -log10(5e-8)]
      rows[[length(rows) + 1L]] <- cbind(data.table(weight_design = design, arm = a, reported = REPORTED[[a]] == design,
        component = "sex", sex_gwas = ref, snp_set = set, reference_n = median(REF[[ref]]$n)), concordance_stats(mm))
    }
  }
}
write_result(rbindlist(rows, fill = TRUE), out_dir, "concordance")
