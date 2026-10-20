#!/usr/bin/env Rscript
# Thresholded comparators: GWAS of seropositivity and of IgG among seropositive participants
#
# What: the usual way to analyse antibody levels is to threshold them. This script runs the two
#   thresholded GWAS on the same 9,145 participants as the G-MoQ antibody analysis, for JCV VP1
#   and MCV VP1: seropositivity (the assay provider's serostatus flag, MFI above 250; logistic
#   REGENIE with Firth correction) and log IgG level among seropositive participants (linear
#   REGENIE). It reports both at the antibody lead and example variants, aligned to the G-MoQ
#   effect allele, with genomic inflation and sample counts. These are the comparators for the
#   membership (-gamma, upper-component membership; gamma = logit P(lower)) and level effects.
# Paper: Fig. 3; Supplementary Note 3
# Input: UK Biobank antibody phenotypes (log MFI and serostatus flags), covariates (age, sex,
#   ten principal components), keep lists, directly genotyped array (REGENIE step 1) and
#   imputed genotypes (step 2); the lead and example variants of
#   results/04_antibodies/decomposition/lead_decomposition.csv (paths in config.yaml, p09b)
# Output: results/04_antibodies/serostatus_gwas/{serostatus_leads,lambda_gc,cohort_counts}.csv;
#   phenotype and keep files in the working directory kept with the UK Biobank data (also read by 03_fut2_threshold.R)
# Runs: on a cluster with UK Biobank data. REGENIE step 2 runs per chromosome (8 threads); the
#   chromosomes are run in parallel on the cluster, here in a loop. A few hours in total.

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
pb <- cfg$p09b
TRAITS <- unlist(cfg$headline$antibody_traits)
stem_of <- function(tr) sub("_log_all$", "", tr)                      # jcv_vp1, mcv_vp1
pathogen_of <- function(tr) sub("_vp1_log_all$", "", tr)              # jcv, mcv
work <- private_dir(cfg, "04_antibodies", "serostatus")
out_dir <- "04_antibodies/serostatus_gwas"
threads <- 8L

## 1. Analysed cohort, serostatus and IgG among seropositives -----------------------------
# Cohort: the Stage 1 sample of the antibody scan, restricted to the keep list of the scan's
# REGENIE (participants complete for every antigen); it must be the 9,145 of the G-MoQ analysis.
# Serostatus: the assay provider's flag <pathogen>_seropos (any assay instance at or above 250
# MFI), as in the antibody analysis manifest.
to01 <- function(v) { v <- toupper(trimws(as.character(v)))
  ifelse(v %in% c("1", "TRUE", "T", "YES"), 1, ifelse(v %in% c("0", "FALSE", "F", "NO"), 0, NA_real_)) }
pheno <- fread(pb$phenotype_file, colClasses = list(character = c("FID", "IID")))
keep_all <- fread(pb$keep_all_complete, header = FALSE, colClasses = "character")[[1]]
cols <- list(); ids_all <- list(); counts <- list()
for (tr in TRAITS) {
  inp <- stage1_inputs(cfg, tr)
  s1 <- stage1_analysed_sample(inp$pheno, inp$trait, inp$covar, inp$keep, inp$loco, inp$bgen_sample)
  ids <- s1[FID %in% keep_all, .(FID, IID)]
  if (nrow(ids) != cfg$cohorts[[tr]]$n) stop(tr, ": analysed cohort has ", nrow(ids), " participants")
  x <- pheno[match(ids$FID, pheno$FID)]
  y <- as.numeric(x[[tr]])
  seropos <- to01(x[[paste0(pathogen_of(tr), "_seropos")]])
  if (anyNA(seropos) || !all(is.finite(y))) stop(tr, ": missing serostatus or antibody level")
  d <- data.table(FID = x$FID, IID = x$IID, s = seropos, lp = ifelse(seropos == 1, y, NA_real_))
  setnames(d, c("s", "lp"), paste0(stem_of(tr), c("_seropos", "_log_seropos")))
  cols[[tr]] <- d; ids_all[[tr]] <- ids
  fwrite(ids[which(seropos == 1)], file.path(work, paste0("keep_", tr, "_seropos.txt")), sep = " ", col.names = FALSE)
  counts[[tr]] <- data.table(trait = tr, n_analysed = nrow(ids), n_stage1_rule = nrow(s1),
                             n_seropositive = sum(seropos == 1), n_seronegative = sum(seropos == 0),
                             seroprevalence = mean(seropos == 1))
}
pheno_file <- file.path(work, "phenotypes_serostatus.txt")
fwrite(Reduce(function(a, b) merge(a, b, by = c("FID", "IID"), all = TRUE), cols), pheno_file, sep = " ", na = "NA", quote = FALSE)
keep_analysed <- file.path(work, "keep_analysed.txt")
fwrite(unique(rbindlist(ids_all)), keep_analysed, sep = " ", col.names = FALSE, quote = FALSE)
counts <- rbindlist(counts)

## 2. REGENIE -----------------------------------------------------------------------------
# Step 1 on the directly genotyped array after plink2 quality control in the tested sample;
# step 2 on imputed HapMap3 variants. Seropositivity: --bt with Firth correction for both
# antigens together; IgG among seropositives: --qt, one antigen at a time.
covar <- paste("--covarFile", pb$covar_file, "--covarCol", paste(unlist(pb$covariates), collapse = ","))
bt_cols <- paste(paste0(stem_of(TRAITS), "_seropos"), collapse = ",")
regenie_step1 <- function(dir, keep, pheno_args, type, prefix) {
  if (file.exists(file.path(dir, paste0(prefix, "_pred.list")))) return(file.path(dir, paste0(prefix, "_pred.list")))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  run_shell(paste("module load", pb$plink_module, "&& plink2 --bfile", pb$step1_bed, "--keep", keep, pb$step1_qc,
                  "--threads", threads, "--write-snplist --write-samples --no-id-header --out", file.path(dir, "qc_pass")))
  run_shell(paste("module load", pb$regenie_module, "&& regenie --step 1 --bed", pb$step1_bed, "--phenoFile", pheno_file,
                  pheno_args, covar, "--keep", keep, "--extract", file.path(dir, "qc_pass.snplist"), type,
                  "--bsize 1000 --lowmem --lowmem-prefix", file.path(dir, paste0("tmp_", prefix)),
                  "--threads", threads, "--out", file.path(dir, prefix)))
  file.path(dir, paste0(prefix, "_pred.list"))
}
regenie_step2 <- function(dir, keep, pred, pheno_args, type, prefix) for (chr in 1:22) {
  cc <- sprintf("%02d", chr); od <- file.path(dir, paste0("chr", cc))
  if (length(Sys.glob(file.path(od, paste0(prefix, "_chr", cc, "_*.regenie"))))) next
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  run_shell(paste("module load", pb$regenie_module, "&& regenie --step 2 --bgen", file.path(pb$bgen_dir, paste0("data.chr", cc, ".bgen")),
                  "--sample", pb$bgen_sample, "--pred", pred, "--phenoFile", pheno_file, pheno_args, covar,
                  "--keep", keep, "--extract", gsub("NN", cc, pb$extract_template), type, "--bsize", pb$step2_bsize,
                  "--threads", threads, "--out", file.path(od, paste0(prefix, "_chr", cc))))
}
pred_bt <- regenie_step1(file.path(work, "step1_bt"), keep_analysed, paste("--phenoColList", bt_cols), "--bt", "seropos_step1")
regenie_step2(file.path(work, "step2_bt"), keep_analysed, pred_bt, paste("--phenoColList", bt_cols), pb$bt_flags, "seropos")
for (tr in TRAITS) {
  keep_pos <- file.path(work, paste0("keep_", tr, "_seropos.txt"))
  pheno_args <- paste("--phenoCol", paste0(stem_of(tr), "_log_seropos"))
  pred_qt <- regenie_step1(file.path(work, paste0("step1_qt_", tr)), keep_pos, pheno_args, "--qt", "log_seropos_step1")
  regenie_step2(file.path(work, paste0("step2_qt_", tr)), keep_pos, pred_qt, pheno_args, "--qt", "log_seropos")
}

## 3. Results at the lead and example variants --------------------------------------------
LD <- as.data.table(read_result("04_antibodies/decomposition", "lead_decomposition"))
SNPS <- unique(LD[trait %in% TRAITS, .(trait, snp_id, chr, bp, effect_allele, other_allele, is_example)])
ANALYSES <- data.table(analysis = c("seropositivity", "seropositive_igg"),
                       beta_scale = c("log_odds_seropositive", "log_mfi_among_seropositives"))
regenie_file <- function(tr, analysis, chr) {
  st <- stem_of(tr); cc <- sprintf("%02d", chr)
  if (analysis == "seropositivity") file.path(work, "step2_bt", paste0("chr", cc), sprintf("seropos_chr%s_%s_seropos.regenie", cc, st))
  else file.path(work, paste0("step2_qt_", tr), paste0("chr", cc), sprintf("log_seropos_chr%s_%s_log_seropos.regenie", cc, st))
}
read_regenie <- function(f) {
  hdr <- names(fread(f, nrows = 0L))
  x <- fread(f, select = intersect(c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "N_CASES",
                                     "N_CONTROLS", "TEST", "BETA", "SE", "LOG10P"), hdr), showProgress = FALSE)
  if ("TEST" %in% names(x)) x <- x[TEST == "ADD"]
  x
}
leads <- list(); lam <- list(); n_regenie <- list()
for (tr in TRAITS) {
  ct <- counts[trait == tr]
  for (a in seq_len(nrow(ANALYSES))) {
    an <- ANALYSES$analysis[a]; binary <- an == "seropositivity"
    G <- rbindlist(lapply(1:22, function(chr) read_regenie(regenie_file(tr, an, chr))), fill = TRUE)
    G[, chisq := qchisq(-LOG10P * log(10), 1, lower.tail = FALSE, log.p = TRUE)]
    n_med <- stats::median(G$N, na.rm = TRUE)
    if (n_med != if (binary) ct$n_analysed else ct$n_seropositive) stop(tr, " / ", an, ": REGENIE N differs from the cohort")
    n_regenie[[paste(tr, an)]] <- n_med
    lam[[length(lam) + 1L]] <- data.table(trait = tr, analysis = an, n_snps = nrow(G),
      lambda_gc = stats::median(G$chisq, na.rm = TRUE) / qchisq(0.5, 1),
      n_gws = sum(G$LOG10P >= -log10(5e-8), na.rm = TRUE), n_p_lt_1e5 = sum(G$LOG10P >= 5, na.rm = TRUE),
      n_median = n_med, n_cases = if (binary) ct$n_seropositive else NA_integer_,
      n_controls = if (binary) ct$n_seronegative else NA_integer_)
    L <- SNPS[trait == tr]
    m <- match(L$snp_id, G$ID)
    flip <- allele_sign(L$effect_allele, L$other_allele, G$ALLELE1[m], G$ALLELE0[m])
    has_cc <- all(c("N_CASES", "N_CONTROLS") %in% names(G))
    leads[[length(leads) + 1L]] <- data.table(trait = tr, analysis = an, snp_id = L$snp_id, chr = L$chr, bp = L$bp,
      effect_allele = L$effect_allele, other_allele = L$other_allele,
      beta = G$BETA[m] * flip, se = G$SE[m], mlog10p = G$LOG10P[m], n = G$N[m],
      n_cases = if (!binary) NA_real_ else if (has_cc) G$N_CASES[m] else ct$n_seropositive,
      n_controls = if (!binary) NA_real_ else if (has_cc) G$N_CONTROLS[m] else ct$n_seronegative,
      a1freq = ifelse(flip == 1, G$A1FREQ[m], 1 - G$A1FREQ[m]), allele_flipped = flip == -1,
      is_example = L$snp_id %in% unlist(cfg$p09$example_loci), beta_scale = ANALYSES$beta_scale[a])
    rm(G); invisible(gc())
  }
}
counts[, `:=`(n_regenie_seropositivity = unlist(n_regenie[paste(trait, "seropositivity")]),
              n_regenie_seropositive_igg = unlist(n_regenie[paste(trait, "seropositive_igg")]))]
write_result(rbindlist(leads), out_dir, "serostatus_leads")
write_result(rbindlist(lam), out_dir, "lambda_gc")
write_result(counts, out_dir, "cohort_counts")
