#!/usr/bin/env Rscript
# Recovery of sex-specific testosterone effects without using sex
#
# What: the testosterone mixture is fitted blind to sex, so its lower component should be the
#   women and its upper component the men. This script compares each G-MoQ component with the
#   GWAS that does use sex, for the unshrunk phenotype and the three shrunk arms, under both
#   weight designs:
#     female = beta_left  against the female-only REGENIE GWAS,
#     male   = beta_right against the male-only REGENIE GWAS,
#     sex    = -gamma     against the GWAS of chromosomal sex (02_sex_gwas.R). gamma is logit
#            P(lower = female component), so -gamma is membership of the upper (male)
#            component, which is how the paper reports it; the sex GWAS codes male = 1.
#   Comparisons: (A) per-variant concordance (correlation and slope over all variants and at
#   the reference's genome-wide significant variants); (B) genetic correlation by LD score
#   regression; (D) internal Mendelian randomisation in both directions (inverse-variance
#   weighted, fixed effect, through the origin); (F) one summary row per arm and component;
#   (L) genomic inflation of G-MoQ and REGENIE with the LDSC intercepts.
#   Two female/male references are used: 'combined' (sex-stratified REGENIE without LOCO
#   predictions, the primary per-variant comparison) and 'paper_matched' (same phenotype file,
#   covariates and LOCO predictions as Stage 1; the LDSC and Mendelian randomisation reference).
# Paper: Fig. 2a,b; Supplementary Fig. S8; Supplementary Table S3; Supplementary Note 2
# Input: per-variant Stage 2 results (02_mixtures/02_genomewide_stage2.R); sex-stratified
#   REGENIE GWAS; the sex GWAS; REGENIE lead variants; LD scores and 1000 Genomes European
#   genotypes for clumping (paths in config.yaml, section p08)
# Output: results/03_testosterone/recovery/{concordance, concordance_leads, ldsc_rg,
#   internal_mr, internal_mr_qc, table_s4_inputs, lambda_gc}.csv
# Runs: on a cluster with UK Biobank data, in about 3 hours (LD score regression dominates)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
p8 <- cfg$p08
ARMS <- p8$arms
N_ANALYSED <- cfg$cohorts$testosterone$n
GWS <- -log10(5e-8)
COMPONENTS <- data.table(component = c("female", "male", "sex"),
                         param = c("beta_left", "beta_right", "gamma"),
                         sign = c(1, 1, p8$gamma_sign_vs_sex))    # -1 for gamma: see header
out_dir <- "03_testosterone/recovery"

## 1. Reference GWAS ----------------------------------------------------------------------
read_combined <- function(arm, sex) {
  x <- fread(file.path(p8$regenie_combined, sprintf("testosterone_%s_%s_all.regenie", sex, arm)),
             select = c("ID", "ALLELE0", "ALLELE1", "N", "BETA", "SE", "LOG10P"), showProgress = FALSE)
  data.table(snp_id = x$ID, a1 = x$ALLELE1, a2 = x$ALLELE0, n = x$N, beta = x$BETA, se = x$SE, mlog10p = x$LOG10P)
}
read_ldsc_format <- function(path) {
  x <- fread(path, select = c("SNP", "A1", "A2", "N", "P", "BETA", "SE"), showProgress = FALSE)
  data.table(snp_id = x$SNP, a1 = x$A1, a2 = x$A2, n = x$N, beta = x$BETA, se = x$SE,
             mlog10p = -log10(pmax(x$P, 1e-320)))
}
REF <- list()
for (a in ARMS) for (s in c("female", "male")) {
  REF[[paste("combined", a, s)]] <- read_combined(a, s)
  REF[[paste("paper_matched", a, s)]] <- read_ldsc_format(file.path(p8$regenie_paper_raw, sprintf("regenie__%s_%s.tsv", a, s)))
}
REF[["sex"]] <- read_ldsc_format(p8$sex_reference)
REF_LAMBDA <- lapply(REF, function(r) lambda_from_z(r$beta / r$se))
ref_key <- function(reference, arm, component) if (reference == "sex") "sex" else paste(reference, arm, component)
references_for <- function(component) if (component == "sex") "sex" else c("combined", "paper_matched")
regenie_leads <- fread(p8$regenie_lead_snps)            # gwas_id, lead_snp (LD-clumped REGENIE leads)

# Inverse-variance weighted Mendelian randomisation slope through the origin (fixed effect),
# after aligning alleles and dropping palindromic variants.
ivw <- function(exposure, outcome) {
  m <- merge(exposure, outcome, by = "snp_id", suffixes = c(".x", ".y"))
  flip <- ifelse(m$a1.x == m$a1.y & m$a2.x == m$a2.y, 1, ifelse(m$a1.x == m$a2.y & m$a2.x == m$a1.y, -1, NA))
  palindromic <- paste0(m$a1.x, m$a2.x) %in% c("AT", "TA", "CG", "GC")
  m[, b_out := beta.y * flip]
  n_pal <- sum(palindromic & !is.na(flip))
  m <- m[!is.na(flip) & !palindromic & is.finite(beta.x) & beta.x != 0 & is.finite(b_out) & se.y > 0]
  if (!nrow(m)) return(list(res = data.table(nsnp = 0L), n_pal = n_pal))
  w <- 1 / m$se.y^2
  b <- sum(w * m$beta.x * m$b_out) / sum(w * m$beta.x^2); s <- sqrt(1 / sum(w * m$beta.x^2))
  list(res = data.table(nsnp = nrow(m), beta = b, se = s, p = 2 * pnorm(-abs(b / s)),
                        ci_low = b - 1.96 * s, ci_high = b + 1.96 * s,
                        beta_correlation = if (nrow(m) > 2) cor(m$beta.x, m$b_out) else NA_real_),
       n_pal = n_pal)
}

OUT <- list()
add <- function(name, x, design) {
  x <- as.data.table(x); x[, weight_design := design]; setcolorder(x, "weight_design")
  OUT[[name]] <<- c(OUT[[name]], list(x))
}

for (design in WEIGHT_DESIGNS) {
  G <- setNames(lapply(ARMS, function(a) read_stage2(cfg, design, a)), ARMS)
  ldsc_dir <- private_dir(cfg, "03_testosterone", "ldsc", design)

  ## 2. (A) Per-variant concordance -----------------------------------------------------
  conc <- list(); lead_rows <- list()
  for (a in ARMS) for (k in seq_len(nrow(COMPONENTS))) {
    cp <- COMPONENTS[k]
    for (rn in references_for(cp$component)) {
      ref <- REF[[ref_key(rn, a, cp$component)]]
      m <- harmonise_with_reference(G[[a]], cp$param, cp$sign, ref)
      for (set in c("all", "reference_p_lt_5e-8")) {
        mm <- if (set == "all") m else m[mlog10p >= GWS]
        conc[[length(conc) + 1L]] <- cbind(data.table(arm = a, component = cp$component, param = cp$param,
          sign_applied = cp$sign, reference = rn, snp_set = set, reference_n = median(ref$n)), concordance_stats(mm))
      }
      # Rows at the reference's LD-clumped lead variants (Fig. 2a,b).
      lead_gwas <- if (cp$component == "sex") "regenie__sex" else sprintf("regenie__%s_%s", a, cp$component)
      ids <- unique(regenie_leads[gwas_id == lead_gwas, lead_snp])
      if (length(ids)) lead_rows[[length(lead_rows) + 1L]] <- m[snp_id %in% ids,
        .(arm = a, component = cp$component, reference = rn, snp_id, chr, bp,
          beta_gmoq = b, se_gmoq = s, beta_regenie = beta_ref, se_regenie = se, mlog10p_regenie = mlog10p)]
    }
  }
  conc <- rbindlist(conc, fill = TRUE)
  add("concordance", conc, design)
  add("concordance_leads", rbindlist(lead_rows), design)

  ## 3. (B) Genetic correlation by LD score regression -----------------------------------
  # Female and male against the paper_matched sex-stratified GWAS; -gamma against the sex GWAS.
  rg_rows <- list(); h2_rows <- list()
  for (a in ARMS) for (k in seq_len(nrow(COMPONENTS))) {
    cp <- COMPONENTS[k]; g <- G[[a]]
    id <- sprintf("gmoq_%s_%s", a, cp$param)
    ref <- if (cp$component == "sex") p8$sex_reference_ldsc_id else sprintf("regenie__%s_%s", a, cp$component)
    ref_munged <- if (cp$component == "sex") p8$sex_reference_munged else
      file.path(p8$ldsc_munged_ref, paste0(ref, ".sumstats.gz"))
    munged <- file.path(ldsc_dir, paste0(id, ".sumstats.gz"))
    if (!file.exists(munged)) {
      raw <- data.table(SNP = g$snp_id, CHR = g$chr, BP = g$bp, A1 = g$effect_allele, A2 = g$other_allele,
                        N = N_ANALYSED, P = pmax(10^(-g[[paste0("mlog10p_", cp$param)]]), 1e-300),
                        BETA = g[[cp$param]] * cp$sign, SE = g[[paste0(cp$param, "_se")]])[is.finite(BETA) & is.finite(SE)]
      raw_file <- file.path(ldsc_dir, paste0(id, ".tsv"))
      fwrite(raw, raw_file, sep = "\t")
      ldsc_munge(cfg, raw_file, sub("\\.sumstats\\.gz$", "", munged))
      unlink(raw_file)
    }
    log <- ldsc_rg(cfg, munged, ref_munged, file.path(ldsc_dir, paste0(a, "__", cp$component,
                   if (cp$component == "sex") paste0("__", ref) else "")))
    v <- ldsc_rg_summary(log)
    rg_rows[[length(rg_rows) + 1L]] <- data.table(arm = a, component = cp$component, param = cp$param,
      sign_applied = cp$sign, reference = ref, rg = as.numeric(v$rg), rg_se = as.numeric(v$se),
      z = as.numeric(v$z), p = as.numeric(v$p), h2_ref_obs = as.numeric(v$h2_obs),
      h2_ref_obs_se = as.numeric(v$h2_obs_se), gcov_int = as.numeric(v$gcov_int),
      gcov_int_se = as.numeric(v$gcov_int_se), n_gmoq = N_ANALYSED, n_snps_ldsc = ldsc_n_snps(log))
    h1 <- ldsc_h2_block(log, 1); h2 <- ldsc_h2_block(log, 2)
    if (!is.null(h1)) h2_rows[[length(h2_rows) + 1L]] <- cbind(data.table(arm = a, component = cp$component,
      trait_side = "gmoq", ldsc_trait = id), h1)
    if (!is.null(h2)) h2_rows[[length(h2_rows) + 1L]] <- cbind(data.table(arm = a, component = cp$component,
      trait_side = "regenie_reference", ldsc_trait = ref), h2)
  }
  rg <- rbindlist(rg_rows); h2 <- rbindlist(h2_rows, fill = TRUE)
  add("ldsc_rg", rg, design)

  ## 4. (D) Internal Mendelian randomisation (female and male components) -----------------
  # REGENIE -> G-MoQ: the REGENIE lead variants as instruments, G-MoQ as outcome.
  # G-MoQ -> REGENIE: G-MoQ variants clumped with plink2 (1000 Genomes European; P < 5e-8,
  # secondary P < 1e-4, r2 < 0.001, 10 Mb), REGENIE (paper_matched) as outcome.
  mr_rows <- list(); mr_qc <- list()
  clump_dir <- private_dir(cfg, "03_testosterone", "clump", design)
  for (a in ARMS) for (k in 1:2) {
    cp <- COMPONENTS[k]; g <- G[[a]]
    gm <- data.table(snp_id = g$snp_id, a1 = g$effect_allele, a2 = g$other_allele, beta = g[[cp$param]],
                     se = g[[paste0(cp$param, "_se")]], P = pmax(10^(-g[[paste0("mlog10p_", cp$param)]]), 1e-300))
    clump_in <- file.path(clump_dir, sprintf("%s__%s.tsv", a, cp$param))
    fwrite(gm[, .(SNP = snp_id, P)], clump_in, sep = "\t")
    run_shell(sprintf("module load %s && plink2 --bfile %s --clump %s --clump-id-field SNP --clump-p-field P --clump-p1 %g --clump-p2 %g --clump-r2 %g --clump-kb %d --out %s",
                      p8$plink_module, p8$plink_bfile, clump_in, p8$clump$p1, p8$clump$p2, p8$clump$r2,
                      as.integer(p8$clump$kb), sub("\\.tsv$", "", clump_in)))
    clumps <- sub("\\.tsv$", ".clumps", clump_in)
    gmoq_instruments <- if (file.exists(clumps)) fread(clumps)$ID else character(0)
    reference <- REF[[paste("paper_matched", a, cp$component)]]
    regenie_instruments <- unique(regenie_leads[gwas_id == sprintf("regenie__%s_%s", a, cp$component), lead_snp])
    for (direction in c("REGENIE -> G-MoQ", "G-MoQ -> REGENIE")) {
      if (direction == "REGENIE -> G-MoQ") {
        exposure <- reference[snp_id %in% regenie_instruments]; outcome <- gm; inst <- regenie_instruments
        rule <- "reviewed REGENIE lead_snps.tsv (plink clump p<5e-8, r2 0.001, 10 Mb)"
      } else {
        exposure <- gm[snp_id %in% gmoq_instruments]; outcome <- reference; inst <- gmoq_instruments
        rule <- paste("re-clumped", design, "Stage 2 (plink2, same settings)")
      }
      z <- ivw(exposure[, .(snp_id, a1, a2, beta, se)], outcome[, .(snp_id, a1, a2, beta, se)])
      mr_rows[[length(mr_rows) + 1L]] <- cbind(data.table(arm = a, component = cp$component,
        direction = direction, reference = "paper_matched"), z$res)
      mr_qc[[length(mr_qc) + 1L]] <- data.table(arm = a, component = cp$component, direction = direction,
        n_instruments = length(inst), n_exposure_rows = nrow(exposure), n_palindromic_dropped = z$n_pal,
        n_harmonised = z$res$nsnp, selection_rule = rule)
    }
  }
  mr <- rbindlist(mr_rows, fill = TRUE)
  add("internal_mr", mr, design)
  add("internal_mr_qc", rbindlist(mr_qc), design)

  ## 5. (F) One summary row per arm and component (Supplementary Table S3) ----------------
  s3 <- list()
  for (a in ARMS) for (k in seq_len(nrow(COMPONENTS))) {
    cp <- COMPONENTS[k]; g <- G[[a]]
    row <- data.table(arm = a, component = cp$component, param = cp$param, n_snps = nrow(g),
                      n_gws = sum(g[[paste0("mlog10p_", cp$param)]] >= GWS),
                      lambda_gc = lambda_from_z(g[[cp$param]] / g[[paste0(cp$param, "_se")]]))
    pick <- function(ref, set, col) {
      v <- conc[arm == a & component == cp$component & reference == ref & snp_set == set][[col]]
      if (length(v)) v[1] else NA_real_
    }
    for (ref in c("combined", "paper_matched", "sex")) {
      row[[paste0("r_all_", ref)]] <- pick(ref, "all", "pearson_beta")
      row[[paste0("slope_all_", ref)]] <- pick(ref, "all", "slope_gmoq_on_regenie")
      row[[paste0("r_gws_", ref)]] <- pick(ref, "reference_p_lt_5e-8", "pearson_beta")
      row[[paste0("slope_gws_", ref)]] <- pick(ref, "reference_p_lt_5e-8", "slope_gmoq_on_regenie")
    }
    rr <- rg[arm == a & component == cp$component]
    row[, `:=`(ldsc_rg = if (nrow(rr)) rr$rg[1] else NA_real_, ldsc_rg_se = if (nrow(rr)) rr$rg_se[1] else NA_real_)]
    for (dirn in c("REGENIE -> G-MoQ", "G-MoQ -> REGENIE")) {
      mm <- mr[arm == a & component == cp$component & direction == dirn]
      tag <- if (dirn == "REGENIE -> G-MoQ") "mr_regenie_to_gmoq" else "mr_gmoq_to_regenie"
      row[[paste0(tag, "_beta")]] <- if (nrow(mm)) mm$beta[1] else NA_real_
      row[[paste0(tag, "_se")]] <- if (nrow(mm)) mm$se[1] else NA_real_
      row[[paste0(tag, "_nsnp")]] <- if (nrow(mm)) mm$nsnp[1] else NA_integer_
    }
    s3[[length(s3) + 1L]] <- row
  }
  add("table_s4_inputs", rbindlist(s3, fill = TRUE), design)

  ## 6. (L) Genomic inflation of G-MoQ and REGENIE, with LDSC intercepts ------------------
  # A raised lambda with an LDSC intercept near 1 (small ratio) reflects polygenic signal;
  # an intercept well above 1 would indicate inflation.
  lrows <- list()
  for (a in ARMS) for (k in seq_len(nrow(COMPONENTS))) {
    cp <- COMPONENTS[k]; g <- G[[a]]
    for (rn in references_for(cp$component)) {
      key <- ref_key(rn, a, cp$component)
      cr <- conc[arm == a & component == cp$component & reference == rn & snp_set == "all"]
      hg <- h2[arm == a & component == cp$component & trait_side == "gmoq"]
      hr <- h2[arm == a & component == cp$component & trait_side == "regenie_reference"]
      lrows[[length(lrows) + 1L]] <- data.table(arm = a, component = cp$component, param = cp$param,
        reference = rn, lambda_gc_gmoq_all = lambda_from_z(g[[cp$param]] / g[[paste0(cp$param, "_se")]]),
        n_gmoq = nrow(g), lambda_gc_regenie_all = REF_LAMBDA[[key]], n_regenie = nrow(REF[[key]]),
        lambda_gc_gmoq_common = cr$lambda_gc_gmoq_common, lambda_gc_regenie_common = cr$lambda_gc_regenie_common,
        n_common = cr$n, median_se_ratio_gmoq_over_regenie = cr$median_se_ratio_gmoq_over_regenie,
        ldsc_reference = if (nrow(hr)) hr$ldsc_trait[1] else NA_character_,
        ldsc_gmoq_lambda_gc = hg$ldsc_lambda_gc[1], ldsc_gmoq_mean_chi2 = hg$mean_chi2[1],
        ldsc_gmoq_intercept = hg$intercept[1], ldsc_gmoq_intercept_se = hg$intercept_se[1],
        ldsc_gmoq_ratio = hg$ratio[1], ldsc_gmoq_ratio_se = hg$ratio_se[1], ldsc_gmoq_h2 = hg$h2[1],
        ldsc_ref_lambda_gc = hr$ldsc_lambda_gc[1], ldsc_ref_mean_chi2 = hr$mean_chi2[1],
        ldsc_ref_intercept = hr$intercept[1], ldsc_ref_intercept_se = hr$intercept_se[1],
        ldsc_ref_ratio = hr$ratio[1], ldsc_ref_ratio_se = hr$ratio_se[1])
    }
  }
  add("lambda_gc", rbindlist(lrows, fill = TRUE), design)
  rm(G); invisible(gc())
}

## 7. Write -------------------------------------------------------------------------------
for (nm in names(OUT)) write_result(rbindlist(OUT[[nm]], fill = TRUE), out_dir, nm)
