#!/usr/bin/env Rscript
# Lead variant sets
#
# What: defines the lead variants of every trait by one clumping rule. Candidates are the
#   variants with P < 5e-8 in the trait's additive REGENIE GWAS; plink2 clumps them against
#   the 1000 Genomes European reference (--clump-p1 5e-8 --clump-p2 1e-3 --clump-r2 0.01
#   --clump-kb 10000) and the index variants are the leads ("rule A"). This is the rule
#   that gave the 639 BMI leads, and the script checks that it reproduces them; log BMI
#   uses the same 639. Trait-survey traits have no REGENIE GWAS, so the same clumping is
#   applied to the Stage 1 location effect delta_mu = sum_t c_t beta_t (trapezoid weights
#   c), with SE sqrt(c' Sigma0 c / d_resid) and P = 2 Phi(-|z|). A stricter rule
#   ("rule B": --clump-p2 1e-4 --clump-r2 0.001, as used by the earlier trait survey) is
#   kept as a sensitivity set, and the sex-stratified testosterone leads are read from their
#   REGENIE lead table. Finally the expected FUT2 and HLA loci are looked up among the
#   JCV/MCV leads.
# Paper: lead variants of Figs. 1-6 and Supplementary Figs. S8-S14; Supplementary Table S7
# Input: REGENIE step 2 GWAS (BMI, weight, testosterone x4, JCV/MCV VP1); Stage 1 quantile
#   GWAS of the survey traits and weight; the 1000 Genomes European plink files
# Output: results/01_data/lead_variants/{leads,lead_sets,antibody_expected_loci}.csv
# Runs: on a cluster with the GWAS files, about 30 minutes (traits run 4 at a time)

.repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(.repo, "R", "load.R"))
cfg <- read_config()

plink2 <- cfg$tools$plink2                 # plink2 v2.0.0-a.6.8LM
ld_reference <- cfg$tools$ld_reference     # 1000 Genomes phase 3 European, plink prefix
work <- file.path(tempdir(), "clump")      # summary-statistic and plink files
dir.create(work, showWarnings = FALSE)
cores <- 4

RULES <- list(A = list(p1 = "5e-8", p2 = "1e-3", r2 = "0.01", kb = "10000"),
              B = list(p1 = "5e-8", p2 = "1e-4", r2 = "0.001", kb = "10000"))
GWS <- -log10(5e-8)
u <- tau_grid(cfg)
cw <- trap_weights(u)

## 1. Clumping
# sumstats: data.table(SNP, CHR, BP, P, LOG10P). Returns the index variant IDs.
run_clump <- function(sumstats, tag, rule) {
  pr <- RULES[[rule]]
  ss <- file.path(work, paste0(tag, ".sumstats"))
  fwrite(sumstats, ss, sep = "\t", quote = FALSE)
  prefix <- file.path(work, paste0(tag, "_clump"))
  rc <- system2(plink2, c("--bfile", ld_reference, "--clump", ss, "--clump-id-field", "SNP",
                          "--clump-p-field", "P", "--clump-p1", pr$p1, "--clump-p2", pr$p2,
                          "--clump-r2", pr$r2, "--clump-kb", pr$kb, "--threads", "2", "--out", prefix),
                stdout = paste0(prefix, ".log"), stderr = paste0(prefix, ".log"))
  clumps <- paste0(prefix, ".clumps")
  if (rc != 0 || !file.exists(clumps)) {
    if (nrow(sumstats) == 0) return(character())
    stop("plink2 clumping failed for ", tag)
  }
  cl <- fread(clumps, showProgress = FALSE)
  cl[[if ("ID" %in% names(cl)) "ID" else names(cl)[3]]]
}

## 2. Rule A on a REGENIE GWAS
# The candidate P is written as the original awk pipeline wrote it, exp(-LOG10P log 10)
# printed with %.6g (P below the double range prints as 0), so plink2 orders ties as it did.
awk_p <- function(log10p) sprintf("%.6g", exp(-as.numeric(log10p) * log(10)))
regenie_candidates <- function(files) {
  rbindlist(lapply(files, function(f) {
    cmd <- sprintf(paste0("awk 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; print \"SNP\\tCHR\\tBP\\tLOG10P\"; next} ",
                          "$h[\"LOG10P\"] >= %.15g {print $h[\"ID\"]\"\\t\"$h[\"CHROM\"]\"\\t\"$h[\"GENPOS\"]\"\\t\"$h[\"LOG10P\"]}' %s"),
                   GWS, shQuote(f))
    fread(cmd = cmd, colClasses = "character", showProgress = FALSE)
  }))
}
lead_regenie <- function(trait, files, set, description) {
  cand <- regenie_candidates(files)
  ids <- run_clump(cand[, .(SNP, CHR, BP, P = awk_p(LOG10P), LOG10P)], set, "A")
  leads <- unique(cand[SNP %in% ids, .(snp_id = SNP, chr = as.integer(CHR), bp = as.integer(BP),
                                      log10p = as.numeric(LOG10P))], by = "snp_id")
  leads[, `:=`(lead_set = set, trait = trait, rule = "A_regenie", p_source = description)]
  list(leads = leads, n_candidates = nrow(cand))
}

# Collect lead sets (leads.csv) and one summary row per set (lead_sets.csv).
sets <- list(); summary_rows <- list()
add_set <- function(res, set, trait, rule, primary, source, extra = list()) {
  leads <- res$leads
  leads[, primary := primary]
  sets[[set]] <<- leads
  summary_rows[[set]] <<- as.data.table(c(
    list(lead_set = set, trait = trait, rule = rule, primary = primary, source = source,
         n_candidates = res$n_candidates, n_leads = nrow(leads)),
    RULES[[substr(rule, 1, 1)]], extra))
}

## 3. BMI: the original REGENIE GWAS, and the check that the rule returns the 639 leads
bmi_file <- regenie_files(cfg, "bmi", "leads")$path
bmi <- lead_regenie("bmi", bmi_file, "bmi_A_regenie", "REGENIE combined bmi_all.regenie (--ignore-pred)")
reference_639 <- fread(cfg$leads$bmi$source)
bmi_reproduced <- setequal(bmi$leads$snp_id, reference_639$ID)
add_set(bmi, "bmi_A_regenie", "bmi", "A_regenie", TRUE, bmi_file,
        list(reproduces = bmi_reproduced, n_reference = nrow(reference_639),
             n_only_new = length(setdiff(bmi$leads$snp_id, reference_639$ID)),
             n_only_reference = length(setdiff(reference_639$ID, bmi$leads$snp_id)),
             reference = cfg$leads$bmi$source))
# Log BMI uses the same 639 variants (aligned to the BMI-increasing allele in lead_profiles.R).
log_bmi <- copy(bmi$leads)[, `:=`(lead_set = "log_bmi_A_regenie", trait = "log_bmi",
                                  p_source = "inherited from bmi_A_regenie")]
add_set(list(leads = log_bmi, n_candidates = bmi$n_candidates), "log_bmi_A_regenie", "log_bmi",
        "A_regenie", TRUE, "inherited: same 639 as BMI")

## 4. Weight, testosterone x4 and the antibodies: rule A on their REGENIE GWAS
regenie_description <- c(
  weight = "REGENIE GWAS of height-adjusted weight with PC1-PC10 (16 covariates)",
  testosterone = "sex-combined REGENIE, --ignore-pred", testosterone_shrink_0_2 = "sex-combined REGENIE, --ignore-pred",
  testosterone_shrink_0_5 = "sex-combined REGENIE, --ignore-pred", testosterone_shrink_0_7 = "sex-combined REGENIE, --ignore-pred",
  jcv_vp1_log_all = "44-antigen scan REGENIE", mcv_vp1_log_all = "44-antigen scan REGENIE")
regenie_sets <- parallel::mclapply(names(regenie_description), function(tr) {
  files <- regenie_files(cfg, tr, "leads")[file.exists(path)]$path
  list(trait = tr, files = files,
       res = lead_regenie(tr, files, paste0(tr, "_A_regenie"), regenie_description[[tr]]))
}, mc.cores = cores)
for (x in regenie_sets) {
  if (inherits(x, "try-error")) stop(x)
  add_set(x$res, paste0(x$trait, "_A_regenie"), x$trait, "A_regenie", TRUE,
          if (length(x$files) == 1) x$files else paste0(dirname(x$files[1]), " (22 files)"))
}

## 5. Survey traits and weight: clumping on the Stage 1 location effect delta_mu
stage1_complete <- function(tr) {
  f <- stage1_chr_paths(cfg, tr)
  nrow(f) == 22 && all(f$stage1_exists & file.exists(f$sigma0) & file.exists(f$sigma0_meta))
}
deltamu_leads <- function(tr) {
  files <- stage1_chr_paths(cfg, tr)[stage1_exists == TRUE]
  genome <- rbindlist(lapply(seq_len(nrow(files)), function(j) {
    st <- stage1_taus(files$stage1[j])
    x <- fread(cmd = sprintf("zcat %s", shQuote(files$stage1[j])),
               select = c("snp_id", "chr", "bp", "effect_allele", "other_allele", "d_resid", st$beta_cols),
               showProgress = FALSE)
    Sigma0 <- read_sigma0(files$sigma0[j], st$taus)
    delta_mu <- as.numeric(as.matrix(x[, st$beta_cols, with = FALSE]) %*% cw)
    se <- sqrt(as.numeric(t(cw) %*% Sigma0 %*% cw) / x$d_resid)
    data.table(snp_id = x$snp_id, chr = x$chr, bp = x$bp, delta_mu = delta_mu, delta_mu_se = se,
               z = delta_mu / se)
  }))
  genome[, p := 2 * pnorm(-abs(z))]
  genome <- genome[is.finite(p)]
  out <- list(trait = tr, n_variants = nrow(genome), n_p = sum(genome$p < 5e-8))
  for (rule in c("A", "B")) {
    set <- paste0(tr, "_", rule, "_deltamu")
    ids <- run_clump(genome[, .(SNP = snp_id, CHR = chr, BP = bp, P = p, LOG10P = -log10(p))], set, rule)
    leads <- genome[snp_id %in% ids, .(snp_id, chr, bp, log10p = -log10(p))]
    leads[, `:=`(lead_set = set, trait = tr, rule = paste0(rule, "_deltamu"),
                 p_source = "Stage 1 delta_mu, pooled Sigma0")]
    out[[rule]] <- leads
  }
  out
}
traits <- stage1_traits(cfg)
survey <- Filter(stage1_complete, traits[group == "survey"]$trait)
deltamu_sets <- parallel::mclapply(c(survey, "weight"), deltamu_leads, mc.cores = cores, mc.preschedule = FALSE)
earlier_survey_leads <- file.path(dirname(cfg$cohorts$survey$root), "leads")
for (x in deltamu_sets) {
  if (inherits(x, "try-error")) stop(x)
  for (rule in c("A", "B")) {
    extra <- list(n_variants = x$n_variants)
    # Rule B is the earlier trait-survey rule: check it returns the earlier survey leads.
    reference <- file.path(earlier_survey_leads, paste0(x$trait, "_leads.csv"))
    if (rule == "B" && file.exists(reference)) {
      o <- fread(reference)
      extra <- c(extra, list(reproduces = setequal(x$B$snp_id, o$snp_id), n_reference = nrow(o),
                             n_only_new = length(setdiff(x$B$snp_id, o$snp_id)),
                             n_only_reference = length(setdiff(o$snp_id, x$B$snp_id)), reference = reference))
    }
    # The primary weight set is the REGENIE one; its delta_mu sets are sensitivity sets.
    add_set(list(leads = x[[rule]], n_candidates = x$n_p), paste0(x$trait, "_", rule, "_deltamu"),
            x$trait, paste0(rule, "_deltamu"), rule == "A" && x$trait != "weight", "Stage 1 delta_mu", extra)
  }
}

## 6. Sex-stratified testosterone leads (read from the REGENIE lead table)
sexstrat_file <- file.path(cfg$leads$testosterone$source_dir, "lead_snps.tsv")
sexstrat <- fread(sexstrat_file)[source == "regenie" & phenotype == "sex_strat"]
for (pm in unique(sexstrat$param)) {
  base <- sub("_(female|male)$", "", pm); sex <- sub(".*_", "", pm)
  set <- paste0(base, "_", sex, "_sexstrat_B")
  leads <- data.table(snp_id = sexstrat[param == pm]$lead_snp, chr = NA_integer_, bp = NA_integer_,
                      log10p = NA_real_, lead_set = set, trait = base, rule = "sexstrat_B",
                      p_source = paste("sex-stratified REGENIE", sex, "(lead_snps.tsv)"))
  add_set(list(leads = leads, n_candidates = NA), set, base, "sexstrat_B", FALSE, sexstrat_file)
}

## 7. Expected antibody loci among the JCV/MCV leads
# FUT2 (rs681343, rs676388) and HLA (rs9271055, rs3130287): exact lead, or the nearest lead
# within 500 kb (positions from the reference panel).
expected <- data.table(locus = c("FUT2", "FUT2", "HLA", "HLA"),
                       snp_id = c("rs681343", "rs676388", "rs9271055", "rs3130287"))
bim <- fread(paste0(ld_reference, ".bim"), header = FALSE, select = c(1, 2, 4), col.names = c("chr", "snp_id", "bp"))
expected <- merge(expected, bim[snp_id %in% expected$snp_id], by = "snp_id", all.x = TRUE, sort = FALSE)
antibody_loci <- rbindlist(lapply(c("jcv_vp1_log_all", "mcv_vp1_log_all"), function(tr) {
  L <- sets[[paste0(tr, "_A_regenie")]]
  rbindlist(lapply(seq_len(nrow(expected)), function(i) {
    e <- expected[i]
    near <- L[chr == e$chr & abs(bp - e$bp) <= 500000][order(abs(bp - e$bp))]
    data.table(trait = tr, locus = e$locus, expected_snp = e$snp_id, expected_chr = e$chr,
               expected_bp_1kg = e$bp, is_lead = e$snp_id %in% L$snp_id,
               nearest_lead = if (nrow(near)) near$snp_id[1] else NA_character_,
               nearest_lead_bp = if (nrow(near)) near$bp[1] else NA_integer_,
               nearest_lead_log10p = if (nrow(near)) near$log10p[1] else NA_real_,
               n_leads_within_500kb = nrow(near))
  }))
}))

## 8. Write
leads <- rbindlist(sets, fill = TRUE)
setcolorder(leads, c("lead_set", "trait", "rule", "primary", "snp_id", "chr", "bp", "log10p", "p_source"))
lead_sets <- rbindlist(summary_rows, fill = TRUE)
print(lead_sets[, .(lead_set, n_candidates, n_leads, primary)])
write_result(leads, "01_data/lead_variants", "leads")
write_result(lead_sets, "01_data/lead_variants", "lead_sets")
write_result(antibody_loci, "01_data/lead_variants", "antibody_expected_loci")
stopifnot("the clumping rule must reproduce the 639 BMI leads" = bmi_reproduced)
