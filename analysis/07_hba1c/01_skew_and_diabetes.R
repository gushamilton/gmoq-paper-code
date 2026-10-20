#!/usr/bin/env Rscript
# Are the log-HbA1c leads with skew beyond shift + spread diabetes loci?
#
# What: Among the 472 log-HbA1c lead variants, 69 have skew beyond the shift + spread frame
#   (third L-moment of the frame residual, two-sided Bonferroni within trait;
#   06_traits/01_moments_and_principal_components.R). Here each lead is tested for association
#   with diabetes (logistic regression on dosage with age, sex and ten genetic PCs, in the same
#   participants) and with HbA1c among participants without diabetes, and the same SNPs are read
#   from the no-diabetes Stage 1 GWAS (converted exactly to the log scale) to see whether the skew
#   survives there. Skew-beyond leads are compared with the rest (enrichment tests). Finally the
#   six Fig. 5c leads are chosen by a rule fixed in advance: the three diabetes-associated leads
#   with the largest positive skew beyond the frame (one per 1 Mb), each matched on mean-effect z
#   to a lead without skew.
#   Diabetes is defined without HbA1c: doctor-diagnosed diabetes (UK Biobank field 2443) other
#   than gestational only (field 4041), or insulin use (fields 6177/6153), the definition used to
#   exclude participants from the no-diabetes Stage 1 run.
# Paper: Fig. 5c (lead selection, profiles); Supplementary Fig. S11; Results (HbA1c paragraph);
#   Methods "Moment decomposition and principal components" (exemplar rule)
# Input: UK Biobank HbA1c phenotypes, keep lists (all participants; no diabetes), diabetes case
#   list, Stage 1 covariates, imputed BGEN files; the no-diabetes HbA1c Stage 1 GWAS; aggregate
#   results/01_data/{lead_profiles,baseline_quantiles} and results/06_traits/moments/lead_moments.csv
# Output: results/07_hba1c/skew_diabetes/{per_lead,tests,profiles,examples,sample_counts}.csv
# Runs: on a cluster with UK Biobank data (needs bgenix and plink2; about two hours, mostly
#   genotype extraction over the 22 chromosomes)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
work <- file.path(cfg$bp$private_root, "07_hba1c", "skew_diabetes")   # participant-level files stay here
for (d in c("keep", "pheno", "snps", "geno", "glm", "ld")) dir.create(file.path(work, d), recursive = TRUE, showWarnings = FALSE)
u <- tau_grid(cfg); n_tau <- length(u); cw <- trap_weights(u)
n_leads_all <- 472L
z_bonf <- qnorm(1 - 0.025 / n_leads_all)                  # two-sided Bonferroni, 0.05 / 472
p_bonf <- 0.05 / n_leads_all
SET_ALL <- "log_hba1c_A_deltamu"; SET_ND <- "hba1c_nodiab_A_deltamu"
window_bp <- 500000L
lead_profiles <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))

## 1. Analysis samples (participant level, private)
# Rebuild the two Stage 1 samples: keep list -> non-missing phenotype -> complete covariates ->
# in the genotype sample file. They must reproduce the Stage 1 N exactly.
read_keep <- function(f) unique(fread(f, header = FALSE, col.names = c("FID", "IID"), colClasses = "character"))
pheno <- fread(cfg$hba1c$phenotypes, select = c("FID", "IID", "hba1c", "log_hba1c"),
               colClasses = list(character = c("FID", "IID")))
covar <- fread(cfg$prs$covar, colClasses = list(character = c("FID", "IID")))
covar_complete <- covar[complete.cases(covar), .(FID, IID)]
bgen_sample <- fread(cfg$ukb$bgen_sample, skip = 2, header = FALSE, select = 1:2, col.names = c("FID", "IID"),
                     colClasses = "character")
diabetes_ids <- unique(fread(cfg$hba1c$diabetes_cases, header = FALSE, colClasses = "character")[[1]])
funnel <- function(k, col, label) {
  f <- data.table(sample = label, stage = "keep", n = nrow(k))
  k <- merge(k, pheno[!is.na(get(col)), .(FID, IID)], by = c("FID", "IID"))
  f <- rbind(f, data.table(sample = label, stage = paste0("non_missing_", col), n = nrow(k)))
  k <- merge(k, covar_complete, by = c("FID", "IID"))
  f <- rbind(f, data.table(sample = label, stage = "complete_covariates", n = nrow(k)))
  k <- merge(k, bgen_sample, by = c("FID", "IID"))
  f <- rbind(f, data.table(sample = label, stage = "in_bgen_sample", n = nrow(k)))
  list(k = k, f = f)
}
all_p <- funnel(read_keep(cfg$hba1c$keep_all), "log_hba1c", "log_hba1c")
no_diab <- funnel(read_keep(cfg$hba1c$keep_no_diabetes), "hba1c", "hba1c_nodiab")
all_p$k[, diab := as.integer(IID %in% diabetes_ids)]
sample_counts <- rbind(all_p$f, no_diab$f,
                       data.table(sample = "log_hba1c", stage = c("diabetes_cases", "diabetes_controls"),
                                  n = c(sum(all_p$k$diab), sum(all_p$k$diab == 0L))))
stopifnot(nrow(all_p$k) == cfg$cohorts$survey$traits$log_hba1c$n,
          nrow(no_diab$k) == cfg$cohorts$survey$traits$hba1c_nodiab$n)
# The no-diabetes sample is exactly the controls of the full sample.
controls <- all_p$k[diab == 0L, .(FID, IID)]
stopifnot(nrow(controls) == nrow(no_diab$k), nrow(fsetdiff(controls, no_diab$k[, .(FID, IID)])) == 0L)
keep_all <- file.path(work, "keep", "keep_log_hba1c.txt"); keep_union <- file.path(work, "keep", "keep_union.txt")
fwrite(all_p$k[, .(FID, IID)], keep_all, sep = " ", col.names = FALSE)
fwrite(unique(rbind(all_p$k[, .(FID, IID)], no_diab$k[, .(FID, IID)])), keep_union, sep = " ", col.names = FALSE)
P <- merge(all_p$k[, .(FID, IID, diab)], pheno, by = c("FID", "IID"), all.x = TRUE)
P[, in_nd := paste(FID, IID) %in% paste(no_diab$k$FID, no_diab$k$IID)]
pheno_file <- file.path(work, "pheno", "pheno.txt")       # diab coded 1/2 for plink2
fwrite(P[, .(FID, IID, diab = diab + 1L, hba1c_nd = fifelse(in_nd, hba1c, NA_real_),
             log_hba1c_nd = fifelse(in_nd, log_hba1c, NA_real_), log_hba1c_all = log_hba1c)],
       pheno_file, sep = "\t", na = "NA", quote = FALSE)
rm(P, pheno, covar); invisible(gc())

## 2. Lead genotypes, regressions and LD (plink2), chromosome by chromosome
lead_snps <- unique(lead_profiles[lead_set %in% c(SET_ALL, SET_ND), .(snp_id, chr)])
plink <- function(...) { status <- system2(cfg$ukb$plink2, c(...), stdout = FALSE); status }
covar_args <- c("--covar", cfg$prs$covar, "--covar-name", "age", "sex", "PC1-PC10", "--covar-variance-standardize")
for (cc in sort(unique(lead_snps$chr))) {
  tag <- sprintf("chr%02d", cc)
  bgen <- file.path(work, "geno", paste0(tag, ".bgen")); pfile <- file.path(work, "geno", tag)
  system2(cfg$ukb$bgenix, c("-g", file.path(cfg$ukb$bgen_dir, sprintf("data.chr%02d.bgen", cc)),
                            "-incl-rsids", lead_snps[chr == cc]$snp_id), stdout = bgen)
  plink("--bgen", bgen, "ref-first", "--sample", cfg$ukb$bgen_sample, "--keep", keep_union,
        "--make-pgen", "--out", pfile)
  unlink(bgen)
  # diabetes ~ dosage + covariates (logistic, Firth fallback), all participants
  plink("--pfile", pfile, "--keep", keep_all, "--pheno", pheno_file, "--pheno-name", "diab", covar_args,
        "--glm", "hide-covar", "firth-fallback", "cols=+a1freq", "--out", file.path(work, "glm", paste0(tag, "_diab")))
  # HbA1c (raw and log) without diabetes, and log HbA1c in everyone (linear)
  plink("--pfile", pfile, "--keep", keep_union, "--pheno", pheno_file,
        "--pheno-name", "hba1c_nd", "log_hba1c_nd", "log_hba1c_all", covar_args,
        "--glm", "hide-covar", "cols=+a1freq", "--out", file.path(work, "glm", paste0(tag, "_lin")))
  # unphased r2 between leads within 500 kb (only used to find nearby no-diabetes leads)
  plink("--pfile", pfile, "--keep", keep_all, "--r2-unphased", "--ld-window-kb", 500, "--ld-window", 999999,
        "--ld-window-r2", 0, "--out", file.path(work, "ld", tag))
}

## 3. The 472 log-HbA1c leads and their moment functionals
moments <- as.data.table(read_result("06_traits/moments", "lead_moments"))
baseline <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))
baseline_quantiles <- function(trait_name) {
  x <- baseline[baseline$trait == trait_name]
  as.numeric(x$q[match(round(u, 3), round(as.numeric(x$tau), 3))])
}
q_all <- baseline_quantiles("log_hba1c")                    # log scale
q_nd_raw <- baseline_quantiles("hba1c_nodiab"); q_nd <- log(q_nd_raw)
leads <- lead_set_matrix(lead_profiles, SET_ALL, u)          # aligned to the HbA1c-raising allele
D <- merge(as.data.table(leads$meta)[, .(snp_id, chr, bp, effect_allele, other_allele, d_resid)],
           moments[arm == "log_hba1c", .(snp_id, delta_mu, delta_mu_se, delta_mu_z, delta_lambda2_z,
                                          skew_beyond = skew_beyond_frame, skew_beyond_se = skew_beyond_frame_se,
                                          skew_beyond_z = skew_beyond_frame_z, lof_p_davies)],
           by = "snp_id", sort = FALSE)
stopifnot(nrow(D) == n_leads_all)
D[, skew_hit := abs(skew_beyond_z) > z_bonf]

## 4. The same SNPs in the no-diabetes Stage 1, converted exactly to log
# beta_log = beta / Q and Sigma_log = D Sigma0 D with D = diag(1 / Q(u)); functionals use the
# no-diabetes baseline on the log scale.
paths_nd <- stage1_chr_paths(cfg, "hba1c_nodiab")[stage1_exists == TRUE]
cols_nd <- stage1_taus(paths_nd$stage1[1]); ci <- match(round(u, 3), round(cols_nd$taus, 3))
nd <- rbindlist(lapply(sort(unique(D$chr)), function(cc) {
  x <- read_stage1(paths_nd[chr == cc]$stage1[1], snps = D[chr == cc]$snp_id)
  x[, snp_id := as.character(snp_id)]; x
}), fill = TRUE)
nd <- nd[match(D$snp_id, nd$snp_id)]
B_nd <- as.matrix(nd[, cols_nd$beta_cols[ci], with = FALSE]); storage.mode(B_nd) <- "double"
allele_sign <- ifelse(nd$effect_allele == D$effect_allele & nd$other_allele == D$other_allele, 1,
                      ifelse(nd$effect_allele == D$other_allele & nd$other_allele == D$effect_allele, -1, NA))
B_nd <- B_nd * allele_sign
to_log <- 1 / q_nd_raw
B_nd_log <- sweep(B_nd, 2, to_log, `*`)
sigma0_nd <- p05_sigma0("hba1c_nodiab", u); sigma0_all <- p05_sigma0("log_hba1c", u)
W_nd <- functional_weights(u, q_nd)[c("delta_mu", "delta_var", "skew_beyond_frame")]
nd_fun <- matrix(NA_real_, nrow(D), 5, dimnames = list(NULL, c("delta_mu", "delta_mu_se", "delta_var_z",
                                                             "skew_beyond", "skew_beyond_se")))
for (k in unique(as.character(D$chr))) {
  i <- which(as.character(D$chr) == k & !is.na(allele_sign))
  if (!length(i)) next
  S <- change_scale(matrix(0, 1, n_tau), sigma0_nd[[k]], to_log)$Sigma0
  fe <- functional_estimate(B_nd_log[i, , drop = FALSE], S, nd$d_resid[i], W_nd)
  nd_fun[i, ] <- cbind(fe$estimate[, "delta_mu"], fe$se[, "delta_mu"],
                       fe$estimate[, "delta_var"] / fe$se[, "delta_var"],
                       fe$estimate[, "skew_beyond_frame"], fe$se[, "skew_beyond_frame"])
}
D[, `:=`(nd_stage1_found = !is.na(allele_sign), nd_d_resid = nd$d_resid,
         nd_delta_mu = nd_fun[, "delta_mu"], nd_delta_mu_se = nd_fun[, "delta_mu_se"],
         nd_delta_mu_z = nd_fun[, "delta_mu"] / nd_fun[, "delta_mu_se"], nd_delta_var_z = nd_fun[, "delta_var_z"],
         nd_skew_beyond = nd_fun[, "skew_beyond"], nd_skew_beyond_se = nd_fun[, "skew_beyond_se"],
         nd_skew_beyond_z = nd_fun[, "skew_beyond"] / nd_fun[, "skew_beyond_se"])]
D[, nd_skew_hit_472 := abs(nd_skew_beyond_z) > z_bonf]

## 5. Regression results, aligned to the HbA1c-raising allele
read_glm <- function(pattern) {
  x <- rbindlist(lapply(list.files(file.path(work, "glm"), pattern = pattern, full.names = TRUE), fread), fill = TRUE)
  setnames(x, names(x), sub("^#", "", names(x)))
  x
}
align <- function(x, prefix) {
  x <- x[TEST == "ADD"]
  b <- if ("BETA" %in% names(x)) x$BETA else log(x$OR)
  se <- if ("SE" %in% names(x)) x$SE else x[["LOG(OR)_SE"]]
  z_col <- intersect(c("T_STAT", "Z_STAT"), names(x))[1]
  other <- ifelse(x$A1 == x$ALT, x$REF, x$ALT)
  y <- data.table(snp_id = x$ID, A1 = x$A1, O = other, b = b, se = se, z = x[[z_col]], p = x$P, n = x$OBS_CT,
                  a1freq = x$A1_FREQ, err = if ("ERRCODE" %in% names(x)) x$ERRCODE else NA_character_,
                  firth = if ("FIRTH?" %in% names(x)) x[["FIRTH?"]] else NA_character_)
  m <- y[match(D$snp_id, y$snp_id)]
  s <- ifelse(m$A1 == D$effect_allele & m$O == D$other_allele, 1,
              ifelse(m$A1 == D$other_allele & m$O == D$effect_allele, -1, NA))
  out <- data.table(b = m$b * s, se = m$se, z = m$z * s, p = m$p, n = m$n,
                    eaf = ifelse(s == 1, m$a1freq, 1 - m$a1freq), err = m$err, firth = m$firth)
  setnames(out, paste0(prefix, "_", names(out)))
  out
}
G_diab <- align(read_glm("_diab\\.diab\\.glm\\.logistic"), "diab")
G_nd <- align(read_glm("_lin\\.log_hba1c_nd\\.glm\\.linear$"), "lin_nd_log")
G_nd_raw <- align(read_glm("_lin\\.hba1c_nd\\.glm\\.linear$"), "lin_nd_raw")
G_all <- align(read_glm("_lin\\.log_hba1c_all\\.glm\\.linear$"), "lin_all_log")
D <- cbind(D, G_diab[, .(diab_b, diab_se, diab_z, diab_p, diab_n, diab_eaf, diab_err, diab_firth)],
           G_nd[, .(lin_nd_log_b, lin_nd_log_se, lin_nd_log_z, lin_nd_log_p, lin_nd_log_n)],
           G_nd_raw[, .(lin_nd_raw_b, lin_nd_raw_se, lin_nd_raw_z, lin_nd_raw_p)],
           G_all[, .(lin_all_log_b, lin_all_log_se, lin_all_log_z, lin_all_log_p, lin_all_log_n)])
D[, ratio_diab_to_nd := diab_b / lin_nd_log_b]   # log OR(diabetes) per unit log-HbA1c effect without diabetes

## 6. Nearby no-diabetes leads: nearest within 500 kb, and the one in highest LD
nd_leads <- merge(moments[arm == "hba1c_nodiab_log", .(snp_id, chr, bp, skew_log_z = skew_beyond_frame_z)],
                  moments[arm == "hba1c_nodiab", .(snp_id, skew_raw_z = skew_beyond_frame_z)], by = "snp_id")
ld_files <- list.files(file.path(work, "ld"), pattern = "\\.vcor$|\\.ld$", full.names = TRUE)
LD <- if (length(ld_files)) {
  x <- rbindlist(lapply(ld_files, fread), fill = TRUE); setnames(x, names(x), sub("^#", "", names(x)))
  r2_col <- intersect(c("UNPHASED_R2", "R2"), names(x))[1]
  rbind(x[, .(a = ID_A, b = ID_B, r2 = get(r2_col))], x[, .(a = ID_B, b = ID_A, r2 = get(r2_col))])
} else data.table(a = character(), b = character(), r2 = numeric())
near <- rbindlist(lapply(seq_len(nrow(D)), function(i) {
  c0 <- nd_leads[chr == D$chr[i] & abs(bp - D$bp[i]) <= window_bp]
  if (!nrow(c0)) return(data.table(nd_lead_within_500kb = FALSE, n_nd_leads_500kb = 0L, nearest_nd_lead = NA_character_,
                                   nearest_nd_dist_bp = NA_integer_, nearest_nd_skew_log_z = NA_real_,
                                   nearest_nd_skew_raw_z = NA_real_, maxr2_nd_lead = NA_character_, maxr2 = NA_real_,
                                   maxr2_nd_skew_log_z = NA_real_, same_snp_is_nd_lead = FALSE))
  c0[, dist := abs(bp - D$bp[i])]
  j <- which.min(c0$dist)
  r <- LD[a == D$snp_id[i] & b %in% c0$snp_id]
  r2 <- ifelse(c0$snp_id == D$snp_id[i], 1, r$r2[match(c0$snp_id, r$b)])
  k <- if (all(is.na(r2))) NA_integer_ else which.max(r2)
  data.table(nd_lead_within_500kb = TRUE, n_nd_leads_500kb = nrow(c0), nearest_nd_lead = c0$snp_id[j],
             nearest_nd_dist_bp = as.integer(c0$dist[j]), nearest_nd_skew_log_z = c0$skew_log_z[j],
             nearest_nd_skew_raw_z = c0$skew_raw_z[j],
             maxr2_nd_lead = if (is.na(k)) NA_character_ else c0$snp_id[k], maxr2 = if (is.na(k)) NA_real_ else r2[k],
             maxr2_nd_skew_log_z = if (is.na(k)) NA_real_ else c0$skew_log_z[k],
             same_snp_is_nd_lead = D$snp_id[i] %in% c0$snp_id)
}))
D <- cbind(D, near)
D[, nd_lead_r2_gt_0.1 := !is.na(maxr2) & maxr2 > 0.1]

## 7. Nearest gene (hg19 UCSC known genes)
D[, `:=`(nearest_gene = NA_character_, gene_distance_bp = NA_integer_)]
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- suppressMessages(GenomicFeatures::genes(txdb))
genes <- genes[as.character(GenomicRanges::seqnames(genes)) %in% paste0("chr", 1:22)]
hits <- GenomicRanges::distanceToNearest(GenomicRanges::GRanges(paste0("chr", D$chr), IRanges::IRanges(D$bp, D$bp)), genes)
qi <- S4Vectors::queryHits(hits); si <- S4Vectors::subjectHits(hits)
symbol <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(S4Vectors::mcols(genes)$gene_id[si]),
                                                 column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
D$nearest_gene[qi] <- unname(symbol); D$gene_distance_bp[qi] <- as.integer(S4Vectors::mcols(hits)$distance)

## 8. Skew-beyond leads against the rest
h <- D$skew_hit; abs_diab_z <- abs(D$diab_z)
tests <- list()
add <- function(q, v, note = "") tests[[length(tests) + 1L]] <<- data.table(quantity = q, value = as.numeric(v), note = note)
add("n_leads", nrow(D)); add("n_skew_hits", sum(h), sprintf("|z| > %.3f (Bonferroni 0.05/472, two-sided)", z_bonf))
add("n_skew_hits_positive", sum(D$skew_beyond_z > z_bonf))
for (grp in c("hit", "nonhit")) {
  s <- if (grp == "hit") h else !h
  add(paste0("median_abs_z_diab_", grp), median(abs_diab_z[s], na.rm = TRUE))
  add(paste0("frac_diab_p_lt_5e-8_", grp), mean(D$diab_p[s] < 5e-8, na.rm = TRUE))
  add(paste0("n_diab_p_lt_5e-8_", grp), sum(D$diab_p[s] < 5e-8, na.rm = TRUE))
  add(paste0("frac_diab_p_lt_bonf472_", grp), mean(D$diab_p[s] < p_bonf, na.rm = TRUE))
  add(paste0("n_diab_p_lt_bonf472_", grp), sum(D$diab_p[s] < p_bonf, na.rm = TRUE))
  add(paste0("frac_diab_z_pos_", grp), mean(D$diab_z[s] > 0, na.rm = TRUE), "HbA1c-increasing allele raises diabetes risk")
  add(paste0("median_abs_delta_mu_z_", grp), median(abs(D$delta_mu_z[s])))
  add(paste0("median_ratio_diab_to_nd_", grp), median(D$ratio_diab_to_nd[s], na.rm = TRUE),
      "log OR(diabetes) / beta(log HbA1c, no diabetes), per allele")
  add(paste0("median_abs_nd_skew_z_", grp), median(abs(D$nd_skew_beyond_z[s]), na.rm = TRUE), "same SNP, no-diabetes Stage 1, log scale")
  add(paste0("n_nd_skew_hit_same_snp_", grp), sum(D$nd_skew_hit_472[s], na.rm = TRUE))
  add(paste0("median_abs_nd_delta_mu_z_", grp), median(abs(D$nd_delta_mu_z[s]), na.rm = TRUE))
  add(paste0("frac_nd_lead_within_500kb_", grp), mean(D$nd_lead_within_500kb[s]))
  add(paste0("frac_nd_lead_r2_gt_0.1_", grp), mean(D$nd_lead_r2_gt_0.1[s]))
}
add("wilcoxon_p_abs_z_diab_hit_vs_nonhit", wilcox.test(abs_diab_z[h], abs_diab_z[!h])$p.value)
ft <- fisher.test(table(h, D$diab_p < 5e-8)); add("fisher_or_diab_gws_hit", ft$estimate); add("fisher_p_diab_gws_hit", ft$p.value)
ft <- fisher.test(table(h, D$diab_p < p_bonf)); add("fisher_or_diab_bonf_hit", ft$estimate); add("fisher_p_diab_bonf_hit", ft$p.value)
ct <- cor.test(D$skew_beyond_z, D$diab_z, method = "spearman", exact = FALSE)
add("spearman_skew_z_vs_diab_z", ct$estimate); add("spearman_p_skew_z_vs_diab_z", ct$p.value)
add("pearson_skew_z_vs_diab_z", cor(D$skew_beyond_z, D$diab_z, use = "complete.obs"))
add("wilcoxon_p_ratio_hit_vs_nonhit", wilcox.test(D$ratio_diab_to_nd[h], D$ratio_diab_to_nd[!h])$p.value)
# Adjusted for the strength of the mean effect (power for both skew and diabetes grows with it).
cf <- summary(glm(h ~ log(abs(delta_mu_z)) + abs(diab_z), family = binomial, data = D))$coefficients
add("logit_hit_abs_diab_z_coef_adj_logmeanz", cf["abs(diab_z)", 1]); add("logit_hit_abs_diab_z_p_adj_logmeanz", cf["abs(diab_z)", 4])
lf <- lm(skew_beyond_z ~ delta_mu_z + diab_z, data = D); lc <- summary(lf)$coefficients
add("lm_skew_z_on_diab_z_adj_meanz_coef", lc["diab_z", 1]); add("lm_skew_z_on_diab_z_adj_meanz_p", lc["diab_z", 4])
add("lm_skew_z_partial_r2_diab_z", summary(lf)$r.squared - summary(lm(skew_beyond_z ~ delta_mu_z, data = D))$r.squared)
# Inverse-variance-weighted slope through the origin of diabetes log OR on the no-diabetes effect.
for (grp in c("hit", "nonhit")) {
  s <- if (grp == "hit") h else !h
  add(paste0("ivw_slope_logOR_per_log_hba1c_nd_", grp),
      coef(lm(diab_b ~ 0 + lin_nd_log_b, data = D[s], weights = 1 / diab_se^2))[1], "IVW through origin")
}
# Agreement checks: Stage 1 location against the linear regressions; the converted no-diabetes
# skew against the moment run of the no-diabetes leads (same SNPs).
add("cor_stage1_delta_mu_vs_plink_all_log_b", cor(D$delta_mu, D$lin_all_log_b, use = "complete.obs"), "Stage 1 delta_mu (x 0.96) vs OLS")
add("cor_nd_stage1_delta_mu_vs_plink_nd_log_b", cor(D$nd_delta_mu, D$lin_nd_log_b, use = "complete.obs"))
add("n_nd_stage1_found", sum(D$nd_stage1_found)); add("n_diab_found", sum(!is.na(D$diab_z)))
add("n_ld_pairs", nrow(LD) / 2)
ov <- merge(D[, .(snp_id, here = nd_skew_beyond_z)], moments[arm == "hba1c_nodiab_log", .(snp_id, mom = skew_beyond_frame_z)], by = "snp_id")
add("n_overlap_with_nd_leads", nrow(ov))
add("max_abs_diff_abs_nd_skew_z_vs_moments_overlap", if (nrow(ov)) max(abs(abs(ov$here) - abs(ov$mom))) else NA,
    "reproduction check of the change of scale")
add("diab_cases_log_hba1c_sample", sample_counts[stage == "diabetes_cases"]$n, "aggregate")
add("diab_controls_log_hba1c_sample", sample_counts[stage == "diabetes_controls"]$n, "aggregate")
tests <- rbindlist(tests)
print(tests)

## 9. Profiles of both analyses on the log scale, with pooled SE and the shift + spread fit
frame_fit <- function(B, q) { X <- cbind(1, q - sum(cw * q) / sum(cw)); B %*% (X %*% solve(crossprod(X), t(X))) }
se_all <- se_nd <- matrix(NA_real_, nrow(D), n_tau)
for (k in unique(as.character(D$chr))) {
  i <- which(as.character(D$chr) == k)
  se_all[i, ] <- sqrt(outer(1 / D$d_resid[i], diag(sigma0_all[[k]])))
  se_nd[i, ] <- sqrt(outer(1 / nd$d_resid[i], diag(change_scale(matrix(0, 1, n_tau), sigma0_nd[[k]], to_log)$Sigma0)))
}
long <- function(B, S, F, label) data.table(snp_id = rep(D$snp_id, n_tau), analysis = label,
                                            tau = rep(round(u, 3), each = nrow(D)), beta = as.vector(B),
                                            se_pooled = as.vector(S), frame_fit = as.vector(F))
profiles <- rbind(long(leads$B, se_all, frame_fit(leads$B, q_all), "all_participants_log"),
                  long(B_nd_log, se_nd, frame_fit(B_nd_log, q_nd), "no_diabetes_log"))

## 10. Fig. 5c leads (rule fixed before looking at the profiles)
candidates <- D[skew_hit & skew_beyond_z > 0 & diab_p < p_bonf][order(-skew_beyond_z)]
hits <- candidates[0]
for (i in seq_len(nrow(candidates))) {                     # one per 1 Mb
  if (nrow(hits) && any(hits$chr == candidates$chr[i] & abs(hits$bp - candidates$bp[i]) < 1e6)) next
  hits <- rbind(hits, candidates[i]); if (nrow(hits) == 3L) break
}
pool <- D[abs(skew_beyond_z) < 1.96]                       # no skew beyond the frame
matched <- pool[0]
for (i in seq_len(nrow(hits))) {                           # closest log|mean z|, without replacement
  pp <- pool[!(snp_id %in% matched$snp_id)]
  matched <- rbind(matched, pp[which.min(abs(log(abs(pp$delta_mu_z)) - log(abs(hits$delta_mu_z[i]))))])
}
examples <- rbind(data.table(role = "skew_hit", rank = seq_len(nrow(hits)), snp_id = hits$snp_id, matched_to = NA_character_),
                  data.table(role = "non_skew_matched", rank = seq_len(nrow(matched)), snp_id = matched$snp_id,
                             matched_to = hits$snp_id))
examples <- merge(examples, D[, .(snp_id, chr, bp, nearest_gene, gene_distance_bp, delta_mu_z, skew_beyond_z, diab_z,
                                  diab_p, nd_skew_beyond_z, lin_nd_log_z)], by = "snp_id", sort = FALSE)
examples[, rule := paste("skew_hit: top 3 positive Bonferroni skew-beyond hits (log HbA1c) with diabetes p < 0.05/472,",
                         "by skew-beyond z, one per 1 Mb; non_skew_matched: |skew-beyond z| < 1.96, closest log|mean z|",
                         "to each hit, without replacement")]
print(examples)

## 11. Write
out <- "07_hba1c/skew_diabetes"
write_result(D, out, "per_lead")
write_result(tests, out, "tests")
write_result(profiles, out, "profiles")
write_result(examples, out, "examples")
setnames(sample_counts, "n", "n_participants")
write_result(sample_counts, out, "sample_counts")
