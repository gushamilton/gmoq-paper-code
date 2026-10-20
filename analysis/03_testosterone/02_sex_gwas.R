#!/usr/bin/env Rscript
# GWAS of chromosomal sex in the testosterone Stage 1 sample
#
# What: runs a GWAS of chromosomal sex (male = 1) by logistic regression in REGENIE on exactly
#   the 399,710 participants of the testosterone Stage 1 GWAS (199,572 female, 200,138 male).
#   It is the reference for the testosterone membership parameter: -gamma, membership of the
#   upper (male) component, should track the variants whose allele frequencies differ by sex in
#   this sample. The options are those of the full-cohort sex GWAS (step 2 only, no LOCO
#   predictions, age and ten principal components), with the sample restricted by --keep.
#   The result is converted to LD score regression input and munged.
# Paper: Fig. 2a,b; Supplementary Table S3; Methods
# Input: UK Biobank sex phenotype, covariates without sex, keep list, LOCO predictions and
#   imputed genotypes (BGEN) with their sample file; HapMap3 variant lists
# Output: the sex reference files named in config.yaml (p08.sex_reference: LDSC-format
#   summary statistics; p08.sex_reference_munged); participant lists stay in the private
#   working directory
# Runs: on a cluster with UK Biobank data. REGENIE runs per chromosome (8 threads, about 1
#   hour each); the chromosomes are run in parallel on the cluster, here in a loop.

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
p8 <- cfg$p08
work <- private_dir(cfg, "03_testosterone", "sex_gwas")
sex_pheno_file <- file.path(dirname(p8$phenotype_file), "phenotypes_sex_paper.txt")   # FID IID sex

## 1. The testosterone Stage 1 sample ---------------------------------------------------
s <- stage1_analysed_sample(p8$phenotype_file, "testosterone", p8$covar_file, p8$keep_file,
                            p8$loco_list, p8$bgen_sample)
sx <- fread(sex_pheno_file, colClasses = list(character = c("FID", "IID")))
s <- merge(s, sx[!is.na(sex), .(FID, IID, sex)], by = c("FID", "IID"))
if (nrow(s) != cfg$cohorts$testosterone$n || sum(s$sex == 0) != 199572L || sum(s$sex == 1) != 200138L)
  stop("the keep list does not reproduce the Stage 1 sample: ", nrow(s), " (", sum(s$sex == 0), " F / ",
       sum(s$sex == 1), " M)")
keep_file <- file.path(work, "keep_testosterone_stage1.txt")
fwrite(s[, .(FID, IID)], keep_file, sep = " ", col.names = FALSE)

## 2. REGENIE step 2 per chromosome -----------------------------------------------------
for (chr in 1:22) {
  cc <- sprintf("%02d", chr)
  out <- file.path(work, paste0("sex_chr", cc))
  if (file.exists(paste0(out, "_sex.regenie"))) next
  run_shell(paste(
    "module load", cfg$p09b$regenie_module, "&& regenie --step 2",
    "--bgen", file.path(cfg$p09b$bgen_dir, paste0("data.chr", cc, ".bgen")), "--sample", p8$bgen_sample,
    "--phenoFile", sex_pheno_file, "--phenoCol sex", "--covarFile", p8$covar_file,
    "--extract", gsub("NN", cc, cfg$p09b$extract_template), "--bsize 1000 --threads 8",
    "--out", out, "--bt --ignore-pred --keep", keep_file))
}

## 3. LD score regression input --------------------------------------------------------
reg <- rbindlist(lapply(file.path(work, sprintf("sex_chr%02d_sex.regenie", 1:22)), fread))
if (!all(reg$N == cfg$cohorts$testosterone$n)) stop("REGENIE did not analyse the whole sample")
cap_p <- function(p) { p[!is.finite(p)] <- NA_real_; p[p <= 0] <- 1e-300; p[p > 1] <- NA_real_; p }
raw <- reg[, .(SNP = ID, CHR = as.integer(CHROM), BP = as.integer(GENPOS), A1 = ALLELE1, A2 = ALLELE0,
               N = as.integer(N), P = cap_p(10^(-LOG10P)), BETA = BETA, SE = SE)]
if (file.exists(p8$sex_reference)) stop("refusing to overwrite ", p8$sex_reference)
dir.create(dirname(p8$sex_reference), recursive = TRUE, showWarnings = FALSE)
fwrite(raw, p8$sex_reference, sep = "\t")
dir.create(dirname(p8$sex_reference_munged), recursive = TRUE, showWarnings = FALSE)
ldsc_munge(cfg, p8$sex_reference, sub("\\.sumstats\\.gz$", "", p8$sex_reference_munged))
