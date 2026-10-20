#!/usr/bin/env bash
# LD score regression: SNP heritability and genetic correlations
#
# What: converts GWAS summary statistics (REGENIE step 2, or genome-wide Stage 2 parameter
#   estimates) to a common format, munges them onto the HapMap3 alleles and runs LDSC for
#   SNP heritability and for genetic correlations between pairs. Used for the SNP heritability
#   of log BMI and for the genetic correlations of Stage 2 parameters with conventional GWAS
#   (testosterone component effects with sex-stratified testosterone GWAS, the testosterone
#   membership parameter with the GWAS of sex, and the genome-wide BMI models).
# Paper: Fig. 2 (genetic correlations with the sex-stratified reference); Fig. 5d (SNP heritability
#   of log BMI); Methods (genome-wide BMI models)
# Input: REGENIE .regenie files (02_regenie.sh) or Stage 2 per-variant tables
#   (columns snp_id chr bp effect_allele other_allele <parameter> <parameter>_se)
# Output: <out_dir>/raw/<id>.tsv, munged/<id>.sumstats.gz, h2/<id>.log, rg/<id1>__<id2>.log
# Runs: on a cluster. LDSC 2.0.0 (munge_sumstats.py, ldsc.py); 1000 Genomes European LD scores
#   (eur_w_ld_chr) and the HapMap3 SNP list (w_hm3.snplist).
set -euo pipefail

# ---------------------------------------------------------------- user settings --
OUT=/path/to/ldsc
MERGE_ALLELES=/path/to/ldsc/w_hm3.snplist.gz
LD_DIR=/path/to/ldsc/eur_w_ld_chr/                      # used for both --ref-ld-chr and --w-ld-chr
REGENIE_STEP2=/path/to/gwas_inputs/regenie_step2
mkdir -p "$OUT"/{raw,munged,h2,rg}

## 1. Summary statistics to LDSC input format -----------------------------------------------
# Columns SNP CHR BP A1 A2 N P BETA SE, A1 = effect allele.
#   REGENIE: A1 = ALLELE1, A2 = ALLELE0, N per variant, P = 10^-LOG10P.
#   Stage 2: N = the trait's analysed sample size, P from the Wald chi-square (estimate / SE)^2, 1 df.
# P values that underflow to 0 are set to 1e-300; values outside (0, 1] become missing.
regenie_to_raw() {   # id file1.regenie [file2.regenie ...]
  local id=$1; shift
  Rscript -e '
    suppressPackageStartupMessages(library(data.table))
    a <- commandArgs(TRUE); out <- a[1]
    d <- rbindlist(lapply(a[-1], fread), use.names = TRUE, fill = FALSE)
    p <- 10^(-d$LOG10P); p[!is.finite(p)] <- NA; p[p <= 0] <- 1e-300; p[p > 1] <- NA
    fwrite(d[, .(SNP = ID, CHR = as.integer(CHROM), BP = as.integer(GENPOS), A1 = ALLELE1, A2 = ALLELE0,
                 N = as.integer(N), P = p, BETA, SE)], out, sep = "\t")' "$OUT/raw/$id.tsv" "$@"
}
stage2_to_raw() {    # id stage2.tsv.gz parameter N
  Rscript -e '
    suppressPackageStartupMessages(library(data.table))
    a <- commandArgs(TRUE); prm <- a[3]; se <- paste0(prm, "_se")
    d <- fread(cmd = paste("zcat", shQuote(a[2])), select = c("snp_id", "chr", "bp", "effect_allele", "other_allele", prm, se))
    p <- pchisq((d[[prm]] / d[[se]])^2, df = 1, lower.tail = FALSE)
    p[!is.finite(p)] <- NA; p[p <= 0] <- 1e-300; p[p > 1] <- NA
    fwrite(data.table(SNP = d$snp_id, CHR = as.integer(d$chr), BP = as.integer(d$bp), A1 = d$effect_allele,
                      A2 = d$other_allele, N = as.integer(a[4]), P = p, BETA = d[[prm]], SE = d[[se]]),
           a[1], sep = "\t")' "$OUT/raw/$1.tsv" "$2" "$3" "$4"
}

## 2. Munge onto the HapMap3 alleles ---------------------------------------------------------------
# The effect direction is read from BETA (null value 0); SE is ignored by LDSC.
munge() {   # id
  munge_sumstats.py --sumstats "$OUT/raw/$1.tsv" --out "$OUT/munged/$1" --merge-alleles "$MERGE_ALLELES" \
    --snp SNP --a1 A1 --a2 A2 --p P --N-col N --signed-sumstats BETA,0 --ignore SE
}

## 3. SNP heritability and genetic correlation ------------------------------------------------------
h2() { ldsc.py --h2 "$OUT/munged/$1.sumstats.gz" --ref-ld-chr "$LD_DIR" --w-ld-chr "$LD_DIR" --out "$OUT/h2/$1"; }
rg() { ldsc.py --rg "$OUT/munged/$1.sumstats.gz,$OUT/munged/$2.sumstats.gz" --ref-ld-chr "$LD_DIR" --w-ld-chr "$LD_DIR" \
         --out "$OUT/rg/$1__$2"; }

## 4. Runs ---------------------------------------------------------------------------------------------
# REGENIE GWAS of BMI and log BMI: SNP heritability (log BMI is the Fig. 5d value).
for t in bmi log_bmi; do
  regenie_to_raw "regenie__${t}_additive" $(ls "$REGENIE_STEP2/$t/${t}_additive_chr"*"_$t.regenie" | sort)
  munge "regenie__${t}_additive"; h2 "regenie__${t}_additive"
done
# Sex-stratified testosterone GWAS and the GWAS of sex on the testosterone Stage 1 sample.
for t in testosterone testosterone_shrink_0_2 testosterone_shrink_0_5 testosterone_shrink_0_7; do
  for s in female male; do
    regenie_to_raw "regenie__${t}_$s" $(ls "$REGENIE_STEP2/sex_strat/${t}_${s}_chr"*.regenie | sort)
    munge "regenie__${t}_$s"
  done
done
regenie_to_raw regenie__sex_matched $(ls "$REGENIE_STEP2/sex_gwas_matched/sex_chr"*_sex.regenie | sort)
munge regenie__sex_matched

# Stage 2 parameters are converted with stage2_to_raw, munged, and correlated with the matching
# conventional GWAS by rg(), e.g. for each testosterone arm (Stage 2 tables from analysis/02_mixtures):
#   stage2_to_raw gmoq_testosterone_beta_left  <stage2 table> beta_left  <N>; munge ...; rg gmoq_testosterone_beta_left  regenie__testosterone_female
#   stage2_to_raw gmoq_testosterone_beta_right <stage2 table> beta_right <N>; munge ...; rg gmoq_testosterone_beta_right regenie__testosterone_male
# The membership parameter gamma is correlated with regenie__sex_matched after flipping its sign
# (gamma moves mass to the lower, female, component; the sex GWAS codes male = 1). The analysis
# scripts that make these calls are in analysis/03_testosterone and analysis/05_bmi.
