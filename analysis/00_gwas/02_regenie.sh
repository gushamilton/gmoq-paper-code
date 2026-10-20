#!/usr/bin/env bash
# Conventional mean GWAS with REGENIE (step 1 LOCO predictions and step 2 association tests)
#
# What: fits REGENIE step 1 on the directly genotyped array to obtain leave-one-chromosome-out
#   (LOCO) polygenic predictions, which the Stage 1 quantile GWAS includes as a covariate, and
#   runs REGENIE step 2 on the imputed HapMap3 variants as the conventional mean-effect GWAS:
#   BMI and log BMI, height-adjusted weight, sex-stratified testosterone, JC and Merkel cell
#   polyomavirus VP1 antibodies (all-sample levels, serostatus, seropositive-only levels) and
#   the GWAS of chromosomal sex on the testosterone Stage 1 sample. Step 1 is also run for the
#   three survey traits whose Stage 1 used LOCO (height, log CRP, log triglycerides).
# Paper: Methods (comparator GWAS); references for Fig. 2 and Fig. 3; LOCO inputs to every Stage 1 run
# Input: phenotype, covariate and keep files from 01_build_phenotypes.R; UK Biobank directly
#   genotyped array (PLINK bed) and imputed BGEN; HapMap3 variant lists per chromosome
# Output: <out_dir>/regenie_step1/<trait>/ (LOCO predictions) and <out_dir>/regenie_step2/<trait>/
#   (per-chromosome .regenie files)
# Runs: on a cluster with UK Biobank data. REGENIE v4.1.2, PLINK 2.0 (2.00a6.8LM). Chromosomes
#   were run in parallel on the cluster; here they are a loop.
set -euo pipefail

# ---------------------------------------------------------------- user settings --
IN=/path/to/gwas_inputs                         # out_dir of 01_build_phenotypes.R
OUT=/path/to/gwas_inputs                        # REGENIE output root
ARRAY_BED=/path/to/ukb_genotyped/chr1-22_merged # directly genotyped array, autosomes merged
BGEN_DIR=/path/to/ukb_imputed                   # data.chrNN.bgen
SAMPLE=/path/to/ukb_imputed/data.chr1-22.sample
HM3=/path/to/hapmap3/w_hm3_chrNN_only.txt       # HapMap3 rsIDs on chromosome NN
THREADS=8
KEEP=$IN/keep
hm3() { echo "${HM3//NN/$1}"; }
chr2() { printf "%02d" "$1"; }

# ------------------------------------------------------------------ helpers --
# Step 1: variant QC on the analysed sample (MAF >= 1%, MAC >= MIN_MAC, call rate >= 90%,
# HWE p >= 1e-15, sample call rate >= 90%), then the REGENIE ridge fit (five-fold
# cross-validation by default, block size 1000, low-memory mode).
step1() {   # out_prefix pheno_file keep covar_file min_mac  [extra REGENIE arguments...]
  local prefix=$1 pheno=$2 keep=$3 covar=$4 mac=$5; shift 5
  mkdir -p "$(dirname "$prefix")"
  plink2 --bfile "$ARRAY_BED" --keep "$keep" \
    --maf 0.01 --mac "$mac" --geno 0.1 --hwe 1e-15 0 --mind 0.1 \
    --threads "$THREADS" --write-snplist --write-samples --no-id-header --out "${prefix}_qc_pass"
  regenie --step 1 --bed "$ARRAY_BED" --phenoFile "$pheno" --covarFile "$covar" --keep "$keep" \
    --extract "${prefix}_qc_pass.snplist" --bsize 1000 \
    --lowmem --lowmem-prefix "${prefix}_tmp" --threads "$THREADS" --out "${prefix}_step1" "$@"
}
# Step 2: association test on the imputed dosages, HapMap3 variants only.
step2() {   # out_prefix chr pheno_file covar_file  [extra REGENIE arguments...]
  local prefix=$1 c=$2 pheno=$3 covar=$4; shift 4
  mkdir -p "$(dirname "$prefix")"
  regenie --step 2 --bgen "$BGEN_DIR/data.chr$(chr2 "$c").bgen" --sample "$SAMPLE" \
    --phenoFile "$pheno" --covarFile "$covar" --extract "$(hm3 "$(chr2 "$c")")" \
    --threads "$THREADS" --out "$prefix" "$@"
}
S1=$OUT/regenie_step1; S2=$OUT/regenie_step2

## 1. Step 1 for the main traits -----------------------------------------------------------
# BMI and log BMI: covariates age, sex, PC1-PC10.
for t in bmi log_bmi; do
  step1 "$S1/$t/$t" "$IN/phenotypes_paper.txt" "$KEEP/keep_$t.txt" "$IN/covars_with_sex.txt" 100 --phenoCol "$t"
done
# Height-adjusted weight: 16 covariates (sex, age, height terms, PC1-PC10).
step1 "$S1/weight/weight" "$IN/phenotypes_paper.txt" "$KEEP/keep_weight.txt" \
  "$IN/covars_weight_height_model_pcs.txt" 100 --phenoCol weight
# Log testosterone and its shrunk arms: covariates age, PC1-PC10 (no sex), one sample for all arms.
for t in testosterone testosterone_shrink_0.2 testosterone_shrink_0.5 testosterone_shrink_0.7; do
  step1 "$S1/${t//./_}/${t//./_}" "$IN/phenotypes_paper.txt" "$KEEP/keep_testosterone_series.txt" \
    "$IN/covars_no_sex.txt" 100 --phenoCol "$t"
done

## 2. Step 2 for the main traits (quantitative, block size 1000) ----------------------------
for c in $(seq 1 22); do
  cc=$(chr2 "$c")
  for t in bmi log_bmi; do
    step2 "$S2/$t/${t}_additive_chr$cc" "$c" "$IN/phenotypes_paper.txt" "$IN/covars_with_sex.txt" \
      --phenoCol "$t" --keep "$KEEP/keep_$t.txt" --pred "$S1/$t/${t}_step1_pred.list" --bsize 1000 --qt
  done
  step2 "$S2/weight/weight_additive_chr$cc" "$c" "$IN/phenotypes_paper.txt" "$IN/covars_weight_height_model_pcs.txt" \
    --phenoCol weight --keep "$KEEP/keep_weight.txt" --pred "$S1/weight/weight_step1_pred.list" --bsize 1000 --qt
  # Sex-stratified testosterone: each stratum is tested with the LOCO predictions of the
  # sex-combined step 1 for that arm, covariates age and PC1-PC10.
  for t in testosterone testosterone_shrink_0.2 testosterone_shrink_0.5 testosterone_shrink_0.7; do
    tt=${t//./_}
    for s in male female; do
      step2 "$S2/sex_strat/${tt}_${s}_chr$cc" "$c" "$IN/phenotypes_paper.txt" "$IN/covars_no_sex.txt" \
        --phenoCol "$t" --keep "$KEEP/keep_$s.txt" --pred "$S1/$tt/${tt}_step1_pred.list" --bsize 1000 --qt
    done
  done
done

## 3. GWAS of chromosomal sex ----------------------------------------------------------------
# Reference for the mixture membership parameter of testosterone. Binary trait (male = 1),
# covariates age and PC1-PC10, no LOCO predictions (--ignore-pred), restricted to exactly the
# testosterone Stage 1 sample (keep_testosterone_stage1.txt, section 7 of 01_build_phenotypes.R).
# The earlier full-cohort version is the same command without --keep.
for c in $(seq 1 22); do
  step2 "$S2/sex_gwas_matched/sex_chr$(chr2 "$c")" "$c" "$IN/phenotypes_sex.txt" "$IN/covars_no_sex.txt" \
    --phenoCol sex --bsize 1000 --bt --ignore-pred --keep "$KEEP/keep_testosterone_stage1.txt"
done

## 4. Antibodies: all-sample log MFI ----------------------------------------------------------
# One joint step 1 over every all-sample log-MFI trait with >= 100 values (antibody_qt_all_traits.txt,
# alphabetical), on participants with complete covariates. The serology sample is small (about
# 9,600), so the variant QC uses MAC >= 20 instead of 100. Covariates age, sex, PC1-PC10.
QT_ALL=$(paste -sd, "$IN/antibody_qt_all_traits.txt")
AB_COV=(--covarCol age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
step1 "$S1/antibody_qt_all/antibody_qt_all" "$IN/phenotypes_antibody_all.txt" "$KEEP/keep_antibody_all_complete.txt" \
  "$IN/covars_antibody_all.txt" 20 --phenoColList "$QT_ALL" "${AB_COV[@]}" --qt
# Step 2, all traits in one pass per chromosome (block size 400). The step 1 prediction list
# maps each trait to its LOCO file in the same order as --phenoColList.
for c in $(seq 1 22); do
  cc=$(chr2 "$c")
  step2 "$S2/antibody_qt_all/chr$cc/antigen_qt_chr$cc" "$c" "$IN/phenotypes_antibody_all.txt" "$IN/covars_antibody_all.txt" \
    --pred "$S1/antibody_qt_all/antibody_qt_all_step1_pred.list" --phenoColList "$QT_ALL" "${AB_COV[@]}" \
    --keep "$KEEP/keep_antibody_all_complete.txt" --qt --bsize 400
done

## 5. Antibodies: serostatus and seropositive-only levels (JCV VP1, MCV VP1) ---------------
# Run after section 7 of 01_build_phenotypes.R, which restricts both traits to their Stage 1
# sample. Serostatus is binary (Firth-corrected, approximate, for p < 0.01); the seropositive-only
# trait is log MFI among seropositives. Step 1 QC uses MAC >= 100 here.
SERO=$IN/phenotypes_serostatus.txt
step1 "$S1/serostatus/seropos" "$SERO" "$KEEP/keep_antibody_analysed.txt" "$IN/covars_antibody_all.txt" 100 \
  --phenoColList jcv_vp1_seropos,mcv_vp1_seropos "${AB_COV[@]}" --bt
for t in jcv_vp1 mcv_vp1; do
  step1 "$S1/serostatus/${t}_log_seropos" "$SERO" "$KEEP/keep_${t}_log_all_seropos.txt" "$IN/covars_antibody_all.txt" 100 \
    --phenoCol "${t}_log_seropos" "${AB_COV[@]}" --qt
done
for c in $(seq 1 22); do
  cc=$(chr2 "$c")
  step2 "$S2/serostatus/bt/seropos_chr$cc" "$c" "$SERO" "$IN/covars_antibody_all.txt" \
    --pred "$S1/serostatus/seropos_step1_pred.list" --phenoColList jcv_vp1_seropos,mcv_vp1_seropos "${AB_COV[@]}" \
    --keep "$KEEP/keep_antibody_analysed.txt" --bt --firth --approx --pThresh 0.01 --bsize 400
  for t in jcv_vp1 mcv_vp1; do
    step2 "$S2/serostatus/${t}_log_seropos/log_seropos_chr$cc" "$c" "$SERO" "$IN/covars_antibody_all.txt" \
      --pred "$S1/serostatus/${t}_log_seropos_step1_pred.list" --phenoCol "${t}_log_seropos" "${AB_COV[@]}" \
      --keep "$KEEP/keep_${t}_log_all_seropos.txt" --qt --bsize 400
  done
done

## 6. Step 1 for the survey traits analysed with LOCO -----------------------------------------
# Height, log CRP and log triglycerides; covariates age, sex, PC1-PC10; keep = participants with a
# positive value in the covariate-complete frame. The other survey traits were run without LOCO.
step1 "$S1/height/height" "$IN/phenotypes_survey_core.txt" "$KEEP/keep_height.txt" "$IN/covars_with_sex.txt" 100 --phenoCol height
step1 "$S1/log_crp/log_crp" "$IN/phenotypes_survey_core.txt" "$KEEP/keep_crp.txt" "$IN/covars_with_sex.txt" 100 --phenoCol log_crp
step1 "$S1/log_triglycerides/log_triglycerides" "$IN/phenotypes_survey_core.txt" "$KEEP/keep_triglycerides.txt" \
  "$IN/covars_with_sex.txt" 100 --phenoCol log_triglycerides

## 7. Earlier GWAS used only to select lead variants and to train the BMI score ---------------
# The lead variants of BMI (639) and of testosterone and its shrunk arms, and the summary
# statistics behind the BMI polygenic score, come from an earlier REGENIE run whose settings were recorded with its
# output but whose submission script is not part of this repository. Recorded settings:
#   BMI:          --step 2 --qt --ignore-pred (no LOCO), phenotype bmi, covariates age, sex,
#                 PC1-PC10 (12), HapMap3 rsIDs, all chromosomes in one file, no keep file
#                 (N 461,320 used)
#   testosterone (each arm): --step 2 --qt --ignore-pred, sex-combined, covariates age,
#                 PC1-PC10 (11), HapMap3 rsIDs (N 399,710 used)
# Same-sample LOCO runs of these traits are in sections 1-2 above.
# The per-variant concordance of Fig. 2 likewise used sex-stratified testosterone files from that
# earlier run (199,572 female and 200,138 male participants, the testosterone Stage 1 sample),
# which the Methods describe as step 2 without step 1 predictions; the LOCO versions in section 2
# are the reference for the genetic correlations.
