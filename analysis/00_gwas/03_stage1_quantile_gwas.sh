#!/usr/bin/env bash
# Stage 1 of G-MoQ: unconditional quantile GWAS by RIF regression at 45 quantile levels
#
# What: runs the gmoq-stage1 command-line program for every trait and chromosome. For
#   each variant it estimates the effect on the 45 phenotype quantiles u = 0.02, ..., 0.98
#   (recentred influence function regression on the covariates and, where used, the REGENIE
#   LOCO prediction as one more covariate) and writes the per-variant profile beta(u), the
#   residual genotype variance d_resid and the pooled quantile-by-quantile covariance Sigma_0
#   that Stage 2 uses. Settings differ by trait group and are shown explicitly. A last section
#   runs Stage 1 at the 639 BMI lead variants only, for the polygenic-score analyses.
# Paper: Methods (Stage 1); input to every Stage 2 result (Figs. 1-6)
# Input: phenotype, covariate and keep files (01_build_phenotypes.R), REGENIE step 1 LOCO
#   prediction lists (02_regenie.sh), imputed BGEN, HapMap3 variant lists
# Output: <out_dir>/<trait>/chrNN/scan.stage1.tsv.gz, scan.sigma0.tsv.gz, scan.sigma0.meta.json
#   (and rif_nuisance.json, the Stage 1 quantiles and kernel density read by the Stage 2 weights)
# Runs: on a cluster with UK Biobank data. gmoq-stage1 at the package commit given in
#   config.yaml (package.stage1). Chromosomes (and traits) were run in parallel on the cluster;
#   here they are loops. Threads: 4 per job.
set -euo pipefail

# ---------------------------------------------------------------- user settings --
IN=/path/to/gwas_inputs                       # out_dir of 01_build_phenotypes.R and 02_regenie.sh
OUT=/path/to/stage1
BGEN_DIR=/path/to/ukb_imputed                 # data.chrNN.bgen
SAMPLE=/path/to/ukb_imputed/data.chr1-22.sample
HM3=/path/to/hapmap3/w_hm3_chrNN_only.txt     # HapMap3 rsIDs on chromosome NN
LEAD_SNPS=/path/to/bmi_leads/chrNN.txt        # the 639 BMI lead rsIDs on chromosome NN
PREP=/path/to/prs_prep                        # keep files of the polygenic-score analysis (section 6)
PY=python                                     # Python with gmoq installed
THREADS=4
export OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
KEEP=$IN/keep; S1=$IN/regenie_step1
chr2() { printf "%02d" "$1"; }
hm3() { echo "${HM3//NN/$(chr2 "$1")}"; }

# The 45 quantile levels, linspace(0.02, 0.98, 45), passed at full precision (%.17g).
TAUS=$("$PY" -c 'import numpy as np; print(",".join(f"{x:.17g}" for x in np.linspace(.02, .98, 45)))')

stage1() {   # out_dir chr pheno_file pheno_col covar_file keep_file snp_list  [flags...]
  local od=$1 c=$2 pheno=$3 col=$4 covar=$5 keep=$6 snps=$7; shift 7
  mkdir -p "$od"
  "$PY" -m gmoq.cli --bgen "$BGEN_DIR/data.chr$(chr2 "$c").bgen" --sample "$SAMPLE" \
    --pheno "$pheno" --pheno-col "$col" --covar "$covar" --keep "$keep" --snps "$snps" \
    --threads "$THREADS" --out "$od/scan" "$@"
}

# Stage 1 kernel density and quantiles on the matched sample, saved as rif_nuisance.json.
# Same sample rule as Stage 1 (keep -> complete covariates -> LOCO -> BGEN sample) and the
# same Silverman Gaussian kernel. With --maf-out, also writes the variants with in-sample
# dosage MAF > 0.05 (the antibody variant list).
nuisance() {   # out_dir chr pheno_file pheno_col covar_file keep_file loco_list snp_list [--maf-out]
  "$PY" - "$@" <<'PY'
import json, os, sys
import numpy as np
from gmoq import cli
od, c, pheno, col, covar, keep, loco, snps = sys.argv[1:9]
maf_out = "--maf-out" in sys.argv[9:]
taus = np.linspace(.02, .98, 45)
bgen = os.path.join(os.environ["BGEN_DIR"], f"data.chr{int(c):02d}.bgen")
ids, aliases, ys, _, _ = cli.load_analysis_inputs(pheno, [col], covar, keep_file=keep)
_, loco_ok = cli.load_loco_predictions(loco, [col], aliases)
ids = [x for x, ok in zip(ids, loco_ok) if ok]
y = np.asarray(ys[col], dtype=float)[loco_ok]
variants, samples, handle = cli.get_variants_bgen(bgen, os.environ["SAMPLE"], cli.load_snps_list(snps))
valid, indices = cli.match_samples(samples, ids)
y = y[valid]
_, q, f = cli.compute_rif(y, taus, density_method="silverman_gaussian")
os.makedirs(od, exist_ok=True)
json.dump({"n": len(y), "mean": float(y.mean()), "sd": float(y.std(ddof=1)), "taus": taus.tolist(),
           "q_tau": q.tolist(), "fhat": f.tolist()}, open(os.path.join(od, "rif_nuisance.json"), "w"), indent=2, sort_keys=True)
if maf_out:
    keep_ids = []
    for info, _, dosage in cli._decode_variants(variants, bgen, indices, cli._SkipCount()):
        g = dosage[np.isfinite(dosage)]
        p = g.sum() / (2 * len(g))
        if min(p, 1 - p) > .05:
            keep_ids.append(info[0])
    open(os.path.join(od, "maf_gt_005.snps"), "w").write("\n".join(keep_ids) + "\n")
handle.close()
PY
}
export BGEN_DIR SAMPLE

## 1. BMI, log BMI and testosterone (with its shrunk arms) ------------------------------------
# LOCO prediction as a covariate; full jackknife covariance output (--cov-mode full).
COMMON=(--taus "$TAUS" --blocks 1000 --batch-size 100 --cov-dtype float32 --seed 42)
for c in $(seq 1 22); do
  for t in bmi log_bmi; do
    od=$OUT/$t/chr$(chr2 "$c"); loco=$S1/$t/${t}_step1_pred.list
    nuisance "$od" "$c" "$IN/phenotypes_paper.txt" "$t" "$IN/covars_with_sex.txt" "$KEEP/keep_$t.txt" "$loco" "$(hm3 "$c")"
    stage1 "$od" "$c" "$IN/phenotypes_paper.txt" "$t" "$IN/covars_with_sex.txt" "$KEEP/keep_$t.txt" "$(hm3 "$c")" \
      --loco-preds "$loco" --loco-mode covariate --cov-mode full "${COMMON[@]}"
  done
  for t in testosterone testosterone_shrink_0.2 testosterone_shrink_0.5 testosterone_shrink_0.7; do
    tt=${t//./_}; od=$OUT/$tt/chr$(chr2 "$c"); loco=$S1/$tt/${tt}_step1_pred.list
    # Covariates age and PC1-PC10: sex is left out so the two sex components stay in the distribution.
    nuisance "$od" "$c" "$IN/phenotypes_paper.txt" "$t" "$IN/covars_no_sex.txt" "$KEEP/keep_testosterone_series.txt" "$loco" "$(hm3 "$c")"
    stage1 "$od" "$c" "$IN/phenotypes_paper.txt" "$t" "$IN/covars_no_sex.txt" "$KEEP/keep_testosterone_series.txt" "$(hm3 "$c")" \
      --loco-preds "$loco" --loco-mode covariate --cov-mode full "${COMMON[@]}"
  done
done

## 2. Height-adjusted weight --------------------------------------------------------------------
# 16 covariates (sex, age, height terms, PC1-PC10); LOCO as a covariate; pooled covariance only
# (--cov-mode pooled: Sigma_0 and per-variant SEs, no per-variant jackknife covariance file).
for c in $(seq 1 22); do
  od=$OUT/weight/chr$(chr2 "$c"); loco=$S1/weight/weight_step1_pred.list
  nuisance "$od" "$c" "$IN/phenotypes_paper.txt" weight "$IN/covars_weight_height_model_pcs.txt" "$KEEP/keep_weight.txt" "$loco" "$(hm3 "$c")"
  stage1 "$od" "$c" "$IN/phenotypes_paper.txt" weight "$IN/covars_weight_height_model_pcs.txt" "$KEEP/keep_weight.txt" "$(hm3 "$c")" \
    --loco-preds "$loco" --loco-mode covariate --cov-mode pooled "${COMMON[@]}"
done

## 3. JC and Merkel cell polyomavirus VP1 antibodies ---------------------------------------------
# All-sample log MFI. Because the sample is small (9,145), only HapMap3 variants with in-sample
# MAF > 0.05 are scanned (computed from dosages on the matched sample by nuisance --maf-out).
# LOCO from the joint antibody step 1, as a covariate, regression-scaled; pooled covariance only
# (--cov-mode pooled-only). Covariates age, sex, PC1-PC10.
for t in jcv_vp1_log_all mcv_vp1_log_all; do
  pred=$OUT/$t/${t}_step1_pred.list; mkdir -p "$OUT/$t"
  grep "^$t " "$S1/antibody_qt_all/antibody_qt_all_step1_pred.list" > "$pred"   # this trait's LOCO file only
  for c in $(seq 1 22); do
    od=$OUT/$t/chr$(chr2 "$c")
    nuisance "$od" "$c" "$IN/phenotypes_antibody_all.txt" "$t" "$IN/covars_antibody_all.txt" "$KEEP/keep_$t.txt" "$pred" "$(hm3 "$c")" --maf-out
    stage1 "$od" "$c" "$IN/phenotypes_antibody_all.txt" "$t" "$IN/covars_antibody_all.txt" "$KEEP/keep_$t.txt" "$od/maf_gt_005.snps" \
      --loco-preds "$pred" --loco-mode covariate --loco-scale regression --cov-mode pooled-only "${COMMON[@]}"
  done
done

## 4. Trait survey -----------------------------------------------------------------------------------
# Covariates age, sex, PC1-PC10; pooled covariance; HapMap3. LOCO (as a covariate) only for height,
# log CRP and log triglycerides; every other survey trait was run without LOCO (covariates only).
# These runs did not pass --taus or --cov-dtype, so the program defaults applied: the same 45
# levels written to 3 decimals, and int8 storage of the jackknife covariance.
# Phenotype column : keep file : phenotype file (keep_panel = biochemistry/blood-count panel)
SURVEY=(
  "height:keep/keep_height:core"            "log_crp:keep/keep_crp:core"
  "log_triglycerides:keep/keep_triglycerides:core"
  "albumin:keep_panel/keep_albumin:panel"   "fev1_fvc_ratio:keep_panel/keep_fev1_fvc_ratio:panel"
  "hba1c:keep_panel/keep_hba1c_nodiab:panel" "ldl:keep_panel/keep_ldl:panel"
  "log_alt:keep_panel/keep_alt:panel"       "log_bilirubin_total:keep_panel/keep_bilirubin_total:panel"
  "log_calcium:keep_panel/keep_calcium:panel" "log_ggt:keep_panel/keep_ggt:panel"
  "log_hba1c:keep_panel/keep_hba1c:panel"   "log_hdl:keep_panel/keep_hdl:panel"
  "log_igf1:keep_panel/keep_igf1:panel"     "log_lipoprotein_a:keep_panel/keep_lipoprotein_a:panel"
  "log_total_protein:keep_panel/keep_total_protein:panel" "log_urate:keep_panel/keep_urate:panel"
  "mch:keep_panel/keep_mch:panel"           "monocyte_pct:keep_panel/keep_monocyte_pct:panel"
  "phosphate:keep_panel/keep_phosphate:panel"
)
for c in $(seq 1 22); do
  for row in "${SURVEY[@]}"; do
    IFS=: read -r col keep src <<< "$row"
    tag=$col; [[ "$keep" == *hba1c_nodiab ]] && tag=hba1c_nodiab     # raw HbA1c, diabetes excluded
    loco=()
    case $col in height|log_crp|log_triglycerides)
      loco=(--loco-preds "$S1/$col/${col}_step1_pred.list" --loco-mode covariate) ;; esac
    stage1 "$OUT/survey/$tag/chr$(chr2 "$c")" "$c" "$IN/phenotypes_survey_$src.txt" "$col" "$IN/covars_with_sex.txt" \
      "$IN/$keep.txt" "$(hm3 "$c")" ${loco[@]+"${loco[@]}"} --blocks 1000 --batch-size 100 --seed 42 --cov-mode pooled
  done
done

## 5. Leads only: log BMI without LOCO and with LOCO as a covariate --------------------------------
# The polygenic-score analyses use Stage 1 without LOCO, because the LOCO prediction is itself a
# polygenic score and would absorb part of the interaction under study. Both arms are run at the
# 639 BMI lead variants on the BMI Stage 1 sample (keep_bmi_stage1.txt; N 461,320), covariates
# age, sex, PC1-PC10, pooled covariance, program defaults for --taus and --cov-dtype.
LEADS_COMMON=(--blocks 1000 --batch-size 100 --seed 42 --cov-mode pooled)
lead_list() { echo "${LEAD_SNPS//NN/$(chr2 "$1")}"; }
for c in $(seq 1 22); do
  [[ -s "$(lead_list "$c")" ]] || continue                     # chromosomes without a lead
  cc=$(chr2 "$c")
  stage1 "$OUT/leads/log_bmi_none/chr$cc" "$c" "$IN/phenotypes_paper.txt" log_bmi "$IN/covars_with_sex.txt" \
    "$KEEP/keep_bmi_stage1.txt" "$(lead_list "$c")" "${LEADS_COMMON[@]}"
  stage1 "$OUT/leads/log_bmi_covariate/chr$cc" "$c" "$IN/phenotypes_paper.txt" log_bmi "$IN/covars_with_sex.txt" \
    "$KEEP/keep_bmi_stage1.txt" "$(lead_list "$c")" \
    --loco-preds "$S1/log_bmi/log_bmi_step1_pred.list" --loco-mode covariate "${LEADS_COMMON[@]}"

## 6. Leads only: arms of the polygenic-score analysis ----------------------------------------------
# Keep files come from the polygenic-score preparation (analysis/05_bmi): the same BMI sample
# (keep_all), its two sexes (keep_sex0, keep_sex1) and a random split into halves (keep_half1,
# keep_half2), used for the cross-half link test on log and raw BMI.
  stage1 "$OUT/leads/log_bmi_subtract/chr$cc" "$c" "$IN/phenotypes_paper.txt" log_bmi "$IN/covars_with_sex.txt" \
    "$PREP/keep_all.txt" "$(lead_list "$c")" \
    --loco-preds "$S1/log_bmi/log_bmi_step1_pred.list" --loco-mode subtract --loco-scale regression "${LEADS_COMMON[@]}"
  for s in 0 1; do   # within one sex, sex is constant and is dropped from the covariates
    stage1 "$OUT/leads/log_bmi_sex$s/chr$cc" "$c" "$IN/phenotypes_paper.txt" log_bmi "$IN/covars_with_sex.txt" \
      "$PREP/keep_sex$s.txt" "$(lead_list "$c")" --covar-col age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 "${LEADS_COMMON[@]}"
  done
  for h in 1 2; do
    for t in log_bmi bmi; do
      stage1 "$OUT/leads/${t}_h$h/chr$cc" "$c" "$IN/phenotypes_paper.txt" "$t" "$IN/covars_with_sex.txt" \
        "$PREP/keep_half$h.txt" "$(lead_list "$c")" "${LEADS_COMMON[@]}"
    done
  done
done
