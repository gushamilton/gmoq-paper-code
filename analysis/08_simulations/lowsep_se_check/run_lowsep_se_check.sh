#!/usr/bin/env bash
# Standard errors at the empirical mixture separations: run all parts
#
# What: do pooled Stage 1 covariance standard errors agree with 1,000-block jackknife
#   standard errors for the mixture parameters (membership and the two component
#   locations) at the separations of the paper's empirical mixtures? Part A reruns Stage 1
#   for JCV/MCV VP1 on two chromosomes keeping the jackknife payload and fits Stage 2 with
#   both covariances; part B does the same Stage 2 comparison for the four testosterone
#   mixtures from stored Stage 1 output; part C aggregates the ratios and the separations;
#   part D runs the simulation; part E joins the real-data and simulated arms.
# Paper: Supplementary Information, study 7; Supplementary Fig. S4; Supplementary Tables S1, S2
# Input: UK Biobank imputed BGEN and the 44-antigen scan inputs (part A); stored Stage 1
#   output and mixture weights (parts B, C); the G-MoQ packages and the simulation harness
#   in ../code/validation/.
# Output: ../inputs/lowsep_se_check/{empirical_se_ratio,separations,reproducibility_checks,
#   sim_summary,lowsep_se_summary}.csv
# Runs: on a cluster with UK Biobank data; part D about 3 hours on 32 cores
set -euo pipefail

## User-set paths --------------------------------------------------------------------
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=/path/to/lowsep_se_check_output             # per-variant Stage 2 reruns (not shared)
PYTHON=/path/to/gmoq-env/bin/python   # Python with the gmoq package
RSCRIPT=Rscript                                 # R with the gmoq package
STAGE2_LIB=/path/to/R/library/with/gmoq
STAGE2_CHR_SCRIPT=/path/to/run_stage2_chr.R     # per-chromosome Stage 2 driver (G-MoQ scripts/stage2/)
VALIDATION_DIR=$HERE/../code/validation        # simulation harness (run_examples.py, truth.py, stage2.R)
ENV_MANIFEST=/path/to/environment.json          # names the Python, Rscript, R library and Stage 1 build
THREADS=8
WORKERS=32

BGEN_DIR=/path/to/ukb_imputed                   # data.chrNN.bgen
SAMPLE=/path/to/ukb_imputed/data.chr1-22.sample
ANTIBODY_DIR=/path/to/antibody_phenotypes       # phenotypes/covariates/keep files of the 44-antigen scan
SCAN_DIR=/path/to/antigen_scan                  # 44-antigen scan output: <trait>/chrNN/{maf_gt_005.snps,weights/,stage2/}
CORRECTED_DIR=/path/to/testosterone_model_density   # corrected_weights/<trait>/chrNN/weights_skew.rds, stage2_corrected/
OLD_RUN_TEMPLATES='/path/to/earlier_run_%1$s/full/%1$s,/path/to/earlier_run_reused_chromosomes/full/%1$s'

export OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT"
TAUS=$("$PYTHON" -c 'import numpy as np; print(",".join(f"{x:.17g}" for x in np.linspace(.02,.98,45)))')

# Stage 2 twice on the same Stage 1 betas and weights: pooled covariance, then jackknife.
stage2 () {  # trait model chr stage1_file weights_file covariance_prefix out_dir
  local t=$1 m=$2 c=$3 s1=$4 w=$5 pre=$6 od=$7
  for src in pooled jackknife; do
    mkdir -p "$od/$src"
    extra=(); [[ $src == jackknife ]] && extra=(--cov-file "$pre.cov.gz")
    "$RSCRIPT" --vanilla "$STAGE2_CHR_SCRIPT" --trait "$t" --model "$m" --chr "$c" \
      --stage1 "$s1" --weights "$w" --out "$od/$src/result" --cov-source "$src" \
      --estimator ols --threads "$THREADS" --r-lib "$STAGE2_LIB" \
      --sigma0-file "$pre.sigma0.tsv.gz" --cov-meta-file "$pre.cov.meta.json" "${extra[@]}" \
      > "$od/$src/stage2.log" 2>&1
  done
}

## A. JCV/MCV VP1: Stage 1 rerun with the full (jackknife) covariance -----------------
# Same inputs and seed as the 44-antigen scan (variants with MAF > 0.05), so the pooled
# arm reproduces the stored scan exactly; then Stage 2 with the scan's own mixture weights.
for t in jcv_vp1_log_all mcv_vp1_log_all; do
  for c in 6 22; do
    cc=$(printf %02d $c); src="$SCAN_DIR/$t/chr$cc"; w="$OUT/antibody/$t/chr$cc"; mkdir -p "$w"
    "$PYTHON" -m gmoq.cli \
      --bgen "$BGEN_DIR/data.chr$cc.bgen" --sample "$SAMPLE" \
      --pheno "$ANTIBODY_DIR/phenotypes_antibody_all.txt" --pheno-col "$t" \
      --covar "$ANTIBODY_DIR/covars_antibody_all.txt" --keep "$ANTIBODY_DIR/keep_files/keep_$t.txt" \
      --loco-preds "$SCAN_DIR/$t/${t}_step1_pred.list" --snps "$src/maf_gt_005.snps" \
      --taus "$TAUS" --blocks 1000 --batch-size 100 --cov-dtype float32 --cov-mode full \
      --loco-mode covariate --loco-scale regression --seed 42 --threads "$THREADS" \
      --out "$w/scan" > "$w/stage1.log" 2>&1
    stage2 "$t" mixture "$c" "$w/scan.stage1.tsv.gz" "$src/weights/weights_mixture.rds" "$w/scan" "$w/stage2"
  done
done

## B. Testosterone x4 on the model-density weights, chromosomes 19-22 -----------------
# Stage 2 from the stored Stage 1 and its stored full covariance (scan.cov.gz).
for t in testosterone testosterone_shrink_0.2 testosterone_shrink_0.5 testosterone_shrink_0.7; do
  for c in 19 20 21 22; do
    cc=$(printf %02d $c)
    wf="$CORRECTED_DIR/corrected_weights/$t/chr$cc/weights_skew.rds"
    s1=$("$RSCRIPT" --vanilla -e "cat(readRDS('$wf')\$provenance\$stage1_file)")
    pre="${s1%.stage1.tsv.gz}"
    stage2 "$t" skew "$c" "$s1" "$wf" "$pre" "$OUT/testosterone_corrected/$t/chr$cc"
  done
done

## C. Aggregate the real-data ratios and the separations ------------------------------
"$RSCRIPT" --vanilla "$HERE/empirical_se_ratio.R" --run-dir "$OUT" --out-dir "$OUT/aggregate" \
  --old-run-templates "$OLD_RUN_TEMPLATES" --corrected-dir "$CORRECTED_DIR" --scan-dir "$SCAN_DIR"

## D. Simulation: 40 cells x 200 replicates, one thread per worker --------------------
OMP_NUM_THREADS=1 "$PYTHON" "$HERE/sim_lowsep_se.py" --out-dir "$OUT/sim" \
  --validation-dir "$VALIDATION_DIR" --environment-manifest "$ENV_MANIFEST" \
  --workers "$WORKERS" --blocks 1000 --replicates 200

## E. Combined table; copy the aggregates next to this script --------------------------
"$PYTHON" "$HERE/assemble_summary.py" --empirical-dir "$OUT/aggregate" \
  --sim-summary "$OUT/sim/sim_summary.csv" --out "$OUT/lowsep_se_summary.csv"
DEST="$HERE/../inputs/lowsep_se_check"
cp "$OUT"/aggregate/{empirical_se_ratio,separations,reproducibility_checks}.csv \
   "$OUT/sim/sim_summary.csv" "$OUT/lowsep_se_summary.csv" "$DEST/"
