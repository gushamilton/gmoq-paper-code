#!/usr/bin/env bash
# Lead variants: LD clumping of genome-wide significant variants
#
# What: selects independent lead variants from a GWAS by LD clumping in plink2 against the
#   1000 Genomes European reference: candidates at P < 5e-8, index P < 5e-8, secondary
#   P < 1e-3, r2 < 0.01 within 10 Mb. The same rule is used for every lead set in the paper:
#   on the REGENIE GWAS of BMI (639 leads, also used for log BMI and every BMI scale),
#   height-adjusted weight, testosterone and its shrunk arms, and the JCV/MCV VP1 antibodies;
#   and, for survey traits without a REGENIE GWAS, on the P value of the Stage 1 location effect
#   delta_mu = sum_t c_t beta(u_t) (computed and clumped in analysis/01_data/lead_variants.R).
# Paper: Methods (lead variants); lead sets used in Figs. 1-6
# Input: REGENIE .regenie file(s) with columns ID CHROM GENPOS ... LOG10P
# Output: <out_dir>/<name>.sumstats (candidates), <name>_clump.clumps (plink2), <name>_leads.txt
# Runs: on a cluster (reads genome-wide summary statistics). PLINK v2.0.0-a.6.8LM.
set -euo pipefail

# ---------------------------------------------------------------- user settings --
LD_REF=/path/to/1000G/EUR          # 1000 Genomes phase 3 European, PLINK bed/bim/fam prefix
OUT=/path/to/leads
mkdir -p "$OUT"

clump_leads() {   # name file1.regenie [file2.regenie ...]
  local name=$1; shift
  local ss=$OUT/$name.sumstats
  # 1. Candidates: LOG10P >= -log10(5e-8). P is written as awk prints exp(-LOG10P * log(10)),
  #    i.e. %.6g, and underflows to 0 for the strongest signals, which plink2 still treats as
  #    the smallest P values. Columns are located by name in the REGENIE header.
  awk 'BEGIN { OFS = "\t"; thr = -log(5e-8) / log(10); print "SNP", "CHR", "BP", "P", "LOG10P" }
       FNR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
       $h["LOG10P"] >= thr { print $h["ID"], $h["CHROM"], $h["GENPOS"], exp(-$h["LOG10P"] * log(10)), $h["LOG10P"] }' \
    "$@" > "$ss"
  # 2. Clump: index variants P < 5e-8, variants within 10,000 kb with r2 >= 0.01 and P < 1e-3
  #    are assigned to the index variant's clump.
  plink2 --bfile "$LD_REF" --clump "$ss" --clump-id-field SNP --clump-p-field P \
    --clump-p1 5e-8 --clump-p2 1e-3 --clump-r2 0.01 --clump-kb 10000 --out "$OUT/${name}_clump"
  # 3. Leads are the index variants of the clumps (column 3 of the .clumps file is ID).
  awk 'BEGIN { OFS = "\t"; print "ID", "CHROM", "GENPOS", "LOG10P" }
       NR == FNR { if (FNR > 1) keep[$3] = 1; next }
       FNR > 1 && ($1 in keep) { print $1, $2, $3, $5 }' "$OUT/${name}_clump.clumps" "$ss" > "$OUT/${name}_leads.txt"
}

## 1. Lead sets from REGENIE GWAS -------------------------------------------------------------
# BMI: the earlier all-chromosome REGENIE BMI file (02_regenie.sh, section 7) -> 639 leads.
clump_leads bmi /path/to/regenie/bmi_all.regenie
# Height-adjusted weight (16 covariates, LOCO), one file per chromosome.
clump_leads weight /path/to/gwas_inputs/regenie_step2/weight/weight_additive_chr*_weight.regenie
# Testosterone and its shrunk arms: the earlier sex-combined REGENIE files.
for t in testosterone testosterone_shrink_0_2 testosterone_shrink_0_5 testosterone_shrink_0_7; do
  clump_leads "$t" "/path/to/regenie/${t}_all.regenie"
done
# Antibodies: all-sample log MFI, REGENIE step 2 per chromosome.
for t in jcv_vp1_log_all mcv_vp1_log_all; do
  clump_leads "$t" /path/to/gwas_inputs/regenie_step2/antibody_qt_all/chr*/antigen_qt_chr*_"$t".regenie
done

## 2. Survey traits ------------------------------------------------------------------------------
# Same clumping call on the Stage 1 location effect: P = 2 Phi(-|z|), z = delta_mu / SE with
# SE = sqrt(c' Sigma_0 c / d_resid), written at full precision, all variants with a finite P.
# See analysis/01_data/lead_variants.R.
