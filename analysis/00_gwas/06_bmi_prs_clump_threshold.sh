#!/bin/bash
# BMI polygenic score by PLINK clumping and thresholding, one chromosome at a time
#
# What: clumps the REGENIE BMI GWAS against the 1000 Genomes European reference (r2 < 0.1
#   within 250 kb), keeps the index variants with P < 1e-5, and scores every UK Biobank
#   participant per chromosome with their GWAS betas (a "partial score"). For a lead variant on
#   chromosome c, the interaction analysis (analysis/05_bmi/07_prs_scores.R) uses the sum of
#   the other 21 partial scores, standardised, so the score contains no variant from the lead's
#   own chromosome. The weights come from a GWAS of the same participants (an in-sample score).
# Paper: Fig. 6; Supplementary Fig. S14; Supplementary Table S6; Methods (BMI x PRS interaction)
# Input: REGENIE BMI GWAS (all chromosomes in one file); 1000 Genomes EUR PLINK files; UK
#   Biobank imputed BGEN and sample file
# Output: <out_dir>/partial_scores/chrNN_partial.sscore (BETA_SUM per participant; stays on
#   the cluster)
# Runs: on a cluster with UK Biobank data; about 10 seconds and under 1 GB per chromosome.
#   PLINK 2.00a6.8LM, bgenix 1.1.7.

set -euo pipefail

# ---------------------------------------------------------------- user settings --
GWAS=/path/to/regenie/bmi_all.regenie        # columns CHROM GENPOS ID ALLELE0 ALLELE1 ... TEST BETA SE CHISQ LOG10P
REF=/path/to/1000G/EUR                        # PLINK prefix, 1000 Genomes phase 3 European
BGEN=/path/to/ukb_imputed                     # data.chrNN.bgen (+ .bgi), data.chr1-22.sample
OUT=/path/to/bmi_score
THRESHOLD=1e-5
MEM=4000                                      # MB; plink2 must be told its budget under a scheduler

mkdir -p "$OUT/partial_scores" "$OUT/work"
for CHR in $(seq 1 22); do
  C=$(printf "%02d" "$CHR"); W="$OUT/work/chr$C"; mkdir -p "$W"

  # 1. Summary statistics for this chromosome (additive test): SNP, effect allele, beta, P.
  awk -v c="$CHR" 'BEGIN{OFS="\t"; print "SNP","A1","BETA","P"}
       NR>1 && $1==c && $9=="ADD" {p=exp(-$13*log(10)); if (p<1e-300) p=1e-300; print $3,$5,$10,p}' "$GWAS" > "$W/ss.tsv"

  # 2. Clump: index variants with P < 1e-3, no two within 250 kb at r2 >= 0.1.
  plink2 --memory $MEM --threads 1 --bfile "$REF" --chr "$CHR" --clump "$W/ss.tsv" \
    --clump-id-field SNP --clump-p-field P --clump-p1 1e-3 --clump-p2 1e-3 --clump-r2 0.1 --clump-kb 250 \
    --out "$W/clump"

  # 3. Weights: index variants below the threshold, with their GWAS beta for the effect allele.
  awk 'NR==1{for(i=1;i<=NF;i++) if($i=="ID") k=i; next} {print $k}' "$W/clump.clumps" > "$W/index.txt"
  awk -v t="$THRESHOLD" 'NR==FNR{keep[$1]=1; next} FNR==1{print "SNP\tA1\tBETA"; next}
       ($1 in keep) && $4<t {print $1"\t"$2"\t"$3}' "$W/index.txt" "$W/ss.tsv" > "$W/weights.tsv"
  tail -n +2 "$W/weights.tsv" | cut -f1 > "$W/rsids.txt"

  # 4. Score every participant: extract the variants, convert, sum dosage x beta.
  if [[ -s "$W/rsids.txt" ]]; then
    bgenix -g "$BGEN/data.chr$C.bgen" -incl-rsids "$W/rsids.txt" > "$W/sub.bgen"
    plink2 --memory $MEM --threads 1 --bgen "$W/sub.bgen" ref-first --sample "$BGEN/data.chr1-22.sample" \
      --make-pgen --out "$W/sub"
    plink2 --memory $MEM --threads 1 --pfile "$W/sub" --rm-dup force-first \
      --score "$W/weights.tsv" 1 2 3 header-read cols=+scoresums --out "$OUT/partial_scores/chr${C}_partial"
    rm -f "$W/sub.bgen" "$W/sub.pgen" "$W/sub.pvar" "$W/sub.psam"
  else   # no variant below the threshold on this chromosome: a zero partial score
    awk 'NR>2{print $1"\t"$2"\t0"}' "$BGEN/data.chr1-22.sample" | sed '1i #FID\tIID\tBETA_SUM' \
      > "$OUT/partial_scores/chr${C}_partial.sscore"
  fi
done
