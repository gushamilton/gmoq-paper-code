#!/bin/bash
# Genetic correlation of testosterone between women and men
#
# What: LD score regression genetic correlation between the female-only and male-only REGENIE
#   GWAS of log testosterone (the 'combined' sex-stratified reference: 199,572 women, 200,138
#   men). It shows that a variant's effect in one sex says almost nothing about its effect in
#   the other (rg = 0.015, SE 0.043), which is why the two mixture components must be estimated
#   separately. Only the LDSC logs are kept.
# Paper: Results (testosterone; rg between the sex-stratified GWAS)
# Input: sex-stratified REGENIE GWAS (config.yaml p08.regenie_combined); HapMap3 alleles and
#   1000 Genomes European LD scores (p08.ldsc_merge_alleles, p08.ldsc_ld_dir)
# Output: results/03_testosterone/checks/{comb_female,comb_male,rg_combined_female_male}.log
# Runs: on a cluster node with LDSC, in about 15 minutes
set -euo pipefail
repo=${GMOQ_REPO:-$(git rev-parse --show-toplevel)}
cfg() { Rscript -e "cat(yaml::read_yaml('$repo/config.yaml')\$$1)"; }
combined=$(cfg 'p08$regenie_combined')
work=$(cfg 'bp$private_root')/gmoq-paper-code/03_testosterone/male_female_rg
out=${GMOQ_OUT:-$repo/results}/03_testosterone/checks
mkdir -p "$work" "$out"; cd "$work"
module load "$(cfg 'p08$ldsc_module')"
for sex in female male; do
  # REGENIE columns: 3 ID, 4 ALLELE0, 5 ALLELE1 (effect allele), 8 N, 10 BETA, 13 LOG10P.
  awk 'BEGIN{OFS="\t"; print "SNP","A1","A2","N","P","BETA"}
       NR>1 && $10!="NA" && $13!="NA" {p=10^(-$13); if(p<1e-300)p=1e-300; print $3,$5,$4,$8,p,$10}' \
      "$combined/testosterone_${sex}_testosterone_all.regenie" > raw_${sex}.tsv
  munge_sumstats.py --sumstats raw_${sex}.tsv --out comb_${sex} --merge-alleles "$(cfg 'p08$ldsc_merge_alleles')" \
    --snp SNP --a1 A1 --a2 A2 --p P --N-col N --signed-sumstats BETA,0 > /dev/null
  rm -f raw_${sex}.tsv
done
ld=$(cfg 'p08$ldsc_ld_dir')
ldsc.py --rg comb_female.sumstats.gz,comb_male.sumstats.gz --ref-ld-chr "$ld" --w-ld-chr "$ld" \
  --out rg_combined_female_male > /dev/null
cp comb_female.log comb_male.log rg_combined_female_male.log "$out/"
