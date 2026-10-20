#!/bin/bash
# Linkage disequilibrium of FUT2-region variants with the secretor variant rs601338
#
# What: r2, D' and phase in the 1000 Genomes European reference between the FUT2 nonsense
#   variant rs601338 (non-secretor allele A) and the JCV FUT2-region lead rs681343, and the
#   FUT2 variants reported by earlier serology GWAS (rs2432132, rs492602, rs12462111), with
#   their allele frequencies. Only the plink logs and frequency files are kept.
# Paper: Supplementary Note 3 (rs681343 with rs601338: r2 = 0.996, D' = 1)
# Input: 1000 Genomes European genotypes (config.yaml p08.plink_bfile)
# Output: results/04_antibodies/fut2_threshold/ld_1000g/{ld_<snp>_rs601338.log, freq_fut2_lit.afreq}
# Runs: on a cluster node with plink2, in under a minute
set -euo pipefail
repo=${GMOQ_REPO:-$(git rev-parse --show-toplevel)}
cfg() { Rscript -e "cat(yaml::read_yaml('$repo/config.yaml')\$$1)"; }
out=${GMOQ_OUT:-$repo/results}/04_antibodies/fut2_threshold/ld_1000g
mkdir -p "$out"; cd "$out"
module load "$(cfg 'p08$plink_module')"
reference=$(cfg 'p08$plink_bfile')
for snp in rs2432132 rs492602 rs12462111 rs681343; do
  plink2 --bfile "$reference" --ld "$snp" rs601338 --out "ld_${snp}_rs601338" || true   # a variant absent from the reference is skipped
done
plink2 --bfile "$reference" --snps rs2432132,rs492602,rs12462111,rs681343,rs601338 --freq --out freq_fut2_lit
