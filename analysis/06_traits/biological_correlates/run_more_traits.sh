#!/bin/bash
# Run inside a BP compute allocation after supplying the inputs listed in README.md.
# Participant-level outputs stay in GMOQ_PRIVATE_ROOT. No allocation is requested here.
set -eo pipefail

ROOT=${GMOQ_PRIVATE_ROOT:?Set GMOQ_PRIVATE_ROOT outside this repository}
PY=${PYTHON:-python3}
: "${BGENIX:?Set BGENIX to the indexed BGEN extraction executable}"
BG=${GMOQ_BGEN_DIR:?Set GMOQ_BGEN_DIR to the imputed BGEN directory}
COVAR=${GMOQ_COVARIATES:?Set GMOQ_COVARIATES to the covariate table}
SAMPLE=${GMOQ_BGEN_SAMPLE:?Set GMOQ_BGEN_SAMPLE to the BGEN sample file}
: "${PLINK_MEMORY_MB:?Set PLINK_MEMORY_MB to memory available to PLINK in this allocation}"
THREADS=${SLURM_CPUS_PER_TASK:-1}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

echo "[aux_more2d] start $(date -Is) host=$(hostname) job=$SLURM_JOB_ID"
"$PY" "$SCRIPT_DIR/run_more_traits.py" --mode prep --root "$ROOT"

set +u
source /etc/profile >/dev/null 2>&1 || true
module load apps/plink2/2.00a68LM
set -u
mkdir -p "$ROOT/geno" "$ROOT/glm"

for chr in $(seq 1 22); do
  c=$(printf %02d "$chr")
  L="$ROOT/snps/union_chr${c}.txt"
  [[ -s "$L" ]] || continue
  "$BGENIX" -g "$BG/data.chr${c}.bgen" -incl-rsids $(tr '\n' ' ' < "$L") > "$ROOT/geno/union_chr${c}.bgen"
  plink2 --bgen "$ROOT/geno/union_chr${c}.bgen" ref-first --sample "$SAMPLE" \
    --keep "$ROOT/keep/keep_union.txt" --make-pgen --threads "$THREADS" --memory "$PLINK_MEMORY_MB" \
    --out "$ROOT/geno/union_chr${c}" > /dev/null
  rm -f "$ROOT/geno/union_chr${c}.bgen"
  plink2 --pfile "$ROOT/geno/union_chr${c}" --keep "$ROOT/keep/keep_union.txt" \
    --export A --threads "$THREADS" --memory "$PLINK_MEMORY_MB" --out "$ROOT/geno/union_chr${c}" > /dev/null
  plink2 --pfile "$ROOT/geno/union_chr${c}" --keep "$ROOT/keep/keep_fev.txt" \
    --pheno "$ROOT/pheno/pheno_copd.txt" --pheno-name copd \
    --covar "$COVAR" --covar-name age sex PC1-PC10 --covar-variance-standardize \
    --glm hide-covar firth-fallback --threads "$THREADS" --memory "$PLINK_MEMORY_MB" \
    --out "$ROOT/glm/fev_chr${c}" > /dev/null
  plink2 --pfile "$ROOT/geno/union_chr${c}" --keep "$ROOT/keep/keep_urate.txt" \
    --pheno "$ROOT/pheno/pheno_gout.txt" --pheno-name gout \
    --covar "$COVAR" --covar-name age sex PC1-PC10 --covar-variance-standardize \
    --glm hide-covar firth-fallback --threads "$THREADS" --memory "$PLINK_MEMORY_MB" \
    --out "$ROOT/glm/ura_chr${c}" > /dev/null
  echo "[aux_more2d] chr${c} done $(date -Is)"
done

"$PY" "$SCRIPT_DIR/run_more_traits.py" --mode analyze --root "$ROOT"
echo "[aux_more2d] ALL DONE $(date -Is)"
