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

echo "[aux_hidbio] start $(date -Is) host=$(hostname) job=$SLURM_JOB_ID"
"$PY" "$SCRIPT_DIR/run_hidden_biology.py" --mode prep --root "$ROOT"

set +u
source /etc/profile >/dev/null 2>&1 || true
module load apps/plink2/2.00a68LM
set -u
mkdir -p "$ROOT/geno" "$ROOT/glm"

# Extract ALT leads and run liver-disease GLMs + dosage export
for chr in $(seq 1 22); do
  c=$(printf %02d "$chr")
  L="$ROOT/snps/alt_chr${c}.txt"
  if [[ -s "$L" ]]; then
    "$BGENIX" -g "$BG/data.chr${c}.bgen" -incl-rsids $(tr '\n' ' ' < "$L") > "$ROOT/geno/alt_chr${c}.bgen"
    plink2 --bgen "$ROOT/geno/alt_chr${c}.bgen" ref-first --sample "$SAMPLE" \
      --keep "$ROOT/keep/keep_alt.txt" --make-pgen --threads "$THREADS" --memory "$PLINK_MEMORY_MB" \
      --out "$ROOT/geno/alt_chr${c}" > /dev/null
    rm -f "$ROOT/geno/alt_chr${c}.bgen"
    plink2 --pfile "$ROOT/geno/alt_chr${c}" --keep "$ROOT/keep/keep_alt.txt" \
      --pheno "$ROOT/pheno/pheno_alt.txt" --pheno-name masld_k760 fo_k76 cirrhosis masld_any masld_nonalc \
      --covar "$COVAR" --covar-name age sex PC1-PC10 --covar-variance-standardize \
      --glm hide-covar firth-fallback cols=+a1freq --threads "$THREADS" --memory "$PLINK_MEMORY_MB" \
      --out "$ROOT/glm/alt_chr${c}" > /dev/null
    plink2 --pfile "$ROOT/geno/alt_chr${c}" --keep "$ROOT/keep/keep_alt.txt" \
      --export A --threads "$THREADS" --memory "$PLINK_MEMORY_MB" --out "$ROOT/geno/alt_chr${c}" > /dev/null
    echo "[aux_hidbio] ALT chr${c} done $(date -Is)"
  fi

  L2="$ROOT/snps/ldl_chr${c}.txt"
  if [[ -s "$L2" ]]; then
    "$BGENIX" -g "$BG/data.chr${c}.bgen" -incl-rsids $(tr '\n' ' ' < "$L2") > "$ROOT/geno/ldl_chr${c}.bgen"
    plink2 --bgen "$ROOT/geno/ldl_chr${c}.bgen" ref-first --sample "$SAMPLE" \
      --keep "$ROOT/keep/keep_ldl.txt" --make-pgen --threads "$THREADS" --memory "$PLINK_MEMORY_MB" \
      --out "$ROOT/geno/ldl_chr${c}" > /dev/null
    rm -f "$ROOT/geno/ldl_chr${c}.bgen"
    plink2 --pfile "$ROOT/geno/ldl_chr${c}" --keep "$ROOT/keep/keep_ldl.txt" \
      --pheno "$ROOT/pheno/pheno_ldl.txt" --pheno-name statin \
      --covar "$COVAR" --covar-name age sex PC1-PC10 --covar-variance-standardize \
      --glm hide-covar firth-fallback cols=+a1freq --threads "$THREADS" --memory "$PLINK_MEMORY_MB" \
      --out "$ROOT/glm/ldl_chr${c}" > /dev/null
    plink2 --pfile "$ROOT/geno/ldl_chr${c}" --keep "$ROOT/keep/keep_ldl.txt" \
      --export A --threads "$THREADS" --memory "$PLINK_MEMORY_MB" --out "$ROOT/geno/ldl_chr${c}" > /dev/null
    echo "[aux_hidbio] LDL chr${c} done $(date -Is)"
  fi
done

"$PY" "$SCRIPT_DIR/run_hidden_biology.py" --mode analyze --root "$ROOT"
echo "[aux_hidbio] ALL DONE $(date -Is)"
