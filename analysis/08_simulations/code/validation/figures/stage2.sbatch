#!/usr/bin/env bash
# Cluster job for simulation study 3 (Stage 2 recovery of simple simulated effects)
#
# What: loads the environment made by ../bootstrap.sh, runs simulate_stage2.py for four
#   scenarios x <SNPS> variants at sample size <N>, then draws working plots.
# Paper: Supplementary Fig. S3; Supplementary Table S2 (study 3)
# Input: none (simulated data)
# Output: <ROOT>/.local/stage2-recovery-n<N>-m<SNPS>/{results.csv,config.json,
#   summary.json,figures/}
# Runs: sbatch --time=16:00:00 validation/figures/stage2.sbatch 400000 1000
#   The paper's run (N = 400,000, 1,000 variants per scenario) took 14.5 hours on one core,
#   so it needs a longer --time than the default below. The #SBATCH --chdir and --output
#   lines hold the authors' cluster paths; edit them, and the user settings below.
#SBATCH --job-name=gmoq-stage2-recovery
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --chdir=/user/work/fh6520/Gmoq_paper-figures-20260917
#SBATCH --output=/user/work/fh6520/Gmoq_paper-figures-20260917/.local/stage2-%j.out
set -Eeuo pipefail
: "${SLURM_JOB_ID:?Slurm required}"

## User settings
# Directory that contains validation/ (here: analysis/08_simulations/code).
ROOT=${GMOQ_SIM_ROOT:-/user/work/fh6520/Gmoq_paper-figures-20260917}
# Environment made by validation/bootstrap.sh.
BOOT=${GMOQ_SIM_ENV:-/user/work/fh6520/Gmoq_paper-validation-20260917/.local/validation-bootstrap-20260917-02}

N="${1:-50000}"
SNPS="${2:-1000}"
cd "$ROOT"
OUT="$ROOT/.local/stage2-recovery-n${N}-m${SNPS}"
source "$BOOT/activate.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$GMOQ_STAGE1_PYTHON" validation/figures/simulate_stage2.py --n "$N" --snps "$SNPS" \
  --environment-manifest "$GMOQ_VALIDATION_ENV/environment.json" --out-dir "$OUT"
"$GMOQ_RSCRIPT_BIN" --vanilla validation/figures/plot_effect_recovery.R "$OUT/results.csv" "$OUT/figures"
