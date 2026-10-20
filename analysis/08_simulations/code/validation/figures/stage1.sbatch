#!/usr/bin/env bash
# Cluster job for simulation studies 1 and 2 (Stage 1 against ordinary RIF-OLS)
#
# What: loads the environment made by ../bootstrap.sh, pins every numerical library to one
#   thread and runs run_stage1.py, which simulates the four workloads, times the three
#   methods and exports the comparison tables.
# Paper: Supplementary Figs. S1, S2; Supplementary Table S2 (studies 1 and 2)
# Input: none (simulated data)
# Output: a new run directory under <ROOT>/.local/ (see run_stage1.py); the log in
#   <ROOT>/.local/stage1-figures-log-<attempt>/
# Runs: sbatch validation/figures/stage1.sbatch 05 400000   (about 25 minutes, one core)
#   The #SBATCH --chdir and --output lines hold the authors' cluster paths; edit them, and
#   the two user settings below, for your cluster.
#SBATCH --job-name=gmoq-s1s2
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --chdir=/user/work/fh6520/Gmoq_paper-figures-20260917
#SBATCH --output=/user/work/fh6520/Gmoq_paper-figures-20260917/.local/slurm-%j.out
set -Eeuo pipefail
: "${SLURM_JOB_ID:?Slurm required}"

## User settings
# Directory that contains validation/ (here: analysis/08_simulations/code).
ROOT=${GMOQ_SIM_ROOT:-/user/work/fh6520/Gmoq_paper-figures-20260917}
# Environment made by validation/bootstrap.sh.
BOOT=${GMOQ_SIM_ENV:-/user/work/fh6520/Gmoq_paper-validation-20260917/.local/validation-bootstrap-20260917-02}

ATTEMPT="${1:-01}"     # two-digit label of the run directory
N_LARGE="${2:-50000}"  # sample size of the two large workloads; the paper used 400000
[[ "$ATTEMPT" =~ ^[0-9]{2}$ ]] || exit 2
[[ "$N_LARGE" == 50000 || "$N_LARGE" == 400000 ]] || exit 2
cd "$ROOT"
LOG="$ROOT/.local/stage1-figures-log-$ATTEMPT"
mkdir "$LOG"
trap 'rc=$?; printf "%s\n" "$rc" > "$LOG/exit_code.txt"; exit "$rc"' EXIT
source "$BOOT/activate.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$GMOQ_STAGE1_PYTHON" validation/figures/run_stage1.py --attempt "$ATTEMPT" --n-large "$N_LARGE" > "$LOG/run.log" 2>&1
