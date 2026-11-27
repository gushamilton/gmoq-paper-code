#!/usr/bin/env bash
# One paired cell per array task. Resource requests are supplied at submission
# from the pilot; the script does not submit jobs or choose resource padding.
#SBATCH --job-name=paired-speed
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
set -eo pipefail
case "$(hostname)" in bp1-compute*) ;; *) echo "requires a BP compute node" >&2; exit 2 ;; esac
for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    if [[ -f "$f" ]]; then source "$f"; break; fi
done
module load git/2.45.1-pqk5
module load languages/R/4.5.1
set -u
MODE=${1:?pilot or full}
OUT=${2:?output root}
SCRIPT_DIR=${PAIRED_RUNNER_DIR:?absolute directory containing paired_speed.py}
[[ -f "$SCRIPT_DIR/paired_speed.py" ]] || { echo "missing runner in PAIRED_RUNNER_DIR" >&2; exit 2; }
PAPER_ROOT=${PAPER_ROOT:-/user/work/fh6520/Gmoq_paper-speedprof-6bf1af5}
PIN_JSON=${PIN_JSON:-/user/work/fh6520/gmoq-paper-v1/Rlib/gmoq-166635d4/PIN.json}
PYTHON_BIN=${PYTHON_BIN:-/user/home/fh6520/miniforge3/envs/fungwas-stage1/bin/python}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export R_LIBS_USER=${R_LIBS_USER:-/user/work/fh6520/r_packages}
IDX=${SLURM_ARRAY_TASK_ID:?array task ID required}
NS=(10000 50000 400000)
case "$MODE" in
    pilot)
        [[ "$IDX" -ge 0 && "$IDX" -le 2 ]] || exit 2
        N=${NS[$IDX]}; REP=1; M=3
        ;;
    full)
        [[ "$IDX" -ge 0 && "$IDX" -le 8 ]] || exit 2
        N=${NS[$((IDX / 3))]}; REP=$((IDX % 3 + 1)); M=100
        ;;
    *) echo "mode must be pilot or full" >&2; exit 2 ;;
esac
ARGS=(cell --paper-root "$PAPER_ROOT" --pin-json "$PIN_JSON" --out "$OUT" --n "$N" --rep "$REP" --variants "$M")
if [[ -n "${METHODS:-}" ]]; then ARGS+=(--methods "$METHODS"); fi
if [[ -n "${TIMEOUT_SECONDS:-}" ]]; then ARGS+=(--timeout-seconds "$TIMEOUT_SECONDS"); fi
exec "$PYTHON_BIN" "$SCRIPT_DIR/paired_speed.py" "${ARGS[@]}"
