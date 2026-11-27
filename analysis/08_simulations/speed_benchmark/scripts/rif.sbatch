#!/bin/bash
#SBATCH --job-name=sim-speed-rif
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
set -eo pipefail
case "$(hostname)" in bp1-compute*) ;; *) exit 2 ;; esac
source /etc/profile.d/modules.sh
module load git/2.45.1-pqk5
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py=/user/home/fh6520/miniforge3/envs/fungwas-stage1/bin/python
pin=/user/work/fh6520/gmoq-paper-v1/Rlib/gmoq-166635d4/PIN.json
case "$1" in
  smoke) args=(--preset smoke) ;;
  n400000_m1|n400000_m1000) args=(--workload "$1") ;;
  *) exit 2 ;;
esac
exec "$py" scripts/rif-speed/run_stage1.py --pin-json "$pin" --out "$PWD/rif-$1-$SLURM_JOB_ID" "${args[@]}"
