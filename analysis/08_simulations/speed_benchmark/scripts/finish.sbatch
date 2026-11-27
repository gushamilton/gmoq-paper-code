#!/bin/bash
#SBATCH --job-name=sim-speed-summary
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
set -euo pipefail
case "$(hostname)" in bp1-compute*) ;; *) exit 2 ;; esac
py=/user/home/fh6520/miniforge3/envs/fungwas-stage1/bin/python
components=${1:?completed component directory}
paired=${2:-$PWD/scripts/paired-speed}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$py" scripts/component_summary.py "$components" results/component_costs.csv
"$py" scripts/summarize_timings.py "$PWD" "$PWD/results"
"$py" "$paired/paired_speed.py" collate --out "$PWD/mixture-full" --reps 3 --variants 100
