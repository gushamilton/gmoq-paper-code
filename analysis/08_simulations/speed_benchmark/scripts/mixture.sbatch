#!/bin/bash
#SBATCH --job-name=sim-speed-mixture
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
set -euo pipefail
export PAIRED_RUNNER_DIR="${PAIRED_RUNNER_DIR:-$PWD/scripts/paired-speed}"
export R_LIBS_USER=/user/work/fh6520/r_packages
export R_PROFILE_USER=$(mktemp)
trap 'rm -f "$R_PROFILE_USER"' EXIT
echo ".libPaths(c('/user/work/fh6520/gmoq-paper-v1/Rlib/gmoq-166635d4','/user/work/fh6520/r_packages',.libPaths()))" > "$R_PROFILE_USER"
# Invoke by the actual path so the runner can locate its Python and R sources.
bash "$PAIRED_RUNNER_DIR/paired_speed.sbatch" "$1" "${2:-$PWD/mixture-$1}"
