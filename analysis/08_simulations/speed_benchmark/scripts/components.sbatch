#!/bin/bash
#SBATCH --job-name=sim-speed-profile
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:12:00
set -eo pipefail
case "$(hostname)" in bp1-compute*) ;; *) exit 2 ;; esac
source /etc/profile.d/modules.sh
module load git/2.45.1-pqk5 languages/R/4.5.1
set -u
export PAPER_REPO=/user/work/fh6520/Gmoq_paper-speedprof-6bf1af5
lib=/user/work/fh6520/gmoq-paper-v1/Rlib/gmoq-166635d4
py=/user/home/fh6520/miniforge3/envs/fungwas-stage1/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export R_PROFILE_USER=$(mktemp)
trap 'rm -f "$R_PROFILE_USER"' EXIT
echo ".libPaths(c('$lib', '/user/work/fh6520/r_packages', .libPaths()))" > "$R_PROFILE_USER"
run_dir=$PWD
mode=${1:?smoke or full}
case "$mode" in
  smoke) args=(--n 5000 --reps 1 --ml-snps 1) ;;
  full) args=(--n 10000,50000,400000 --reps 3 --ml-snps 10,10,5) ;;
  *) exit 2 ;;
esac
out="$run_dir/components-$mode-$SLURM_JOB_ID"
cd "$PAPER_REPO"
scripts=pipeline/aux/speed_profile
hostname
git rev-parse HEAD
"$py" "$scripts/speed_profile_gmoq.py" --out "$out" --pin-json "$lib/PIN.json" "${args[@]}"
Rscript "$scripts/stage2_profile.R" "$out" "$lib" gmoq166635d4 gmoq
"$py" - "$out" "$mode" <<'PY'
import csv, json, math, pathlib, sys
root = pathlib.Path(sys.argv[1])
with (root / 'python_timings.csv').open() as f:
    rows = list(csv.DictReader(f))
with (root / 'stage2_timings_gmoq166635d4.csv').open() as f:
    rrows = list(csv.DictReader(f))
ns = {5000} if sys.argv[2] == 'smoke' else {10000, 50000, 400000}
checks = [r for r in rows if r['component']=='check_max_rel_diff_binned_vs_dense']
assert {int(r['n']) for r in checks} == ns
assert all(float(r['seconds']) < 1e-4 for r in checks), checks
fits = [r for r in rows if r['component']=='ml_fit_one_snp']
assert len(fits) == (1 if sys.argv[2]=='smoke' else 25) and all(r['note']=='ok' for r in fits), fits
s2 = [r for r in rrows if r['component']=='r_param_gwas_bench_args' and float(r['M'])==100000]
assert {int(float(r['n'])) for r in s2} == ns and all(float(r['seconds'])>0 for r in s2), s2
pin = json.loads((root/'python_env.json').read_text())
assert pin['package_commit']=='166635d49527e27f8a714247d551c17093e02a55'
print('PASS',sys.argv[2],len(fits),'ML fits, Stage 1 checks', [r['seconds'] for r in checks], 'Stage 2 calibrated', sorted(ns))
PY
