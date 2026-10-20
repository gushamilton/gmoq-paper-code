#!/bin/bash
# Cost of G-MoQ and of direct mixture maximum likelihood: run all three parts
#
# What: runs the timing study behind Fig. 2c and Supplementary Fig. S5d on one core of a
#   cluster node (all numerical libraries limited to one thread): part 1 times Stage 1 and
#   direct mixture maximum likelihood on the benchmark's simulated panels; part 2 times Stage 2
#   alone with the Stage 2 package the benchmark used (labelled pin9c47e7b) and with the one
#   used for the paper's analyses (Stage 2 commit 9fee67d); part 3 computes
#   the cost per variant. Simulated data only.
# Paper: Fig. 2c; Supplementary Fig. S5d
# Input: the mixture benchmark code and its validated software environment (their locations
#   follow from config.yaml: p08.timings and bp.env_bootstrap)
# Output: results/03_testosterone/speed/{python_component_timings,stage2_timings,per_snp_costs}.csv
# Runs: on a cluster node, single-threaded, in about 1.5 hours (24 GB memory)
set -euo pipefail
repo=${GMOQ_REPO:-$(git rev-parse --show-toplevel)}
here=$repo/analysis/03_testosterone
cfg() { Rscript -e "cat(yaml::read_yaml('$repo/config.yaml')\$$1)"; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
benchmark_root=$(dirname "$(dirname "$(dirname "$(cfg 'p08$timings')")")")   # checkout holding validation/
environment=$(dirname "$(cfg 'bp$env_bootstrap')")                               # activate.sh, environment.json, Rlib/
python=$(cfg 'package$stage1$python')
run=$(cfg 'bp$private_root')/gmoq-paper-code/03_testosterone/speed_$(date +%Y%m%d%H%M%S)

"$python" "$here/06_speed_profile.py" --out "$run" --benchmark-root "$benchmark_root" \
  --environment "$environment" --rscript "$(command -v Rscript)"
Rscript "$here/06_speed_stage2.R" "$run" "$environment/Rlib" pin9c47e7b
Rscript "$here/06_speed_stage2.R" "$run" "$(cfg 'package$stage2$lib')" pin9fee67d
Rscript "$here/06_speed_costs.R" "$run"
