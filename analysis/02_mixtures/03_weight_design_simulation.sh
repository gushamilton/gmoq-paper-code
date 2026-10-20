#!/bin/bash
# Weight-design simulation: model against hybrid weights for a membership effect
#
# What: runs the simulation of Supplementary Note 2. For each of ten scenarios (nine
#   two-component skew-normal mixtures and one shifted-lognormal control), 400 replicates are
#   simulated and passed through Stage 1 (03_simulate_stage1.py); the mixture is then fitted
#   to each replicate and Stage 2 is run under the model, hybrid, kernel-density-only and
#   oracle weight designs (03_stage2_designs.R); finally the replicates are summarised.
#   Replicates are processed in batches of 25. The batches are independent (each replicate
#   has its own seed) and are run in parallel on a cluster; here they run one after another.
# Paper: Supplementary Note 2 (Supplementary Note Table 1); Supplementary Fig. S7b
# Input: 03_weight_design_scenarios.json (scenarios, effects, seeds); synthetic data only
# Output: results/02_mixtures/weight_design/sim/*.csv
# Runs: on a cluster; about 4,000 CPU-minutes in total (12 GB memory per batch at N = 400,000)
set -euo pipefail
repo=${GMOQ_REPO:-$(git rev-parse --show-toplevel)}
here=$repo/analysis/02_mixtures
scenarios=$here/03_weight_design_scenarios.json
cfg() { Rscript -e "cat(yaml::read_yaml('$repo/config.yaml')\$$1)"; }
package_root=$(cfg 'package$stage1$root')     # pinned G-MoQ checkout (Stage 1 Python package)
python=$(cfg 'package$stage1$python')
work=${GMOQ_SIM_WORK:-$(cfg 'bp$private_root')/gmoq-paper-code/02_mixtures/weight_design_sim}
mkdir -p "$work/stage1" "$work/stage2"

for scenario in $(python3 -c "import json,sys; print(' '.join(json.load(open(sys.argv[1]))['scenarios']))" "$scenarios"); do
  for start in $(seq 0 25 375); do
    end=$((start + 25))
    [[ -e "$work/stage2/$(printf '%s_%04d.rds' "$scenario" "$start")" ]] && continue
    "$python" "$here/03_simulate_stage1.py" --config "$scenarios" --scenario "$scenario" \
      --rep-start "$start" --rep-end "$end" --out "$work/stage1" --package-root "$package_root"
    Rscript "$here/03_stage2_designs.R" --task --scenario "$scenario" \
      --rep-start "$start" --rep-end "$end" --in "$work/stage1" --out "$work/stage2"
  done
done
Rscript "$here/03_stage2_designs.R" --summary --in "$work/stage2"
