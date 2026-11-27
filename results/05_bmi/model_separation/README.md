# Equal-rank BMI model comparison

These saved aggregate results calibrate the OLS residual-score difference
`Q_B - Q_A` by Gaussian simulation at each lead's observed truncated mean effect,
using the chromosome-specific pooled covariance. `pA_obs` is its lower-tail
probability under a pure A profile, and `pB_obs` its upper-tail probability under
a pure B profile. A nominal 0.05 rule is applied per lead, without a
multiple-comparison adjustment.

For log BMI and `additive` versus `log_additive`, 570 of 639 leads reject the
raw-scale additive shift (89.2%); one lead rejects log-additivity. The per-lead
values are in `lead_pair_power.csv`; summaries across all leads and strength
bins are in `pair_power_summary.csv` and `pair_power_by_strength.csv`.
`q_crosscheck.csv` records the maximum relative difference from the saved
model scores for each model and scale (below 1e-12).

The original calibration used 20,000 draws and seed 20260928. The producing
numerical operations are in `analysis/05_bmi/03_model_separation.R`. It reads
only the released aggregate profiles, model weights and pooled covariance. Its
three equal-rank model pairs and draw order are preserved. The unrelated
component-genetic-correlation and nested-model calculations are excluded. No
chi-square reference is claimed for a nested unweighted-OLS score difference.

The 20-draw validation run checks data plumbing and the exact residual-score
cross-check; it is not the reported simulation calibration. A full 20,000-draw
rerun was unnecessary to recover the saved calibrated per-lead outputs.
