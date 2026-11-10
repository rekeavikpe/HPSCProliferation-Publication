# Cross-validation error tables

Run from the repository root after the balanced-nine cross-validation has finished:

```bash
python3 scripts/report_cross_validation.py results/cross_validation_balanced9_lac100_kg100
```

This reads saved `state_metrics.csv`, `predictions.csv`, `fit_summary.csv`, and `folds.csv`. It does not run Julia, refit parameters, or change saved fits. All 390 state/overall MRPE entries are independently checked against saved prediction rows before reporting.

Outputs under `results/cross_validation_balanced9_lac100_kg100/tables/`:
- `index.html`: readable tables for Models 6 and 0, calibration and prediction, every fold, with the 100 mL-only reference and fold averages.
- `errors_by_run.csv`: per-run state-specific and overall MRPE, experimental noise, eligible observation counts, training BIC, and stopping status.
- `average_errors_by_design.csv`: equal-fold averages and observation-count-weighted averages, explicitly distinguished.
- `provenance.json`: input hashes and metric definitions.

MRPE is the mean absolute relative error against eligible condition-day means (day > 0, nonzero observed mean). Overall pools all eligible state-day observations. It is not a signed bias or a mean across optimizer starts. The 100 optimization starts provide candidate fits; these tables evaluate the selected fit from each calibration run/fold.

LOCO predicts a held-out 100 mL condition. The reference predicts all nine 500 mL cultures. These assess different prediction tasks. In partial-500 evaluation, compare the partial-500 fit with the 100 mL-only reference on exactly the same five held-out cultures for each split. The nine overlapping splits are not independent validation experiments; do not interpret the fold variation as an independent-replicate confidence interval.

Failed stopping flags remain visible. The reporting script expects the balanced-nine design and fails if files are incomplete or saved MRPEs disagree with their predictions.
