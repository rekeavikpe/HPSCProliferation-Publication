# Pooled replicate cross-validation and exclusion sensitivity

## Design

Only ESI-017 is pooled. There are 46 independent bioreactor replicates: 31 in the original five 100 mL conditions, nine at 500 mL, and three for each of the 20% and 50% media-exchange protocols. Cell lines are not pooled with ESI-017.

Three folds, seed 20261014, keep complete replicate histories together. Each of eight groups is represented in every held-out fold. Three folds were chosen because each media-exchange group has only three replicates. Both Model 0 and publication Model 6 (storage model 5) use identical folds. All active parameters are fitted using calibration data only. The original condition-mean Gaussian likelihood, floors, bounds and hybrid optimizer are retained. Media-exchange schedules and day-zero assumptions match the existing media-exchange evaluation. Day-zero measurements specify initial conditions, not fitted outcomes. Prediction error excludes day zero and zero denominators. Replicates lacking measurements for a state do not contribute to that state's error summary. The effective replicate count is reported.

For every fold, compare a full calibration fit with a trimmed calibration fit against EXACTLY the same held-out replicates. Rank calibration replicates within each experimental group by final-day measured viable cell density. Remove one lowest and one highest, breaking ties by series ID, only when at least four calibration replicates and all terminal cell counts exist. Thus at least two calibration replicates remain. No held-out outcome determines exclusions. Media-exchange calibration groups have two replicates and are never trimmed. This is an exclusion sensitivity analysis, not a declaration that extreme measurements are invalid.

Also fit full and trimmed original 100 mL datasets and predict all nine 500 mL replicates without removing any validation observations. Existing completed fits are preserved. This additional reference is rerun under the same locked optimizer settings as the sensitivity analysis.

## Computation

Eight fits per model (three full folds, three trimmed folds, full-100 reference and trimmed-100 reference), two models and 100 hybrid starts = 1,600 starts. Default population budget is 20,000 evaluations followed by 3,000 local iterations. Parallelism is across starts, with 42 threads by default. Every completed start is checked and saved. Run the same command again after interruption to resume. Do not change source files or budgets within an existing output directory.

```bash
export PATH="$HOME/.juliaup/bin:$PATH"
mkdir -p logs
nohup env JULIA_NUM_THREADS=42 bash scripts/run_pooled_validation.sh \
  > logs/pooled_validation.log 2>&1 < /dev/null &
```

Progress:
```bash
tail -n 30 logs/pooled_validation.log
find results/pooled_validation3_lac100_kg100 -mindepth 4 -maxdepth 4 -name COMPLETE | wc -l
```

Preparation only (no fitting):
```bash
HPSC_CV_PREPARE_ONLY=1 bash scripts/run_pooled_validation.sh
```

If preparation used reduced budgets, use a DIFFERENT output directory for the final run. Optional smoke budgets: HPSC_CV_STARTS=1 HPSC_CV_EVALUATIONS=100 HPSC_CV_LOCAL_ITERATIONS=5. These do not produce scientific results. HPSC_JULIA can specify an executable without juliaup.

## Results

`results/pooled_validation3_lac100_kg100/` contains:

- `folds.csv`, per-fit `split.csv`: exact included/held-out replicate identities.
- `*_trimming_decisions.csv`: rankings, terminal counts and exclusions or reasons for retaining a group.
- `replicate_prediction_errors.csv`: MRPE for each excluded replicate and state, plus overall.
- `replicate_error_summary.csv`: mean, sample SD and sample variance of replicate MRPEs, by experiment/model/group/state and pooled across groups. The pooled statistic weights replicates equally. Overall within a replicate pools eligible measured state/day points.
- `fold_error_summary.csv`: mean and variance across fold means. These are descriptive, not confidence intervals, because calibration sets overlap.
- `replicate_predictions.csv`: auditable measured and predicted values.
- `state_metrics.csv`: companion errors against condition/day means, for continuity with previous reporting. These are NOT the same estimand as errors against individual replicates.
- Per-fit parameter tables, optimizer diagnostics and saved stages.

Replicate prediction uses that replicate's measured day-zero values where the established protocol does so. Day-zero fallback, when needed, uses calibration values. No held-out post-day-zero measurements enter calibration or noise weighting. The same measured starting conditions are used for paired full/trimmed predictions.

This pooled analysis evaluates prediction of excluded replicates among represented experimental conditions. It does not replace the separate 100 mL-to-500 mL assessment, or establish prediction of an entirely new condition. Excluded extreme calibration replicates remain in the raw data and are held out normally in their own fold.
