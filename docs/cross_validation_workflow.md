# Condition holdout and partial 500 mL calibration

This workflow compares publication Model 0 (baseline) with publication Model 6
(storage model 5, legacy model 6). It uses the established Lagrangian P50 exposure,
likelihood, bounds, and hybrid optimizer. Every kinetic parameter is re-estimated
inside each fit. Fixed values and parameter starts from the full-data fit are not
imported. K_G and K_Lac retain the specified upper bounds of 100 mM.

## Experiments

1. Five leave-one-condition-out fits: use four 100 mL conditions to predict the
   omitted condition. Publication conditions 1–5 correspond to source conditions
   1, 2, 3, 5, 6. All days/states of each culture stay together.
2. Nine balanced partial-500 splits: fit all 31 small-vessel cultures plus four
   of the nine 500 mL cultures; predict the remaining five. A seeded random order
   (20261007) and rotating four-culture subsets give every 500 mL culture four
   training and five prediction appearances. These are nine of the 126 possible
   subsets, not exhaustive cross-validation. Folds overlap and are not independent.
3. A reference fit to all 100 mL cultures only. Evaluate that same fit on the
   exact held-out cultures in every partial-500 split for a paired comparison.

The default is 15 fit groups × two models × 100 starts = 3,000 optimizer runs.
Each run uses 20,000 population objective evaluations and up to 3,000 local
iterations. Completion is not a guarantee of convergence or a global optimum.

These are retrospective robustness experiments: the conditions and larger-vessel
results have already been examined during model development. Do not describe them
as prospective independent validation or use the outcomes to silently reselect
models. No additional profile-likelihood fits are part of this workflow.

## Scoring and leakage controls

Fitting means and SD-based weights are calculated from training cultures only.
Prediction NLL reuses training condition/day/state weights where available. For
an unseen condition or missing training time/state, its weight is the RMS of the
training sigma values for that state. This rule is fixed in advance; held-out
means and SDs do not set scoring weights. Report MRPE as the primary prediction
measure and NLL as a secondary diagnostic. Do not compare BIC across folds with
different training datasets. BIC is a within-fold comparison between the models.

Day-zero glucose from the prediction cultures supplies their measured initial
condition; it is not used for kinetic estimation. Fresh-medium glucose is assumed
equal to initial glucose, as in the existing workflow. If day-zero glucose is
missing, a training value from that condition is used; absent such a value, the
workflow stops rather than borrowing a held-out post-day-zero measurement.

State and overall MRPE use nonzero observed means after day zero. Experimental
noise uses SD/absolute mean only when at least two cultures were measured. Overall
metrics pool condition/day/state means; they are not an unweighted average of the
four state errors. Missing measurements are omitted.

## Mac: commit and push

```bash
cd ~/Desktop/HPSCProliferation-R2
git add population/CrossValidation.jl population/run_cross_validation.jl \
  population/test_cross_validation.jl figures/render_cross_validation.jl \
  scripts/run_cross_validation.sh docs/cross_validation_workflow.md
git diff --cached --stat
git commit -m "Add condition holdout and partial 500 mL cross-validation"
git push origin main
```

These explicit paths exclude unrelated local changes.

## Server: launch

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
export PATH="$HOME/.juliaup/bin:$PATH"
julia +1.12.4 --version
mkdir -p logs
nohup env JULIA_NUM_THREADS=42 HPSC_CV_STARTS=100 \
  bash scripts/run_cross_validation.sh \
  > logs/cross_validation.log 2>&1 < /dev/null &
echo $! > logs/cross_validation.pid
```

Use this on an allocated compute session if required by your lab's server policy.
Do not run another 42-thread analysis at the same time unless resources permit it.

## Monitor and resume

```bash
tail -n 30 logs/cross_validation.log
pgrep -af '[j]ulia.*run_cross_validation'
find results/cross_validation_balanced9_lac100_kg100 \
  -mindepth 4 -maxdepth 4 -name COMPLETE | wc -l
```

The last count targets 3,000 for the default run. Each verified optimizer run is
saved separately. Repeating the launch command with unchanged settings resumes
incomplete runs. Use `>> logs/cross_validation.log` instead of `>` on a restart to
preserve earlier log text. A per-directory lock prevents concurrent launches.
Source/data hashes and settings are checked before reuse; if they changed, use a
new output directory. Do not pull code updates while this workflow is running.

## Outputs

`results/cross_validation_balanced9_lac100_kg100/`:

- `folds.csv`: planned folds and culture counts.
- `<fold>/split.csv`: individual culture assignments.
- `<fold>/held_weight_sources.csv`: prediction weights and their origin.
- `<fold>/model*/run*/COMPLETE`: verified optimizer checkpoints.
- `fit_summary.csv`: training BIC, NLL, boundaries (in per-model fit tables),
  prediction NLL, and selected-fit stopping status.
- `state_metrics.csv`: state-specific and overall MRPE and experimental noise,
  including matched 100 mL-only reference predictions.
- `predictions.csv`: measured means, SDs, sample counts, weights and predictions.
- `figures/leave_condition_out.*` and `figures/partial500_comparison.*`:
  PDF, EPS, PNG, JPEG, SVG, TIFF (raster exports at 600 dpi).
- `COMPLETE`: written after numerical evaluation; the shell reports completion
  after rendering. If plotting fails, rerun the renderer alone.

```bash
julia +1.12.4 --project=figures figures/render_cross_validation.jl \
  results/cross_validation_balanced9_lac100_kg100
```

## Optional exhaustive design

`HPSC_CV_DESIGN=all126` uses all 126 four-culture subsets: 132 groups including
LOCO and reference, or 26,400 optimizer runs. Give this a distinct directory:

```bash
HPSC_CV_DESIGN=all126 bash scripts/run_cross_validation.sh \
  results/cross_validation_all126_lac100_kg100
```

## Local workflow checks

```bash
julia +1.12.4 --project=population population/test_cross_validation.jl
```

A reduced-budget smoke test checks execution and file generation only; its fits
must not be used as scientific results. Use a separate output directory and
`HPSC_CV_STARTS=1 HPSC_CV_EVALUATIONS=100 HPSC_CV_LOCAL_ITERATIONS=2`.

Local verification on 7 October 2026: 119 split/weight/metric checks passed;
a one-start reduced-budget run completed all 15 groups for both models; a second
invocation reused all 30 checkpoints; 390 metric rows were generated with finite
MRPE values; both figure layouts exported successfully in all six formats.
The local workflow-check results are not publication fits.
