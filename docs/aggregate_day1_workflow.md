# Model 6 with aggregate size initialized at day 1

This is David's additional P50-shear/P50-EDR fit. It is separate from the running
P50-shear/P75-EDR fit and does not change that job or its checkpoints.

## What changes

- Keep the original 19 training and 12 selection cultures, Model 6 equations,
  all 11 fitted parameters, parameter bounds, P50 exposure functions and feeds.
- Cells, glucose and lactate begin at their original day 0 values. Only those
  three equations are integrated before day 1. Their size-dependent terms use
  the existing 15 µm input held fixed. The aggregate-size equation is inactive.
  This fixed input is a modeling assumption, not a measured aggregate diameter
  or a prediction of first-day aggregation.
- The CFD size grid starts at 200 µm. The existing nearest-endpoint rule is
  retained, so the 15 µm input uses the 200 µm P50 exposure values. No single-cell
  CFD values are inferred. The aggregate-size growth limitation is retained
  with 15 µm as its input before day 1.
- At exactly 24 h, start the aggregate equation with the measured day 1 mean
  of the cultures in that cohort. Carry the predicted cell density, glucose
  and lactate forward without resetting them to measurements. Preserve the
  absolute culture age in the death term and the original feeding schedule.
  No aggregate-size prediction is made or plotted before day 1.
- Training initialization uses only training cultures. Selection initialization
  uses only selection cultures. External experiments use their own day 1 size.
  There is no parameter fitting on selection or validation data.
- Exclude day 1 aggregate measurements from the fit and comparison scores because
  they supply initial conditions. Keep all previously scored days for the other
  three states. Score both the original model and the new fit on the same mask.
  Compare the new selection BIC with this matched original score, not directly
  with the earlier 1458.73 score, which included day 1 aggregate observations.
- Evaluate training, selection, all 100 mL data descriptively, 500 mL,
  20%/50% media exchange and direct transfer to the two iPSC lines. External
  predictions are conditional on the measured day 1 aggregate size.

The earlier October 20 day 1 analysis initialized all four states at day 1.
This new fit does not do that. It does not establish that aggregates are absent
or that all first-day biology is captured by the fixed-size approximation.

## Search and output

Run 100 independent hybrid starts, each with 20,000 population evaluations and
up to 3,000 local iterations. Select the optimum by training likelihood.
Use the same existing fixed Gaussian noise weights for retained observations.
Selection weights remain training-derived. Overall MRPE pools eligible
condition-day-state observations, following the existing implementation.

Output: `results/aggregate_day1_fixed15_p50shear_p50edr`.

This is a new output root and a new analysis identifier. Earlier aggregate-D1
checkpoints that evolved aggregate size before D1 are incompatible and must not
be reused. They remain untouched in their original folder. A manifest check
rejects reuse under changed equations or settings.

- `aggregate_day1_fit/selected_parameters.csv`: selected training estimate.
- `selection_scores.csv`: matched-observation selection NLL and BIC for both arms.
- `comparison.csv`: state-specific and overall MRPE with differences.
- `predictions.csv`: predicted and measured values, with an explicit scoring role.
- `curves.csv`: trajectories for plotting. The new aggregate curve begins at D1.
- `initial_states.csv`: all four states immediately before/after D1 initialization.
- `reference_original_scores.csv`: unchanged original scores for reference only.
- `index.html`: comparison and assumptions.
- `figures/`: PNG/PDF comparison plots for each subset after the optional rendering step.

## Mac

```bash
cd "$HOME/Desktop/PhD Work/Code/HPSCProliferation-R2"
git add revision/AggregateDayOne.jl revision/run_aggregate_day1.jl \
  scripts/run_aggregate_day1.sh tests/test_aggregate_day1.jl \
  figures/render_aggregate_day1.jl docs/aggregate_day1_workflow.md
git diff --cached --stat
git commit -m "Hold size input at 15 um before day 1 and start aggregate ODE at day 1"
git push origin main
```

Only these new files are needed. Existing modified figure scripts are unrelated
and should not be included in this commit.

## Server in parallel

Do not stop or restart the mixed-exposure job. This script defaults to eight
threads while that job is running. Increase the new job's thread count only if
the server has spare allocated cores. Two concurrent jobs share resources.

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
mkdir -p logs
nohup env JULIA_NUM_THREADS=8 HPSC_AGGDAY1_STARTS=100 \
  HPSC_AGGDAY1_EVALUATIONS=20000 HPSC_AGGDAY1_ITERATIONS=3000 \
  HPSC_AGGDAY1_OUTPUT=results/aggregate_day1_fixed15_p50shear_p50edr \
  bash scripts/run_aggregate_day1.sh \
  >> logs/aggregate_day1_fixed15.log 2>&1 < /dev/null &
echo $!
```

The script validates the design, fits and evaluates. It saves trajectories for
local plotting with Arial. Set `HPSC_AGGDAY1_RENDER=1` to also render on the server
if Arial is installed there. Otherwise download the results and render locally:

```bash
julia +1.12.4 --project=figures figures/render_aggregate_day1.jl PATH_TO_DOWNLOADED_RESULTS
```

The existing Julia 1.12.4 population and figures environments are reused. Its output lock is
independent of the mixed-exposure lock. Only one job may write this new output.

## Progress and resume

```bash
cat results/aggregate_day1_fixed15_p50shear_p50edr/progress.log
tail -n 30 logs/aggregate_day1_fixed15.log
find results/aggregate_day1_fixed15_p50shear_p50edr/aggregate_day1_fit \
  -mindepth 2 -maxdepth 2 -type f -name COMPLETE | wc -l
pgrep -af '[j]ulia.*run_aggregate_day1'
```

The count is completed starts, not the largest start identifier in the log.
After all starts finish, evaluation writes `COMPLETE`. Optional plotting writes
`FIGURES_COMPLETE`. A same-settings restart after this job exits verifies and
reuses completed starts. It does not restart the mixed-exposure job.

## Collect

```bash
cd ~/HPSCProliferation-R2
test -f results/aggregate_day1_fixed15_p50shear_p50edr/COMPLETE && \
tar -czf aggregate_day1_fixed15_for_review.tar.gz \
  --exclude='run[0-9]*' \
  results/aggregate_day1_fixed15_p50shear_p50edr logs/aggregate_day1_fixed15.log
ls -lh aggregate_day1_fixed15_for_review.tar.gz
```

Download that archive for review. Neither this run nor the mixed-exposure run
automatically changes the submission documents.
