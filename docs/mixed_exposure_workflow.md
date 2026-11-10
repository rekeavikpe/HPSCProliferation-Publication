# Model 6 with P50 shear and P75 EDR

This run tests David's requested combination of Lagrangian P50 shear stress and
Lagrangian P75 EDR against the current Model 6 with P50 for both quantities.
It does not change the manuscript or replace existing results automatically.

## Analysis

- Reuse the current P50/P50 Model 6 training estimate and verify its parameters,
  training likelihood and original replicate assignments against the saved split.
- Fit the P50/P75 version to the same 19 training cultures with all 11 Model 6
  parameters free, the same bounds and the same fixed Gaussian noise weights.
- Run 100 independent hybrid searches, each with 20,000 population evaluations
  followed by up to 3,000 local iterations. Choose the saved optimum by training
  likelihood, never by selection or validation error.
- Apply the existing size interpolation and endpoint rules separately to P50
  shear and P75 EDR. Both inputs follow the predicted aggregate diameter.
- Evaluate selection BIC and NLL on the same 12 selection cultures using the
  training-derived noise weights. No new noise parameters are fitted.
- Compare each state's MRPE and overall MRPE on training, selection, all 100 mL
  data (descriptive only), 500 mL validation and 20%/50% media-exchange experiments.
- Include direct predictions of both iPSC lines using the ESI-017 parameters.
  These are not new cell-line recalibrations. The existing interpolation between
  40 and 60 RPM is preserved for the 50 RPM iPSC experiment.
- Initial conditions, feeding schedules, eligible observations and MRPE
  aggregation match the current implementation. Overall MRPE pools eligible
  condition-day-state observations. It is not an average over optimization runs.

Only the mixed-exposure Model 6 requires fitting. Other models, profiles and
cross-validation folds are not rerun.

## Mac: commit only the new analysis files

```bash
cd "$HOME/Desktop/PhD Work/Code/HPSCProliferation-R2"
git add revision/MixedExposure.jl revision/run_mixed_exposure.jl \
  scripts/run_mixed_exposure.sh tests/test_mixed_exposure.jl \
  docs/mixed_exposure_workflow.md
git diff --cached --stat
git commit -m "Compare Model 6 with P50 shear and P75 EDR"
git push origin main
```

## Server: start the comparison

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
mkdir -p logs
nohup env JULIA_NUM_THREADS=42 HPSC_MIXED_STARTS=100 \
  bash scripts/run_mixed_exposure.sh \
  >> logs/mixed_exposure.log 2>&1 < /dev/null &
echo $!
```

The script uses Julia 1.12.4. It verifies the baseline scores and data allocation
before fitting. Results go into a new directory and existing analyses are not
modified. Do not launch two copies. On Linux the output lock prevents this.

## Progress and resume

```bash
tail -n 30 logs/mixed_exposure.log
cat results/mixed_exposure_p50shear_p75edr/progress.log
find results/mixed_exposure_p50shear_p75edr/mixed_fit \
  -mindepth 2 -maxdepth 2 -name COMPLETE | wc -l
pgrep -af '[j]ulia.*run_mixed_exposure'
```

The count reaches 100. After that the script evaluates all data subsets and
writes the top-level `COMPLETE` file. A silent interval can occur while each
hybrid search is still optimizing. Repeating the same launch after the process
has exited resumes verified completed starts. Keep the same settings, source
files and output path. A changed source or budget is rejected rather than mixed
with earlier checkpoints.

## Results and download

- `index.html`: selection BIC and MRPE comparison with differences.
- `selection_scores.csv`: selection NLL, BIC, parameter count and sample count.
- `comparison.csv`: state-specific and overall MRPE, including each condition.
- `reference_predictions.csv` and `mixed_predictions.csv`: paired predictions
  and measurements for checking errors.
- `reference_curves.csv` and `mixed_curves.csv`: saved trajectories for plotting.
- `mixed_fit/selected_parameters.csv`: training-selected mixed-exposure estimate.
- `manifest.toml` and `split.csv`: settings, source hashes and culture allocation.

Negative MRPE differences favour mixed exposure. Check every requested subset
and state rather than describing an improvement in one quantity as universal.

```bash
cd ~/HPSCProliferation-R2
test -f results/mixed_exposure_p50shear_p75edr/COMPLETE && \
tar -czf mixed_exposure_for_review.tar.gz \
  --exclude='run[0-9]*' \
  results/mixed_exposure_p50shear_p75edr logs/mixed_exposure.log
ls -lh mixed_exposure_for_review.tar.gz
```

Download that archive from the server for interpretation. The excluded search
checkpoints remain on the server. The archive includes the selected fit and
comparison data.
