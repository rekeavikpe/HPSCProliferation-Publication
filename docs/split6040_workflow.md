# Random culture-level 60/40 comparison

The primary selection criterion is held-out negative log-likelihood (NLL), not held-out BIC. Each model's parameter set is selected using calibration data alone; those parameters are scored on the selection cultures without refitting. Training BIC is reported separately. Within 1e-6 selection NLL, prefer fewer fitted parameters, then lower training BIC. The rule and random seed are fixed before production fitting.

The original repository split seed, 20261003, assigns 19 of 31 recorded 100-mL cultures to calibration and 12 to selection. All days and states from one culture stay together. Conditions 1,2,3,5,6 have calibration/selection counts 2/2, 5/3, 4/2, 4/3 and 4/2. The saved split uses source condition IDs; these correspond to publication Conditions 1–5.

## Statistical scope

Condition–day–state means and SDs for calibration use only calibration cultures. Selection means are scored using the matching calibration SDs and configured floors, never selection SDs. Selection SDs remain descriptive experimental noise. Each partition uses its own measured day-zero glucose to initialize predictions, with a calibration fallback if unavailable. No post-day-zero selection outcomes are used to estimate parameters or noise weights. Fresh glucose retains the documented initial-glucose assumption.

The Gaussian score is a working diagonal likelihood on means, weighted by replicate SD rather than SEM. It does not model repeated-day covariance. All candidate kinetic parameters are free during selection; old all-data parameter estimates are not used as starts or fixed values. Bounds remain K_G ≤100 mM and K_Lac ≤100 mM. No 500-mL observations enter this workflow. Because all these data have already been inspected during model development, describe this as retrospective internal selection, not a newly untouched validation study.

## Server workflow

After the updated files are committed and pushed from the Desktop repository, use the server terminal:

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
mkdir -p logs
nohup env JULIA_NUM_THREADS=42 bash scripts/run_split6040.sh > logs/split6040.log 2>&1 < /dev/null &
echo $!
```

The script runs 100 independent population-search + local-refinement runs for each of the seven models, including baseline Model 0. Per run: population size 50, 20,000 population evaluations and up to 3,000 local iterations. Jobs run across Julia threads; BLAS uses one thread. Linux `flock` prevents two copies of this wrapper using the same output directory. `nohup` disconnects the job from the SSH session; follow the lab scheduler's policies if required.

```bash
tail -n 30 logs/split6040.log
find results/split6040_hybrid_lac100_kg100 -mindepth 3 -maxdepth 3 -name COMPLETE | wc -l
```

There are 700 per-run completion markers when fitting finishes. If terminated, repeat the same launch command with the same code and settings. Completed candidates are verified and reused; only incomplete runs are repeated. Changed source, data or settings are rejected on resume. Use a new result directory for a genuinely different analysis. Do not launch multiple direct Julia commands against the same directory.

## Outputs

Everything is under `results/split6040_hybrid_lac100_kg100/`:

- `split.csv`, `split_counts.csv`: reproducible culture assignment and counts.
- `protocol.toml`, `source_hashes.csv`: analysis settings and provenance; protocol budgets override the historical defaults in the copied calibration config.
- `calibration_observations.csv`, `selection_observations.csv`: separate means, SDs and scoring weights.
- `training_ranking.csv`: calibration NLL/BIC and stopping diagnostics; its model field retains storage IDs.
- `selection_ranking.csv`: primary predictive ranking with explicit publication and storage IDs.
- `selection_lock.csv`: selected model and training-parameter hash.
- `state_errors.csv`, `predictions.csv`: calibration/selection MRPEs, scaled errors and individual predicted means.
- `figures/split_model_comparison.*`: (a) training BIC distributions; (b) held-out NLL differences. PDF/EPS/SVG and 600-dpi PNG/JPEG/TIFF.
- `figures/supplementary_convergence.csv`: optimizer-status counts separate from the main figure.

Existing figure exporters and profile scripts that assume all-100-mL calibration must not be pointed at this split directory. This workflow intentionally does not refit on all data or evaluate external cohorts. Inspect the production ranking, then refit the chosen structure on all 100-mL cultures and derive any fixed values from that fit before profiles and external evaluations.

## Local checks

The all-seven-model short workflow uses one start, 100 population evaluations and 20 local iterations. These results only test execution, export and resumption; they are not manuscript evidence. The production comparison has not been run locally.

```bash
julia +1.12.4 --project=population population/test_split_selection.jl
```
