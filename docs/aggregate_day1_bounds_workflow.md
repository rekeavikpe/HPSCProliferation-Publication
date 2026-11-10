# Aggregate day 1: targeted K_G and K_eps boundary checks

This is a separate diagnostic analysis. It does not replace the completed 100-start fit, modify either existing exposure run or edit any manuscript, response or supplementary document.

## Design

- Reuse the completed aggregate-day-1 fit and eight distinct saved training solutions ranked by training likelihood. No new 100-start population search is requested.
- Preserve the 19 training cultures, 12 selection cultures, Model 6 equations, P50/P50 inputs, feeds, likelihood scales and all other parameter bounds. The fixed 15 µm input before day 1 and measured cohort-specific day 1 aggregate initialization remain unchanged.
- Reoptimize all 11 parameters in three nested search boxes. K_G upper limits are 100, 1,000 and 10,000 mM. K_eps upper limits are 5, 50 and 500 m²/s³. The 10-fold and 100-fold extensions test sensitivity to numerical search limits. They are not literature-derived biological ranges or proposed final bounds. Other limits, including mu_max, stay unchanged and are reported if reached.
- Use bounded least-squares refinement followed by a bounded Nelder–Mead check in normalized log coordinates. Each original seed is retained if refinement is worse. Additional deterministic seeds test the newly opened ranges. Actual optimizer status is recorded, including non-convergence.
- At a sparse set of fixed K_G values and fixed K_eps values, optimize the other ten parameters using up to four distinct saved/continuation starts. Include the fitted value and original upper limits. These are training-likelihood diagnostics, not a full confidence-interval calculation. The other target parameter remains free.
- For K_G scan initialization, also try a mu_max starting value that preserves the seed's glucose growth factor at the median initial training glucose. This supplies feasible starting trajectories when fixing K_G far below its fitted value. It is only an initial guess. mu_max remains free, its original bounds are retained and all ten nuisance parameters are optimized. Invalid starting trajectories are screened before optimization.
- Refine the lowest training-likelihood candidate found by these checks without fixing either target. Choose using training data only. Evaluate every bound case and the best diagnostic candidate on the same eight datasets. Do not refit external data. iPSC results remain direct transfers, not line-specific recalibrations.
- Score identical observations in all comparisons. Day 1 aggregate observations initialize the new model and are excluded from both models' scores. Day 1 observations of the other states remain scored. Selection has 111 scored observations and 11 fitted parameters.

The original model is still the day-0 prediction. Matching scored observations does not remove the additional measured day-1 information supplied to the D1 model. Interpret the new aggregate predictions as predictions after measured D1 initialization.

## Preservation and resume

Input: `results/aggregate_day1_fixed15_p50shear_p50edr`.

Output: `results/aggregate_day1_boundary_checks`.

The launcher rejects overlapping paths, including symbolic links, before writing output. It does not delete, move or write to any completed fit. All files in the original D1 result tree are hashed at startup and verified again at completion. The original calculation's source/data hashes and likelihood are checked. A separate manifest locks the diagnostic settings and inputs. Each optimizer start has a checkpoint. Restarting the same command reuses and numerically verifies completed checkpoints. Changed inputs/settings cause an error rather than silent reuse. On Linux a file lock prevents two jobs from writing this output concurrently.

## Run on the server

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
mkdir -p logs
nohup env JULIA_NUM_THREADS=8 \
  bash scripts/run_aggregate_day1_bounds.sh \
  >> logs/aggregate_day1_boundary_checks.log 2>&1 < /dev/null &
echo $!
```

There is no need to stop or rerun either completed fit. Do not change the diagnostic source or budgets while this job is running. Restart the same command only after the previous process has exited. Do not delete its output directory.

Optional settings before the first launch: `HPSC_AGGBOUND_SEEDS=8`, `HPSC_AGGBOUND_SCAN_SEEDS=4`, `HPSC_AGGBOUND_LM_ITERATIONS=200`, `HPSC_AGGBOUND_POLISH_ITERATIONS=3000`. Paths can be overridden with `HPSC_AGGBOUND_REFERENCE`, `HPSC_AGGBOUND_SPLIT`, `HPSC_AGGBOUND_SOURCE` and `HPSC_AGGBOUND_OUTPUT`. The default settings should be used for this comparison.

## Progress

```bash
cd ~/HPSCProliferation-R2
tail -n 25 logs/aggregate_day1_boundary_checks.log
tail -n 10 results/aggregate_day1_boundary_checks/progress.log
pgrep -af '[j]ulia.*run_aggregate_day1_bounds'
```

The log reports each completed optimizer start, each bound case and each fixed-value scan point. It subsequently evaluates all datasets and writes a final `COMPLETE` marker. Completion means the calculations finished, not that the parameters have been established as identifiable.

## Results to review

- `audit.csv`: likelihood improvement, convergence flags, all active bounds and verification that the original result tree is unchanged.
- `training_candidates.csv`: every free-bound candidate and scan candidate, with training NLL and parameter values.
- `likelihood_scans.csv`: fixed values, optimized training NLL, change from the best candidate and nuisance parameters at bounds. A sparse scan does not establish an exact confidence limit. A non-converged point must not be interpreted as evidence of non-identifiability.
- `selection_scores.csv` and `state_errors.csv`: selection BIC and state-specific MRPE for all eight datasets. These results never determine the fitting seed or optimization outcome.
- `predictions.csv`, `curves.csv` and `initial_states.csv`: evaluated points, plot-ready curves and day-1 initialization records.
- `index.html`: readable numerical report. No manuscript or submission figure is replaced.

Review whether the objective continues improving beyond the original limits, whether it becomes insensitive to a target and whether another parameter reaches a bound. If this model is adopted, any published profile intervals must be computed for this new model and final parameter domain. This workflow does not reuse the current manuscript's intervals or claim to have recomputed them.

## Return completed results

```bash
cd ~/HPSCProliferation-R2
test -f results/aggregate_day1_boundary_checks/COMPLETE && \
tar -czf aggregate_day1_boundary_checks_for_review.tar.gz \
  results/aggregate_day1_boundary_checks \
  logs/aggregate_day1_boundary_checks.log
ls -lh aggregate_day1_boundary_checks_for_review.tar.gz
```

The archive is in the repository's top-level folder on the server. Download it to the Mac and share it for review before changing the documents.
