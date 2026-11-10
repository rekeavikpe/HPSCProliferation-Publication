# Targeted profile optimization verification

This workflow preserves the 60% training cohort and all previously selected parameters, bounds, likelihood settings and prediction results. It recomputes the same 25 profiles identified in `config/training60_profile_repairs.csv`, using `ProfileLikelihood.jl`.

Each fixed-parameter optimization starts from the preceding profile solution, the profile MLE and the three nearest saved nuisance solutions from the previous full-range repair. BOBYQA optimizes each start. The lowest objective is retained. Near the confidence region, a bounded BFGS refinement provides an additional check. A guard rejects solutions worse than their supplied seeds. The optimizer log records the number of starts, the best starting objective and the spread between optimized starts.

Confidence limits use the nearest adjacent threshold bracket on each side of the MLE and the package profile spline. Where the profile returns above the threshold, the complete curve must be reviewed for disconnected regions or numerical problems. Reaching a bound without crossing does not establish a finite confidence limit.

Run `bash scripts/run_training60_profile_multistart.sh` on the server after pulling. Inputs are `results/training60_final/separate_ablations/profile_inputs` and `results/training60_final/profile_repairs`. Outputs are separate in `results/training60_final/profile_multistart`.

The workflow is resumable and rejects changed source hashes. Do not edit its source while it is running. A completion marker means computation finished, not that all intervals are approved. Review `audit.csv`, profile curves, nuisance solutions and optimizer records. A new likelihood maximum or remaining near-threshold dips requires investigation before finalizing the intervals. Original calibration and prediction outputs are not modified.
