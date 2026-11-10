# David’s October 20 revision

Run `bash scripts/run_october20_revision.sh`. Default output is
`results/david_comments_20261020`. This is a separate, resumable analysis directory.
The original fits and manuscript are not overwritten.

## Analysis contract

- Original publication models: 0, 1, 3, 4, 5, 6, 7. Original conditions: 1–6,
  including condition 4 = 500 mL. Storage folder IDs are not publication numbers.
- The existing 19 training and 12 selection culture IDs are checked. No new split.
- P75: Model 0 is exposure independent and reused. The existing Model 6 P75 fit
  is reused after checking its 100 starts, both split memberships and recomputed
  training likelihood. The other five candidates receive 100 starts each. If the
  saved Model 6 P75 fit is absent it is also fitted. Selected parameters minimize training NLL. The 12 selection cultures
  are used to score those parameters. No eleven-summary refit.
- Day 1: all eight formulations, 100 starts each. Measured day 1 means for all four
  states initialise each partition separately. None of those initial observations
  contribute to the objective or selection score. Feeding retains its absolute
  calendar and the death term retains culture age. The 500 mL day 1 comparison is
  labelled as day-1-initialised, not the original day-zero external prediction.
  Using measured day 1 states is an explicit initialisation assumption. No day 2–5
  selection measurements enter training. Bounds and active parameter sets remain
  those of the respective formulations. The early death parameter may become weakly
  identifiable after excluding day zero to day one.
- Cells: initial free set is mu_max, K_d, K_oG, Y_GLac, alpha. The other six
  parameters are fixed at their ESI-017 training values. The three cultures per line
  remain pooled for recalibration. These are in-sample fits, not 60/40 validation.
  All requested parameters receive ProfileLikelihood.jl profiles using the existing
  cell-line weighted-residual objective. Its residual scales are assumed, not
  independently estimated cell-line measurement variances. Profile intervals are
  conditional on this error model and the fixed ESI-017 parameters.
- No cell parameter is automatically dropped. Both lines retain the requested five
  parameters. Unbounded intervals produce `parameter_reduction_review_required`.
  Numerical issues produce `profile_numerical_review_required`. The report lists
  affected parameters explicitly. Any reduced model requires a separately justified
  choice and refitting/reprofiling. An unfinished review is not a completed
  identifiability claim.
- Existing ablations supply selection BIC, 500 mL likelihood scores, parameter-bound
  checks and actual response terms. No refit to 500 mL. The external Gaussian score
  uses the existing external-data SD/floor convention, exported separately from
  selection scores. Held-out `2*NLL + k*log(n)` is a predictive BIC-form score.
- Overall MRPE here pools eligible state/day relative errors, preserving the
  established calculation. It is not silently changed to an equal-state mean.

## Comment-to-output map (Word comment IDs)

| IDs | Deliverable |
|---|---|
| 0 | p75/selection_scores.csv. R1 exposure tables/figure to be restored during document editing. No eleven-summary rerun. |
| 1,4,20 | statistics/culture_values.csv, tests.csv. Culture-level day-5 aggregate sizes and within-culture day5/day1 cell fold expansion. Three-RPM Kruskal–Wallis and separate two-vessel rank comparison. No unsupported causal claims or inherited p-values. |
| 2,3 | Original numbering in new reports. Global document/legacy-figure relabel is a writing-stage check. |
| 5 | Code audit: src/Proliferation.jl simulate evaluates exposure(u[4]) at ODE evaluations. linear_exposure interpolates CFD size bins and clamps below/above endpoints. No measured-size time series drives the primary ODE. |
| 6,16 | reuse/model6_errors.csv, training and selection per-state/per-condition. |
| 7–10 | oxygen/condition_summary.csv, aggregate_core_sensitivity.csv, scaling_assumptions.json and figures/oxygen_contours.pdf. |
| 11 | Methods: general size-dependent growth limitation. Do not claim strongest post-day1 inhibition at 60 RPM without fold results. |
| 12 | cells/status.json, each final round’s errors.csv, parameter_policy.csv, profiles/audit.csv and landscape curves. |
| 13–15 | Parameter source/units/growth-rate explanation checked during writing. mu_max has units h^-1, K_d h^-1, alpha um, K_Agg um^3, K_G and K_Lac mM, K_oG mL/cell, Y_GX mM*mL/cell, Y_GLac dimensionless, K_tau Pa, K_eps m^2/s^3 for Model 6. Routed channel units follow their input. |
| 17,18 | reuse/ablation_scores.csv, ablation_parameters.csv, ablation_response_terms.csv, figures/ablation_selection_bic.pdf. Equal rounded MRPEs alone do not establish identical mechanisms. |
| 19 | reuse/fold_predictions.csv, statistics/culture_values.csv and day1/scores.csv. |
| 21 | reuse/shear_inhibition.csv recomputed along current Model 6 trajectories. |

## Oxygen assumptions and sources

Lamont and Scott (1970), https://doi.org/10.1002/aic.690160403, motivates
kL proportional to sqrt(D) times (surface dissipation/nu)^(1/4).
Neto (2022), https://doi.org/10.3390/bioengineering9080386, supplies the distinct
Sherwood correlation used to anchor 100 mL at 37 C. They are not the same model.
At fixed RPM, ideal geometric similarity and constant power number give average
EDR proportional to V^(2/3) and specific interface area proportional to V^(-1/3).
Assuming surface EDR follows that average gives kLa ratio = 5^(-1/6).
The alternative in oxygen/scaling_scenarios.csv also uses the available volume-average CFD EDR ratio as an explicitly labelled proxy sensitivity. Neither approach substitutes a measured surface EDR. Actual 500 mL geometry and surface dissipation have not verified the scaling assumption.
The new calculation therefore reports a scaling scenario, not a measured value.
Do not call the earlier geometry discrepancy a unit-conversion error without evidence.

There is no single Da field for six different peak cell densities. The one-panel
figure overlays each condition’s Da=0.1, 0.5 and 1 contours and nominal point on common linear
kLa and uptake axes. Values at the points appear in its legend and companion table.
Core oxygen calculations state packing, diffusivity, diameter and uptake assumptions.
They do not establish that all aggregates are oxygen sufficient.

## Resume, progress and collection

Rerun the identical launch command to resume. Source hashes and budgets must match.
A changed source requires a fresh output directory. Never delete completed fits to
force a different scientific conclusion. Logs are in the output logs directory and
progress.log records each finished stage. `COMPUTATIONS_COMPLETE` means the scripts
finished, not that all profile/assumption reviews passed. See review_status.json.

Collect with `python3 scripts/collect_october20.py` when the suite finishes.

## Updating an already running pre-correction suite

After pulling this correction, run `python3 scripts/october20_resume.py stop`.
This stops only the October 20 suite process tree. Then run
`python3 scripts/october20_resume.py migrate` and relaunch the usual command.
Use `--out PATH` for a nondefault output directory on both helper calls.
Migration permits only this exact driver/report correction and the new helper.
All other source hashes, saved inputs and budgets must match. It preserves P75,
day 1, statistics, oxygen and first-round cell checkpoints. Superseded reduced
rounds and downstream reports are moved to a sibling archive. It invalidates
only cells, figures and report completion markers. Do not delete the manifest.
