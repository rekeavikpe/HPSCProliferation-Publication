# Model implementation and numbering

The public scripts accept models 0–6. The internal equation kernel retains historical IDs to permit direct verification against archived results. Public IDs 0, 1, 2, 3, 4, 5, 6 map to historical IDs 0, 1, 3, 4, 5, 6, 7 respectively. **Model numbers in the equation notes below are historical IDs.**

The Julia module implements the manuscript's four-state equations with shared growth and metabolic calculations. Models 0, 1, 3, 4, 5, 6 and 7 are the distinct candidates. Model 2 is retained as an alias for reproducibility but excluded from `UNIQUE_MODELS`, because its formula is algebraically identical to Model 1.

## Equations and units

State order is viable cells/mL, glucose mM, lactate mM and aggregate diameter µm. Time is hours. The equations follow the manuscript baseline equations and model-variant table, checked against `hPSC_Proliferation_Paper_v18_clean.tex` lines 127–216 and the archived implementation.

| Parameter | Units required by the implemented equations |
|---|---|
| K_Agg | µm³ |
| K_d, mu_max, delta_max | h⁻¹ |
| K_G, K_Lac | mM |
| K_oG | mL/cell |
| Y_GLac | mol glucose/mol lactate |
| Y_GX | mM·mL/cell, equivalent to µmol/cell |
| alpha | µm |
| K_tau | Pa |
| K_eps | m²/s³ |
| gamma | dimensionless |

Model 4 retains the archived interpretation: shear changes cell growth, while metabolism and the aggregate growth contribution use the unsheared growth rate. This detail should be explicit in the final model table. Model 3 uses the stated rational approximation, not an exponential function.

The historical state-evaluation floors (1 cell/mL, 0.001 mM glucose, zero lactate, 1 µm diameter) and a 12-hour death adaptation time are explicit compatibility choices. They are not measurements or newly inferred biological parameters. The floors affect RHS evaluation, not saved observations, and do not guarantee positivity for arbitrary parameter sets. The reference trajectories were checked for nonnegative predictions. Final fitting must reject or handle physically invalid trajectories explicitly.

## Feeding and observations

The solver integrates each continuous interval up to the next exchange. It saves observations at that endpoint first, applies the exchange once, then starts the next interval at the same absolute time. There is no artificial 0.02-day delay, and the adaptation clock is not restarted.

For exchange fraction q, glucose becomes `(1-q)*G + q*G_fresh`, and lactate becomes `(1-q)*Lac` for lactate-free medium. Cells and aggregate diameter are retained. These are the existing model's idealized mixing/retention assumptions. The daily sample-volume removal and possible losses due to centrifugation are not newly modelled. Condition 6's handling is recorded in feeding metadata.

## Exposure input

The new CFD importer extracts 100 P50/P75 records directly from the verified rebuilt workbook, with cell references. The module uses linear interpolation in diameter and nearest-endpoint values outside the available size grid. This avoids silently extrapolating beyond the CFD grid but remains a modelling assumption to report, especially for early small aggregates. The publication exposure-summary comparison is a separate pending stage; this repository's main reproduction command uses P50.

## Likelihood

The module supplies one normalized fixed-sigma Gaussian NLL and BIC formula. It requires explicit positive sigma values and an explicit parameter count. It does not yet select a residual model, infer independence across days, choose SD versus SEM, apply uncertainty floors, or decide whether to fit replicate records versus condition means. Those calibration decisions must be recorded before fitting. Censored values require an explicit treatment; they cannot be supplied as measured zeros.

## Verification

- 35 packaged assertions cover zero-exposure limits, metabolic balance, Model 1/2 equivalence, Model 4 behavior, exact dilution, pre-exchange observation order, an analytic death trajectory across an exchange, exposure interpolation and normalized scoring.
- A one-time migration check compared 800 RHS evaluations against the archived models across all eight historical identifiers. All matched to relative tolerance 1e-12. This audit used the preserved archive; it is not a dependency of the clean module.
- Six Model-6 P50 reference simulations apply 19 exchanges. Each is rerun at tighter tolerances and must agree within 1e-6 after scaling by `max(abs(prediction),1)`.

The saved reference parameters are from the earlier conditional P50 diagnostic run, not a new fit. Historical starting values (14 mM glucose, except 15 mM for Condition 4; zero lactate; 15 µm diameter) are used only to check implementation continuity. Their policy must be reconciled with day-zero measurements before final calibration.

Fresh calibration, model selection and regenerated paper figures have not been performed in this stage.
