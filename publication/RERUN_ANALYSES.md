# Optional analysis reruns

The archived outputs already contain the final results. `reproduce.py` checks them without optimization. The commands here are for new analyses and can take hours or days. Run from the release root with Julia 1.12.4 and instantiated environments. Use new output folders, never `publication/cache/`.

## Primary calibration and selection

The final workflow uses the same seeded 19/12 whole-culture split. Population searches use 100 starts followed by the archived refinement routine. Model 6 remains the fit from the 19 training cultures for validation.

```bash
export JULIA_NUM_THREADS=42 OPENBLAS_NUM_THREADS=1
julia --project=population population/run_split_selection.jl config/calibration_lac100_kg100.toml results/new_split
julia --project=population population/refine_split_selection.jl results/new_split results/new_split_refined
julia --project=population revision/prepare_training60.jl results/new_split_refined results/new_training/reference_fit
```

`prepare_training60.jl` verifies selected values and the 19/12 cohort, and copies them into an evaluation adapter; it does not refit. Some legacy reports use intermediate numbering. Use the final mapping in `publication/model_numbers.csv`. Do not follow the historical suggestion to refit on all small-vessel data before validation.

## Frozen validation

```bash
julia --project=. scripts/predict_media_exchange.jl publication/cache/training60/reference_fit results/new_media
julia --project=. scripts/predict_cell_lines.jl publication/cache/training60/reference_fit results/new_cell_transfer
```

For 500 mL predictions, `scripts/predict_scale.jl FIT` writes a `scale_transfer` directory inside FIT. Use a copy of the reference fit as FIT to preserve the archive. `publication/verify_primary.jl` already re-evaluates both reported 500 mL models without changing the fit directory.

## P75, statistical and oxygen analyses

```bash
julia --project=population revision/run_october20.jl p75 publication/cache/training60/reference_fit publication/cache/training60/reference_fit results/new_p75
julia --project=revision/statistics revision/october_statistics.jl results/new_statistics
python publication/fetch_reference.py
julia --project=population scripts/oxygen_sensitivity.jl results/new_oxygen_reference
python3 scripts/october_oxygen.py results/new_oxygen_reference results/new_oxygen
```

The oxygen scaling pass is required: the intermediate correlation reference contains the earlier assumed 500 mL coefficient. The final scaling pass gives 1.18 h⁻¹. These remain analytical estimates with the documented assumptions. The P75 default is 100 starts with 20,000 population evaluations and 3,000 local iterations; see saved `manifest.toml` and the `HPSC_REVIEW_*` environment settings in `revision/October20.jl`.

## Nine Model 6 profiles

```bash
HPSC_PROFILE_FULL_RANGE=0 julia --threads=9 --project=profiles profiles/run_fixed_two.jl publication/cache/training60/reference_fit results/new_nine_profiles publication/cache/profile_seeds_historical 200 all
python3 scripts/audit_fixed_two_profiles.py results/new_nine_profiles
```

This fixes K_Lac and Y_GX and profiles the other nine parameters. It uses isolated per-parameter seed banks. The historical seed-directory intervals are not publication results. The saved final profiles include the successful recovery/reoptimization diagnostics.

## Two fixed-parameter diagnostics

```bash
julia --threads=42 --project=profiles profiles/run_fixed_parameter_diagnostics.jl publication/cache/training60/reference_fit publication/cache/profile_seeds_historical results/new_two_diagnostics
```

This profiles K_Lac and Y_GX separately, with the other ten parameters free, and writes point-level checkpoints. Inspect the numerical audit before interpreting limits.

## iPSC recalibrations

For each line (`585A1` and `SCTi003-A`), use a new output folder:

```bash
julia --threads=42 --project=population revision/fit_october_cells.jl publication/cache/training60/reference_fit 585A1 results/new_585A1_three mu_max,K_oG,alpha
julia --project=profiles profiles/run_october_cells.jl results/new_585A1_three results/new_585A1_three/profiles
julia --project=population revision/export_october_cells.jl publication/cache/training60/reference_fit results/new_585A1_three
```

For Recalibration 1, the parameter list is `mu_max,K_d,K_oG,Y_GLac,alpha`. The final three-parameter intervals are in `cache/ipsc/`; the original five-parameter profiles are in `cache/october20/cells/`. All line-specific errors are in-sample.

## Leave-one-condition-out cross-validation

`population/run_cross_validation.jl config/calibration_lac100_kg100.toml OUTPUT` contains the original prespecified LOCO and additional historical holdout designs. See `population/CrossValidation.jl::plans` for definitions. The accepted-paper subset consists of the five `loco_condition*` folds and manuscript Model 6 (storage 5); `publication/verify_cross_validation.jl` replays exactly those five from saved parameters. The whole historical driver runs additional analyses too; it is not necessary to reproduce Table S9.

All reruns should retain the supplied data, bounds, likelihood and weights unless a new analysis is intended. Examine convergence and scoring before replacing a saved result.
