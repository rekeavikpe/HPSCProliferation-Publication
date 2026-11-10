# Identifiability workflow — original package approach restored

**Numbering update (5 October 2026):** Historical examples and CLI arguments below retain storage IDs. Storage Model 5 is now publication **Model 6**; storage Model 6 is publication **Model 5**. See [the crosswalk](model_numbering.md). Commands remain valid.

The production practical-identifiability workflow uses **ProfileLikelihood.jl 0.4.1**, following the original `scripts/22_profile_likelihood_all11.jl` and `04_profile_likelihood.jl`:

1. Construct `LikelihoodProblem` with the saved 100-mL likelihood and bounds.
2. Refine the saved calibration estimate using `mle` with Fminbox/BFGS and BackTracking.
3. Use `construct_profile_ranges` and `ProfileLikelihood.profile` with NLopt `LN_BOBYQA`, previous-solution starts (`:prev`), and nominal 95% package confidence intervals.
4. Export the package's profile values, nuisance parameters and confidence intervals; restore the package plotting object from these CSVs and call `ProfileLikelihood.plot_profiles` with the original manuscript settings, using CairoMakie as its backend.

The original bound-safe restart method is retained as an explicit, version-pinned compatibility override in `profiles/PackageWorkflow.jl`. Package installations are not edited. The small BOBYQA wrapper only records solver return codes; it delegates optimization directly to NLopt through Optimization.jl. There is no custom profile search engine in this workflow.

## Documented improvements and differences

Parameters use normalized logarithmic coordinates to accommodate their different orders of magnitude. Values and package interval endpoints are converted back to physical units for reporting. The mathematical model, saved bounds and likelihood remain those of the 100-mL calibration. The MLE is checked against the saved objective before profiling. Profile objectives are independently recomputed from the saved nuisance estimates. Source hashes and an MLE checkpoint protect resume consistency; completed parameters are skipped. An interrupted parameter is restarted from the package MLE.

The original resolution setting of 100 is retained: in ProfileLikelihood.jl this means 100 candidate positions **per direction**, not exactly 100 final points. The package stops at the likelihood threshold and may refine short branches. Parameter profiles and their two directions can run concurrently, but 42 requested threads do not guarantee 42 simultaneous useful tasks.

The separate `profiles` environment pins DiffEqBase 6.155.0 because the package's SimpleNonlinearSolve 1 dependency fails with newer DiffEqBase releases. LineSearches 7.3.0 is also pinned for compatibility with the resolved Optim version. Its manifest pins all dependencies. The saved calibration likelihood must match within 1e-6 before a run proceeds. The main calibration and figures environments are unchanged.

## Run on the server

From the repository root, after pulling the updated files:

```bash
julia +1.12.4 --project=profiles -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --threads=4 --project=profiles profiles/test.jl
julia +1.12.4 --threads=42 --project=profiles profiles/run.jl results/calibration100 results/profilelikelihood100 5 100
julia +1.12.4 --project=profile_figures -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=profile_figures profile_figures/render.jl results/profilelikelihood100 figures/output/profilelikelihood100
```

Repeat the production command to resume. Use a new output directory for changed settings or source code. The `--smoke` option is for software testing with reduced iteration limits only and must not be used for final inference.

Read `mle_summary.csv`, `profile_confidence_intervals.csv` and each parameter's `profile.csv`. A package interval endpoint at a search bound is not a finite confidence limit. `unconfirmed_points` and `better_minimum` require review. A returned interval alone does not establish convergence or calibrated coverage under parameter bounds and correlated repeated observations. Fixing parameters later produces conditional inference; fixed values must come solely from the 100-mL analysis.

## Earlier custom profiles

The previous custom implementation did not follow the requested package workflow. Its output is retained only as an exploratory cross-check. `scripts/profile_all.jl` now stops with a migration message, preventing accidental continuation as the production analysis. The custom profiles are not relabeled as package results.

## Structural analysis and subsequent evaluations

All-seven-model structural analysis remains in `structural/`; see `docs/structural_identifiability.md` for completed results and qualifications. No original structural script was found in the supplied source. The rational-segment analysis and nominal-initial-state rank checks are distinct from practical identifiability.

After reviewing package profiles and fixing the parameter policy, evaluate the media-exchange protocols and additional cell lines. Resolve cohort/protocol mapping first. Any parameters fitted to an evaluation cohort must be reported as adaptation, rather than unchanged-parameter validation.

## Verification on 4 October 2026

Seven analytic package-workflow checks passed, including known Gaussian profile intervals, objective values, solver-code recording and bound-safe starts. A reduced-iteration Model 5 K_Lac integration run completed: the package MLE returned Success and matched the saved NLL within 6e-10. All nine exported profile objectives were recomputed successfully, and the test correctly retained MaxIters flags. Resuming this completed run succeeded without refitting. This smoke run is a software check, not an inferential result. The full 11-parameter, resolution-100 package run was started separately and was not complete when this note was written.

## Native manuscript plotting

The figure now comes directly from `ProfileLikelihood.plot_profiles`, using the original all-11 script: 3 rows, 4 columns, mathematical parameter labels in original order, `show_mles=true`, `shade_ci=true`, package-default splines and linear physical axes, original dimensions and typography, and four x-axis ticks. The separate `profile_figures` environment provides the Makie version supported by the package without modifying the fitting environment. The plot-only adapter restores exported package curves, MLEs and interval endpoints; it does not run optimization or recompute intervals. Native CI titles need the qualifications in the companion CSV when bounds are reached or nuisance fits are unconfirmed. PNG/JPEG/TIFF export at 600 dpi; PDF/EPS/SVG also exported. The prior custom renderer now stops with migration instructions.

Julia 1.12 compatibility: Makie 0.21 uses the removed internal `TypeName.mt` when naming plot functions. The plotting entry point substitutes the public `nameof` equivalent for that naming helper only. `profile_figures/LocalPreferences.toml` disables CairoMakie’s precompile drawing workload so this runtime compatibility fix can load. No installed package source is edited; ProfileLikelihood plotting routines and numerical fitting remain unchanged.

## Conditional analysis fixing K_Lac and Y_GX

Use the saved 100-mL calibration estimates for both parameters, leaving nine free parameters. Values are read from the chosen source fit, not copied from the original combined-data manuscript. Re-estimation and profiling still use ProfileLikelihood.jl, Fminbox/BFGS and BOBYQA. The all-free results remain separate.

```bash
julia +1.12.4 --threads=42 --project=profiles profiles/run.jl results/calibration100 results/profilelikelihood_fixed2 5 100 all --fix-lactate-yield
julia +1.12.4 --project=. profiles/compare_fixed.jl results/profilelikelihood_fixed2
julia +1.12.4 --project=profile_figures profile_figures/render.jl results/profilelikelihood_fixed2 figures/output/profilelikelihood_fixed2
```

The output includes `fixed_parameters.csv`, the conditional `calibration.toml`, the re-estimated `mle_parameters.csv`, and native package profiles in a 3-by-3 layout. These are conditional intervals, not evidence that the fixed parameters became identifiable. Because the two values were estimated from the same training data, this exercise does not justify reducing the original BIC parameter count from 11 to 9. The comparison script evaluates the already fitted parameters at 500 mL only after calibration; that cohort never enters re-estimation.

Local conditional run (4 October 2026): K_Lac = 10000.00000000001 mM and Y_GX = 3.2126404034728867e-6 mmol/cell, read from the local saved 100-mL calibration. The nine-parameter MLE returned Success (NLL 613.0855788295045); all nine profile exports completed and their retained points reported successful stopping. MRPE changes were below 0.000016 percentage points. Y_GLac's package interval narrowed to approximately 0.420–0.496. K_G remained above the likelihood threshold at its upper search bound of 400 (normalized log likelihood −0.655); mu_max likewise remained above threshold at 2 (−1.613). The package's upper endpoints at the MLEs for these two parameters are not valid finite upper confidence limits. This conditional analysis has not resolved those upper uncertainties.

## K_G upper bound reduced to 40 mM

The requested comparison keeps K_Lac and Y_GX fixed at the source 100-mL MLEs and changes only K_G's allowed range from [0.1,400] to [0.1,40] mM. The saved source calibration is preserved. The override is recorded in the new run's settings and calibration files; an out-of-range source MLE is rejected rather than silently clipped.

```bash
julia +1.12.4 --threads=42 --project=profiles profiles/run.jl results/calibration100 results/profilelikelihood_fixed2_KG40 5 100 all --fix-lactate-yield --kg-upper=40
julia +1.12.4 --project=. profiles/compare_fixed.jl results/profilelikelihood_fixed2_KG40
julia +1.12.4 --project=profile_figures profile_figures/render.jl results/profilelikelihood_fixed2_KG40 figures/output/profilelikelihood_fixed2_KG40
```

Interpret any tightened intervals conditional on this upper bound. A profile ending above the likelihood threshold at 40 mM is bound-limited, not an identified upper confidence limit. This run does not change the optimizer or fix K_Agg.

Completed local KG40 comparison: the refit returned Success with NLL 613.0855788344838. All nine profile exports completed without unconfirmed retained points. The package interval for mu_max is approximately 0.0940–0.2684 h^-1, conditional on K_G <= 40 mM. K_G still has no upper likelihood-threshold crossing within that range; its package upper endpoint at the MLE is a fallback and must not be reported as a finite confidence limit.

## Inactive lactate term with conditional yield

Run `bash scripts/run_inactive_profiles.sh results/calibration100` from the repository root. It uses Julia 1.12.4 and 42 threads by default (`JULIA_NUM_THREADS` overrides this). The saved 100-mL-only calibration provides the initial parameters and fixed Y_GX; no 500-mL data enter fitting or profiling. The source likelihood is verified before modifying the formulation.

The active comparison reoptimizes all parameters except Y_GX. The inactive comparison sets the lactate inhibition factor exactly to one, removes K_Lac from the free vector, and profiles the remaining nine parameters using ProfileLikelihood.jl. K_Lac=1 in its parameter file is an unused placeholder, not an estimated biological constant. No new K_G bound is imposed. Original-source bounds remain in force.

Outputs are `results/inactive_lactate_yield_active_fit` and `results/inactive_lactate_yield`. Native package profiles are in `figures`; prediction sensitivity plots are in `prediction_figures`. Medium-exchange jumps are retained. CSVs include fit likelihoods, per-state MRPE and the accepted profile vectors. The envelope is the pointwise range over converged sampled profile-path vectors within the nominal one-parameter likelihood reference. It is conditional on Y_GX, not a simultaneous 95% prediction interval, and excludes measurement noise and model discrepancy. A substantially better likelihood found during profiling stops envelope construction and requires refitting. Ordinary BIC is not reported for this plug-in fixed-yield comparison, because that yield was estimated from the same data.

To resume, repeat the command with unchanged sources and inputs. Source changes require a fresh output prefix, supplied as the second argument. Do not run two processes against the same prefix. A smoke run is not a publication analysis.

## Additional fixed-parameter diagnostic

To test one additional parameter at the current inactive-model MLE, append `--conditional-from=results/inactive_lactate_yield --fix-at-mle=K_Agg` to the inactive profile command and use a new output directory. The parent configuration, source calibration, convergence and likelihood are verified before its parameter value is fixed. `K_G` is an alternative diagnostic directly targeting the K_G–mu_max trade-off. Neither is fixed at a chosen boundary or an old combined-data estimate.

For example:

```sh
julia +1.12.4 --threads=42 --project=profiles profiles/run.jl results/calibration100 results/inactive_fixed_KAgg 5 100 all --inactive-lactate --fix-yield --conditional-from=results/inactive_lactate_yield --fix-at-mle=K_Agg
julia +1.12.4 --project=profile_figures profile_figures/render.jl results/inactive_fixed_KAgg results/inactive_fixed_KAgg/figures
```

These eight-parameter profiles condition on an additional estimate obtained from the same data. Closed conditional intervals do not establish identifiability of the original nine-parameter fit. Retain the parent profiles and every tested fixation as the diagnostic record. Use `profiles/compare_additional_fixed.jl` to check likelihoods, both threshold crossings, optimizer status and nuisance-parameter boundaries for all compared suites. The native eight-panel layout follows the original fixed-three script (two rows, four columns).

## Corrected fixed-three analysis: original lactate inhibition retained

`bash scripts/run_fixed_three_profiles.sh results/calibration100 results/active_fixed_three` fits and profiles eight free parameters while holding K_G, K_Lac and Y_GX at their exact saved values in the source 100-mL calibration. It uses `profiles/run.jl ... --fix-three`. This preserves the original finite lactate-inhibition factor K_Lac/(L + K_Lac); it does not use the inactive model or its MLEs. The default is 42 threads, controlled by JULIA_NUM_THREADS. No bounds are changed. Outputs include the native two-by-four package figure and verification tables under `review/`.

The finite K_Lac estimate in the source fit is at its upper bound. Fixing that constrained MLE must be disclosed; it is not an independently known or practically identified constant. Remaining intervals condition on all three fixed estimates. Retain the active nine-parameter profiles (K_Lac and Y_GX fixed, K_G free) to document the K_G–mu_max trade-off. Keep inactive-inhibition results only as a separately labelled sensitivity analysis. Do not describe this corrected fit as inactive inhibition.

## Historical lactate-value sensitivity

```sh
julia +1.12.4 --threads=42 --project=profiles profiles/historical_lactate.jl results/calibration100 results/historical_lactate100 8
julia +1.12.4 --project=profile_figures profile_figures/historical_lactate.jl results/historical_lactate100
```

This holds only K_Lac=30.40495 and reoptimizes the other ten parameters with eight deterministic starts using ProfileLikelihood.jl's MLE interface and NLopt BOBYQA. All bounds and 100-mL observations remain unchanged. The source all-free likelihood is verified. Starts, objective values, convergence, parameter boundaries, statewise MRPE and standardized residual errors are saved. The 500-mL prediction is retrospective and does not influence optimization. Because the historical value was obtained from combined data, this is a sensitivity analysis, not an independently specified constant or leakage-free validation. The likelihood penalty and statewise error changes—not visual smoothness—determine the empirical cost of that assumption. No BIC advantage is claimed by treating a previously estimated constant as independently known.

## Candidate 100-mM lactate-bound rerun (4 October)

Run the new all-free, seven-model comparison into a separate directory:

```sh
julia +1.12.4 --threads=42 --project=. scripts/calibrate.jl results/calibration100_lac100 100 config/calibration_lac100.toml
```

This sets 100 initial starts per model and retains the original finite lactate-inhibition equation. The 100-mM upper bound is an operational constraint, not a demonstrated biological maximum. A seven-model workflow test completed with one initial start and two refinement starts per model; its output is `results/lac100_workflow_smoke` and is not the final comparison.

Inspect the final ranking and optimizer diagnostics before choosing the model and fixed values. If Model 5 remains selected, derive K_G, K_Lac and Y_GX from this new calibration and run:

```sh
bash scripts/run_fixed_three_profiles.sh results/calibration100_lac100 results/active_fixed_three_lac100
```

Do not run that Model-5-specific command automatically if another model is selected. Conditional profile intervals do not include uncertainty in the three fixed estimates. Keep the unrestricted profiles and broader-bound sensitivity as supporting results. Closure must be verified from actual likelihood crossings, not from plotting limits.

## Updated candidate: K_G and K_Lac capped at 100 mM

The separate configuration below supersedes the lactate-only candidate for the next production run. The earlier configuration and saved results remain unchanged. K_G is the glucose half-saturation parameter, not the glucose concentration state. The choice of 100 rather than 50 mM is the less restrictive of the requested options; it does not establish practical identifiability.

```sh
julia +1.12.4 --threads=42 --project=. scripts/calibrate.jl results/calibration100_lac100_kg100 100 config/calibration_lac100_kg100.toml
```

Use this new output directory even if the lactate-only run has already started. Changing a configuration file does not change an already running process. The earlier seven-model smoke test used K_G <= 400 and does not validate numerical convergence under this revised bound. The new config has been checked for correct bounds and no fixed parameters; final convergence and ranking require the new run.

## Population-search comparison

`population/run.jl` is a separate diagnostic using BlackBoxOptim 0.6.12's original `adaptive_de_rand_1_bin_radiuslimited` method. It keeps the current 100-mL cohort, fixed-SD likelihood, bounds, normalized-log coordinates and local refinement (300 least-squares iterations, 3000 Nelder–Mead iterations). No selected parameters seed the population. Each independent search uses 50 population members and a declared evaluation budget, rather than a machine-dependent wall-time budget. The primary calibration, selection lock and profiles are not overwritten.

The first comparison is five searches each for Models 0 and 5, with 20,000 evaluations per search. Seeds are 20261004 + 1000 × model ID + run ID. Each task seeds its own random-number generator; searches run in parallel across Julia tasks. The pinned environment is separate from the existing calibration environment.

```sh
julia +1.12.4 --project=population -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --threads=42 --project=population population/run.jl config/calibration_lac100_kg100.toml results/population_pilot100 5 20000 0,5 results/calibration100_lac100_kg100/figure_data/optimizer_bic.csv
```

The optional reference CSV supplies comparison scores only. It never initializes or stops optimization; it must come from the matching calibration cohort and likelihood. Every population and local stage is retained. `runs.csv` reports individual outcomes, `summary.csv` gives model-level spread, and any failed runs appear in `failures.csv`. Saved objectives are reevaluated before completion. Reaching the population evaluation budget is not claimed to demonstrate convergence. Five repeats are a diagnostic, not a replacement for the final 100-run BIC distributions.

Completed local population diagnostic: five repeats for all seven models (35 fits), in `results/population_pilot100` and `results/population_other_models100`. All 35 finished within 2 BIC units of the matching best saved reference. Models 0/1/2/4/5 had BIC ranges below 0.00001, Model 3 about 0.740 and Model 6 about 0.451. No failed runs were omitted. The assembled standard-format output is `results/population_pilot_all7`; all 133 likelihood, stage and ODE-precision checks passed.

Stopping flags must be read at stage level: the strict-minimum stage succeeded in 4/35 runs; 9/35 had a successful stage within 1e-6 NLL of their minimum. The existing assembly policy therefore selects a successfully stopped Model 5 fit at BIC 1280.666532778122. Models 1, 2, 3 and 4 retain unresolved stopping flags. These facts support objective repeatability, not a global-optimum proof or resolved parameter identifiability. `results/population_comparison/index.html` reports all outcomes.

The server helper `bash scripts/run_population100.sh` prepares 100 population searches per model, preserves original results, converts output to the standard calibration layout and verifies all saved likelihoods. It has not been run at the 100-repeat scale. Resolve the remaining local stopping diagnostics before presenting that production analysis as finalized.

## Production hybrid requested by the author

The original `src/calibration.jl::calibrate_model` uses differential evolution followed directly by Nelder–Mead, with clipping in `scaled_to_params`. Its returned success flag is unconditional when no exception occurs; the new implementation records `Optim.converged` instead. `population/HybridRefinement.jl` restores the direct local stage and clipped-coordinate objective. The current 100-mL likelihood, bounds and log coordinates remain unchanged. Population size 50, a 20,000-evaluation search limit, and a 3,000-iteration local cap are explicit revisions to the historical settings. Local optimization can stop earlier.

The 35 saved population candidates were first checked with direct original-style Nelder–Mead at 3,000 iterations. All five baseline and Model 5 runs met the stopping criterion; selected runs for Models 1, 2, 4 and 6 did too. Model 3 required further diagnostics. This test is separate from the earlier population → least-squares → Nelder–Mead diagnostic.

After synchronizing the repository, run on the server:

```sh
JULIA_NUM_THREADS=42 bash scripts/run_population100.sh
```

The helper runs 100 independent searches for each model, retains every stage, assembles the selected fits, verifies likelihoods and tighter ODE predictions, generates 500-mL predictions, and exports the calibration and BIC figures. Output: `results/calibration100_hybrid_lac100_kg100`, with figures in its `figures` subfolder. Individual search results: `results/hybrid100_lac100_kg100_search`. Existing output folders cause a safe stop; they are not overwritten. The optional reference scores never seed or terminate fitting.

The longer Model 3 diagnostic (12,000 iterations) did not resolve its stopping flag and changed BIC negligibly. Production retains the verified 3,000-iteration cap rather than silently loosening the stopping criterion or extending all fits. The test output preserves the longer-budget metadata and source snapshot.

Verification of the restored original hybrid: 35 direct-local repeats at the production 3,000-iteration cap; all five Model 5 and baseline runs converged. Best runs for Models 1, 2, 4 and 6 also converged. The five longer Model 3 checks (12,000 iterations) still did not meet their stopping rule. The independent seven-model end-to-end run passed all 42 likelihood/stage/ODE-precision checks and generated predictions and plot inputs. Its longer local cap is recorded in its own metadata; the production helper retains 3,000. Report: `results/original_hybrid_refinement/index.html`.

### Frozen media-exchange evaluation

After confirming the hybrid calibration ranking, run:

```bash
bash scripts/run_media_exchange.sh \
  results/calibration100_hybrid_lac100_kg100 \
  results/media_exchange_hybrid_lac100_kg100
```

This reads the locked Model 5 estimates without re-estimation and also evaluates the saved Model 0 baseline. It uses the collaborator workbook's six labeled vessels: L85/L87/L89 (20%) and L84/L86/L88 (50%). The tracked measurement table includes workbook cell references and source hashes. Regenerate it with `python3 tools/import_media_exchange.py` only when the original workbook is available under `data/source` (requires openpyxl).

Both protocols exchange 50% on day 2, followed by their assigned fraction on days 3–5. Predictions at observation times are pre-exchange; plots include both sides of each exchange. Initial glucose is the measured protocol mean at day zero (20%: 14.4333333333 mM; 50%: 14.3666666667 mM). Fresh glucose is assumed equal to this value, consistent with the calibration convention; it is not a separately measured feed concentration. Initial aggregate size retains the calibration assumption, excluding workbook day-zero zero placeholders. This replaces rounded cell counts and the common 14 mM initial glucose in the archived validation CSV with source workbook values.

The runner verifies the selection hash, nonnegative finite trajectories, medium-exchange balance, pre-exchange sampling, and stability under tighter integration tolerances. It saves predictions, trajectories, experimental mean/sample SD, MRPE/RMSE/bias by state and protocol, exchange events, assumptions, and hashes. The figure uses the original blue/red 2×2 arrangement and exports PDF/SVG/PNG/JPEG/TIFF/EPS (raster formats at 600 dpi).

These experiments have previously been inspected; describe the evaluation as frozen-parameter protocol transfer rather than a newly blinded test. The plotted SD is experimental dispersion, not a predictive interval. The local `media_exchange_workflow_check` uses the local workflow-check calibration; final server results must be generated from `calibration100_hybrid_lac100_kg100`.

### MRPE and experimental variability across evaluations

After syncing the repository, run `bash scripts/run_evaluation_metrics.sh` on the server. This reuses existing complete media-exchange predictions, generates frozen cell-line predictions if absent, updates the media legend, and writes `results/prediction_metrics_hybrid_lac100_kg100/index.html` and `prediction_metrics.csv`. No kinetic optimization is launched. Existing evaluation parameter hashes must match the supplied calibration; mismatches stop reporting rather than combine different fits.

The table includes all seven calibration models, each original calibration condition separately, the 500-mL transfer, each media-exchange protocol, and both frozen cell-line transfers. Each has per-state rows and an all-state pooled-point row. MRPE is the mean absolute relative error against the observed means. Experimental noise is the mean replicate coefficient of variation (sample SD / absolute mean), following the original definition; it is not variance or an acceptance threshold. Days greater than zero are used. Zero means are excluded from relative measures and singleton measurements do not receive invented SDs. Counts and a matched-subset MRPE are reported. Actual replicate SDs are used, never the fitting sigma floors. Statewise RMSE and bias retain physical units; these are not pooled across states.

The cell-line importer selects workbook groups explicitly labeled no dextran sulfate: 585A1 C9/C13/C21 and SCTi003-A SCT3/SCT5/SCT7. All 108 post-day-zero measurements match the archived cell-line analysis CSV (cell counts and metabolites exactly; aggregate sizes within 0.000061% rounding). This preserves the original analyzed cohorts; it does not resolve the conflicting 585A1 manuscript Methods description. Initial glucose uses the measured day-zero line mean, with the inherited assumption that fresh medium has the same concentration. The original 585A1 50-rpm exposure interpolation between 40 and 60 rpm is retained and exported for inspection. Cell-line-specific recalibration remains separate and must be labeled as fitting rather than validation.

The evaluation helper also renders `results/cell_lines_frozen_hybrid_lac100_kg100/figures/cell_lines_frozen.png` (plus PDF/SVG/JPEG/TIFF/EPS). Its four-state, two-column layout retains the original cell-line colors and distinguishes frozen Model 5 predictions (gray dashed lines) from experimental means and SD. It contains no recalibrated curves. Existing predictions are reused when rerunning the helper; figures and metric summaries are regenerated.

### Separate iPSC recalibration

Run `JULIA_NUM_THREADS=42 bash scripts/run_cell_line_recalibration.sh` after the frozen evaluation. This launches 20 independent hybrid searches for each line (40 jobs), each with a population of 30, 20,000 global objective evaluations, and up to 3,000 Nelder–Mead iterations. The seeded evaluation budget replaces the archived 45-second global search for reproducibility. Results go to `results/cell_lines_recalibrated_hybrid_lac100_kg100`; existing results are never overwritten.

This preserves the archived cell-line method's seven free parameters and bounds: mu_max [0.05,0.5], alpha [10,200], K_G [1,50], K_oG [1e-7,1e-4], Y_GLac [0.1,1.5], K_d [0.01,0.5], K_Agg [1e6,5e8]. Only K_oG and K_Agg use logarithmic search coordinates, as in the original. K_Lac, Y_GX, K_tau and K_eps retain the exact saved ESI-017 values. Each line has its own biological parameter estimates. The objective retains the original state-scaled weighted SSR, with scales [500000 cells/mL, 5 mM, 5 mM, 100 µm] and weights [1,1,1,1.5]; it is not the main calibration Gaussian likelihood and must not be used for an unqualified BIC comparison.

Initial conditions, feeding, exposure interpolation, and observations are taken directly from the verified frozen evaluation. Statewise predictions are checked against a tighter integration, and the four fixed parameters are checked for exact equality. Every run's scores, actual convergence flag, boundary flags, seed, and parameter vectors are retained. The selected fit is the minimum SSR, preferring a converged solution only within 1e-8 SSR of that minimum. No global-optimum or parameter-identifiability guarantee is implied.

The comparison figure is `figures/cell_lines_comparison.png` (and five companion formats): frozen ESI-017 predictions are gray dashed, recalibrated curves are solid in the cell-line colors, and experimental means/SD are separate markers. `metrics/prediction_metrics.csv` labels frozen prediction and in-sample recalibration separately. Inspect `summary.csv` for convergence and boundaries before manuscript use.
