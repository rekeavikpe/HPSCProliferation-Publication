# Training-only final workflow

The original culture split is retained: 19 cultures for training and 12 cultures for selection. All days from a culture remain together. There is no subsequent fit to all 31 cultures.

## Parameter and data sources

- `reference_fit` copies the saved selected parameters of all seven candidates from `split6040_refined_hybrid_lac100_kg100`. It verifies the split, observations, noise weights and training likelihoods first. It does not optimize parameters.
- Selection BIC is `2 × selection NLL + k × log(n_selection)`. Selection MRPE uses the same held-out 40%. The likelihood weights come from training. This is a held-out BIC-form score, not a new maximum-likelihood fit to the selection data.
- Each boxplot value evaluates one saved training optimization run on the selection data. The selected parameter vector for each model remains the best training fit. The code never chooses a run by its selection score. Model 6 must remain the selected formulation or the workflow stops for review.
- The five training-condition curves, 500 mL predictions and media-exchange predictions use the copied training parameters. The fit grid contains training-culture observations, not all-culture means. Its sixth row contains the independent 500 mL observations and predictions.
- Existing 60% exposure-summary fits are verified against their saved training and selection likelihoods then evaluated at 500 mL. No exposure optimization is repeated.
- Cell-line recalibrations inherit unestimated parameters from the training-only reference. Both the five-metabolic-parameter and seven-parameter comparisons are retained, alongside the other approved comparisons.
- The six separate-response ablations are fitted to the 19 training cultures with 100 hybrid starts each. Their held-out BIC distributions and MRPEs are evaluated on the same 12 selection cultures. Model 0 and Model 6 are reused. The 500 mL outcomes do not select fits.
- All-active and fixed-three profiles use `ProfileLikelihood.jl` and the training observations only. Each formulation's fixed K_G, K_Lac and Y_GX values come from its own training fit. The profile optimizer may refine the training optimum. This does not introduce selection or validation observations. Closed intervals are not guaranteed.

Original all-100-mL results are preserved as historical outputs. They are not inputs to this workflow. Measured-data statistics, structural identifiability, oxygen calculations and the separately designed pooled cross-validation do not require rerunning for this change.

## Run on the server

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
mkdir -p logs
nohup env JULIA_NUM_THREADS=42 \
  bash scripts/run_training60_revision.sh \
  > logs/training60_final.log 2>&1 < /dev/null &
echo $!
```

Defaults:

```text
Saved split: results/split6040_refined_hybrid_lac100_kg100
Saved exposure fits: results/david_revision_20261014/exposure_selection
New results: results/training60_final
```

The script accepts these three paths as positional arguments in that order. It runs prediction exports and figures first, followed by cell-line refits, ablations and both profile sets. The latter stages are computationally expensive. Repeating the same command resumes completed stages and fitting checkpoints. Do not launch a second copy. A workflow lock prevents concurrent runs on Linux. Changes to inputs, code or budgets require a new output folder so different analyses are not mixed.

Progress:

```bash
cat results/training60_final/progress.log
tail -n 40 logs/training60_final.log
pgrep -af '[j]ulia'
```

## Results

- `reference_fit/selection_scores.csv`: seven-model selection BIC
- `reference_fit/selection_state_errors.csv`: selection MRPEs
- `reference_fit/selection_optimizer_bic.csv`: held-out BIC for each saved training run
- `reference_fit/model*/selected_parameters.csv`: unchanged training parameter vectors
- `figures/model_selection_40.*`: selection BIC and MRPE panels
- `figures/model6_six_by_four.*`: five training conditions plus 500 mL validation
- `figures/exposure_selection_and_transfer.*`: exposure selection and external prediction comparison
- `prediction_metrics/prediction_metrics.csv`: training and external prediction errors and experimental noise
- `media/figures` and `cell_transfer/figures`: transfer figures
- `cell_lines`: updated line-specific refits and their error tables
- `figures/cell_lines_*.{pdf,png,...}`: cell-line comparison figures
- `separate_ablations/selection_ranking.csv`: held-out ablation BIC
- `separate_ablations/selection_state_errors.csv`: held-out ablation MRPEs
- `separate_ablations/figures`: actual held-out BIC boxplots with 500 mL MRPEs
- `separate_ablations/profiles/all-active` and `profiles/conditional`: profile results for all eight formulations

The legacy storage ID 5 corresponds to publication Model 6. Storage ID 6 corresponds to publication Model 5. Figure labels use publication numbering.

## Validation

The adapter and prediction paths were checked locally against the archived final split. All seven selected training likelihoods were reproduced. The 22-check cohort test verifies 19/12 separation, unchanged external observations and preservation of training observations after fixing parameters or changing hydrodynamic routing. Prediction checks passed for media exchange (170 checks) and cell-line transfer (102 checks). Low-budget smoke runs completed for all six ablations, both cell-line refit suites and one ProfileLikelihood.jl parameter. Full 100-start optimization and full profile grids run on the server. A local smoke test is not a substitute for those final runs.
