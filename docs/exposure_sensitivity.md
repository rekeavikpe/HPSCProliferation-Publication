# Hydrodynamic exposure-summary sensitivity

Run from the repository root after the final Model 5 calibration:

```bash
JULIA_NUM_THREADS=42 bash scripts/run_exposure_sensitivity.sh
```

Defaults: source `results/calibration100_hybrid_lac100_kg100`, destination `results/exposure_sensitivity_hybrid_lac100_kg100`. Two optional arguments override these paths. Existing destinations are rejected. This is a new analysis, not a rerun or replacement of the primary calibration.

The comparison holds public Model 5 (legacy Model 6) fixed and estimates all 11 parameters separately for each of 11 supplied input summaries:

- Lagrangian mean, P25, P50, P75, P90 and maximum.
- Eulerian volume average, P50, P75, P90 and maximum.

All scenarios use the same five 100-mL conditions, 116 condition–day–state means, measured initial glucose, feed schedules, likelihood and sigma floors, K_G/K_Lac upper bounds of 100 mM, and normalized logarithmic parameter coordinates. Each receives 100 independent population searches (50 candidates, 20,000 evaluation budget) followed by up to 3,000 iterations of the same direct Nelder–Mead refinement used in the primary hybrid calibration. The seed schedule matches the primary Model 5 calibration and is identical across exposure summaries. Every scenario is independently searched; no prior best fit is used to initialize or stop searches. Parameter counts remain 11 for each scenario. The reference P50 BIC is recorded only for numerical comparison.

Lagrangian inputs depend on predicted aggregate size through the existing linear interpolation and endpoint clamping. Eulerian inputs are constant at each vessel/RPM. The tracked source table contains workbook cells. Its importer verifies P50/P75 against the primary inputs and Eulerian volume averages against the separate volume-average workbook. Blank, descriptive and unidentified workbook rows are excluded. Shear and EDR summaries are marginal statistics; paired maxima need not occur simultaneously. No raw trajectory distribution or CFD weighting definition is reconstructed.

`tools/import_exposure_summaries.py` regenerates the 350 input records from the three original CFD workbooks under `data/source` (requires openpyxl). Imported values and workbook hashes are tracked, so server fitting does not require the workbooks or Python Excel dependencies.

The pipeline verifies identical training observations and exact P50 lookup equivalence before fitting. It retains all optimization runs, scores, parameters, status and boundary flags. Selection uses minimum achieved training NLL, preferring a converged result only within 1e-6 NLL of that minimum. Failed starts are retained and prevent finalization. Selected likelihoods are checked, and tighter solver predictions must agree within 1e-4 of the fitting sigma. The script does not consume 500-mL, media-exchange or cell-line outcomes.

Outputs include `ranking.csv`, `runs.csv`, `predictions.csv`, `state_errors.csv`, source hashes, selected parameter vectors, an HTML review page and six-format figures for BIC and statewise MRPE. Experimental CV uses measured sample SD and is distinct from the fitting uncertainty floor. BIC boxplots represent optimizer variability. Some starts or selected fits may not meet stopping criteria; these are retained and explicitly marked.

Interpretation is conditional on Model 5. This does not establish that its model-structure ranking is unchanged for every exposure summary. Any switch of primary exposure/model needs an explicitly defined further analysis. Keep P50 as the declared primary analysis while evaluating this sensitivity; report alternatives that perform better or worse. The Gaussian score remains a working likelihood without repeated-day covariance. No global optimum or practical identifiability is guaranteed.

The local `results/exposure_sensitivity_workflow_check` uses two starts per scenario and the local workflow-check calibration as its configuration/reference. It tests the full implementation; it is not the final 100-start server result.
