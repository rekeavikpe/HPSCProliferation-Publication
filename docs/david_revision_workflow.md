# Complete analysis run for David's review

Run `bash scripts/run_david_revision.sh` from the repository. The default output is `results/david_revision_20261014`.

## What runs

1. Verify the saved all-100 mL observations and likelihoods. Re-score the saved seven 60/40 candidate fits against their original split. Export calibration BIC, parameter counts, held-out likelihood, per-state and overall MRPE, and Experimental Noise (%). These are prediction errors on withheld cultures, not a relabeled training BIC. The split contains 19 calibration and 12 selection cultures. Each partition is summarized by condition, day and state, so both can contain 116 mean observations without sharing cultures.
2. Experimental RPM comparisons at 100 mL and 80,000 cells/mL. The culture is the experimental unit. Use days 1–4 averages and day 5 separately. Report two-sided Welch tests, unadjusted mean-difference confidence intervals and Holm-adjusted p-values over 24 contrasts (four states × two endpoints × three RPM pairs). Require complete measurements for each endpoint. Flag unmatched handling. Do not attribute centrifugation differences to RPM.
3. Pooled three-fold cross-validation of Models 0 and 6 over all 46 ESI-017 cultures (100 mL, 500 mL, 20% and 50% protocols). Whole cultures stay together. Includes training-only exclusion sensitivity and the original 100-to-500 mL reference, with and without exclusions. Eight scenarios × two models × 100 starts = 1,600 starts. Details: `docs/pooled_validation_workflow.md`.
4. Eleven exposure representations, each fitted on the original 60% calibration cultures and evaluated on the same withheld 40%. 100 hybrid starts per representation = 1,100. No all-data fit anchors these searches. The existing all-100 mL/500 mL sensitivity remains a separate result. This stage does not silently switch the main model or exposure representation.
5. Six routed ablations with **separate** cell/aggregate response parameters in the two both-equation variants, as David requested. 100 hybrid starts each = 600. Baseline and full Model 6 reuse the verified all-100 mL fits. The existing shared-parameter results remain intact. Outputs include actual calibration BIC boxplots, 500 mL MRPEs and noise. No 500 mL refitting here.
6. Both iPSC lines: unchanged Model 0 and Model 6 transfer predictions, plus five refitting strategies per line. Each refit has 100 hybrid starts = 1,000 total. All strategies use the original cell-line state-scaled weighted squared-error objective, the same data, and the final ESI-017 bounds. Seven-parameter refitting preserves the parameter set, not the old narrower bounds.
7. Updated figures and tables, including separate cell-line landscape figures to avoid seven overlapping curves. PDF, SVG, PNG, JPEG, EPS and TIFF are exported, with 600 dpi raster formats.
8. ProfileLikelihood.jl profiles for all eight ablation formulations, both all-active and with K_G, K_Lac and Y_GX fixed at each formulation's own all-100 mL estimates. Native package plotting is retained. These fixed-three profiles describe uncertainty with those values held fixed. Open profiles remain open if the data do not constrain them.

The default is **4,300 hybrid searches plus profiles**. It runs stages sequentially, with 42 Julia threads available within fitting stages. This avoids launching several 42-thread jobs simultaneously. Runtime depends on optimization and profiling. The driver is not an assurance of bounded profiles or a preferred ranking.

### Cell-line parameter sets

| Strategy | Parameters re-estimated |
|---|---|
| Baseline refit | All nine Model 0 parameters |
| Seven-parameter refit | mu_max, alpha, K_G, K_oG, Y_GLac, K_d, K_Agg |
| Five metabolic parameters | K_G, K_Lac, K_oG, Y_GX, Y_GLac |
| Hydrodynamic-only refit | K_tau, K_eps |
| All Model 6 parameters | All eleven parameters |

Other parameters retain the corresponding all-100 mL estimates. Both transfer comparisons retain all parameters. Refit errors are in-sample. They are not independent cell-line validation. The author confirmed dextran sulfate was used for 585A1 experimentally, superseding the raw CSV label. No explicit dextran-sulfate effect is fitted.

## Launch on the server

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
export PATH="$HOME/.juliaup/bin:$PATH"
julia +1.12.4 --version
mkdir -p logs
nohup env JULIA_NUM_THREADS=42 \
  bash scripts/run_david_revision.sh \
  > logs/david_revision.log 2>&1 < /dev/null &
echo $!
```

Default source folders:
- `results/calibration100_hybrid_lac100_kg100`
- `results/split6040_refined_hybrid_lac100_kg100`

Override only if the corresponding final results are elsewhere, using `HPSC_FINAL_FIT` and `HPSC_SPLIT_FIT`. The script checks required files before starting. It does not download server data.

## Progress and recovery

```bash
tail -f logs/david_revision.log
cat results/david_revision_20261014/progress.log
pgrep -af '[j]ulia'
```

Ctrl+C exits `tail` only. With nohup, disconnecting the laptop does not normally terminate the server process. Server shutdowns, memory limits or administrator actions still can.

If interrupted, rerun the **same launch command** with the same source files, settings and output folder. Completed hybrid searches and individual profiles are reused. Do not edit source code or pull changes during a run. Hash checks reject changed analysis inputs instead of mixing results. Incomplete searches are repeated. Do not run duplicate suites in the same folder.

A suite-level `COMPLETE` file appears only after the profile stage and all other scheduled stages finish. Per-stage logs are under `results/david_revision_20261014/logs`. The HTML report is generated before profiles, so it can be inspected while profiles continue.

## Writing outputs

- `selection/ranking.csv`, `selection/state_errors.csv`
- `exposure_selection/ranking.csv`, `exposure_selection/state_errors.csv`
- `rpm_statistics/welch_holm.csv`, `rpm_statistics/culture_endpoints.csv`
- `pooled_validation/replicate_error_summary.csv`, `fold_error_summary.csv`, `replicate_prediction_errors.csv`
- `separate_ablations/ranking.csv`, `state_errors.csv`, `figures/`, `profiles/`
- `cell_lines/fits.csv`, `state_errors.csv`, `parameter_policy.csv`, `curves.csv`
- `figures/` and `index.html`

Overall MRPE in the selection/cell-line tables pools eligible condition-day-state observations. Pooled CV additionally reports averages and sample variances across culture-level MRPEs. These estimands are explicitly separate. Experimental Noise is mean sample SD divided by the observed mean, expressed as a percentage. The code preserves diagnostics in its result files.

Oxygen CFD remains with Pouyan. Its results are not manufactured or approximated by this launch. Existing structural-identifiability, media-exchange and final model-profile results remain available and are not needlessly re-estimated. Results must be reviewed before claims and response text are finalized.

## Collect results for writing

After the suite finishes:

```bash
python3 scripts/collect_david_revision.py results/david_revision_20261014
```

Download `results/david_revision_20261014_writing.zip` from the VS Code Explorer. It contains tables, report, PDF/PNG/SVG figures and profile summaries, without the individual optimizer-run directories. The full result folder remains on the server.
