# Final results after 60/40 selection

Run `bash scripts/run_final_results.sh` after pulling the updated repository. This is a reporting and verification job, not another calibration. It does not run ablations or oxygen sensitivity.

Default inputs under `results/`:

- `split6040_refined_hybrid_lac100_kg100`: completed 60/40 selection and refinement.
- `calibration100_hybrid_lac100_kg100`: matching all-100-mL final fit.
- `active_fixed_three_hybrid_lac100_kg100`: native ProfileLikelihood.jl results.
- `media_exchange_hybrid_lac100_kg100`.
- `cell_lines_frozen_hybrid_lac100_kg100`.
- `cell_lines_recalibrated_hybrid_lac100_kg100`.

The script verifies disjoint culture assignments, 19/12 counts, selected-model parameter hash, refined split NLL and training BIC from saved predictions, and recomputes split MRPE/Experimental Noise. It verifies all-data parameter hashes used by the profile and external-evaluation suites, and that the three fixed profile values match the all-data fit. A changed model or parameter provenance stops assembly.

The final output is `results/final_publication_results/index.html`. It links to publication-numbered parameter/metric tables and regenerated figures. Model 6 is stored as model5; reporting labels use Model 6. Five calibration conditions are numbered 1–5; original data identifiers stay intact.

The combined figure uses split-training BIC on the left and predictions from the all-data fit on the right. Split held-out NLL determines model selection. Supplementary convergence information remains attached. The all-data calibration renderer runs with HPSC_FIT_FIGURES_ONLY=1 so older all-data BIC figures are not generated in this collection. Conditional profiles are rendered from saved ProfileLikelihood.jl data without optimization. Media and cell-line plots use verified copies; numerical source results are preserved. The fit's derived figure_data directory is refreshed.

The final COMPLETE marker is written only after all renderers succeed. TABLES_VERIFIED means numerical tables passed but does not mean figure generation finished. Rerunning the same command refreshes reporting artifacts, not fits.

Tables:

- tables/split_selection.csv: refined held-out likelihood, training BIC and convergence.
- tables/split_metrics.csv: refined training/selection statewise and overall MRPE/Experimental Noise.
- tables/evaluation_metrics.csv: all-data calibration, 500-mL, media exchange, cell-line transfer and recalibration.
- tables/all100_parameters.csv, tables/fixed_parameters.csv, tables/conditional_intervals.csv.

Figures: combined/, calibration/, profiles/, media/figures/, cell_recalibration/figures/. Source hashes are in source_manifest.csv. The profiles condition on data-derived fixed values; their uncertainty is not propagated. Cell-line recalibration is in-sample and the 585A1 no-DS Methods discrepancy remains explicitly flagged.

Local verification uses the short-budget split fixture with the uploaded final all-data bundle. It is a workflow check, not a substitute for assembling the production split results on the server.

The separate combined_cfd/cfd_bic_and_500mL figure uses three columns, retaining one-third of the width for CFD, 26% for BIC and the remaining width for predictions: the original CFD image (a–b) on the left, split BIC (c) in the centre, and four 500-mL state panels (d–g) in a 2 × 2 grid on the right. The ordinary combined figure is retained. Both cell-line orientations have a–h labels in reading order. The embedded CFD JPEG retains its original resolution; exporting at higher DPI does not add CFD detail.
