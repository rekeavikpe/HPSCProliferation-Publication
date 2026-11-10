# Verification record

Local platform: Julia 1.12.4, Apple arm64, 10 logical CPUs. Checked 3 October 2026.

- 94 equation, feeding, optimizer, whole-culture splitting and split-feasibility checks passed.
- 56 reference checks passed: saved likelihood/BIC for all seven models and saved 500-mL predictions for the selected model and baseline.
- 42 short optimizer tasks were executed using one thread and again using 42 threads. All corresponding likelihood values were exactly equal. These are short optimization checks, not 42 independent full calibrations.
- The complete workflow was exercised with a deliberately short seven-model smoke fit, producing a selection lock and scale-transfer outputs. Those smoke-fit scores are not scientific results.
- CairoMakie rendered the four-state scale-transfer figure. PDF, EPS, SVG, PNG, TIFF and JPEG outputs exist. Raster images are 4250 × 3167 pixels with 600-dpi metadata (PNG rounds to 599.9988 dpi). The PNG was visually inspected: labels and error bars are legible; glucose/lactate feeding jumps are present.
- A fresh full calibration completed all 56 starts/refinements. Model 5 was selected again: BIC 1278.4606497690702, compared with 1278.460649793328 previously. The largest BIC discrepancy across models was 0.0001328 (Model 6, unconverged); model order was unchanged.
- All 140 checks of fresh saved stages, selected scores, and tighter-tolerance predictions passed. Across core, reference and full-run checks, 290 assertions passed.
- Selected starts for Models 2 and 6 did not meet their stopping criteria. These flags are retained; this is not a claim of convergence for all models.
- Seven experimental CSVs regenerated from the source workbooks were byte-identical to the curated inputs.

The lab server has not been accessed or benchmarked. A thread-count option does not guarantee linear acceleration. Separate statistical uncertainty from numerical reproducibility. Structural-identifiability proofs and the media/cell-line pipelines remain unfinished; these checks do not certify the manuscript as publication-ready.

The first reproduction run inherited stale optimizer-description fields from the earlier template. The actual seed (20261003) was recorded correctly in run.toml throughout; the effective configuration metadata was corrected afterward, with the original template preserved. No numerical inputs or fitted values were changed. The runner now writes its effective seed/iteration settings and input/source hashes directly, and that metadata path passed the subsequent smoke run.
# Structural and figure extension (3 October 2026)

All seven pooled lower-clamp global symbolic analyses completed. All seven exact finite-field local-rank checks completed at the configured initial states. See `reference/structural/summary.csv` and `docs/structural_identifiability.md` for assumptions and finite ambiguities; these are not practical-identifiability results.

`structural/test.jl`: 840 symbolic-versus-numerical RHS comparisons passed across all six CFD segments and all five training conditions; 2 identifiable/confounded control checks passed. These supplement the earlier numerical regression checks below.

`scripts/verify_figure_data.jl`: 488 checks passed comparing plotted sampling-time values with the numerical predictions and testing both sides of every training medium-exchange event. The single-40-rpm structural design check correctly retains the separate K_tau/mu_max confounding in Model 5's constant-exposure segment.

The manuscript-style renderer was run on the saved fresh reproduction data. Baseline/hydrodynamic training grids, full-start BIC boxplots, best-score BIC, and 500-mL predictions were exported as PDF, SVG, EPS, 600-dpi PNG/TIFF and JPEG convenience copies. Generated raster previews were inspected for layout, legends and exchange jumps. Scientific values were not altered to match the manuscript's historical results. This does not validate the remaining cell-line/media analyses or the lab-server environment.
