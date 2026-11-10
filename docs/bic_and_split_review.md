# BIC display and data-split audit (5 October 2026)

The completed `calibration100_hybrid_lac100_kg100` run used all eligible cultures from the five 100-mL conditions. `population/run.jl` calls `context(cfg=cfg)` without a replicate-ID subset; `Workflow.context` includes all series in conditions 1,2,3,5,6. The 116 likelihood entries are condition–day–state means. The 100 searches per model are optimizer repetitions, not random data splits. Existing figures must not be described as 60/40 results.

## Why the original BIC plot differs

The supplied original `scripts/19_lagrangian_recalibration.jl`, function `plot_bic_academic`, sets `minw = (bic_max - bic_min) * 0.02` and `aq3 = max(q3, q1 + minw)`. Thus the displayed rectangle can be wider than the true interquartile range. Current figures display actual quartiles. Narrow boxes reflect closely agreeing optimizer outcomes.

That original script also calibrates the combined six-condition dataset, minimizes weighted SSR, and computes a separate residual-variance-based BIC score. The revised pipeline calibrates only 100 mL and optimizes its fixed-sigma Gaussian likelihood directly. Bounds, noise floors and search/refinement settings also differ. Removing a redundant model is therefore not the only change, and absolute BICs across these pipelines are not directly comparable.

## Combined figure

Run from the repository root:

```bash
julia +1.12.4 --project=figures figures/render_bic_scale.jl results/calibration100_hybrid_lac100_kg100/figure_data results/calibration100_hybrid_lac100_kg100/combined_figure
```

The left panel contains the six hydrodynamic models (publication IDs 1–6); the existing seven-model figure retains baseline Model 0. The right panels compare Model 6 and Model 0 predictions with 500-mL measurements. Exports include PDF, EPS, SVG and 600-dpi PNG, JPEG and TIFF. `bic_box_statistics.csv` gives the actual quartiles. This figure currently represents all-100-mL calibration.

## Required workflow for the requested 60/40 analysis

This split analysis has not yet been run. Lock a seeded random assignment of whole bioreactor cultures within each condition before fitting; every day and state from a culture stays in one partition. Account for missing metabolite coverage when defining strata, without using outcome values or fit quality to choose the split. There are 31 recorded 100-mL series, so a near-60/40 allocation is 19 calibration and 12 selection cultures. Audit culture identities and partition coverage before freezing the manifest.

Fit every candidate on the calibration partition using the same starts and settings. Derive observation means, noise weights and calibration initial conditions from that partition only. Evaluate the saved fits on the selection partition without reoptimizing. Report training BIC separately from predictive selection scores; BIC is based on maximized training likelihood and is not a holdout score. Declare the primary predictive selection criterion before obtaining split results.

After model selection, refit the chosen model on all 100-mL data, derive fixed parameter values from that fit, and perform conditional profiles and the external evaluations. The already-inspected 500-mL data remain a retrospective transfer evaluation, not a newly untouched test set. Preserve the previous all-data results and their provenance.
