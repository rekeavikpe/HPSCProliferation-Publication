# Accepted-paper analysis map

All cache paths below are relative to `publication/cache/`. The table text in `paper_tables.json` is a comparison reference, not an input to fitting.

| Paper item | Numerical source | Reproduction/check |
|---|---|---|
| Table 1; Figure 2; Table S6 | `training60/reference_fit/`; `october20/reuse/candidate_selection_optimizer_bic.csv` | `verify_primary.jl`; `plotters/figure2.py`; `report_tables.py` |
| Table 2; Figure S5 | Primary Model 6 estimate; `nine_profiles/` | `verify_nine_profiles.jl`; `plotters/nine_profiles.jl` |
| Table 3; Figure 3 | `training60/reference_fit/`; `data/processed/measurements.csv` | `verify_primary.jl`; `plotters/figure3.py` recomputes all marker means/SDs from current cultures |
| Table 3; Figure 4 | `training60/media/` | `scripts/predict_media_exchange.jl`; `plotters/figure4.py` |
| Table S1 | `training60/reference_fit/split.csv`; `data/processed/series.csv` and other measurements | `report_tables.py` |
| Tables S2/S3; Figures S1–S3 | `data/processed/exposure_summary_inputs.csv`; `october20/reuse/exposure_rpm_sweep.csv` | `plotters/exposure_correlation.py`, `plotters/exposure_sweeps.py`; `report_tables.py` |
| Table S4; statistical comparisons | `october20/statistics/`; current measurements | `revision/october_statistics.jl`; `report_tables.py` |
| Table S5; Figure S4 numerical values | `october20/oxygen/`; `october20/oxygen_reference/` | `scripts/oxygen_sensitivity.jl`, `scripts/october_oxygen.py`; `report_tables.py` |
| Table S7 | `october20/p75/` | `verify_additional.jl` replays all P75 scores |
| Table S8 | Selection subset of primary Model 6 | `verify_primary.jl`; `report_tables.py` |
| Table S9 | `cross_validation/loco_condition*/` | `verify_cross_validation.jl` |
| Table S10; iPSC direct predictions | `training60/cell_transfer/`, `october20/cells/*/round1/`, `ipsc/*/fixed_K_d_Y_GLac/` | `scripts/predict_cell_lines.jl`; `verify_additional.jl` |
| Table S11; Figures S7/S8 | `ipsc/*/fixed_K_d_Y_GLac/` including `profiles/` | `verify_additional.jl`; `plotters/ipsc_fits.py`, `plotters/ipsc_profiles.py` |
| Figure S6 | `two_diagnostics/` | `verify_two_diagnostics.jl`; `plotters/two_diagnostics.jl` |
| Figure 1 | Coauthor-supplied CFD artwork | Reference TIFF provided; no CFD simulation project supplied |

`profile_seeds_historical/` contains numerical starting points for optional profile reruns. Its intervals are superseded by `nine_profiles/` and `two_diagnostics/`; do not report its earlier K_G interval.

The P50/P75 mixed-exposure study, aggregate-day-1 fits and eight-formulation ablations are not part of this accepted-paper result set. Historical algorithms for them remain in the source tree, but they are not invoked by `publication/reproduce.py`.
