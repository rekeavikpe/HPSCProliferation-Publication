# Local workflow validation — 14 October 2026

These checks validate the implementation. They are not production scientific results and do not establish optimizer global convergence.

- 181 revision design/scoring checks passed: disjoint 19/12 culture split, unchanged exposure-scoring observations and weights, constant Eulerian inputs, Lagrangian P50 equivalence, cell-line parameter sets, feeding schedules, zero/10% error controls, and separate/shared response-parameter counts.
- 84 pooled-fold and training-only exclusion checks passed.
- 54 existing routed-ablation checks passed with separate response parameters.
- Short hybrid runs completed all eleven exposure scenarios and all five recalibration strategies for each cell line. Baseline and unchanged Model 6 transfer predictions were also exported.
- Exposure and cell-line checkpoint reloads completed without repeating saved searches. Source/settings fingerprints prevent mixing different analyses in one output folder.
- All 24 Welch p-values matched independent numerical integration of the Student t density within 1e-9. Holm adjustment was checked independently.
- All 70 cell-line MRPE and experimental-noise values matched calculations from exported trajectories and observations. Dense-curve and sampling-time predictions agreed within the numerical check tolerance.
- ProfileLikelihood.jl completed a small all-active shear-both profile and a fixed-three EDR-both profile. The native ProfileLikelihood.plot_profiles renderer completed with PDF, PNG and other exports.
- Revised heatmaps and both landscape cell-line figures were rendered and visually checked. Experimental colors, matching per-state y limits, panel letters and central legend boxes are retained.
- Shell syntax and the writing-bundle filter passed checks.

Short-run budgets were one population search, 100 evaluations and five local iterations. Profile smoke tests used three points per direction for a single parameter. Those outputs are stored under local `results/david_revision_smoke*` folders and must not be used in the manuscript. The final server driver defaults to 100 searches, 20,000 population evaluations, 3,000 local iterations and profile resolution 100.

The full 4,300-search production suite and all complete profiles have not been run locally. Existing reference fits are checked numerically when the driver starts on the server.
