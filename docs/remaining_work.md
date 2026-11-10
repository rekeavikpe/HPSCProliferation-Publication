# Revision analysis status — 5 October 2026

Model 6 is selected by the completed retrospective 60/40 culture split. The refined split training BIC is 1272.37149246 and selection NLL is 703.22003967. Its advantage over Model 1 on selection NLL is small (0.76022). Split Models 1, 4, 5 and 6 retain unmet stopping flags; repeated refinement stabilized scores but did not resolve those flags. No further optimization is scheduled solely to change those flags.

The matching completed all-100-mL Model 6 fit can be reused. Current Julia code reproduced its saved predictions within 1.1e-9 experimental sigma; the all-data selected fit meets its stopping criterion. Baseline predictions were also reproduced. Do not confuse this all-data convergence result with the split status.

The final ProfileLikelihood.jl suite fixes KG=21.267814603918406 mM, KLac=100.00000000000004 mM and YGX=2.258385491255778e-6 micromol/cell at that fit. Eight saved conditional intervals have both crossings and no unconfirmed-point flags. Their interpretation is conditional on the fixed values. KLac remains a constrained estimate at an operational upper bound, not an independently identified biological constant. All 11 estimated parameters still count in selection BIC.

Final 500-mL, media-exchange and both cell-line transfer/recalibration outputs exist. Parameter hashes agree with the all-data fit, and MRPE/experimental noise have been recomputed. Cell-line recalibration is in-sample and uses the original weighted-SSR objective. It must not be presented as independent validation.

## Still needed before submission

1. Download the final refined split state_errors.csv and checkpoints for the archival result bundle. The last supplied statewise MRPEs predate refinement; do not silently relabel them.
2. Complete separate shear-only and EDR-only ablations if making claims about their individual contributions. Model 0 versus Model 6 establishes only the combined comparison.
3. Audit the six-condition oxygen calculation, correlation and units, then complete uptake/kLa sensitivity. Bulk oxygen calculations cannot exclude intra-aggregate gradients.
4. Reconcile the no-dextran-sulfate 585A1 source cohort with manuscript Methods. Retain experimental-label sensitivities and unknown-batch limitations in aggregate statistics.
5. Report exposure-summary sensitivity as an all-data comparison conditional on Model 6. Eulerian maximum beats Lagrangian P50; neither P50 superiority nor a new held-out exposure-selection study is established. Below-200-um extrapolation is a four-start exploratory sensitivity, with unresolved power-law convergence.
6. Use the refreshed final-parameter vessel contrast: observed day-5 100-minus-500 aggregate-size difference +68.985 um; Model 6 -3.238 um. The empirical aggregate tests have unchanged source hashes and can be reused.
7. Synchronize figures, equations, bibliography, parameter units and publication numbering. Archive final code/configuration/parameter hashes and server outputs. Structural results remain conditional on their documented observation and segment assumptions.

See reviewer_action_matrix.csv for comment-level actions and config/publication_protocol.json for the current protocol. Numerical files retain storage IDs; publication Model 6 is stored as model5.
