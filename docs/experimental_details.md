# Experimental details confirmed by the author

Recorded 2 October 2026. These clarifications supplement the source-data audit without overwriting source workbooks or previous results.

## Authoritative Methods reconciliation

The author directs that the collaborators' manuscript Methods take precedence for experimental protocols. Checked `hPSC_Proliferation_Paper_v18_clean.tex`, lines 117–119. The author also confirms sampling before medium changes, consistent with line 119, which specifies daily 1-mL sampling immediately prior to exchange.

| Main condition | Exchange events before harvest | Harvest |
|---|---|---|
| 1: 100 mL, 20k/mL, 40 rpm | Day 4: 50% only | Day 6 |
| 2: 100 mL, 80k/mL, 40 rpm | Day 2: 50%; days 3–5: 80% | Day 6 |
| 3: 100 mL, 150k/mL, 40 rpm | Day 2: 50%; days 3–5: 80% | Day 6 |
| 4: 500 mL, 80k/mL, 30 rpm | Day 2: 50%; days 3–4: 80% | Day 5 |
| 5: 100 mL, 80k/mL, 30 rpm | Day 2: 50%; days 3–4: 80% | Day 5 |
| 6: 100 mL, 80k/mL, 60 rpm | Day 2: 50%; days 3–5: 80% | Day 6 |

No day-1 exchange is performed. Condition 1 is the important exception to the general schedule. Condition 6 additionally undergoes centrifugation on days 2 and 3 (400 rpm, approximately 36 rcf, 3 minutes). The media-exchange study changes the later 80% exchanges to 20% or 50%, retaining the day-2 50% exchange.

Predictions compared with observations on exchange days must explicitly use the pre-exchange state. The fresh implementation should enforce and test this event-side choice, rather than rely on interpolation defaults.

Methods specify 585A1 exchanges of 80% on days 2 and 3, and SCTi003-A exchanges of 50% on day 2 then 80% on days 3 and 4. The additional 585A1 workbook includes multiple schedule/treatment groups; map the Methods-described subset explicitly rather than apply one schedule to every workbook series.

The Methods also specify at least 80 measured aggregates per replicate. Individual aggregates must not be counted as independent bioreactor replicates. The Figure 3 caption at line 366 says n=3 for every condition; that caption must be reconciled with actual included replicate counts rather than used to discard additional cultures.

- Each labelled series represents a replicate in a separate bioreactor, followed across days. K1–K4 are four replicates of the 20k/mL condition. Experimental batch membership has not been established.
- The Condition 3/C1/day-5 glucose and lactate entries marked “No data” represent unmeasured observations. Fresh source imports must retain them as missing, not use the numerical replacements in the legacy CSV as measured observations.
- R1P1, R2P1 and R2P2 were excluded because glucose/lactate could not be measured. This explains the historical exclusion. It does not invalidate their recorded cell counts or aggregate sizes. Their use in state-specific analyses versus the historical complete-state calibration must be documented separately.
- The correct L29/day-6 aggregate value is unknown. Both 264.20876342953454 and 276.8101654846336 µm remain source candidates. Flag for collaborator clarification. Until resolved, keep the prior RPM-sheet value as an explicitly provisional continuity choice and evaluate the alternative in sensitivity analysis before finalizing affected conclusions. Do not average the candidates or call either verified.
- The L33/L34 versus L34/L35 header discrepancy and reuse of L35 across RPM conditions remain unresolved. Retain source sheet, condition and column as identifiers alongside the original labels. Do not silently merge series on the label alone or assert corrected identities.

## Feeding schedules transcribed from author-supplied screenshots

The first screenshot specifies 100-mL cultures starting at 8×10^4 cells/mL and 40 rpm. All three groups receive 50% exchange on day 2. On days 3, 4 and 5, the groups receive 20%, 50% or 80% exchange according to their group. Harvest is day 6. No day-1 exchange is shown. This screenshot alone does not establish schedules for other seeding densities or RPMs.

The second screenshot specifies 40×10^6 cells in 500 mL at 30 rpm, equivalent to 8×10^4 cells/mL. Exchange is 50% on day 2, then 80% on days 3 and 4, followed by harvest on day 5. No day-1 exchange is shown.

Neither screenshot specifies sampling order or exact clock times. Sampling order is now resolved by the manuscript Methods and author confirmation above: before exchange. Exact within-day clock times are not provided; conversion of protocol days to hours remains a documented model convention.

An experimental batch means cultures initiated together from the same starting cell preparation, typically on the same date. Separate vessels can belong to the same batch. Batch information, if recorded, helps distinguish vessel variation from shared starting-material or day effects; it must not be invented if unavailable.
