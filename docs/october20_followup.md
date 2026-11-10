# Saved day-1 predictions and iPSC parameter reductions

Run `bash scripts/run_october20_followup.sh` after the original suite completes.
Default output: `results/david_october20_followup`. Original results are read only.

Day-1 evaluation uses the saved eight parameter vectors, retaining partition-specific
initial states, absolute culture age, original feeding times and post-day-1 scoring.
Exports training, selection and day-1-initialized 500 mL trajectories, observations,
state-specific errors and day-5/final-day fold expansion. Dense trajectories are
checked against the scoring solver. Figures show all eight formulations and a
readable three-model comparison (Model 0, EDR on aggregate size, Model 6).
This stage performs no optimization.

For 585A1, compare every nonempty subset of fixing K_d, Y_GLac and alpha
at the ESI-017 reference values (7 candidates). For SCTi003-A compare every
nonempty subset of fixing K_d and Y_GLac (3 candidates). This design is derived
from the original unbounded profiles, with no sequential dropping order.
Previously bounded parameters remain free in this first reduction comparison.
Each candidate has 100 hybrid starts and profiles of every remaining parameter.
All three cultures per line are used, as in the approved original line fits.
The six parameters already fixed stay fixed. Fixing does not delete model terms.

The existing five-parameter results are retained as the reference. Report fit
likelihood, fit BIC (with the existing cell-line weighting), per-state MRPE and
profile status for every candidate. This is not held-out selection or validation.
No candidate is automatically accepted. Choose only after considering boundedness,
loss of fit and biological justification for fixing parameters at ESI-017 values.
If none has bounded profiles, report that outcome rather than claiming success
or broadening bounds merely to make profiles cross.

Identical reruns resume completed candidate fits and profiles. The follow-up
manifest rejects changed source, data, reference estimates or budgets. Use a new
output directory if any of those change. Do not remove the old manifest.

Progress: `tail -f logs/david_october20_followup.log` and
`cat results/david_october20_followup/progress.log`.

Collect when complete:

```
python3 scripts/collect_october20_followup.py
```

The review archive omits per-start directories and binary profile caches. Keep
them on the server for the separate final reproducibility archive.
