# Model numbering for the revised manuscript

| Publication model | Original model | Saved storage ID |
|---|---|---|
| 0 | 0 (baseline) | 0 |
| 1 | 1 | 1 |
| 2 | 3 | 2 |
| 3 | 4 | 3 |
| 4 | 5 | 4 |
| 5 | 7 (dual interaction) | 6 |
| 6 | 6 (growth-coupled EDR) | 5 |

Original Model 2 is redundant with original Model 1 and is omitted. Original Model 6 retains its manuscript number. This is a labeling change only: equations, likelihoods, parameter estimates, and predictions are unchanged.

Saved results, checkpoint folders, and numerical CLI arguments retain their storage IDs to allow interrupted runs to resume safely. Thus `model5` and CLI argument `5` refer to publication Model 6. Existing CSV fields named `model` or `publication_model` in archived computational outputs also retain the earlier storage numbering. Do not manually rename these folders or change fitted-result identifiers.

Updated renderers translate storage IDs into publication numbers. To export explicitly numbered publication tables without changing source data, run:

```bash
python3 scripts/export_publication_tables.py results/calibration100_hybrid_lac100_kg100 results/publication_tables/calibration
```

The exporter copies CSV tables recursively, renames the old ID field to `storage_model`, and adds `publication_model`. These are reporting tables, not inputs for fitting or figure-generation scripts. Renderers continue to read the original computational outputs.
