# Model numbering in the accepted manuscript

| Saved storage ID | Accepted manuscript model |
|---|---|
| 0 | 0 |
| 1 | 1 |
| 2 | 3 |
| 3 | 4 |
| 4 | 5 |
| 5 | 6 |
| 6 | 7 |

Original Model 2 duplicates Model 1 and is omitted. Model equations use the original manuscript identifiers internally, through `Workflow.legacy_model`.

`src/PublicationModels.jl` and `scripts/export_publication_tables.py` preserve an intermediate reporting convention used by historical workflows. Do not apply their labels to the final manuscript. `publication/` uses the accepted mapping above explicitly. Cached CSV identifiers remain unchanged; changing them would invalidate joins and stored source hashes. In LOCO folder names only, `loco_condition4` means manuscript Condition 5 and `loco_condition5` means manuscript Condition 6. The underlying `condition` column retains the original condition number.
