# Curated data used by the publication

Numerical CSVs preserve the exact inputs and source hashes used for the saved analyses. Cells are viable cells/mL. `series.csv` and `measurements.csv` retain manuscript condition IDs (Condition 4 is 500 mL).

## Correction to historical dextran-sulfate metadata

The `dextran_sulfate=absent` flag in the original `cell_line_measurements.csv` import and the no-DS wording in `cell_line_source.json` are superseded for **585A1**: the author confirmed that those cultures contained dextran sulfate. The final manuscript reports this. `publication/culture_metadata_corrections.csv` supplies the corrected annotation. This correction changes neither the measured values nor culture selection. The model has no explicit dextran-sulfate term. Original numerical input files remain byte-for-byte intact to preserve saved-run checksums; use the corrected annotation when describing the cultures. Historical cached cohort labels containing “no DS” do not establish absence of dextran sulfate experimentally.

The exposure inputs are summaries from the rebuilt particle workbook and the supplied Eulerian summaries. The raw CFD projects and trajectory TXT exports are not included. Source workbook names, cells and hashes are retained in the CSV/JSON records; source workbooks are not required to run the model.

License: CC BY 4.0 for the data; see ../../LICENSE-DATA.txt.
