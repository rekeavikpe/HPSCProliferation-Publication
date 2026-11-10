# Release changes relative to Git commit cd7e985

- Added final accepted-paper cached results and reference figures, including nine Model 6 profiles and the two diagnostic profiles.
- Added portable versions of the final figure scripts that were previously outside the Git repository. Figure 3 uses all seven Condition 5 cultures, with available-measurement counts by day/state.
- Added reproducible score/likelihood checks, table exports, integrity hashes and an analysis map.
- Replaced the landing documentation with the accepted-paper workflow and corrected its model/condition numbering guide. Earlier guides are preserved under `docs/historical/`.
- Clarified the author-confirmed dextran-sulfate annotation for 585A1, and removed the stale no-DS label from new prediction exports. Original numerical CSVs and their hashes are unchanged.
- Replaced bundled third-party Neto source files with attributed, pinned source links and a checksum-verified download of the numerical table required by the oxygen check.
- Omitted unused historical cross-validation plot exports; all reported numerical cross-validation results and final paper figures are retained.
- Declared Apache 2.0 for code and the team-approved CC BY 4.0 license for data/results.

No model equations, calibration data, saved parameter values, profile values, cohort assignments or scientific conclusions were changed for packaging. The original Git working directory was not modified. The publication additions are recorded in the GitHub `v1.0.2` tag together with the archived source, data and result caches.

## Version 1.0.1

Citation-only update: added the accepted-paper and code/data citations to README.md, a BibTeX file, and the paper as the preferred citation in CITATION.cff. Updated version metadata. No code, input data, cached results, figures or numerical verification outputs changed from v1.0.0.

## Version 1.0.2

Added the author-reserved Zenodo DOI 10.5281/zenodo.24117723 to the README, software citation and deposit metadata. This DOI identifies the code/data record, not the journal article. No scientific code, numerical inputs, cached results, plots or analysis settings changed from v1.0.1.
