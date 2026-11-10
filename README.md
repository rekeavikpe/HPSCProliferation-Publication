# hPSC proliferation in vertical-wheel bioreactors — publication release 1.0.2

Code, curated data and cached numerical results supporting **In Silico Modeling of Shear Stress and Energy Dissipation Rate Effects on Human Pluripotent Stem Cell Proliferation in Vertical-Wheel Bioreactors**, accepted in *Stem Cells Translational Medicine* (SCTM-26-0107.R2).

## How to cite

**Code/data DOI:** [10.5281/zenodo.24117723](https://doi.org/10.5281/zenodo.24117723)

If you use the models, code, data or results in this repository, please cite **both the associated paper and the version of the code/data release you used**.

### Associated paper

Avikpe, F. R., Alibhai, F. J., Romero, D. A., Mostofinejad, A., Bauer, J. E. S., Montague, E. C., Laflamme, M. A., & Amon, C. H. (in press). In Silico Modeling of Shear Stress and Energy Dissipation Rate Effects on Human Pluripotent Stem Cell Proliferation in Vertical-Wheel Bioreactors. *Stem Cells Translational Medicine*. Accepted for publication.

### Code and data

Avikpe, F. R., Alibhai, F. J., Romero, D. A., Mostofinejad, A., Bauer, J. E. S., Montague, E. C., Laflamme, M. A., & Amon, C. H. (2026). *Code, data and cached results for In Silico Modeling of Shear Stress and Energy Dissipation Rate Effects on Human Pluripotent Stem Cell Proliferation in Vertical-Wheel Bioreactors* (Version 1.0.2) [Computer software and data]. Zenodo. https://doi.org/10.5281/zenodo.24117723

The DOI above identifies this code/data record. The accepted-paper citation is separate; the journal article DOI will be added when available.

Bibliography files are provided in [CITATION.bib](CITATION.bib) and [CITATION.cff](CITATION.cff). GitHub’s “Cite this repository” option provides the associated-paper citation.

## Start here

This release contains the final 19-culture calibration, 12-culture selection, frozen 500 mL and feeding-protocol predictions, P75 sensitivity, five leave-one-condition-out folds, iPSC recalibrations, nine-parameter Model 6 profiles, and the two fixed-parameter diagnostic profiles. It includes the corrected Figure 3 observations and final plotting scripts.

The biological equations, saved estimates, cohort assignments and likelihood definition have not been changed for this release. The code is based on Git commit `cd7e98567ee3be0e4b3884982d3dbe83246b7f0d`, with publication packaging and portable plotting/verification entry points added. See `publication/RELEASE_CHANGES.md`. The matching GitHub publication release is tagged [`v1.0.2`](https://github.com/rekeavikpe/HPSCProliferation-R2/tree/v1.0.2); this tag fixes the source, data and caches used by the Zenodo ZIP.

**Use `publication/` for the accepted paper.** Older workflows and their documentation are retained for research history. They include superseded model labels, all-data calibrations and exploratory analyses that were not included in the accepted paper. They are not the release entry points. In particular, the accepted primary model was never refitted to all 31 small-vessel cultures after selection.

## Quick verification and reproduction

Requirements: Julia **1.12.4** and Python **3.12**. Julia environments are pinned by the included `Manifest.toml` files; Python plotting dependencies are in `publication/requirements.txt`.

From the extracted release root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r publication/requirements.txt
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=population -e 'using Pkg; Pkg.instantiate()'
julia --project=profiles -e 'using Pkg; Pkg.instantiate()'
julia --project=profile_figures -e 'using Pkg; Pkg.instantiate()'
julia --project=revision/statistics -e 'using Pkg; Pkg.instantiate()'
python publication/reproduce.py
```

This verifies file integrity, replays the saved fits and profile likelihoods through the model, checks frozen predictions, and regenerates reporting tables. It does **not** run optimization. First-time package installation/precompilation can take longer than subsequent runs. The oxygen check downloads one checksum-verified numerical table from Neto’s pinned public repository; internet access is required for that first retrieval. No credentials are required. Upstream files without a declared license are not relicensed or bundled.

To also regenerate the plots:

```bash
python publication/render_figures.py --profiles
```

The original S7/S8 plots used Matplotlib 3.10.8, whereas the other Python figures used 3.11.2. To reproduce their exact original text layout, install `publication/requirements_ipsc.txt` in a second environment and pass its Python path as `--ipsc-python /path/to/that/python`. A single environment preserves the numerical curves but can change text spacing.

For the paper's Arial typography, install Arial locally. The font is not redistributed. For Julia plots, `HPSC_FONT` may specify the installed Arial font file; otherwise the font name `Arial` is used. A substitute font can change appearance but not numerical results. `HPSC_JULIA` can name a Julia executable if `julia` is not on PATH.

## What is where

| Folder/file | Contents |
|---|---|
| `src/` | Model equations, feeding events, likelihood and calibration utilities |
| `data/processed/` | Curated run-level measurements, culture identifiers, feeding schedules, and exported CFD exposure summaries |
| `population/`, `profiles/`, `revision/` | Parameter estimation, selection, profile algorithms and analysis routines |
| `publication/cache/` | Saved parameters, predictions, observations, profile points, nuisance parameters and diagnostics |
| `publication/plotters/` | Portable final figure scripts |
| `publication/paper_figures/` | Figures from the final submission package, retained as reference artwork |
| `publication/paper_tables.json` | Table captions and rounded values transcribed programmatically from the final submitted table files, for comparison |
| `publication/release_verification/` | Verification results recorded when this archive was built |
| `publication/generated/` | Regenerated outputs; safe to recreate |
| `publication/ANALYSIS_MAP.md` | Maps every figure/table to the numerical inputs and code |
| `publication/RERUN_ANALYSES.md` | Optional optimization commands; not needed to use the saved results |
| `SHA256SUMS.txt` | Integrity manifest for release files |

## Scientific conventions

- `cells` means **viable cells/mL**, not total cells. Day-zero nominal inocula and censored/missing observations remain explicitly identified in the curated data.
- Complete cultures, not individual observations, are assigned to calibration (19) or selection (12). Condition 4 is the nine-culture **500 mL** dataset. Conditions 1, 2, 3, 5 and 6 are 100 mL. Condition 5 has seven cultures; one ends on day 4, so six contribute day-5 measurements. Metabolite availability is lower for this condition.
- Storage IDs `0,1,2,3,4,5,6` correspond to manuscript Models `0,1,3,4,5,6,7`. Thus `model5` is manuscript **Model 6**. Some original CSV fields called `model` or `publication_model` use earlier numbering. Consult `publication/model_numbers.csv` and the analysis map before joining files. No cache identifiers have been silently renamed.
- Selection likelihoods use the fixed calibration noise scales, as implemented in `population/SplitSelection.jl`. Selection SDs are descriptive; they do not replace these scales. The full Gaussian likelihood includes its normalization terms. BIC uses `2*NLL + k*log(n)`.
- MRPE is `100*mean(abs(prediction-observed_mean)/abs(observed_mean))` over eligible condition–day–state means. Overall MRPE pools these means. Experimental variability is `100*mean(sample_SD/abs(observed_mean))`, without likelihood floors. Initial values are not scored. Missing measurements are not imputed.
- Model 6's nine confidence intervals are conditional on `K_Lac` and `Y_GX` fixed at the calibration estimates; `K_G` is profiled. The additional two diagnostic profiles let the other ten parameters vary. Their unresolved limits are reported as unresolved within the tested bounds, not as finite confidence limits.
- The iPSC recalibrations are fits to the three runs of each line, not independent validations. Recalibration 1 fits five parameters; Recalibration 2 fits three. The 585A1 cultures contained dextran sulfate; its effect is not an explicit model term.
- The 500 mL oxygen-transfer coefficient is the analytical scaling estimate **1.179809899 h⁻¹**, not a measurement or near-surface CFD calculation. `cache/october20/oxygen/` contains this final estimate. `oxygen_reference/` retains the earlier correlation reproduction, whose inherited 500 mL value of 1.8 h⁻¹ is not the accepted-paper estimate.

## Scope of reproducibility

Saved fits, predictions, intervals and plots can be checked without rerunning optimization. Selected fits and score distributions are included; not every historical optimizer trajectory/checkpoint is included. The full fitting algorithms, fixed random seeds and settings are provided. Nonlinear optimization can depend on platform and numerical libraries; identical final iterates are not guaranteed on a new system.

The package starts from the supplied CFD summary exports. Raw CFD simulation projects and particle trajectory TXT exports were not available for this release. The spreadsheet rebuild itself cannot be independently reconstructed from these summaries. Figure 1 (CFD artwork) and Figure S4 (Damköhler ternary artwork) were supplied by coauthors; their original drawing projects are not included. Their final image files are provided, and the oxygen calculations and tabular inputs are reproducible. These limits do not require an external workbook to run the biological model.

Original settings and source records may contain historical computer paths. They are documentation of the saved run; the publication verification and figure scripts resolve inputs inside this archive. Some historical fitting/import scripts do require original source locations, as their documentation indicates.

## License

Code: Apache License 2.0 (`LICENSE`). Curated data, cached results, figures and publication reporting outputs: CC BY 4.0 (`LICENSE-DATA.txt`). These apply to different parts of the archive; third-party dependencies retain their own licenses. See “How to cite” above for the paper and code/data references.
