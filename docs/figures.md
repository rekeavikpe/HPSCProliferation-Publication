# Figure export specification

Checked 3 October 2026 against the journal's current author guidelines:
https://academic.oup.com/stcltm/pages/author-guidelines

- Submit each complete multipanel figure as one file, with panel letters at the upper left.
- Keep figure titles and full captions in the manuscript, not in the image.
- Raster photographs: at least 350 dpi at final physical size.
- Raster line art: at least 600 dpi; fine or complex line art: 1200 dpi.
- Vector graphics are preferred for plots; EPS and PDF are listed examples. Vector elements have no intrinsic dpi.
- JPEG/JPG figures are not accepted. We export a JPEG only for convenient viewing.
- Include alt text beneath each figure legend.

This code renders using CairoMakie with physical point dimensions (72 points per inch), then exports PDF/SVG and converts the PDF with pinned Ghostscript to EPS and 600-dpi PNG, TIFF and JPEG. For a 510-point-wide plot, 600 dpi gives 4250 pixels. Change the export dpi to 1200 for especially fine line art. Physical size and font size must be checked at the final layout size, not only by zooming into a large image.

EPS does not support PDF transparency in the same way. The light grids, SD bars and boxplot fills use transparency; Ghostscript handles conversion, which can rasterize affected content. Prefer the native PDF for vector submission. Check the EPS at final size before submission. The saved caption distinguishes biological SD from uncertainty in the prediction.

The journal limits ordinary Research Articles to seven main figures and tables combined; check article type and editor instructions before assigning the final figure/table layout. Extra analyses belong in the supplement.

CairoMakie documentation:
https://docs.makie.org/stable/explanations/backends/cairomakie
https://docs.makie.org/stable/how-to/match-figure-size-font-sizes-and-dpi

## Manuscript design audit

The current reference is Figures 2, 3 and 5 from the supplied R1 resubmission, cross-checked against the original Julia plotting code. The manuscript condition palette is blue, orange, green, magenta, sky blue and vermilion, with circle, square, diamond, hexagon, cross and x-cross markers respectively. Colors stay attached to their original experimental conditions when numbering changes. The 500-mL condition is excluded from the five-column calibration grid and now uses the requested dark teal in standalone plots.

The exported calibration grids preserve the original row order (cells, glucose, lactate, diameter), faint grid, black frame, condition-colored means/SD and curves, and framed bottom legend. The two baseline/hydrodynamic grids use identical axis limits. Arial is used locally; the server must have the font or receive an explicitly chosen substitute. The 500-mL comparison follows the manuscript's four-panel geometry with the requested dark-teal colour; a gray dashed baseline is added for this new comparison.

The BIC display preserves horizontal colored boxes and jittered optimizer points. Historical duplicate Model 2 is omitted; colors follow the surviving hypotheses. Unlike the old plotting helper, the new code does not inflate narrow boxes: quartiles remain their actual values. Only initial starts enter the boxes, with stopping flags; dependent refinements are saved in a separate CSV. The bundled snapshot has six initial runs; the 100-run command produces 100 initial runs per model. Full-range poor starts are retained; a separate best-achieved ΔBIC figure makes the ranking legible. Neither plot constitutes an uncertainty interval.

These are design-matched revision figures, not pixel-identical copies: the training cohort, condition numbering, numerical values and number of optimizer starts have changed. Final manuscript assembly must set panel letters and numbering consistently with the chosen figure grouping.

The renderer also writes `training_selection_zoom`: a companion focused view for distributions compressed by extreme optimizer outcomes. All runs still determine the quartiles and whiskers, and the full-range figure retains every point. The zoom reports the number of off-scale points; per-model counts are in `training_selection_zoom_counts.csv`. Its automatic upper bound uses the pooled 90th percentile and includes every model's best initial fit, with padding. Set `HPSC_BIC_ZOOM_MAX` to choose a different focus bound, and document it in the caption. Zooming changes neither scores nor selection.

## 4 October requested design revision

The controlling design references are the manuscript's actual Lagrangian figure generators in original `scripts/19_lagrangian_recalibration.jl` (`plot_bic_academic` and `forward_grid`), not just the generic plotting helper. Model 0 is restored to the top of the horizontal BIC display. The forward grid uses five columns, source-style combined line/marker legend and proportional typography. The original helper's artificial minimum box width is not copied: real quartiles must remain real. Final boxes require the new 100-start results.

The standalone 500-mL plots now use dark teal (#005A70), with the requested legend label “Experimental data.” Bars continue to represent between-vessel SD, stated in captions. Existing-result renders are design previews, not final results for the proposed 100-mM lactate cap. See publication_results_plan.html and config/publication_protocol.json for the remaining six-figure and supplementary work.

## Spacious calibration layout and BIC audit

The calibration master is now 900 × 680 points rather than 540 × 450, with 1-point curves, 5.2-point markers and 0.95-point SD bars with 5.5-point caps. Axis text is enlarged and panel gaps increased. These are large master figures; inspect text at the eventual journal placement size. The BIC distribution canvas is 720 × 440 points with smaller points; no BIC values or run exclusions have changed.

The original Lagrangian script used 100 time-limited differential-evolution searches (20 seconds each, `adaptive_de_rand_1_bin_radiuslimited`), with physical parameter bounds. Its success flag meant finite BIC, not confirmed optimizer convergence. The revised calibration uses normalized-log starts, bounded local least squares and Nelder–Mead polishing. Dataset and bounds also differ. Therefore the two 100-run distributions are not equivalent optimization experiments. A wider distribution cannot be attributed to plotting alone. Audit the server's individual outcomes before claiming a cause or changing the calibration algorithm.

```sh
python3 scripts/summarize_bic_runs.py results/calibration100_lac100_kg100/figure_data
```

This includes only the initial starts, reports convergence and objective spread, and leaves all points in the figures. The local spacious-layout preview uses the earlier six-start reference data, not the server's final 100-start results.

The subsequent label revision uses “Aggregate size (µm)” and single-line y-axis labels throughout the calibration grid, with 10-point y labels and an 11-point condition legend. The spacious master dimensions and thinner data/curve styling are retained.
