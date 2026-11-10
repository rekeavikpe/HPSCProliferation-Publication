# Controlled original-script comparison

This diagnostic uses the numerical definitions extracted verbatim from `original_fixed3.jl`, the supplied `21_profile_likelihood_lagrangian_fixKAgg.jl`. It retains the original CSV, original CFD workbook, Vern7 with 1e-8 tolerances, feed times, prescribed initial/fresh glucose, physical parameter coordinates, likelihood weighting and BOBYQA MLE starts (1, 0.5, 1.5, 2 times the original start). Full original source and input hashes are saved. Comments in that script contain stale free-parameter counts; the executable parameter vectors determine the dimensions.

The baseline includes the original K_G profile window of 2–30 mM. A second baseline uses the full 0.1–40 range. Subsequent single-factor branches compare against the full-range 100-mL-only original-data analysis: changed fixed values; release of K_Agg; K_G upper bound 400; wider shear/EDR parameter bounds; or lower-endpoint exposure clamping. The parent and changed factor are recorded in design.csv. Only K_G and mu_max are profiled in this targeted diagnostic; all other free parameters are reoptimized.

These are retrospective diagnostics, including deliberately reusing the original combined-data fixed values to reproduce the old analysis. They are not leakage-free validation results and must not replace the main 100-mL-only analysis. Do not compare raw NLL values across different datasets or observation counts. The original likelihood omits the Gaussian constant used by the revised likelihood.

A dataset branch here removes the 500-mL condition from the original dataset. It does not simultaneously substitute the subsequently curated dataset or change feeding assumptions. Those are separate differences. Native package plot titles remain subject to interval and convergence review; return-code recording does not prove global optimality.

## Running the comparison

From the repository root (Julia 1.12.4):

```sh
julia +1.12.4 --project=controlled -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --threads=42 --project=controlled controlled/run.jl results/controlled_comparison
julia +1.12.4 --project=profile_figures -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=profile_figures controlled/render.jl results/controlled_comparison
python3 controlled/report.py results/controlled_comparison
```

The default replacement fixed values are the saved local 100-mL-only MLEs included in `inputs/new_100mL_selected_parameters.csv`. To use the corresponding server fit, append `all results/calibration100/model5/selected_parameters.csv` to the fitting command, choosing a fresh output directory. A case name can replace `all` to run one case. Changed code or inputs require a fresh output directory; matching completed cases are skipped.

H, I and J replace only K_Lac, Y_GX or K_Agg respectively. They share B_100mL as their parent, so their effects are isolated. The original fixed values came from combined data and are used only for this diagnostic reconstruction. They must not be adopted as new 100-mL-only estimates.

The native plotting command exports PNG, JPEG, TIFF, PDF, SVG and EPS. The report checks computed threshold crossings independently of package confidence-interval titles. Saved profiles are checked against direct likelihood reevaluation. Optimizer return codes are recorded, but successful termination does not certify a global maximum.
