# Hydrodynamic input assignment comparison

This exploratory comparison uses all 100 mL calibration data and predicts the 500 mL condition without refitting. It does not replace the completed 60/40 model-selection exercise. The original four-version comparison is retained separately.

## Eight formulations

| Formulation | Cell-growth inhibition input | Aggregate-loss input | Active parameters | Conditional profile parameters |
|---|---|---|---:|---:|
| Model 0 | None | None | 9 | 6 |
| Shear: cell | Shear | None | 10 | 7 |
| Shear: aggregate | None | Shear | 10 | 7 |
| Shear: both | Shear | Shear | 10 | 7 |
| EDR: cell | EDR | None | 10 | 7 |
| EDR: aggregate | None | EDR | 10 | 7 |
| EDR: both | EDR | EDR | 10 | 7 |
| Model 6 | Shear | EDR | 11 | 8 |

Shear: both shares one K_tau between the two terms; EDR: both shares one K_epsilon. Model 6 retains separate shear and EDR parameters. The shared versions have 10 active parameters (7 conditional profiles). Inactive response parameters are excluded from estimation and BIC. Model 6 is stored as model5 in the existing calibration output; publication numbering is unchanged.

With aggregate size A, glucose G, lactate L and cell density X, the common equations are

\[
\mu_b=\mu_{\max}\frac{G}{G+K_G}\frac{K_{\mathrm{Lac}}}{L+K_{\mathrm{Lac}}}\frac{K_{\mathrm{Agg}}}{A^3+K_{\mathrm{Agg}}},\qquad
\mu=\mu_b\frac{K_X}{h_X+K_X},
\]
\[
\frac{dX}{dt}=\left(\mu-K_d e^{-t/12}\right)X,\qquad
\frac{dA}{dt}=\mu\left(\alpha-\frac{h_A}{K_{A,\mathrm{hyd}}}A\right).
\]

Time in these model equations is in hours. An absent cell input sets the inhibition factor to one; an absent aggregate input removes the hydrodynamic loss term. The metabolite equations and feeding events are unchanged. The table specifies h_X and h_A. The aggregate hydrodynamic response parameter is distinct from K_Agg in the growth expression.

“Cell” and “aggregate” identify direct input assignment. Cell-growth inhibition still affects aggregate growth through the shared growth rate, and state coupling can affect other trajectories. The swapped-input and both-term versions are alternative formulations, not simply deletions of terms from Model 6. Their fitted response parameters have units appropriate to their assigned exposure. Internal storage retains K_tau and K_eps; in shared both-term variants K_eps is tied to K_tau during every likelihood evaluation and is not independently estimated or profiled. Profile labels identify the shared exposure.

## Estimation and profiles

Six routed variants each receive 100 independent population searches followed by local refinement. The saved baseline and Model 6 fits and their original optimization distributions are reused. All variants use the same source data, likelihood, feeding assumptions and exposure lookup treatment. BIC uses each formulation's active parameter count. Optimization-run boxplots describe numerical repeatability, not sampling uncertainty.

Both profile suites use ProfileLikelihood.jl and its native plot_profiles function:

- **All-active:** profile every active parameter, reoptimizing the others.
- **Conditional:** fix K_G, K_Lac and Y_GX at that formulation's own saved all-100-mL estimates; profile the remaining parameters while reoptimizing the other unfixed parameters.

There are 16 suites in total. Fixed K_Lac retains the original lactate inhibition expression. Fixed estimates at bounds remain bound-constrained estimates. Conditional intervals do not propagate uncertainty in the fixed values. Unresolved intervals and unsuccessful optimizer stops remain reported. Inspect mle_summary.csv for any improvement over the source fit before using final profile results in the manuscript.

## Server execution

On the Mac, stage this workflow and push:

```bash
cd ~/Desktop/HPSCProliferation-R2
git add src/Calibration.jl ablations/RoutedExposure.jl ablations/test_routing.jl docs/ablation_eight_workflow.md figures/AblationPlots.jl figures/render_remaining_analyses.jl figures/render_routed_ablations.jl population/run_routed_ablations.jl profiles/run_routed.jl profile_figures/render.jl scripts/run_routed_ablations.sh scripts/run_routed_profiles.sh
git diff --cached --stat
git commit -m "Add eight exposure assignments and both profile likelihood modes"
git push origin main
```

On the server:

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
export PATH="$HOME/.juliaup/bin:$PATH"
mkdir -p logs
nohup env JULIA_NUM_THREADS=42 bash -c 'bash scripts/run_routed_ablations.sh && bash scripts/run_routed_profiles.sh' > logs/ablation_eight_shared.log 2>&1 < /dev/null &
```

This performs the six calibrations, assembles the eight-version comparison, then runs both profile suites for every version. It uses completed-run and completed-parameter checkpoints. Repeat the same command after interruption to resume; changed inputs or budgets require a new output directory. A currently interrupted optimization or parameter profile may need to restart, while completed checkpoints are reused.

Progress:

```bash
tail -n 30 logs/ablation_eight_shared.log
find results/ablation_eight_shared_hybrid_lac100_kg100 -mindepth 3 -maxdepth 3 -name COMPLETE | wc -l
```

The second command counts completed calibration starts (600 expected). Profile progress appears in the log afterward.

To regenerate the existing four-version figure without recalibration:

```bash
julia +1.12.4 --project=figures figures/render_remaining_analyses.jl ablation results/ablation_hybrid_lac100_kg100
```

If the source calibration was moved, set HPSC_ABLATION_SOURCE_FIT to its directory for this rendering command. Its figure_data/optimizer_bic.csv supplies the baseline and Model 6 optimization distributions.

## Outputs

Under results/ablation_eight_shared_hybrid_lac100_kg100:

- index.html and figures/ablation_eight_comparison.*: comparison figure, with actual BIC vertical boxplots and 500 mL state-specific MRPEs.
- ranking.csv, runs.csv, optimizer_bic.csv: fitted likelihoods, BIC and optimizer results.
- state_errors.csv: state-specific and overall prediction metrics and experimental noise.
- predictions.csv: predictions and corresponding experimental summaries.
- profiles/all-active/SCENARIO/figures/: all-active native package profile plots.
- profiles/conditional/SCENARIO/figures/: conditional native package profile plots.
- Each profile directory contains fixed_parameters.csv, mle_summary.csv and profile diagnostics.

Local workflow checks use deliberately short optimization budgets and selected profile parameters. They verify execution, equation limits, parameter accounting and rendering; their fitted values are not production scientific results.

The earlier separate-parameter outputs remain in results/ablation_eight_hybrid_lac100_kg100. Do not pull changed source into a currently running chained workflow: later stages can otherwise load different code. Stop or finish that workflow before updating. The shared comparison uses a new directory and does not overwrite earlier results.
