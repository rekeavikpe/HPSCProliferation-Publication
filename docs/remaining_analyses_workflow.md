# Remaining analyses

Run `bash scripts/run_remaining_analyses.sh` from the repository root. Julia 1.12.4 and the existing environments are used. The script restores the Juliaup launcher path and limits BLAS to one thread. `JULIA_NUM_THREADS` defaults to 42.

## Hydrodynamic ablations

The primary model selection remains unchanged. This is an all-100-mL sensitivity conditional on publication Model 6. Two variants are fitted with 100 independent differential-evolution searches (50 population members, 20,000 evaluations) followed by up to 3,000 Nelder–Mead iterations per start:

- Shear only: set EDR to zero in the exposure lookup; omit the now-inactive K_eps parameter from estimation.
- EDR only: set shear to zero; omit the now-inactive K_tau parameter from estimation.

Both variants retain all nine common biological parameters, including finite lactate inhibition, and therefore have ten free parameters. A fixed numerical placeholder for the inactive hydrodynamic coefficient cancels from the equations exactly. It is not a biological parameter estimate. The common bounds, observations, SD weights, feeds and initial conditions match the final all-100-mL calibration. No 500-mL data enter optimization or selection. The existing baseline and full fits are reused and evaluated with the same kernel; their parameter counts remain nine and eleven.

The output contains training BIC, actual optimizer status, boundaries, selected parameters, statewise/overall calibration and 500-mL MRPE, and Experimental Noise. Model removal measures conditional contribution within this formulation; it does not establish an independently measured causal effect. Results are not a replacement 60/40 selection exercise.

Per-start checkpoints are saved in `results/ablation_hybrid_lac100_kg100/{shear_only,edr_only}/run*/`. Restart with the same command: completed starts have their saved likelihood checked before reuse. Incomplete directories are archived before retry. Changes to inputs, source files or budgets require a new output directory. Linux flock prevents duplicate wrapper runs.

## Oxygen sensitivity

The short oxygen analysis runs first and writes `results/oxygen_sensitivity_final/`. Its main quantity is the kLa required to balance bulk oxygen demand at the maximum observed condition-day mean density. A second column uses the maximum individual measurement, making the peak definition explicit. The six source conditions include the 500-mL evaluation condition.

The reference oxygen uptake is 17.7 amol/cell/s from TC1133 aggregates at day 7, not an ESI-017 measurement. Uptake multipliers 0.5/1/2 and kLa 0.1–10 per hour are exploratory scenarios, not measured uncertainty intervals. Assumptions are editable in `config/oxygen_sensitivity.toml` and copied into every result folder.

The original Neto implementation needs correction: G is D*A/V, and Table 1 gives specific area 16.4/m at 100 mL, whereas the previous script used V/Vref and 40/m. The audit retains the inherited 0.05-m diameter and water properties explicitly as unverified assumptions. It reports both previous and corrected conditional kLa values, together with correlation-range checks. The 500-mL value of 1.8/h remains an assumption. Do not present those conditional kLa values as measured vessel properties. Culture geometry and medium/headspace applicability remain limitations.

Da = q*X/(kLa*Cstar) is a bulk-capacity ratio using the maximal driving force. It neither measures actual dissolved oxygen nor excludes intra-aggregate oxygen gradients. No automatic oxygen-replete classification is used.

Sources: [Neto et al., equations 6–7 and Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC9405225/); [Kim et al., section 3.3](https://pmc.ncbi.nlm.nih.gov/articles/PMC12024368/).

## Outputs

- `results/ablation_hybrid_lac100_kg100/index.html`
- `results/ablation_hybrid_lac100_kg100/ranking.csv`
- `results/ablation_hybrid_lac100_kg100/state_errors.csv`
- `results/ablation_hybrid_lac100_kg100/figures/ablation_comparison.*`
- `results/oxygen_sensitivity_final/index.html`
- `results/oxygen_sensitivity_final/condition_summary.csv`
- `results/oxygen_sensitivity_final/correlation_audit.csv`
- `results/oxygen_sensitivity_final/figures/oxygen_sensitivity.*`

Figures use CairoMakie and the existing manuscript style/export helper. Local short-budget workflow results under `results/ablation_workflow_check` are software checks only, not final fits. The six oxygen unit/threshold checks and 46 ablation-limit/count checks passed locally.
