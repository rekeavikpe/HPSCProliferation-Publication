# Structural identifiability

The repository includes and has executed a structural analysis for all seven public models (0–6). This is distinct from practical identifiability with sparse, noisy measurements.

## Equations and experimental design

`structural/models.jl` builds a pooled system for the five 100-mL conditions with shared kinetic parameters and their actual 30-, 40- and 60-rpm P50 CFD values. Four biological states are treated as continuously and exactly observed. It preserves the implementation's rational shear response, Hill response, metabolic-growth distinction, aggregate equations and coupled shear/EDR response. The exponential death transient is represented by a known clock, E′ = −E/12, with E(0) = 1 in the specialized analysis.

CFD values are parsed as exact decimal rationals. The code supports the lower endpoint clamp, all four affine interpolation segments, and the upper clamp. Shear/EDR are functions of predicted diameter within a segment, not freely varying independent inputs. Numerical evaluation floors and medium-exchange resets are outside these rational segment equations. The analyzed early-time segment is valid at the configured 15-µm initial diameter, before the first exchange. Positive trajectories remain in this segment for a nonzero initial interval.

The symbolic equations were compared with the actual numerical fitting kernel at all six exposure segments and all five conditions for every model: 840 derivative comparisons passed. Identifiable and confounded control systems also passed (2 checks).

## Completed results

StructuralIdentifiability.jl's global analysis was run for the pooled lower-clamp segment with generic initial states and a probability threshold of 0.999 per analysis. It returned:

| Model | Parameters | Globally identifiable in generic segment analysis | Locally, but not globally, identifiable |
|---|---:|---:|---|
| 0 | 9 | 9 | None |
| 1 | 12 | 12 | None |
| 2 | 12 | 12 | None |
| 3 | 12 | 12 | None |
| 4 | 12 | 11 | K_tau |
| 5 | 13 | 10 | K_tau, gamma, mu_max |
| 6 | 11 | 11 | None |

Model 4 uses K_tau squared, giving a sign symmetry over the unrestricted algebraic parameter domain. The fitted positive domain excludes the negative sign counterpart. Model 5's finite ambiguities require further investigation before any positive-domain global uniqueness claim. These generic-segment classifications are not a global-identifiability proof for the entire hybrid, piecewise model at specialized initial states.

A separate exact finite-field sensitivity-rank calculation fixes the actual nominal initial cell densities, condition-specific measured initial glucose means, zero lactate, diameter 15 µm and clock E = 1. All seven parameter Jacobians have full row rank: 9, 12, 12, 12, 12, 13 and 11 respectively. This is a nonzero-minor witness for generic **local** identifiability at those fixed initial states, conditional on ideal observation of all four states and the specified pooled design. The saved finite-field matrices, parameter witnesses, initial states, prime and precision make the calculation reproducible. No numerical singular-value threshold is used. A failed random rank check would be inconclusive, not proof of nonidentifiability.

An identifiable early-time prefix is sufficient for local identifiability of the complete observation record under those ideal assumptions; it does not require claiming global identifiability across unknown segment transitions. The result does not imply that daily biological observations determine each parameter precisely. Boundary estimates and flat profile likelihoods must still be reported. It also does not establish identifiability of a single RPM condition or of a dataset lacking glucose/lactate observations.

As a design check, Model 6 was also analyzed at only 40 rpm in the constant lower-clamp segment. K_tau and mu_max were nonidentifiable separately, unlike the pooled result. At one constant shear value they enter growth through a single product. This directly illustrates why the different operating conditions matter to the structural conclusion. It is not a statement about the full single-condition model after entering size-dependent segments.

## Reproduce

From the repository root:

```sh
julia --project=structural -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/export_structural_fixture.jl
julia --project=structural structural/test.jl
julia --project=structural structural/run_all.jl results/structural_new
```

For an additional segment or design check:

```sh
julia --project=structural structural/run.jl 5 pooled 1 local results/structural_extra
julia --project=structural structural/run.jl 5 single40 0 global results/structural_extra
```

CLI model IDs and saved CSV IDs are storage IDs: storage 5 is publication Model 6, and storage 6 is publication Model 5. The table above uses publication numbering. Segments 0 and 5 are endpoint clamps; 1–4 interpolate 200–400, 400–600, 600–800 and 800–1000 µm. Existing model-result directories are rejected. The finite-field calculation uses internal APIs from the pinned package version; run its tests before changing that environment. Symbolic work uses separate processes rather than the calibration's thread pool. The nominal initial-state values match the current configuration; revise the specialization if the initial-state policy changes.

Saved results are under `reference/structural/`. Local working outputs are under `results/structural/`.

Methods and definitions: [StructuralIdentifiability.jl documentation](https://docs.sciml.ai/StructuralIdentifiability/stable/tutorials/identifiability/), including the cited local-observability and differential-elimination algorithms. The package's generic-known-initial-condition option does not specialize initial conditions to particular numerical values, which is why the additional rank calculation is supplied.
