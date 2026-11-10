# Targeted convergence and lactate checks

These commands preserve the original calibration, ranking and selection lock. Use a new output directory for each command. They require the full five-condition 100-mL calibration, not a split fit. Their recorded observation identities, means and SDs must match the source fit.

```sh
julia +1.12.4 --threads=42 --project=. scripts/refine_selected.jl results/calibration100 results/refinement100
julia +1.12.4 --threads=42 --project=. scripts/lactate_limit.jl results/calibration100 results/lactate_limit100
julia +1.12.4 --threads=42 --project=. scripts/profile_multistart.jl results/calibration100 results/profile100_K_Lac 5 K_Lac 21
```

Run these sequentially on an allocated server. Threading occurs within each command; running three 42-thread processes together would oversubscribe the allocation. The first two commands have only eight and four new optimization jobs, respectively, so cannot use all 42 threads at once. The profile has many more independent jobs.

The refinement command analyzes Models 1 and 2 using the saved optimum plus three seeded perturbations. Original fits remain candidates so a failed refinement cannot discard a better original score. It uses up to 600 projected least-squares iterations and 12,000 Nelder–Mead iterations. It retains the original convergence tolerances; it does not relax them to manufacture a passing flag. The lowest-score converged stage within 1e-6 NLL of the achieved minimum is preferred; both scores are saved. Invalid starts are recorded with FAILED files and omitted from valid-candidate selection; they must be reported, not presented as successful runs.

Additional numerical diagnostics are available with `scripts/diagnose_precision.jl FIT_DIRECTORY OUTPUT_DIRECTORY` (ODE tolerances and finite-difference step sizes) and `scripts/refine_simplex.jl FIT_DIRECTORY OUTPUT_DIRECTORY` (two narrower simplex restarts at each saved competitor fit). The latter retains the 1e-8 stopping tolerance and does not prove global optimality even if its local stopping test passes.

The lactate-limit command sets the inhibition factor exactly to one while retaining lactate production and medium-exchange dynamics. It reoptimizes the ten remaining Model 5 parameters. K_Lac is stored as a dummy fixed value of 1 solely to preserve the parameter-file interface; it does not enter the equations in this branch. Its BIC counts ten estimated kinetic parameters. This is an exploratory nested-limit sensitivity, not a prespecified additional winner in the original seven-model comparison. A smaller BIC does not establish that lactate cannot inhibit growth biologically.

The multistart profile fixes one parameter at each grid point and reoptimizes all other parameters. It first uses three starts per point, borrowing the earlier single-start profile if present, followed by starts from the neighboring grid points in both directions. Endpoint points have one distinct neighboring start; interior points have two. Saved columns include stopping flags and the number of starts attaining the best achieved NLL within 1e-5. Such agreement is useful but does not prove a global nuisance optimum. Grid resolution, active bounds and the working likelihood limit interval interpretation. No confidence intervals are generated automatically. Negative profile likelihood differences indicate an improved minimum and must be investigated.

An interrupted profile can be resumed by adding `--resume` after the grid-point argument. The script validates the configuration, grid and recorded scientific source hashes, preserves completed starts, and records the resumed helper/script versions. Earlier FAILED records are retained under distinct filenames when retried. Run only one process per output directory.

Once the profile and limit commands finish, render their comparison with:

```sh
julia +1.12.4 --project=figures figures/render_targeted.jl results/profile100_K_Lac results/lactate_limit100 figures/targeted_output
```

For another Model 5 parameter, replace K_Lac and give a new directory. The other parameter names are K_Agg, K_d, K_G, K_oG, Y_GLac, Y_GX, alpha, mu_max, K_tau and K_eps. Their profiles must be completed and checked before claiming all Model 5 parameters are practically identifiable or selecting fixed values.

Local sensitivity results use the local reproduction calibration. The server's supplied 100-run ranking agrees closely, but its parameter files and diagnostics have not been transferred; execute these commands on the server to assess its exact saved solutions.
