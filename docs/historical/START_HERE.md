> **Current training-only workflow (16 October 2026):** use
> `scripts/run_training60_revision.sh` and [the training-only workflow guide](docs/training60_workflow.md).
> This retains the original 60/40 culture split and does not refit on all 100 mL data.
> Earlier all-100-mL workflows below remain available for reproducing historical results.

> **4 October update:** the production calibration now uses the original-style population-search → Nelder–Mead hybrid. After syncing the repository, run `JULIA_NUM_THREADS=42 bash scripts/run_population100.sh`. It saves results and figures under `results/calibration100_hybrid_lac100_kg100`. The local-only calibration instructions below describe the earlier diagnostic workflow. See `docs/identifiability_workflow.md` for settings and convergence qualifications.

Model numbering for the revision: the growth-coupled EDR model is **Model 6** (original Model 6); the dual-interaction model is **Model 5** (original Model 7). Saved analysis folders and CLI model arguments retain their storage IDs: `model5` stores publication Model 6, and `model6` stores publication Model 5. See [the numbering crosswalk](docs/model_numbering.md) before reusing saved outputs.


# Start the 100-run analysis

This is a runnable development snapshot, not the final publication release. Use Julia 1.12.4. The raw workbooks are not required: curated input CSVs and reference outputs are included. See README.md for GitHub and VS Code Remote–SSH instructions.

Open this directory in VS Code and run the following from its terminal:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. scripts/test.jl
julia --startup-file=no --project=. scripts/verify_reference.jl
```

On an allocated lab-server machine with at least 42 CPUs:

```sh
julia --startup-file=no --threads=42 --project=. scripts/calibrate.jl results/calibration100 100
```

This fits all seven models (0–6) to all five 100-mL conditions: **100 initial runs per model, 700 total**, followed by two refinements per model. The first initial vector is the bounds midpoint; the remaining 99 are seeded randomized initial vectors. These are optimizer starts using the same biological dataset, not 100 biological samples or bootstrap replicates. Each run saves both optimization stages and diagnostics. The BIC boxplots use the 100 initial runs; the 14 dependent refinements are exported separately. Best-achieved model ranking includes refinement results.

The full 100-run job has not yet been executed and verified here. The included completed reference calibration used six initial runs plus two refinements per model. The fitting kernel and prior 42-thread scheduling checks have passed; runtime and memory use on your server remain unmeasured. Use `--threads=auto` on a laptop. A Slurm template is provided in server/slurm.sbatch; set your lab's required resource fields before submitting.

After calibration finishes:

```sh
julia --project=. scripts/verify_run.jl results/calibration100
julia --project=. scripts/predict_scale.jl results/calibration100
julia --project=. scripts/export_figure_data.jl results/calibration100 results/calibration100/figure_data
julia --project=figures -e 'using Pkg; Pkg.instantiate()'
julia --project=figures figures/render.jl results/calibration100/figure_data
```

Figures are saved under figures/output. Rendering overwrites previous figure files, so copy that directory when archiving a run. Arial and Arial Bold are used; provide HPSC_FONT and HPSC_BOLD_FONT paths if needed on the server. PDF, EPS, SVG, 600-dpi PNG/TIFF and JPEG convenience copies are produced.

Run the structural analysis separately:

```sh
julia --project=structural -e 'using Pkg; Pkg.instantiate()'
julia --project=structural structural/test.jl
julia --project=structural structural/run_all.jl results/structural_new
```

The completed structural results and assumptions are in reference/structural and docs/structural_identifiability.md. Do not interpret structural identifiability as practical certainty about fitted parameter values.

The calibration command above uses all 100-mL data; it does not perform the optional random 60/40 analysis. Final profiles, split analyses and other-cell-line/media-protocol evaluations are still being developed. Keep this snapshot and its run outputs together. Existing calibration directories are rejected; choose a new directory for each run. Cross-platform numerical agreement is checked by tolerances, not promised to be bitwise identical.
