# hPSC proliferation: R2 analysis repository

Model numbering for the revision: the growth-coupled EDR model is **Model 6** (original Model 6); the dual-interaction model is **Model 5** (original Model 7). Saved analysis folders and CLI model arguments retain their storage IDs: `model5` stores publication Model 6, and `model6` stores publication Model 5. See [the numbering crosswalk](docs/model_numbering.md) before reusing saved outputs.


Mechanistic models of viable cells, glucose, lactate and aggregate diameter in vertical-wheel bioreactors. Julia 1.12.4; public model IDs 0–6. Model 0 omits hydrodynamic terms; the algebraically duplicate historical Model 2 is excluded.

**Status: verified analysis software under development, not a completed publication release.** The included reference results reproduce the earlier 100-mL calibration / retrospective 500-mL evaluation. Structural analysis for Models 0–6 and manuscript-style figures are included, with their assumptions documented. Cross-validation, final profiles and evaluations of other cell lines/media protocols remain unfinished. See [verification](docs/verification.md) and [remaining work](docs/remaining_work.md).

For the requested 100-run analysis, start with [START_HERE.md](START_HERE.md). The supplied reference figures came from six initial runs plus two refinements, not 100. New boxplots include initial runs only; refinement scores are exported separately.

## Start in VS Code

Open this folder in VS Code, install the Julia extension, and use a terminal from the repository root. Install Julia 1.12.4 on each machine. The committed manifests pin package versions; do not run `Pkg.update()` when reproducing a release.

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. scripts/test.jl
julia --startup-file=no --project=. scripts/verify_reference.jl
```

The reference check recomputes likelihoods and 500-mL predictions from saved parameter vectors. It is deliberately separate from reoptimizing the models.

## Structural identifiability

```sh
julia --project=structural -e 'using Pkg; Pkg.instantiate()'
julia --project=structural structural/test.jl
julia --project=structural structural/run_all.jl results/structural_new
```

The equations preserve the actual CFD feedback and pool the five 100-mL conditions. Completed checks support generic local identifiability at the configured initial states for all seven models. Global classifications have qualifications for Models 4 and 6 and for the segment/initial-state assumptions. Read [the structural analysis](docs/structural_identifiability.md) before using these results in the manuscript. Structural identifiability does not establish practical identifiability.

## Production hybrid calibration

The current revision workflow uses the original hybrid structure: differential evolution followed directly by Nelder–Mead, with actual stopping flags retained. It fits the current fixed-SD likelihood on the five 100-mL conditions, with K_G and K_Lac capped at 100 mM.

```sh
JULIA_NUM_THREADS=42 bash scripts/run_population100.sh
```

This runs 100 searches per model, verifies the selected fits, and generates the calibration, BIC and 500-mL figures in `results/calibration100_hybrid_lac100_kg100/figures`. Population searches use 20,000 evaluations and local refinement has a 3,000-iteration cap. Existing result directories are preserved. See [the method and checks](docs/identifiability_workflow.md). Model 3's stopping flag remains a reported limitation; more local iterations did not resolve it in the targeted check.

## Earlier local-only diagnostic calibration

```sh
julia --startup-file=no --threads=42 --project=. scripts/calibrate.jl results/calibration
```

Use `--threads=auto` on your laptop. On a server, request at least as many allocated CPUs as the Julia thread count. The script distributes 42 independent starts across threads, then runs 14 refinements. It sets BLAS to one thread to avoid nested oversubscription. More threads than independent jobs do not accelerate that stage. To use a larger prespecified multistart budget, add a starts-per-model argument, e.g. `scripts/calibrate.jl results/calibration12 12` schedules 84 initial jobs. Changing the start budget defines a new analysis run; record it rather than silently comparing it with the six-start reference. Seeds depend on model/start, not worker identity. Files are written separately per start; no shared output file is mutated by workers.

Every completed start retains its initial vector, final vectors, objective scores, stopping flags and timing. Existing run directories are rejected to prevent accidental mixing of configurations. Automatic resume is not yet implemented: completed start folders are checkpoints for diagnosis, not a promise of automatic restart. Use a new output directory for another run.

Verify a completed full-data run with `julia --project=. scripts/verify_run.jl results/calibration`.

The final `ranking.csv` is based on achieved training BIC. Read convergence flags and active bounds before interpreting its ranking; neither a stopping flag nor multistart agreement proves a global optimum or identifiability.

## Random whole-culture 60/40 splits

```sh
julia --project=. scripts/make_splits.jl
julia --threads=42 --project=. scripts/cross_validate.jl 1 results/split1
```

The first command records ten predetermined random splits. Repeat the second with split numbers 2–10 and distinct directories to run the other splits. Each split is stratified by the five 100-mL operating conditions. Every culture's days and states stay together. Fractions are approximate because numbers of cultures are integers. Split-seed selection never uses observed model performance. Splits share observations, so their variation is not an independent-sample confidence interval.

Calibration means, SD weights and initial glucose are computed from the calibration cultures only. Selection-set day-zero glucose is an allowed measured input; if no such measurement exists, use and document the calibration-set initial glucose for that condition. Endpoints absent from a selection subset are not imputed or scored. No selection outcomes are used for fitting or initial guesses. Held-out errors are evaluated on condition/day/state means, using calibration-set scales; they are not held-out BIC. These scripts support assessing sensitivity to splitting; the final manuscript selection rule remains to be agreed. Unknown batches limit claims of independence.

## Profiles and scale transfer

The additional convergence refinements, multistart profiles and exact inactive-lactate sensitivity are documented in [targeted checks](docs/targeted_checks.md). They write separate output directories and preserve the original calibration. The single-start profile command below is retained for preliminary diagnosis; use the multistart command before interpreting weak-parameter profiles.

After a full calibration (not a split fit):

```sh
julia --threads=42 --project=. scripts/profile.jl results/calibration 5 K_Lac 21
julia --project=. scripts/predict_scale.jl results/calibration
```

The profile reoptimizes all other free parameters at each fixed grid value. The first implementation uses one MLE-based nuisance start per point; unresolved minima require additional starts before confidence intervals are reported. It profiles within the saved bounds. It does not turn a flat or truncated profile into a finite interval. Check negative likelihood differences, which indicate an improved minimum relative to the saved fit.

The prediction script requires the saved selection lock and verifies the selected-parameter hash. It uses the corresponding 500-mL CFD inputs and protocol without refitting. These are retrospective evaluations because the data and model family were previously examined.

## CairoMakie figures

Graphics have a separate pinned environment so graphics installation cannot change the fitting dependencies.

```sh
julia --project=figures -e 'using Pkg; Pkg.instantiate()'
julia --project=figures figures/render.jl
# To plot a new calibration after running predict_scale.jl:
julia --project=. scripts/export_figure_data.jl results/calibration results/calibration/figure_data
julia --project=figures figures/render.jl results/calibration/figure_data
```

PDF/SVG are drawn with CairoMakie. A pinned Ghostscript binary converts PDF to EPS, 600-dpi PNG/TIFF and high-quality JPEG. JPEG is a convenience copy and is **not accepted by the journal for figure submission**. Figure dimensions are in points at physical publication size, not an arbitrary screen size. Feed jumps are retained. Legends/captions and alt text are separate from the image. See [figure specifications](docs/figures.md).

The default renderer follows manuscript Figures 2, 3 and 5: condition colors/markers, boxed axes, light grids, Arial, a four-state-by-condition calibration layout, horizontal BIC boxplots, and a four-panel prediction layout. It exports both baseline and hydrodynamic calibration plots with shared axis limits, plus full-start and best-score BIC views. Five training columns replace the original six because 500 mL is evaluated separately. Numerical results are fresh; the old plotted values are not reused. Arial must be available on the rendering machine; `HPSC_FONT=/path/to/font.ttf` and `HPSC_BOLD_FONT=/path/to/bold.ttf` select explicitly supplied fonts. Font substitution changes the appearance and requires a new visual check. The older `render_diagnostic.jl` is retained only as a diagnostic alternative.

## Desktop → GitHub → lab server

This folder is separate from your existing `Desktop/HPSCProliferation` repository. Nothing has been pushed. The local source workbooks are excluded from Git; curated CSVs and saved reference values are included. Check data-sharing arrangements before making a remote public.

For a new repository, create an empty **private** GitHub repository named `HPSCProliferation-R2` (no generated README), then run:

```sh
cd ~/Desktop/HPSCProliferation-R2
# This working folder has already been initialized; git init is only needed for a new copy.
git init -b main
git status --short
git add .
git diff --cached --stat
git commit -m "Prepare reproducible R2 model analysis"
git remote add origin https://github.com/rekeavikpe/HPSCProliferation-R2.git
git push -u origin main
```

The account name above comes from the remote in your existing Desktop repository; the new repository must be created before pushing. To use a different name/account, change the URL. No force push is needed. Keep data/source and generated results out of a routine code push. Reference CSVs are deliberately versioned for verification.

In VS Code, use Remote–SSH to connect to the lab server, then clone and open the new repository there:

```sh
git clone https://github.com/rekeavikpe/HPSCProliferation-R2.git
cd HPSCProliferation-R2
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/test.jl
julia --project=. scripts/verify_reference.jl
```

If Slurm is used, edit `server/slurm.sbatch` for your lab's account, partition, memory allocation and Julia module before `sbatch server/slurm.sbatch`. If you have an allocated machine and run directly, use `HPSC_THREADS=42 bash server/run.sh`. Do not run intensive fitting on a shared login node. No server configuration or connection has been assumed or tested.

After changes on your Desktop: `git add`, `git commit`, `git push`; on the server: `git pull --ff-only`. Keep run directories named and separate. Archive the exact commit, configuration and output directory for each analysis. Do not expect bitwise agreement between different CPU architectures; use the numerical verification tolerances.

## Layout

- `src/`: model equations, solver, optimization and replicate splitting.
- `scripts/`: runnable analysis commands.
- `config/`: likelihood, bounds and solver settings.
- `data/processed/`: curated measurements with source-cell provenance.
- `data/conditions.csv`: proposed manuscript-to-legacy condition mapping; source IDs are unchanged.
- `reference/`: saved results for reproducibility checks.
- `figures/`: separate CairoMakie environment and export code.
- `structural/`: pinned symbolic environment, equations, checks and nominal-initial-state rank calculations.
- `test/`: analytic and numerical checks.
- `docs/`: protocol, verification, journal requirements and unresolved work.

The Gaussian score currently uses condition means and between-culture SDs with floors; temporal covariance is not modeled. BIC therefore uses a working likelihood. All kinetic parameters in this reproduction are free. Fixing parameters after estimating them from these same data does not automatically reduce the fitted parameter count.

## Optional workbook import audit

All fitting/plotting commands above use Julia and the included curated CSVs. The source-cell audit can also be rebuilt from the locally retained workbooks with Python 3 and `openpyxl==3.1.5`:

```sh
python3 -m venv .venv
.venv/bin/pip install -r tools/requirements.txt
.venv/bin/python tools/import_data.py --out results/reimport
```

The importer checks source workbook hashes before extraction. Seven regenerated experimental CSVs were verified byte-for-byte against the curated files. `tools/import_exposure.py` rebuilds the CFD table in `data/processed`; use it only when deliberately refreshing that tracked table. Workbooks are ignored by Git and must be transferred separately if the server needs this optional import audit.

### Complete identifiability analysis

See [the identifiability workflow](docs/identifiability_workflow.md) for all-seven-model structural checks and the restored **ProfileLikelihood.jl** workflow, following the original Fminbox/BFGS → BOBYQA approach. The earlier custom profile suite is retired.

The requested random 60/40 culture-level comparison is documented in [the split workflow](docs/split6040_workflow.md). Existing all-100-mL results remain separate.
