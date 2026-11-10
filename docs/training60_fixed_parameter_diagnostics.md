# K_Lac and Y_GX diagnostic profiles

This separate run supplies the two panels remaining in the supplementary
diagnostic figure after K_G moves into the nine-parameter figure. It does not
change the primary calibration, predictions, MRPE, BIC or the completed nine
conditional profiles.

Both diagnostic parameters are allowed to vary. At each fixed grid value of
K_Lac, the other ten parameters are optimized, including Y_GX and K_G. At each
fixed grid value of Y_GX, the other ten are optimized, including K_Lac and K_G.
The primary 19-culture training data, likelihood, exposure inputs and bounds are
unchanged. No selection or external validation observations enter the profiles.
The existing converged all-active profile reference is checked against the
primary training likelihood and reused. No reference fit is rerun and no
parameter file used for predictions is replaced.

The saved all-active ProfileLikelihood grids provide coordinates and starting
guesses, not accepted likelihood results. Every distinct grid point is
reoptimized using isolated multistart BOBYQA and the same near-interval BFGS
check used in the nine-parameter recovery. Starts include the saved neighboring
solutions, the reference and an interpolated neighboring solution. Coordinates
near a bound also have an interior starting copy, which improves exploration
without changing the parameter bounds. Grid coordinates
within 1e-9 in normalized log space of the reference or a bound are merged; the
old K_Lac file had many effectively identical points at its upper bound.

Independent grid points run across up to JULIA_NUM_THREADS workers. This differs
from the sequential stepping within each profile in the nine-parameter runner.
Each point writes its parameters and diagnostics before writing COMPLETE. A
restart with identical sources and settings skips those saved points. Work in
an interrupted optimizer call cannot be retained. The grid extends at the
original edge spacing if a corrected curve has not yet crossed the confidence
threshold or reached the tested bound.

The final audit flags sharp dips, disconnected confidence regions, optimizer
failures or a likelihood better than the reference. A confidence limit that is
not reached within the tested bounds is reported as unresolved. The procedure
does not force a crossing or assume either parameter must be non-identifiable.
Original result files are never replaced.

## Server

```bash
cd ~/HPSCProliferation-R2
if git pull --ff-only origin main &&
   test -f scripts/run_training60_fixed_parameter_diagnostics.sh; then
  mkdir -p logs
  nohup env JULIA_NUM_THREADS=42 \
    bash scripts/run_training60_fixed_parameter_diagnostics.sh \
    >> logs/training60_fixed_parameter_diagnostics.log 2>&1 < /dev/null &
  echo "Started launcher PID: $!"
fi
```

Output: `results/training60_final/profile_fixed_parameter_diagnostics`.

```bash
tail -n 25 logs/training60_fixed_parameter_diagnostics.log
pgrep -af '[j]ulia.*run_fixed_parameter_diagnostics'
find results/training60_final/profile_fixed_parameter_diagnostics \
  -type f -path '*/points/*/COMPLETE' 2>/dev/null | wc -l
```

The count is the number of saved profile points across both curves. The log
states the number of points in each batch; extra batches are added if a curve
needs a wider range to reach its threshold or the existing parameter bound.

On completion, send the numerical output and log for inspection and local Arial
figure rendering. A server font substitution is avoided by not rendering there.

```bash
tar -czf training60_fixed_parameter_diagnostics_for_review.tar.gz \
  results/training60_final/profile_fixed_parameter_diagnostics \
  logs/training60_fixed_parameter_diagnostics.log
```
