# Model 6 profiles with two fixed parameters

This analysis tests the uncertainty of nine Model 6 parameters using the same 19 training cultures and P50 shear/P50 EDR inputs as the current manuscript. Only `K_Lac` and `Y_GX` are fixed, at their saved training estimates (approximately 100 mM and 3.961389e-6 µmol/cell). `K_G` is profiled and is re-optimized when profiling every other parameter. The 12 selection cultures and external validation datasets are not used to estimate these profiles.

The nine profiled parameters are `mu_max`, `K_d`, `K_G`, `K_oG`, `Y_GLac`, `alpha`, `K_Agg`, `K_tau` and `K_eps`. The fixed values were estimated from the training data, so these profiles still exclude uncertainty in `K_Lac` and `Y_GX`.

The reference is `results/training60_final/reference_fit`. The script preserves its data allocation, original day-zero initialization, likelihood weights and parameter bounds. Model 6 has storage ID 5 in this repository. This is separate from the aggregate-day-one and mixed-exposure analyses.

The reference optimum is checked/refined under the two fixed values before profiling. Any likelihood improvement is reported in `mle_summary.csv`. It does not replace the saved manuscript fit or its predictions. At each profile point, the existing multistart optimizer uses the previous solution, the new profile reference and three nearby nuisance solutions from the saved all-active Model 6 profiles. Those saved solutions are starting guesses only. Every likelihood evaluation uses the new two-fixed-parameter configuration. A BFGS refinement checks the best BOBYQA solution near the confidence region.

Both directions are explored over the existing bounds with 200 grid values per direction. Boundary-limited intervals, optimizer concerns, disconnected regions and a better likelihood than the reference are flagged in `audit.csv`. Completion alone does not establish identifiability. Review the curves and audit before updating any documents.

## Run on the server

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
mkdir -p logs
nohup env JULIA_NUM_THREADS=8 \
  bash scripts/run_training60_fixed_two_profiles.sh \
  >> logs/training60_fixed_two_profiles.log 2>&1 < /dev/null &
echo $!
```

Use 8 threads alongside an existing run. The script has its own output lock and does not stop other jobs. If more cores are available, the thread count may be increased. Keep the same inputs, output and resolution when resuming.

Output: `results/training60_final/profile_fixed_two`. Original eight-parameter profiles, all-active profiles and calibration results are read-only. Each completed parameter has a `COMPLETE` marker and is skipped when the exact same command is resumed. Changed source files, settings or saved seeds cause a resume error instead of silently mixing results.

```bash
tail -n 30 logs/training60_fixed_two_profiles.log
find results/training60_final/profile_fixed_two \
  -mindepth 2 -maxdepth 2 -type f -name COMPLETE | wc -l
pgrep -af '[j]ulia.*run_fixed_two'
```

The count is out of nine parameters, not 100 starts. Point-level progress appears after every 20 nuisance optimizations. The script exports a 3×3 figure using the existing ProfileLikelihood.jl plotting style and requests Arial. If Arial is not installed on the server, Makie reports a fallback font. The saved results can be rendered again locally with Arial without repeating any optimization. A finished suite has `audit.csv` and a final `FINISHED` line in the log.

## Share completed results

```bash
cd ~/HPSCProliferation-R2
tar -czf training60_fixed_two_profiles_for_review.tar.gz \
  --exclude='*.tiff' --exclude='*.eps' --exclude='*.jpeg' --exclude='*.svg' --exclude='*.bin' \
  results/training60_final/profile_fixed_two \
  logs/training60_fixed_two_profiles.log
```

The archive retains numeric profiles, nuisance solutions, diagnostics, PDF and PNG figures. It does not modify the manuscript submission folder or ZIP.

## Recover the October 22 interval-pass failure

The failure `KeyError: key 20 not found` occurred in the starting-guess selector.
The previous runner captured a shared `old` profile table and a shared `fixed`
coordinate in closures used by parallel jobs. Each parameter now owns its seed
bank and each call owns its fixed coordinate. The likelihood, training cultures,
fixed parameters, bounds and optimization budgets are unchanged.

After the interval process has exited, run:

```bash
cd ~/HPSCProliferation-R2
git pull --ff-only origin main
nohup env JULIA_NUM_THREADS=9 \
  bash scripts/run_training60_fixed_two_recovery.sh \
  >> logs/training60_fixed_two_recovery.log 2>&1 < /dev/null &
echo $!
```

This creates `results/training60_final/profile_fixed_two_intervals_recovered`.
It never edits the failed interval folder, the full-range folder or the original
training fit, and it does not stop the separate full-range process. It rejects
changed data, seeds, likelihood code or settings. Only the specific known
defective runner checksum is allowed to migrate to the corrected runner.

The recovery copies every completed profile and the numerical reference MLE.
It reconstructs the likelihood under the current code, verifies every saved
profile likelihood and reoptimizes the points bracketing each 95% limit plus
any unconfirmed optimizer points near the interval. If the improvement exceeds
1e-5 log-likelihood units or a recheck does not converge, only that profile is
recomputed. Its previous files remain in `before_seed_recheck`. Missing profiles
run with isolated starting guesses. `recovery.json` records the hashes of the
original copied files. The old binary closure is not imported.

The already-running full-range job still has the old code in memory. Its results
are not independent verification of the corrected seed selection and must not
be used to replace recovered results without checking them. The recovery pass
remains a threshold scan rather than a search for distant disconnected regions.

```bash
tail -n 30 logs/training60_fixed_two_recovery.log
python3 - <<'PY'
from pathlib import Path
root = Path('results/training60_final/profile_fixed_two_intervals_recovered')
names = ['K_Agg', 'K_d', 'K_G', 'K_oG', 'Y_GLac', 'alpha', 'mu_max', 'K_tau', 'K_eps']
done = [n for n in names if (root/n/'COMPLETE').is_file()
        and not (root/n/'RECOVERED_NEEDS_CHECK').exists()
        and not (root/n/'RECOMPUTE_REQUIRED').exists()]
print(f'Checked or newly completed profiles: {len(done)}/9')
print('Done:', ', '.join(done) or 'none yet')
print('Remaining:', ', '.join(n for n in names if n not in done))
PY
```

Resume the same recovery command if interrupted. Completed checks and newly
completed profiles are skipped, subject to the same source and setting guards.
Review `audit.csv`, `recovery_recheck.csv` files and curves before inference.
The recovery does not change any manuscript or submission document.
