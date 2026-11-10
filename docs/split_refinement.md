# Targeted refinement of the 60/40 comparison

Refine publication Models 6, 1, 4 and 5 (storage IDs 5,1,4,6) using calibration cultures only. This addresses the unmet optimizer stopping flags in their selected production fits. Six local runs per model give 24 jobs:

- Restart from the saved selected training fit.
- Restart from the best candidates in two other original runs, chosen by training NLL.
- Use three seeded small perturbations of the selected fit in normalized logarithmic coordinates.

If a test fixture has fewer original runs, use seeded perturbations for the missing alternatives. The production dataset has 100 original runs per model. The local iteration budget is increased to 12,000; the original Nelder–Mead method, parameter bounds and convergence tolerance are unchanged. The workflow records actual optimizer status; it does not guarantee convergence.

The original result directory is read only. Copy its candidates into `results/split6040_refined_hybrid_lac100_kg100`, add run101–run106 for each target model, then reselect training fits using the existing training-only rule. Only afterward recompute the held-out scores using the unchanged selection criterion. The original 100 starts remain the boxplot distribution; targeted refinements are not counted as independent population searches.

After pushing and pulling the updated files, run on the server:

```bash
cd ~/HPSCProliferation-R2
mkdir -p logs
nohup env JULIA_NUM_THREADS=42 bash scripts/run_split_refinement.sh \
  > logs/split_refinement.log 2>&1 < /dev/null &
```

Watch with `tail -f logs/split_refinement.log`. Ctrl+C exits the watcher only. Repeat the same launch command after an interruption to reuse verified completed candidates. Do not change code, data, settings or source outputs while resuming; fingerprints detect such changes.

Main outputs:

- `refinement_comparison.csv`: old/new training NLL, BIC, stopping flags and held-out NLL.
- `selection_ranking.csv`: updated predictive ranking.
- `state_errors.csv`: updated calibration/selection MRPE by state.
- `figures/split_model_comparison.pdf`: labeled BIC and predictive-score panels.
- `model*/run10*/local/diagnostics.toml`: iteration counts and actual convergence status for refinements.

If selection still matches the existing all-100-mL prediction model, regenerate the combined figure using:

```bash
bash scripts/render_split_combined.sh results/split6040_refined_hybrid_lac100_kg100
```

The exporter verifies the model and settings before reusing the all-100-mL predictions. This refinement does not itself recalibrate on all 100-mL cultures or rerun profiles/external predictions.
