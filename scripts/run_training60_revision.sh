#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
unset HPSC_ALLOW_SUMMARY_ONLY
export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-42}
split=${1:-results/split6040_refined_hybrid_lac100_kg100}
exposure=${2:-results/david_revision_20261014/exposure_selection}
out=${3:-results/training60_final}
fit="$out/reference_fit"
mkdir -p "$out" logs
if command -v flock >/dev/null 2>&1; then
  exec 9>"$out/workflow.lock"
  flock -n 9 || { echo "This workflow is already running." >&2; exit 1; }
fi
# Prevent reuse of completed stages after changing the source fits or analysis code.
J=(julia +1.12.4)
for project in . population figures profiles profile_figures; do
  "${J[@]}" --project="$project" -e 'using Pkg; Pkg.instantiate()'
done
python3 revision/lock_training60.py "$split" "$exposure" "$out"
run_stage() {
  local name=$1
  shift
  if [[ -f "$out/$name.DONE" ]]; then echo "Already complete: $name"; return; fi
  printf '%s START %s\n' "$(date -Is)" "$name" | tee -a "$out/progress.log"
  "$@"
  touch "$out/$name.DONE"
  printf '%s FINISHED %s\n' "$(date -Is)" "$name" | tee -a "$out/progress.log"
}
# Always verify saved parameters, split membership and likelihood before resuming.
"${J[@]}" --project=population revision/prepare_training60.jl "$split" "$fit"
"${J[@]}" --project=. revision/test_training60.jl "$fit"
# Interrupted prediction-only folders are retained before deterministic regeneration.
predict_scale() {
  [[ ! -d "$fit/scale_transfer" ]] || mv "$fit/scale_transfer" "$fit/scale_transfer.interrupted.$(date +%s)"
  "${J[@]}" --project=. scripts/predict_scale.jl "$fit"
}
predict_media() {
  [[ ! -d "$out/media" ]] || mv "$out/media" "$out/media.interrupted.$(date +%s)"
  "${J[@]}" --project=. scripts/predict_media_exchange.jl "$fit" "$out/media"
}
predict_cells() {
  [[ ! -d "$out/cell_transfer" ]] || mv "$out/cell_transfer" "$out/cell_transfer.interrupted.$(date +%s)"
  "${J[@]}" --project=. scripts/predict_cell_lines.jl "$fit" "$out/cell_transfer"
}
run_stage scale predict_scale
run_stage figure_data "${J[@]}" --project=. scripts/export_figure_data.jl "$fit" "$fit/figure_data"
run_stage media predict_media
run_stage cell_transfer predict_cells
run_stage metrics python3 scripts/report_prediction_metrics.py "$fit" --media "$out/media" --cell-lines "$out/cell_transfer" --out "$out/prediction_metrics"
run_stage exposure "${J[@]}" --project=. revision/evaluate_training60_exposure.jl "$fit" "$exposure" "$out/exposure_transfer"
run_stage main_figures "${J[@]}" --project=figures figures/render_training60.jl "$out" main
run_stage media_figures "${J[@]}" --project=figures figures/render_media_exchange.jl "$out/media"
run_stage transfer_figures "${J[@]}" --project=figures figures/render_cell_lines.jl "$out/cell_transfer"
# These fits change because their inherited reference parameters now come from 60%.
run_stage cell_refits "${J[@]}" --threads="$JULIA_NUM_THREADS" --project=population revision/run_david_analyses.jl cells "$fit" "$split" "$out/cell_lines"
run_stage cell_figures "${J[@]}" --project=figures figures/render_training60.jl "$out" cells
# Preserve David's separate-response ablation design. No new model-selection search.
export HPSC_ABLATION_SHARED=false
run_stage ablations "${J[@]}" --threads="$JULIA_NUM_THREADS" --project=population population/run_routed_ablations.jl "$fit" "$out/separate_ablations"
run_stage ablation_selection "${J[@]}" --project=. revision/evaluate_training60_sensitivity.jl "$fit" "$out/separate_ablations" "$split"
run_stage ablation_figures "${J[@]}" --project=figures figures/render_routed_ablations.jl "$out/separate_ablations"
run_stage profiles bash scripts/run_routed_profiles.sh "$out/separate_ablations"
printf 'Completed training-only workflow. Original results preserved.\n' > "$out/COMPLETE"
echo "Results: $out"
