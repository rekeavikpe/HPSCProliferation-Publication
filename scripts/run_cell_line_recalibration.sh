#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
fit=${1:-results/calibration100_hybrid_lac100_kg100}
frozen=${2:-results/cell_lines_frozen_hybrid_lac100_kg100}
out=${3:-results/cell_lines_recalibrated_hybrid_lac100_kg100}
threads=${JULIA_NUM_THREADS:-42}
if [[ ! -e "$frozen" ]]; then
  julia +1.12.4 --project=. scripts/predict_cell_lines.jl "$fit" "$frozen"
fi
julia +1.12.4 --threads="$threads" --project=population population/recalibrate_cell_lines.jl "$fit" "$frozen" "$out" 20 20000
julia +1.12.4 --project=figures figures/render_cell_lines.jl "$frozen" "$out"
metric_args=("$fit" --cell-lines "$out" --out "$out/metrics")
media=${4:-results/media_exchange_hybrid_lac100_kg100}
if [[ -f "$media/COMPLETE" ]]; then
  metric_args+=(--media "$media")
fi
python3 scripts/report_prediction_metrics.py "${metric_args[@]}"
