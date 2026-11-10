#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
fit=${1:-results/calibration100_hybrid_lac100_kg100}
media=${2:-results/media_exchange_hybrid_lac100_kg100}
cells=${3:-results/cell_lines_frozen_hybrid_lac100_kg100}
out=${4:-results/prediction_metrics_hybrid_lac100_kg100}
if [[ ! -e "$media" ]]; then
  julia +1.12.4 --project=. scripts/predict_media_exchange.jl "$fit" "$media"
fi
[[ -f "$media/COMPLETE" ]] || { echo "Incomplete media evaluation: $media"; exit 1; }
if [[ ! -e "$cells" ]]; then
  julia +1.12.4 --project=. scripts/predict_cell_lines.jl "$fit" "$cells"
fi
[[ -f "$cells/COMPLETE" ]] || { echo "Incomplete cell-line evaluation: $cells"; exit 1; }
python3 scripts/report_prediction_metrics.py "$fit" --media "$media" --cell-lines "$cells" --out "$out"
julia +1.12.4 --project=figures figures/render_media_exchange.jl "$media"
julia +1.12.4 --project=figures figures/render_cell_lines.jl "$cells"
