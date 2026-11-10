#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
split=${1:-results/split6040_hybrid_lac100_kg100}
full=${2:-results/calibration100_hybrid_lac100_kg100/figure_data}
julia +1.12.4 --project=figures figures/export_split_combined.jl \
  "$split" "$full" "$split/combined_inputs"
julia +1.12.4 --project=figures figures/render_bic_scale.jl \
  "$split/combined_inputs" "$split/combined_figure"
echo "Combined figure saved to $split/combined_figure/bic_and_500mL.pdf"
