#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
threads=${JULIA_NUM_THREADS:-42}
reference=${1:-results/calibration100_lac100_kg100/figure_data/optimizer_bic.csv}
search=${2:-results/hybrid100_lac100_kg100_search}
fit=${3:-results/calibration100_hybrid_lac100_kg100}
test -f "$reference" || { echo "Missing reference BIC file: $reference" >&2; exit 1; }
test ! -e "$search" && test ! -e "$fit" || { echo "Choose new result directories." >&2; exit 1; }
julia +1.12.4 --project=. -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=population -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=figures -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --threads="$threads" --project=population population/run.jl \
  config/calibration_lac100_kg100.toml "$search" 100 20000 0,1,2,3,4,5,6 "$reference" --original-hybrid
julia +1.12.4 --project=population population/finalize.jl "$fit" "$search"
julia +1.12.4 --project=. scripts/verify_run.jl "$fit"
julia +1.12.4 --project=. scripts/predict_scale.jl "$fit"
julia +1.12.4 --project=. scripts/export_figure_data.jl "$fit" "$fit/figure_data"
julia +1.12.4 --project=figures figures/render.jl "$fit/figure_data" "$fit/figures"
echo "Calibration saved to $fit. Inspect ranking.csv and stopping diagnostics before interpretation."
