#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
fit=${1:-results/calibration100_hybrid_lac100_kg100}
out=${2:-results/media_exchange_hybrid_lac100_kg100}
julia +1.12.4 --project=. scripts/predict_media_exchange.jl "$fit" "$out"
julia +1.12.4 --project=figures figures/render_media_exchange.jl "$out"
