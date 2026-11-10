#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
fit=${1:-results/exposure_sensitivity_hybrid_lac100_kg100}
out=${2:-results/exposure_transfer_hybrid_lac100_kg100}
julia +1.12.4 --project=. scripts/evaluate_exposure_transfer.jl "$fit" "$out"
julia +1.12.4 --project=figures figures/render_exposure_transfer.jl "$out"
