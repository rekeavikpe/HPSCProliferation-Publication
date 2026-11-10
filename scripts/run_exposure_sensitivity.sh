#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
fit=${1:-results/calibration100_hybrid_lac100_kg100}
out=${2:-results/exposure_sensitivity_hybrid_lac100_kg100}
threads=${JULIA_NUM_THREADS:-42}
resume=${3:-}
extra=()
if [[ -n "$resume" ]]; then
  [[ "$resume" == "--resume" ]] || { echo "Third argument must be --resume" >&2; exit 2; }
  extra=(--resume)
fi
# Prevent two copies of this wrapper from writing to the same result directory.
if command -v flock >/dev/null 2>&1; then
  mkdir -p "$(dirname "$out")"
  exec 9>"${out}.lock"
  flock -n 9 || { echo "This output is already locked by another run" >&2; exit 1; }
fi
julia +1.12.4 --threads="$threads" --project=population population/exposure_sensitivity.jl "$fit" "$out" 100 20000 "${extra[@]}"
julia +1.12.4 --project=figures figures/render_exposure_sensitivity.jl "$out"
