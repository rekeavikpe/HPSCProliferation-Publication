#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
threads=${JULIA_NUM_THREADS:-42}
out=${1:-results/split6040_hybrid_lac100_kg100}
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
  exec 9>"${out}.lock"
  flock -n 9 || { echo "This result directory is already running." >&2; exit 1; }
fi
export OPENBLAS_NUM_THREADS=1
julia +1.12.4 --project=population -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --threads="$threads" --project=population population/run_split_selection.jl \
  config/calibration_lac100_kg100.toml "$out"
julia +1.12.4 --project=figures -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=figures figures/render_split_selection.jl "$out"
echo "Completed split comparison: $out/selection_ranking.csv"
echo "Next: inspect the selected model, then refit it on all 100-mL cultures before final profiles and external evaluations."
