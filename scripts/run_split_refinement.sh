#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
source=${1:-results/split6040_hybrid_lac100_kg100}
out=${2:-results/split6040_refined_hybrid_lac100_kg100}
threads=${JULIA_NUM_THREADS:-42}
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
  exec 9>"${out}.lock"
  flock -n 9 || { echo "This refinement output is already running." >&2; exit 1; }
fi
julia +1.12.4 --threads="$threads" --project=population population/refine_split_selection.jl "$source" "$out"
julia +1.12.4 --project=figures figures/render_split_selection.jl "$out"
echo "Inspect $out/refinement_comparison.csv and $out/selection_ranking.csv"
