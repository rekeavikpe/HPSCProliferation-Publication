#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
fit=${1:-results/calibration100_hybrid_lac100_kg100}
out=${2:-results/ablation_eight_shared_hybrid_lac100_kg100}
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
  exec 9>"${out}.lock"
  flock -n 9 || { echo "Eight-version calibration is already running." >&2; exit 1; }
fi
julia +1.12.4 --project=population -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=figures -e 'using Pkg; Pkg.instantiate()'
if [[ ! -f "$fit/figure_data/optimizer_bic.csv" ]]; then
  julia +1.12.4 --project=. scripts/export_figure_data.jl "$fit" "$fit/figure_data"
fi
julia +1.12.4 --threads="${JULIA_NUM_THREADS:-42}" --project=population population/run_routed_ablations.jl "$fit" "$out"
julia +1.12.4 --project=figures figures/render_routed_ablations.jl "$out"
echo "Calibration and comparison complete: $out/index.html"
