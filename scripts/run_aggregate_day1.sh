#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
# Conservative default while the independent mixed-exposure fit is running.
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS=1
export HPSC_AGGDAY1_STARTS="${HPSC_AGGDAY1_STARTS:-100}"
export HPSC_AGGDAY1_EVALUATIONS="${HPSC_AGGDAY1_EVALUATIONS:-20000}"
export HPSC_AGGDAY1_ITERATIONS="${HPSC_AGGDAY1_ITERATIONS:-3000}"
fit="${HPSC_AGGDAY1_REFERENCE:-results/training60_final/reference_fit}"
split="${HPSC_AGGDAY1_SPLIT:-results/split6040_refined_hybrid_lac100_kg100}"
out="${HPSC_AGGDAY1_OUTPUT:-results/aggregate_day1_fixed15_p50shear_p50edr}"
mkdir -p "$out"
if command -v flock >/dev/null 2>&1; then
    exec 9>"$out.launch.lock"
    flock -n 9 || { echo "This aggregate-day-1 output already has a running job" >&2; exit 1; }
fi
read -r -a julia_command <<< "${HPSC_JULIA_COMMAND:-julia +1.12.4}"
# Reuse the already installed environment without modifying the shared manifest.
echo "Checking the fixed-15-µm/day-1 design and preserved data split."
"${julia_command[@]}" --project=population tests/test_aggregate_day1.jl "$fit" "$split"
echo "Starting or resuming $HPSC_AGGDAY1_STARTS aggregate-day-1 starts with $JULIA_NUM_THREADS threads. Output: $out"
"${julia_command[@]}" --project=population revision/run_aggregate_day1.jl "$fit" "$split" "$out"
if [[ "${HPSC_AGGDAY1_RENDER:-0}" == "1" ]]; then
    "${julia_command[@]}" --project=figures figures/render_aggregate_day1.jl "$out"
else
    echo "Numerical comparison complete. Saved curves are ready for local plotting with Arial."
fi
