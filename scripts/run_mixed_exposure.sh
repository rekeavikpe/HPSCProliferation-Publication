#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-42}"
export OPENBLAS_NUM_THREADS=1
export HPSC_MIXED_STARTS="${HPSC_MIXED_STARTS:-100}"
export HPSC_MIXED_EVALUATIONS="${HPSC_MIXED_EVALUATIONS:-20000}"
export HPSC_MIXED_ITERATIONS="${HPSC_MIXED_ITERATIONS:-3000}"
fit="${HPSC_MIXED_REFERENCE:-results/training60_final/reference_fit}"
split="${HPSC_MIXED_SPLIT:-results/split6040_refined_hybrid_lac100_kg100}"
out="${HPSC_MIXED_OUTPUT:-results/mixed_exposure_p50shear_p75edr}"
mkdir -p "$out"
if command -v flock >/dev/null 2>&1; then
    exec 9>"$out.launch.lock"
    flock -n 9 || { echo "This output already has a running job" >&2; exit 1; }
fi
read -r -a julia_command <<< "${HPSC_JULIA_COMMAND:-julia +1.12.4}"
"${julia_command[@]}" --project=population -e 'using Pkg; Pkg.instantiate()'
"${julia_command[@]}" --project=population tests/test_mixed_exposure.jl "$fit" "$split"
"${julia_command[@]}" --project=population revision/run_mixed_exposure.jl "$fit" "$split" "$out"
