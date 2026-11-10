#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
threads=${JULIA_NUM_THREADS:-42}
out=${1:-results/cross_validation_balanced9_lac100_kg100}
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
  exec 9>"${out}.lock"
  flock -n 9 || { echo "This result directory is already running." >&2; exit 1; }
fi
export OPENBLAS_NUM_THREADS=1
# Explicit executable override also supports installations without juliaup.
if [[ -n "${HPSC_JULIA:-}" ]]; then
  julia_cmd=("$HPSC_JULIA")
else
  export PATH="$HOME/.juliaup/bin:$PATH"
  julia_cmd=(julia +1.12.4)
fi
"${julia_cmd[@]}" --version
"${julia_cmd[@]}" --project=population -e 'using Pkg; Pkg.instantiate()'
"${julia_cmd[@]}" --threads="$threads" --project=population population/run_cross_validation.jl \
  config/calibration_lac100_kg100.toml "$out"
if [[ "${HPSC_CV_PREPARE_ONLY:-0}" == "1" ]]; then exit 0; fi
"${julia_cmd[@]}" --project=figures -e 'using Pkg; Pkg.instantiate()'
"${julia_cmd[@]}" --project=figures figures/render_cross_validation.jl "$out"
echo "Completed: $out/state_metrics.csv and $out/figures"
