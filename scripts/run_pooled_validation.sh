#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
out=${1:-results/pooled_validation3_lac100_kg100}
export OPENBLAS_NUM_THREADS=1
if [[ -n "${HPSC_JULIA:-}" ]]; then julia_cmd=("$HPSC_JULIA"); else
  export PATH="$HOME/.juliaup/bin:$PATH"
  julia_cmd=(julia +1.12.4)
fi
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
  exec 9>"${out}.lock"
  flock -n 9 || { echo "Already running: $out" >&2; exit 1; }
fi
"${julia_cmd[@]}" --project=population -e 'using Pkg; Pkg.instantiate()'
"${julia_cmd[@]}" --project=population population/test_pooled_validation.jl
"${julia_cmd[@]}" --threads="${JULIA_NUM_THREADS:-42}" --project=population population/run_pooled_validation.jl config/calibration_lac100_kg100.toml "$out"
