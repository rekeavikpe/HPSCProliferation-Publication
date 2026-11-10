#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS=1
fit="${HPSC_AGGBOUND_REFERENCE:-results/training60_final/reference_fit}"
split="${HPSC_AGGBOUND_SPLIT:-results/split6040_refined_hybrid_lac100_kg100}"
source_fit="${HPSC_AGGBOUND_SOURCE:-results/aggregate_day1_fixed15_p50shear_p50edr}"
out="${HPSC_AGGBOUND_OUTPUT:-results/aggregate_day1_boundary_checks}"
# Reject overlapping or symlinked input/output paths before creating a directory or lock.
python3 - "$out" "$fit" "$split" "$source_fit" results/mixed_exposure_p50shear_p75edr <<'PY'
from pathlib import Path
import sys
out = Path(sys.argv[1]).resolve()
for arg in sys.argv[2:]:
    source = Path(arg).resolve()
    if out == source or out in source.parents or source in out.parents:
        raise SystemExit(f"Output overlaps a protected input: {source}")
PY
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
    exec 9>"$out.launch.lock"
    flock -n 9 || { echo "Boundary checks are already running for this output" >&2; exit 1; }
fi
read -r -a julia_command <<< "${HPSC_JULIA_COMMAND:-julia +1.12.4}"
"${julia_command[@]}" --project=population tests/test_aggregate_day1_bounds.jl "$fit" "$split" "$source_fit"
"${julia_command[@]}" --project=population revision/run_aggregate_day1_bounds.jl "$fit" "$split" "$source_fit" "$out"
