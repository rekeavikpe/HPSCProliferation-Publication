#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS=1
export HPSC_PROFILE_FULL_RANGE=1
export HPSC_PROFILE_POINT_PROGRESS=1
fit="${HPSC_FIXED2_REFERENCE:-results/training60_final/reference_fit}"
seeds="${HPSC_FIXED2_SEEDS:-results/training60_final/separate_ablations/profiles/all-active/full}"
out="${HPSC_FIXED2_OUTPUT:-results/training60_final/profile_fixed_two}"
resolution="${HPSC_PROFILE_RESOLUTION:-200}"
python3 - "$out" "$fit" "$seeds" <<'PY'
from pathlib import Path
import sys
out=Path(sys.argv[1]).resolve()
for value in sys.argv[2:]:
    source=Path(value).resolve()
    if out==source or out in source.parents or source in out.parents:
        raise SystemExit(f"Output overlaps protected input: {source}")
PY
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
    exec 9>"$out.launch.lock"
    flock -n 9 || { echo "Nine-parameter profiles are already running for this output" >&2; exit 1; }
fi
read -r -a julia_command <<< "${HPSC_JULIA_COMMAND:-julia +1.12.4}"
"${julia_command[@]}" --project=profiles tests/test_training60_fixed_two.jl "$fit"
echo "$(date -Iseconds) START nine-parameter Model 6 profiles"
"${julia_command[@]}" --project=profiles profiles/run_fixed_two.jl \
    "$fit" "$out" "$seeds" "$resolution" all
python3 scripts/audit_fixed_two_profiles.py "$out"
"${julia_command[@]}" --project=profile_figures profiles/render_fixed_two.jl "$out" "$out/figures"
echo "$(date -Iseconds) FINISHED nine profiles, audit and figures. Review audit.csv before interpretation."
