#!/usr/bin/env bash
# Run the existing nine-parameter workflow to the 95% threshold on each side.
# This uses a separate output and leaves the full-range scan running untouched.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-9}"
export OPENBLAS_NUM_THREADS=1
export HPSC_PROFILE_FULL_RANGE=0
export HPSC_PROFILE_POINT_PROGRESS=1
fit="${HPSC_FIXED2_REFERENCE:-results/training60_final/reference_fit}"
seeds="${HPSC_FIXED2_SEEDS:-results/training60_final/separate_ablations/profiles/all-active/full}"
out="${HPSC_FIXED2_INTERVAL_OUTPUT:-results/training60_final/profile_fixed_two_intervals}"
full="${HPSC_FIXED2_FULL_OUTPUT:-results/training60_final/profile_fixed_two}"
resolution="${HPSC_PROFILE_RESOLUTION:-200}"
python3 - "$out" "$fit" "$seeds" "$full" <<'PY'
from pathlib import Path
import sys
out = Path(sys.argv[1]).resolve()
for value in sys.argv[2:]:
    protected = Path(value).resolve()
    if out == protected or out in protected.parents or protected in out.parents:
        raise SystemExit(f"Output overlaps protected input or full scan: {protected}")
PY
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
    exec 9>"$out.launch.lock"
    flock -n 9 || { echo "The interval scan is already running for this output" >&2; exit 1; }
fi
read -r -a julia_command <<< "${HPSC_JULIA_COMMAND:-julia +1.12.4}"
HPSC_PROFILE_POINT_PROGRESS=0 "${julia_command[@]}" --project=profiles tests/test_training60_fixed_two.jl "$fit"
HPSC_PROFILE_POINT_PROGRESS=0 "${julia_command[@]}" --project=profiles tests/test_fixed_two_seed_isolation.jl
HPSC_PROFILE_POINT_PROGRESS=0 "${julia_command[@]}" --project=profiles tests/test_fixed_two_interval_scan.jl
echo "$(date -Iseconds) START nine-parameter Model 6 confidence-interval scan"
"${julia_command[@]}" --project=profiles profiles/run_fixed_two.jl \
    "$fit" "$out" "$seeds" "$resolution" all
python3 scripts/audit_fixed_two_profiles.py "$out"
"${julia_command[@]}" --project=profile_figures profiles/render_fixed_two.jl "$out" "$out/figures"
echo "$(date -Iseconds) FINISHED interval scan, audit and figures. Review audit.csv and compare the full-range scan before final interpretation."
