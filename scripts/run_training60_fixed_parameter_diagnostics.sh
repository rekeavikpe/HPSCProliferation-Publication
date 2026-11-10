#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-42}"
export OPENBLAS_NUM_THREADS=1
export HPSC_PROFILE_POINT_PROGRESS=0
fit="${HPSC_DIAGNOSTIC_REFERENCE:-results/training60_final/reference_fit}"
seeds="${HPSC_DIAGNOSTIC_SEEDS:-results/training60_final/separate_ablations/profiles/all-active/full}"
out="${HPSC_DIAGNOSTIC_OUTPUT:-results/training60_final/profile_fixed_parameter_diagnostics}"
read -r -a julia_command <<< "${HPSC_JULIA_COMMAND:-julia +1.12.4}"
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
    exec 9>"$out.launch.lock"
    flock -n 9 || { echo "Diagnostics are already running for this output" >&2; exit 1; }
fi
"${julia_command[@]}" --project=profiles tests/test_fixed_parameter_diagnostics.jl "$fit" "$seeds"
echo "$(date -Iseconds) START K_Lac and Y_GX diagnostics; primary fits and nine-parameter results are preserved"
"${julia_command[@]}" --project=profiles profiles/run_fixed_parameter_diagnostics.jl "$fit" "$seeds" "$out"
if [[ -f "$out/SMOKE_COMPLETE" ]]; then
    echo "Smoke test only; no publication profiles produced."
    exit 0
fi
cat "$out/audit.csv"
if [[ -f "$out/NUMERICAL_REVIEW_REQUIRED" ]]; then
    echo "Numerical review required. Read audit.csv; completed point checkpoints are preserved." >&2
    exit 2
fi
echo "$(date -Iseconds) FINISHED two diagnostic curves and numerical audit. Figures will be rendered locally with Arial."
