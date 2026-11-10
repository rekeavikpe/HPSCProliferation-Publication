#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
root=${1:-results/ablation_eight_shared_hybrid_lac100_kg100}
[[ -f "$root/COMPLETE" ]] || { echo "Complete eight-version calibration first." >&2; exit 1; }
if command -v flock >/dev/null 2>&1; then
  exec 9>"${root}.profiles.lock"
  flock -n 9 || { echo "Profiles already running." >&2; exit 1; }
fi
julia +1.12.4 --project=profiles -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=profile_figures -e 'using Pkg; Pkg.instantiate()'
for mode in all-active conditional; do
  for scenario in baseline shear_cell shear_aggregate shear_both edr_cell edr_aggregate edr_both full; do
    input="$root/profile_inputs/$scenario"
    out="$root/profiles/$mode/$scenario"
    julia +1.12.4 --threads="${JULIA_NUM_THREADS:-42}" --project=profiles profiles/run_routed.jl "$input" "$out" "$mode" "${HPSC_PROFILE_RESOLUTION:-100}" all
    julia +1.12.4 --project=profile_figures profile_figures/render.jl "$out" "$out/figures"
  done
done
echo "Both profile sets complete: $root/profiles"
