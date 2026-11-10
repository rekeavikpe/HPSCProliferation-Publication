#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
export HPSC_PROFILE_FULL_RANGE=1
export HPSC_PROFILE_BRANCH_SEEDS=1
export HPSC_PROFILE_POINT_PROGRESS=1
root=${1:-results/training60_final}
out="$root/profile_alpha_branch"
mkdir -p "$out"
if command -v flock >/dev/null 2>&1; then
  exec 9>"$out/run.lock"
  flock -n 9 || { echo "Alpha branch check already running" >&2; exit 1; }
fi
export HPSC_PROFILE_SEED_DIR="$root/profile_multistart/all-active/shear_aggregate"
test -f "$HPSC_PROFILE_SEED_DIR/alpha/nuisance_parameters.csv"
julia +1.12.4 --project=profiles -e 'using Pkg; Pkg.instantiate()'
echo "$(date -Iseconds) START shear_aggregate alpha branch check" | tee -a "$out/progress.log"
julia +1.12.4 --threads="${JULIA_NUM_THREADS:-42}" --project=profiles profiles/run_routed.jl \
 "$root/separate_ablations/profile_inputs/shear_aggregate" "$out/all-active/shear_aggregate" all-active 400 alpha
julia +1.12.4 --project=profile_figures profile_figures/render.jl \
 "$out/all-active/shear_aggregate" "$out/all-active/shear_aggregate/figures"
python3 scripts/audit_profile_repairs.py "$out"
echo "$(date -Iseconds) FINISHED alpha branch check. Inspect saved curve and interval." | tee -a "$out/progress.log"
