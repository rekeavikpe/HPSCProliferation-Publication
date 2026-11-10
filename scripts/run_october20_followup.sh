#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-42} OPENBLAS_NUM_THREADS=1
export HPSC_REVIEW_STARTS=${HPSC_REVIEW_STARTS:-100} HPSC_REVIEW_EVALUATIONS=${HPSC_REVIEW_EVALUATIONS:-20000} HPSC_REVIEW_ITERATIONS=${HPSC_REVIEW_ITERATIONS:-3000}
export HPSC_PROFILE_RESOLUTION=${HPSC_PROFILE_RESOLUTION:-100} HPSC_PROFILE_MAXITERS=${HPSC_PROFILE_MAXITERS:-10000} HPSC_PROFILE_POINT_PROGRESS=1
export HPSC_JULIA_COMMAND=${HPSC_JULIA_COMMAND:-'julia +1.12.4'}
read -r -a J <<< "$HPSC_JULIA_COMMAND"
source=results/david_comments_20261020
out=${1:-results/david_october20_followup}
fit=results/training60_final/reference_fit
split=results/split6040_refined_hybrid_lac100_kg100
[[ -f "$source/COMPUTATIONS_COMPLETE" ]] || { echo 'Original suite is incomplete'; exit 1; }
mkdir -p "$out" logs
if command -v flock >/dev/null 2>&1; then exec 9>"$out/launch.lock"; flock -n 9 || exit 1; fi
python3 scripts/lock_october20_followup.py "$source" "$out"
stage() {
 local name=$1; shift
 if [[ -f "$out/$name.DONE" ]]; then echo "Already completed: $name"; return; fi
 echo "$(date -Iseconds) START $name" | tee -a "$out/progress.log"
 "$@"
 touch "$out/$name.DONE"
 echo "$(date -Iseconds) FINISHED $name" | tee -a "$out/progress.log"
}
mkdir -p "$out/day1"
stage day1_export "${J[@]}" --project=population revision/export_day1_saved.jl "$fit" "$split" "$source" "$out/day1"
stage day1_figures "${J[@]}" --project=figures figures/render_day1_saved.jl "$out/day1"
stage cell_reductions python3 scripts/run_cell_reductions.py "$fit" "$source" "$out/cell_reductions"
stage reduction_figures "${J[@]}" --project=figures figures/render_cell_reductions.jl "$out/cell_reductions"
touch "$out/COMPUTATIONS_COMPLETE"
echo "Results: $out"
