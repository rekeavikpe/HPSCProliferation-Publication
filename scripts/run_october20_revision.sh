#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-42}
export HPSC_REVIEW_STARTS=${HPSC_REVIEW_STARTS:-100}
export HPSC_REVIEW_EVALUATIONS=${HPSC_REVIEW_EVALUATIONS:-20000}
export HPSC_REVIEW_ITERATIONS=${HPSC_REVIEW_ITERATIONS:-3000}
export HPSC_PROFILE_RESOLUTION=${HPSC_PROFILE_RESOLUTION:-100}
export HPSC_PROFILE_POINT_PROGRESS=1 HPSC_PROFILE_MAXITERS=10000
fit=${HPSC_FINAL_FIT:-results/training60_final/reference_fit}
split=${HPSC_SPLIT_FIT:-results/split6040_refined_hybrid_lac100_kg100}
out=${1:-results/david_comments_20261020}
export HPSC_JULIA_COMMAND=${HPSC_JULIA_COMMAND:-'julia +1.12.4'}
read -r -a J <<< "$HPSC_JULIA_COMMAND"
mkdir -p "$out/logs" logs
if command -v flock >/dev/null 2>&1; then
  exec 9>"$out/launch.lock"
  flock -n 9 || { echo 'This suite is already running.' >&2; exit 1; }
fi
printf '%s START preparation\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" | tee -a "$out/progress.log"
python3 scripts/october20_manifest.py "$fit" "$split" "$out"
for project in population profiles figures revision/statistics; do
  "${J[@]}" --project="$project" -e 'using Pkg; Pkg.instantiate()'
done
printf '%s FINISHED preparation\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" | tee -a "$out/progress.log"
stage() {
  local name=$1; shift
  if [[ -f "$out/$name.DONE" ]]; then echo "Already completed: $name"; return; fi
  printf '%s START %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$name" | tee -a "$out/progress.log"
  if "$@" 2>&1 | tee -a "$out/logs/$name.log"; then
    touch "$out/$name.DONE"
    printf '%s FINISHED %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$name" | tee -a "$out/progress.log"
  else
    printf '%s FAILED %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$name" | tee -a "$out/progress.log"
    return 1
  fi
}
stage checks "${J[@]}" --project=population revision/test_october20.jl "$fit" "$split"
stage reuse "${J[@]}" --project=population revision/run_october20.jl reuse "$fit" "$split" "$out/reuse"
stage statistics "${J[@]}" --project=revision/statistics revision/october_statistics.jl "$out/statistics"
stage oxygen_reference "${J[@]}" --project=population scripts/oxygen_sensitivity.jl "$out/oxygen_reference"
stage oxygen python3 scripts/october_oxygen.py "$out/oxygen_reference" "$out/oxygen"
stage p75 "${J[@]}" --project=population revision/run_october20.jl p75 "$fit" "$split" "$out/p75"
stage day1 "${J[@]}" --project=population revision/run_october20.jl day1 "$fit" "$split" "$out/day1"
stage cells python3 scripts/run_october_cells.py "$fit" "$out/cells"
stage figures "${J[@]}" --project=figures figures/render_october20.jl "$out"
stage report python3 scripts/report_october20.py "$out"
printf 'Computations finished. See index.html and cells/status.json for interpretation and review flags.\n' > "$out/COMPUTATIONS_COMPLETE"
echo "Results: $out/index.html"
