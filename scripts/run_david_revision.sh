#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH" OPENBLAS_NUM_THREADS=1
if [[ -n "${HPSC_JULIA:-}" ]]; then J=("$HPSC_JULIA"); else J=(julia +1.12.4); fi
fit=${HPSC_FINAL_FIT:-results/calibration100_hybrid_lac100_kg100}
split=${HPSC_SPLIT_FIT:-results/split6040_refined_hybrid_lac100_kg100}
out=${1:-results/david_revision_20261014}
export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-42}
export HPSC_REVIEW_STARTS=${HPSC_REVIEW_STARTS:-100}
export HPSC_REVIEW_EVALUATIONS=${HPSC_REVIEW_EVALUATIONS:-20000}
export HPSC_REVIEW_ITERATIONS=${HPSC_REVIEW_ITERATIONS:-3000}
export HPSC_CV_STARTS=${HPSC_CV_STARTS:-$HPSC_REVIEW_STARTS}
export HPSC_CV_EVALUATIONS=${HPSC_CV_EVALUATIONS:-$HPSC_REVIEW_EVALUATIONS}
export HPSC_CV_LOCAL_ITERATIONS=${HPSC_CV_LOCAL_ITERATIONS:-$HPSC_REVIEW_ITERATIONS}
unset HPSC_CV_PREPARE_ONLY HPSC_CV_FOLDS
export HPSC_CV_STORAGE_IDS=0,5 HPSC_CV_DESIGN=pooled3
export HPSC_ABLATION_SHARED=false
export HPSC_ABLATION_STARTS=${HPSC_ABLATION_STARTS:-$HPSC_REVIEW_STARTS}
export HPSC_ABLATION_EVALUATIONS=${HPSC_ABLATION_EVALUATIONS:-$HPSC_REVIEW_EVALUATIONS}
export HPSC_ABLATION_ITERATIONS=${HPSC_ABLATION_ITERATIONS:-$HPSC_REVIEW_ITERATIONS}
mkdir -p "$out/logs"
if command -v flock >/dev/null 2>&1; then
 exec 9>"$out/launch.lock"
 flock -n 9 || { echo "This analysis suite is already running." >&2; exit 1; }
fi
rm -f "$out/COMPLETE"
# Fail before starting costly work if the saved reference inputs are unavailable.
for f in "$fit/calibration.toml" "$fit/observations.csv" "$fit/ranking.csv" "$fit/model0/selected_parameters.csv" "$fit/model5/selected_parameters.csv" "$split/split.csv"; do
 [[ -f "$f" ]] || { echo "Missing required saved result: $f" >&2; exit 1; }
done
for id in 0 1 2 3 4 5 6; do
 [[ -f "$split/model$id/selected_parameters.csv" && -f "$split/model$id/summary.csv" ]] || { echo "Missing split fit model$id" >&2; exit 1; }
done
for project in population figures profiles profile_figures revision/statistics; do
 "${J[@]}" --project="$project" -e 'using Pkg; Pkg.instantiate()'
done
run_stage() {
 local name=$1; shift
 printf '%s START %s\n' "$(date +%Y-%m-%dT%H:%M:%S)" "$name" | tee -a "$out/progress.log"
 "$@" 2>&1 | tee -a "$out/logs/$name.log"
 printf '%s FINISHED %s\n' "$(date +%Y-%m-%dT%H:%M:%S)" "$name" | tee -a "$out/progress.log"
}
run_stage reference_audit "${J[@]}" --project=population revision/preflight.jl "$fit" "$split"
run_stage checks "${J[@]}" --project=population revision/test_revision.jl
run_stage selection "${J[@]}" --project=population revision/run_david_analyses.jl selection "$fit" "$split" "$out/selection"
run_stage rpm_statistics "${J[@]}" --project=revision/statistics revision/rpm_statistics.jl "$out/rpm_statistics"
run_stage pooled_validation bash scripts/run_pooled_validation.sh "$out/pooled_validation"
run_stage exposure_selection "${J[@]}" --threads="$JULIA_NUM_THREADS" --project=population revision/run_david_analyses.jl exposure "$fit" "$split" "$out/exposure_selection"
if [[ ! -f "$fit/figure_data/optimizer_bic.csv" ]]; then
 run_stage reference_exports "${J[@]}" --project=. scripts/export_figure_data.jl "$fit" "$fit/figure_data"
fi
run_stage separate_ablations "${J[@]}" --threads="$JULIA_NUM_THREADS" --project=population population/run_routed_ablations.jl "$fit" "$out/separate_ablations"
run_stage ablation_figures "${J[@]}" --project=figures figures/render_routed_ablations.jl "$out/separate_ablations"
run_stage cell_lines "${J[@]}" --threads="$JULIA_NUM_THREADS" --project=population revision/run_david_analyses.jl cells "$fit" "$split" "$out/cell_lines"
run_stage summary_figures "${J[@]}" --project=figures figures/render_david_revision.jl "$out"
# Both profile sets are included. Each individual profile is independently resumable.
for mode in all-active conditional; do
 for scenario in baseline shear_cell shear_aggregate shear_both edr_cell edr_aggregate edr_both full; do
  input="$out/separate_ablations/profile_inputs/$scenario"
  dest="$out/separate_ablations/profiles/$mode/$scenario"
  run_stage "profile_${mode}_${scenario}" "${J[@]}" --threads="$JULIA_NUM_THREADS" --project=profiles profiles/run_routed.jl "$input" "$dest" "$mode" "${HPSC_PROFILE_RESOLUTION:-100}" all
  run_stage "profile_figure_${mode}_${scenario}" "${J[@]}" --project=profile_figures profile_figures/render.jl "$dest" "$dest/figures"
 done
done
printf 'All scheduled analyses and plots completed. Oxygen CFD awaits Pouyan.\n' > "$out/COMPLETE"
echo "Results: $out/index.html"
