#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
split=${1:-results/split6040_refined_hybrid_lac100_kg100}
root=${2:-results}
out=${3:-results/final_publication_results}
fit="$root/calibration100_hybrid_lac100_kg100"
profiles="$root/active_fixed_three_hybrid_lac100_kg100"
mkdir -p "$(dirname "$out")"
if command -v flock >/dev/null 2>&1; then
  exec 9>"${out}.lock"
  flock -n 9 || { echo "Final assembly already running." >&2; exit 1; }
fi
python3 scripts/assemble_final_results.py "$split" "$root" "$out"
julia +1.12.4 --project=figures -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=profile_figures -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=. scripts/export_figure_data.jl "$fit" "$fit/figure_data"
julia +1.12.4 --project=figures figures/export_split_combined.jl "$split" "$fit/figure_data" "$out/combined_inputs"
julia +1.12.4 --project=figures figures/render_bic_scale.jl "$out/combined_inputs" "$out/combined"
julia +1.12.4 --project=figures figures/render_bic_scale.jl "$out/combined_inputs" "$out/combined_cfd" figures/assets/Figure1_original.jpg
env HPSC_FIT_FIGURES_ONLY=1 julia +1.12.4 --project=figures figures/render.jl "$fit/figure_data" "$out/calibration"
julia +1.12.4 --project=profile_figures profile_figures/render.jl "$profiles" "$out/profiles"
for file in "$out"/profiles/model5_profiles.*; do
  [[ -f "$file" ]] || continue
  mv "$file" "${file/model5_profiles/model6_profiles}"
done
julia +1.12.4 --project=figures figures/render_media_exchange.jl "$out/media"
julia +1.12.4 --project=figures figures/render_cell_lines.jl "$out/cell_transfer" "$out/cell_recalibration"
printf 'Verified tables and regenerated figures complete. No optimization performed.\n' > "$out/COMPLETE"
echo "Final results: $out/index.html"
