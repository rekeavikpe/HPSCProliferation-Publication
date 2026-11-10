#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
out=${1:-results/david_revision_20261014}
python3 scripts/finalize_david_tables.py "$out"
python3 scripts/correct_david_reporting.py "$out"
if [[ ! -f "$out/cell_lines/COMPLETE" ]]; then
 echo "Available tables exported. Cell-line fits are pending. Full figure regeneration deferred."
 exit 0
fi
export PATH="$HOME/.juliaup/bin:$PATH"
if [[ -n "${HPSC_JULIA:-}" ]]; then J=("$HPSC_JULIA"); else J=(julia +1.12.4); fi
"${J[@]}" --project=figures figures/render_david_revision_final.jl "$out"
oxygen=${HPSC_OXYGEN_RESULTS:-results/oxygen_sensitivity_neto_properties}
if [[ -f "$oxygen/COMPLETE" ]]; then
 "${J[@]}" --project=figures figures/render_remaining_analyses_readable.jl oxygen "$oxygen"
else
 echo "Oxygen result folder unavailable. Existing oxygen figure not regenerated."
fi
