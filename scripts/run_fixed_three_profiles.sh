#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
fit=${1:-results/calibration100}
out=${2:-results/active_fixed_three}
threads=${JULIA_NUM_THREADS:-42}
julia +1.12.4 --project=profiles -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=profile_figures -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --threads="$threads" --project=profiles profiles/run.jl "$fit" "$out" 5 100 all --fix-three
julia +1.12.4 --project=profile_figures profile_figures/render.jl "$out" "$out/figures"
julia +1.12.4 --project=profiles profiles/compare_additional_fixed.jl "$out/review" "$out"
python3 profiles/report_fixed_three.py "$out"
