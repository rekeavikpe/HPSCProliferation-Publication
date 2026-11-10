#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
fit=${1:-results/calibration100}
threads=${JULIA_NUM_THREADS:-42}
out=${2:-results/inactive_lactate_yield}
julia +1.12.4 --project=profiles -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --project=profile_figures -e 'using Pkg; Pkg.instantiate()'
julia +1.12.4 --threads="$threads" --project=profiles profiles/run.jl "$fit" "${out}_active_fit" 5 100 all --fix-yield --mle-only
julia +1.12.4 --threads="$threads" --project=profiles profiles/run.jl "$fit" "$out" 5 100 all --inactive-lactate --fix-yield
julia +1.12.4 --project=profile_figures profile_figures/render.jl "$out" "$out/figures"
julia +1.12.4 --threads="$threads" --project=profiles profiles/prediction_sensitivity.jl "${out}_active_fit" "$out"
julia +1.12.4 --project=profile_figures profile_figures/prediction_sensitivity.jl "$out"
