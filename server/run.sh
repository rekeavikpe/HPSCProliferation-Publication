#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
threads="${SLURM_CPUS_PER_TASK:-${HPSC_THREADS:-42}}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
exec julia --startup-file=no --threads="$threads" --project=. scripts/calibrate.jl "${1:-results/calibration100}" "${2:-100}"
