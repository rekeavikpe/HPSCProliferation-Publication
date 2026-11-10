#!/usr/bin/env bash
#SBATCH --job-name=hpsc-calibration
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=42
#SBATCH --time=12:00:00
#SBATCH --output=hpsc-%j.log
set -euo pipefail
# Add the account, partition, memory request and Julia module required by your lab.
# Submit from the repository root after instantiating the pinned environment.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
srun julia --startup-file=no --threads="$SLURM_CPUS_PER_TASK" --project=. scripts/calibrate.jl "results/calibration100-${SLURM_JOB_ID}" 100
