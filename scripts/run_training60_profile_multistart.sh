#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.juliaup/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
export HPSC_PROFILE_FULL_RANGE=1
root=${1:-results/training60_final}
out="$root/profile_multistart"
mkdir -p "$out"
if command -v flock >/dev/null 2>&1; then
  exec 9>"$out/run.lock"
  flock -n 9 || { echo "Profile repair is already running" >&2; exit 1; }
fi
julia +1.12.4 --project=profiles -e 'using Pkg; Pkg.instantiate()'
python3 - "$root" > "$out/jobs.tsv" <<'PYTHON'
import csv,sys
from pathlib import Path
from collections import defaultdict
root=Path(sys.argv[1]);jobs=defaultdict(list)
for row in csv.DictReader(open('config/training60_profile_repairs.csv')):
    old=root/'separate_ablations/profiles'/row['mode']/row['scenario']/row['parameter']/'profile.csv'
    if not old.is_file(): raise SystemExit(f'Missing source profile: {old}')
    jobs[row['mode'],row['scenario']].append(row['parameter'])
for (mode,scenario),names in jobs.items(): print(mode,scenario,','.join(names),sep='\t')
PYTHON
while IFS=$'\t' read -r mode scenario parameters; do
  echo "$(date -Iseconds) START $mode $scenario $parameters" | tee -a "$out/progress.log"
  export HPSC_PROFILE_SEED_DIR="$root/profile_repairs/$mode/$scenario"
  julia +1.12.4 --threads="${JULIA_NUM_THREADS:-42}" --project=profiles profiles/run_routed.jl \
    "$root/separate_ablations/profile_inputs/$scenario" "$out/$mode/$scenario" "$mode" "${HPSC_PROFILE_RESOLUTION:-200}" "$parameters"
  julia +1.12.4 --project=profile_figures profile_figures/render.jl "$out/$mode/$scenario" "$out/$mode/$scenario/figures"
  echo "$(date -Iseconds) FINISHED $mode $scenario" | tee -a "$out/progress.log"
done < "$out/jobs.tsv"
python3 scripts/audit_profile_repairs.py "$out"
echo "Targeted computations finished. Interval and likelihood verification still required."
