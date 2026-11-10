#!/usr/bin/env bash
# Recover only the failed interval pass. Never stops or writes the full scan.
set -euo pipefail
cd "$(dirname "$0")/.."
failed="${HPSC_FIXED2_FAILED_OUTPUT:-results/training60_final/profile_fixed_two_intervals}"
out="${HPSC_FIXED2_INTERVAL_OUTPUT:-results/training60_final/profile_fixed_two_intervals_recovered}"
if [[ -d "$out" ]]; then
    python3 - "$failed" "$out" <<'PY'
import json, sys
from pathlib import Path
source, target = (Path(p).resolve() for p in sys.argv[1:])
record = json.loads((target / "recovery.json").read_text())
if Path(record["source"]).resolve() != source:
    raise SystemExit("Existing recovery output belongs to a different source")
print("Resuming the existing recovery output. Exact source and setting checks follow.")
PY
else
    python3 scripts/recover_fixed_two_intervals.py "$failed" "$out"
fi
export HPSC_FIXED2_INTERVAL_OUTPUT="$out"
exec bash scripts/run_training60_fixed_two_intervals.sh
