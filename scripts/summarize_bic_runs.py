"""Summarize initial optimizer outcomes without dropping poor fits."""
import csv
import math
import sys
from pathlib import Path

source = Path(sys.argv[1])
rows = list(csv.DictReader((source / "optimizer_bic.csv").open()))
# The exported source metadata records the number of independent initial starts.
import re
match = re.search(r"^initial_starts_per_model\s*=\s*(\d+)", (source / "source.toml").read_text(), re.M)
if match is None:
    raise ValueError("Missing initial start count")
nstarts = int(match.group(1))
rows = [r for r in rows if int(r["run"].removeprefix("run")) <= nstarts]
writer = csv.writer(sys.stdout)
writer.writerow(["model", "initial_runs", "finite_runs", "stopping_met", "best_bic", "median_bic", "worst_bic", "within_2_of_model_best", "more_than_10_above_model_best"])
from statistics import median
for model in range(7):
    group = [r for r in rows if int(r["model"]) == model]
    values = [float(r["bic"]) for r in group if math.isfinite(float(r["bic"]))]
    if not values:
        raise ValueError(f"No finite BIC for model {model}")
    best = min(values)
    writer.writerow([model, len(group), len(values), sum(r["stopping_met"].lower() == "true" for r in group), best, median(values), max(values), sum(v-best <= 2 for v in values), sum(v-best > 10 for v in values)])
