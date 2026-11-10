"""Review nine-parameter profiles without treating COMPLETE as numerical approval."""
import csv
import math
import sys
from pathlib import Path

THRESHOLD = -1.920729410347062
EXPECTED = {"K_Agg", "K_d", "K_G", "K_oG", "Y_GLac", "alpha", "mu_max", "K_tau", "K_eps"}


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def audit(root):
    intervals = read_csv(root / "profile_confidence_intervals.csv")
    if len(intervals) != 9 or {r["parameter"] for r in intervals} != EXPECTED:
        raise ValueError("Expected exactly nine profiles, including K_G")
    fixed = read_csv(root / "fixed_parameters.csv")
    if {r["parameter"] for r in fixed} != {"K_Lac", "Y_GX"}:
        raise ValueError("Unexpected fixed parameters")
    summary = read_csv(root / "mle_summary.csv")[0]
    if int(summary["free_parameters"]) != 9:
        raise ValueError("MLE did not use nine free parameters")
    rows = []
    for interval in intervals:
        name = interval["parameter"]
        points = sorted(read_csv(root / name / "profile.csv"), key=lambda p: float(p["parameter_value"]))
        nuisance = read_csv(root / name / "nuisance_parameters.csv")
        if {p["parameter"] for p in nuisance} != EXPECTED:
            raise ValueError(f"Wrong nuisance coordinate set for {name}")
        if not (root / name / "COMPLETE").is_file():
            raise ValueError(f"Incomplete profile: {name}")
        if any((root / name / marker).exists() for marker in ("RECOVERED_NEEDS_CHECK", "RECOMPUTE_REQUIRED")):
            raise ValueError(f"Recovery verification is pending: {name}")
        mle = float(interval["mle"])
        values = [float(p["profile_loglik"]) for p in points]
        if not all(math.isfinite(v) for v in values):
            raise ValueError(f"Nonfinite profile likelihood: {name}")
        issues = []
        if summary["stopping_met"] != "true":
            issues.append("reference optimization requires review")
        if summary["smoke_test"] == "true":
            issues.append("smoke test only")
        for side in ("lower", "upper"):
            branch = [p for p in points if (float(p["parameter_value"]) < mle if side == "lower" else float(p["parameter_value"]) > mle)]
            crossed = any(float(p["profile_loglik"]) <= THRESHOLD for p in branch)
            bound = any((float(p["normalized_parameter"]) <= 1e-7 if side == "lower" else float(p["normalized_parameter"]) >= 1 - 1e-7) for p in branch)
            endpoint = float(interval["package_" + side])
            if not crossed:
                issues.append(side + (" bound-limited" if bound else " incomplete"))
            if crossed and (math.isclose(endpoint, mle, rel_tol=1e-8) or not min(float(p["parameter_value"]) for p in points) <= endpoint <= max(float(p["parameter_value"]) for p in points)):
                issues.append(side + " invalid interval endpoint")
        inside = [v >= THRESHOLD for v in values]
        components = sum(v and (i == 0 or not inside[i - 1]) for i, v in enumerate(inside))
        if components > 1:
            issues.append(f"{components} confidence-region segments: inspect curve")
        dips = sum(min(a, c) - b > 0.2 and max(a, c) > -3 for a, b, c in zip(values, values[1:], values[2:]))
        if dips:
            issues.append(f"{dips} sharp dips near confidence region: inspect curve")
        if max(values) > 1e-5:
            issues.append("better likelihood than reference: investigate MLE")
        if any(p["stopping_met"] != "true" for p in points):
            issues.append("nuisance optimization requires review")
        rows.append(dict(parameter=name, mle=mle, lower=interval["package_lower"], upper=interval["package_upper"], points=len(points), status=" | ".join(issues) or "thresholds crossed on both sides"))
    with (root / "audit.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Audited {len(rows)} profiles. Fixed only K_Lac and Y_GX. K_G is free.")
    print(f"Reference likelihood improvement over saved training fit: {float(summary['saved_nll']) - float(summary['nll']):.8g}")
    for row in rows:
        print(f"{row['parameter']}: {row['status']}")
    return rows


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Supply the completed profile directory")
    audit(Path(sys.argv[1]))
