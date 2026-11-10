"""Copy the failed interval run into a new, guarded recovery output.

Never changes the source run. Only the known defective runner hash may differ.
Completed profiles are preserved but must pass the runner's seed rechecks.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
try:
    import tomllib
except ModuleNotFoundError:
    tomllib = None

REPO = Path(__file__).resolve().parents[1]
LEGACY_RUNNER_SHA = "09c65d795fcca26514d036bfa1cc0f42abbb65bafbf92984a8da7f6ecfb66832"
RUNNER = "profiles/run_fixed_two.jl"
HELPER = "profiles/FixedTwoSeeds.jl"
PARAMETERS = {"K_Agg", "K_d", "K_G", "K_oG", "Y_GLac", "alpha", "mu_max", "K_tau", "K_eps"}
ROOT_FILES = ("settings.toml", "calibration.toml", "fixed_parameters.csv", "mle_parameters.csv", "mle_summary.csv")
PROFILE_FILES = ("COMPLETE", "profile.csv", "nuisance_parameters.csv", "optimizer_calls.csv", "interval.csv")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def overlapping(a, b):
    return a == b or a in b.parents or b in a.parents


def ensure_stopped(source):
    # Look at the output argument, not a generic Julia name. The separate
    # full-range job is explicitly allowed to continue.
    listing = subprocess.run(["ps", "-eo", "pid=,args="], check=True, capture_output=True, text=True).stdout
    for line in listing.splitlines():
        fields = line.split()
        if len(fields) < 6 or "julia" not in Path(fields[1]).name:
            continue
        for position, token in enumerate(fields):
            if token.endswith("profiles/run_fixed_two.jl") and len(fields) > position + 2:
                candidate = Path(fields[position + 2])
                if candidate.resolve() == source:
                    raise ValueError(f"Interval process {fields[0]} is still active. Do not recover a changing source.")


def recover(source, destination, repo=REPO):
    source, destination, repo = source.resolve(), destination.resolve(), repo.resolve()
    if overlapping(source, destination):
        raise ValueError("Recovery output must be separate from the failed run")
    if destination.exists():
        raise ValueError("Recovery output already exists. Resume it with run_training60_fixed_two_recovery.sh instead.")
    manifest = read_csv(source / "source_hashes.csv")
    paths = [row["path"] for row in manifest]
    if len(paths) != len(set(paths)):
        raise ValueError("Duplicate source hashes")
    runner = repo / RUNNER
    updated = []
    found = False
    for row in manifest:
        path = Path(row["path"])
        if path.resolve() == runner:
            if row["sha256"] != LEGACY_RUNNER_SHA:
                raise ValueError("Recovery only supports the known failed runner version")
            found = True
            updated.append(dict(path=str(path), sha256=digest(path)))
            updated.append(dict(path=str(path.parent / "FixedTwoSeeds.jl"), sha256=digest(repo / HELPER)))
        else:
            if digest(path) != row["sha256"]:
                raise ValueError(f"Protected input changed: {path}")
            updated.append(row)
            if overlapping(destination, path):
                raise ValueError(f"Recovery output overlaps input: {path}")
    if not found:
        raise ValueError("Original runner hash not found")
    # The Julia runner also checks every setting exactly on resume. Python's
    # TOML check is supplemental because some server Python versions lack it.
    if tomllib:
        settings = tomllib.loads((source / "settings.toml").read_text())
        if settings["full_range_repair"] or settings["profile_scope"] != "fixed-two":
            raise ValueError("Expected the fixed-two confidence-interval run")
    completed = sorted(p for p in PARAMETERS if (source / p / "COMPLETE").is_file())
    selected = [source / f for f in ROOT_FILES] + [source / "source_hashes.csv"]
    selected += [source / p / f for p in completed for f in PROFILE_FILES]
    snapshot = {str(path.relative_to(source)): digest(path) for path in selected}
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=destination.name + ".preparing-", dir=destination.parent))
    try:
        for f in ROOT_FILES:
            shutil.copy2(source / f, staging / f)
        shutil.copy2(source / "source_hashes.csv", staging / "source_hashes_before_recovery.csv")
        for p in completed:
            (staging / p).mkdir()
            for f in PROFILE_FILES:
                shutil.copy2(source / p / f, staging / p / f)
            (staging / p / "RECOVERED_NEEDS_CHECK").write_text("Verify saved likelihoods and re-optimize confidence-limit bracket points before reuse.\n")
        with (staging / "source_hashes.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["path", "sha256"])
            writer.writeheader()
            writer.writerows(updated)
        record = dict(source=str(source), completed_profiles=completed, source_snapshot=snapshot,
                      reason="Isolate each parameter's seed coordinates and vectors from parallel jobs",
                      original_runner_sha256=LEGACY_RUNNER_SHA, corrected_runner_sha256=digest(runner))
        (staging / "recovery.json").write_text(json.dumps(record, indent=2) + "\n")
        for relative, expected in snapshot.items():
            if digest(source / relative) != expected:
                raise ValueError("Source changed while making recovery copy")
            copy = staging / ("source_hashes_before_recovery.csv" if relative == "source_hashes.csv" else relative)
            if digest(copy) != expected:
                raise ValueError(f"Recovery copy differs: {relative}")
        os.rename(staging, destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    print(f"Preserved {len(completed)} completed profiles: {', '.join(completed) or 'none'}")
    print(f"Original run unchanged: {source}")
    print(f"Recovery output: {destination}")
    print("Completed profiles will be checked before reuse. Missing profiles will run with isolated seed banks.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: recover_fixed_two_intervals.py FAILED_INTERVAL_OUTPUT NEW_OUTPUT")
    try:
        source, destination = map(Path, sys.argv[1:])
        ensure_stopped(source.resolve())
        recover(source, destination)
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(str(error))
