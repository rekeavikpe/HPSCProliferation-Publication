"""Import the six ESI-017 conditions without replacing missing observations."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import mean, stdev

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = "mTesR1_VWBR_100mL_vessel_Data_for_Amon_Lab_March_10_2025.xlsx"
STATES = {
    "cells": (4, "cells/mL", "ExperimentalCells_per_mL"),
    "glucose": (24, "mM", "ExperimentalGlucose_mM"),
    "lactate": (34, "mM", "ExperimentalLactate_mM"),
    "aggregate": (42, "um", "AggregateSize_um"),
}


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def classify(value):
    if value is None or value == "No data":
        return "missing", "", ""
    if isinstance(value, (int, float)) and math.isfinite(value):
        return "measured", value, ""
    if isinstance(value, str) and value.startswith("<"):
        return "left_censored", "", float(value[1:])
    if isinstance(value, str) and value.startswith("HARVESTED"):
        return "not_observed_after_harvest", "", ""
    raise ValueError(f"Unrecognized source value: {value!r}")


def verify_sources(root):
    for entry in json.loads((root / "config/source_manifest.json").read_text()):
        path = root / entry["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"Source changed: {entry['path']}")


def extract(root):
    verify_sources(root)
    workbook = openpyxl.load_workbook(root / "data/source" / WORKBOOK, data_only=False)
    conditions = json.loads((root / "config/main_conditions.json").read_text())
    records, series, feeding = [], [], []
    for block in conditions:
        c, sheet = block["condition"], workbook[block["sheet"]]
        columns = range(block["first_column"], block["last_column"] + 1)
        labels = [sheet.cell(3, col).value for col in columns]
        if labels != block["expected_labels"]:
            raise ValueError(f"Source headers changed in condition {c}")
        for col, label in zip(columns, labels):
            letter = openpyxl.utils.get_column_letter(col)
            # Condition and source column prevent merging reused or uncertain labels.
            series_id = f"ESI017_C{c}_{letter}"
            legacy = not (c == 5 and label in {"R1P1", "R2P1", "R2P2"})
            identity_issue = ""
            if c == 2 and label in {"L34", "L35"}:
                identity_issue = "density_sheet_header_conflict"
            if label == "L35":
                identity_issue += ";label_reused_across_conditions"
            series.append(dict(series_id=series_id, source_label=label, condition=c,
                source_sheet=sheet.title, source_column=letter,
                volume_ml=block["volume_ml"], rpm=block["rpm"],
                seeding_cells_ml=block["seeding_cells_ml"],
                separate_bioreactor="author_confirmed", batch_id="",
                legacy_calibration_series=legacy, identity_issue=identity_issue.strip(";")))
            for day in range(block["last_day"] + 1):
                for state, (offset, unit, _) in STATES.items():
                    if state == "aggregate" and day == 0:
                        continue  # No aggregate day-0 block exists in this source.
                    row = (41 + day) if state == "aggregate" and c == 4 else offset + day
                    day_column = {1: 2, 2: 11, 3: 18, 4: 2, 5: 2, 6: 24}[c]
                    if not (c == 4 and state == "aggregate") and sheet.cell(row, day_column).value != day:
                        raise ValueError(f"Day label mismatch: {sheet.title}, row {row}")
                    cell = sheet.cell(row, col)
                    status, value, limit = classify(cell.value)
                    if day == 0 and state == "cells" and status == "measured":
                        status = "nominal_inoculum"
                    if c == 5 and label == "R2P1" and day == 5:
                        status = "not_observed_after_harvest"
                    issue = "L29_day6_conflicting_value" if (c, label, day, state) == (2, "L29", 6, "aggregate") else ""
                    records.append(dict(series_id=series_id, source_label=label, condition=c,
                        day=day, time_h=day * 24, state=state, value=value, unit=unit,
                        status=status, censor_limit=limit, source_file=WORKBOOK,
                        source_sheet=sheet.title, source_cell=cell.coordinate,
                        source_raw="" if cell.value is None else cell.value,
                        day_basis="positional_legacy_alignment" if c == 4 and state == "aggregate" else "source_day",
                        sample_phase="before_exchange", legacy_calibration_series=legacy,
                        issue=issue))
        events = [(4, .5)] if c == 1 else [(d, .5 if d == 2 else .8) for d in range(2, block["last_day"])]
        for day, fraction in events:
            feeding.append(dict(condition=c, day=day, time_h=day * 24,
                exchange_fraction=fraction, observation_side="before_exchange",
                time_basis="protocol_day_in_hours", source="manuscript Methods, lines 117-119",
                additional_handling="centrifugation" if c == 6 and day in (2, 3) else "gravity_settling"))
    return records, series, feeding


def summarize(records, legacy_only):
    groups = {}
    for row in records:
        if row["day"] == 0 or row["status"] != "measured":
            continue
        if legacy_only and not row["legacy_calibration_series"]:
            continue
        groups.setdefault((row["condition"], row["day"], row["state"]), []).append(row)
    result = []
    for (condition, day, state), rows in sorted(groups.items()):
        values = [r["value"] for r in rows]
        result.append(dict(condition=condition, day=day, state=state,
            mean=mean(values), sd=stdev(values) if len(values) > 1 else "",
            n_vessels=len(values), provisional_value_count=sum(bool(r["issue"]) for r in rows)))
    return result


def reconcile(root, records):
    source = {}
    for row in records:
        if row["day"] > 0 and row["legacy_calibration_series"]:
            source.setdefault((row["condition"], row["day"], row["series_id"]), {})[row["state"]] = row
    with (root / "data/source/Combined100mLand500mLDataAllRPM.csv").open() as stream:
        old = list(csv.DictReader(stream))
    differences, mapping = [], []
    for line, row in enumerate(old, 2):
        c, day = int(row["ConditionSet"]), int(row["t (days)"])
        candidates = [(key, values) for key, values in source.items() if key[:2] == (c, day)]
        def score(candidate):
            values = candidate[1]
            return sum(values[s]["status"] == "measured" and
                abs(float(row[meta[2]]) - values[s]["value"]) <=
                (.00501 if s == "aggregate" else .51 if s == "cells" else .000051)
                for s, meta in STATES.items())
        best_score = max(map(score, candidates))
        best = [candidate for candidate in candidates if score(candidate) == best_score]
        if len(best) != 1 or best_score < 2:
            raise ValueError(f"Ambiguous legacy match at row {line}")
        key, values = best[0]
        mapping.append(dict(csv_row=line, condition=c, day=day, series_id=key[2],
            source_label=values["cells"]["source_label"], matching_states=best_score,
            mapping_basis="unique_best_value_match"))
        for state, meta in STATES.items():
            new = values[state]
            previous = float(row[meta[2]])
            if new["status"] != "measured":
                change = "remove_unmeasured_replacement"
            elif previous != new["value"]:
                change = "restore_source_precision_or_transcription"
            else:
                continue
            differences.append(dict(csv_row=line, series_id=key[2], day=day, state=state,
                previous_value=previous, imported_value=new["value"], change=change,
                source_sheet=new["source_sheet"], source_cell=new["source_cell"]))
    keys = [(r["series_id"], r["day"]) for r in mapping]
    if len(set(keys)) != len(keys):
        raise ValueError("Multiple legacy rows mapped to the same vessel/day")
    return mapping, differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "data/processed")
    args = parser.parse_args()
    records, series, feeding = extract(ROOT)
    mapping, differences = reconcile(ROOT, records)
    args.out.mkdir(parents=True, exist_ok=True)
    for name, rows in {
        "measurements.csv": records, "series.csv": series, "feeding.csv": feeding,
        "summary_legacy_cohort.csv": summarize(records, True),
        "summary_available_measurements.csv": summarize(records, False),
        "legacy_mapping.csv": mapping, "legacy_changes.csv": differences,
    }.items():
        write_csv(args.out / name, rows)
    print(f"Imported {len(records)} measurement slots from {len(series)} vessel series; reconciled {len(mapping)} legacy rows.")


if __name__ == "__main__":
    main()
