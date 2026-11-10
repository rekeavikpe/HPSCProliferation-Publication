"""Export reporting CSVs with explicit model-number provenance; never modify fit inputs."""
import argparse
import csv
from pathlib import Path

MAPPING = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 5}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    if not source.is_dir() or output == source or source in output.parents:
        parser.error("Use an existing source and an output outside the source tree.")
    count = 0
    for path in sorted(source.rglob("*.csv")):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            if not fields:
                continue
            rows = list(reader)
        if "storage_model" in fields:
            raise ValueError(f"Already mapped reporting table: {path}")
        key = "model" if "model" in fields else "publication_model" if "publication_model" in fields else None
        if key:
            if key == "model" and "publication_model" in fields:
                raise ValueError(f"Ambiguous ID schema: {path}")
            fields = ["storage_model" if f == key else f for f in fields] + ["publication_model"]
            for row in rows:
                original = row.pop(key)
                row["storage_model"] = original
                row["publication_model"] = MAPPING[int(original)]
        dest = output / path.relative_to(source)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        count += 1
    output.mkdir(parents=True, exist_ok=True)
    (output / "README.txt").write_text("Publication reporting tables only. storage_model retains source IDs; publication_model uses the revised manuscript numbering. Do not use these exports as computational inputs. Source: " + str(source) + "\n")
    print(f"Exported {count} tables to {output}")


if __name__ == "__main__":
    main()
