#!/usr/bin/env python3
"""Merge independently saved CESR orbit-response rho-sweep chunks."""

from __future__ import annotations

import argparse
import csv
import json
import tomllib
from pathlib import Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = arguments()
    root = args.root.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    chunks = sorted(path for path in root.glob("chunk_*") if path.is_dir())
    if len(chunks) != 5:
        raise RuntimeError(f"Expected five chunk directories in {root}, found {len(chunks)}")

    trial_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    metadata: list[dict[str, object]] = []
    baseline_written = False
    for chunk in chunks:
        trials = read_csv(chunk / "rho_sweep_trial_errors.csv")
        for row in trials:
            is_baseline = row["scenario"] == "baseline"
            if is_baseline and baseline_written:
                continue
            row["sample_id"] = len(trial_rows)
            row["source_chunk"] = chunk.name
            trial_rows.append(row)
            baseline_written = baseline_written or is_baseline

        summaries = read_csv(chunk / "rho_sweep_summary.csv")
        for row in summaries:
            if float(row["rho"]) == 0.0:
                continue
            row["source_chunk"] = chunk.name
            summary_rows.append(row)
        with (chunk / "rho_sweep_metadata.toml").open("rb") as stream:
            metadata.append(tomllib.load(stream))

    scenario_order = {"all": 0, "horizontal": 1, "vertical": 2}
    summary_rows.sort(key=lambda row: (scenario_order[str(row["scenario"])], float(row["rho"])))
    trial_fields = list(trial_rows[0])
    summary_fields = list(summary_rows[0])
    write_csv(output / "rho_sweep_trial_errors.csv", trial_rows, trial_fields)
    write_csv(output / "rho_sweep_summary.csv", summary_rows, summary_fields)

    combined = {
        "format": "cesr-orbit-response-rho-sweep-chunked-v1",
        "chunks": [chunk.name for chunk in chunks],
        "chunk_metadata": metadata,
        "total_unique_samples": len(trial_rows),
        "positive_rho_values": sorted({float(row["rho"]) for row in summary_rows}),
        "trials_per_positive_rho_scenario": 600,
        "converged_count": sum(int(item["converged_count"]) - 1 for item in metadata) + 1,
        "failed_count": sum(int(item["failed_count"]) for item in metadata),
        "fallback_count": sum(int(item["fallback_count"]) for item in metadata),
        "exact_physics_seconds_sum": sum(float(item["exact_physics_seconds"]) for item in metadata),
        "response_evaluation_seconds_sum": sum(float(item["response_evaluation_seconds"]) for item in metadata),
        "wall_time_is_not_sum_because_four_chunks_ran_concurrently": True,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "rho_sweep_metadata.json").write_text(
        json.dumps(combined, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Merged {len(trial_rows)} unique samples into {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
