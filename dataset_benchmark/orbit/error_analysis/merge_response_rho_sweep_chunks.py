#!/usr/bin/env python3
"""Merge independently saved CESR orbit-response rho-sweep chunks."""

from __future__ import annotations

import argparse
import csv
import json
import tomllib
from pathlib import Path


# These fields define the physical/solver contract of a merge.  Chunk-specific
# rho values, dates, and timing are intentionally excluded; the chunks are
# expected to partition those dimensions.
METADATA_AGREEMENT_KEYS = (
    "ring_id",
    "lattice_path",
    "branch",
    "rf_on",
    "scibmad_version",
    "control_count",
    "control_names",
    "observable_count",
    "observable_labels",
    "detector_count",
    "detector_count_per_plane",
    "input_csv",
    "detector_response_csv",
    "closed_orbit_response_csv",
    "base_kick_rad",
    "trials_per_rho_scenario",
    "seed",
    "exact_reference",
    "reltol",
    "abstol",
    "maxiter",
    # The response matrices are a paired, method-scoped input to every
    # chunk.  Comparing only the ring and output dimensions would allow a
    # mixed GTPSA/finite-difference merge, or two independently generated
    # response pairs, to look like one experiment.
    "response_method",
    "response_step_rad",
    "response_controls_per_batch",
    "response_pair_id",
    "response_reltol",
    "response_abstol",
    "response_maxiter",
    "rho_definition",
    "direction_distribution",
    "direction_reuse",
    "scenarios",
)


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


def _canonical(value: object) -> str:
    """Make TOML scalar/list/table values comparable without hashing."""

    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def validate_metadata(metadata: list[tuple[str, dict[str, object]]]) -> None:
    """Reject a merge whose chunks use different rings or solver contracts."""

    if not metadata:
        raise RuntimeError("No chunk metadata was loaded")
    reference_name, reference = metadata[0]
    for key in METADATA_AGREEMENT_KEYS:
        if key not in reference:
            raise RuntimeError(f"{reference_name} metadata is missing required merge field '{key}'")
    pair_id = str(reference["response_pair_id"])
    if not pair_id:
        raise RuntimeError(
            f"{reference_name} metadata has an empty response_pair_id; "
            "latest-ring chunks must use one published response pair"
        )
    for chunk_name, current in metadata[1:]:
        for key in METADATA_AGREEMENT_KEYS:
            if key not in current:
                raise RuntimeError(f"{chunk_name} metadata is missing required merge field '{key}'")
            if _canonical(current[key]) != _canonical(reference[key]):
                raise RuntimeError(
                    f"Cannot merge {chunk_name}: metadata field '{key}' differs from "
                    f"{reference_name}"
                )
        if not str(current["response_pair_id"]):
            raise RuntimeError(f"{chunk_name} metadata has an empty response_pair_id")


def _metadata_rhos(metadata: dict[str, object]) -> set[float]:
    values = metadata.get("rho_values")
    if not isinstance(values, list):
        raise RuntimeError("Chunk metadata must contain a rho_values list")
    rhos = {float(value) for value in values}
    if 0.0 not in rhos:
        raise RuntimeError("Each rho chunk must include the shared rho=0 baseline")
    if any(value < 0.0 for value in rhos):
        raise RuntimeError("Chunk metadata contains a negative rho")
    return rhos


def _validate_chunk_rows(
    chunk_name: str,
    rows: list[dict[str, str]],
    metadata: dict[str, object],
) -> tuple[set[tuple[str, float, int]], list[dict[str, str]]]:
    """Validate row keys before de-duplicating the shared baseline.

    A merge is only meaningful when every chunk is a disjoint partition of
    scenario/rho/trial cells.  The old merger silently accepted overlapping
    cells and silently kept whichever chunk happened to sort first.
    """

    if not rows:
        raise RuntimeError(f"{chunk_name} has no rho-sweep trial rows")
    expected_rhos = _metadata_rhos(metadata)
    expected_scenarios = {str(value) for value in metadata.get("scenarios", [])}
    if not expected_scenarios:
        raise RuntimeError(f"{chunk_name} metadata has no scenarios")
    expected_trials = int(metadata.get("trials_per_rho_scenario", 0))
    if expected_trials < 1:
        raise RuntimeError(f"{chunk_name} metadata has an invalid trial count")

    keys: set[tuple[str, float, int]] = set()
    baselines: list[dict[str, str]] = []
    grouped: dict[tuple[str, float], list[int]] = {}
    for row_number, row in enumerate(rows, start=2):
        try:
            scenario = str(row["scenario"])
            rho = float(row["rho"])
            trial_id = int(row["trial_id"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RuntimeError(f"{chunk_name} row {row_number} has invalid scenario/rho/trial fields") from exc
        if scenario == "baseline":
            if rho != 0.0 or trial_id != 0:
                raise RuntimeError(f"{chunk_name} baseline row must have rho=0 and trial_id=0")
            baselines.append(row)
            continue
        if scenario not in expected_scenarios:
            raise RuntimeError(f"{chunk_name} contains unknown scenario '{scenario}'")
        if rho <= 0.0 or rho not in expected_rhos:
            raise RuntimeError(f"{chunk_name} row {row_number} uses rho={rho} absent from metadata")
        if not 1 <= trial_id <= expected_trials:
            raise RuntimeError(
                f"{chunk_name} row {row_number} has trial_id={trial_id}; "
                f"expected 1..{expected_trials}"
            )
        key = (scenario, rho, trial_id)
        if key in keys:
            raise RuntimeError(f"{chunk_name} contains duplicate cell {key}")
        keys.add(key)
        grouped.setdefault((scenario, rho), []).append(trial_id)

    if len(baselines) != 1:
        raise RuntimeError(f"{chunk_name} must contain exactly one shared baseline row")
    expected_positive = expected_rhos - {0.0}
    for scenario in expected_scenarios:
        for rho in expected_positive:
            cell = (scenario, rho)
            trial_ids = sorted(grouped.get(cell, []))
            if trial_ids != list(range(1, expected_trials + 1)):
                raise RuntimeError(
                    f"{chunk_name} does not contain the complete trial cell "
                    f"{scenario}, rho={rho:g}"
                )
    return keys, baselines


def main() -> int:
    args = arguments()
    root = args.root.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    chunks = sorted(path for path in root.glob("chunk_*") if path.is_dir())
    if not chunks:
        raise RuntimeError(f"No chunk directories found in {root}")

    trial_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    metadata: list[tuple[str, dict[str, object]]] = []
    baseline_written = False
    all_trial_keys: set[tuple[str, float, int]] = set()
    for chunk in chunks:
        trials = read_csv(chunk / "rho_sweep_trial_errors.csv")
        metadata_path = chunk / "rho_sweep_metadata.toml"
        if not metadata_path.is_file():
            raise RuntimeError(f"Missing metadata sidecar in {chunk}")
        with metadata_path.open("rb") as stream:
            chunk_metadata = tomllib.load(stream)
        chunk_keys, _ = _validate_chunk_rows(
            chunk.name,
            trials,
            chunk_metadata,
        )
        overlap = all_trial_keys.intersection(chunk_keys)
        if overlap:
            raise RuntimeError(
                f"Chunk {chunk.name} overlaps earlier positive-rho cells: "
                f"{sorted(overlap)[:3]}"
            )
        all_trial_keys.update(chunk_keys)
        for row in trials:
            is_baseline = row["scenario"] == "baseline"
            if is_baseline and baseline_written:
                continue
            # Keep the original per-chunk id for auditability while assigning
            # a unique merged id below.
            row["source_sample_id"] = row["sample_id"]
            row["sample_id"] = len(trial_rows)
            row["source_chunk"] = chunk.name
            trial_rows.append(row)
            baseline_written = baseline_written or is_baseline

        summaries = read_csv(chunk / "rho_sweep_summary.csv")
        summary_keys = {
            (str(row["scenario"]), float(row["rho"]))
            for row in summaries
            if float(row["rho"]) > 0.0
        }
        expected_summary_keys = {
            (scenario, rho)
            for scenario, rho, _ in chunk_keys
        }
        if summary_keys != expected_summary_keys:
            raise RuntimeError(
                f"{chunk.name} summary cells do not match its trial rows"
            )
        for row in summaries:
            if float(row["rho"]) == 0.0:
                continue
            row["source_chunk"] = chunk.name
            summary_rows.append(row)
        metadata.append((chunk.name, chunk_metadata))

    validate_metadata(metadata)
    metadata_values = [item for _, item in metadata]
    if not baseline_written:
        raise RuntimeError("No shared rho=0 baseline was found")
    if not summary_rows:
        raise RuntimeError("No positive-rho summary rows were found")

    scenario_order = {"all": 0, "horizontal": 1, "vertical": 2}
    summary_rows.sort(key=lambda row: (scenario_order[str(row["scenario"])], float(row["rho"])))
    trial_fields = list(trial_rows[0])
    summary_fields = list(summary_rows[0])
    write_csv(output / "rho_sweep_trial_errors.csv", trial_rows, trial_fields)
    write_csv(output / "rho_sweep_summary.csv", summary_rows, summary_fields)

    trials_per_rho = sorted({
        int(item.get("trials_per_rho_scenario", item.get("trials_per_rho", 0)))
        for item in metadata_values
        if int(item.get("trials_per_rho_scenario", item.get("trials_per_rho", 0))) > 0
    })
    reference = metadata_values[0]
    positive_rho_values = sorted({rho for _, rho, _ in all_trial_keys})
    converged_count = sum(
        str(item["converged"]).lower() == "true"
        for item in trial_rows
    )
    combined = {
        "format": "cesr-orbit-response-rho-sweep-chunked-v1",
        "chunks": [chunk.name for chunk in chunks],
        "chunk_metadata": metadata_values,
        "ring_id": reference["ring_id"],
        "lattice_path": reference["lattice_path"],
        "branch": reference["branch"],
        "rf_on": reference["rf_on"],
        "scibmad_version": reference["scibmad_version"],
        "control_count": reference["control_count"],
        "control_names": reference["control_names"],
        "observable_count": reference["observable_count"],
        "observable_labels": reference["observable_labels"],
        "detector_count": reference["detector_count"],
        "detector_count_per_plane": reference["detector_count_per_plane"],
        "input_csv": reference["input_csv"],
        "detector_response_csv": reference["detector_response_csv"],
        "closed_orbit_response_csv": reference["closed_orbit_response_csv"],
        "response_method": reference["response_method"],
        "response_step_rad": reference["response_step_rad"],
        "response_controls_per_batch": reference["response_controls_per_batch"],
        "response_pair_id": reference["response_pair_id"],
        "response_reltol": reference["response_reltol"],
        "response_abstol": reference["response_abstol"],
        "response_maxiter": reference["response_maxiter"],
        "total_unique_samples": len(trial_rows),
        "positive_rho_values": positive_rho_values,
        "trials_per_positive_rho_scenario": trials_per_rho[0] if len(trials_per_rho) == 1 else None,
        "converged_count": converged_count,
        "failed_count": len(trial_rows) - converged_count,
        "fallback_count": sum(int(item["fallback_count"]) for item in metadata_values),
        "fallback_count_is_chunk_sum": True,
        "exact_physics_seconds_sum": sum(float(item["exact_physics_seconds"]) for item in metadata_values),
        "response_evaluation_seconds_sum": sum(float(item["response_evaluation_seconds"]) for item in metadata_values),
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
