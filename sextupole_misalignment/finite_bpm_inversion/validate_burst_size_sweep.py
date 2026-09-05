#!/usr/bin/env python3
"""Validate the paired full-error burst-size sweep."""

from __future__ import annotations

import csv
import inspect
import json
import tomllib
from pathlib import Path

import numpy as np

import analyze_burst_size_sweep as sweep
import analyze_state_space_bpm_gtpsa_inverse as maintained


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results" / "burst_size_sweep"
SOURCE = maintained.SCAN_ROOT / maintained.DEFAULT_CASE
LATENT = maintained.SCAN_ROOT / "paired_latents"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    metadata = json.loads((RESULTS / "analysis_metadata.json").read_text())
    burst_sizes = np.load(RESULTS / "burst_sizes.npy")
    expected_bursts = np.asarray(metadata["burst_sizes"], dtype=int)
    require(np.array_equal(burst_sizes, expected_bursts), "Burst inventory changed")
    nb = len(burst_sizes)
    na = int(metadata["stochastic_augmentations"])
    nm = int(metadata["machine_count"])
    nt = int(metadata["target_count"])
    vector_shape = (nb, na, nm, nt, 2)
    scalar_shape = (nb, na, nm, nt)
    arrays = {}
    for name in (
        "unfiltered_relative_center_estimates",
        "unfiltered_absolute_offset_estimates",
        "filtered_relative_center_estimates",
        "filtered_absolute_offset_estimates",
    ):
        values = np.load(RESULTS / f"{name}.npy")
        require(values.shape == vector_shape, f"Unexpected {name} shape")
        require(np.all(np.isfinite(values)), f"Non-finite {name}")
        arrays[name] = values
    for name in (
        "unfiltered_bpm_state_error_rms_m",
        "filtered_bpm_state_error_rms_m",
    ):
        values = np.load(RESULTS / f"{name}.npy")
        require(values.shape == scalar_shape, f"Unexpected {name} shape")
        require(np.all(np.isfinite(values)), f"Non-finite {name}")

    protocol = rows(RESULTS / "protocol_summary.csv")
    require(len(protocol) == nb, "Protocol row count changed")
    visits = []
    for row, burst in zip(protocol, burst_sizes):
        b = int(burst)
        require(int(row["burst_size"]) == b, "Protocol burst order changed")
        require(int(row["total_signal_turns"]) == 24_576, "Signal turns changed")
        require(
            int(row["total_protocol_acquisitions"]) == 24_732,
            "Protocol acquisition count changed",
        )
        require(int(row["reference_cycle_count"]) == 13, "Reference cycles changed")
        require(int(row["reference_event_count"]) == 156, "Reference count changed")
        require(
            int(row["idealized_signal_state_visits"]) == 24_576 // b,
            "Signal visit scaling changed",
        )
        visits.append(int(row["physical_state_visits_excluding_calibration"]))
    require(all(left > right for left, right in zip(visits[:-1], visits[1:])),
            "State visits do not fall monotonically")

    with (SOURCE / "scan_metadata.toml").open("rb") as stream:
        source_metadata = tomllib.load(stream)
    require(source_metadata["baseline_orbit_correction_applied"] is True,
            "Full-error source lacks orbit correction")
    require(source_metadata["baseline_gtpsa_response_model"] == "nominal",
            "Full-error source changed GTPSA semantics")
    exact_reference = np.asarray(
        np.load(SOURCE / "reference_target_orbits.npy", mmap_mode="r")
    )[:nm, :nt]
    latent_offsets = np.asarray(
        np.load(LATENT / "sextupole_offsets.npy", mmap_mode="r")
    )[:nm, :nt]
    relative_truth = latent_offsets - exact_reference
    summary = rows(RESULTS / "burst_summary.csv")
    require(len(summary) == nb, "Summary row count changed")
    maximum_summary_difference = 0.0
    for index, row in enumerate(summary):
        error_sets = {
            "unfiltered_relative": arrays["unfiltered_relative_center_estimates"][index]
            - relative_truth[None],
            "unfiltered_absolute": arrays["unfiltered_absolute_offset_estimates"][index]
            - latent_offsets[None],
            "filtered_relative": arrays["filtered_relative_center_estimates"][index]
            - relative_truth[None],
            "filtered_absolute": arrays["filtered_absolute_offset_estimates"][index]
            - latent_offsets[None],
        }
        for prefix, errors in error_sets.items():
            calculated = sweep.error_metrics(errors)
            for key, value in calculated.items():
                saved = float(row[f"{prefix}_{key}"])
                maximum_summary_difference = max(
                    maximum_summary_difference, abs(saved - value)
                )
                require(
                    np.isclose(saved, value, rtol=2.0e-11, atol=2.0e-11),
                    f"Saved {prefix}_{key} changed for burst {burst_sizes[index]}",
                )

    production_difference = max(
        metadata["production_b1_max_abs_differences_m"].values(), default=0.0
    )
    require(production_difference <= 2.0e-15, "B=1 no longer reproduces production")
    source_text = inspect.getsource(sweep.main)
    persist_index = source_text.index("Persist every machine-facing estimate")
    truth_index = source_text.index("Evaluation-only boundary")
    require(persist_index < truth_index, "Truth boundary moved ahead of persistence")
    require((RESULTS / "burst_size_tradeoff.png").stat().st_size > 10_000,
            "Tradeoff figure is missing or empty")

    report = {
        "format": "cesr-full-error-burst-size-sweep-validation-v1",
        "status": "PASS",
        "burst_sizes": burst_sizes.tolist(),
        "array_shape": list(vector_shape),
        "maximum_summary_recalculation_difference": maximum_summary_difference,
        "production_b1_max_abs_difference_m": production_difference,
        "truth_persistence_order_check": "PASS",
        "protocol_invariants_check": "PASS",
    }
    (RESULTS / "VALIDATION.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
