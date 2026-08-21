#!/usr/bin/env python3
"""Validate the two maintained end-to-end targeted-inversion smoke cases."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
EXPECTED_STATES = 15
EXPECTED_OBSERVATIONS_PER_STATE = 1190
EXPECTED_ESTIMATES = 8


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(source))


def validate_case(name: str) -> dict[str, object]:
    root = RESULTS / name
    states = rows(root / "scan_states.csv")
    observations = rows(root / "scan_observations.csv")
    estimates = rows(root / "inversion" / "offset_estimates.csv")
    if len(states) != EXPECTED_STATES:
        raise ValueError(f"{name}: expected {EXPECTED_STATES} states, found {len(states)}")
    if len(observations) != EXPECTED_STATES * EXPECTED_OBSERVATIONS_PER_STATE:
        raise ValueError(
            f"{name}: expected {EXPECTED_STATES * EXPECTED_OBSERVATIONS_PER_STATE} "
            f"observations, found {len(observations)}"
        )
    if len(estimates) != EXPECTED_ESTIMATES:
        raise ValueError(
            f"{name}: expected {EXPECTED_ESTIMATES} estimates, found {len(estimates)}"
        )
    for row in observations:
        if not math.isfinite(float(row["observable_readback"])):
            raise ValueError(f"{name}: non-finite observation: {row}")
    for row in estimates:
        for field in (
            "estimated_x_offset_m",
            "estimated_y_offset_m",
            "error_x_m",
            "error_y_m",
            "absolute_error_2d_m",
            "whitened_residual_rms",
        ):
            if not math.isfinite(float(row[field])):
                raise ValueError(f"{name}: non-finite {field}: {row}")
    orbit_modes = next(
        row for row in estimates
        if row["observable_ablation"] == "orbit_only"
        and row["input_representation"] == "modes2"
    )
    return {
        "case": name,
        "states": len(states),
        "observation_rows": len(observations),
        "estimate_rows": len(estimates),
        "truth_um": {
            "x": 1e6 * float(orbit_modes["true_x_offset_m"]),
            "y": 1e6 * float(orbit_modes["true_y_offset_m"]),
        },
        "orbit_modes_estimate_um": {
            "x": 1e6 * float(orbit_modes["estimated_x_offset_m"]),
            "y": 1e6 * float(orbit_modes["estimated_y_offset_m"]),
        },
        "orbit_modes_error_2d_um": 1e6 * float(orbit_modes["absolute_error_2d_m"]),
    }


def main() -> None:
    results = [validate_case(name) for name in ("smoke_exact", "smoke_background")]
    output = RESULTS / "smoke_validation_summary.json"
    output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    for result in results:
        print(
            f"{result['case']}: {result['states']} states, "
            f"{result['observation_rows']} finite observations, "
            f"orbit-mode 2D error {result['orbit_modes_error_2d_um']:.3f} um"
        )
    print(f"Summary: {output}")


if __name__ == "__main__":
    main()

