#!/usr/bin/env python3
"""Validate the maintained relative local-orbit predictor outputs."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results" / "local_orbit_predictors"
METHODS = ("command_only", "two_sided_transport", "global_map")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    summary = rows(RESULTS / "summary.csv")
    targets = rows(RESULTS / "per_target_summary.csv")
    predictions = rows(RESULTS / "per_prediction_errors.csv")
    neighbors = rows(RESULTS / "two_sided_neighbors.csv")
    cv = rows(RESULTS / "global_map_cv.csv")

    assert [row["method"] for row in summary] == list(METHODS)
    assert len(targets) == 3 * 76
    assert len(predictions) == 3 * 76 * 8 * 4
    assert len(neighbors) == 76
    assert sum(row["selected"] == "True" for row in cv) == 1
    assert len({row["target"] for row in targets}) == 76
    assert max(float(row["momentum_block_condition"]) for row in neighbors) < 5.0

    for expected in summary:
        selected = [row for row in predictions if row["method"] == expected["method"]]
        errors = np.asarray(
            [(float(row["error_x_um"]), float(row["error_y_um"])) for row in selected]
        )
        radial = np.linalg.norm(errors, axis=1)
        recomputed = {
            "x_rmse_um": np.sqrt(np.mean(errors[:, 0] ** 2)),
            "y_rmse_um": np.sqrt(np.mean(errors[:, 1] ** 2)),
            "rmse_2d_um": np.sqrt(np.mean(radial**2)),
            "median_2d_um": np.median(radial),
            "p90_2d_um": np.percentile(radial, 90),
            "p99_2d_um": np.percentile(radial, 99),
            "max_2d_um": np.max(radial),
        }
        for key, value in recomputed.items():
            assert np.isclose(value, float(expected[key]), rtol=1e-12, atol=1e-12), (
                expected["method"],
                key,
                value,
                expected[key],
            )

    print("validated local-orbit predictor results")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
