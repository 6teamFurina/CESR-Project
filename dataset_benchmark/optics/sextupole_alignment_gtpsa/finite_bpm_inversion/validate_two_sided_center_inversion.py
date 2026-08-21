#!/usr/bin/env python3
"""Validate maintained two-sided-BPM center-inversion outputs."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results" / "two_sided_center_inversion"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    summary = rows(RESULTS / "summary.csv")
    realizations = rows(RESULTS / "per_realization_fits.csv")
    targets = rows(RESULTS / "per_target_summary.csv")
    oracle = rows(RESULTS / "oracle_difference_summary.csv")
    estimates = np.load(RESULTS / "relative_center_estimates.npy")

    assert len(summary) == 1
    assert summary[0]["method"] == "two_sided_transport"
    assert len(realizations) == 76 * 8
    assert len(targets) == 76
    assert estimates.shape == (76, 8, 2)
    assert np.all(np.isfinite(estimates))
    assert len({row["target"] for row in realizations}) == 76

    errors = np.asarray(
        [(float(row["error_x_um"]), float(row["error_y_um"])) for row in realizations]
    )
    radial = np.linalg.norm(errors, axis=1)
    expected = summary[0]
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
        assert np.isclose(value, float(expected[key]), rtol=1e-12, atol=1e-12)

    oracle_difference = np.asarray(
        [float(row["error_vector_difference_from_oracle_um"]) for row in realizations]
    )
    oracle_expected = oracle[0]
    oracle_recomputed = {
        "rms_2d_um": np.sqrt(np.mean(oracle_difference**2)),
        "median_2d_um": np.median(oracle_difference),
        "p90_2d_um": np.percentile(oracle_difference, 90),
        "max_2d_um": np.max(oracle_difference),
    }
    for key, value in oracle_recomputed.items():
        assert np.isclose(value, float(oracle_expected[key]), rtol=1e-12, atol=1e-12)

    print("validated two-sided-BPM center inversion")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
