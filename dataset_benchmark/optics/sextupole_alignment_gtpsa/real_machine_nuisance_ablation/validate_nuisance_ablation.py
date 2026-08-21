#!/usr/bin/env python3
"""Independently validate the maintained nuisance-ablation summaries."""

from __future__ import annotations

import csv
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PHYSICAL_ROOT = HERE / "results" / "physical_scans"
RESULTS = HERE / "results" / "analysis"
EXPECTED_CASES = (
    "baseline",
    "bpm_gain",
    "corrector_gain",
    "k2_calibration",
    "quadrupole_strength",
    "quadrupole_roll",
    "quadrupole_misalignment",
    "time_drift",
    "bpm_noise",
)
PHYSICAL_CASES = (
    "baseline",
    "corrector_gain",
    "k2_calibration",
    "quadrupole_strength",
    "quadrupole_roll",
    "quadrupole_misalignment",
    "time_drift",
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    summary = rows(RESULTS / "summary.csv")
    realizations = rows(RESULTS / "per_realization_fits.csv")
    targets = rows(RESULTS / "per_target_summary.csv")
    assert tuple(row["case"] for row in summary) == EXPECTED_CASES

    with (PHYSICAL_ROOT / "baseline" / "scan_metadata.toml").open("rb") as stream:
        baseline = tomllib.load(stream)
    nt = int(baseline["target_count"])
    nr = int(baseline["realization_count_per_target"])
    assert len(realizations) == len(EXPECTED_CASES) * nt * nr
    assert len(targets) == len(EXPECTED_CASES) * nt

    baseline_truth = np.load(PHYSICAL_ROOT / "baseline" / "target_truth.npy")
    baseline_sext = np.load(PHYSICAL_ROOT / "baseline" / "latent_sextupole_offsets.npy")
    for case_name in PHYSICAL_CASES:
        source = PHYSICAL_ROOT / case_name
        with (source / "scan_metadata.toml").open("rb") as stream:
            metadata = tomllib.load(stream)
        assert metadata["nuisance_case"] == case_name
        assert "latest_cesr_scibmad_repaired.jl" in str(metadata["lattice"])
        assert np.array_equal(np.load(source / "target_truth.npy"), baseline_truth)
        assert np.array_equal(
            np.load(source / "latent_sextupole_offsets.npy"), baseline_sext
        )

    baseline_errors: np.ndarray | None = None
    for expected in summary:
        case_name = expected["case"]
        selected = [row for row in realizations if row["case"] == case_name]
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
            assert np.isclose(value, float(expected[key]), rtol=1e-12, atol=1e-12)
        if case_name == "baseline":
            baseline_errors = errors
            expected_incremental = 0.0
        else:
            assert baseline_errors is not None
            expected_incremental = np.sqrt(
                np.mean(np.sum((errors - baseline_errors) ** 2, axis=1))
            )
        assert np.isclose(
            expected_incremental,
            float(expected["incremental_error_vector_rms_um"]),
            rtol=1e-12,
            atol=1e-12,
        )
        estimates = np.load(RESULTS / f"{case_name}_relative_center_estimates.npy")
        predicted = np.load(RESULTS / f"{case_name}_predicted_local_orbits.npy")
        assert estimates.shape == (nt, nr, 2)
        assert predicted.shape == (nt, nr, 5, 2)
        assert np.all(np.isfinite(estimates))
        assert np.all(np.isfinite(predicted))

    print("validated real-machine nuisance ablation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
