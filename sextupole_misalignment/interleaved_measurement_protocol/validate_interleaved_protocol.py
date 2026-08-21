#!/usr/bin/env python3
"""Independently validate the interleaved-protocol result tables."""

from __future__ import annotations

import csv
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PHYSICAL_ROOT = HERE.parent / "real_machine_nuisance_ablation" / "results" / "physical_scans"
RESULTS = HERE / "results"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    summary = rows(RESULTS / "summary.csv")
    realizations = rows(RESULTS / "per_realization_fits.csv")
    validation = rows(RESULTS / "drift_response_validation.csv")
    optimizer_validation = rows(RESULTS / "clean_optimizer_validation.csv")
    assert len(validation) == 1
    assert len(optimizer_validation) == 1
    assert summary
    case_names = tuple(row["case"] for row in summary)
    assert len(case_names) == len(set(case_names))
    expected_nuisances = {"bpm_noise", "random_walk_drift", "combined"}
    assert {row["nuisance"] for row in summary} == expected_nuisances

    with (PHYSICAL_ROOT / "baseline" / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    assert "latest_cesr_scibmad_repaired.jl" in str(metadata["lattice"])
    assert np.max(np.abs(np.load(PHYSICAL_ROOT / "baseline" / "latent_quadrupole_offsets.npy"))) == 0.0
    assert np.any(np.load(PHYSICAL_ROOT / "baseline" / "latent_sextupole_offsets.npy") != 0.0)

    archives = np.load(RESULTS / "relative_center_estimates.npz")
    assert set(archives.files) == {"clean", *case_names}
    fit_count = int(summary[0]["fit_count"])
    assert len(realizations) == len(summary) * fit_count

    for expected in summary:
        selected = [row for row in realizations if row["case"] == expected["case"]]
        assert len(selected) == fit_count
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
        boundary_fraction = np.mean(
            np.any(np.abs(archives[expected["case"]]) >= 1.49e-3, axis=-1)
        )
        assert np.isclose(
            boundary_fraction,
            float(expected["fit_boundary_fraction"]),
            rtol=1e-12,
            atol=1e-12,
        )
        repeats = int(expected["repeats_per_nonzero_point"])
        if expected["protocol"] == "blocked":
            assert int(expected["acquisitions_per_scan"]) == 15 * repeats
        else:
            assert int(expected["acquisitions_per_scan"]) == 5 * (4 * repeats + 1)
        estimate = archives[expected["case"]]
        assert estimate.shape == (
            int(metadata["target_count"]),
            int(metadata["realization_count_per_target"]),
            2,
        )
        assert np.all(np.isfinite(estimate))

    clean = archives["clean"]
    assert np.all(np.isfinite(clean))
    print("validated interleaved and repeated acquisition protocol")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
