#!/usr/bin/env python3
"""Regression checks for the maintained paired P0--P3 benchmark."""

from __future__ import annotations

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results" / "paired_benchmark"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(source))


def main() -> None:
    estimates = {row["method"] + ":" + row["observable_view"]: row for row in rows(RESULTS / "p0_p2_offset_estimates.csv")}
    expected = {
        "P0_nominal_mixed_GTPSA:orbit_only": (20.0, 40.0),
        "P1_background_conditioned_mixed:orbit_only": (5.0, 10.0),
        "P0_nominal_mixed_GTPSA:orbit_phase_coupling_tune": (140.0, 165.0),
        "P1_background_conditioned_mixed:orbit_phase_coupling_tune": (1.0, 4.0),
        "P2a_two_local_dipole_kicks:orbit_only": (5.0, 10.0),
        "P2b_four_local_kicks:orbit_phase_coupling_tune": (1.0, 4.0),
    }
    for key, (lower_um, upper_um) in expected.items():
        error_um = 1e6 * float(estimates[key]["absolute_error_2d_m"])
        assert lower_um <= error_um <= upper_um, (key, error_um)
    p3 = rows(RESULTS / "p3_exact_history.csv")
    assert int(p3[0]["iteration"]) == 0
    assert 1e6 * float(p3[0]["absolute_error_2d_m"]) < 4.0
    assert 1e6 * float(p3[-1]["absolute_error_2d_m"]) < 1e-3
    assert float(p3[-1]["weighted_residual_rms"]) < 1e-10
    print("Paired P0--P3 benchmark validation passed.")


if __name__ == "__main__":
    main()
