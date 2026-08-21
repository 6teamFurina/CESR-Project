#!/usr/bin/env python3
"""Regression checks for the unknown-background linear/nonlinear P2b test."""

from __future__ import annotations

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE / "results" / "unknown_background_benchmark"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(source))


def final_error_um(path: Path) -> float:
    return 1e6 * float(rows(path)[-1]["absolute_error_2d_m"])


def main() -> None:
    unknown_linear = rows(ROOT / "linear" / "p0_p2_offset_estimates.csv")
    p1 = next(row for row in unknown_linear if row["method"] == "P1_nominal_bump_conditioned_mixed" and row["observable_view"] == "orbit_phase_coupling_tune")
    p2 = next(row for row in unknown_linear if row["method"] == "P2b_four_local_kicks")
    p1_error = 1e6 * float(p1["absolute_error_2d_m"])
    p2_error = 1e6 * float(p2["absolute_error_2d_m"])
    nonlinear_error = final_error_um(ROOT / "nonlinear_p2b" / "nonlinear_p2b_history.csv")
    closure_linear = rows(ROOT / "closure_linear" / "p0_p2_offset_estimates.csv")
    closure_p2 = next(row for row in closure_linear if row["method"] == "P2b_four_local_kicks")
    closure_linear_error = 1e6 * float(closure_p2["absolute_error_2d_m"])
    closure_nonlinear_error = final_error_um(ROOT / "closure_nonlinear" / "nonlinear_p2b_history.csv")
    assert 140.0 < p1_error < 180.0
    assert 140.0 < p2_error < 180.0
    assert 140.0 < nonlinear_error < 180.0
    assert abs(nonlinear_error - p2_error) < 2.0
    assert 1.0 < closure_linear_error < 4.0
    assert closure_nonlinear_error < 1.0
    calibration = rows(ROOT / "nonlinear_p2b" / "reconstructed_sources.csv")
    assert max(float(row["calibration_source_rms"]) for row in calibration) < 1e-8
    print("Unknown-background nonlinear P2b validation passed.")


if __name__ == "__main__":
    main()
