#!/usr/bin/env python3
"""Validate all repaired-lattice exact-11 scan and inversion outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import tomllib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--scan-root", type=Path, default=HERE / "results" / "scans")
    result.add_argument("--aggregate-dir", type=Path, default=HERE / "results" / "aggregate")
    return result


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    args = parser().parse_args()
    scan_root = args.scan_root.expanduser().resolve()
    scan_dirs = [
        path for path in scan_root.iterdir()
        if path.is_dir() and (path / "scan_metadata.toml").exists()
    ]
    assert len(scan_dirs) == 76
    for scan_dir in scan_dirs:
        metadata = tomllib.loads((scan_dir / "scan_metadata.toml").read_text(encoding="utf-8"))
        scenarios = [
            line.strip()
            for line in (scan_dir / "scenario_labels.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        observations = np.load(scan_dir / "bpm_orbits.npy")
        target_orbits = np.load(scan_dir / "target_orbits.npy")
        expected = (len(scenarios), 11, 9, 5, 111, 2)
        assert observations.shape == expected
        assert target_orbits.shape == (len(scenarios), 11, 9, 5, 2)
        assert np.isfinite(observations).all()
        assert np.isfinite(target_orbits).all()
        triplets = rows(scan_dir / "quadratic_bump_triplet_inversion.csv")
        assert len(triplets) == 10
        assert sorted(int(row["information_rank"]) for row in triplets) == list(range(1, 11))
        assert sorted(int(row["precision_rank"]) for row in triplets) == list(range(1, 11))
        assert sorted(int(row["error_rank"]) for row in triplets) == list(range(1, 11))
        assert all(
            math.isfinite(float(row[key]))
            for row in triplets
            for key in (
                "position_error_um",
                "quadratic_center_information_logdet",
                "matched_nominal_information_gain_logdet",
                "matched_nominal_precision_improvement_worst_axis",
                "predicted_worst_axis_sigma_um",
                "center_equation_condition_number",
            )
        )
    aggregate_dir = args.aggregate_dir.expanduser().resolve()
    aggregate = rows(aggregate_dir / "exact11_validation_by_sextupole.csv")
    summary = json.loads(
        (aggregate_dir / "exact11_validation_summary.json").read_text(encoding="utf-8")
    )
    assert len(aggregate) == 76
    assert summary["target_count"] == 76
    linearity = rows(aggregate_dir / "k2_linearity_by_sextupole.csv")
    linearity_summary = json.loads(
        (aggregate_dir / "k2_linearity_summary.json").read_text(encoding="utf-8")
    )
    assert len(linearity) == 76
    assert linearity_summary["target_count"] == 76
    matched_path = aggregate_dir / "matched_nominal_response_triplets.csv"
    matched_summary_path = aggregate_dir / "matched_nominal_response_summary.json"
    if matched_path.exists() or matched_summary_path.exists():
        assert matched_path.exists() and matched_summary_path.exists()
        matched = rows(matched_path)
        matched_summary = json.loads(matched_summary_path.read_text(encoding="utf-8"))
        assert len(matched) == 760
        assert matched_summary["triplet_count"] == 760
    print(
        json.dumps(
            {
                "targets": len(scan_dirs),
                "triplets": 760,
                "response_matches_exact_information": summary[
                    "response_selection_matches_exact_information_count"
                ],
                "response_matches_exact_precision": summary[
                    "response_selection_matches_exact_precision_count"
                ],
                "response_matches_exact_error": summary[
                    "response_selection_matches_exact_error_count"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
