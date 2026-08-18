#!/usr/bin/env python3
"""Measure five-point K2 linearity for every completed exact-11 scan."""

from __future__ import annotations

import argparse
import csv
import json
import tomllib
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--scan-root", type=Path, default=HERE / "results" / "scans")
    result.add_argument(
        "--output-dir", type=Path, default=HERE / "results" / "aggregate"
    )
    return result


def distribution(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "p90": float(np.percentile(array, 90.0)),
        "maximum": float(np.max(array)),
    }


def relative_l2(numerator: np.ndarray, denominator: np.ndarray) -> float:
    denominator_norm = float(np.linalg.norm(denominator))
    if denominator_norm == 0.0:
        return 0.0
    return float(np.linalg.norm(numerator) / denominator_norm)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parser().parse_args()
    scan_root = args.scan_root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    rows: list[dict[str, Any]] = []
    for scan_dir in sorted(path for path in scan_root.iterdir() if path.is_dir()):
        metadata_path = scan_dir / "scan_metadata.toml"
        orbit_path = scan_dir / "bpm_orbits.npy"
        if not metadata_path.exists() or not orbit_path.exists():
            continue
        metadata = tomllib.loads(metadata_path.read_text(encoding="utf-8"))
        labels = [
            line.strip()
            for line in (scan_dir / "scenario_labels.txt").read_text(
                encoding="utf-8"
            ).splitlines()
            if line.strip()
        ]
        truth_index = labels.index("truth")
        observations = np.load(orbit_path, mmap_mode="r")[truth_index]
        levels = (
            np.asarray(metadata["k2_levels"], dtype=float)
            * float(metadata["k2_step_m3"])
        )
        if not np.allclose(levels, [-0.02, -0.01, 0.0, 0.01, 0.02]):
            raise ValueError(f"Unexpected K2 grid for {scan_dir.name}: {levels}")
        zero_index = int(np.flatnonzero(levels == 0.0)[0])
        signal = observations - observations[:, :, zero_index : zero_index + 1]
        flattened = np.moveaxis(observations, 2, 0).reshape(5, -1)
        signal_flattened = np.moveaxis(signal, 2, 0).reshape(5, -1)

        linear_design = np.column_stack([np.ones(5), levels])
        quadratic_design = np.column_stack([np.ones(5), levels, levels**2])
        linear_fit = linear_design @ np.linalg.lstsq(
            linear_design, flattened, rcond=None
        )[0]
        quadratic_fit = quadratic_design @ np.linalg.lstsq(
            quadratic_design, flattened, rcond=None
        )[0]
        inner_slope = (flattened[3] - flattened[1]) / (levels[3] - levels[1])
        outer_slope = (flattened[4] - flattened[0]) / (levels[4] - levels[0])
        five_point_slope = (levels @ flattened) / float(levels @ levels)

        rows.append(
            {
                "sextupole": metadata["target_sextupole"],
                "linear_fit_relative_l2_residual": relative_l2(
                    flattened - linear_fit, signal_flattened
                ),
                "quadratic_fit_relative_l2_residual": relative_l2(
                    flattened - quadratic_fit, signal_flattened
                ),
                "inner_outer_slope_relative_l2_difference": relative_l2(
                    outer_slope - inner_slope, five_point_slope
                ),
                "five_point_inner_slope_relative_l2_difference": relative_l2(
                    five_point_slope - inner_slope, five_point_slope
                ),
                "five_point_outer_slope_relative_l2_difference": relative_l2(
                    five_point_slope - outer_slope, five_point_slope
                ),
                "k2_signal_rms_m": float(np.sqrt(np.mean(signal_flattened**2))),
            }
        )

    if len(rows) != 76:
        raise RuntimeError(f"Expected 76 completed targets, found {len(rows)}")
    rows.sort(key=lambda row: str(row["sextupole"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "k2_linearity_by_sextupole.csv", rows)
    metric_names = [
        "linear_fit_relative_l2_residual",
        "quadratic_fit_relative_l2_residual",
        "inner_outer_slope_relative_l2_difference",
        "five_point_inner_slope_relative_l2_difference",
        "five_point_outer_slope_relative_l2_difference",
        "k2_signal_rms_m",
    ]
    summary = {
        "format": "cesr-repaired-lattice-exact-11-k2-linearity-v1",
        "target_count": len(rows),
        "k2_grid_m3": [-0.02, -0.01, 0.0, 0.01, 0.02],
        "metrics": {
            name: distribution([float(row[name]) for row in rows])
            for name in metric_names
        },
        "interpretation": (
            "The inner-versus-outer symmetric slope difference isolates odd "
            "higher-order K2 contamination; the linear-fit residual also includes "
            "even K2 curvature that cancels from a symmetric slope."
        ),
    }
    (output_dir / "k2_linearity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
