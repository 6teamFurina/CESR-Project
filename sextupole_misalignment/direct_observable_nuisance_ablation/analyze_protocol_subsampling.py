#!/usr/bin/env python3
"""Subsample K2/bump grids and refit with a shared thin-sextupole source model."""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def k2_slope(values: np.ndarray, delta: np.ndarray) -> np.ndarray:
    centered = delta - delta.mean()
    return np.einsum("bkc,k->bc", values, centered) / np.dot(centered, centered)


def source_matrix(xy: np.ndarray, center: np.ndarray) -> np.ndarray:
    x = xy[:, 0] - center[0]
    y = xy[:, 1] - center[1]
    return np.column_stack((0.5 * (x * x - y * y), x * y))


def fit_center(slopes: np.ndarray, xy: np.ndarray) -> np.ndarray:
    scale = np.sqrt(np.mean(slopes * slopes, axis=0))
    positive = scale[scale > 0]
    floor = np.median(positive) * 1e-8 if positive.size else 1.0
    normalized = slopes / np.maximum(scale, floor)

    def residual(center: np.ndarray) -> np.ndarray:
        source = source_matrix(xy, center)
        propagation = np.linalg.lstsq(source, normalized, rcond=1e-12)[0]
        return (source @ propagation - normalized).ravel()

    starts = [
        np.zeros(2),
        np.mean(xy, axis=0),
        np.array([xy[:, 0].min(), 0.0]),
        np.array([xy[:, 0].max(), 0.0]),
        np.array([0.0, xy[:, 1].min()]),
        np.array([0.0, xy[:, 1].max()]),
    ]
    solutions = [
        least_squares(
            residual, start, bounds=(-1.5e-3, 1.5e-3),
            xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=500,
        ) for start in starts
    ]
    best = min(solutions, key=lambda result: np.dot(result.fun, result.fun))
    return best.x


def summarize(estimates: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    error = np.linalg.norm(estimates - truth, axis=1) * 1e6
    return {
        "rmse_2d_um": float(np.sqrt(np.mean(error**2))),
        "median_2d_um": float(np.median(error)),
        "p90_2d_um": float(np.percentile(error, 90)),
        "max_2d_um": float(np.max(error)),
    }


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir", type=Path,
        default=here / "results" / "sex_09aw_paired_pilot",
    )
    args = parser.parse_args()
    source = args.input_dir.resolve()
    with (source / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    orbit = np.load(source / "bpm_orbits.npy")
    target_orbit = np.load(source / "target_orbits.npy")
    truth = np.load(source / "target_truth.npy")
    orbit = orbit.reshape(*orbit.shape[:3], -1)
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    full_delta = levels * float(metadata["k2_step_m3"])
    nominal = int(np.flatnonzero(levels == 0)[0])
    local_xy = target_orbit[:, :, nominal, :]

    # Julia grid ordering: x outer, y inner.
    bump_sets = {
        "9_grid": np.arange(9),
        "5_cross": np.array([1, 3, 4, 5, 7]),
        "5_corners_center": np.array([0, 2, 4, 6, 8]),
    }
    k2_sets = {
        "5_full": np.arange(5),
        "3_inner": np.array([1, 2, 3]),       # -K, 0, +K
        "3_outer": np.array([0, 2, 4]),       # -2K, 0, +2K
    }
    rows = []
    for bump_name, bump_indices in bump_sets.items():
        for k2_name, k2_indices in k2_sets.items():
            estimates = np.zeros_like(truth)
            for realization in range(orbit.shape[0]):
                selected = orbit[realization][bump_indices][:, k2_indices, :]
                slopes = k2_slope(selected, full_delta[k2_indices])
                estimates[realization] = fit_center(
                    slopes, local_xy[realization, bump_indices],
                )
            row = {
                "bump_protocol": bump_name,
                "bump_count": len(bump_indices),
                "k2_protocol": k2_name,
                "k2_count": len(k2_indices),
                **summarize(estimates, truth),
            }
            rows.append(row)
            print(row)
    output = source / "orbit_protocol_subsampling.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
