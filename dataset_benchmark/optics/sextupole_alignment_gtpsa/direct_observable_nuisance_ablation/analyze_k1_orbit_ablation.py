#!/usr/bin/env python3
"""Fit the paired orbit-only K1 ablation with a common sextupole center."""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from analyze_protocol_subsampling import k2_slope, source_matrix


def fit_conditions(slopes: np.ndarray, xy: np.ndarray, selected: np.ndarray) -> np.ndarray:
    normalized = []
    for condition in selected:
        values = slopes[condition]
        scale = np.sqrt(np.mean(values * values, axis=0))
        positive = scale[scale > 0]
        floor = np.median(positive) * 1e-8 if positive.size else 1.0
        normalized.append(values / np.maximum(scale, floor))

    def residual(center: np.ndarray) -> np.ndarray:
        pieces = []
        for values, condition in zip(normalized, selected):
            source = source_matrix(xy[condition], center)
            propagation = np.linalg.lstsq(source, values, rcond=1e-12)[0]
            pieces.append((source @ propagation - values).ravel())
        return np.concatenate(pieces)

    starts = [np.zeros(2), *[xy[0, i] for i in range(xy.shape[1])]]
    solutions = [
        least_squares(
            residual, start, bounds=(-1.5e-3, 1.5e-3),
            xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=500,
        ) for start in starts
    ]
    return min(solutions, key=lambda result: np.dot(result.fun, result.fun)).x


def metrics(estimates: np.ndarray, truth: np.ndarray):
    error = np.linalg.norm(estimates - truth, axis=1) * 1e6
    return (
        float(np.sqrt(np.mean(error**2))), float(np.median(error)),
        float(np.percentile(error, 90)), float(np.max(error)), error,
    )


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir", type=Path,
        default=here / "results" / "sex_09aw_k1_orbit_ablation",
    )
    args = parser.parse_args()
    source = args.input_dir.resolve()
    with (source / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    orbit = np.load(source / "bpm_orbits.npy")
    target_orbit = np.load(source / "target_orbits.npy")
    truth = np.load(source / "target_truth.npy")
    names = (source / "condition_names.txt").read_text().splitlines()
    orbit = orbit.reshape(*orbit.shape[:4], -1)
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta = levels * float(metadata["k2_step_m3"])
    nominal_k2 = int(np.flatnonzero(levels == 0)[0])
    nr, nc = orbit.shape[:2]
    slopes = np.empty((nr, nc, orbit.shape[2], orbit.shape[-1]))
    for realization in range(nr):
        for condition in range(nc):
            slopes[realization, condition] = k2_slope(orbit[realization, condition], delta)
    local_xy = target_orbit[:, :, :, nominal_k2, :]

    protocols = {
        "no_k1_scan": np.array([0]),
        "all_7_k1_conditions": np.arange(nc),
    }
    estimates = {}
    rows = []
    for label, selected in protocols.items():
        estimate = np.zeros_like(truth)
        for realization in range(nr):
            estimate[realization] = fit_conditions(
                slopes[realization], local_xy[realization], selected,
            )
        estimates[label] = estimate
        rmse, median, p90, maximum, error = metrics(estimate, truth)
        rows.append({
            "protocol": label,
            "condition_count": len(selected),
            "rmse_2d_um": rmse,
            "median_2d_um": median,
            "p90_2d_um": p90,
            "max_2d_um": maximum,
            "paired_better_than_no_k1_count": 0 if label == "no_k1_scan" else int(np.sum(
                error < metrics(estimates["no_k1_scan"], truth)[4]
            )),
        })
    with (source / "k1_orbit_ablation_summary.csv").open(
        "w", newline="", encoding="utf-8",
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("conditions:", names)
    for row in rows:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
