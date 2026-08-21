#!/usr/bin/env python3
"""Invert measured BPM K2 slopes with direct latest-lattice GTPSA maps."""

from __future__ import annotations

import argparse
import csv
import math
import tomllib
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from analyze_inverses import (
    channel_scales,
    polynomial_k2_slope,
    read_rows,
    summarize_vectors,
    write_rows,
)


HERE = Path(__file__).resolve().parent


def evaluate_k2_slope_map(
    derivatives: np.ndarray,
    monomials: list[tuple[int, int]],
    effective_offset_m: np.ndarray,
    retained_offset_order: int,
) -> np.ndarray:
    x, y = effective_offset_m
    result = np.zeros(derivatives.shape[1])
    for index, (x_power, y_power) in enumerate(monomials):
        if x_power + y_power > retained_offset_order:
            continue
        result += (
            derivatives[index]
            * x**x_power
            * y**y_power
            / (math.factorial(x_power) * math.factorial(y_power))
        )
    return result


def fit_gtpsa_center(
    measured_slopes: np.ndarray,
    local_xy_m: np.ndarray,
    derivatives: np.ndarray,
    monomials: list[tuple[int, int]],
    retained_offset_order: int,
    zero_bump: int,
    initial_m: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    observed = measured_slopes - measured_slopes[zero_bump : zero_bump + 1]
    scale = channel_scales(observed)

    def residual(center_m: np.ndarray) -> np.ndarray:
        reference = evaluate_k2_slope_map(
            derivatives,
            monomials,
            center_m - local_xy_m[zero_bump],
            retained_offset_order,
        )
        predictions = np.asarray(
            [
                evaluate_k2_slope_map(
                    derivatives,
                    monomials,
                    center_m - local_xy,
                    retained_offset_order,
                )
                - reference
                for local_xy in local_xy_m
            ]
        )
        return ((predictions - observed) / scale[None, :]).ravel()

    starts = [
        initial_m,
        np.zeros(2),
        np.array([5.0e-4, 0.0]),
        np.array([-5.0e-4, 0.0]),
        np.array([0.0, 5.0e-4]),
        np.array([0.0, -5.0e-4]),
    ]
    solutions = [
        least_squares(
            residual,
            np.clip(start, -1.5e-3, 1.5e-3),
            bounds=(-1.5e-3, 1.5e-3),
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
            max_nfev=500,
        )
        for start in starts
    ]
    best = min(solutions, key=lambda item: float(np.dot(item.fun, item.fun)))
    jacobian_singular = np.linalg.svd(best.jac, compute_uv=False)
    condition = (
        float(jacobian_singular[0] / jacobian_singular[-1])
        if jacobian_singular[-1] > 0
        else np.inf
    )
    residual_rms = float(np.linalg.norm(best.fun) / np.sqrt(best.fun.size))
    return best.x, condition, residual_rms


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-dir", type=Path, default=HERE / "results" / "exact_scans")
    parser.add_argument("--local-analysis-dir", type=Path, default=HERE / "results" / "analysis")
    parser.add_argument("--map-dir", type=Path, default=HERE / "results" / "gtpsa_maps")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results" / "gtpsa_map_analysis")
    args = parser.parse_args()
    scan_dir = args.scan_dir.resolve()
    local_dir = args.local_analysis_dir.resolve()
    map_dir = args.map_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with (scan_dir / "scan_metadata.toml").open("rb") as stream:
        scan_metadata = tomllib.load(stream)
    with (map_dir / "map_metadata.toml").open("rb") as stream:
        map_metadata = tomllib.load(stream)
    scan_targets = (scan_dir / "target_names.txt").read_text(encoding="utf-8").splitlines()
    map_targets = (map_dir / "target_names.txt").read_text(encoding="utf-8").splitlines()
    scan_bpms = (scan_dir / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    map_bpms = (map_dir / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    if scan_targets != map_targets or scan_bpms != map_bpms:
        raise ValueError("Scan and GTPSA map inventories differ")
    target_names = scan_targets
    levels = np.asarray(scan_metadata["k2_levels"], dtype=float)
    delta_k2 = levels * float(scan_metadata["k2_step_m3"])
    nominal_k2 = int(np.flatnonzero(levels == 0.0)[0])
    bump_rows = read_rows(scan_dir / "bump_points.csv")
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    zero_bump = int(np.flatnonzero(np.all(bump_commands == 0.0, axis=1))[0])
    monomial_rows = read_rows(map_dir / "offset_monomials.csv")
    monomials = [
        (int(row["x_offset_power"]), int(row["y_offset_power"]))
        for row in monomial_rows
    ]
    derivatives = np.load(map_dir / "k2_offset_derivatives.npy")
    predicted_local = np.load(local_dir / "predicted_relative_local_orbits.npy")
    bpm_orbits = np.load(scan_dir / "bpm_orbits.npy", mmap_mode="r")
    nt, nr, nb, nk, nd, planes = bpm_orbits.shape
    if derivatives.shape != (nt, len(monomials), 2 * nd):
        raise ValueError(f"Unexpected GTPSA derivative shape: {derivatives.shape}")
    bpm_flat = np.asarray(bpm_orbits, dtype=float).reshape(nt, nr, nb, nk, 2 * nd)

    estimates = {
        f"gtpsa_k2_offset_order{order}_predicted_local": np.zeros((nt, nr, 2))
        for order in (1, 2, 3)
    }
    diagnostics: list[dict[str, object]] = []
    for target in range(nt):
        for realization in range(nr):
            slopes = polynomial_k2_slope(
                bpm_flat[target, realization], delta_k2, min(4, nk - 1)
            )
            local_xy = predicted_local[target, realization, :, nominal_k2, :]
            initial = np.zeros(2)
            for order in (1, 2, 3):
                method = f"gtpsa_k2_offset_order{order}_predicted_local"
                result = fit_gtpsa_center(
                    slopes,
                    local_xy,
                    derivatives[target],
                    monomials,
                    order,
                    zero_bump,
                    initial,
                )
                estimates[method][target, realization] = result[0]
                initial = result[0]
                diagnostics.append(
                    {
                        "target": target_names[target],
                        "target_index": target + 1,
                        "realization": realization + 1,
                        "method": method,
                        "jacobian_condition": result[1],
                        "scaled_residual_rms": result[2],
                    }
                )
        print(f"GTPSA inverse {target + 1}/{nt}: {target_names[target]}")
    write_rows(output / "fit_diagnostics.csv", diagnostics)
    for method, estimate in estimates.items():
        np.save(output / f"{method}_relative_center_estimates.npy", estimate)

    # Evaluation-only truth.
    target_orbits = np.load(scan_dir / "target_orbits.npy", mmap_mode="r")
    target_truth = np.load(scan_dir / "target_truth.npy")
    nominal_centers = np.load(scan_dir / "nominal_target_centers.npy")
    exact_reference = np.asarray(
        target_orbits[:, :, zero_bump, nominal_k2, :], dtype=float
    )
    relative_truth = nominal_centers[:, None, :] + target_truth - exact_reference
    summary_rows = []
    per_case_rows = []
    for method, estimate in estimates.items():
        errors = estimate - relative_truth
        summary_rows.append({"method": method, **summarize_vectors(errors)})
        for target, name in enumerate(target_names):
            for realization in range(nr):
                vector_um = errors[target, realization] * 1e6
                per_case_rows.append(
                    {
                        "method": method,
                        "target": name,
                        "target_index": target + 1,
                        "realization": realization + 1,
                        "truth_x_um": relative_truth[target, realization, 0] * 1e6,
                        "truth_y_um": relative_truth[target, realization, 1] * 1e6,
                        "estimate_x_um": estimate[target, realization, 0] * 1e6,
                        "estimate_y_um": estimate[target, realization, 1] * 1e6,
                        "error_2d_um": np.linalg.norm(vector_um),
                    }
                )
    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "per_case_estimates.csv", per_case_rows)
    lines = "\n".join(
        f"| {row['method']} | {row['rmse_2d_um']:.6f} | {row['median_2d_um']:.6f} | "
        f"{row['p90_2d_um']:.6f} | {row['max_2d_um']:.6f} |"
        for row in summary_rows
    )
    report = f"""# Direct SciBmad/GTPSA K2--offset map inverse

The map is generated directly by SciBmad/GTPSA on the nominal validated latest
lattice.  Measured quartic-in-K2 BPM slopes and BPM-predicted target-local
coordinates come from the common sextupole-misalignment-only exact scan.  The
zero-bump slope is subtracted from both measurement and map, so the inverse is
driven by bump-dependent K2 response rather than an absolute nominal baseline.

- targets / realizations: {nt} / {nr} per target
- descriptor: `{map_metadata['descriptor']}`
- maximum saved offset order: {map_metadata['maximum_offset_order']}
- BPM channels: {2 * nd}

| method | beam-relative 2D RMSE [um] | median [um] | P90 [um] | maximum [um] |
|---|---:|---:|---:|---:|
{lines}

These direct-map methods are more model-dependent than the scan-profiled
physical-source and empirical Taylor-surface methods.  Other 75 sextupole
offsets are hidden and are not used to condition the nominal GTPSA map.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

