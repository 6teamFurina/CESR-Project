#!/usr/bin/env python3
"""Fit a nonlinear target-offset model in the reconstructed four-source space."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import benchmark_physical_inverses as physical
import invert_scan as base


HERE = Path(__file__).resolve().parent
DEFAULT_SCAN = HERE / "results" / "smoke_background"
DEFAULT_MODELS = HERE / "results" / "unknown_background_benchmark" / "nominal_models"
DEFAULT_LINEAR = HERE / "results" / "unknown_background_benchmark" / "linear"
DEFAULT_OUTPUT = HERE / "results" / "unknown_background_benchmark" / "nonlinear_p2b"
FAMILIES = {"orbit", "phase", "coupling", "tune"}
SOURCE_COLUMNS = ("d_kn0l", "d_ks0l", "d_kn1l", "d_ks1l")


def calibration_groups(rows: list[dict[str, str]]) -> dict[int, dict[int, list[dict[str, str]]]]:
    grouped: dict[int, dict[int, list[dict[str, str]]]] = {}
    for row in rows:
        grouped.setdefault(int(row["bump_index"]), {}).setdefault(
            int(row["calibration_index"]), []
        ).append(row)
    return grouped


def features(x: float, y: float) -> np.ndarray:
    return np.asarray([x, y, 0.5 * x * x, x * y, 0.5 * y * y])


def feature_jacobian(x: float, y: float) -> np.ndarray:
    return np.asarray([
        [1.0, 0.0],
        [0.0, 1.0],
        [x, 0.0],
        [y, x],
        [0.0, y],
    ])


def source_inverse(source: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u, singular, vt = np.linalg.svd(source, full_matrices=False)
    threshold = np.finfo(float).eps * max(source.shape) * singular[0]
    rank = int(np.sum(singular > threshold))
    return u[:, :rank], singular[:rank], vt[:rank]


def reconstruct_q(u: np.ndarray, singular: np.ndarray, vt: np.ndarray, right: np.ndarray) -> np.ndarray:
    return vt.T @ ((u.T @ right) / singular)


def build_bump_models(
    slopes: dict[int, dict[tuple[str, str, str], float]],
    conditioned: dict[int, dict[tuple[str, str, str], dict[str, str]]],
    sources: dict[int, dict[tuple[str, str, str], dict[str, str]]],
    calibrations: dict[int, dict[int, list[dict[str, str]]]],
    scales: dict[str, float],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    models: list[dict[str, object]] = []
    source_rows: list[dict[str, object]] = []
    for bump in sorted(slopes):
        response_rows = conditioned[bump]
        source_rows_by_key = sources[bump]
        keys = [
            key for key in sorted(response_rows)
            if base.observable_family(key[2]) in FAMILIES
        ]
        scale = np.asarray([scales[key[2]] for key in keys])
        source = np.asarray([
            [float(source_rows_by_key[key][column]) for column in SOURCE_COLUMNS]
            for key in keys
        ]) / scale[:, None]
        u, singular, vt = source_inverse(source)
        reference = np.asarray([
            float(response_rows[key]["reference_k2_slope"]) for key in keys
        ])
        measured = np.asarray([slopes[bump][key] for key in keys])
        q_measured = reconstruct_q(u, singular, vt, (measured - reference) / scale)

        point_features, point_q = [], []
        for calibration_index in sorted(calibrations[bump]):
            rows = {physical.row_key(row): row for row in calibrations[bump][calibration_index]}
            first = next(iter(rows.values()))
            x = float(first["calibration_x_offset_m"])
            y = float(first["calibration_y_offset_m"])
            calibration_slope = np.asarray([float(rows[key]["k2_slope"]) for key in keys])
            q_calibration = reconstruct_q(
                u, singular, vt, (calibration_slope - reference) / scale
            )
            point_features.append(features(x, y))
            point_q.append(q_calibration)
        design = np.vstack(point_features)
        q_values = np.vstack(point_q)
        coefficients = np.linalg.lstsq(design, q_values, rcond=None)[0]
        calibration_residual = float(
            np.linalg.norm(design @ coefficients - q_values) / np.sqrt(q_values.size)
        )
        whitener = singular[:, None] * vt
        models.append({
            "bump": bump,
            "q_measured": q_measured,
            "coefficients": coefficients,
            "whitener": whitener,
            "calibration_residual": calibration_residual,
        })
        source_rows.append({
            "bump_index": bump,
            "source_rank": len(singular),
            "calibration_source_rms": calibration_residual,
            **{f"measured_{name}_per_k2": float(value) for name, value in zip(SOURCE_COLUMNS, q_measured)},
        })
    return models, source_rows


def fit_offset(models: list[dict[str, object]], initial: np.ndarray, iterations: int) -> list[dict[str, float | int]]:
    estimate = initial.copy()
    history: list[dict[str, float | int]] = []
    for iteration in range(iterations + 1):
        design_blocks, right_blocks, residual_blocks = [], [], []
        for model in models:
            coefficients = np.asarray(model["coefficients"])
            whitener = np.asarray(model["whitener"])
            q_measured = np.asarray(model["q_measured"])
            q_predicted = features(*estimate) @ coefficients
            q_jacobian = coefficients.T @ feature_jacobian(*estimate)
            residual = q_measured - q_predicted
            design_blocks.append(whitener @ q_jacobian)
            right_blocks.append(whitener @ residual)
            residual_blocks.append(whitener @ residual)
        design = np.vstack(design_blocks)
        right = np.concatenate(right_blocks)
        residual_rms = float(np.linalg.norm(np.concatenate(residual_blocks)) / np.sqrt(len(right)))
        history.append({
            "iteration": iteration,
            "estimated_x_offset_m": float(estimate[0]),
            "estimated_y_offset_m": float(estimate[1]),
            "source_fit_residual_rms": residual_rms,
        })
        if iteration == iterations:
            break
        update = np.linalg.lstsq(design, right, rcond=None)[0]
        estimate += np.clip(update, -2e-4, 2e-4)
    return history


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-dir", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--models-dir", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--linear-dir", type=Path, default=DEFAULT_LINEAR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--iterations", type=int, default=6)
    args = parser.parse_args()
    scan_dir, models_dir, linear_dir, output_dir = map(
        Path.resolve, (args.scan_dir, args.models_dir, args.linear_dir, args.output_dir)
    )
    observations = base.read_rows(scan_dir / "scan_observations.csv")
    states = base.read_rows(scan_dir / "scan_states.csv")
    slopes, _ = base.fit_slopes(observations)
    scales = base.load_scales(base.RESPONSE_MAP / "local_response_svd_scales.csv")
    conditioned = physical.conditioned_by_bump(
        base.read_rows(models_dir / "conditioned_mixed_response.csv")
    )
    sources = physical.conditioned_by_bump(
        base.read_rows(models_dir / "local_source_response.csv")
    )
    calibrations = calibration_groups(
        base.read_rows(models_dir / "nonlinear_offset_calibration_slopes.csv")
    )
    models, source_rows = build_bump_models(slopes, conditioned, sources, calibrations, scales)
    linear_rows = base.read_rows(linear_dir / "p0_p2_offset_estimates.csv")
    linear_p2b = next(row for row in linear_rows if row["method"] == "P2b_four_local_kicks")
    initial = np.asarray([
        float(linear_p2b["estimated_x_offset_m"]),
        float(linear_p2b["estimated_y_offset_m"]),
    ])
    truth = np.asarray([
        float(states[0]["true_x_offset_m"]), float(states[0]["true_y_offset_m"])
    ])
    history = fit_offset(models, initial, args.iterations)
    for row in history:
        estimate = np.asarray([row["estimated_x_offset_m"], row["estimated_y_offset_m"]])
        error = estimate - truth
        row.update({
            "true_x_offset_m": float(truth[0]),
            "true_y_offset_m": float(truth[1]),
            "error_x_m": float(error[0]),
            "error_y_m": float(error[1]),
            "absolute_error_2d_m": float(np.linalg.norm(error)),
        })
    output_dir.mkdir(parents=True, exist_ok=True)
    base.write_csv(output_dir / "nonlinear_p2b_history.csv", history)
    base.write_csv(output_dir / "reconstructed_sources.csv", source_rows)
    summary = {
        "format": "cesr-unknown-background-nonlinear-p2b-v1",
        "scan_dir": str(scan_dir),
        "model_background": "nominal; saved other-sextupole offsets are hidden",
        "nonlinear_model": "per-bump quadratic target-offset polynomial in reconstructed four-source space",
        "input_axis": "K2 slopes from the maintained three-point scan",
        "initial_linear_p2b": linear_p2b,
        "final_nonlinear_p2b": history[-1],
        "interpretation_boundary": (
            "This is a first nonlinear second-stage test; it is not yet the per-finite-K2 thick-sextupole fit."
        ),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Linear P2b error: {1e6*float(linear_p2b['absolute_error_2d_m']):.3f} um")
    for row in history:
        print(
            f"Nonlinear P2b iteration {int(row['iteration'])}: "
            f"x={1e6*float(row['estimated_x_offset_m']):+.3f} um "
            f"y={1e6*float(row['estimated_y_offset_m']):+.3f} um "
            f"error={1e6*float(row['absolute_error_2d_m']):.3f} um"
        )


if __name__ == "__main__":
    main()
