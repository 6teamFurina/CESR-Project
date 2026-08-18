#!/usr/bin/env python3
"""Recover the sextupole center directly from the exact 3x3 bump curvature."""

from __future__ import annotations

import argparse
import csv
import json
import math
import tomllib
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--scan-dir", type=Path, required=True)
    result.add_argument("--coordinate-scale-m", type=float, default=5.0e-4)
    result.add_argument(
        "--bpm-noise-m",
        type=float,
        default=5.0e-6,
        help="Independent raw orbit noise per BPM plane and K2 state.",
    )
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parser().parse_args()
    scan_dir = args.scan_dir.expanduser().resolve()
    metadata = tomllib.loads((scan_dir / "scan_metadata.toml").read_text(encoding="utf-8"))
    observations = np.load(scan_dir / "bpm_orbits.npy")
    target_orbits = np.load(scan_dir / "target_orbits.npy")
    scenarios = [
        line.strip()
        for line in (scan_dir / "scenario_labels.txt").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    truth_index = scenarios.index("truth")
    conditions = read_rows(scan_dir / "k1_conditions.csv")
    candidates = list(metadata["candidate_quadrupoles"])
    levels = np.asarray(metadata["k2_levels"], dtype=float) * float(metadata["k2_step_m3"])
    k2_zero = int(np.flatnonzero(levels == 0.0)[0])
    slopes = np.einsum("k,cbkmp->cbmp", levels, observations[truth_index]) / float(
        levels @ levels
    )
    local_orbits = target_orbits[truth_index, :, :, k2_zero, :] / args.coordinate_scale_m

    condition_by_candidate: dict[str, tuple[int, int]] = {}
    for candidate in candidates:
        indices = tuple(
            index
            for index, row in enumerate(conditions)
            if row["quadrupole"] == candidate
        )
        if len(indices) != 2:
            raise ValueError(f"Expected two K1 conditions for {candidate}")
        condition_by_candidate[candidate] = indices

    truth = np.asarray(
        [
            float(metadata.get("baseline_x_offset_m", 0.0))
            + float(metadata["true_x_offset_m"]),
            float(metadata.get("baseline_y_offset_m", 0.0))
            + float(metadata["true_y_offset_m"]),
        ]
    )
    truth_scaled = truth / args.coordinate_scale_m
    slope_variance = args.bpm_noise_m**2 / float(levels @ levels)
    equations_by_condition: list[dict[str, Any]] = []
    for condition_index in range(11):
        x = local_orbits[condition_index, :, 0]
        y = local_orbits[condition_index, :, 1]
        polynomial = np.column_stack(
            [np.ones_like(x), x, y, x * x, x * y, y * y]
        )
        coefficients, _, rank, _ = np.linalg.lstsq(
            polynomial, slopes[condition_index].reshape(9, -1), rcond=None
        )
        if rank != 6:
            raise RuntimeError(f"Rank-deficient 3x3 bump polynomial for condition {condition_index}")
        linear_x = coefficients[1]
        linear_y = coefficients[2]
        quadratic_xx = coefficients[3]
        quadratic_xy = coefficients[4]
        quadratic_yy = coefficients[5]
        equation = np.empty((2 * linear_x.size, 2))
        rhs = np.empty(2 * linear_x.size)
        equation[0::2, 0] = 2.0 * quadratic_xx
        equation[0::2, 1] = quadratic_xy
        equation[1::2, 0] = quadratic_xy
        equation[1::2, 1] = 2.0 * quadratic_yy
        rhs[0::2] = -linear_x
        rhs[1::2] = -linear_y
        fitted = polynomial @ coefficients
        residual_rms = float(
            np.sqrt(np.mean((fitted - slopes[condition_index].reshape(9, -1)) ** 2))
        )
        coefficient_covariance = slope_variance * np.linalg.inv(polynomial.T @ polynomial)
        residual_projection = np.asarray(
            [
                [0.0, 1.0, 0.0, 2.0 * truth_scaled[0], truth_scaled[1], 0.0],
                [0.0, 0.0, 1.0, 0.0, truth_scaled[0], 2.0 * truth_scaled[1]],
            ]
        )
        residual_covariance = (
            residual_projection @ coefficient_covariance @ residual_projection.T
        )
        residual_precision = np.linalg.inv(residual_covariance)
        information = np.zeros((2, 2))
        normal_rhs = np.zeros(2)
        for output_index in range(linear_x.size):
            output_slice = slice(2 * output_index, 2 * output_index + 2)
            output_equation = equation[output_slice]
            output_rhs = rhs[output_slice]
            information += output_equation.T @ residual_precision @ output_equation
            normal_rhs += output_equation.T @ residual_precision @ output_rhs
        equations_by_condition.append(
            {
                "equation": equation,
                "rhs": rhs,
                "residual_rms": residual_rms,
                "information": information,
                "normal_rhs": normal_rhs,
            }
        )
    rows: list[dict[str, Any]] = []
    repeated_nominal_information = 7.0 * equations_by_condition[0]["information"]
    repeated_nominal_sign, repeated_nominal_logdet = np.linalg.slogdet(
        repeated_nominal_information
    )
    if repeated_nominal_sign <= 0:
        raise RuntimeError("Non-positive repeated-nominal center information")
    repeated_nominal_covariance = (
        args.coordinate_scale_m**2 * np.linalg.inv(repeated_nominal_information)
    )
    repeated_nominal_worst_axis_sigma = math.sqrt(
        float(np.max(np.linalg.eigvalsh(repeated_nominal_covariance)))
    )
    for triplet in combinations(candidates, 3):
        selected_conditions = [0]
        for candidate in triplet:
            selected_conditions.extend(condition_by_candidate[candidate])
        equation = np.vstack(
            [equations_by_condition[index]["equation"] for index in selected_conditions]
        )
        rhs = np.concatenate(
            [equations_by_condition[index]["rhs"] for index in selected_conditions]
        )
        unweighted_center_scaled, _, rank, singular_values = np.linalg.lstsq(
            equation, rhs, rcond=None
        )
        if rank != 2:
            raise RuntimeError(f"Rank-deficient center equations for {triplet}")
        information = sum(
            (
                equations_by_condition[index]["information"]
                for index in selected_conditions
            ),
            start=np.zeros((2, 2)),
        )
        normal_rhs = sum(
            (
                equations_by_condition[index]["normal_rhs"]
                for index in selected_conditions
            ),
            start=np.zeros(2),
        )
        center_scaled = np.linalg.solve(information, normal_rhs)
        center = center_scaled * args.coordinate_scale_m
        unweighted_center = unweighted_center_scaled * args.coordinate_scale_m
        error = center - truth
        sign, logdet = np.linalg.slogdet(information)
        if sign <= 0:
            raise RuntimeError(f"Non-positive center information for {triplet}")
        center_covariance = (
            args.coordinate_scale_m**2 * np.linalg.inv(information)
        )
        center_standard_deviation = np.sqrt(np.diag(center_covariance))
        worst_axis_standard_deviation = math.sqrt(
            float(np.max(np.linalg.eigvalsh(center_covariance)))
        )
        rows.append(
            {
                "target_sextupole": metadata["target_sextupole"],
                "quadrupole_1": triplet[0],
                "quadrupole_2": triplet[1],
                "quadrupole_3": triplet[2],
                "estimated_x_center_um": 1.0e6 * center[0],
                "estimated_y_center_um": 1.0e6 * center[1],
                "x_bias_um": 1.0e6 * error[0],
                "y_bias_um": 1.0e6 * error[1],
                "position_error_um": 1.0e6 * float(np.linalg.norm(error)),
                "unweighted_x_center_um": 1.0e6 * unweighted_center[0],
                "unweighted_y_center_um": 1.0e6 * unweighted_center[1],
                "quadratic_center_information_logdet": float(logdet),
                "matched_nominal_information_gain_logdet": float(
                    logdet - repeated_nominal_logdet
                ),
                "predicted_x_sigma_um": 1.0e6 * center_standard_deviation[0],
                "predicted_y_sigma_um": 1.0e6 * center_standard_deviation[1],
                "predicted_worst_axis_sigma_um": 1.0e6
                * worst_axis_standard_deviation,
                "matched_nominal_precision_improvement_worst_axis": float(
                    repeated_nominal_worst_axis_sigma
                    / worst_axis_standard_deviation
                ),
                "center_equation_condition_number": float(
                    singular_values[0] / singular_values[-1]
                ),
                "maximum_bump_polynomial_residual_rms": max(
                    equations_by_condition[index]["residual_rms"]
                    for index in selected_conditions
                ),
            }
        )

    information_order = sorted(
        rows,
        key=lambda row: float(row["quadratic_center_information_logdet"]),
        reverse=True,
    )
    error_order = sorted(rows, key=lambda row: float(row["position_error_um"]))
    precision_order = sorted(
        rows, key=lambda row: float(row["predicted_worst_axis_sigma_um"])
    )
    for rank, row in enumerate(information_order, start=1):
        row["information_rank"] = rank
    for rank, row in enumerate(error_order, start=1):
        row["error_rank"] = rank
    for rank, row in enumerate(precision_order, start=1):
        row["precision_rank"] = rank
    rows.sort(key=lambda row: int(row["information_rank"]))
    write_csv(scan_dir / "quadratic_bump_triplet_inversion.csv", rows)

    chosen = information_order[0]
    summary = {
        "format": "cesr-repaired-lattice-exact-quadratic-bump-center-v1",
        "target_sextupole": metadata["target_sextupole"],
        "candidate_quadrupoles": candidates,
        "truth_center_um": (1.0e6 * truth).tolist(),
        "triplet_count": 10,
        "coordinate_scale_m": args.coordinate_scale_m,
        "assumed_independent_bpm_noise_m": args.bpm_noise_m,
        "best_information_triplet": [
            chosen["quadrupole_1"],
            chosen["quadrupole_2"],
            chosen["quadrupole_3"],
        ],
        "best_information_position_error_um": chosen["position_error_um"],
        "best_information_predicted_worst_axis_sigma_um": chosen[
            "predicted_worst_axis_sigma_um"
        ],
        "best_precision_triplet": [
            precision_order[0]["quadrupole_1"],
            precision_order[0]["quadrupole_2"],
            precision_order[0]["quadrupole_3"],
        ],
        "minimum_predicted_worst_axis_sigma_um": precision_order[0][
            "predicted_worst_axis_sigma_um"
        ],
        "best_error_triplet": [
            error_order[0]["quadrupole_1"],
            error_order[0]["quadrupole_2"],
            error_order[0]["quadrupole_3"],
        ],
        "minimum_position_error_um": error_order[0]["position_error_um"],
        "limitations": [
            "The direct quadratic center fit is noise-free and does not yet marginalize other-sextupole or calibration nuisance.",
            "Only BPM closed-orbit K2 slopes are used; direct phase/coupling observables are deferred.",
            "Reported precision assumes independent 5 um raw BPM-plane noise in every K2 state and propagates it through the slope and bump-polynomial fits.",
            "Corrector-generated bumps are model-based and not machine-approved knobs.",
        ],
    }
    (scan_dir / "quadratic_bump_triplet_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
