#!/usr/bin/env python3
"""Invert exact 11-condition scans and compare all three-of-five triplets."""

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

HERE = Path(__file__).resolve().parent


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--scan-dir", type=Path, required=True)
    result.add_argument("--bpm-noise-m", type=float, default=5.0e-6)
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parser().parse_args()
    scan_dir = args.scan_dir.expanduser().resolve()
    metadata = tomllib.loads((scan_dir / "scan_metadata.toml").read_text(encoding="utf-8"))
    observations = np.load(scan_dir / "bpm_orbits.npy")
    scenarios = [
        line.strip()
        for line in (scan_dir / "scenario_labels.txt").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    conditions = read_rows(scan_dir / "k1_conditions.csv")
    expected_shape = (6, 11, 9, 5, int(metadata["bpm_count"]), 2)
    if observations.shape != expected_shape:
        raise ValueError(f"Unexpected observation shape {observations.shape}, expected {expected_shape}")
    scenario_index = {name: index for index, name in enumerate(scenarios)}
    levels = np.asarray(metadata["k2_levels"], dtype=float) * float(metadata["k2_step_m3"])
    slope = np.einsum("k,scbkmp->scbmp", levels, observations) / float(levels @ levels)
    zero = slope[scenario_index["zero"]]
    truth = slope[scenario_index["truth"]]
    fd_step = float(metadata["offset_fd_step_m"])
    jacobian_x = (
        slope[scenario_index["x_plus"]] - slope[scenario_index["x_minus"]]
    ) / (2.0 * fd_step)
    jacobian_y = (
        slope[scenario_index["y_plus"]] - slope[scenario_index["y_minus"]]
    ) / (2.0 * fd_step)
    residual = truth - zero
    candidates = list(metadata["candidate_quadrupoles"])
    condition_by_candidate: dict[str, tuple[int, int]] = {}
    for candidate in candidates:
        indices = [
            index
            for index, row in enumerate(conditions)
            if row["quadrupole"] == candidate
        ]
        if len(indices) != 2:
            raise ValueError(f"Expected two conditions for {candidate}")
        condition_by_candidate[candidate] = (indices[0], indices[1])

    truth_center = np.asarray(
        [float(metadata["true_x_offset_m"]), float(metadata["true_y_offset_m"])]
    )
    slope_noise = args.bpm_noise_m / math.sqrt(float(levels @ levels))
    rows: list[dict[str, Any]] = []
    for triplet in combinations(candidates, 3):
        selected_conditions = [0]
        for candidate in triplet:
            selected_conditions.extend(condition_by_candidate[candidate])
        y = residual[selected_conditions].reshape(-1)
        design = np.column_stack(
            [
                jacobian_x[selected_conditions].reshape(-1),
                jacobian_y[selected_conditions].reshape(-1),
            ]
        )
        estimate, _, rank, singular_values = np.linalg.lstsq(design, y, rcond=None)
        information = (design.T @ design) / slope_noise**2
        covariance = np.linalg.inv(information)
        sigmas = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        sign, logdet = np.linalg.slogdet(information)
        if sign <= 0 or rank != 2:
            raise RuntimeError(f"Rank-deficient triplet {triplet}")
        error = estimate - truth_center
        rows.append(
            {
                "target_sextupole": metadata["target_sextupole"],
                "quadrupole_1": triplet[0],
                "quadrupole_2": triplet[1],
                "quadrupole_3": triplet[2],
                "estimated_x_offset_um": 1.0e6 * estimate[0],
                "estimated_y_offset_um": 1.0e6 * estimate[1],
                "x_bias_um": 1.0e6 * error[0],
                "y_bias_um": 1.0e6 * error[1],
                "position_error_um": 1.0e6 * float(np.linalg.norm(error)),
                "information_logdet": float(logdet),
                "predicted_sigma_x_um": 1.0e6 * sigmas[0],
                "predicted_sigma_y_um": 1.0e6 * sigmas[1],
                "predicted_worst_axis_sigma_um": 1.0e6 * float(np.max(sigmas)),
                "design_condition_number": float(singular_values[0] / singular_values[-1]),
            }
        )

    info_order = sorted(rows, key=lambda row: float(row["information_logdet"]), reverse=True)
    error_order = sorted(rows, key=lambda row: float(row["position_error_um"]))
    precision_order = sorted(rows, key=lambda row: float(row["predicted_worst_axis_sigma_um"]))
    for rank, row in enumerate(info_order, start=1):
        row["information_rank"] = rank
    for rank, row in enumerate(error_order, start=1):
        row["error_rank"] = rank
    for rank, row in enumerate(precision_order, start=1):
        row["precision_rank"] = rank
    rows.sort(key=lambda row: int(row["information_rank"]))
    write_csv(scan_dir / "exact_triplet_inversion.csv", rows)

    chosen = info_order[0]
    summary = {
        "format": "cesr-repaired-lattice-exact-11-triplet-inversion-v1",
        "target_sextupole": metadata["target_sextupole"],
        "candidate_quadrupoles": candidates,
        "truth_offset_um": (1.0e6 * truth_center).tolist(),
        "triplet_count": len(rows),
        "bpm_noise_m": args.bpm_noise_m,
        "slope_noise_per_readback": slope_noise,
        "best_information_triplet": [
            chosen["quadrupole_1"],
            chosen["quadrupole_2"],
            chosen["quadrupole_3"],
        ],
        "best_information_position_error_um": chosen["position_error_um"],
        "best_error_triplet": [
            error_order[0]["quadrupole_1"],
            error_order[0]["quadrupole_2"],
            error_order[0]["quadrupole_3"],
        ],
        "minimum_position_error_um": error_order[0]["position_error_um"],
        "best_precision_triplet": [
            precision_order[0]["quadrupole_1"],
            precision_order[0]["quadrupole_2"],
            precision_order[0]["quadrupole_3"],
        ],
        "best_precision_worst_axis_sigma_um": precision_order[0][
            "predicted_worst_axis_sigma_um"
        ],
        "limitations": [
            "Nominal target-offset scenario only; other-sextupole and calibration nuisances are not yet included.",
            "Closed orbit is the only observable group in this pilot; direct phase/coupling ablation is deferred.",
            "The inverse is a central-finite-difference local linear model evaluated against an exact finite-offset scan.",
        ],
    }
    (scan_dir / "exact_triplet_inversion_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
