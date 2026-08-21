#!/usr/bin/env python3
"""Aggregate repaired-lattice exact-11 triplet validation results."""

from __future__ import annotations

import argparse
import csv
import json
import math
import tomllib
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_SELECTION = (
    HERE.parent
    / "results"
    / "scibmad_latest"
    / "selection"
    / "best_triplets_by_sextupole.csv"
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--scan-root", type=Path, default=HERE / "results" / "scans")
    result.add_argument("--selection-csv", type=Path, default=DEFAULT_SELECTION)
    result.add_argument("--output-dir", type=Path, default=HERE / "results" / "aggregate")
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def triplet_set(row: dict[str, Any], prefix: str) -> frozenset[str]:
    return frozenset(str(row[f"{prefix}_{index}"]) for index in range(1, 4))


def distribution(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "maximum": float(np.max(array)),
    }


def main() -> int:
    args = parser().parse_args()
    scan_root = args.scan_root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    selected_rows = read_rows(args.selection_csv.expanduser().resolve())
    selected_by_target = {row["sextupole"]: row for row in selected_rows}
    rows: list[dict[str, Any]] = []

    for scan_dir in sorted(path for path in scan_root.iterdir() if path.is_dir()):
        metadata_path = scan_dir / "scan_metadata.toml"
        result_path = scan_dir / "quadratic_bump_triplet_inversion.csv"
        if not metadata_path.exists() or not result_path.exists():
            continue
        metadata = tomllib.loads(metadata_path.read_text(encoding="utf-8"))
        target = str(metadata["target_sextupole"])
        triplets = read_rows(result_path)
        if len(triplets) != 10:
            raise ValueError(f"Expected ten triplets for {target}, found {len(triplets)}")
        exact_info = min(triplets, key=lambda row: int(row["information_rank"]))
        exact_precision = min(triplets, key=lambda row: int(row["precision_rank"]))
        exact_error = min(triplets, key=lambda row: int(row["error_rank"]))
        response_selection = selected_by_target[target]
        response_prefix = (
            "selected_quadrupole"
            if "selected_quadrupole_1" in response_selection
            else "candidate"
        )
        response_set = frozenset(
            response_selection[f"{response_prefix}_{index}"] for index in range(1, 4)
        )
        response_exact = next(
            row
            for row in triplets
            if frozenset(row[f"quadrupole_{index}"] for index in range(1, 4))
            == response_set
        )
        exact_info_set = frozenset(
            exact_info[f"quadrupole_{index}"] for index in range(1, 4)
        )
        exact_error_set = frozenset(
            exact_error[f"quadrupole_{index}"] for index in range(1, 4)
        )
        exact_precision_set = frozenset(
            exact_precision[f"quadrupole_{index}"] for index in range(1, 4)
        )

        target_orbits = np.load(scan_dir / "target_orbits.npy")
        scenario_labels = [
            line.strip()
            for line in (scan_dir / "scenario_labels.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        truth_index = scenario_labels.index("truth")
        levels = np.asarray(metadata["k2_levels"], dtype=float)
        zero_k2 = int(np.flatnonzero(levels == 0.0)[0])
        nominal_target_orbit = target_orbits[truth_index, 0, :, zero_k2, :]
        center_index = 4
        realized_bumps = nominal_target_orbit - nominal_target_orbit[center_index]
        bump_commands = read_rows(scan_dir / "bump_points.csv")
        commands = np.asarray(
            [
                [float(row["bump_x_command_m"]), float(row["bump_y_command_m"])]
                for row in bump_commands
            ]
        )
        bump_error = realized_bumps - commands

        rows.append(
            {
                "sextupole": target,
                "response_selected_1": response_selection[f"{response_prefix}_1"],
                "response_selected_2": response_selection[f"{response_prefix}_2"],
                "response_selected_3": response_selection[f"{response_prefix}_3"],
                "exact_information_1": exact_info["quadrupole_1"],
                "exact_information_2": exact_info["quadrupole_2"],
                "exact_information_3": exact_info["quadrupole_3"],
                "exact_error_1": exact_error["quadrupole_1"],
                "exact_error_2": exact_error["quadrupole_2"],
                "exact_error_3": exact_error["quadrupole_3"],
                "exact_precision_1": exact_precision["quadrupole_1"],
                "exact_precision_2": exact_precision["quadrupole_2"],
                "exact_precision_3": exact_precision["quadrupole_3"],
                "response_matches_exact_information": int(response_set == exact_info_set),
                "response_matches_exact_precision": int(response_set == exact_precision_set),
                "response_matches_exact_error": int(response_set == exact_error_set),
                "response_triplet_exact_information_rank": int(response_exact["information_rank"]),
                "response_triplet_exact_precision_rank": int(response_exact["precision_rank"]),
                "response_triplet_exact_error_rank": int(response_exact["error_rank"]),
                "response_triplet_predicted_worst_axis_sigma_um": float(
                    response_exact["predicted_worst_axis_sigma_um"]
                ),
                "best_exact_predicted_worst_axis_sigma_um": float(
                    exact_precision["predicted_worst_axis_sigma_um"]
                ),
                "predicted_worst_axis_sigma_penalty_um": float(
                    response_exact["predicted_worst_axis_sigma_um"]
                )
                - float(exact_precision["predicted_worst_axis_sigma_um"]),
                "response_triplet_position_error_um": float(response_exact["position_error_um"]),
                "best_exact_position_error_um": float(exact_error["position_error_um"]),
                "position_error_penalty_um": float(response_exact["position_error_um"])
                - float(exact_error["position_error_um"]),
                "exact_information_logdet_gap": float(
                    exact_info["quadratic_center_information_logdet"]
                )
                - float(response_exact["quadratic_center_information_logdet"]),
                "max_abs_realized_bump_error_um": 1.0e6 * float(np.max(np.abs(bump_error))),
                "max_bump_polynomial_residual_rms": float(
                    response_exact["maximum_bump_polynomial_residual_rms"]
                ),
                "response_triplet_matched_nominal_information_gain_logdet": float(
                    response_exact["matched_nominal_information_gain_logdet"]
                ),
                "response_triplet_matched_nominal_precision_improvement": float(
                    response_exact[
                        "matched_nominal_precision_improvement_worst_axis"
                    ]
                ),
                "exact_information_triplet_matched_nominal_information_gain_logdet": float(
                    exact_info["matched_nominal_information_gain_logdet"]
                ),
                "exact_precision_triplet_matched_nominal_precision_improvement": float(
                    exact_precision[
                        "matched_nominal_precision_improvement_worst_axis"
                    ]
                ),
            }
        )

    if len(rows) != 76:
        raise RuntimeError(f"Expected 76 completed targets, found {len(rows)}")
    rows.sort(key=lambda row: row["sextupole"])
    write_csv(output_dir / "exact11_validation_by_sextupole.csv", rows)
    summary = {
        "format": "cesr-repaired-lattice-exact-11-validation-summary-v1",
        "target_count": len(rows),
        "response_selection_matches_exact_information_count": sum(
            int(row["response_matches_exact_information"]) for row in rows
        ),
        "response_selection_matches_exact_error_count": sum(
            int(row["response_matches_exact_error"]) for row in rows
        ),
        "response_selection_matches_exact_precision_count": sum(
            int(row["response_matches_exact_precision"]) for row in rows
        ),
        "exact_information_matches_exact_precision_count": sum(
            int(
                frozenset(row[f"exact_information_{index}"] for index in range(1, 4))
                == frozenset(row[f"exact_precision_{index}"] for index in range(1, 4))
            )
            for row in rows
        ),
        "response_triplet_exact_information_rank": distribution(
            [float(row["response_triplet_exact_information_rank"]) for row in rows]
        ),
        "response_triplet_position_error_um": distribution(
            [float(row["response_triplet_position_error_um"]) for row in rows]
        ),
        "response_triplet_exact_precision_rank": distribution(
            [float(row["response_triplet_exact_precision_rank"]) for row in rows]
        ),
        "response_triplet_predicted_worst_axis_sigma_um": distribution(
            [float(row["response_triplet_predicted_worst_axis_sigma_um"]) for row in rows]
        ),
        "predicted_worst_axis_sigma_penalty_um": distribution(
            [float(row["predicted_worst_axis_sigma_penalty_um"]) for row in rows]
        ),
        "position_error_penalty_um": distribution(
            [float(row["position_error_penalty_um"]) for row in rows]
        ),
        "exact_information_logdet_gap": distribution(
            [float(row["exact_information_logdet_gap"]) for row in rows]
        ),
        "response_triplet_matched_nominal_information_gain_logdet": distribution(
            [
                float(
                    row[
                        "response_triplet_matched_nominal_information_gain_logdet"
                    ]
                )
                for row in rows
            ]
        ),
        "response_triplet_matched_nominal_precision_improvement": distribution(
            [
                float(
                    row[
                        "response_triplet_matched_nominal_precision_improvement"
                    ]
                )
                for row in rows
            ]
        ),
        "exact_information_triplet_matched_nominal_information_gain_logdet": distribution(
            [
                float(
                    row[
                        "exact_information_triplet_matched_nominal_information_gain_logdet"
                    ]
                )
                for row in rows
            ]
        ),
        "exact_precision_triplet_matched_nominal_precision_improvement": distribution(
            [
                float(
                    row[
                        "exact_precision_triplet_matched_nominal_precision_improvement"
                    ]
                )
                for row in rows
            ]
        ),
        "max_abs_realized_bump_error_um": distribution(
            [float(row["max_abs_realized_bump_error_um"]) for row in rows]
        ),
        "limitations": [
            "One nominal target-offset realization per sextupole; other-sextupole and calibration nuisance are not yet included.",
            "Closed-orbit K2 slopes only; direct phase/coupling observable ablation is deferred.",
            "Exact information and precision propagate a provisional independent 5 um raw BPM-plane noise model through K2-slope and bump-polynomial fits; use recovery penalty and later nuisance tests before replacing a response-level triplet.",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "exact11_validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
