#!/usr/bin/env python3
"""Choose the best three-quadrupole subset from each retained five-candidate set.

The available design consists of the nominal optics condition and symmetric
K1 +/- conditions for each of the five retained quadrupoles (11 conditions in
total).  All ten three-of-five subsets are evaluated with the same
nuisance-marginalized information model.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

from analyze_affinity import (
    LATEST_RESULTS,
    covariance_metrics,
    load_target_bundle,
    marginalized_information,
    nuisance_inverse,
    slope_noise,
    target_bundles,
    write_csv,
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--response-dir", type=Path, default=LATEST_RESULTS / "responses")
    result.add_argument("--selection-dir", type=Path, default=LATEST_RESULTS / "selection")
    result.add_argument("--nuisance-rms-m", type=float, default=3.0e-4)
    result.add_argument("--k2-step-m3", type=float, default=0.01)
    result.add_argument("--k2-levels", default="-2,-1,0,1,2")
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def positive_logdet(matrix: np.ndarray, label: str) -> float:
    sign, value = np.linalg.slogdet(matrix)
    if sign <= 0.0 or not math.isfinite(float(value)):
        raise RuntimeError(f"Non-positive or non-finite information determinant: {label}")
    return float(value)


def distribution(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "maximum": float(np.max(array)),
    }


def main() -> int:
    args = parser().parse_args()
    response_dir = args.response_dir.expanduser().resolve()
    selection_dir = args.selection_dir.expanduser().resolve()
    levels = [float(value) for value in args.k2_levels.split(",")]
    if not math.isfinite(args.nuisance_rms_m) or args.nuisance_rms_m <= 0.0:
        raise ValueError("--nuisance-rms-m must be positive and finite")

    five_rows = read_rows(selection_dir / "greedy_quadrupole_selection.csv")
    five_by_target: dict[str, list[dict[str, str]]] = {}
    for row in five_rows:
        five_by_target.setdefault(row["sextupole"], []).append(row)
    for target in five_by_target:
        five_by_target[target].sort(key=lambda row: int(row["greedy_step"]))

    all_triplets: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    changed_count = 0

    for target_path in target_bundles(response_dir):
        target, labels, nominal, nuisance, candidates, plus, minus = load_target_bundle(
            target_path, response_dir
        )
        retained = [row["quadrupole"] for row in five_by_target[target]]
        if len(retained) != 5 or len(set(retained)) != 5:
            raise RuntimeError(f"Expected five unique retained candidates for {target}")
        candidate_index = {name: index for index, name in enumerate(candidates)}
        if any(name not in candidate_index for name in retained):
            raise RuntimeError(f"Retained candidate missing from response bundle for {target}")

        sigma = slope_noise(labels, args.k2_step_m3, levels)
        whitened_nuisance = nuisance / sigma[:, None]
        baseline_inverse_nuisance = nuisance_inverse(
            whitened_nuisance, 1, args.nuisance_rms_m
        )
        baseline_information = marginalized_information(
            [nominal], whitened_nuisance, sigma, baseline_inverse_nuisance
        )
        baseline_logdet = positive_logdet(baseline_information, f"{target}/baseline")
        baseline_sigma, baseline_worst_axis, baseline_worst_direction = covariance_metrics(
            baseline_information
        )

        # Every three-knob trial uses nominal plus six one-at-a-time K1 blocks.
        triplet_inverse_nuisance = nuisance_inverse(
            whitened_nuisance, 7, args.nuisance_rms_m
        )
        target_triplets: list[dict[str, Any]] = []
        for triplet in combinations(retained, 3):
            blocks = [nominal]
            for name in triplet:
                index = candidate_index[name]
                blocks.extend((plus[index], minus[index]))
            information = marginalized_information(
                blocks, whitened_nuisance, sigma, triplet_inverse_nuisance
            )
            sigmas, worst_axis, worst_direction = covariance_metrics(information)
            logdet = positive_logdet(information, f"{target}/{'|'.join(triplet)}")
            target_triplets.append(
                {
                    "sextupole": target,
                    "sextupole_s_m": float(five_by_target[target][0]["sextupole_s_m"]),
                    "quadrupole_1": triplet[0],
                    "quadrupole_2": triplet[1],
                    "quadrupole_3": triplet[2],
                    "information_gain_logdet": logdet - baseline_logdet,
                    "precision_improvement_worst_axis": baseline_worst_axis / worst_axis,
                    "precision_improvement_worst_direction": baseline_worst_direction
                    / worst_direction,
                    "sigma_x_um": 1.0e6 * float(sigmas[0]),
                    "sigma_y_um": 1.0e6 * float(sigmas[1]),
                }
            )

        information_order = sorted(
            target_triplets,
            key=lambda row: (
                float(row["information_gain_logdet"]),
                float(row["precision_improvement_worst_axis"]),
                str(row["quadrupole_1"]),
                str(row["quadrupole_2"]),
                str(row["quadrupole_3"]),
            ),
            reverse=True,
        )
        precision_order = sorted(
            target_triplets,
            key=lambda row: (
                float(row["precision_improvement_worst_axis"]),
                float(row["information_gain_logdet"]),
            ),
            reverse=True,
        )
        information_rank = {id(row): rank for rank, row in enumerate(information_order, start=1)}
        precision_rank = {id(row): rank for rank, row in enumerate(precision_order, start=1)}
        for row in target_triplets:
            row["information_rank"] = information_rank[id(row)]
            row["precision_rank"] = precision_rank[id(row)]
            row["selected"] = int(row is information_order[0])
        all_triplets.extend(target_triplets)

        chosen = information_order[0]
        greedy_prefix = set(retained[:3])
        chosen_set = {
            str(chosen["quadrupole_1"]),
            str(chosen["quadrupole_2"]),
            str(chosen["quadrupole_3"]),
        }
        changed = chosen_set != greedy_prefix
        changed_count += int(changed)
        greedy_row = next(
            row
            for row in target_triplets
            if {str(row["quadrupole_1"]), str(row["quadrupole_2"]), str(row["quadrupole_3"])}
            == greedy_prefix
        )
        best_precision = precision_order[0]
        summaries.append(
            {
                "sextupole": target,
                "sextupole_s_m": chosen["sextupole_s_m"],
                "candidate_1": retained[0],
                "candidate_2": retained[1],
                "candidate_3": retained[2],
                "candidate_4": retained[3],
                "candidate_5": retained[4],
                "selected_quadrupole_1": chosen["quadrupole_1"],
                "selected_quadrupole_2": chosen["quadrupole_2"],
                "selected_quadrupole_3": chosen["quadrupole_3"],
                "selected_information_gain_logdet": chosen["information_gain_logdet"],
                "selected_precision_improvement_worst_axis": chosen[
                    "precision_improvement_worst_axis"
                ],
                "selected_precision_rank": chosen["precision_rank"],
                "greedy_prefix_changed": int(changed),
                "greedy_prefix_information_gain_logdet": greedy_row[
                    "information_gain_logdet"
                ],
                "information_gain_over_greedy_prefix": float(
                    chosen["information_gain_logdet"]
                )
                - float(greedy_row["information_gain_logdet"]),
                "best_precision_quadrupole_1": best_precision["quadrupole_1"],
                "best_precision_quadrupole_2": best_precision["quadrupole_2"],
                "best_precision_quadrupole_3": best_precision["quadrupole_3"],
                "best_precision_improvement_worst_axis": best_precision[
                    "precision_improvement_worst_axis"
                ],
            }
        )

    all_triplets.sort(
        key=lambda row: (float(row["sextupole_s_m"]), int(row["information_rank"]))
    )
    summaries.sort(key=lambda row: float(row["sextupole_s_m"]))
    write_csv(selection_dir / "triplet_combinations.csv", all_triplets)
    write_csv(selection_dir / "best_triplets_by_sextupole.csv", summaries)

    chosen_information = [
        float(row["selected_information_gain_logdet"]) for row in summaries
    ]
    chosen_precision = [
        float(row["selected_precision_improvement_worst_axis"]) for row in summaries
    ]
    gain_over_prefix = [
        float(row["information_gain_over_greedy_prefix"]) for row in summaries
    ]
    information_precision_disagreement = sum(
        {
            row["selected_quadrupole_1"],
            row["selected_quadrupole_2"],
            row["selected_quadrupole_3"],
        }
        != {
            row["best_precision_quadrupole_1"],
            row["best_precision_quadrupole_2"],
            row["best_precision_quadrupole_3"],
        }
        for row in summaries
    )
    summary = {
        "target_count": len(summaries),
        "retained_candidate_count_per_target": 5,
        "available_k1_condition_count": 11,
        "triplet_count_per_target": 10,
        "selected_triplet_size": 3,
        "targets_changed_from_greedy_prefix": changed_count,
        "information_vs_precision_triplet_disagreement_count": information_precision_disagreement,
        "selected_information_gain_logdet": distribution(chosen_information),
        "selected_precision_improvement_worst_axis": distribution(chosen_precision),
        "information_gain_over_greedy_prefix": distribution(gain_over_prefix),
        "selected_triplet_union_quadrupole_count": len(
            {
                str(row[key])
                for row in summaries
                for key in (
                    "selected_quadrupole_1",
                    "selected_quadrupole_2",
                    "selected_quadrupole_3",
                )
            }
        ),
    }
    (selection_dir / "triplet_selection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    metadata = {
        "format": "cesr-sextupole-five-to-three-quadrupole-selection-v1",
        "engine": "SciBmad/GTPSA repaired-lattice response bundles",
        "selection_objective": "maximum nuisance-marginalized log determinant among all 3-of-5 subsets; worst-axis precision is the tie-breaker",
        "available_condition_design": "nominal plus K1 +/- for each retained quadrupole (11 one-at-a-time conditions)",
        "evaluated_condition_count_per_triplet": 7,
        "nuisance_rms_m": args.nuisance_rms_m,
        "k2_step_m3": args.k2_step_m3,
        "k2_levels": levels,
        "limitations": [
            "This is the response-level 11-condition selection, not yet the exact 3 x 3 bump by five-point K2 finite-amplitude scan.",
            "Other-sextupole nuisance responses are evaluated at nominal optics and reused in every K1 block.",
            "Combined K1 states are not part of this one-at-a-time design; they belong to the subsequent three-knob 3^3 validation.",
        ],
    }
    (selection_dir / "triplet_selection_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Targets: {len(summaries)}")
    print(f"Triplets evaluated: {len(all_triplets)}")
    print(f"Changed from greedy prefix: {changed_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
