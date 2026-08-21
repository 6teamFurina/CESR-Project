#!/usr/bin/env python3
"""Select five complementary quadrupoles per sextupole from SciBmad responses.

The first three greedy selections are marked as the provisional operational
set.  All five remain available for the finite-amplitude and observable
ablation study.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
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
    result.add_argument("--affinity-dir", type=Path, default=LATEST_RESULTS / "affinity")
    result.add_argument("--output-dir", type=Path, default=LATEST_RESULTS / "selection")
    result.add_argument("--retain-count", type=int, default=5)
    result.add_argument("--operational-count", type=int, default=3)
    result.add_argument("--nuisance-rms-m", type=float, default=3.0e-4)
    result.add_argument("--k2-step-m3", type=float, default=0.01)
    result.add_argument("--k2-levels", default="-2,-1,0,1,2")
    result.add_argument(
        "--exclude-quadrupoles",
        default="",
        help="Comma-separated quadrupoles excluded by an operational audit.",
    )
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def positive_logdet(matrix: np.ndarray, label: str) -> float:
    sign, value = np.linalg.slogdet(matrix)
    if sign <= 0.0 or not math.isfinite(float(value)):
        raise RuntimeError(f"Non-positive or non-finite information determinant: {label}")
    return float(value)


def main() -> int:
    args = parser().parse_args()
    if args.retain_count < 1:
        raise ValueError("--retain-count must be positive")
    if not 1 <= args.operational_count <= args.retain_count:
        raise ValueError("--operational-count must lie in [1, retain-count]")
    if not math.isfinite(args.nuisance_rms_m) or args.nuisance_rms_m <= 0.0:
        raise ValueError("--nuisance-rms-m must be positive and finite")

    levels = [float(value) for value in args.k2_levels.split(",")]
    excluded_quadrupoles = {
        value.strip().upper()
        for value in args.exclude_quadrupoles.split(",")
        if value.strip()
    }
    response_dir = args.response_dir.expanduser().resolve()
    affinity_dir = args.affinity_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    score_rows = read_rows(affinity_dir / "quadrupole_affinity_scores.csv")
    score_lookup = {
        (row["sextupole"], row["quadrupole"]): row for row in score_rows
    }
    sextupole_s = {
        row["sextupole"]: float(row["sextupole_s_m"]) for row in score_rows
    }

    selection_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    target_paths = target_bundles(response_dir)
    if not target_paths:
        raise FileNotFoundError(f"No target response files under {response_dir / 'targets'}")

    for target_path in target_paths:
        target, labels, nominal, nuisance, candidates, plus, minus = load_target_bundle(
            target_path, response_dir
        )
        retained_indices = [
            index
            for index, candidate in enumerate(candidates)
            if candidate not in excluded_quadrupoles
        ]
        candidates = [candidates[index] for index in retained_indices]
        plus = plus[retained_indices]
        minus = minus[retained_indices]
        if args.retain_count > len(candidates):
            raise ValueError(
                f"{target} has {len(candidates)} candidates, fewer than retain-count={args.retain_count}"
            )
        sigma = slope_noise(labels, args.k2_step_m3, levels)
        whitened_nuisance = nuisance / sigma[:, None]
        inverse_nuisance_one = nuisance_inverse(
            whitened_nuisance, 1, args.nuisance_rms_m
        )
        baseline_information = marginalized_information(
            [nominal], whitened_nuisance, sigma, inverse_nuisance_one
        )
        baseline_sigma, baseline_worst_axis, baseline_worst_direction = covariance_metrics(
            baseline_information
        )
        baseline_logdet = positive_logdet(baseline_information, f"{target}/baseline")

        candidate_index = {name: index for index, name in enumerate(candidates)}
        selected: list[str] = []
        previous_logdet = baseline_logdet
        target_rows: list[dict[str, Any]] = []

        for step in range(1, args.retain_count + 1):
            # Every trial at this step has nominal plus two blocks per selected
            # quadrupole, so the reused-nuisance information is candidate-independent.
            trial_block_count = 1 + 2 * step
            inverse_nuisance = nuisance_inverse(
                whitened_nuisance, trial_block_count, args.nuisance_rms_m
            )
            trials: list[tuple[float, float, str, np.ndarray, np.ndarray, float, float]] = []
            for candidate in candidates:
                if candidate in selected:
                    continue
                trial_names = selected + [candidate]
                target_blocks = [nominal]
                for name in trial_names:
                    index = candidate_index[name]
                    target_blocks.extend((plus[index], minus[index]))
                information = marginalized_information(
                    target_blocks, whitened_nuisance, sigma, inverse_nuisance
                )
                sigmas, worst_axis, worst_direction = covariance_metrics(information)
                logdet = positive_logdet(information, f"{target}/{'+'.join(trial_names)}")
                precision = baseline_worst_axis / worst_axis
                worst_direction_precision = baseline_worst_direction / worst_direction
                trials.append(
                    (
                        logdet,
                        precision,
                        candidate,
                        information,
                        sigmas,
                        worst_axis,
                        worst_direction_precision,
                    )
                )

            # Information gain is primary; worst-axis precision is a deterministic
            # tie-breaker and remains explicit for the later operational guardrail.
            (
                chosen_logdet,
                chosen_precision,
                chosen,
                _chosen_information,
                chosen_sigmas,
                _chosen_worst_axis,
                chosen_worst_direction_precision,
            ) = max(trials, key=lambda trial: (trial[0], trial[1], trial[2]))
            selected.append(chosen)
            single = score_lookup[(target, chosen)]
            row = {
                "sextupole": target,
                "sextupole_s_m": sextupole_s[target],
                "greedy_step": step,
                "quadrupole": chosen,
                "provisional_operational": int(step <= args.operational_count),
                "cumulative_information_gain_logdet": chosen_logdet - baseline_logdet,
                "incremental_information_gain_logdet": chosen_logdet - previous_logdet,
                "cumulative_precision_improvement_worst_axis": chosen_precision,
                "cumulative_precision_improvement_worst_direction": chosen_worst_direction_precision,
                "joint_sigma_x_um": 1.0e6 * float(chosen_sigmas[0]),
                "joint_sigma_y_um": 1.0e6 * float(chosen_sigmas[1]),
                "single_information_gain_logdet": float(single["information_gain_logdet"]),
                "single_precision_improvement_worst_axis": float(
                    single["precision_improvement_worst_axis"]
                ),
                "single_screen_rank": int(single["screen_rank"]),
                "single_optics_leverage": float(single["optics_leverage"]),
            }
            target_rows.append(row)
            selection_rows.append(row)
            previous_logdet = chosen_logdet

        by_step = {int(row["greedy_step"]): row for row in target_rows}
        summary: dict[str, Any] = {
            "sextupole": target,
            "sextupole_s_m": sextupole_s[target],
            "baseline_sigma_x_um": 1.0e6 * float(baseline_sigma[0]),
            "baseline_sigma_y_um": 1.0e6 * float(baseline_sigma[1]),
        }
        for step, quadrupole in enumerate(selected, start=1):
            summary[f"candidate_{step}"] = quadrupole
        operational = by_step[args.operational_count]
        retained = by_step[args.retain_count]
        summary.update(
            {
                "operational_count": args.operational_count,
                "operational_information_gain_logdet": operational[
                    "cumulative_information_gain_logdet"
                ],
                "operational_precision_improvement_worst_axis": operational[
                    "cumulative_precision_improvement_worst_axis"
                ],
                "retained_count": args.retain_count,
                "retained_information_gain_logdet": retained[
                    "cumulative_information_gain_logdet"
                ],
                "retained_precision_improvement_worst_axis": retained[
                    "cumulative_precision_improvement_worst_axis"
                ],
            }
        )
        summary_rows.append(summary)

    selection_rows.sort(key=lambda row: (float(row["sextupole_s_m"]), int(row["greedy_step"])))
    summary_rows.sort(key=lambda row: float(row["sextupole_s_m"]))
    write_csv(output_dir / "greedy_quadrupole_selection.csv", selection_rows)
    write_csv(output_dir / "quadrupole_sets_by_sextupole.csv", summary_rows)

    operational_information = np.asarray(
        [float(row["operational_information_gain_logdet"]) for row in summary_rows]
    )
    retained_information = np.asarray(
        [float(row["retained_information_gain_logdet"]) for row in summary_rows]
    )
    operational_precision = np.asarray(
        [float(row["operational_precision_improvement_worst_axis"]) for row in summary_rows]
    )
    retained_precision = np.asarray(
        [float(row["retained_precision_improvement_worst_axis"]) for row in summary_rows]
    )

    def distribution(values: np.ndarray) -> dict[str, float]:
        return {
            "minimum": float(np.min(values)),
            "median": float(np.median(values)),
            "maximum": float(np.max(values)),
        }

    operational_rows = [
        row for row in selection_rows if int(row["greedy_step"]) <= args.operational_count
    ]
    operational_counts = Counter(str(row["quadrupole"]) for row in operational_rows)
    summary = {
        "target_count": len(summary_rows),
        "selection_row_count": len(selection_rows),
        "operational_union_quadrupole_count": len(operational_counts),
        "retained_union_quadrupole_count": len(
            {str(row["quadrupole"]) for row in selection_rows}
        ),
        "operational_information_gain_logdet": distribution(operational_information),
        "retained_information_gain_logdet": distribution(retained_information),
        "operational_fraction_of_retained_information_gain": distribution(
            operational_information / retained_information
        ),
        "operational_precision_improvement_worst_axis": distribution(
            operational_precision
        ),
        "retained_precision_improvement_worst_axis": distribution(retained_precision),
        "most_frequent_operational_quadrupoles": [
            {"quadrupole": name, "target_count": count}
            for name, count in operational_counts.most_common(10)
        ],
    }
    (output_dir / "selection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    metadata = {
        "format": "cesr-sextupole-greedy-quadrupole-selection-v1",
        "engine": "SciBmad/GTPSA repaired-lattice response bundles",
        "response_directory": str(response_dir),
        "affinity_directory": str(affinity_dir),
        "target_count": len(summary_rows),
        "screened_candidates_per_target": len(score_rows) // len(summary_rows),
        "retained_candidates_per_target": args.retain_count,
        "provisional_operational_candidates_per_target": args.operational_count,
        "selection_objective": "greedy maximum cumulative nuisance-marginalized log determinant; worst-axis precision is the tie-breaker and retained diagnostic",
        "condition_model": "nominal plus symmetric one-at-a-time K1 +/- blocks for every selected quadrupole",
        "nuisance_rms_m": args.nuisance_rms_m,
        "k2_step_m3": args.k2_step_m3,
        "k2_levels": levels,
        "excluded_quadrupoles": sorted(excluded_quadrupoles),
        "limitations": [
            "This selection reuses nominal other-sextupole nuisance responses in every K1 block.",
            "It uses the nominal-launch trajectory dictionary and does not yet contain the exact 3 x 3 orbit-bump grid.",
            "The first three candidates are provisional until finite-amplitude, observable-ablation, and operational-knob checks are complete.",
        ],
    }
    (output_dir / "selection_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Targets: {len(summary_rows)}")
    print(f"Selections: {len(selection_rows)}")
    print(f"Retained per target: {args.retain_count}")
    print(f"Provisional operational per target: {args.operational_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
