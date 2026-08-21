#!/usr/bin/env python3
"""Compare each seven-condition triplet with seven repeated nominal blocks."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
AFFINITY = HERE.parent
sys.path.insert(0, str(AFFINITY))

from analyze_affinity import (  # noqa: E402
    covariance_metrics,
    load_target_bundle,
    marginalized_information,
    nuisance_inverse,
    slope_noise,
    target_bundles,
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--response-dir", type=Path, default=AFFINITY / "results" / "scibmad_latest" / "responses"
    )
    result.add_argument(
        "--candidate-csv",
        type=Path,
        default=AFFINITY
        / "results"
        / "scibmad_latest"
        / "selection"
        / "quadrupole_sets_by_sextupole.csv",
    )
    result.add_argument("--nuisance-rms-m", type=float, default=3.0e-4)
    result.add_argument("--k2-step-m3", type=float, default=0.01)
    result.add_argument("--output-dir", type=Path, default=HERE / "results" / "aggregate")
    return result


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def positive_logdet(matrix: np.ndarray) -> float:
    sign, value = np.linalg.slogdet(matrix)
    if sign <= 0:
        raise ValueError("Information matrix is not positive definite")
    return float(value)


def distribution(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "p90": float(np.percentile(array, 90.0)),
        "maximum": float(np.max(array)),
    }


def main() -> int:
    args = parser().parse_args()
    candidate_rows = rows(args.candidate_csv.expanduser().resolve())
    candidate_by_target = {
        row["sextupole"]: [row[f"candidate_{index}"] for index in range(1, 6)]
        for row in candidate_rows
    }
    output_rows: list[dict[str, Any]] = []
    levels = [-2.0, -1.0, 0.0, 1.0, 2.0]
    for target_path in target_bundles(args.response_dir.expanduser().resolve()):
        target, labels, nominal, nuisance, candidates, plus, minus = load_target_bundle(
            target_path, args.response_dir.expanduser().resolve()
        )
        retained = candidate_by_target[target]
        candidate_index = {name: index for index, name in enumerate(candidates)}
        sigma = slope_noise(labels, args.k2_step_m3, levels)
        whitened_nuisance = nuisance / sigma[:, None]
        inverse_nuisance = nuisance_inverse(
            whitened_nuisance, 7, args.nuisance_rms_m
        )
        repeated_information = marginalized_information(
            [nominal] * 7, whitened_nuisance, sigma, inverse_nuisance
        )
        repeated_logdet = positive_logdet(repeated_information)
        _, repeated_worst_axis, _ = covariance_metrics(repeated_information)
        for triplet in combinations(retained, 3):
            blocks = [nominal]
            for name in triplet:
                index = candidate_index[name]
                blocks.extend((plus[index], minus[index]))
            information = marginalized_information(
                blocks, whitened_nuisance, sigma, inverse_nuisance
            )
            _, worst_axis, _ = covariance_metrics(information)
            output_rows.append(
                {
                    "sextupole": target,
                    "quadrupole_1": triplet[0],
                    "quadrupole_2": triplet[1],
                    "quadrupole_3": triplet[2],
                    "matched_nominal_information_gain_logdet": positive_logdet(information)
                    - repeated_logdet,
                    "matched_nominal_precision_improvement_worst_axis": repeated_worst_axis
                    / worst_axis,
                }
            )
    if len(output_rows) != 760:
        raise RuntimeError(f"Expected 760 triplets, found {len(output_rows)}")
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "matched_nominal_response_triplets.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    information = [
        float(row["matched_nominal_information_gain_logdet"]) for row in output_rows
    ]
    precision = [
        float(row["matched_nominal_precision_improvement_worst_axis"])
        for row in output_rows
    ]
    best_by_target = []
    for target in sorted(candidate_by_target):
        target_rows = [row for row in output_rows if row["sextupole"] == target]
        best_by_target.append(
            max(
                target_rows,
                key=lambda row: float(
                    row["matched_nominal_information_gain_logdet"]
                ),
            )
        )
    summary = {
        "format": "cesr-repaired-lattice-matched-nominal-triplet-gain-v1",
        "target_count": 76,
        "triplet_count": 760,
        "comparison": "one nominal plus six K1-conditioned blocks versus seven repeated nominal blocks",
        "matched_nominal_information_gain_logdet": distribution(information),
        "matched_nominal_precision_improvement_worst_axis": distribution(precision),
        "best_triplet_matched_nominal_information_gain_logdet": distribution(
            [
                float(row["matched_nominal_information_gain_logdet"])
                for row in best_by_target
            ]
        ),
        "best_triplet_matched_nominal_precision_improvement_worst_axis": distribution(
            [
                float(row["matched_nominal_precision_improvement_worst_axis"])
                for row in best_by_target
            ]
        ),
    }
    (output_dir / "matched_nominal_response_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
