#!/usr/bin/env python3
"""Validate exhaustive three-of-five quadrupole selection outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path

from analyze_affinity import LATEST_RESULTS


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--selection-dir", type=Path, default=LATEST_RESULTS / "selection")
    return result


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    args = parser().parse_args()
    selection_dir = args.selection_dir.expanduser().resolve()
    triplets = rows(selection_dir / "triplet_combinations.csv")
    best = rows(selection_dir / "best_triplets_by_sextupole.csv")
    summary = json.loads(
        (selection_dir / "triplet_selection_summary.json").read_text(encoding="utf-8")
    )
    metadata = json.loads(
        (selection_dir / "triplet_selection_metadata.json").read_text(encoding="utf-8")
    )

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in triplets:
        grouped[row["sextupole"]].append(row)
    best_by_target = {row["sextupole"]: row for row in best}

    assert len(grouped) == 76
    assert len(best) == 76
    assert len(triplets) == 76 * 10
    assert summary["target_count"] == 76
    assert summary["available_k1_condition_count"] == 11
    assert summary["triplet_count_per_target"] == 10
    assert metadata["evaluated_condition_count_per_triplet"] == 7

    for target, group in grouped.items():
        summary_row = best_by_target[target]
        retained = [summary_row[f"candidate_{index}"] for index in range(1, 6)]
        expected = {frozenset(item) for item in combinations(retained, 3)}
        actual = {
            frozenset((row["quadrupole_1"], row["quadrupole_2"], row["quadrupole_3"]))
            for row in group
        }
        assert expected == actual
        assert sorted(int(row["information_rank"]) for row in group) == list(range(1, 11))
        assert sorted(int(row["precision_rank"]) for row in group) == list(range(1, 11))
        selected = [row for row in group if int(row["selected"]) == 1]
        assert len(selected) == 1
        selected_row = selected[0]
        assert int(selected_row["information_rank"]) == 1
        selected_names = {
            selected_row["quadrupole_1"],
            selected_row["quadrupole_2"],
            selected_row["quadrupole_3"],
        }
        summary_names = {
            summary_row["selected_quadrupole_1"],
            summary_row["selected_quadrupole_2"],
            summary_row["selected_quadrupole_3"],
        }
        assert selected_names == summary_names
        values = [
            float(row[key])
            for row in group
            for key in (
                "information_gain_logdet",
                "precision_improvement_worst_axis",
                "precision_improvement_worst_direction",
                "sigma_x_um",
                "sigma_y_um",
            )
        ]
        assert all(math.isfinite(value) and value > 0.0 for value in values)
        assert float(summary_row["information_gain_over_greedy_prefix"]) >= -1.0e-12

    print(
        json.dumps(
            {
                "targets": len(grouped),
                "triplets": len(triplets),
                "triplets_per_target": 10,
                "available_k1_conditions": 11,
                "changed_from_greedy_prefix": summary[
                    "targets_changed_from_greedy_prefix"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
