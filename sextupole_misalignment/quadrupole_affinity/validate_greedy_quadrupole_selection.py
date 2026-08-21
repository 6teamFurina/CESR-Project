#!/usr/bin/env python3
"""Validate the five-candidate, three-operational quadrupole selection."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from analyze_affinity import LATEST_RESULTS


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--affinity-dir", type=Path, default=LATEST_RESULTS / "affinity")
    result.add_argument("--selection-dir", type=Path, default=LATEST_RESULTS / "selection")
    return result


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    args = parser().parse_args()
    affinity_dir = args.affinity_dir.expanduser().resolve()
    selection_dir = args.selection_dir.expanduser().resolve()
    scores = rows(affinity_dir / "quadrupole_affinity_scores.csv")
    selected = rows(selection_dir / "greedy_quadrupole_selection.csv")
    sets = rows(selection_dir / "quadrupole_sets_by_sextupole.csv")
    metadata = json.loads((selection_dir / "selection_metadata.json").read_text(encoding="utf-8"))
    summary = json.loads((selection_dir / "selection_summary.json").read_text(encoding="utf-8"))

    score_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    selection_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in scores:
        score_groups[row["sextupole"]].append(row)
    for row in selected:
        selection_groups[row["sextupole"]].append(row)

    assert len(score_groups) == 76
    assert len(sets) == 76
    assert len(selected) == 76 * 5
    assert metadata["target_count"] == 76
    assert metadata["retained_candidates_per_target"] == 5
    assert metadata["provisional_operational_candidates_per_target"] == 3
    assert summary["target_count"] == 76
    assert summary["selection_row_count"] == 380

    for target, group in selection_groups.items():
        ordered = sorted(group, key=lambda row: int(row["greedy_step"]))
        assert [int(row["greedy_step"]) for row in ordered] == [1, 2, 3, 4, 5]
        assert len({row["quadrupole"] for row in ordered}) == 5
        assert [int(row["provisional_operational"]) for row in ordered] == [1, 1, 1, 0, 0]
        candidate_names = {row["quadrupole"] for row in score_groups[target]}
        assert all(row["quadrupole"] in candidate_names for row in ordered)

        information = [float(row["cumulative_information_gain_logdet"]) for row in ordered]
        precision = [
            float(row["cumulative_precision_improvement_worst_axis"]) for row in ordered
        ]
        assert all(math.isfinite(value) for value in information + precision)
        assert all(right >= left for left, right in zip(information, information[1:]))
        assert all(right >= left for left, right in zip(precision, precision[1:]))

        best_single = max(
            score_groups[target], key=lambda row: float(row["information_gain_logdet"])
        )
        assert ordered[0]["quadrupole"] == best_single["quadrupole"]

    print(
        json.dumps(
            {
                "targets": len(selection_groups),
                "selection_rows": len(selected),
                "retained_per_target": 5,
                "provisional_operational_per_target": 3,
                "operational_union_quadrupoles": summary[
                    "operational_union_quadrupole_count"
                ],
                "retained_union_quadrupoles": summary[
                    "retained_union_quadrupole_count"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
