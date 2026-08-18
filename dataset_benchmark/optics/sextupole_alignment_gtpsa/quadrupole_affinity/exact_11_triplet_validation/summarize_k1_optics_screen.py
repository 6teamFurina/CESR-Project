#!/usr/bin/env python3
"""Summarize a SciBmad K1 optics screen and its impact on five-candidate pools."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--screen-csv", type=Path, required=True)
    result.add_argument("--candidate-csv", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    return result


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, values: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def main() -> int:
    args = parser().parse_args()
    screen_rows = rows(args.screen_csv.expanduser().resolve())
    representative: dict[str, dict[str, str]] = {}
    for row in screen_rows:
        representative.setdefault(row["quadrupole"], row)
    if len(representative) != 113:
        raise RuntimeError(f"Expected 113 quadrupoles, found {len(representative)}")
    unsafe = [
        {
            "quadrupole": name,
            "max_abs_tune_shift": float(row["max_abs_tune_shift"]),
            "max_detector_beta_beating": float(row["max_detector_beta_beating"]),
        }
        for name, row in representative.items()
        if int(row["allowed"]) == 0
    ]
    unsafe.sort(key=lambda row: str(row["quadrupole"]))
    unsafe_names = {str(row["quadrupole"]) for row in unsafe}

    impact_rows: list[dict[str, Any]] = []
    for row in rows(args.candidate_csv.expanduser().resolve()):
        candidates = [row[f"candidate_{index}"] for index in range(1, 6)]
        rejected = [name for name in candidates if name in unsafe_names]
        impact_rows.append(
            {
                "sextupole": row["sextupole"],
                "safe_original_candidate_count": 5 - len(rejected),
                "unsafe_original_candidates": ";".join(rejected),
            }
        )
    safe_counts = Counter(
        int(row["safe_original_candidate_count"]) for row in impact_rows
    )
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "k1_optics_unsafe_quadrupoles.csv", unsafe)
    write_csv(output_dir / "k1_optics_original_pool_impact.csv", impact_rows)
    summary = {
        "format": "cesr-repaired-lattice-k1-optics-screen-summary-v1",
        "quadrupole_count": len(representative),
        "allowed_count": len(representative) - len(unsafe),
        "rejected_count": len(unsafe),
        "rejected_quadrupoles": unsafe,
        "original_five_candidate_safe_count_distribution": {
            str(count): safe_counts[count] for count in sorted(safe_counts)
        },
        "all_targets_retain_at_least_three_original_candidates": all(
            int(row["safe_original_candidate_count"]) >= 3 for row in impact_rows
        ),
    }
    (output_dir / "k1_optics_screen_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
