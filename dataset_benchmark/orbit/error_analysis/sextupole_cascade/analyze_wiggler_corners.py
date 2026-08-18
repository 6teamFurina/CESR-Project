#!/usr/bin/env python3

"""Vector inclusion-exclusion analysis of K2/wiggler ablation corners."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def norm(values):
    return math.sqrt(sum(value * value for value in values))


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def percentile(values, probability):
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vectors", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir or args.vectors.parent
    data = defaultdict(list)
    with args.vectors.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            key = (float(row["lambda2"]), float(row["wiggler_scale"]), int(row["trial"]))
            data[key].append(float(row["c3_y_m"]))
    trials = sorted({key[2] for key in data})
    rows = []
    concatenated = {label: [] for label in ("residual", "sextupole", "wiggler", "interaction", "baseline")}
    for trial in trials:
        c00 = data[(0.0, 0.0, trial)]
        c10 = data[(1.0, 0.0, trial)]
        c01 = data[(0.0, 1.0, trial)]
        c11 = data[(1.0, 1.0, trial)]
        components = {
            "residual": c00,
            "sextupole": [a - b for a, b in zip(c10, c00)],
            "wiggler": [a - b for a, b in zip(c01, c00)],
            "interaction": [a - b - c + d for a, b, c, d in zip(c11, c10, c01, c00)],
        }
        denominator = norm(c11)
        row = {"trial": trial}
        for label, component in components.items():
            row[f"{label}_norm_fraction"] = norm(component) / denominator
            row[f"{label}_signed_projection"] = dot(component, c11) / denominator**2
            concatenated[label].extend(component)
        concatenated["baseline"].extend(c11)
        rows.append(row)

    with (output_dir / "wiggler_corner_direction_attribution.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    baseline = concatenated["baseline"]
    baseline_norm = norm(baseline)
    lines = [
        "# K2/wiggler cubic-vector corner attribution",
        "",
        "The four-corner inclusion-exclusion decomposition is",
        "",
        "`C11 = C00 + (C10-C00) + (C01-C00) + (C11-C10-C01+C00)`.",
        "",
        "| Component | Global norm / nominal | Global signed projection | Direction median norm [P10, P90] |",
        "|---|---:|---:|---:|",
    ]
    for label in ("residual", "sextupole", "wiggler", "interaction"):
        component = concatenated[label]
        fractions = [row[f"{label}_norm_fraction"] for row in rows]
        lines.append(
            f"| {label} | {norm(component)/baseline_norm:.6f} | "
            f"{dot(component, baseline)/baseline_norm**2:.6f} | "
            f"{percentile(fractions, 0.5):.6f} "
            f"[{percentile(fractions, 0.1):.6f}, {percentile(fractions, 0.9):.6f}] |"
        )
    lines.extend(
        [
            "",
            "Norm fractions do not add because the signed detector vectors interfere. Signed projections add to one (up to numerical extraction error) and quantify alignment with the nominal response.",
        ]
    )
    (output_dir / "WIGGLER_CORNER_RESULTS.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
