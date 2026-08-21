#!/usr/bin/env python3
"""Compare nominal Bmad finite-difference responses with saved SciBmad GTPSA maps."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

HERE = Path(__file__).resolve().parent
ALIGNMENT_DIR = HERE.parent


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--response-dir", type=Path, default=HERE / "results" / "responses")
    result.add_argument(
        "--gtpsa-coefficients",
        type=Path,
        default=ALIGNMENT_DIR
        / "archived_methods"
        / "response_map"
        / "results"
        / "full"
        / "alignment_coefficients.csv",
    )
    result.add_argument("--output-dir", type=Path, default=HERE / "results" / "affinity")
    return result


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def observation_label(row: dict[str, str]) -> str:
    source = "ring" if row["observation_scope"] == "ring" else row["observation_name"]
    return f"{source}:{row['observable']}"


def main() -> int:
    args = parser().parse_args()
    response_dir = args.response_dir.expanduser().resolve()
    coefficient_path = args.gtpsa_coefficients.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    gtpsa: dict[str, dict[str, tuple[float, float]]] = {}
    with coefficient_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            gtpsa.setdefault(row["sextupole"], {})[observation_label(row)] = (
                float(row["d2_k2_x"]),
                float(row["d2_k2_y"]),
            )

    rows: list[dict[str, Any]] = []
    target_files = sorted((response_dir / "targets").glob("*_responses.npz"))
    for path in target_files:
        target = path.name.removesuffix("_responses.npz").upper()
        with np.load(path) as saved:
            labels = [str(value) for value in saved["observation_labels"]]
            bmad = np.asarray(saved["target_response_nominal"], dtype=float)
        reference = np.asarray([gtpsa[target][label] for label in labels], dtype=float)
        families = sorted({label.split(":", 1)[1] for label in labels})
        for family in families:
            indices = [index for index, label in enumerate(labels) if label.endswith(f":{family}")]
            calculated = bmad[indices].ravel()
            expected = reference[indices].ravel()
            denominator = float(np.linalg.norm(expected))
            calculated_norm = float(np.linalg.norm(calculated))
            dot_denominator = denominator * calculated_norm
            cosine = (
                float(np.dot(calculated, expected) / dot_denominator)
                if dot_denominator > 0.0
                else math.nan
            )
            rows.append(
                {
                    "sextupole": target,
                    "observable": family,
                    "value_count": len(calculated),
                    "cosine_similarity": cosine,
                    "norm_ratio_bmad_to_scibmad": (
                        calculated_norm / denominator if denominator > 0.0 else math.nan
                    ),
                    "relative_l2_difference": (
                        float(np.linalg.norm(calculated - expected) / denominator)
                        if denominator > 0.0
                        else math.nan
                    ),
                }
            )
    write_csv(output_dir / "nominal_bmad_vs_scibmad.csv", rows)
    finite_cosines = [float(row["cosine_similarity"]) for row in rows if math.isfinite(float(row["cosine_similarity"]))]
    finite_differences = [float(row["relative_l2_difference"]) for row in rows if math.isfinite(float(row["relative_l2_difference"]))]
    summary = {
        "target_count": len(target_files),
        "family_comparisons": len(rows),
        "minimum_cosine_similarity": min(finite_cosines),
        "maximum_relative_l2_difference": max(finite_differences),
        "interpretation": "Independent nominal Bmad finite differences versus saved SciBmad GTPSA mixed K2-offset derivatives.",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "nominal_response_validation.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
