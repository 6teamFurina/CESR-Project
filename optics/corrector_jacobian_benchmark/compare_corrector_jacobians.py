#!/usr/bin/env python3
"""Compare labeled Bmad and SciBmad corrector-optics Jacobian matrices."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--results-dir", type=Path, default=HERE / "results")
    return result


def read_matrix(path: Path) -> tuple[list[str], list[str], list[list[float]]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    columns = rows[0][1:]
    labels = [row[0] for row in rows[1:]]
    values = [[float(value) for value in row[1:]] for row in rows[1:]]
    if any(not math.isfinite(value) for row in values for value in row):
        raise RuntimeError(f"Non-finite value in {path}")
    return labels, columns, values


def matrix_metrics(
    reference: list[list[float]], candidate: list[list[float]]
) -> dict[str, float | None]:
    ref = [value for row in reference for value in row]
    got = [value for row in candidate for value in row]
    difference = [b - a for a, b in zip(ref, got)]
    ref_norm = math.sqrt(sum(value * value for value in ref))
    got_norm = math.sqrt(sum(value * value for value in got))
    diff_norm = math.sqrt(sum(value * value for value in difference))
    dot = sum(a * b for a, b in zip(ref, got))
    return {
        "relative_frobenius": (
            diff_norm / ref_norm if ref_norm else (0.0 if not got_norm else None)
        ),
        "cosine_correlation": (
            dot / (ref_norm * got_norm)
            if ref_norm and got_norm
            else (1.0 if not ref_norm and not got_norm else None)
        ),
        "maximum_absolute_difference": max(map(abs, difference)),
        "reference_maximum_absolute": max(map(abs, ref)),
    }


def select_rows(values: list[list[float]], indices: list[int]) -> list[list[float]]:
    return [values[index] for index in indices]


def phase_origin_aligned(
    labels: list[str], values: list[list[float]]
) -> list[list[float]]:
    result = [row[:] for row in values]
    for quantity in ("phi_1", "phi_2", "phi_3"):
        indices = [
            index for index, label in enumerate(labels)
            if label.rsplit(":", 1)[1] == quantity
        ]
        if not indices:
            continue
        origin = result[indices[0]][:]
        for index in indices:
            result[index] = [
                value - offset for value, offset in zip(result[index], origin)
            ]
    return result


def compare_pair(bmad_path: Path, scibmad_path: Path) -> dict[str, object]:
    bmad_labels, bmad_columns, bmad = read_matrix(bmad_path)
    sci_labels, sci_columns, sci = read_matrix(scibmad_path)
    if bmad_columns != sci_columns:
        raise RuntimeError(f"Control labels differ: {bmad_path.name}")
    if bmad_labels != sci_labels:
        # SciBmad's coasting tune vector has an additional slip row. Compare
        # the common labeled rows without discarding that extra saved result.
        sci_by_label = dict(zip(sci_labels, sci))
        if len(sci_by_label) != len(sci_labels) or any(
            label not in sci_by_label for label in bmad_labels
        ):
            raise RuntimeError(f"Observable labels differ: {bmad_path.name}")
        sci = [sci_by_label[label] for label in bmad_labels]
    result: dict[str, object] = {"overall": matrix_metrics(bmad, sci)}
    if "detector_optics" in bmad_path.name:
        quantities = [label.rsplit(":", 1)[1] for label in bmad_labels]
        result["by_quantity"] = {
            quantity: matrix_metrics(
                select_rows(bmad, [index for index, item in enumerate(quantities) if item == quantity]),
                select_rows(sci, [index for index, item in enumerate(quantities) if item == quantity]),
            )
            for quantity in dict.fromkeys(quantities)
        }
        phase_aligned_bmad = phase_origin_aligned(bmad_labels, bmad)
        phase_aligned_sci = phase_origin_aligned(bmad_labels, sci)
        result["phase_origin_aligned_overall"] = matrix_metrics(
            phase_aligned_bmad,
            phase_aligned_sci,
        )
        orbit_indices = [
            index for index, quantity in enumerate(quantities)
            if quantity in ("orbit_x", "orbit_px", "orbit_y", "orbit_py")
        ]
        result["transverse_detector_orbit"] = matrix_metrics(
            select_rows(bmad, orbit_indices),
            select_rows(sci, orbit_indices),
        )
    return result


def main() -> int:
    args = parser().parse_args()
    root = args.results_dir.resolve()
    matrices = {
        "detector_optics": "detector_optics_jacobian.csv",
        "closed_orbit_6x119": "closed_orbit_jacobian_6x119.csv",
        "ring_tunes": "ring_tune_jacobian.csv",
    }
    results = {
        name: compare_pair(
            root / "bmad" / f"bmad_{suffix}",
            root / "scibmad" / f"scibmad_{suffix}",
        )
        for name, suffix in matrices.items()
    }
    output = root / "comparison.json"
    output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for name, detail in results.items():
        metrics = detail["overall"]
        print(
            f"{name}: relative Frobenius={metrics['relative_frobenius']:.6e}, "
            f"correlation={metrics['cosine_correlation']:.12f}, "
            f"max abs difference={metrics['maximum_absolute_difference']:.6e}"
        )
    print(f"Comparison: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
