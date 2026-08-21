#!/usr/bin/env python3
"""Compare the labeled Bmad and SciBmad extended optics Jacobians."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from compare_corrector_jacobians import (
    matrix_metrics,
    phase_origin_aligned,
    read_matrix,
    select_rows,
)

HERE = Path(__file__).resolve().parent


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--case-dir",
        type=Path,
        default=HERE / "results" / "extended" / "correctors_quads",
    )
    return result


def select_columns(values: list[list[float]], indices: list[int]) -> list[list[float]]:
    return [[row[index] for index in indices] for row in values]


def family_name(label: str) -> str:
    prefix = label.split(":", 1)[0]
    return {"COR": "corrector", "K1": "quadrupole", "K2": "sextupole"}[prefix]


def compare_pair(bmad_path: Path, scibmad_path: Path) -> dict[str, object]:
    bmad_labels, bmad_columns, bmad = read_matrix(bmad_path)
    sci_labels, sci_columns, sci = read_matrix(scibmad_path)
    if bmad_columns != sci_columns:
        raise RuntimeError(f"Parameter labels differ: {bmad_path.name}")
    if bmad_labels != sci_labels:
        sci_by_label = dict(zip(sci_labels, sci))
        if len(sci_by_label) != len(sci_labels) or any(
            label not in sci_by_label for label in bmad_labels
        ):
            raise RuntimeError(f"Observable labels differ: {bmad_path.name}")
        sci = [sci_by_label[label] for label in bmad_labels]

    family_indices: dict[str, list[int]] = {}
    for index, label in enumerate(bmad_columns):
        family_indices.setdefault(family_name(label), []).append(index)

    result: dict[str, object] = {
        "shape": [len(bmad_labels), len(bmad_columns)],
        "overall": matrix_metrics(bmad, sci),
        "by_parameter_family": {
            family: matrix_metrics(
                select_columns(bmad, indices), select_columns(sci, indices)
            )
            for family, indices in family_indices.items()
        },
    }

    if "detector_optics" in bmad_path.name:
        quantities = [label.rsplit(":", 1)[1] for label in bmad_labels]
        quantity_indices = {
            quantity: [
                index for index, item in enumerate(quantities) if item == quantity
            ]
            for quantity in dict.fromkeys(quantities)
        }
        result["by_quantity"] = {
            quantity: matrix_metrics(
                select_rows(bmad, indices), select_rows(sci, indices)
            )
            for quantity, indices in quantity_indices.items()
        }
        result["by_quantity_and_parameter_family"] = {
            quantity: {
                family: matrix_metrics(
                    select_columns(select_rows(bmad, row_indices), column_indices),
                    select_columns(select_rows(sci, row_indices), column_indices),
                )
                for family, column_indices in family_indices.items()
            }
            for quantity, row_indices in quantity_indices.items()
        }
        aligned_bmad = phase_origin_aligned(bmad_labels, bmad)
        aligned_sci = phase_origin_aligned(bmad_labels, sci)
        result["phase_origin_aligned_overall"] = matrix_metrics(
            aligned_bmad, aligned_sci
        )
        orbit_indices = [
            index
            for index, quantity in enumerate(quantities)
            if quantity in ("orbit_x", "orbit_px", "orbit_y", "orbit_py")
        ]
        result["transverse_detector_orbit"] = matrix_metrics(
            select_rows(bmad, orbit_indices), select_rows(sci, orbit_indices)
        )
        result["transverse_detector_orbit_by_parameter_family"] = {
            family: matrix_metrics(
                select_columns(select_rows(bmad, orbit_indices), indices),
                select_columns(select_rows(sci, orbit_indices), indices),
            )
            for family, indices in family_indices.items()
        }
    return result


def main() -> int:
    root = parser().parse_args().case_dir.resolve()
    matrices = {
        "detector_optics": "detector_optics_jacobian.csv",
        "closed_orbit": "closed_orbit_jacobian.csv",
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
    output.write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    for matrix_name, detail in results.items():
        print(matrix_name)
        for family, metrics in detail["by_parameter_family"].items():
            print(
                f"  {family}: relative Frobenius="
                f"{metrics['relative_frobenius']:.6e}, correlation="
                f"{metrics['cosine_correlation']:.12f}"
            )
    print(f"Comparison: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
