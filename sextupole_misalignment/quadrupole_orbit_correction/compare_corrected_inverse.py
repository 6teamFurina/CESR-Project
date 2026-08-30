#!/usr/bin/env python3
"""Compare zero-offset, uncorrected-offset, and corrected-offset inverses."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SEQUENTIAL = HERE.parent / "sequential_joint_inverse"
WITHOUT = "without_quadrupole_misalignment"
UNCORRECTED = "with_quadrupole_misalignment"
CORRECTED = "with_quadrupole_misalignment_corrected"
MODELS = (
    "physics_gls",
    "shared_target_local_ridge",
    "shared_joint_ridge",
    "shared_joint_random_feature",
)
LEARNED_MODELS = MODELS[1:]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def same_distribution(
    rows: list[dict[str, str]], case: str, model: str
) -> dict[str, str]:
    candidates = [
        row
        for row in rows
        if row["evaluation_case"] == case
        and row["model"] == model
        and (
            row["training_case"] == "fixed_nominal_physics"
            or row["training_case"] == case
        )
    ]
    if len(candidates) != 1:
        raise ValueError(f"Expected one same-distribution row for {case}/{model}")
    return candidates[0]


def metric(row: dict[str, str], name: str) -> float:
    value = float(row[name])
    if not np.isfinite(value):
        raise ValueError(f"Non-finite {name}")
    return value


def best_learned(
    rows: list[dict[str, str]], case: str
) -> dict[str, str]:
    return min(
        (same_distribution(rows, case, model) for model in LEARNED_MODELS),
        key=lambda row: metric(row, "rmse_2d_um"),
    )


def ood_joint(rows: list[dict[str, str]], evaluation: str) -> dict[str, str]:
    candidates = [
        row
        for row in rows
        if row["training_case"] == WITHOUT
        and row["evaluation_case"] == evaluation
        and row["model"] == "shared_joint_ridge"
    ]
    if len(candidates) != 1:
        raise ValueError(f"Expected one zero-trained joint row for {evaluation}")
    return candidates[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--uncorrected-analysis",
        type=Path,
        default=SEQUENTIAL / "results" / "joint_inverse_analysis",
    )
    parser.add_argument(
        "--corrected-analysis",
        type=Path,
        default=SEQUENTIAL / "results" / "joint_inverse_analysis_corrected",
    )
    args = parser.parse_args()
    uncorrected_dir = args.uncorrected_analysis.resolve()
    corrected_dir = args.corrected_analysis.resolve()
    uncorrected_rows = read_rows(uncorrected_dir / "summary.csv")
    corrected_rows = read_rows(corrected_dir / "summary.csv")

    comparison_rows: list[dict[str, object]] = []
    for model in MODELS:
        zero = same_distribution(corrected_rows, WITHOUT, model)
        uncorrected = same_distribution(uncorrected_rows, UNCORRECTED, model)
        corrected = same_distribution(corrected_rows, CORRECTED, model)
        zero_rmse = metric(zero, "rmse_2d_um")
        uncorrected_rmse = metric(uncorrected, "rmse_2d_um")
        corrected_rmse = metric(corrected, "rmse_2d_um")
        comparison_rows.append(
            {
                "model": model,
                "zero_offset_rmse_2d_um": zero_rmse,
                "uncorrected_offset_rmse_2d_um": uncorrected_rmse,
                "corrected_offset_rmse_2d_um": corrected_rmse,
                "corrected_rmse_reduction_percent_vs_uncorrected": 100.0
                * (1.0 - corrected_rmse / uncorrected_rmse),
                "corrected_rmse_excess_percent_vs_zero_offset": 100.0
                * (corrected_rmse / zero_rmse - 1.0),
                "zero_offset_p99_2d_um": metric(zero, "p99_2d_um"),
                "uncorrected_offset_p99_2d_um": metric(uncorrected, "p99_2d_um"),
                "corrected_offset_p99_2d_um": metric(corrected, "p99_2d_um"),
                "zero_offset_worst_target_rmse_2d_um": metric(
                    zero, "worst_target_rmse_2d_um"
                ),
                "uncorrected_offset_worst_target_rmse_2d_um": metric(
                    uncorrected, "worst_target_rmse_2d_um"
                ),
                "corrected_offset_worst_target_rmse_2d_um": metric(
                    corrected, "worst_target_rmse_2d_um"
                ),
            }
        )
    write_rows(corrected_dir / "protocol_comparison.csv", comparison_rows)

    best_zero = best_learned(corrected_rows, WITHOUT)
    best_uncorrected = best_learned(uncorrected_rows, UNCORRECTED)
    best_corrected = best_learned(corrected_rows, CORRECTED)
    zero_rmse = metric(best_zero, "rmse_2d_um")
    uncorrected_rmse = metric(best_uncorrected, "rmse_2d_um")
    corrected_rmse = metric(best_corrected, "rmse_2d_um")
    removed_excess = 100.0 * (
        1.0 - (corrected_rmse - zero_rmse) / (uncorrected_rmse - zero_rmse)
    )
    uncorrected_ood = metric(ood_joint(uncorrected_rows, UNCORRECTED), "rmse_2d_um")
    corrected_ood = metric(ood_joint(corrected_rows, CORRECTED), "rmse_2d_um")
    diagnostics = json.loads(
        (corrected_dir / "diagnostics.json").read_text(encoding="utf-8")
    )
    summary = {
        "format": "cesr-corrected-sequential-inverse-comparison-v1",
        "zero_offset_best_model": best_zero["model"],
        "zero_offset_best_rmse_2d_um": zero_rmse,
        "uncorrected_offset_best_model": best_uncorrected["model"],
        "uncorrected_offset_best_rmse_2d_um": uncorrected_rmse,
        "corrected_offset_best_model": best_corrected["model"],
        "corrected_offset_best_rmse_2d_um": corrected_rmse,
        "corrected_best_rmse_reduction_percent_vs_uncorrected": 100.0
        * (1.0 - corrected_rmse / uncorrected_rmse),
        "corrected_best_rmse_excess_percent_vs_zero_offset": 100.0
        * (corrected_rmse / zero_rmse - 1.0),
        "removed_excess_rmse_percent": removed_excess,
        "zero_trained_joint_eval_uncorrected_rmse_2d_um": uncorrected_ood,
        "zero_trained_joint_eval_corrected_rmse_2d_um": corrected_ood,
        "precorrection_bpm_change_rms_um": diagnostics[
            "precorrection_reference_bpm_change_rms_um"
        ],
        "postcorrection_bpm_change_rms_um": diagnostics[
            "paired_reference_bpm_change_rms_um"
        ],
        "precorrection_target_change_rms_2d_um": diagnostics[
            "precorrection_beam_relative_truth_change_rms_2d_um"
        ],
        "postcorrection_target_change_rms_2d_um": diagnostics[
            "paired_beam_relative_truth_change_rms_2d_um"
        ],
        "postcorrection_outside_bump_radius_fraction": diagnostics[
            "with_truth_fraction_outside_bump_radius"
        ],
    }
    (corrected_dir / "protocol_comparison.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    table_rows = []
    for row in comparison_rows:
        table_rows.append(
            f"| {row['model']} | {row['zero_offset_rmse_2d_um']:.3f} | "
            f"{row['uncorrected_offset_rmse_2d_um']:.3f} | "
            f"{row['corrected_offset_rmse_2d_um']:.3f} | "
            f"{row['corrected_rmse_reduction_percent_vs_uncorrected']:.2f}% | "
            f"{row['corrected_rmse_excess_percent_vs_zero_offset']:+.2f}% |"
        )
    report = f"""# Fixed-baseline orbit-correction protocol comparison

All rows use the same deterministic 16-machine SciBmad ensemble, 10/3/3
machine split, measurement seeds, noise/drift augmentation, and inverse
definitions.  The only physical protocol change between the two quadrupole-
offset columns is whether the BPM-reference baseline correction is applied and
held fixed during every sextupole scan.

| inverse | zero-offset RMSE [um] | uncorrected offset RMSE [um] | corrected offset RMSE [um] | corrected reduction | corrected excess vs zero |
|---|---:|---:|---:|---:|---:|
{chr(10).join(table_rows)}

The best learned uncorrected result is
`{best_uncorrected['model']}` at `{uncorrected_rmse:.3f} um`.  After fixed
baseline correction, the best learned result is
`{best_corrected['model']}` at `{corrected_rmse:.3f} um`, a
`{summary['corrected_best_rmse_reduction_percent_vs_uncorrected']:.2f}%`
reduction.  It is only
`{summary['corrected_best_rmse_excess_percent_vs_zero_offset']:.2f}%` above the
best zero-offset result and removes `{removed_excess:.2f}%` of the excess RMSE
attributed to the uncorrected quadrupole-drift protocol.

The zero-offset-trained joint ridge gives `{uncorrected_ood:.3f} um` when
evaluated on uncorrected drift but `{corrected_ood:.3f} um` after the baseline
correction.  Thus correction removes most of the operational distribution
shift before the inverse is asked to estimate sextupole centers.

The corrected result still fails the strict tail gate: its best learned P99 is
`{metric(best_corrected, 'p99_2d_um'):.3f} um` and worst-target RMSE is
`{metric(best_corrected, 'worst_target_rmse_2d_um'):.3f} um`, both above 50
micrometers.  The result establishes the workflow benefit, not final CESR
precision or hardware safety.
"""
    (corrected_dir / "CORRECTION_COMPARISON.md").write_text(
        report, encoding="utf-8"
    )
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
