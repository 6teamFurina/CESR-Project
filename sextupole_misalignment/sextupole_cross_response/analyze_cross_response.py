#!/usr/bin/env python3
"""Analyze locality and effective rank of the all-sextupole response matrices."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "results" / "raw"
DEFAULT_OUTPUT = HERE / "results" / "analysis"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def locality(values: np.ndarray, target: int) -> dict[str, float]:
    """Summarize one target's nonnegative radial response over observations."""
    values = np.asarray(values, dtype=float)
    energy = values**2
    total = float(np.sum(energy))
    diagonal = float(values[target])
    off = np.delete(values, target)
    sorted_energy = np.sort(energy)[::-1]
    cumulative = np.cumsum(sorted_energy) / total if total > 0.0 else np.zeros_like(energy)

    def count_for(fraction: float) -> int:
        return int(np.searchsorted(cumulative, fraction, side="left") + 1) if total > 0 else 0

    participation = total**2 / float(np.sum(energy**2)) if np.sum(energy**2) > 0 else 0.0
    return {
        "target_amplitude": diagonal,
        "off_target_rms": float(np.sqrt(np.mean(off**2))),
        "off_target_max": float(np.max(off)),
        "offdiag_l2_over_target": (
            float(np.linalg.norm(off) / diagonal) if diagonal > 0 else float("inf")
        ),
        "target_energy_fraction": float(energy[target] / total) if total > 0 else 0.0,
        "participation_count": participation,
        "observations_for_50pct_energy": count_for(0.50),
        "observations_for_90pct_energy": count_for(0.90),
        "observations_for_99pct_energy": count_for(0.99),
    }


def matrix_svd(name: str, matrix: np.ndarray) -> tuple[list[dict[str, object]], dict[str, object]]:
    singular = np.linalg.svd(matrix, compute_uv=False)
    squared = singular**2
    fraction = squared / np.sum(squared)
    cumulative = np.cumsum(fraction)
    positive = fraction[fraction > 0]
    effective_rank = float(np.exp(-np.sum(positive * np.log(positive))))
    tolerance = max(matrix.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.sum(singular > tolerance))

    def modes_for(level: float) -> int:
        return int(np.searchsorted(cumulative, level, side="left") + 1)

    rows = [
        {
            "matrix": name,
            "mode": index + 1,
            "singular_value": float(value),
            "energy_fraction": float(fraction[index]),
            "cumulative_energy_fraction": float(cumulative[index]),
        }
        for index, value in enumerate(singular)
    ]
    summary = {
        "matrix": name,
        "rows": int(matrix.shape[0]),
        "columns": int(matrix.shape[1]),
        "numerical_rank": rank,
        "effective_rank": effective_rank,
        "modes_for_90pct_energy": modes_for(0.90),
        "modes_for_99pct_energy": modes_for(0.99),
        "sigma_max": float(singular[0]),
        "sigma_min_retained": float(singular[rank - 1]) if rank else 0.0,
        "condition_number_retained": (
            float(singular[0] / singular[rank - 1]) if rank else float("inf")
        ),
    }
    return rows, summary


def normalized_heatmap(values: np.ndarray, names: list[str], title: str, path: Path) -> None:
    radial = np.asarray(values, dtype=float)
    scale = np.linalg.norm(radial, axis=1, keepdims=True)
    normalized = np.divide(radial, scale, out=np.zeros_like(radial), where=scale > 0)
    image = np.log10(np.maximum(normalized, 1.0e-5))
    fig, axis = plt.subplots(figsize=(10.5, 8.5), constrained_layout=True)
    display = axis.imshow(image, aspect="auto", origin="lower", cmap="viridis", vmin=-5, vmax=0)
    ticks = np.arange(0, len(names), 5)
    axis.set_xticks(ticks, [names[index] for index in ticks], rotation=90, fontsize=6)
    axis.set_yticks(ticks, [names[index] for index in ticks], fontsize=6)
    axis.set_xlabel("Observation sextupole")
    axis.set_ylabel("Excited sextupole")
    axis.set_title(title)
    colorbar = fig.colorbar(display, ax=axis)
    colorbar.set_label("log10(row-L2-normalized radial response)")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def aggregate_locality(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for response in sorted({str(row["response"]) for row in rows}):
        selected = [row for row in rows if row["response"] == response]
        item: dict[str, object] = {"response": response, "target_count": len(selected)}
        for key in (
            "target_amplitude",
            "off_target_rms",
            "off_target_max",
            "offdiag_l2_over_target",
            "target_energy_fraction",
            "participation_count",
            "observations_for_50pct_energy",
            "observations_for_90pct_energy",
            "observations_for_99pct_energy",
        ):
            values = np.asarray([float(row[key]) for row in selected])
            finite = values[np.isfinite(values)]
            item[f"{key}_median"] = float(np.median(finite)) if finite.size else float("inf")
            item[f"{key}_p90"] = float(np.percentile(finite, 90)) if finite.size else float("inf")
            item[f"{key}_max"] = float(np.max(finite)) if finite.size else float("inf")
        result.append(item)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source = args.input_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    inventory = read_rows(source / "target_inventory.csv")
    names = [row["target"] for row in inventory]
    nt = len(names)
    kick = np.load(source / "periodic_kick_response.npy")
    bump = np.load(source / "bump_response.npy")
    sext_source = np.load(source / "sextupole_source_response.npy")
    design = np.load(source / "alignment_design.npy")
    if kick.shape != (nt, nt, 2, 2) or bump.shape != kick.shape or sext_source.shape != kick.shape:
        raise ValueError("Unexpected propagation matrix shape")
    if design.shape != (nt, 2, nt, 2, 2):
        raise ValueError(f"Unexpected alignment-design shape: {design.shape}")

    response_fields: list[tuple[str, np.ndarray]] = []
    for axis, label in enumerate(("x", "y")):
        response_fields.append((f"bump_{label}", np.linalg.norm(bump[..., axis], axis=2)))
    for source_axis, label in enumerate(("normal", "skew")):
        response_fields.append(
            (
                f"source_{label}",
                np.linalg.norm(sext_source[..., source_axis], axis=2),
            )
        )
    for bump_axis, bump_label in enumerate(("bx", "by")):
        for center_axis, center_label in enumerate(("cx", "cy")):
            response_fields.append(
                (
                    f"alignment_{bump_label}_{center_label}",
                    np.linalg.norm(design[:, bump_axis, :, :, center_axis], axis=2),
                )
            )

    locality_rows: list[dict[str, object]] = []
    for response_name, radial in response_fields:
        for target in range(nt):
            locality_rows.append(
                {
                    "response": response_name,
                    "target_index": target + 1,
                    "target": names[target],
                    **locality(radial[target], target),
                }
            )
    write_rows(output / "per_target_locality.csv", locality_rows)
    aggregate = aggregate_locality(locality_rows)
    write_rows(output / "aggregate_locality.csv", aggregate)

    bump_matrix = np.transpose(bump, (1, 2, 0, 3)).reshape(2 * nt, 2 * nt)
    kick_matrix = np.transpose(kick, (1, 2, 0, 3)).reshape(2 * nt, 2 * nt)
    source_matrix = np.transpose(sext_source, (1, 2, 0, 3)).reshape(2 * nt, 2 * nt)
    shared_alignment_template = np.transpose(design, (1, 2, 3, 0, 4)).reshape(
        4 * nt, 2 * nt
    )
    np.save(output / "bump_matrix_152x152.npy", bump_matrix)
    np.save(output / "periodic_kick_matrix_152x152.npy", kick_matrix)
    np.save(output / "sextupole_source_matrix_152x152.npy", source_matrix)
    np.save(
        output / "shared_alignment_template_matrix_304x152.npy",
        shared_alignment_template,
    )

    target_design_rows: list[dict[str, object]] = []
    target_design_singular_values: list[np.ndarray] = []
    for target in range(nt):
        target_matrix = design[target].reshape(4 * nt, 2)
        singular = np.linalg.svd(target_matrix, compute_uv=False)
        target_design_singular_values.append(singular)
        tolerance = max(target_matrix.shape) * np.finfo(float).eps * singular[0]
        gram = target_matrix.T @ target_matrix
        column_cosine = float(gram[0, 1] / np.sqrt(gram[0, 0] * gram[1, 1]))
        target_design_rows.append(
            {
                "target_index": target + 1,
                "target": names[target],
                "rows": int(target_matrix.shape[0]),
                "columns": int(target_matrix.shape[1]),
                "numerical_rank": int(np.sum(singular > tolerance)),
                "sigma_1": float(singular[0]),
                "sigma_2": float(singular[1]),
                "condition_number": float(singular[0] / singular[1]),
                "column_cosine": column_cosine,
            }
        )
    write_rows(output / "per_target_design_svd.csv", target_design_rows)

    singular_rows: list[dict[str, object]] = []
    singular_summary: list[dict[str, object]] = []
    for name, matrix in (
        ("bump_152x152", bump_matrix),
        ("periodic_kick_152x152", kick_matrix),
        ("sextupole_source_152x152", source_matrix),
        ("shared_alignment_template_304x152", shared_alignment_template),
    ):
        rows, summary = matrix_svd(name, matrix)
        singular_rows.extend(rows)
        singular_summary.append(summary)

    # The physical joint design for separate target scans is block diagonal,
    # with one 304-by-2 block per target.  Its singular values are the union of
    # the small-block values; allocating a dense 23104-by-152 matrix would be
    # unnecessary.
    block_name = "separate_scan_block_design_23104x152"
    block_singular = np.sort(np.concatenate(target_design_singular_values))[::-1]
    block_squared = block_singular**2
    block_fraction = block_squared / np.sum(block_squared)
    block_cumulative = np.cumsum(block_fraction)
    block_tolerance = (4 * nt * nt) * np.finfo(float).eps * block_singular[0]
    block_rank = int(np.sum(block_singular > block_tolerance))
    block_positive = block_fraction[block_fraction > 0]
    block_effective_rank = float(
        np.exp(-np.sum(block_positive * np.log(block_positive)))
    )
    singular_rows.extend(
        {
            "matrix": block_name,
            "mode": index + 1,
            "singular_value": float(value),
            "energy_fraction": float(block_fraction[index]),
            "cumulative_energy_fraction": float(block_cumulative[index]),
        }
        for index, value in enumerate(block_singular)
    )
    singular_summary.append(
        {
            "matrix": block_name,
            "rows": 4 * nt * nt,
            "columns": 2 * nt,
            "numerical_rank": block_rank,
            "effective_rank": block_effective_rank,
            "modes_for_90pct_energy": int(
                np.searchsorted(block_cumulative, 0.90, side="left") + 1
            ),
            "modes_for_99pct_energy": int(
                np.searchsorted(block_cumulative, 0.99, side="left") + 1
            ),
            "sigma_max": float(block_singular[0]),
            "sigma_min_retained": float(block_singular[block_rank - 1]),
            "condition_number_retained": float(
                block_singular[0] / block_singular[block_rank - 1]
            ),
        }
    )
    write_rows(output / "singular_values.csv", singular_rows)
    write_rows(output / "svd_summary.csv", singular_summary)

    bump_combined = np.sqrt(np.sum(bump**2, axis=(2, 3)))
    source_combined = np.sqrt(np.sum(sext_source**2, axis=(2, 3)))
    alignment_combined = np.sqrt(np.sum(design**2, axis=(1, 3, 4)))
    normalized_heatmap(
        bump_combined,
        names,
        "Corrector-bump propagation at all sextupoles",
        output / "bump_cross_response_heatmap.png",
    )
    normalized_heatmap(
        source_combined,
        names,
        "Periodic local sextupole-source propagation",
        output / "source_cross_response_heatmap.png",
    )
    normalized_heatmap(
        alignment_combined,
        names,
        "K2-bump-center alignment-signal propagation",
        output / "alignment_cross_response_heatmap.png",
    )

    target_conditions = np.asarray(
        [float(row["condition_number"]) for row in target_design_rows]
    )
    summary_json = {
        "format": "cesr-sextupole-cross-response-analysis-v1",
        "target_count": nt,
        "aggregate_locality": aggregate,
        "svd": singular_summary,
        "per_target_design": {
            "minimum_rank": int(
                min(int(row["numerical_rank"]) for row in target_design_rows)
            ),
            "maximum_rank": int(
                max(int(row["numerical_rank"]) for row in target_design_rows)
            ),
            "condition_number_median": float(np.median(target_conditions)),
            "condition_number_maximum": float(np.max(target_conditions)),
        },
    }
    (output / "summary.json").write_text(
        json.dumps(summary_json, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    by_name = {row["response"]: row for row in aggregate}
    bump_x = by_name["bump_x"]
    bump_y = by_name["bump_y"]
    align_rows = [
        row for row in aggregate if str(row["response"]).startswith("alignment_")
    ]
    align_energy = np.asarray(
        [float(row["target_energy_fraction_median"]) for row in align_rows]
    )
    align_participation = np.asarray(
        [float(row["participation_count_median"]) for row in align_rows]
    )
    svd_by_name = {row["matrix"]: row for row in singular_summary}
    report = f"""# All-sextupole GTPSA cross-response result

The latest repaired SciBmad lattice was linearized with GTPSA at the nominal
RF-on closed orbit.  The calculation covers all {nt} active normal sextupoles
as both excitation and observation locations.

## Compact derivative construction

No descriptor containing all sextupole K2, bump, and offset parameters was
formed.  The response uses an order-1 GTPSA periodic local-kick map, one
first-derivative corrector calculation, and the exact local normal-sextupole
polynomial.  The resulting selected K2--bump--center derivative is saved as
`alignment_design.npy` with axes
`target, bump_axis, observation_sextupole, output_plane, center_axis`.

## Locality

- x-bump median target-only energy fraction:
  `{100.0 * float(bump_x['target_energy_fraction_median']):.3f}%`;
  median participation count:
  `{float(bump_x['participation_count_median']):.2f}` sextupoles.
- y-bump median target-only energy fraction:
  `{100.0 * float(bump_y['target_energy_fraction_median']):.3f}%`;
  median participation count:
  `{float(bump_y['participation_count_median']):.2f}` sextupoles.
- Per unit target bump, the median off-target radial-orbit RMS is
  `{float(bump_x['off_target_rms_median']):.3f}` for x commands and
  `{float(bump_y['off_target_rms_median']):.3f}` for y commands.  For a
  0.5 mm command these are `{0.5 * float(bump_x['off_target_rms_median']):.3f}`
  and `{0.5 * float(bump_y['off_target_rms_median']):.3f} mm`; the largest
  individual off-target responses in the complete matrix are
  `{0.5 * float(bump_x['off_target_max_max']):.3f}` and
  `{0.5 * float(bump_y['off_target_max_max']):.3f} mm`.
- Across the four K2--bump--center channels, the median target-only energy
  fraction ranges from `{100.0 * float(np.min(align_energy)):.3f}%` to
  `{100.0 * float(np.max(align_energy)):.3f}%`; the median participation count
  ranges from `{float(np.min(align_participation)):.2f}` to
  `{float(np.max(align_participation)):.2f}` sextupoles.

These matrices are not block-local under a target-only definition.  Small
absolute orbit at a distant location must not be interpreted as negligible
until the response is whitened by the intended measurement covariance.

## Effective rank

| matrix | numerical rank | effective rank | modes for 90% energy | modes for 99% energy |
|---|---:|---:|---:|---:|
"""
    for name in (
        "bump_152x152",
        "periodic_kick_152x152",
        "sextupole_source_152x152",
        "shared_alignment_template_304x152",
        "separate_scan_block_design_23104x152",
    ):
        row = svd_by_name[name]
        report += (
            f"| `{name}` | {int(row['numerical_rank'])} | "
            f"{float(row['effective_rank']):.3f} | "
            f"{int(row['modes_for_90pct_energy'])} | "
            f"{int(row['modes_for_99pct_energy'])} |\n"
        )
    report += f"""

The shared `304 x 152` template matrix puts every target/center template into
one common channel coordinate system.  Its compact spectrum is evidence for a
candidate shared observation basis, not for a 20-dimensional joint inverse.
The physical design for separate one-target-at-a-time scans retains the target
axis and is block diagonal: it has full column rank {block_rank}, effective
rank `{block_effective_rank:.3f}`, and unwhitened condition number
`{block_singular[0] / block_singular[-1]:.3f}`.  Every individual `304 x 2`
target block has rank two; its condition-number median / maximum is
`{np.median(target_conditions):.12f} / {np.max(target_conditions):.12f}`.
These Euclidean nominal values are structural checks, not experimental
position precision.  A truncated shared basis must be tested against
covariance-whitened center-recovery error rather than selected by response
energy alone.

## Interpretation boundary

The source factorization is a nominal, first-order periodic-response result.
It uses an integrated thin normal-sextupole source placed at the element entry.
Finite-length effects, finite bump/K2 amplitudes, misaligned-background
relinearization, BPM noise, and machine-operating limits are outside this
matrix and require separate exact SciBmad validation.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
