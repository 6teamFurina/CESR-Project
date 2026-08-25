#!/usr/bin/env python3
"""Compare paired exact finite-amplitude scans with the GTPSA alignment design."""

from __future__ import annotations

import argparse
import csv
import json
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, values: list[dict[str, object]]) -> None:
    if not values:
        raise ValueError("Cannot write an empty validation table")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=HERE / "results" / "raw")
    parser.add_argument(
        "--exact-dir", type=Path, default=HERE / "results" / "exact_validation"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=HERE / "results" / "exact_validation_analysis"
    )
    args = parser.parse_args()
    raw = args.raw_dir.resolve()
    exact = args.exact_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with (exact / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    orbits = np.load(exact / "exact_sextupole_orbits.npy")
    centers = np.load(exact / "scenario_centers.npy")
    design = np.load(raw / "alignment_design.npy")
    selected = rows(exact / "selected_targets.csv")
    states = rows(exact / "states.csv")
    nt = len(selected)
    if orbits.shape != (nt, 2, 9, 76, 2) or centers.shape != (nt, 2, 2):
        raise ValueError(f"Unexpected exact-validation arrays: {orbits.shape}, {centers.shape}")

    state_index = {row["state"]: int(row["state_index"]) - 1 for row in states}
    bump = float(metadata["bump_amplitude_m"])
    kstep = float(metadata["k2_step_m3"])
    gradients = np.zeros((nt, 2, 2, 76, 2))
    for target in range(nt):
        for scenario in range(2):
            for axis, label in enumerate(("x", "y")):
                minus_b = (
                    orbits[target, scenario, state_index[f"{label}_b-1_k1"]]
                    - orbits[target, scenario, state_index[f"{label}_b-1_k-1"]]
                ) / (2.0 * kstep)
                plus_b = (
                    orbits[target, scenario, state_index[f"{label}_b1_k1"]]
                    - orbits[target, scenario, state_index[f"{label}_b1_k-1"]]
                ) / (2.0 * kstep)
                gradients[target, scenario, axis] = (plus_b - minus_b) / (2.0 * bump)

    exact_increment = gradients[:, 1] - gradients[:, 0]
    predicted = np.zeros_like(exact_increment)
    rows_out: list[dict[str, object]] = []
    zero = state_index["zero"]
    for selected_index, row in enumerate(selected):
        inventory_index = int(row["inventory_index"]) - 1
        local_orbit = orbits[selected_index, :, zero, inventory_index, :]
        relative_center = centers[selected_index] - local_orbit
        center_increment = relative_center[1] - relative_center[0]
        predicted[selected_index] = np.einsum(
            "aopc,c->aop", design[inventory_index], center_increment
        )
        for axis, axis_name in enumerate(("x", "y")):
            residual = predicted[selected_index, axis] - exact_increment[selected_index, axis]
            exact_norm = float(np.linalg.norm(exact_increment[selected_index, axis]))
            residual_norm = float(np.linalg.norm(residual))
            # This is the residual in the raw orbit difference between the
            # +delta-K2 and -delta-K2 states at one signed bump.  The complete
            # four-corner odd/odd contrast is exactly twice this value.
            k2_pair_difference = (
                np.linalg.norm(residual, axis=-1) * bump * (2.0 * kstep) * 1e9
            )
            rows_out.append(
                {
                    "target": row["target"],
                    "bump_axis": axis_name,
                    "center_increment_x_um": float(center_increment[0] * 1e6),
                    "center_increment_y_um": float(center_increment[1] * 1e6),
                    "exact_gradient_l2": exact_norm,
                    "residual_gradient_l2": residual_norm,
                    "relative_l2_residual": residual_norm / exact_norm,
                    "k2_pair_orbit_difference_rmse_nm": float(
                        np.sqrt(np.mean(k2_pair_difference**2))
                    ),
                    "k2_pair_orbit_difference_max_nm": float(
                        np.max(k2_pair_difference)
                    ),
                }
            )

    np.save(output / "exact_alignment_gradient_increment.npy", exact_increment)
    np.save(output / "predicted_alignment_gradient_increment.npy", predicted)
    write_rows(output / "per_target_validation.csv", rows_out)
    exact_flat = exact_increment.reshape(-1)
    predicted_flat = predicted.reshape(-1)
    residual_flat = (predicted - exact_increment).reshape(-1)
    aggregate_relative = float(np.linalg.norm(residual_flat) / np.linalg.norm(exact_flat))
    cosine_similarity = float(
        np.dot(predicted_flat, exact_flat)
        / (np.linalg.norm(predicted_flat) * np.linalg.norm(exact_flat))
    )
    fitted_predicted_scale = float(
        np.dot(predicted_flat, exact_flat) / np.dot(predicted_flat, predicted_flat)
    )
    k2_pair_difference_all = (
        np.linalg.norm(predicted - exact_increment, axis=-1) * bump * (2.0 * kstep) * 1e9
    )
    summary = {
        "format": "cesr-sextupole-cross-response-exact-validation-analysis-v1",
        "target_count": nt,
        "channel_block_count": len(rows_out),
        "aggregate_relative_l2_residual": aggregate_relative,
        "cosine_similarity": cosine_similarity,
        "fitted_scale_multiplying_prediction": fitted_predicted_scale,
        "k2_pair_orbit_difference_rmse_nm": float(
            np.sqrt(np.mean(k2_pair_difference_all**2))
        ),
        "k2_pair_orbit_difference_p90_nm": float(
            np.percentile(k2_pair_difference_all, 90)
        ),
        "k2_pair_orbit_difference_max_nm": float(np.max(k2_pair_difference_all)),
        "four_corner_contrast_rmse_nm": float(
            2.0 * np.sqrt(np.mean(k2_pair_difference_all**2))
        ),
        "maximum_block_relative_l2_residual": float(
            max(float(row["relative_l2_residual"]) for row in rows_out)
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    figure, axis = plt.subplots(figsize=(6.2, 6.0), constrained_layout=True)
    axis.scatter(exact_flat, predicted_flat, s=8, alpha=0.35, linewidths=0)
    limit = 1.04 * max(float(np.max(np.abs(exact_flat))), float(np.max(np.abs(predicted_flat))))
    axis.plot([-limit, limit], [-limit, limit], color="black", linewidth=1.0, linestyle="--")
    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Exact paired-scan gradient increment [m^3]")
    axis.set_ylabel("GTPSA-factorized prediction [m^3]")
    axis.set_title("Selected finite-amplitude cross-response validation")
    axis.text(
        0.03,
        0.97,
        f"relative L2 residual = {aggregate_relative:.3%}\n"
        f"cosine similarity = {cosine_similarity:.6f}",
        transform=axis.transAxes,
        va="top",
    )
    figure.savefig(output / "exact_vs_predicted.png", dpi=180)
    plt.close(figure)

    report = f"""# Exact finite-amplitude cross-response validation

The paired exact SciBmad check uses {nt} selected targets, aligned and
misaligned machine states, signed {1e3 * bump:.3f} mm bumps, and
`delta K2 = +/-{kstep:.3f} m^-3`.  The aligned K2-odd/bump-odd gradient is
subtracted from the matched misaligned gradient before comparison with the
nominal GTPSA alignment design.

- aggregate relative L2 residual: `{aggregate_relative:.6e}`;
- cosine similarity: `{cosine_similarity:.9f}`;
- fitted scale multiplying the GTPSA prediction:
  `{fitted_predicted_scale:.9f}`;
- residual in the raw `+delta K2` minus `-delta K2` orbit difference at one
  signed bump, RMS:
  `{summary['k2_pair_orbit_difference_rmse_nm']:.6f} nm`;
- the same two-state residual, P90:
  `{summary['k2_pair_orbit_difference_p90_nm']:.6f} nm`;
- the same two-state residual, maximum:
  `{summary['k2_pair_orbit_difference_max_nm']:.6f} nm`;
- full four-corner odd/odd contrast residual RMS:
  `{summary['four_corner_contrast_rmse_nm']:.6f} nm`;
- maximum target/axis block relative L2 residual:
  `{summary['maximum_block_relative_l2_residual']:.6e}`.

This validates the compact source factorization only for the selected
finite-amplitude, single-target-offset cases.  It does not validate a
misaligned all-magnet background or measured-machine covariance.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
