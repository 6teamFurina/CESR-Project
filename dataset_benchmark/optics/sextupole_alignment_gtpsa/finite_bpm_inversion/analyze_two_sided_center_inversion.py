#!/usr/bin/env python3
"""Invert all sextupole beam-relative centers using two-sided BPM local orbit.

The fit consumes full-ring BPM K2 slopes and relative local coordinates already
predicted from the nearest upstream/downstream BPM pair. Exact target-local
orbit and target offset are loaded only after all 608 fits for evaluation.
"""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_command_space_finite_bpm import fit_center, k2_slope


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = (
    HERE.parent
    / "direct_observable_nuisance_ablation"
    / "results"
    / "all_76_orbit_protocol"
)
DEFAULT_PREDICTIONS = (
    HERE / "results" / "local_orbit_predictors" / "two_sided_transport_local_orbits.npy"
)
DEFAULT_ORACLE_FITS = (
    DEFAULT_SOURCE / "per_realization_fits.csv"
)


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(errors_m: np.ndarray) -> dict[str, float]:
    radial_um = np.linalg.norm(errors_m, axis=-1) * 1e6
    return {
        "x_rmse_um": float(np.sqrt(np.mean(errors_m[..., 0] ** 2)) * 1e6),
        "y_rmse_um": float(np.sqrt(np.mean(errors_m[..., 1] ** 2)) * 1e6),
        "rmse_2d_um": float(np.sqrt(np.mean(radial_um**2))),
        "median_2d_um": float(np.median(radial_um)),
        "p90_2d_um": float(np.percentile(radial_um, 90)),
        "p99_2d_um": float(np.percentile(radial_um, 99)),
        "max_2d_um": float(np.max(radial_um)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--predicted-local-orbits", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--oracle-fits", type=Path, default=DEFAULT_ORACLE_FITS)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=HERE / "results" / "two_sided_center_inversion",
    )
    args = parser.parse_args()
    source = args.input_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with (source / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta_k2 = levels * float(metadata["k2_step_m3"])
    nominal_candidates = np.flatnonzero(levels == 0.0)
    if nominal_candidates.size != 1:
        raise ValueError("Expected exactly one nominal K2 level")
    nominal_k2 = int(nominal_candidates[0])
    bump_rows: list[dict[str, str]]
    with (source / "bump_points.csv").open(newline="", encoding="utf-8") as stream:
        bump_rows = list(csv.DictReader(stream))
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    zero_candidates = np.flatnonzero(np.all(bump_commands == 0.0, axis=1))
    if zero_candidates.size != 1:
        raise ValueError("Expected exactly one zero bump")
    zero_bump = int(zero_candidates[0])

    predicted_local = np.load(args.predicted_local_orbits.resolve())
    bpm_orbits = np.load(source / "bpm_orbits.npy", mmap_mode="r")
    nt, nr, nb, nk, nd, planes = bpm_orbits.shape
    if (nt, nb, nk, nd, planes) != (
        len(target_names),
        len(bump_commands),
        len(levels),
        len(bpm_names),
        2,
    ):
        raise ValueError("BPM tensor and metadata dimensions do not agree")
    if predicted_local.shape != (nt, nr, nb, 2):
        raise ValueError(f"Unexpected predicted local-orbit shape: {predicted_local.shape}")
    if np.max(np.abs(predicted_local[:, :, zero_bump, :])) > 1e-15:
        raise ValueError("Predicted zero-bump relative orbit is not zero")

    # Machine-facing inverse: only BPM scan tensors and the already predicted
    # local coordinates are used. All 111 x/y BPM channels contribute to the
    # K2 slope. Truth is intentionally unavailable in this loop.
    estimates = np.zeros((nt, nr, 2), dtype=float)
    for target in range(nt):
        for realization in range(nr):
            scan = np.asarray(bpm_orbits[target, realization], dtype=float).reshape(
                nb, nk, 2 * nd
            )
            slopes = k2_slope(scan, delta_k2)
            estimates[target, realization] = fit_center(
                slopes, predicted_local[target, realization]
            )
        print(f"two-sided center inversion {target + 1}/{nt}: {target_names[target]}")
    np.save(output / "relative_center_estimates.npy", estimates)

    # Evaluation-only truth begins here.
    target_truth = np.load(source / "target_truth.npy")
    target_orbits = np.load(source / "target_orbits.npy", mmap_mode="r")
    zero_orbit = np.asarray(
        target_orbits[:, :, zero_bump, nominal_k2, :], dtype=float
    )
    relative_truth = target_truth - zero_orbit
    exact_relative_local = (
        np.asarray(target_orbits[:, :, :, nominal_k2, :], dtype=float)
        - zero_orbit[:, :, None, :]
    )
    local_errors = predicted_local - exact_relative_local
    nonzero_bumps = np.arange(nb) != zero_bump
    per_case_local_rmse_um = np.sqrt(
        np.mean(
            np.sum(local_errors[:, :, nonzero_bumps, :] ** 2, axis=-1),
            axis=-1,
        )
    ) * 1e6

    center_errors = estimates - relative_truth
    radial_center_um = np.linalg.norm(center_errors, axis=-1) * 1e6
    overall = summarize(center_errors)
    correlation = float(
        np.corrcoef(per_case_local_rmse_um.ravel(), radial_center_um.ravel())[0, 1]
    )

    with args.oracle_fits.resolve().open(newline="", encoding="utf-8") as stream:
        oracle_rows = list(csv.DictReader(stream))
    if len(oracle_rows) != nt * nr:
        raise ValueError("Oracle fit count does not match target/realization inventory")
    oracle_errors_um = np.zeros((nt, nr, 2), dtype=float)
    for row in oracle_rows:
        target = int(row["target_index"]) - 1
        realization = int(row["realization"]) - 1
        if row["target"] != target_names[target]:
            raise ValueError("Oracle target inventory does not match")
        oracle_errors_um[target, realization, 0] = (
            float(row["estimate_x_um"]) - float(row["truth_x_um"])
        )
        oracle_errors_um[target, realization, 1] = (
            float(row["estimate_y_um"]) - float(row["truth_y_um"])
        )
    oracle_difference_um = center_errors * 1e6 - oracle_errors_um
    oracle_difference_radial_um = np.linalg.norm(oracle_difference_um, axis=-1)
    oracle_difference_summary = {
        "rms_2d_um": float(np.sqrt(np.mean(oracle_difference_radial_um**2))),
        "median_2d_um": float(np.median(oracle_difference_radial_um)),
        "p90_2d_um": float(np.percentile(oracle_difference_radial_um, 90)),
        "max_2d_um": float(np.max(oracle_difference_radial_um)),
    }

    realization_rows: list[dict[str, object]] = []
    target_rows: list[dict[str, object]] = []
    for target, name in enumerate(target_names):
        for realization in range(nr):
            realization_rows.append(
                {
                    "target": name,
                    "target_index": target + 1,
                    "realization": realization + 1,
                    "relative_truth_x_um": relative_truth[target, realization, 0] * 1e6,
                    "relative_truth_y_um": relative_truth[target, realization, 1] * 1e6,
                    "estimate_x_um": estimates[target, realization, 0] * 1e6,
                    "estimate_y_um": estimates[target, realization, 1] * 1e6,
                    "error_x_um": center_errors[target, realization, 0] * 1e6,
                    "error_y_um": center_errors[target, realization, 1] * 1e6,
                    "error_2d_um": radial_center_um[target, realization],
                    "local_orbit_prediction_rmse_2d_um": per_case_local_rmse_um[
                        target, realization
                    ],
                    "oracle_error_x_um": oracle_errors_um[target, realization, 0],
                    "oracle_error_y_um": oracle_errors_um[target, realization, 1],
                    "error_vector_difference_from_oracle_um": oracle_difference_radial_um[
                        target, realization
                    ],
                }
            )
        target_rows.append(
            {
                "target": name,
                "target_index": target + 1,
                **summarize(center_errors[target]),
                "local_orbit_prediction_rmse_2d_um": float(
                    np.sqrt(np.mean(per_case_local_rmse_um[target] ** 2))
                ),
            }
        )
    write_rows(output / "per_realization_fits.csv", realization_rows)
    write_rows(output / "per_target_summary.csv", target_rows)
    write_rows(output / "summary.csv", [{"method": "two_sided_transport", **overall}])
    write_rows(
        output / "oracle_difference_summary.csv",
        [{"comparison": "two_sided_minus_oracle_error_vector", **oracle_difference_summary}],
    )

    target_rmse = np.asarray([float(row["rmse_2d_um"]) for row in target_rows])
    worst_indices = np.argsort(target_rmse)[::-1][:10]
    worst_lines = "\n".join(
        f"| {rank} | {target_names[index]} | {target_rmse[index]:.3f} | "
        f"{float(target_rows[index]['p90_2d_um']):.3f} | "
        f"{float(target_rows[index]['max_2d_um']):.3f} |"
        for rank, index in enumerate(worst_indices, start=1)
    )

    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.plot(np.arange(nt) + 1, target_rmse, color="#4472c4", marker="o", ms=3, lw=1)
    ax.axhline(
        np.median(target_rmse),
        color="#ed7d31",
        ls="--",
        label=f"target median {np.median(target_rmse):.3f} um",
    )
    ax.set_xlabel("Sextupole inventory index")
    ax.set_ylabel("Beam-relative center 2D RMSE [micrometers]")
    ax.set_title("Two-sided-BPM local orbit followed by K2-slope center inversion")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "per_target_center_rmse.png", dpi=180)
    plt.close(fig)

    report = f"""# Two-sided-BPM sextupole-center inversion

This end-to-end experiment replaces exact target-local coordinates with the
relative local orbit predicted from the nearest upstream/downstream BPM pair.
The maintained K2-slope center fit then consumes all 111 BPM x/y response
channels. Exact target-local orbit and target offset are evaluation-only.

- targets / latent realizations: {nt} / {nr} per target
- completed center fits: {nt * nr}
- bump protocol: five-point axial cross at `0.5 mm`
- K2 protocol: `(-0.02, 0, +0.02) m^-3`
- target offset: independent x/y uniform within `+/-350 micrometers`
- other 75 sextupoles: independent x/y `300 micrometer` RMS offsets
- all 113 quadrupoles: independent physical strength errors within `+/-1%`
- BPM noise/offset/gain errors and missing channels: none
- internal target orbit used by the two-sided inverse: **no**
- aggregate x/y/2D center RMSE: **{overall['x_rmse_um']:.3f} / {overall['y_rmse_um']:.3f} / {overall['rmse_2d_um']:.3f} micrometers**
- aggregate median / P90 / P99 / maximum: **{overall['median_2d_um']:.3f} / {overall['p90_2d_um']:.3f} / {overall['p99_2d_um']:.3f} / {overall['max_2d_um']:.3f} micrometers**
- per-target RMSE median / P90 / maximum: **{np.median(target_rmse):.3f} / {np.percentile(target_rmse, 90):.3f} / {np.max(target_rmse):.3f} micrometers**
- correlation between per-case local-orbit prediction RMSE and center error:
  **{correlation:.6f}**
- two-sided versus oracle center-error-vector RMS / median / P90 / maximum:
  **{oracle_difference_summary['rms_2d_um']:.3f} / {oracle_difference_summary['median_2d_um']:.3f} / {oracle_difference_summary['p90_2d_um']:.3f} / {oracle_difference_summary['max_2d_um']:.3f} micrometers**

For context, the same frozen tensor gave `13.913 micrometers` beam-relative
2D RMSE when commanded bump coordinates were used directly, while the frozen
oracle-local-orbit fit gave `5.870 micrometers`. The oracle value is the same
under an exact common coordinate translation because the source fit depends
only on local orbit minus center.

## Ten largest per-target center RMSE values

| rank | target | RMSE [micrometers] | P90 [micrometers] | max [micrometers] |
|---:|---|---:|---:|---:|
{worst_lines}

The result remains a noise-free SciBmad study. It tests propagation of the
two-sided local-orbit estimate through the physical center inverse, not real-
machine BPM or corrector calibration accuracy.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
