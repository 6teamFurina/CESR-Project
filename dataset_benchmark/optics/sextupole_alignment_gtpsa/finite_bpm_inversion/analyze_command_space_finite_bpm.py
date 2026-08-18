#!/usr/bin/env python3
"""Finite-BPM, command-space beam-relative sextupole-center benchmark.

The inverse consumes only BPM closed-orbit tensors and commanded bump points.
Exact target-local orbits are loaded after all fits and are used solely to
define evaluation truth.  Consequently this script does not estimate the
absolute mechanical offset; it estimates center relative to the zero-bump
beam orbit in commanded-bump coordinates.
"""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = (
    HERE.parent
    / "direct_observable_nuisance_ablation"
    / "results"
    / "all_76_orbit_protocol"
)


def k2_slope(values: np.ndarray, delta: np.ndarray) -> np.ndarray:
    centered = delta - delta.mean()
    return np.einsum("bkc,k->bc", values, centered) / np.dot(centered, centered)


def source_matrix(xy: np.ndarray, center: np.ndarray) -> np.ndarray:
    x = xy[:, 0] - center[0]
    y = xy[:, 1] - center[1]
    return np.column_stack((0.5 * (x * x - y * y), x * y))


def fit_center(slopes: np.ndarray, xy_command: np.ndarray) -> np.ndarray:
    """Fit relative center while profiling out two response vectors."""
    from scipy.optimize import least_squares

    scale = np.sqrt(np.mean(slopes * slopes, axis=0))
    positive = scale[scale > 0]
    floor = np.median(positive) * 1e-8 if positive.size else 1.0
    normalized = slopes / np.maximum(scale, floor)

    def residual(center: np.ndarray) -> np.ndarray:
        source = source_matrix(xy_command, center)
        propagation = np.linalg.lstsq(source, normalized, rcond=1e-12)[0]
        return (source @ propagation - normalized).ravel()

    starts = [
        np.zeros(2),
        np.mean(xy_command, axis=0),
        np.array([xy_command[:, 0].min(), 0.0]),
        np.array([xy_command[:, 0].max(), 0.0]),
        np.array([0.0, xy_command[:, 1].min()]),
        np.array([0.0, xy_command[:, 1].max()]),
    ]
    solutions = [
        least_squares(
            residual,
            start,
            bounds=(-1.5e-3, 1.5e-3),
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
            max_nfev=500,
        )
        for start in starts
    ]
    return min(solutions, key=lambda result: np.dot(result.fun, result.fun)).x


def uniform_indices(total: int, requested: int) -> np.ndarray:
    if requested >= total:
        return np.arange(total, dtype=int)
    # Midpoints of equal ring sectors avoid privileging either end of the
    # inventory while remaining deterministic and independent of fit data.
    indices = np.floor((np.arange(requested) + 0.5) * total / requested).astype(int)
    if np.unique(indices).size != requested:
        raise RuntimeError("Uniform BPM selection produced duplicate indices")
    return indices


def load_bump_commands(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    rows.sort(key=lambda row: int(row["bump_index"]))
    return np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in rows
        ],
        dtype=float,
    )


def summarize(error_um: np.ndarray) -> dict[str, float]:
    return {
        "rmse_2d_um": float(np.sqrt(np.mean(error_um**2))),
        "median_2d_um": float(np.median(error_um)),
        "p90_2d_um": float(np.percentile(error_um, 90)),
        "p99_2d_um": float(np.percentile(error_um, 99)),
        "max_2d_um": float(np.max(error_um)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument(
        "--bpm-counts",
        default="1,2,4,8,16,32,64,111",
        help="Comma-separated deterministic uniform BPM subset sizes",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=HERE / "results" / "command_space_uniform_bpm",
    )
    args = parser.parse_args()

    source = args.input_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (source / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)

    # Fit inputs: measured BPM tensor and known commands only.
    bpm_orbits = np.load(source / "bpm_orbits.npy", mmap_mode="r")
    bump_commands = load_bump_commands(source / "bump_points.csv")
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta_k2 = levels * float(metadata["k2_step_m3"])
    nominal_k2 = int(np.flatnonzero(levels == 0)[0])
    zero_bump_candidates = np.flatnonzero(np.all(bump_commands == 0.0, axis=1))
    if zero_bump_candidates.size != 1:
        raise ValueError("Expected exactly one zero-command bump state")
    zero_bump = int(zero_bump_candidates[0])

    nt, nr, nb, nk, bpm_total, planes = bpm_orbits.shape
    if (nb, nk, bpm_total, planes) != (
        len(bump_commands), len(levels), len(bpm_names), 2
    ):
        raise ValueError("BPM tensor and metadata dimensions do not agree")

    requested_counts = [int(item) for item in args.bpm_counts.split(",")]
    counts = sorted(set(min(bpm_total, count) for count in requested_counts if count > 0))
    if not counts:
        raise ValueError("At least one positive BPM count is required")

    estimates_by_count: dict[int, np.ndarray] = {}
    selected_by_count: dict[int, np.ndarray] = {}
    for count in counts:
        indices = uniform_indices(bpm_total, count)
        selected_by_count[count] = indices
        estimates = np.zeros((nt, nr, 2), dtype=float)
        for target in range(nt):
            for realization in range(nr):
                selected = np.asarray(
                    bpm_orbits[target, realization, :, :, indices, :], dtype=float
                )
                # Advanced indexing puts the BPM axis first; restore
                # (bump, K2, BPM, plane), then flatten observable channels.
                if selected.shape[:2] != (nb, nk):
                    selected = np.moveaxis(selected, 0, 2)
                selected = selected.reshape(nb, nk, -1)
                slopes = k2_slope(selected, delta_k2)
                estimates[target, realization] = fit_center(slopes, bump_commands)
        estimates_by_count[count] = estimates
        print(f"completed command-space fits with {count} BPMs")

    # Evaluation-only truth is loaded after the fits.  It is never passed to
    # fit_center or used to choose a BPM subset.
    target_truth = np.load(source / "target_truth.npy")
    target_orbits = np.load(source / "target_orbits.npy", mmap_mode="r")
    zero_bump_orbit = np.asarray(
        target_orbits[:, :, zero_bump, nominal_k2, :], dtype=float
    )
    relative_truth = target_truth - zero_bump_orbit
    actual_bump_displacement = (
        np.asarray(target_orbits[:, :, :, nominal_k2, :], dtype=float)
        - zero_bump_orbit[:, :, None, :]
    )
    bump_mapping_error = actual_bump_displacement - bump_commands[None, None, :, :]
    nonzero_bumps = np.arange(nb) != zero_bump
    bump_mapping_radial_um = (
        np.linalg.norm(bump_mapping_error[:, :, nonzero_bumps, :], axis=-1) * 1e6
    )
    bump_mapping_rms_um = float(np.sqrt(np.mean(bump_mapping_radial_um**2)))
    bump_mapping_median_um = float(np.median(bump_mapping_radial_um))
    bump_mapping_p90_um = float(np.percentile(bump_mapping_radial_um, 90))

    summary_rows = []
    realization_rows = []
    for count in counts:
        estimates = estimates_by_count[count]
        vector_error = estimates - relative_truth
        radial_um = np.linalg.norm(vector_error, axis=2) * 1e6
        row = {
            "bpm_count": count,
            "channel_count": 2 * count,
            **summarize(radial_um),
        }
        summary_rows.append(row)
        for target in range(nt):
            for realization in range(nr):
                realization_rows.append(
                    {
                        "bpm_count": count,
                        "target": target_names[target],
                        "target_index": target + 1,
                        "realization": realization + 1,
                        "relative_truth_x_um": relative_truth[target, realization, 0] * 1e6,
                        "relative_truth_y_um": relative_truth[target, realization, 1] * 1e6,
                        "estimate_x_um": estimates[target, realization, 0] * 1e6,
                        "estimate_y_um": estimates[target, realization, 1] * 1e6,
                        "error_2d_um": radial_um[target, realization],
                    }
                )
        print(row)

    for filename, rows in (
        ("summary.csv", summary_rows),
        ("per_realization_fits.csv", realization_rows),
    ):
        with (output / filename).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    with (output / "selected_bpms.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["bpm_count", "bpm_index", "bpm"])
        writer.writeheader()
        for count in counts:
            for index in selected_by_count[count]:
                writer.writerow(
                    {"bpm_count": count, "bpm_index": index + 1, "bpm": bpm_names[index]}
                )

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    x = np.asarray([row["bpm_count"] for row in summary_rows])
    y = np.asarray([row["rmse_2d_um"] for row in summary_rows])
    ax.plot(x, y, color="#4472C4", marker="o", lw=1.8)
    ax.set_xscale("log", base=2)
    ax.set_xticks(x, [str(value) for value in x])
    ax.set_xlabel("Uniformly retained BPMs")
    ax.set_ylabel("Beam-relative center 2D RMSE [micrometers]")
    ax.set_title("Command-space sextupole-center inversion")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "rmse_vs_bpm_count.png", dpi=180)
    plt.close(fig)

    best = summary_rows[-1]
    table_lines = "\n".join(
        f"| {row['bpm_count']} | {row['rmse_2d_um']:.3f} | "
        f"{row['median_2d_um']:.3f} | {row['p90_2d_um']:.3f} |"
        for row in summary_rows
    )
    report = f"""# Command-space finite-BPM baseline

This experiment fits the beam-relative sextupole center using only selected
BPM closed-orbit channels and the five known bump commands. Exact internal
target orbit is excluded from every fit and from BPM selection. It is loaded
only afterward to evaluate the relative truth `c_s - z_s0`.

- targets / latent realizations: {nt} / {nr} per target
- total fits: {nt * nr * len(counts)}
- K2 protocol: three outer points
- BPM selection: deterministic ring-uniform subsets
- BPM noise/offset/gain errors: none
- internal target orbit used by inverse: **no**
- absolute mechanical offset estimated: **no**
- {best['bpm_count']}-BPM beam-relative 2D RMSE:
  **{best['rmse_2d_um']:.3f} micrometers**
- nominal-command versus actual local-bump displacement 2D
  RMS / median / P90: **{bump_mapping_rms_um:.3f} /
  {bump_mapping_median_um:.3f} / {bump_mapping_p90_um:.3f} micrometers**

| retained BPMs | RMSE [micrometers] | median [micrometers] | P90 [micrometers] |
|---:|---:|---:|---:|
{table_lines}

The result measures command-space recovery of the beam-centering displacement.
It does not yet include reconstruction of the nominal local orbit required to
recover an absolute mechanical offset. The nearly flat result from 8 through
111 BPMs, together with the local-bump mapping diagnostic above, makes bump
calibration/model mismatch the next variable to test before optimizing BPM
placement. See `summary.csv` for the complete BPM-count ablation.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
