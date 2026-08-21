#!/usr/bin/env python3
"""Measurement-level observable ablation for the paired SciBmad scan."""

from __future__ import annotations

import argparse
import csv
import math
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir", type=Path,
        default=here / "results" / "sex_09aw_paired_pilot",
    )
    parser.add_argument("--turns", type=int, default=256)
    parser.add_argument("--fft-size", type=int, default=8192)
    parser.add_argument("--energy-delta", type=float, default=1.0e-3)
    return parser.parse_args()


def estimate_frequency(signal: np.ndarray, fft_size: int) -> float:
    signal = np.asarray(signal, dtype=float)
    signal = signal - signal.mean()
    window = np.hanning(signal.size)
    spectrum = np.abs(np.fft.rfft(signal * window, n=fft_size)) ** 2
    frequency = np.fft.rfftfreq(fft_size)
    valid = (frequency >= 0.05) & (frequency <= 0.49)
    candidates = np.flatnonzero(valid)
    peak = candidates[np.argmax(spectrum[valid])]
    if 0 < peak < spectrum.size - 1:
        y0, y1, y2 = np.log(np.maximum(spectrum[peak - 1 : peak + 2], 1e-300))
        denominator = y0 - 2 * y1 + y2
        correction = 0.0 if abs(denominator) < 1e-20 else 0.5 * (y0 - y2) / denominator
    else:
        correction = 0.0
    return float((peak + correction) / fft_size)


def matrix_trajectory(matrix: np.ndarray, launch: np.ndarray, turns: int) -> np.ndarray:
    result = np.empty((turns, 4), dtype=float)
    state = launch.copy()
    for turn in range(turns):
        result[turn] = state
        state = matrix @ state
    return result


def harmonic_coefficients(signals: np.ndarray, tune: float) -> np.ndarray:
    turns = np.arange(signals.shape[0])
    kernel = np.exp(-2j * np.pi * tune * turns)
    return np.einsum("t,tbc->bc", kernel, signals) / signals.shape[0]


def direct_linear_readbacks(transport: np.ndarray, one_turn: np.ndarray, turns: int, fft_size: int):
    # Use the actual raw quantities available from a launch/angle scan.  These
    # avoid treating simulator Twiss beta or phase as if it were a readback.
    # Same-plane angle and position responses are the inputs to phase/beta
    # reconstruction; cross-plane responses are the coupling spectral lines.
    phase = np.concatenate((transport[:, 0, 1], transport[:, 1, 3]))
    beta = np.concatenate((transport[:, 0, 0], transport[:, 1, 2]))
    coupling = np.concatenate((
        transport[:, 0, 2], transport[:, 0, 3],
        transport[:, 1, 0], transport[:, 1, 1],
    ))

    # Tune remains a genuine TBT frequency fit.
    matrix = one_turn[:4, :4]
    tx = matrix_trajectory(matrix, np.array([1e-4, 0.0, 0.0, 0.0]), turns)
    ty = matrix_trajectory(matrix, np.array([0.0, 0.0, 1e-4, 0.0]), turns)
    qx = estimate_frequency(tx[:, 0], fft_size)
    qy = estimate_frequency(ty[:, 2], fft_size)
    return phase, beta, coupling, np.array([qx, qy])


def two_tunes_from_matrix(matrix: np.ndarray, turns: int, fft_size: int) -> np.ndarray:
    tx = matrix_trajectory(matrix, np.array([1e-4, 0.0, 0.0, 0.0]), turns)
    ty = matrix_trajectory(matrix, np.array([0.0, 0.0, 1e-4, 0.0]), turns)
    return np.array([
        estimate_frequency(tx[:, 0], fft_size),
        estimate_frequency(ty[:, 2], fft_size),
    ])


def fixed_energy_readbacks(
    transport: np.ndarray,
    one_turn: np.ndarray,
    hessian: np.ndarray,
    delta: float,
    turns: int,
    fft_size: int,
):
    a = one_turn[:4, :4]
    b = one_turn[:4, 5]
    dispersion_orbits = []
    tunes = []
    for sign in (-1.0, 1.0):
        d = sign * delta
        closed = np.linalg.lstsq(np.eye(4) - a, b * d, rcond=1e-12)[0]
        orbit = np.einsum("boc,c->bo", transport[:, :, :4], closed)
        orbit += transport[:, :, 5] * d
        dispersion_orbits.append(orbit)
        point = np.zeros(6)
        point[:4] = closed
        point[5] = d
        shifted = a + np.einsum("ijk,k->ij", hessian[:, :4, :], point)
        tunes.append(two_tunes_from_matrix(shifted, turns, fft_size))
    dispersion = (dispersion_orbits[1] - dispersion_orbits[0]) / 2.0
    chromaticity = (tunes[1] - tunes[0]) / (2.0 * delta)
    return dispersion.ravel(), chromaticity


def k2_slopes(values: np.ndarray, k2_delta: np.ndarray) -> np.ndarray:
    centered = k2_delta - k2_delta.mean()
    return np.einsum("rbkc,k->rbc", values, centered) / np.dot(centered, centered)


def retain_sensitive(a: np.ndarray, b: np.ndarray, relative_floor: float = 1e-8):
    sensitivity = np.linalg.norm(a, axis=1)
    if not np.any(np.isfinite(sensitivity)):
        return a[:0], b[:0]
    threshold = max(np.nanmax(sensitivity) * relative_floor, np.finfo(float).tiny)
    keep = np.isfinite(b) & np.all(np.isfinite(a), axis=1) & (sensitivity > threshold)
    return a[keep], b[keep]


def orbit_equations(slopes: np.ndarray, xy: np.ndarray):
    x, y = xy[:, 0], xy[:, 1]
    design = np.column_stack((np.ones_like(x), x, y, 0.5 * x * x, x * y, 0.5 * y * y))
    coefficient = np.linalg.lstsq(design, slopes, rcond=1e-12)[0]
    fitted = design @ coefficient
    scale = np.sqrt(np.mean(slopes * slopes, axis=0))
    scale = np.maximum(scale, np.nanmedian(scale[scale > 0]) * 1e-6)
    rows, rhs = [], []
    for channel in range(slopes.shape[1]):
        h = np.array([
            [coefficient[3, channel], coefficient[4, channel]],
            [coefficient[4, channel], coefficient[5, channel]],
        ]) / scale[channel]
        g = coefficient[1:3, channel] / scale[channel]
        rows.extend(h)
        rhs.extend(-g)
    a, b = retain_sensitive(np.asarray(rows), np.asarray(rhs))
    residual_fraction = float(np.linalg.norm(slopes - fitted) / max(np.linalg.norm(slopes), 1e-30))
    return a, b, residual_fraction


def zero_crossing_equations(slopes: np.ndarray, xy: np.ndarray):
    design = np.column_stack((np.ones(xy.shape[0]), xy))
    coefficient = np.linalg.lstsq(design, slopes, rcond=1e-12)[0]
    fitted = design @ coefficient
    scale = np.sqrt(np.mean(slopes * slopes, axis=0))
    positive = scale[scale > 0]
    floor = np.nanmedian(positive) * 1e-6 if positive.size else 1.0
    scale = np.maximum(scale, floor)
    a = (coefficient[1:3, :] / scale).T
    b = -coefficient[0, :] / scale
    a, b = retain_sensitive(a, b)
    residual_fraction = float(np.linalg.norm(slopes - fitted) / max(np.linalg.norm(slopes), 1e-30))
    return a, b, residual_fraction


def normalized_block(a: np.ndarray, b: np.ndarray):
    norm = np.linalg.norm(a)
    return (a / norm, b / norm) if norm > 0 else (a, b)


def solve_center(blocks):
    usable = [normalized_block(a, b) for a, b in blocks if a.shape[0] >= 2]
    if not usable:
        return np.full(2, np.nan), np.inf
    a = np.vstack([item[0] for item in usable])
    b = np.concatenate([item[1] for item in usable])
    center = np.linalg.lstsq(a, b, rcond=1e-12)[0]
    singular = np.linalg.svd(a, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else np.inf
    return center, condition


def summarize(name: str, estimates: np.ndarray, truth: np.ndarray):
    vector_error = estimates - truth
    radial_um = np.linalg.norm(vector_error, axis=1) * 1e6
    return {
        "model": name,
        "rmse_2d_um": float(np.sqrt(np.mean(radial_um**2))),
        "median_2d_um": float(np.median(radial_um)),
        "p90_2d_um": float(np.percentile(radial_um, 90)),
        "max_2d_um": float(np.max(radial_um)),
        "bias_x_um": float(np.mean(vector_error[:, 0]) * 1e6),
        "bias_y_um": float(np.mean(vector_error[:, 1]) * 1e6),
    }


def main() -> int:
    args = parse_args()
    source = args.input_dir.resolve()
    with (source / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    orbit = np.load(source / "bpm_orbits.npy")
    target_orbit = np.load(source / "target_orbits.npy")
    transport = np.load(source / "bpm_jacobians.npy")
    one_turn = np.load(source / "one_turn_jacobians.npy")
    hessian = np.load(source / "one_turn_hessians.npy")
    orm = np.load(source / "corrector_orm.npy")
    truth = np.load(source / "target_truth.npy")
    nr, nb, nk = orbit.shape[:3]
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    k2_delta = levels * float(metadata["k2_step_m3"])
    nominal_k2 = int(np.flatnonzero(levels == 0)[0])
    local_xy = target_orbit[:, :, nominal_k2, :]

    feature_lists = {key: [] for key in (
        "phase", "beta", "coupling", "tune", "dispersion", "chromaticity",
    )}
    for realization in range(nr):
        realization_features = {key: [] for key in feature_lists}
        for bump in range(nb):
            bump_features = {key: [] for key in feature_lists}
            for k2 in range(nk):
                phase, beta, coupling, tune = direct_linear_readbacks(
                    transport[realization, bump, k2], one_turn[realization, bump, k2],
                    args.turns, args.fft_size,
                )
                dispersion, chromaticity = fixed_energy_readbacks(
                    transport[realization, bump, k2], one_turn[realization, bump, k2],
                    hessian[realization, bump, k2], args.energy_delta,
                    args.turns, args.fft_size,
                )
                for key, value in (
                    ("phase", phase), ("beta", beta), ("coupling", coupling),
                    ("tune", tune), ("dispersion", dispersion),
                    ("chromaticity", chromaticity),
                ):
                    bump_features[key].append(value)
            for key in feature_lists:
                realization_features[key].append(np.asarray(bump_features[key]))
        for key in feature_lists:
            feature_lists[key].append(np.asarray(realization_features[key]))
        print(f"measurement extraction {realization + 1}/{nr}")
    features = {key: np.asarray(value) for key, value in feature_lists.items()}
    features["orbit_response"] = orm.reshape(nr, nb, nk, -1)
    features["orbit"] = orbit.reshape(nr, nb, nk, -1)

    slopes = {key: k2_slopes(value, k2_delta) for key, value in features.items()}
    direct_groups = [
        "phase", "beta", "coupling", "tune", "dispersion", "chromaticity",
        "orbit_response",
    ]
    model_names = ["orbit_only"] + [f"orbit_plus_{name}" for name in direct_groups]
    model_names += [
        "orbit_plus_feeddown_direct",
        "orbit_plus_direct_except_chromaticity",
        "orbit_plus_all_direct",
    ]
    estimates = {name: np.zeros((nr, 2)) for name in model_names}
    condition_numbers = {name: np.zeros(nr) for name in model_names}
    fit_rows = []

    for realization in range(nr):
        blocks = {}
        orbit_block = orbit_equations(slopes["orbit"][realization], local_xy[realization])
        blocks["orbit"] = orbit_block[:2]
        fit_rows.append((realization + 1, "orbit", orbit_block[2]))
        for group in direct_groups:
            result = zero_crossing_equations(slopes[group][realization], local_xy[realization])
            blocks[group] = result[:2]
            fit_rows.append((realization + 1, group, result[2]))

        configurations = {"orbit_only": [blocks["orbit"]]}
        configurations.update({
            f"orbit_plus_{group}": [blocks["orbit"], blocks[group]]
            for group in direct_groups
        })
        configurations["orbit_plus_direct_except_chromaticity"] = [
            blocks["orbit"],
            *[blocks[group] for group in direct_groups if group != "chromaticity"],
        ]
        configurations["orbit_plus_feeddown_direct"] = [
            blocks["orbit"],
            *[blocks[group] for group in (
                "phase", "beta", "coupling", "tune", "orbit_response",
            )],
        ]
        configurations["orbit_plus_all_direct"] = [
            blocks["orbit"], *[blocks[group] for group in direct_groups],
        ]
        for name, selected_blocks in configurations.items():
            estimates[name][realization], condition_numbers[name][realization] = solve_center(selected_blocks)

    summaries = [summarize(name, estimates[name], truth) for name in model_names]
    for row in summaries:
        row["median_condition_number"] = float(np.median(condition_numbers[row["model"]]))
        row["paired_better_than_orbit_count"] = int(np.sum(
            np.linalg.norm(estimates[row["model"]] - truth, axis=1)
            < np.linalg.norm(estimates["orbit_only"] - truth, axis=1)
        ))

    with (source / "fit_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    with (source / "surface_fit_residuals.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("realization", "observable_group", "relative_surface_residual"))
        writer.writerows(fit_rows)

    rows = []
    for realization in range(nr):
        row = {
            "realization": realization + 1,
            "truth_x_um": truth[realization, 0] * 1e6,
            "truth_y_um": truth[realization, 1] * 1e6,
        }
        for name in model_names:
            row[f"{name}_x_um"] = estimates[name][realization, 0] * 1e6
            row[f"{name}_y_um"] = estimates[name][realization, 1] * 1e6
            row[f"{name}_error_um"] = np.linalg.norm(estimates[name][realization] - truth[realization]) * 1e6
        rows.append(row)
    with (source / "per_realization_fits.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    plot_names = [
        "orbit_only", "orbit_plus_phase", "orbit_plus_beta", "orbit_plus_coupling",
        "orbit_plus_tune", "orbit_plus_dispersion", "orbit_plus_chromaticity",
        "orbit_plus_orbit_response", "orbit_plus_feeddown_direct",
        "orbit_plus_direct_except_chromaticity",
        "orbit_plus_all_direct",
    ]
    labels = [
        "Orbit only", "+ phase", "+ beta", "+ coupling", "+ tune", "+ dispersion",
        "+ chromaticity", "+ ORM", "+ feed-down set", "+ all except chrom.",
        "+ all direct",
    ]
    rmse = [next(row["rmse_2d_um"] for row in summaries if row["model"] == name) for name in plot_names]
    fig, ax = plt.subplots(figsize=(11, 5.6))
    colors = ["#4472C4"] + ["#70AD47"] * 7 + ["#5B9BD5", "#ED7D31", "#C00000"]
    ax.bar(np.arange(len(rmse)), rmse, color=colors)
    ax.set_xticks(np.arange(len(rmse)), labels, rotation=35, ha="right")
    ax.set_ylabel("2D center RMSE [µm]")
    ax.set_yscale("log")
    ax.set_title("Paired nuisance ablation: 75 sextupole offsets + quadrupole errors ≤1%")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(source / "observable_ablation_rmse.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharex=True, sharey=True)
    limits = np.array([-450.0, 450.0])
    for ax, name, title in zip(
        axes, ("orbit_only", "orbit_plus_all_direct"), ("Orbit only", "Orbit + all direct"),
    ):
        ax.plot(limits, limits, color="0.65", lw=1)
        ax.scatter(truth[:, 0] * 1e6, estimates[name][:, 0] * 1e6, label="x", marker="o")
        ax.scatter(truth[:, 1] * 1e6, estimates[name][:, 1] * 1e6, label="y", marker="s")
        ax.set_title(title)
        ax.set_xlabel("True offset [µm]")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("Fitted offset [µm]")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(source / "truth_vs_fit.png", dpi=180)
    plt.close(fig)

    orbit_summary = next(row for row in summaries if row["model"] == "orbit_only")
    all_summary = next(row for row in summaries if row["model"] == "orbit_plus_all_direct")
    without_chrom = next(row for row in summaries if row["model"] == "orbit_plus_direct_except_chromaticity")
    feeddown = next(row for row in summaries if row["model"] == "orbit_plus_feeddown_direct")
    report = f"""# Direct-observable nuisance ablation

This paired SciBmad pilot uses {nr} complete scan tensors for `{metadata['target_sextupole']}`.
Every tensor contains independent unknown offsets on the other 75 sextupoles
(Gaussian RMS {metadata['other_sextupole_offset_rms_m'] * 1e6:.0f} µm per plane) and
independent, fixed quadrupole strength errors uniformly bounded by
±{metadata['quadrupole_fraction'] * 100:.1f}%. The inverse never receives either
nuisance truth. Target local-orbit coordinates are treated as exact.

The direct readbacks are generated from measurement processes: fixed
launch/angle BPM trajectory differences (the raw inputs to phase, beta and
cross-plane coupling reconstruction), TBT tune spectra,
a fixed beam-energy ±delta probe (dispersion orbit difference and tune shift), and
finite differences of two actual correctors (ORM). The energy probe uses
`delta={args.energy_delta:g}` rather than retuning the harmon-master RF frequency.

| model | 2D RMSE [µm] | median [µm] | P90 [µm] | paired wins / {nr} |
|---|---:|---:|---:|---:|
| orbit only | {orbit_summary['rmse_2d_um']:.3f} | {orbit_summary['median_2d_um']:.3f} | {orbit_summary['p90_2d_um']:.3f} | — |
| orbit + feed-down direct readbacks | {feeddown['rmse_2d_um']:.3f} | {feeddown['median_2d_um']:.3f} | {feeddown['p90_2d_um']:.3f} | {feeddown['paired_better_than_orbit_count']} |
| orbit + all direct except chromaticity | {without_chrom['rmse_2d_um']:.3f} | {without_chrom['median_2d_um']:.3f} | {without_chrom['p90_2d_um']:.3f} | {without_chrom['paired_better_than_orbit_count']} |
| orbit + all direct | {all_summary['rmse_2d_um']:.3f} | {all_summary['median_2d_um']:.3f} | {all_summary['p90_2d_um']:.3f} | {all_summary['paired_better_than_orbit_count']} |

These are noise-free, structurally block-normalized physics fits, not predicted
machine precision. Each observable block receives equal total weight because a
measured covariance has not yet been supplied. Chromaticity is shown separately
because a centered sextupole has an intrinsic chromatic response, so its K2 slope
does not obey the same zero-at-center relation as ordinary linear feed-down.
"""
    (source / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
