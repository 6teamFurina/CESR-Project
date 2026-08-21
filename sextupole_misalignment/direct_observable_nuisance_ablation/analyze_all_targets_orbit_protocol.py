#!/usr/bin/env python3
"""Invert and summarize the 76-target economical orbit protocol."""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_protocol_subsampling import fit_center, k2_slope


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir", type=Path,
        default=here / "results" / "all_76_orbit_protocol",
    )
    args = parser.parse_args()
    source = args.input_dir.resolve()
    with (source / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    orbit = np.load(source / "bpm_orbits.npy")
    target_orbit = np.load(source / "target_orbits.npy")
    truth = np.load(source / "target_truth.npy")
    latent_sext = np.load(source / "latent_sextupole_offsets.npy")
    latent_quad = np.load(source / "latent_quadrupole_relative_errors.npy")
    names = (source / "target_names.txt").read_text().splitlines()
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta = levels * float(metadata["k2_step_m3"])
    nominal = int(np.flatnonzero(levels == 0)[0])
    nt, nr = orbit.shape[:2]
    orbit = orbit.reshape(nt, nr, orbit.shape[2], orbit.shape[3], -1)
    estimates = np.zeros_like(truth)

    for target in range(nt):
        for realization in range(nr):
            slopes = k2_slope(orbit[target, realization], delta)
            xy = target_orbit[target, realization, :, nominal, :]
            estimates[target, realization] = fit_center(slopes, xy)
        print(f"inversion {target + 1}/{nt}: {names[target]}")

    vector_error = estimates - truth
    radial_um = np.linalg.norm(vector_error, axis=2) * 1e6
    per_realization_rows = []
    per_target_rows = []
    for target, name in enumerate(names):
        for realization in range(nr):
            per_realization_rows.append({
                "target": name,
                "target_index": target + 1,
                "realization": realization + 1,
                "truth_x_um": truth[target, realization, 0] * 1e6,
                "truth_y_um": truth[target, realization, 1] * 1e6,
                "estimate_x_um": estimates[target, realization, 0] * 1e6,
                "estimate_y_um": estimates[target, realization, 1] * 1e6,
                "error_2d_um": radial_um[target, realization],
            })
        errors = radial_um[target]
        per_target_rows.append({
            "target": name,
            "target_index": target + 1,
            "rmse_2d_um": float(np.sqrt(np.mean(errors**2))),
            "median_2d_um": float(np.median(errors)),
            "p90_2d_um": float(np.percentile(errors, 90)),
            "max_2d_um": float(np.max(errors)),
            "bias_x_um": float(np.mean(vector_error[target, :, 0]) * 1e6),
            "bias_y_um": float(np.mean(vector_error[target, :, 1]) * 1e6),
        })

    for filename, rows in (
        ("per_realization_fits.csv", per_realization_rows),
        ("per_target_summary.csv", per_target_rows),
    ):
        with (source / filename).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    target_rmse = np.asarray([row["rmse_2d_um"] for row in per_target_rows])
    aggregate_rmse = float(np.sqrt(np.mean(radial_um**2)))
    aggregate_median = float(np.median(radial_um))
    aggregate_p90 = float(np.percentile(radial_um, 90))
    aggregate_p99 = float(np.percentile(radial_um, 99))
    aggregate_max = float(np.max(radial_um))
    worst_indices = np.argsort(target_rmse)[::-1][:10]

    other_mask = np.ones((nt, nt), dtype=bool)
    np.fill_diagonal(other_mask, False)
    other_offsets = latent_sext[other_mask[:, None, :].repeat(nr, axis=1)].reshape(-1, 2)
    other_x_rms = float(np.sqrt(np.mean(other_offsets[:, 0] ** 2)) * 1e6)
    other_y_rms = float(np.sqrt(np.mean(other_offsets[:, 1] ** 2)) * 1e6)

    fig, ax = plt.subplots(figsize=(14, 5.5))
    x = np.arange(nt)
    ax.plot(x, target_rmse, color="#4472C4", marker="o", ms=3, lw=1)
    ax.axhline(np.median(target_rmse), color="#ED7D31", ls="--", lw=1.5, label="target median")
    ax.set_xticks(x[::2], [names[i] for i in x[::2]], rotation=90, fontsize=7)
    ax.set_ylabel("Per-target 2D RMSE [µm]")
    ax.set_title("All 76 sextupoles: 5 axial bumps × 3 K2 × nominal K1")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(source / "per_target_rmse.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.hist(radial_um.ravel(), bins=30, color="#4472C4", edgecolor="white")
    ax.axvline(aggregate_median, color="#ED7D31", ls="--", label=f"median {aggregate_median:.2f} µm")
    ax.axvline(aggregate_p90, color="#C00000", ls=":", label=f"P90 {aggregate_p90:.2f} µm")
    ax.set_xlabel("2D center error [µm]")
    ax.set_ylabel("Realizations")
    ax.set_title("Distribution over 76 targets × 8 nuisance realizations")
    ax.legend()
    fig.tight_layout()
    fig.savefig(source / "all_realization_error_histogram.png", dpi=180)
    plt.close(fig)

    worst_lines = "\n".join(
        f"| {rank} | {names[index]} | {target_rmse[index]:.3f} | "
        f"{per_target_rows[index]['p90_2d_um']:.3f} | {per_target_rows[index]['max_2d_um']:.3f} |"
        for rank, index in enumerate(worst_indices, start=1)
    )
    report = f"""# All-76 economical orbit protocol

The latest repaired SciBmad CESR lattice was evaluated for all {nt} active
normal sextupoles. Each target has {nr} independent latent realizations and
uses five axial-cross bumps, three K2 levels `(-2K,0,+2K)`, and nominal K1
commands. Every tensor contains unknown offsets on the other 75 sextupoles and
independent physical quadrupole strength errors bounded by ±1%. Target-local
orbit coordinates are exact and no BPM noise is added.

- exact SciBmad states: {metadata['total_state_count']}
- aggregate fits: {nt * nr}
- aggregate 2D RMSE: **{aggregate_rmse:.3f} µm**
- aggregate median / P90 / P99: **{aggregate_median:.3f} / {aggregate_p90:.3f} / {aggregate_p99:.3f} µm**
- maximum realization error: **{aggregate_max:.3f} µm**
- per-target RMSE median / P90 / maximum: **{np.median(target_rmse):.3f} / {np.percentile(target_rmse, 90):.3f} / {np.max(target_rmse):.3f} µm**
- realized other-sextupole x/y RMS: **{other_x_rms:.3f} / {other_y_rms:.3f} µm**
- realized quadrupole-error range: **{latent_quad.min()*100:.4f}% to {latent_quad.max()*100:.4f}%**
- SciBmad generation wall time: **{metadata['calculation_wall_seconds']:.1f} s**

## Ten largest per-target RMSE values

| rank | target | RMSE [µm] | P90 [µm] | max [µm] |
|---:|---|---:|---:|---:|
{worst_lines}

These values measure the present noise-free shared thin-sextupole source fit.
They do not include target-local-orbit uncertainty, BPM noise, missing BPMs, or
measured covariance, and therefore are not predicted machine accuracy.
"""
    (source / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
