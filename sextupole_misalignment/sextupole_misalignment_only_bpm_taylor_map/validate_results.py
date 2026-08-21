#!/usr/bin/env python3
"""Independently validate the saved scan and inverse benchmark artifacts."""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def vector_summary(errors_m: np.ndarray) -> dict[str, float]:
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


def assert_close(saved: dict[str, str], calculated: dict[str, float], prefix: str = "") -> None:
    for key, value in calculated.items():
        saved_key = f"{prefix}{key}"
        if not np.isclose(float(saved[saved_key]), value, rtol=2e-12, atol=2e-12):
            raise AssertionError(
                f"Summary mismatch for {saved_key}: saved={saved[saved_key]}, calculated={value}"
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-dir", type=Path, default=HERE / "results" / "exact_scans")
    parser.add_argument("--analysis-dir", type=Path, default=HERE / "results" / "analysis")
    args = parser.parse_args()
    scan = args.scan_dir.resolve()
    analysis = args.analysis_dir.resolve()

    with (scan / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    if metadata["only_machine_error"] != "fixed x/y misalignment of all 76 active normal sextupoles":
        raise AssertionError("The scan does not declare the required isolated error model")
    expected_omissions = {
        "BPM noise/offset/gain/roll/missing channels",
        "time drift",
        "corrector calibration error",
        "target-K2 calibration error",
        "quadrupole strength/roll/misalignment error",
        "additional RF or lattice errors",
    }
    if set(metadata["omitted_errors"]) != expected_omissions:
        raise AssertionError("The declared omitted-error inventory changed")

    target_names = (scan / "target_names.txt").read_text(encoding="utf-8").splitlines()
    sextupole_names = (scan / "sextupole_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_names = (scan / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    nt, ns, nd = len(target_names), len(sextupole_names), len(bpm_names)
    nr = int(metadata["realization_count_per_target"])
    nb = int(metadata["bump_count"])
    nk = int(metadata["k2_count"])
    if (nt, ns, nd, nr, nb, nk) != (76, 76, 111, 1, 25, 5):
        raise AssertionError(f"Unexpected benchmark dimensions: {(nt, ns, nd, nr, nb, nk)}")

    bpm = np.load(scan / "bpm_orbits.npy", mmap_mode="r")
    target_orbits = np.load(scan / "target_orbits.npy", mmap_mode="r")
    target_truth = np.load(scan / "target_truth.npy")
    latent = np.load(scan / "latent_sextupole_offsets.npy")
    nominal_centers = np.load(scan / "nominal_target_centers.npy")
    if bpm.shape != (nt, nr, nb, nk, nd, 2):
        raise AssertionError(f"Unexpected BPM tensor shape: {bpm.shape}")
    if target_orbits.shape != (nt, nr, nb, nk, 2):
        raise AssertionError(f"Unexpected target-orbit tensor shape: {target_orbits.shape}")
    if target_truth.shape != (nt, nr, 2) or latent.shape != (nt, nr, ns, 2):
        raise AssertionError("Unexpected truth/latent tensor shape")
    if not all(np.all(np.isfinite(array)) for array in (bpm, target_orbits, target_truth, latent)):
        raise AssertionError("The exact scan contains non-finite values")
    inventory = rows(scan / "target_inventory.csv")
    for target, row in enumerate(inventory):
        sextupole = int(row["sextupole_inventory_index"]) - 1
        if row["target"] != target_names[target]:
            raise AssertionError("Target inventory order changed")
        if not np.array_equal(target_truth[target], latent[target, :, sextupole, :]):
            raise AssertionError(f"Target truth does not match latent offsets for {row['target']}")

    bump_rows = rows(scan / "bump_points.csv")
    bump = np.asarray(
        [(float(row["bump_x_command_m"]), float(row["bump_y_command_m"])) for row in bump_rows]
    )
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    zero_bump = np.flatnonzero(np.all(bump == 0.0, axis=1))
    zero_k2 = np.flatnonzero(levels == 0.0)
    if zero_bump.size != 1 or zero_k2.size != 1:
        raise AssertionError("Zero-bump or zero-K2 state is not unique")
    ib, ik = int(zero_bump[0]), int(zero_k2[0])

    predicted = np.load(analysis / "predicted_relative_local_orbits.npy")
    predicted_reference = np.load(analysis / "predicted_reference_absolute_orbits.npy")
    exact_reference = np.asarray(target_orbits[:, :, ib, ik, :])
    exact_relative = np.asarray(target_orbits) - exact_reference[:, :, None, None, :]
    if predicted.shape != exact_relative.shape or predicted_reference.shape != exact_reference.shape:
        raise AssertionError("Predicted local-orbit output shape changed")

    nonzero_bumps = np.arange(nb) != ib
    calculated_local = {
        "relative_local_orbit_nominal_k2_nonzero_bumps": vector_summary(
            (predicted - exact_relative)[:, :, nonzero_bumps, ik, :]
        ),
        "relative_local_orbit_all_states": vector_summary(predicted - exact_relative),
        "absolute_zero_bump_reference_orbit": vector_summary(
            predicted_reference - exact_reference
        ),
    }
    saved_local = {row["quantity"]: row for row in rows(analysis / "local_orbit_summary.csv")}
    for name, summary in calculated_local.items():
        assert_close(saved_local[name], summary)

    true_total_center = nominal_centers[:, None, :] + target_truth
    relative_truth = true_total_center - exact_reference
    saved_summary = {row["method"]: row for row in rows(analysis / "summary.csv")}
    expected_methods = {
        "fd_linear_source_predicted",
        "fd_quartic_source_predicted",
        "quadratic_o_derivative_predicted",
        "chain_rule_o_derivative_predicted",
        "o_taylor_order3_nominal_local",
        "o_taylor_order4_nominal_local",
        "o_taylor_order5_nominal_local",
        "o_taylor_order4_all_state_local",
        "fd_quartic_source_oracle_local",
        "o_taylor_order4_oracle_all_state_local",
    }
    if set(saved_summary) != expected_methods:
        raise AssertionError("Inverse method inventory changed")
    for method, saved in saved_summary.items():
        estimate = np.load(analysis / f"{method}_relative_center_estimates.npy")
        if estimate.shape != relative_truth.shape or not np.all(np.isfinite(estimate)):
            raise AssertionError(f"Invalid estimate tensor for {method}")
        assert_close(saved, vector_summary(estimate - relative_truth), "relative_")
        absolute_increment = estimate + predicted_reference - nominal_centers[:, None, :]
        assert_close(
            saved,
            vector_summary(absolute_increment - target_truth),
            "absolute_increment_",
        )

    diagnostics = rows(analysis / "fit_diagnostics.csv")
    if len(diagnostics) != nt * nr * 6:
        raise AssertionError(f"Unexpected diagnostic row count: {len(diagnostics)}")
    expected_ranks = {
        "o_taylor_order3_nominal_local": 10,
        "o_taylor_order4_nominal_local": 20,
        "o_taylor_order5_nominal_local": 34,
        "o_taylor_order4_all_state_local": 20,
    }
    for method, expected_rank in expected_ranks.items():
        subset = [row for row in diagnostics if row["method"] == method]
        ranks = {int(row["taylor_design_rank"]) for row in subset}
        if len(subset) != nt * nr or ranks != {expected_rank}:
            raise AssertionError(f"Taylor design rank failure for {method}: {ranks}")
    derivative_rows = [
        row for row in diagnostics if row["method"].endswith("o_derivative_predicted")
    ]
    if any(
        not np.isfinite(float(row["inverse_condition"]))
        or int(row["retained_rows"]) < 2
        for row in derivative_rows
    ):
        raise AssertionError("Derivative inverse contains an invalid solve")

    print("PASS: isolated-error declaration and omitted-error inventory")
    print(f"PASS: exact SciBmad tensor dimensions ({nt} targets, {nt*nr*nb*nk} states)")
    print("PASS: target truth equals the corresponding latent sextupole offset")
    print("PASS: local-orbit and all inverse summaries independently reproduce the CSV values")
    print("PASS: all derivative solves are finite and every Taylor design has full rank")

    twiss_map_dir = HERE / "results" / "smoke_gtpsa_maps"
    fixed_map_dir = HERE / "results" / "smoke_gtpsa_fixed_point_maps"
    if twiss_map_dir.is_dir() and fixed_map_dir.is_dir():
        twiss_map = np.load(twiss_map_dir / "k2_offset_derivatives.npy")
        fixed_map = np.load(fixed_map_dir / "k2_offset_derivatives.npy")
        relative_difference = float(
            np.linalg.norm(fixed_map - twiss_map) / np.linalg.norm(twiss_map)
        )
        timing = rows(fixed_map_dir / "map_timings.csv")[0]
        final_residual = float(timing["final_fixed_point_residual"])
        if relative_difference > 1e-11 or final_residual > 1e-10:
            raise AssertionError(
                "The direct fixed-point GTPSA smoke map does not reproduce periodic Twiss"
            )
        print(
            "PASS: fixed-point and periodic-Twiss direct GTPSA maps agree "
            f"(relative L2={relative_difference:.3e}, residual={final_residual:.3e})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
