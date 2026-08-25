#!/usr/bin/env python3
"""Independently validate cross-response shapes and analytic composition."""

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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=HERE / "results" / "raw")
    parser.add_argument("--analysis-dir", type=Path, default=HERE / "results" / "analysis")
    parser.add_argument(
        "--exact-dir", type=Path, default=HERE / "results" / "exact_validation"
    )
    parser.add_argument(
        "--exact-analysis-dir",
        type=Path,
        default=HERE / "results" / "exact_validation_analysis",
    )
    args = parser.parse_args()
    raw = args.raw_dir.resolve()
    analysis = args.analysis_dir.resolve()
    exact = args.exact_dir.resolve()
    exact_analysis = args.exact_analysis_dir.resolve()

    with (raw / "response_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    if metadata["format"] != "cesr-sextupole-cross-response-v1":
        raise AssertionError("Unexpected raw-response format")
    if int(metadata["target_count"]) != 76 or int(metadata["observation_count"]) != 76:
        raise AssertionError("The response does not cover all 76 sextupoles")

    inventory = rows(raw / "target_inventory.csv")
    controls = rows(raw / "control_inventory.csv")
    if len(inventory) != 76 or len(controls) != int(metadata["corrector_count"]):
        raise AssertionError("Inventory count mismatch")

    kick = np.load(raw / "periodic_kick_response.npy")
    bump = np.load(raw / "bump_response.npy")
    source = np.load(raw / "sextupole_source_response.npy")
    local = np.load(raw / "local_bump_jacobian.npy")
    design = np.load(raw / "alignment_design.npy")
    response = np.load(raw / "target_control_response.npy")
    knobs = np.load(raw / "bump_knobs.npy")
    if kick.shape != (76, 76, 2, 2):
        raise AssertionError(f"Unexpected kick-response shape: {kick.shape}")
    if bump.shape != kick.shape or source.shape != kick.shape:
        raise AssertionError("Propagation matrices do not share one shape")
    if local.shape != (76, 2, 2) or design.shape != (76, 2, 76, 2, 2):
        raise AssertionError("Local-bump or alignment-design shape mismatch")
    if response.shape != (76, 2, len(controls)) or knobs.shape != (76, 2, len(controls)):
        raise AssertionError("Corrector response or bump-knob shape mismatch")
    if not all(np.all(np.isfinite(value)) for value in (kick, bump, source, local, design)):
        raise AssertionError("A raw response contains non-finite values")

    recomputed_bump = np.einsum("opc,tac->topa", response, knobs)
    np.testing.assert_allclose(recomputed_bump, bump, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(local, bump[np.arange(76), np.arange(76)], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        local,
        np.broadcast_to(np.eye(2), local.shape),
        rtol=0.0,
        atol=1e-10,
    )

    recomputed_design = np.zeros_like(design)
    for target in range(76):
        normal = source[target, :, :, 0]
        skew = source[target, :, :, 1]
        for bump_axis in range(2):
            dx, dy = local[target, :, bump_axis]
            recomputed_design[target, bump_axis, :, :, 0] = -dx * normal - dy * skew
            recomputed_design[target, bump_axis, :, :, 1] = +dy * normal - dx * skew
    np.testing.assert_allclose(recomputed_design, design, rtol=0.0, atol=0.0)

    if len(rows(raw / "cross_response_long.csv")) != 76 * 76 * 2:
        raise AssertionError("Long response table row count changed")
    expected_analysis_shapes = {
        "bump_matrix_152x152.npy": (152, 152),
        "periodic_kick_matrix_152x152.npy": (152, 152),
        "sextupole_source_matrix_152x152.npy": (152, 152),
        "shared_alignment_template_matrix_304x152.npy": (304, 152),
    }
    for filename, shape in expected_analysis_shapes.items():
        matrix = np.load(analysis / filename)
        if matrix.shape != shape or not np.all(np.isfinite(matrix)):
            raise AssertionError(f"Invalid analyzed matrix {filename}: {matrix.shape}")
    shared_template = np.load(
        analysis / "shared_alignment_template_matrix_304x152.npy"
    )
    expected_shared_template = np.transpose(design, (1, 2, 3, 0, 4)).reshape(
        304, 152
    )
    np.testing.assert_allclose(
        shared_template, expected_shared_template, rtol=0.0, atol=0.0
    )
    if len(rows(analysis / "per_target_locality.csv")) != 8 * 76:
        raise AssertionError("Per-target locality row count changed")
    if len(rows(analysis / "aggregate_locality.csv")) != 8:
        raise AssertionError("Aggregate locality row count changed")
    target_svd = rows(analysis / "per_target_design_svd.csv")
    if len(target_svd) != 76:
        raise AssertionError("Per-target design SVD row count changed")
    if {int(row["numerical_rank"]) for row in target_svd} != {2}:
        raise AssertionError("At least one per-target alignment design lost rank")
    svd_summary = rows(analysis / "svd_summary.csv")
    if len(svd_summary) != 5:
        raise AssertionError("SVD summary row count changed")
    block_summary = next(
        row
        for row in svd_summary
        if row["matrix"] == "separate_scan_block_design_23104x152"
    )
    block_singular = np.concatenate(
        [np.linalg.svd(design[target].reshape(304, 2), compute_uv=False) for target in range(76)]
    )
    if int(block_summary["numerical_rank"]) != 152:
        raise AssertionError("Separate-scan block design lost rank")
    np.testing.assert_allclose(
        float(block_summary["condition_number_retained"]),
        float(np.max(block_singular) / np.min(block_singular)),
        rtol=2e-15,
        atol=0.0,
    )

    with (exact / "scan_metadata.toml").open("rb") as stream:
        exact_metadata = tomllib.load(stream)
    if exact_metadata["format"] != "cesr-sextupole-cross-response-exact-validation-v1":
        raise AssertionError("Unexpected exact-validation format")
    if exact_metadata["lattice"] != metadata["lattice"]:
        raise AssertionError("Raw and exact responses used different lattices")
    exact_orbits = np.load(exact / "exact_sextupole_orbits.npy")
    centers = np.load(exact / "scenario_centers.npy")
    selected = rows(exact / "selected_targets.csv")
    states = rows(exact / "states.csv")
    selected_count = len(selected)
    if exact_orbits.shape != (selected_count, 2, 9, 76, 2):
        raise AssertionError(f"Unexpected exact-orbit shape: {exact_orbits.shape}")
    if centers.shape != (selected_count, 2, 2):
        raise AssertionError(f"Unexpected center shape: {centers.shape}")
    if not np.all(np.isfinite(exact_orbits)) or not np.all(np.isfinite(centers)):
        raise AssertionError("Exact validation contains non-finite values")
    if len(states) != 9 or len({row["state"] for row in states}) != 9:
        raise AssertionError("Exact validation does not contain nine unique scan states")

    state_index = {row["state"]: int(row["state_index"]) - 1 for row in states}
    bump_amplitude = float(exact_metadata["bump_amplitude_m"])
    k2_step = float(exact_metadata["k2_step_m3"])
    gradients = np.zeros((selected_count, 2, 2, 76, 2))
    for target in range(selected_count):
        for scenario in range(2):
            for axis, label in enumerate(("x", "y")):
                for bump_sign in (-1, 1):
                    positive_k2 = exact_orbits[
                        target, scenario, state_index[f"{label}_b{bump_sign}_k1"]
                    ]
                    negative_k2 = exact_orbits[
                        target, scenario, state_index[f"{label}_b{bump_sign}_k-1"]
                    ]
                    k2_slope = (positive_k2 - negative_k2) / (2.0 * k2_step)
                    gradients[target, scenario, axis] += (
                        bump_sign * k2_slope / (2.0 * bump_amplitude)
                    )
    exact_increment = gradients[:, 1] - gradients[:, 0]
    saved_exact = np.load(exact_analysis / "exact_alignment_gradient_increment.npy")
    np.testing.assert_allclose(saved_exact, exact_increment, rtol=2e-15, atol=1e-18)

    zero_index = state_index["zero"]
    predicted = np.zeros_like(exact_increment)
    for target, row in enumerate(selected):
        inventory_index = int(row["inventory_index"]) - 1
        local_orbit = exact_orbits[target, :, zero_index, inventory_index, :]
        center_increment = (centers[target, 1] - local_orbit[1]) - (
            centers[target, 0] - local_orbit[0]
        )
        predicted[target] = np.einsum(
            "aopc,c->aop", design[inventory_index], center_increment
        )
    saved_predicted = np.load(
        exact_analysis / "predicted_alignment_gradient_increment.npy"
    )
    np.testing.assert_allclose(saved_predicted, predicted, rtol=2e-15, atol=1e-18)
    residual_ratio = float(
        np.linalg.norm(predicted - exact_increment) / np.linalg.norm(exact_increment)
    )
    cosine = float(
        np.vdot(predicted, exact_increment)
        / (np.linalg.norm(predicted) * np.linalg.norm(exact_increment))
    )
    if residual_ratio >= 0.10 or cosine <= 0.99:
        raise AssertionError(
            f"Exact validation outside qualification: residual={residual_ratio}, cosine={cosine}"
        )
    print("Cross-response validation: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
