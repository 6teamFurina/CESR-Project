#!/usr/bin/env python3
"""Independently validate the sequential BPM/GTPSA local-orbit inverse."""

from __future__ import annotations

import csv
import inspect
import json
import sys
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
STUDY_ROOT = HERE.parent
RESULTS = HERE / "results" / "sequential_bpm_gtpsa_inverse"
SCAN_ROOT = STUDY_ROOT / "sequential_joint_inverse" / "results" / "exact_joint_machines"
CASE = "with_quadrupole_misalignment_gtpsa_noisy_corrected"
MODEL_DIR = HERE / "results" / "local_orbit_model"
KNOBS = (
    STUDY_ROOT
    / "quadrupole_affinity"
    / "exact_11_triplet_validation"
    / "results"
    / "bump_knobs"
    / "local_bump_knobs.csv"
)

sys.path.insert(0, str(HERE))
import analyze_sequential_bpm_gtpsa_inverse as analysis  # noqa: E402


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def assert_close(actual: float, expected: float, label: str) -> None:
    if not np.isclose(actual, expected, rtol=2.0e-10, atol=2.0e-10):
        raise AssertionError(f"{label}: {actual} != {expected}")


def verify_summary(
    saved: dict[str, str],
    prefix: str,
    errors_m: np.ndarray,
) -> None:
    calculated = analysis.summarize_vectors(errors_m)
    for key, value in calculated.items():
        assert_close(float(saved[f"{prefix}_{key}"]), value, f"{prefix}_{key}")


def main() -> int:
    metadata = json.loads((RESULTS / "analysis_metadata.json").read_text(encoding="utf-8"))
    source = SCAN_ROOT / CASE
    with (source / "scan_metadata.toml").open("rb") as stream:
        scan_metadata = tomllib.load(stream)
    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    bump_rows = rows(source / "bump_points.csv")
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    zero_bump = int(np.flatnonzero(np.all(bump_commands == 0.0, axis=1))[0])
    delta_k2 = np.asarray(scan_metadata["k2_delta_m3"], dtype=float)
    nominal_k2 = int(np.flatnonzero(delta_k2 == 0.0)[0])
    machine_count = int(scan_metadata["machine_count"])
    target_count = len(target_names)
    detector_count = len(bpm_names)
    bump_count = len(bump_commands)

    # Recompute both local-orbit products from only observable BPM arrays and
    # the GTPSA model. This validation path never needs target-local truth.
    model = analysis.load_orbit_model(
        MODEL_DIR,
        KNOBS,
        target_names,
        bpm_names,
        bump_commands,
    )
    physical_bpm = np.asarray(np.load(source / "bpm_orbits.npy", mmap_mode="r"))
    bpm_gains = np.asarray(np.load(SCAN_ROOT / "paired_latents" / "bpm_gain_errors.npy"))
    measured = physical_bpm * (1.0 + bpm_gains[:, None, None, None, :, :])
    reference = measured[:, :, zero_bump, nominal_k2]
    nominal = measured[:, :, :, nominal_k2]
    observed_relative = nominal - reference[:, :, None]
    residual = observed_relative.reshape(
        machine_count, target_count, bump_count, 2 * detector_count
    ) - model.model_bpm_bumps[None]
    expected_relative = np.broadcast_to(
        model.model_target_bumps[None],
        (machine_count, target_count, bump_count, 2),
    ).copy()
    expected_reference = np.broadcast_to(
        model.nominal_target_orbits[None],
        (machine_count, target_count, 2),
    ).copy()
    reference_residual = reference.reshape(
        machine_count, target_count, 2 * detector_count
    ) - model.nominal_bpm_orbits.reshape(-1)
    for target, row in enumerate(model.neighbor_rows):
        upstream = int(row["upstream_bpm_index"]) - 1
        downstream = int(row["downstream_bpm_index"]) - 1
        channels = np.asarray(
            (2 * upstream, 2 * upstream + 1, 2 * downstream, 2 * downstream + 1)
        )
        transport = model.two_sided_maps[target]
        expected_relative[:, target] += (
            np.take(residual[:, target], channels, axis=-1) @ transport.T
        )
        expected_reference[:, target] += (
            np.take(reference_residual[:, target], channels, axis=-1) @ transport.T
        )
    saved_relative_local = np.load(
        RESULTS / "deterministic_predicted_relative_local_orbits.npy"
    )
    saved_reference = np.load(
        RESULTS / "deterministic_predicted_reference_absolute_orbits.npy"
    )
    if np.max(np.abs(saved_relative_local - expected_relative)) > 2.0e-15:
        raise AssertionError("Saved relative local-orbit reconstruction changed")
    if np.max(np.abs(saved_reference - expected_reference)) > 2.0e-15:
        raise AssertionError("Saved absolute reference-orbit reconstruction changed")
    if np.max(np.abs(saved_relative_local[:, :, zero_bump])) > 2.0e-15:
        raise AssertionError("Zero-bump relative local orbit is not zero")

    # Evaluation-only truth is introduced only for metric reconstruction.
    exact_target = np.asarray(np.load(source / "target_orbits.npy", mmap_mode="r"))
    exact_reference = np.asarray(np.load(source / "reference_target_orbits.npy"))
    offsets = np.asarray(np.load(SCAN_ROOT / "paired_latents" / "sextupole_offsets.npy"))
    exact_relative = exact_target[:, :, :, nominal_k2] - exact_reference[:, :, None]
    relative_truth = offsets - exact_reference
    nonzero = np.arange(bump_count) != zero_bump

    local_rows = {
        (row["acquisition"], row["quantity"]): row
        for row in rows(RESULTS / "local_orbit_summary.csv")
    }
    calculated_local = analysis.summarize_vectors(
        saved_relative_local[:, :, nonzero] - exact_relative[:, :, nonzero]
    )
    calculated_reference = analysis.summarize_vectors(saved_reference - exact_reference)
    for key, value in calculated_local.items():
        assert_close(
            float(
                local_rows[
                    ("deterministic_static_readback", "relative_local_orbit_nonzero_bumps")
                ][key]
            ),
            value,
            f"deterministic local {key}",
        )
    for key, value in calculated_reference.items():
        assert_close(
            float(
                local_rows[
                    ("deterministic_static_readback", "absolute_reference_orbit")
                ][key]
            ),
            value,
            f"deterministic reference {key}",
        )

    center_rows = {
        (row["acquisition"], row["method"]): row
        for row in rows(RESULTS / "center_summary.csv")
    }
    deterministic_relative = np.load(
        RESULTS / "deterministic_relative_center_estimates.npy"
    )
    deterministic_absolute = np.load(
        RESULTS / "deterministic_absolute_offset_estimates.npy"
    )
    deterministic_row = center_rows[
        ("deterministic_static_readback", "bpm_gtpsa_two_sided")
    ]
    verify_summary(
        deterministic_row,
        "relative",
        deterministic_relative - relative_truth,
    )
    verify_summary(
        deterministic_row,
        "absolute",
        deterministic_absolute - offsets,
    )
    oracle_relative = np.load(RESULTS / "oracle_relative_center_estimates.npy")
    oracle_absolute = np.load(RESULTS / "oracle_absolute_offset_estimates.npy")
    oracle_row = center_rows[
        ("deterministic_static_readback", "exact_local_orbit_oracle")
    ]
    verify_summary(oracle_row, "relative", oracle_relative - relative_truth)
    verify_summary(oracle_row, "absolute", oracle_absolute - offsets)

    stochastic_relative_path = RESULTS / "stochastic_relative_center_estimates.npy"
    if stochastic_relative_path.exists():
        indices = np.load(RESULTS / "stochastic_machine_indices_zero_based.npy")
        stochastic_relative = np.load(stochastic_relative_path)
        stochastic_absolute = np.load(RESULTS / "stochastic_absolute_offset_estimates.npy")
        stochastic_row = center_rows[("stochastic_15_state_means", "bpm_gtpsa_two_sided")]
        verify_summary(
            stochastic_row,
            "relative",
            stochastic_relative - relative_truth[indices][None],
        )
        verify_summary(
            stochastic_row,
            "absolute",
            stochastic_absolute - offsets[indices][None],
        )

    # Structural leakage guard: the callable machine-facing inverse must not
    # mention any evaluation artifact, and the main program must persist its
    # output before the first evaluation-only load.
    machine_source = inspect.getsource(analysis.machine_facing_inverse)
    forbidden = ("target_orbits.npy", "reference_target_orbits.npy", "sextupole_offsets.npy")
    if any(name in machine_source for name in forbidden):
        raise AssertionError("Machine-facing inverse references evaluation truth")
    complete_source = Path(analysis.__file__).read_text(encoding="utf-8")
    save_position = complete_source.index(
        'save_machine_facing_result(output, "deterministic", deterministic)'
    )
    truth_position = complete_source.index(
        'np.load(source / "target_orbits.npy", mmap_mode="r")'
    )
    if not save_position < truth_position:
        raise AssertionError("Evaluation truth is loaded before machine-facing persistence")

    validation = {
        "status": "PASS",
        "format": metadata["format"],
        "case": metadata["case"],
        "machine_count": machine_count,
        "target_count": target_count,
        "bpm_count": detector_count,
        "deterministic_fit_count": machine_count * target_count,
        "stochastic_fit_count": int(
            np.prod(np.load(stochastic_relative_path).shape[:-1])
            if stochastic_relative_path.exists()
            else 0
        ),
        "maximum_recomputed_relative_local_difference_m": float(
            np.max(np.abs(saved_relative_local - expected_relative))
        ),
        "maximum_recomputed_reference_difference_m": float(
            np.max(np.abs(saved_reference - expected_reference))
        ),
        "truth_leakage_check": "PASS",
        "deterministic_relative_center_rmse_um": float(
            deterministic_row["relative_rmse_2d_um"]
        ),
        "deterministic_absolute_offset_rmse_um": float(
            deterministic_row["absolute_rmse_2d_um"]
        ),
    }
    (RESULTS / "VALIDATION.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
