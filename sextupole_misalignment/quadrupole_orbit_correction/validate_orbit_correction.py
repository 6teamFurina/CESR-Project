#!/usr/bin/env python3
"""Validate paired SciBmad quadrupole-offset orbit-correction artifacts."""

from __future__ import annotations

import argparse
import csv
import math
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_SOURCE = HERE / "results" / "orbit_correction_50um"


def finite(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    values = np.load(path)
    if values.shape != shape:
        raise AssertionError(f"{path.name}: expected {shape}, found {values.shape}")
    if not np.all(np.isfinite(values)):
        raise AssertionError(f"{path.name}: non-finite values")
    return np.asarray(values, dtype=float)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()
    source = args.source.resolve()

    with (source / "metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    methods = (source / "method_names.txt").read_text().splitlines()
    machine_count = int(metadata["machine_count"])
    method_count = len(methods)
    bpm_count = int(metadata["bpm_count"])
    target_count = int(metadata["target_count"])
    corrector_count = int(metadata["corrector_count"])
    quadrupole_count = int(metadata["quadrupole_count"])
    history_count = int(metadata["maximum_iterations"]) + 1
    scan_radius = float(metadata["scan_radius_m"])
    if methods != list(metadata["response_methods"]):
        raise AssertionError("Method names do not match metadata")

    reference = finite(source / "reference_bpm_readbacks.npy", (machine_count, bpm_count, 2))
    uncorrected = finite(source / "uncorrected_bpm_readbacks.npy", (machine_count, bpm_count, 2))
    corrected = finite(
        source / "corrected_bpm_readbacks.npy",
        (method_count, machine_count, bpm_count, 2),
    )
    reference_target = finite(
        source / "reference_target_orbits.npy", (machine_count, target_count, 2)
    )
    uncorrected_target = finite(
        source / "uncorrected_target_orbits.npy", (machine_count, target_count, 2)
    )
    corrected_target = finite(
        source / "corrected_target_orbits.npy",
        (method_count, machine_count, target_count, 2),
    )
    commands = finite(
        source / "corrector_commands.npy",
        (method_count, machine_count, corrector_count),
    )
    singular = finite(
        source / "response_singular_values.npy",
        (method_count, machine_count, corrector_count),
    )
    offsets = finite(
        source / "quadrupole_offsets_m.npy",
        (machine_count, quadrupole_count, 2),
    )
    sextupole_offsets = finite(
        source / "sextupole_offsets_m.npy",
        (machine_count, target_count, 2),
    )
    history = np.load(source / "bpm_rms_history_m.npy")
    if history.shape != (method_count, machine_count, history_count):
        raise AssertionError("BPM history shape mismatch")
    if not np.all(np.isfinite(history[..., 0])):
        raise AssertionError("Initial BPM history entries must be finite")
    if np.any(np.diff(singular, axis=-1) > 0):
        raise AssertionError("Response singular values are not descending")
    if not np.any(np.abs(commands) > 0):
        raise AssertionError("All corrector commands are zero")

    requested_rms = float(metadata["quadrupole_alignment_rms_m_per_plane"])
    realized = np.sqrt(np.mean(offsets**2, axis=(0, 1)))
    if np.any(np.abs(realized / requested_rms - 1.0) > 0.15):
        raise AssertionError(f"Unexpected realized quadrupole RMS: {realized}")

    before = float(np.sqrt(np.mean((uncorrected - reference) ** 2)))
    if before <= 0:
        raise AssertionError("Uncorrected paired BPM difference is zero")
    for index, method in enumerate(methods):
        after = float(np.sqrt(np.mean((corrected[index] - reference) ** 2)))
        target_before = float(
            np.sqrt(np.mean(np.sum((uncorrected_target - reference_target) ** 2, axis=-1)))
        )
        target_after = float(
            np.sqrt(
                np.mean(
                    np.sum((corrected_target[index] - reference_target) ** 2, axis=-1)
                )
            )
        )
        if not after < before:
            raise AssertionError(f"{method}: BPM correction did not improve")
        if not target_after < target_before:
            raise AssertionError(f"{method}: target trajectory did not improve")
        final_history = np.nanmin(history[index], axis=-1)
        saved_residual = np.sqrt(
            np.mean((corrected[index] - reference) ** 2, axis=(1, 2))
        )
        if not np.allclose(
            final_history, saved_residual, rtol=2e-8, atol=1e-14
        ):
            raise AssertionError(
                f"{method}: saved history does not match corrected BPM residual"
            )
        print(
            f"{method}: BPM {1e6*before:.6f} -> {1e6*after:.6f} um; "
            f"target 2D {1e6*target_before:.6f} -> {1e6*target_after:.6f} um"
        )

    with (source / "aggregate.csv").open(newline="") as stream:
        aggregate = list(csv.DictReader(stream))
    if [row["method"] for row in aggregate] != methods:
        raise AssertionError("Aggregate method order mismatch")
    reference_truth = sextupole_offsets - reference_target
    uncorrected_truth = sextupole_offsets - uncorrected_target
    reference_outside = float(
        np.mean(np.linalg.norm(reference_truth, axis=-1) > scan_radius)
    )
    uncorrected_outside = float(
        np.mean(np.linalg.norm(uncorrected_truth, axis=-1) > scan_radius)
    )
    for index, row in enumerate(aggregate):
        for key, value in row.items():
            if key not in {"method"} and not math.isfinite(float(value)):
                raise AssertionError(f"Non-finite aggregate field {key}")
        corrected_truth = sextupole_offsets - corrected_target[index]
        corrected_outside = float(
            np.mean(np.linalg.norm(corrected_truth, axis=-1) > scan_radius)
        )
        expected = {
            "before_bpm_rms_um": 1e6 * before,
            "after_bpm_rms_um": 1e6
            * float(np.sqrt(np.mean((corrected[index] - reference) ** 2))),
            "command_rms": float(np.sqrt(np.mean(commands[index] ** 2))),
            "max_abs_command": float(np.max(np.abs(commands[index]))),
            "reference_truth_outside_scan_radius_fraction": reference_outside,
            "uncorrected_truth_outside_scan_radius_fraction": uncorrected_outside,
            "corrected_truth_outside_scan_radius_fraction": corrected_outside,
        }
        for key, expected_value in expected.items():
            if not math.isclose(
                float(row[key]), expected_value, rel_tol=2e-10, abs_tol=1e-13
            ):
                raise AssertionError(f"{row['method']}: inconsistent aggregate {key}")
        if corrected_outside > uncorrected_outside:
            raise AssertionError(f"{row['method']}: scan-range coverage became worse")

    comparison_path = source / "response_comparison.csv"
    if {"reference_orm", "current_orm"}.issubset(methods):
        if not comparison_path.is_file():
            raise AssertionError("Missing response_comparison.csv")
        with comparison_path.open(newline="") as stream:
            comparison = list(csv.DictReader(stream))
        if len(comparison) != machine_count:
            raise AssertionError("Response-comparison machine count mismatch")
        reference_index = methods.index("reference_orm")
        current_index = methods.index("current_orm")
        for machine, row in enumerate(comparison):
            for key, value in row.items():
                if key != "machine" and not math.isfinite(float(value)):
                    raise AssertionError(f"Non-finite response-comparison field {key}")
            measured_delta_um = 1e6 * float(
                np.sqrt(
                    np.mean(
                        (
                            corrected[reference_index, machine]
                            - corrected[current_index, machine]
                        )
                        ** 2
                    )
                )
            )
            command_delta = float(
                np.sqrt(
                    np.mean(
                        (
                            commands[reference_index, machine]
                            - commands[current_index, machine]
                        )
                        ** 2
                    )
                )
            )
            if not math.isclose(
                float(row["corrected_bpm_method_difference_rms_um"]),
                measured_delta_um,
                rel_tol=2e-10,
                abs_tol=1e-12,
            ):
                raise AssertionError("Inconsistent corrected-BPM method difference")
            if not math.isclose(
                float(row["corrector_command_method_difference_rms"]),
                command_delta,
                rel_tol=2e-10,
                abs_tol=1e-14,
            ):
                raise AssertionError("Inconsistent command method difference")

    if not (source / "SUMMARY.md").is_file():
        raise AssertionError("Missing SUMMARY.md")
    print("Orbit-correction validation: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
