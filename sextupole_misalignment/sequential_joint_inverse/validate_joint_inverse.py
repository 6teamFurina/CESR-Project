#!/usr/bin/env python3
"""Validate raw machine-indexed scans and joint-inverse result artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WITHOUT = "without_quadrupole_misalignment"
WITH = "with_quadrupole_misalignment"
WITH_CORRECTED = "with_quadrupole_misalignment_corrected"
WITH_GTPSA_NOISY_CORRECTED = "with_quadrupole_misalignment_gtpsa_noisy_corrected"
CASES = (WITHOUT, WITH)
MODELS = (
    "physics_gls",
    "shared_target_local_ridge",
    "shared_joint_ridge",
    "shared_joint_random_feature",
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def load_metadata(path: Path) -> dict[str, object]:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def finite_file(path: Path, shape: tuple[int, ...] | None = None) -> np.ndarray:
    values = np.load(path, mmap_mode="r")
    if shape is not None and values.shape != shape:
        raise AssertionError(f"{path.name} shape {values.shape} != {shape}")
    if not np.all(np.isfinite(values)):
        raise AssertionError(f"Non-finite values in {path}")
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scan-root", type=Path, default=HERE / "results" / "exact_joint_machines"
    )
    parser.add_argument(
        "--analysis-dir", type=Path, default=HERE / "results" / "joint_inverse_analysis"
    )
    parser.add_argument(
        "--comparison-case",
        choices=(WITH, WITH_CORRECTED, WITH_GTPSA_NOISY_CORRECTED),
        default=WITH,
    )
    parser.add_argument("--expected-machines", type=int, default=16)
    parser.add_argument("--expected-targets", type=int, default=76)
    parser.add_argument("--expected-seed", type=int, default=20260823)
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()
    scan_root = args.scan_root.resolve()
    analysis = args.analysis_dir.resolve()
    comparison_case = args.comparison_case
    cases = (WITHOUT, comparison_case)
    latent_root = scan_root / "paired_latents"
    latent_metadata = load_metadata(latent_root / "latent_metadata.toml")
    if latent_metadata["format"] != "cesr-sequential-joint-paired-latents-v1":
        raise AssertionError("Wrong paired-latent metadata format")
    if int(latent_metadata["random_seed_base"]) != args.expected_seed:
        raise AssertionError("Unexpected paired-latent base seed")
    target_names = (latent_root / "target_names.txt").read_text(
        encoding="utf-8"
    ).splitlines()
    bpm_names = (latent_root / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    quadrupole_names = (latent_root / "quadrupole_names.txt").read_text(
        encoding="utf-8"
    ).splitlines()
    machine_count = args.expected_machines
    target_count = args.expected_targets
    if len(target_names) != 76 or target_count > len(target_names):
        raise AssertionError("Full paired latent target inventory is not 76")
    selected_names = target_names[:target_count]
    bpm_count = len(bpm_names)
    quadrupole_count = len(quadrupole_names)
    if bpm_count != 111:
        raise AssertionError(f"Expected 111 measurable BPMs, found {bpm_count}")

    sextupole_offsets = finite_file(
        latent_root / "sextupole_offsets.npy", (machine_count, 76, 2)
    )
    finite_file(latent_root / "corrector_gain_errors.npy")
    finite_file(latent_root / "k2_gain_errors.npy", (machine_count, 76))
    finite_file(
        latent_root / "quadrupole_relative_errors.npy",
        (machine_count, quadrupole_count),
    )
    finite_file(
        latent_root / "quadrupole_rolls.npy", (machine_count, quadrupole_count)
    )
    alignment_normals = finite_file(
        latent_root / "quadrupole_alignment_standard_normals.npy",
        (machine_count, quadrupole_count, 2),
    )
    bpm_gain_errors = finite_file(
        latent_root / "bpm_gain_errors.npy", (machine_count, bpm_count, 2)
    )
    finite_file(latent_root / "drift_directions.npy", (machine_count, 76, 2))

    case_metadata: dict[str, dict[str, object]] = {}
    references: dict[str, np.ndarray] = {}
    reference_targets: dict[str, np.ndarray] = {}
    corrected_validation: dict[str, object] = {}
    maximum_zero_bpm_error = 0.0
    maximum_zero_target_error = 0.0
    for case in cases:
        source = scan_root / case
        metadata = load_metadata(source / "scan_metadata.toml")
        case_metadata[case] = metadata
        expected_format = {
            WITHOUT: "cesr-sequential-joint-machine-scan-v1",
            WITH: "cesr-sequential-joint-machine-scan-v1",
            WITH_CORRECTED: "cesr-sequential-joint-machine-scan-v2",
            WITH_GTPSA_NOISY_CORRECTED: "cesr-sequential-joint-machine-scan-v3",
        }[case]
        if metadata["format"] != expected_format:
            raise AssertionError(f"Wrong scan format in {case}")
        if metadata["case"] != case:
            raise AssertionError(f"Wrong case label in {case}")
        if "SciBmad" not in str(metadata["engine"]):
            raise AssertionError("Primary engine is not SciBmad")
        if not str(metadata["lattice"]).endswith(
            "Latest_Lattice\\latest_cesr_scibmad_repaired.jl"
        ) and not str(metadata["lattice"]).endswith(
            "Latest_Lattice/latest_cesr_scibmad_repaired.jl"
        ):
            raise AssertionError("Scan does not use the repaired latest lattice")
        if int(metadata["machine_count"]) != machine_count:
            raise AssertionError(f"Machine count changed in {case}")
        if int(metadata["target_count"]) != target_count:
            raise AssertionError(f"Target count changed in {case}")
        expected_atomic_unit = (
            "all latent errors and baseline correction fixed across all target scans"
            if case in (WITH_CORRECTED, WITH_GTPSA_NOISY_CORRECTED)
            else "all latent errors fixed across all target scans"
        )
        if str(metadata["machine_atomic_unit"]) != expected_atomic_unit:
            raise AssertionError("Machine atomic-unit provenance changed")
        if (source / "target_names.txt").read_text(
            encoding="utf-8"
        ).splitlines() != selected_names:
            raise AssertionError(f"Target inventory mismatch in {case}")
        shape = (machine_count, target_count, 5, 3, bpm_count, 2)
        bpm = finite_file(source / "bpm_orbits.npy", shape)
        finite_file(source / "drift_bpm_orbits.npy", shape)
        target = finite_file(
            source / "target_orbits.npy", (machine_count, target_count, 5, 3, 2)
        )
        finite_file(
            source / "drift_target_orbits.npy",
            (machine_count, target_count, 5, 3, 2),
        )
        reference_bpm = finite_file(
            source / "reference_bpm_orbits.npy", (machine_count, bpm_count, 2)
        )
        reference_target = finite_file(
            source / "reference_target_orbits.npy", (machine_count, target_count, 2)
        )
        finite_file(source / "scan_seconds.npy", (machine_count, target_count))
        maximum_zero_bpm_error = max(
            maximum_zero_bpm_error,
            float(np.max(np.abs(bpm[:, :, 2, 1] - reference_bpm[:, None]))),
        )
        maximum_zero_target_error = max(
            maximum_zero_target_error,
            float(np.max(np.abs(target[:, :, 2, 1] - reference_target))),
        )
        references[case] = np.asarray(reference_bpm)
        reference_targets[case] = np.asarray(reference_target)
    if maximum_zero_bpm_error > 1.0e-11 or maximum_zero_target_error > 1.0e-11:
        raise AssertionError(
            "Zero-excitation scan states do not reproduce independent references"
        )
    if (
        float(case_metadata[comparison_case]["quadrupole_alignment_rms_m_per_plane"])
        != 5.0e-5
    ):
        raise AssertionError("Primary quadrupole alignment RMS is not 50 um/plane")
    if not np.any(references[comparison_case] != references[WITHOUT]):
        raise AssertionError("Paired quadrupole-alignment case has no physical effect")
    realized_alignment_rms = np.sqrt(np.mean((5.0e-5 * alignment_normals) ** 2, axis=(0, 1)))
    if np.any((realized_alignment_rms < 2.5e-5) | (realized_alignment_rms > 7.5e-5)):
        raise AssertionError("Finite quadrupole alignment draw is implausible")

    if comparison_case in (WITH_CORRECTED, WITH_GTPSA_NOISY_CORRECTED):
        corrected_source = scan_root / comparison_case
        corrected_metadata = case_metadata[comparison_case]
        if corrected_metadata.get("baseline_orbit_correction_applied") is not True:
            raise AssertionError("Corrected case does not declare baseline correction")
        baseline_count = int(corrected_metadata["baseline_corrector_count"])
        local_count = int(corrected_metadata["local_bump_corrector_count"])
        if baseline_count != 103 or local_count != 62:
            raise AssertionError("Unexpected baseline/local corrector registry size")
        zero_bpm = finite_file(
            corrected_source / "zero_offset_reference_bpm_orbits.npy",
            (machine_count, bpm_count, 2),
        )
        zero_target = finite_file(
            corrected_source / "zero_offset_reference_target_orbits.npy",
            (machine_count, target_count, 2),
        )
        uncorrected_bpm = finite_file(
            corrected_source / "uncorrected_bpm_orbits.npy",
            (machine_count, bpm_count, 2),
        )
        uncorrected_target = finite_file(
            corrected_source / "uncorrected_target_orbits.npy",
            (machine_count, target_count, 2),
        )
        baseline_commands = finite_file(
            corrected_source / "baseline_corrector_commands.npy",
            (machine_count, baseline_count),
        )
        finite_file(
            corrected_source / "baseline_local_corrector_fields.npy",
            (machine_count, local_count),
        )
        finite_file(
            corrected_source / "baseline_corrector_gain_errors.npy",
            (machine_count, baseline_count),
        )
        singular = finite_file(
            corrected_source / "baseline_response_singular_values.npy",
            (machine_count, baseline_count),
        )
        if np.any(np.diff(singular, axis=-1) > 0.0):
            raise AssertionError("Baseline response singular values are not descending")
        histories = np.load(corrected_source / "baseline_correction_history_m.npy")
        expected_history = (
            machine_count,
            int(corrected_metadata["baseline_maximum_iterations"]) + 1,
        )
        if histories.shape != expected_history or not np.all(
            np.isfinite(histories[:, 0])
        ):
            raise AssertionError("Invalid baseline-correction history")
        if comparison_case == WITH_GTPSA_NOISY_CORRECTED:
            if corrected_metadata.get("baseline_response_method") != "reference_gtpsa_orm":
                raise AssertionError("Noisy corrected case does not use a reference GTPSA ORM")
            if float(corrected_metadata["baseline_bpm_noise_rms_m_per_read"]) != 5.0e-6:
                raise AssertionError("Unexpected correction BPM per-read noise")
            if int(corrected_metadata["baseline_measurement_repeats"]) != 3072:
                raise AssertionError("Unexpected correction BPM repeat count")
            physical_histories = np.load(
                corrected_source / "baseline_physical_correction_history_m.npy"
            )
            if physical_histories.shape != expected_history or not np.all(
                np.isfinite(physical_histories[:, 0])
            ):
                raise AssertionError("Invalid physical correction history")
            reference_noise = finite_file(
                corrected_source / "baseline_reference_bpm_noise_m.npy",
                (machine_count, bpm_count, 2),
            )
            finite_file(
                corrected_source / "baseline_final_correction_bpm_noise_m.npy",
                (machine_count, bpm_count, 2),
            )
            finite_file(
                corrected_source / "baseline_validation_bpm_noise_m.npy",
                (machine_count, bpm_count, 2),
            )
            finite_file(
                corrected_source / "baseline_validation_residual_rms_m.npy",
                (machine_count,),
            )
            gtpsa_closure = finite_file(
                corrected_source / "baseline_response_closure_norm.npy",
                (machine_count,),
            )
            gtpsa_fd_relative = finite_file(
                corrected_source / "baseline_gtpsa_vs_fd_relative_l2.npy",
                (machine_count,),
            )
            gtpsa_fd_max_abs = finite_file(
                corrected_source / "baseline_gtpsa_vs_fd_max_abs.npy",
                (machine_count,),
            )
            expected_mean_noise = 5.0e-6 / np.sqrt(3072)
            realized_mean_noise = float(np.sqrt(np.mean(reference_noise**2)))
            if not 0.8 * expected_mean_noise < realized_mean_noise < 1.2 * expected_mean_noise:
                raise AssertionError("Reference BPM mean-noise draw is implausible")
            if np.max(gtpsa_closure) > 1.0e-12:
                raise AssertionError("GTPSA fixed-point response closure is too large")
            if np.max(gtpsa_fd_relative) > 1.0e-6:
                raise AssertionError("GTPSA and finite-difference ORM disagree")
            if np.max(gtpsa_fd_max_abs) > 1.0e-3:
                raise AssertionError("GTPSA and finite-difference ORM absolute error is too large")
        if not np.any(np.abs(baseline_commands) > 0.0):
            raise AssertionError("All saved baseline corrector commands are zero")

        corrected_bpm = references[comparison_case]
        corrected_target = finite_file(
            corrected_source / "reference_target_orbits.npy",
            (machine_count, target_count, 2),
        )
        gains = 1.0 + np.asarray(bpm_gain_errors)
        before_bpm = float(
            np.sqrt(np.mean(((uncorrected_bpm - zero_bpm) * gains) ** 2))
        )
        after_bpm = float(
            np.sqrt(np.mean(((corrected_bpm - zero_bpm) * gains) ** 2))
        )
        before_target = float(
            np.sqrt(
                np.mean(np.sum((uncorrected_target - zero_target) ** 2, axis=-1))
            )
        )
        after_target = float(
            np.sqrt(
                np.mean(np.sum((corrected_target - zero_target) ** 2, axis=-1))
            )
        )
        if not after_bpm < before_bpm:
            raise AssertionError("Baseline correction did not improve BPM restoration")
        if not after_target < before_target:
            raise AssertionError("Baseline correction did not improve target orbit")
        corrected_truth = np.asarray(sextupole_offsets)[:, :target_count] - corrected_target
        uncorrected_truth = (
            np.asarray(sextupole_offsets)[:, :target_count] - uncorrected_target
        )
        radius = float(corrected_metadata["bump_amplitude_m"])
        if np.mean(np.linalg.norm(corrected_truth, axis=-1) > radius) > np.mean(
            np.linalg.norm(uncorrected_truth, axis=-1) > radius
        ):
            raise AssertionError("Baseline correction worsened scan-radius coverage")
        paired_zero_bpm_error = float(
            np.max(np.abs(zero_bpm - references[WITHOUT]))
        )
        paired_zero_target_error = float(
            np.max(np.abs(zero_target - reference_targets[WITHOUT]))
        )
        if paired_zero_bpm_error > 1.0e-11 or paired_zero_target_error > 1.0e-11:
            raise AssertionError("Corrected-case zero-offset reference is not paired")

        correction_root = (
            HERE.parent
            / "quadrupole_orbit_correction"
            / "results"
            / "orbit_correction_50um"
        )
        command_artifact_error = None
        bpm_artifact_error = None
        target_artifact_error = None
        if (
            comparison_case == WITH_CORRECTED
            and machine_count == 16
            and target_count == 76
            and correction_root.is_dir()
        ):
            prior_commands = finite_file(
                correction_root / "corrector_commands.npy", (2, 16, 103)
            )[0]
            prior_bpm = finite_file(
                correction_root / "corrected_bpm_orbits.npy", (2, 16, 111, 2)
            )[0]
            prior_target = finite_file(
                correction_root / "corrected_target_orbits.npy", (2, 16, 76, 2)
            )[0]
            command_artifact_error = float(
                np.max(np.abs(baseline_commands - prior_commands))
            )
            bpm_artifact_error = float(np.max(np.abs(corrected_bpm - prior_bpm)))
            target_artifact_error = float(
                np.max(np.abs(corrected_target - prior_target))
            )
            if max(
                command_artifact_error, bpm_artifact_error, target_artifact_error
            ) > 1.0e-14:
                raise AssertionError(
                    "Corrected scans disagree with standalone correction artifact"
                )
        corrected_validation = {
            "baseline_bpm_rms_before_um": 1.0e6 * before_bpm,
            "baseline_bpm_rms_after_um": 1.0e6 * after_bpm,
            "baseline_target_2d_rms_before_um": 1.0e6 * before_target,
            "baseline_target_2d_rms_after_um": 1.0e6 * after_target,
            "corrected_truth_outside_bump_radius_fraction": float(
                np.mean(np.linalg.norm(corrected_truth, axis=-1) > radius)
            ),
            "paired_zero_bpm_max_error_m": paired_zero_bpm_error,
            "paired_zero_target_max_error_m": paired_zero_target_error,
            "standalone_command_max_error": command_artifact_error,
            "standalone_corrected_bpm_max_error_m": bpm_artifact_error,
            "standalone_corrected_target_max_error_m": target_artifact_error,
        }
        if comparison_case == WITH_GTPSA_NOISY_CORRECTED:
            corrected_validation.update(
                {
                    "baseline_response_method": corrected_metadata[
                        "baseline_response_method"
                    ],
                    "bpm_noise_rms_um_per_read": 1.0e6
                    * float(corrected_metadata["baseline_bpm_noise_rms_m_per_read"]),
                    "measurement_repeats": int(
                        corrected_metadata["baseline_measurement_repeats"]
                    ),
                    "bpm_mean_noise_rms_um": 1.0e6 * realized_mean_noise,
                    "gtpsa_vs_finite_difference_relative_l2_max": float(
                        np.max(gtpsa_fd_relative)
                    ),
                    "gtpsa_response_closure_norm_max": float(
                        np.max(gtpsa_closure)
                    ),
                }
            )

    analysis_metadata = json.loads(
        (analysis / "analysis_metadata.json").read_text(encoding="utf-8")
    )
    if analysis_metadata["format"] != "cesr-sequential-joint-inverse-analysis-v1":
        raise AssertionError("Wrong analysis format")
    if int(analysis_metadata["machine_count"]) != machine_count:
        raise AssertionError("Analysis machine count changed")
    if int(analysis_metadata["target_count"]) != target_count:
        raise AssertionError("Analysis target count changed")
    if analysis_metadata.get("cases", list(cases)) != list(cases):
        raise AssertionError("Analysis case provenance changed")
    train = set(map(int, analysis_metadata["train_indices_zero_based"]))
    validation = set(map(int, analysis_metadata["validation_indices_zero_based"]))
    test = set(map(int, analysis_metadata["test_indices_zero_based"]))
    if train & validation or train & test or validation & test:
        raise AssertionError("Machine split leakage detected")
    if train | validation | test != set(range(machine_count)):
        raise AssertionError("Machine splits do not cover the exact ensemble")
    if not train or not validation or not test:
        raise AssertionError("At least one machine split is empty")

    summary_rows = read_rows(analysis / "summary.csv")
    if len(summary_rows) != 14:
        raise AssertionError(f"Expected 14 evaluation rows, found {len(summary_rows)}")
    expected_keys = {
        ("fixed_nominal_physics", case, "physics_gls") for case in cases
    }
    expected_keys |= {
        (training, evaluation, model)
        for training in cases
        for evaluation in cases
        for model in MODELS[1:]
    }
    actual_keys = {
        (row["training_case"], row["evaluation_case"], row["model"])
        for row in summary_rows
    }
    if actual_keys != expected_keys:
        raise AssertionError("Evaluation matrix is incomplete or duplicated")
    numeric_fields = (
        "rmse_2d_um",
        "p90_2d_um",
        "p99_2d_um",
        "worst_target_rmse_2d_um",
        "fraction_below_50um",
    )
    for row in summary_rows:
        values = np.asarray([float(row[field]) for field in numeric_fields])
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise AssertionError("Invalid summary metric")
        if float(row["fraction_below_50um"]) > 1.0:
            raise AssertionError("Fraction metric exceeds one")
        for field in (
            "preferred_30um_aggregate_gate",
            "hard_50um_rmse_p99_all_target_gate",
        ):
            if row.get(field) not in {"True", "False"}:
                raise AssertionError(f"Missing Boolean acceptance field: {field}")

    target_rows = read_rows(analysis / "per_target_summary.csv")
    if len(target_rows) != len(summary_rows) * target_count:
        raise AssertionError("Per-target summary does not cover every evaluation row")
    if set(row["target"] for row in target_rows) != set(selected_names):
        raise AssertionError("Per-target analysis inventory changed")
    selection_rows = read_rows(analysis / "model_selection.csv")
    if not selection_rows:
        raise AssertionError("Model selection table is empty")
    if set(row["training_case"] for row in selection_rows) != set(cases):
        raise AssertionError("Model selection does not cover both training distributions")

    predictions = np.load(analysis / "held_out_predictions.npz")
    if not predictions.files:
        raise AssertionError("Prediction bundle is empty")
    for key in predictions.files:
        values = predictions[key]
        if not np.all(np.isfinite(values)):
            raise AssertionError(f"Non-finite held-out prediction: {key}")
    required_files = [
        analysis / "SUMMARY.md",
        analysis / "diagnostics.json",
        analysis / "held_out_model_comparison.png",
        analysis / "quadrupole_drift_domain_diagnostic.png",
    ]
    for training_case in cases:
        if training_case == WITHOUT:
            abbreviation = "no_quad_align"
        elif training_case == WITH_CORRECTED:
            abbreviation = "quad_align_50um_corrected"
        elif training_case == WITH_GTPSA_NOISY_CORRECTED:
            abbreviation = "quad_align_50um_gtpsa_noisy_corrected"
        else:
            abbreviation = "quad_align_50um"
        required_files.extend(
            [
                analysis / f"model_{abbreviation}_pipeline.npz",
                analysis / f"model_{abbreviation}_local_ridge.npz",
                analysis / f"model_{abbreviation}_joint_ridge.npz",
                analysis / f"model_{abbreviation}_joint_random_feature.npz",
            ]
        )
    for path in required_files:
        if not path.is_file() or path.stat().st_size == 0:
            raise AssertionError(f"Missing or empty artifact: {path}")

    validation_report = {
        "status": "PASS",
        "machine_count": machine_count,
        "target_count": target_count,
        "bpm_count": bpm_count,
        "quadrupole_count": quadrupole_count,
        "cases": list(cases),
        "maximum_zero_bpm_reference_error_m": maximum_zero_bpm_error,
        "maximum_zero_target_reference_error_m": maximum_zero_target_error,
        "realized_quadrupole_alignment_rms_um_per_plane": (
            realized_alignment_rms * 1.0e6
        ).tolist(),
        "evaluation_row_count": len(summary_rows),
        "per_target_row_count": len(target_rows),
        "machine_splits": {
            "train": sorted(train),
            "validation": sorted(validation),
            "test": sorted(test),
        },
        "corrected_protocol": corrected_validation or None,
    }
    if args.write_report:
        (analysis / "VALIDATION.json").write_text(
            json.dumps(validation_report, indent=2) + "\n", encoding="utf-8"
        )
    print(json.dumps(validation_report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
