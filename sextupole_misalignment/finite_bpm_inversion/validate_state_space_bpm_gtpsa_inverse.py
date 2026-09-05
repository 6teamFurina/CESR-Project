#!/usr/bin/env python3
"""Independent checks for the full-error state-space BPM/GTPSA inverse."""

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
SCAN_ROOT = STUDY_ROOT / "sequential_joint_inverse" / "results" / "exact_joint_machines"
CASE = "with_all_errors_gtpsa_nominal_corrected"
SOURCE = SCAN_ROOT / CASE
LATENT = SCAN_ROOT / "paired_latents"
MODEL = HERE / "results" / "local_orbit_model"
RESULTS = HERE / "results" / "state_space_sequential_bpm_gtpsa_inverse"

sys.path.insert(0, str(HERE))
import analyze_state_space_bpm_gtpsa_inverse as analysis  # noqa: E402


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def assert_summary(
    row: dict[str, str],
    relative_errors: np.ndarray,
    absolute_errors: np.ndarray,
) -> None:
    relative = analysis.base.summarize_vectors(relative_errors)
    absolute = analysis.base.summarize_vectors(absolute_errors)
    for prefix, calculated in (("relative", relative), ("absolute", absolute)):
        for key, value in calculated.items():
            saved = float(row[f"{prefix}_{key}"])
            if not np.isclose(saved, value, rtol=2.0e-11, atol=2.0e-11):
                raise AssertionError(
                    f"Saved {prefix}_{key} {saved} != recomputed {value}"
                )


def dense_state_space_check() -> float:
    """Compare the eigenspace implementation with direct dense conditioning."""
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in rows(SOURCE / "bump_points.csv")
        ]
    )
    with (SOURCE / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    delta_k2 = np.asarray(metadata["k2_delta_m3"], dtype=float)
    protocol = analysis.build_protocol(bump_commands, delta_k2, 16, 4)
    response = np.asarray(
        [
            [1.0, 0.3],
            [0.4, -0.8],
            [-0.2, 1.1],
            [0.7, 0.5],
            [-1.0, 0.2],
            [0.1, -0.6],
        ]
    )
    scalar_step = (10.0e-6) ** 2 / (16 * 8 - 1)
    operator = analysis.build_state_space_operator(
        response,
        protocol,
        5.0e-6,
        32,
        scalar_step,
    )
    rng = np.random.default_rng(71023)
    observations = rng.normal(size=(3, len(protocol.reference_times), 2)) * 2.0e-6
    compact = analysis.hidden_state_filtered_averages(observations, operator)

    reference = protocol.reference_times
    process_step = 0.5 * scalar_step
    process = process_step * np.minimum.outer(reference, reference)
    average_cross = np.asarray(
        [
            process_step
            * np.mean(np.minimum(times[:, None], reference[None, :]), axis=0)
            for times in protocol.signal_times
        ]
    )
    same = (
        protocol.reference_types[:, None] == protocol.reference_types[None, :]
    ).astype(float)
    state_identity = np.eye(2)
    observation_covariance = (
        np.kron(process, state_identity)
        + np.kron(
            np.eye(len(reference)), operator.projected_noise_covariance
        )
        + np.kron(
            same, operator.projected_noise_covariance / 32
        )
    )
    cross = np.kron(average_cross, state_identity)
    direct = np.asarray(
        [
            (
                cross
                @ np.linalg.solve(
                    observation_covariance,
                    observation.reshape(-1),
                )
            ).reshape(8, 2)
            for observation in observations
        ]
    )
    difference = float(np.max(np.abs(compact - direct)))
    if difference > 2.0e-15:
        raise AssertionError(f"Dense state-space difference is {difference} m")
    return difference


def analytic_profiled_jacobian_check() -> float:
    """Check the exact optimizer Jacobian against an independent FD audit."""
    local_orbits = np.asarray(
        [
            (-1.2e-3, 0.1e-3),
            (0.2e-3, -1.1e-3),
            (0.0, 0.0),
            (0.1e-3, 1.3e-3),
            (1.4e-3, -0.2e-3),
        ]
    )
    center = np.asarray((0.23e-3, -0.17e-3))
    normalized = np.random.default_rng(71024).normal(size=(5, 12))
    _, jacobian = analysis.base.profiled_residual_and_jacobian(
        normalized, local_orbits, center
    )
    step = 1.0e-9
    finite_difference = np.column_stack(
        [
            (
                analysis.base.profiled_residual_and_jacobian(
                    normalized,
                    local_orbits,
                    center + np.eye(2)[parameter] * step,
                )[0]
                - analysis.base.profiled_residual_and_jacobian(
                    normalized,
                    local_orbits,
                    center - np.eye(2)[parameter] * step,
                )[0]
            )
            / (2.0 * step)
            for parameter in range(2)
        ]
    )
    relative_difference = float(
        np.linalg.norm(jacobian - finite_difference)
        / max(np.linalg.norm(finite_difference), 1.0e-30)
    )
    require(
        relative_difference <= 1.0e-8,
        f"Profiled analytic Jacobian difference is {relative_difference}",
    )
    return relative_difference


def main() -> int:
    with (SOURCE / "scan_metadata.toml").open("rb") as stream:
        scan = tomllib.load(stream)
    result_metadata = json.loads(
        (RESULTS / "analysis_metadata.json").read_text(encoding="utf-8")
    )
    machine_count = int(scan["machine_count"])
    target_count = int(scan["target_count"])
    bpm_count = int(scan["bpm_count"])
    augmentations = int(result_metadata["stochastic_augmentations"])
    require((machine_count, target_count, bpm_count) == (16, 76, 111), "Scope changed")
    require(scan["baseline_orbit_correction_applied"] is True, "No orbit correction")
    require(scan["baseline_response_method"] == "reference_gtpsa_orm", "ORM is not GTPSA")
    require(scan["baseline_gtpsa_response_model"] == "nominal", "ORM knows latent errors")
    require(scan["baseline_gtpsa_validation_enabled"] is False, "FD ORM was reapplied")
    require(
        "inverse consumes these saved readbacks and receives no gain realization"
        in scan["observable_bpm_readback_semantics"],
        "Observable BPM readbacks were not materialized by the forward generator",
    )
    require(scan["target_scan_parallelism"] == "threads", "Target scan was not threaded")
    require(int(scan["target_scan_worker_count"]) >= 2, "Only one scan worker was used")
    require(
        int(scan["target_scan_worker_count"]) == int(scan["julia_thread_count"]),
        "The target scan did not use every available Julia CPU thread",
    )
    require(
        int(scan["blas_thread_count_during_scan"]) == 1,
        "Nested BLAS threading was not disabled during target scans",
    )
    require(scan["thread_equivalence_checked"] is True, "No serial/thread equivalence check")
    thread_fields = (
        "thread_equivalence_bpm_max_abs_m",
        "thread_equivalence_drift_bpm_max_abs_m",
        "thread_equivalence_target_max_abs_m",
        "thread_equivalence_drift_target_max_abs_m",
    )
    for field in thread_fields:
        require(float(scan[field]) <= 1.0e-13, f"{field} is too large")
    orm = np.load(SOURCE / "baseline_nominal_gtpsa_orm.npy")
    require(orm.shape == (2 * bpm_count, 103), "Saved nominal ORM shape changed")
    require(np.all(np.isfinite(orm)), "Saved nominal ORM is non-finite")

    observable = np.load(SOURCE / "observable_bpm_readbacks.npy", mmap_mode="r")
    drift_observable = np.load(
        SOURCE / "observable_drift_bpm_readbacks.npy", mmap_mode="r"
    )
    require(
        observable.shape == (16, 76, 5, 3, 111, 2),
        "Observable BPM tensor shape changed",
    )
    require(drift_observable.shape == observable.shape, "Drift BPM shape changed")
    physical_bpm = np.load(SOURCE / "bpm_orbits.npy", mmap_mode="r")
    drift_physical_bpm = np.load(SOURCE / "drift_bpm_orbits.npy", mmap_mode="r")
    bpm_gain_errors = np.load(LATENT / "bpm_gain_errors.npy", mmap_mode="r")
    gain_factors = 1.0 + bpm_gain_errors[:, None, None, None, :, :]
    observable_difference = float(
        np.max(np.abs(observable - physical_bpm * gain_factors))
    )
    drift_observable_difference = float(
        np.max(np.abs(drift_observable - drift_physical_bpm * gain_factors))
    )
    require(observable_difference == 0.0, "Observable BPM materialization changed")
    require(
        drift_observable_difference == 0.0,
        "Drift BPM materialization changed",
    )

    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in rows(SOURCE / "bump_points.csv")
        ]
    )
    delta_k2 = np.asarray(scan["k2_delta_m3"], dtype=float)
    zero_bump = int(np.flatnonzero(np.all(bump_commands == 0.0, axis=1))[0])
    zero_k2 = int(np.flatnonzero(delta_k2 == 0.0)[0])
    reference_bpm = np.load(SOURCE / "reference_bpm_orbits.npy", mmap_mode="r")
    reference_target = np.load(
        SOURCE / "reference_target_orbits.npy", mmap_mode="r"
    )
    physical_target = np.load(SOURCE / "target_orbits.npy", mmap_mode="r")
    zero_bpm_closure = float(
        np.max(
            np.abs(
                physical_bpm[:, :, zero_bump, zero_k2]
                - reference_bpm[:, None]
            )
        )
    )
    zero_target_closure = float(
        np.max(
            np.abs(
                physical_target[:, :, zero_bump, zero_k2] - reference_target
            )
        )
    )
    require(zero_bpm_closure <= 1.0e-13, "Corrected BPM scan baseline changed")
    require(
        zero_target_closure <= 1.0e-13,
        "Corrected target scan baseline changed",
    )
    require(float(scan["baseline_bpm_rms_after_um"]) < float(scan["baseline_bpm_rms_before_um"]), "BPM correction failed")
    require(float(scan["baseline_target_2d_rms_after_um"]) < float(scan["baseline_target_2d_rms_before_um"]), "Target orbit correction failed")

    relative = np.load(RESULTS / "filtered_fixed_template_relative_center_estimates.npy")
    absolute = np.load(RESULTS / "filtered_fixed_template_absolute_offset_estimates.npy")
    unfiltered_relative = np.load(
        RESULTS / "unfiltered_fixed_template_relative_center_estimates.npy"
    )
    require(relative.shape == (augmentations, 16, 76, 2), "Filtered estimate shape changed")
    require(absolute.shape == relative.shape, "Absolute estimate shape changed")
    require(
        np.max(np.abs(relative - unfiltered_relative)) > 0.0,
        "Hidden-state correction was not applied to center inputs",
    )

    exact_reference = np.load(SOURCE / "reference_target_orbits.npy")
    offsets = np.load(LATENT / "sextupole_offsets.npy")
    relative_truth = offsets - exact_reference
    center_rows = {
        (row["acquisition"], row["method"]): row
        for row in rows(RESULTS / "center_summary.csv")
    }
    filtered_row = center_rows[
        (
            "periodic_reference_time_series",
            "state_space_filtered_fixed_gtpsa_template",
        )
    ]
    assert_summary(
        filtered_row,
        relative - relative_truth[None],
        absolute - offsets[None],
    )
    require(float(filtered_row["relative_rmse_2d_um"]) < 30.0, "Relative RMSE gate failed")
    require(float(filtered_row["absolute_rmse_2d_um"]) < 30.0, "Absolute RMSE gate failed")

    unfiltered_bpm = np.load(RESULTS / "unfiltered_bpm_state_error_rms_m.npy")
    filtered_bpm = np.load(RESULTS / "filtered_bpm_state_error_rms_m.npy")
    unfiltered_bpm_rms = float(np.sqrt(np.mean(unfiltered_bpm**2)))
    filtered_bpm_rms = float(np.sqrt(np.mean(filtered_bpm**2)))
    require(filtered_bpm_rms < 0.2 * unfiltered_bpm_rms, "State filter did not suppress drift")

    # Structural estimator boundary.
    machine_functions = (
        analysis.hidden_state_filtered_averages,
        analysis.reconstruct_target_local_orbits,
        analysis.fixed_template_centers,
        analysis.fit_noise_aware_profiled_center,
        analysis.run_machine_facing_target,
    )
    forbidden = (
        "target_orbits.npy",
        "sextupole_offsets.npy",
        "bpm_gain_errors.npy",
        "corrector_gain_errors.npy",
        "k2_gain_errors.npy",
        "quadrupole_alignment_standard_normals.npy",
        "forward_drift_response",
        "drift_directions",
    )
    for function in machine_functions:
        source = inspect.getsource(function)
        hits = [token for token in forbidden if token in source]
        require(not hits, f"{function.__name__} leaks {hits}")
    gtpsa_source = (
        STUDY_ROOT
        / "quadrupole_orbit_correction"
        / "gtpsa_noisy_response.jl"
    ).read_text(encoding="utf-8")
    nominal_start = gtpsa_source.index("function nominal_gtpsa_bpm_orm(")
    nominal_end = gtpsa_source.index(
        '"""Apply fixed BPM and corrector gains', nominal_start
    )
    nominal_block = gtpsa_source[nominal_start:nominal_end]
    require(
        "function nominal_gtpsa_bpm_orm(control_names, bpm_names, nominal_orbit)"
        in nominal_block,
        "Nominal GTPSA ORM acquired latent-machine arguments",
    )
    nominal_forbidden = (
        "bpm_gain_errors",
        "physical_corrector_gains",
        "latents",
        "quadrupole",
        "sextupole",
    )
    nominal_hits = [token for token in nominal_forbidden if token in nominal_block]
    require(not nominal_hits, f"Nominal GTPSA ORM leaks {nominal_hits}")
    generator_source = (
        STUDY_ROOT
        / "quadrupole_orbit_correction"
        / "generate_corrected_joint_machine_scans.jl"
    ).read_text(encoding="utf-8")
    defaults_start = generator_source.index(
        "function main_corrected_scans(args=ARGS)"
    )
    defaults_end = generator_source.index(
        "options = parse_exact11_options(defaults, args)", defaults_start
    )
    defaults_block = generator_source[defaults_start:defaults_end]
    require(
        '"baseline-response-method" => "gtpsa"' in defaults_block,
        "Generic correction generator no longer defaults to GTPSA",
    )
    require(
        '"gtpsa-response-model" => "nominal"' in defaults_block,
        "Generic correction generator no longer defaults to the nominal model",
    )
    require(
        '"corrected-case-name" => GTPSA_NOMINAL_CORRECTED_CASE'
        in defaults_block,
        "Generic correction output is mislabeled as the historical FD case",
    )
    require(
        '"correction-bpm-noise-rms-m" => "5.0e-6"' in defaults_block
        and '"correction-measurement-repeats" => "3072"' in defaults_block,
        "Generic correction acquisition no longer matches the maintained case",
    )
    profiled_source = inspect.getsource(analysis.fit_noise_aware_profiled_center)
    require(
        "jac=jacobian" in profiled_source,
        "Profiled optimizer fell back to SciPy numerical differencing",
    )
    complete_source = Path(analysis.__file__).read_text(encoding="utf-8")
    persistence = complete_source.index(
        "# Persist every machine-facing product before evaluation truth is opened."
    )
    truth = complete_source.index(
        'np.load(source / "target_orbits.npy", mmap_mode="r")'
    )
    require(persistence < truth, "Evaluation truth is opened before persistence")
    pre_evaluation_source = complete_source[:persistence]
    pre_evaluation_forbidden = (
        'scan_root / "paired_latents"',
        '"bpm_gain_errors.npy"',
        '"sextupole_offsets.npy"',
    )
    pre_evaluation_hits = [
        token for token in pre_evaluation_forbidden if token in pre_evaluation_source
    ]
    require(
        not pre_evaluation_hits,
        f"The inverse process opens latent artifacts before evaluation: {pre_evaluation_hits}",
    )
    dense_difference = dense_state_space_check()
    profiled_jacobian_difference = analytic_profiled_jacobian_check()

    validation = {
        "status": "PASS",
        "format": result_metadata["format"],
        "case": CASE,
        "machine_count": machine_count,
        "target_count": target_count,
        "bpm_count": bpm_count,
        "stochastic_augmentations": augmentations,
        "filtered_fit_count": int(np.prod(relative.shape[:-1])),
        "nominal_gtpsa_orm_shape": list(orm.shape),
        "finite_difference_orm_reapplied": False,
        "generic_correction_default": "nominal_gtpsa",
        "profiled_optimizer_jacobian": "exact_analytic_variable_projection",
        "profiled_jacobian_fd_audit_relative_difference": profiled_jacobian_difference,
        "thread_worker_count": int(scan["target_scan_worker_count"]),
        "maximum_thread_vs_serial_difference_m": max(float(scan[field]) for field in thread_fields),
        "dense_state_space_max_difference_m": dense_difference,
        "truth_leakage_check": "PASS",
        "observable_bpm_materialization_max_difference_m": observable_difference,
        "drift_observable_bpm_materialization_max_difference_m": drift_observable_difference,
        "zero_state_bpm_closure_max_difference_m": zero_bpm_closure,
        "zero_state_target_closure_max_difference_m": zero_target_closure,
        "unfiltered_bpm_time_error_rms_um": unfiltered_bpm_rms * 1.0e6,
        "filtered_bpm_time_error_rms_um": filtered_bpm_rms * 1.0e6,
        "filtered_relative_center_rmse_um": float(filtered_row["relative_rmse_2d_um"]),
        "filtered_absolute_offset_rmse_um": float(filtered_row["absolute_rmse_2d_um"]),
        "filtered_absolute_offset_p99_um": float(filtered_row["absolute_p99_2d_um"]),
        "aggregate_rmse_below_30_um": float(
            filtered_row["absolute_rmse_2d_um"]
        )
        < 30.0,
        "strict_absolute_p99_below_50_um": float(
            filtered_row["absolute_p99_2d_um"]
        )
        < 50.0,
    }
    (RESULTS / "VALIDATION.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
