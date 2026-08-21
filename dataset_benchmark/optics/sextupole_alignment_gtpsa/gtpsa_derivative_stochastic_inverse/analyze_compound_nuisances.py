#!/usr/bin/env python3
"""Revalidate all maintained nuisances except quadrupole misalignment.

The deterministic states are exact latest-lattice SciBmad scans at the current
eight-state excitation amplitudes.  BPM gain is applied to the simulated
readback, while BPM white noise and continuously evolving drift are propagated
through the maintained covariance-matched and state-space inverse.  Static
nuisance cases use component-matched random seeds, allowing the exact compound
error to be decomposed into a linear sum and a nonlinear interaction residual.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import tomllib
from pathlib import Path

import numpy as np

import analyze_stochastic_inverse as base
import analyze_time_series_inverse as series


HERE = Path(__file__).resolve().parent
STUDY_ROOT = HERE.parent
DEFAULT_SCAN = HERE / "results" / "exact_k5_b3"
DEFAULT_MODEL = STUDY_ROOT / "finite_bpm_inversion" / "results" / "local_orbit_model"

BASELINE = "baseline"
COMPONENT_CASES = (
    "bpm_gain",
    "corrector_gain",
    "k2_calibration",
    "quadrupole_strength",
    "quadrupole_roll",
)
PHYSICAL_COMPONENTS = tuple(case for case in COMPONENT_CASES if case != "bpm_gain")
COMPOUND = "combined_without_quadrupole_misalignment"
COMPOUND_DRIFT = "combined_without_quadrupole_misalignment_time_drift"
LABELS = {
    BASELINE: "Reference: sextupole offsets only",
    "bpm_gain": "BPM gain",
    "corrector_gain": "Corrector gain",
    "k2_calibration": "K2 calibration gain",
    "quadrupole_strength": "Quadrupole strength",
    "quadrupole_roll": "Quadrupole roll",
    "linear_sum_prediction": "Linear sum of paired increments",
    COMPOUND: "All static nuisances combined",
    "baseline_with_white_and_drift": "Reference + white noise + drift inverse",
    "compound_with_white_and_drift": "All nuisances except quadrupole misalignment",
}
MAGNITUDES = {
    BASELINE: "none beyond hidden sextupole offsets",
    "bpm_gain": "1% RMS per BPM/plane, fixed in scan",
    "corrector_gain": "1% RMS per corrector, fixed in scan",
    "k2_calibration": "1% RMS intervention gain per target scan",
    "quadrupole_strength": "independent uniform +/-1%",
    "quadrupole_roll": "1 mrad RMS",
    "linear_sum_prediction": "paired vector sum",
    COMPOUND: "all five static nuisance settings above",
    "baseline_with_white_and_drift": "5 um/read white noise + 10 um endpoint drift",
    "compound_with_white_and_drift": (
        "all five static nuisances + 5 um/read white noise + 10 um endpoint drift"
    ),
}


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_metadata(path: Path) -> dict[str, object]:
    with (path / "scan_metadata.toml").open("rb") as stream:
        return tomllib.load(stream)


def apply_bpm_gain(bpm: np.ndarray, gain_errors: np.ndarray) -> np.ndarray:
    expected = bpm.shape[:2] + bpm.shape[-2:]
    if gain_errors.shape != expected:
        raise ValueError(f"BPM gain shape {gain_errors.shape} != {expected}")
    return np.asarray(bpm) * (1.0 + gain_errors[:, :, None, None, :, :])


def deterministic_inverse(
    source: Path,
    measured_bpm: np.ndarray,
    delta_k2: np.ndarray,
    bump_commands: np.ndarray,
    left_inverses: np.ndarray,
    zero_bump: int,
    zero_k2: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gradients, _amplitude = base.parity_gradients(
        measured_bpm, delta_k2, bump_commands
    )
    right = np.concatenate((gradients[:, :, 0], gradients[:, :, 1]), axis=-1)
    estimates = np.einsum("tic,trc->tri", left_inverses, right)
    target_truth = np.load(source / "target_truth.npy")
    target_orbits = np.asarray(np.load(source / "target_orbits.npy", mmap_mode="r"))
    relative_truth = target_truth - target_orbits[:, :, zero_bump, zero_k2]
    return estimates, relative_truth, estimates - relative_truth


def time_series_covariance(
    static_bpm: np.ndarray,
    drift_bpm: np.ndarray,
    drift_halfwidth_m: float,
    left_inverses: np.ndarray,
    states: list[tuple[int, int, int, int]],
    signal_schedule: list[tuple[int, int]],
    reference_schedule: list[tuple[int, int]],
    reference_bumps: list[int],
    zero_k2: int,
    channel_variance: np.ndarray,
    repeats: int,
    normalization: float,
    endpoint_rms_m: float,
    calibration_reads: int,
    reference_interval: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    response = base.recover_drift_response(
        static_bpm, drift_bpm, drift_halfwidth_m
    )
    target_count, realization_count = static_bpm.shape[:2]
    white_target = series.white_center_covariances(
        left_inverses, states, channel_variance, repeats, normalization
    )
    white = np.broadcast_to(
        white_target[:, None], (target_count, realization_count, 2, 2)
    ).copy()
    vectors = series.core_drift_vectors(
        left_inverses, response, states, repeats, normalization
    )
    reference_noise = series.reference_variances(
        response, reference_bumps, zero_k2, channel_variance
    )
    core_reads = repeats * len(signal_schedule)
    step_variance = endpoint_rms_m**2 / max(core_reads - 1, 1)
    filtered_drift = series.functional_drift_covariance(
        vectors,
        reference_schedule,
        repeats,
        step_variance,
        reference_noise,
        calibration_reads,
        signal_schedule,
        reference_interval,
    ).reshape(target_count, realization_count, 2, 2)
    return white + filtered_drift, white, filtered_drift


def summary_row(
    case: str,
    errors: np.ndarray,
    baseline_errors: np.ndarray,
) -> dict[str, object]:
    sampled = errors if errors.ndim == 4 else errors[None]
    increment = sampled - baseline_errors[None]
    return {
        "case": case,
        "label": LABELS[case],
        "nuisance_magnitude": MAGNITUDES[case],
        "fit_count": int(np.prod(errors.shape[:-1])),
        **base.summarize(errors),
        "increment_vs_clean_reference_rms_um": float(
            np.sqrt(np.mean(np.sum(increment * increment, axis=-1))) * 1e6
        ),
    }


def main() -> int:
    started = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument("--physical-root", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--sextupole-length-m", type=float, default=0.272)
    parser.add_argument("--bpm-gain-rms", type=float, default=0.01)
    parser.add_argument("--bpm-gain-seed", type=int, default=20260829)
    parser.add_argument("--bpm-noise-rms-m", type=float, default=5.0e-6)
    parser.add_argument("--drift-endpoint-rms-m", type=float, default=1.0e-5)
    parser.add_argument("--repeats", type=int, default=3072)
    parser.add_argument("--reference-cycle-interval", type=int, default=256)
    parser.add_argument("--reference-calibration-reads", type=int, default=32)
    parser.add_argument("--monte-carlo-seeds", type=int, default=512)
    parser.add_argument("--measurement-seed", type=int, default=20260921)
    parser.add_argument("--required-rmse-um", type=float, default=50.0)
    parser.add_argument("--preferred-rmse-um", type=float, default=30.0)
    parser.add_argument("--relative-error-scale-um", type=float, default=300.0)
    parser.add_argument(
        "--output-dir", type=Path, default=HERE / "results" / "compound_nuisance_analysis"
    )
    args = parser.parse_args()
    physical_root = args.physical_root.resolve()
    required_directories = (
        BASELINE,
        "time_drift",
        *PHYSICAL_COMPONENTS,
        COMPOUND,
        COMPOUND_DRIFT,
    )
    for case in required_directories:
        if not (physical_root / case / "bpm_orbits.npy").is_file():
            raise FileNotFoundError(f"Missing exact physical scan: {case}")

    baseline_source = physical_root / BASELINE
    metadata = load_metadata(baseline_source)
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta_k2 = levels * float(metadata["k2_step_m3"])
    bump_rows = base.read_rows(baseline_source / "bump_points.csv")
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    zero_bump = int(np.flatnonzero(np.all(bump_commands == 0.0, axis=1))[0])
    zero_k2 = int(np.flatnonzero(levels == 0.0)[0])
    bump_amplitude = float(np.max(np.abs(bump_commands)))
    normalization = float(np.ptp(delta_k2)) * 2.0 * bump_amplitude
    baseline_bpm = np.asarray(
        np.load(baseline_source / "bpm_orbits.npy", mmap_mode="r")
    )
    target_count, realization_count, _nb, _nk, bpm_count, plane_count = baseline_bpm.shape
    channel_count = bpm_count * plane_count
    templates = base.source_templates(args.model_dir.resolve(), args.sextupole_length_m)
    design = base.center_design(templates)
    channel_variance = np.full(channel_count, args.bpm_noise_rms_m**2)
    left_inverses = series.covariance_matched_left_inverses(
        design, channel_variance
    )
    states = base.signed_state_indices(bump_commands, delta_k2)
    signal_schedule = series.core_schedule(len(states))
    reference_schedule, reference_bumps = series.interleaved_reference_schedule(states)

    rng = np.random.default_rng(args.bpm_gain_seed)
    bpm_gain_errors = args.bpm_gain_rms * rng.standard_normal(
        (target_count, realization_count, bpm_count, plane_count)
    )
    deterministic_errors: dict[str, np.ndarray] = {}
    estimates_by_case: dict[str, np.ndarray] = {}
    truths_by_case: dict[str, np.ndarray] = {}
    measured_bpm_by_case: dict[str, np.ndarray] = {BASELINE: baseline_bpm}

    for case in (BASELINE, *PHYSICAL_COMPONENTS):
        source = physical_root / case
        bpm = baseline_bpm if case == BASELINE else np.asarray(
            np.load(source / "bpm_orbits.npy", mmap_mode="r")
        )
        measured_bpm_by_case[case] = bpm
        estimates, truth, errors = deterministic_inverse(
            source,
            bpm,
            delta_k2,
            bump_commands,
            left_inverses,
            zero_bump,
            zero_k2,
        )
        estimates_by_case[case] = estimates
        truths_by_case[case] = truth
        deterministic_errors[case] = errors

    bpm_gain_bpm = apply_bpm_gain(baseline_bpm, bpm_gain_errors)
    bpm_gain_estimates, bpm_gain_truth, bpm_gain_error = deterministic_inverse(
        baseline_source,
        bpm_gain_bpm,
        delta_k2,
        bump_commands,
        left_inverses,
        zero_bump,
        zero_k2,
    )
    estimates_by_case["bpm_gain"] = bpm_gain_estimates
    truths_by_case["bpm_gain"] = bpm_gain_truth
    deterministic_errors["bpm_gain"] = bpm_gain_error

    compound_source = physical_root / COMPOUND
    compound_physical_bpm = np.asarray(
        np.load(compound_source / "bpm_orbits.npy", mmap_mode="r")
    )
    compound_bpm = apply_bpm_gain(compound_physical_bpm, bpm_gain_errors)
    compound_estimates, compound_truth, compound_errors = deterministic_inverse(
        compound_source,
        compound_bpm,
        delta_k2,
        bump_commands,
        left_inverses,
        zero_bump,
        zero_k2,
    )
    estimates_by_case[COMPOUND] = compound_estimates
    truths_by_case[COMPOUND] = compound_truth
    deterministic_errors[COMPOUND] = compound_errors

    baseline_errors = deterministic_errors[BASELINE]
    increments = {
        case: deterministic_errors[case] - baseline_errors for case in COMPONENT_CASES
    }
    linear_sum_errors = baseline_errors + sum(increments.values())
    actual_increment = compound_errors - baseline_errors
    linear_increment = linear_sum_errors - baseline_errors
    interaction = compound_errors - linear_sum_errors
    interaction_summary = {
        "quadrature_component_increment_rms_um": float(
            np.sqrt(
                sum(np.mean(np.sum(value * value, axis=-1)) for value in increments.values())
            )
            * 1e6
        ),
        "linear_vector_sum_increment_rms_um": float(
            np.sqrt(np.mean(np.sum(linear_increment * linear_increment, axis=-1))) * 1e6
        ),
        "actual_compound_increment_rms_um": float(
            np.sqrt(np.mean(np.sum(actual_increment * actual_increment, axis=-1))) * 1e6
        ),
        "nonlinear_interaction_rms_um": float(
            np.sqrt(np.mean(np.sum(interaction * interaction, axis=-1))) * 1e6
        ),
        "nonlinear_interaction_p99_um": float(
            np.percentile(np.linalg.norm(interaction, axis=-1) * 1e6, 99)
        ),
    }
    interaction_norm_um = np.linalg.norm(interaction, axis=-1) * 1e6
    interaction_peak = np.unravel_index(
        int(np.argmax(interaction_norm_um)), interaction_norm_um.shape
    )
    interaction_summary.update(
        {
            "nonlinear_interaction_max_um": float(interaction_norm_um[interaction_peak]),
            "nonlinear_interaction_max_target_index": int(interaction_peak[0] + 1),
            "nonlinear_interaction_max_realization_index": int(interaction_peak[1] + 1),
        }
    )
    interaction_summary["interaction_fraction_of_actual_increment"] = float(
        interaction_summary["nonlinear_interaction_rms_um"]
        / max(interaction_summary["actual_compound_increment_rms_um"], 1.0e-30)
    )

    baseline_drift_bpm = np.asarray(
        np.load(physical_root / "time_drift" / "bpm_orbits.npy", mmap_mode="r")
    )
    baseline_total_cov, baseline_white_cov, baseline_drift_cov = time_series_covariance(
        baseline_bpm,
        baseline_drift_bpm,
        float(load_metadata(physical_root / "time_drift")["drift_halfwidth_m"]),
        left_inverses,
        states,
        signal_schedule,
        reference_schedule,
        reference_bumps,
        zero_k2,
        channel_variance,
        args.repeats,
        normalization,
        args.drift_endpoint_rms_m,
        args.reference_calibration_reads,
        args.reference_cycle_interval,
    )
    compound_drift_physical = np.asarray(
        np.load(physical_root / COMPOUND_DRIFT / "bpm_orbits.npy", mmap_mode="r")
    )
    compound_drift_bpm = apply_bpm_gain(compound_drift_physical, bpm_gain_errors)
    compound_total_cov, compound_white_cov, compound_drift_cov = time_series_covariance(
        compound_bpm,
        compound_drift_bpm,
        float(load_metadata(physical_root / COMPOUND_DRIFT)["drift_halfwidth_m"]),
        left_inverses,
        states,
        signal_schedule,
        reference_schedule,
        reference_bumps,
        zero_k2,
        channel_variance,
        args.repeats,
        normalization,
        args.drift_endpoint_rms_m,
        args.reference_calibration_reads,
        args.reference_cycle_interval,
    )
    baseline_stochastic_rmse_um = float(
        np.sqrt(np.mean(np.trace(baseline_total_cov, axis1=-2, axis2=-1))) * 1e6
    )
    compound_white_rmse_um = float(
        np.sqrt(np.mean(np.trace(compound_white_cov, axis1=-2, axis2=-1))) * 1e6
    )
    compound_filtered_drift_rmse_um = float(
        np.sqrt(np.mean(np.trace(compound_drift_cov, axis1=-2, axis2=-1))) * 1e6
    )
    compound_stochastic_rmse_um = float(
        np.sqrt(np.mean(np.trace(compound_total_cov, axis1=-2, axis2=-1))) * 1e6
    )
    baseline_samples = series.sample_errors(
        baseline_errors,
        baseline_total_cov,
        args.monte_carlo_seeds,
        args.measurement_seed,
    )
    compound_samples = series.sample_errors(
        compound_errors,
        compound_total_cov,
        args.monte_carlo_seeds,
        args.measurement_seed + 1,
    )

    summary_rows: list[dict[str, object]] = []
    for case in (BASELINE, *COMPONENT_CASES):
        summary_rows.append(summary_row(case, deterministic_errors[case], baseline_errors))
    summary_rows.append(
        summary_row("linear_sum_prediction", linear_sum_errors, baseline_errors)
    )
    summary_rows.append(summary_row(COMPOUND, compound_errors, baseline_errors))
    summary_rows.append(
        summary_row(
            "baseline_with_white_and_drift", baseline_samples, baseline_errors
        )
    )
    summary_rows.append(
        summary_row(
            "compound_with_white_and_drift", compound_samples, baseline_errors
        )
    )

    target_names = (baseline_source / "target_names.txt").read_text(
        encoding="utf-8"
    ).splitlines()
    interaction_summary["nonlinear_interaction_max_target"] = target_names[
        interaction_peak[0]
    ]
    per_target_rows: list[dict[str, object]] = []
    for case, values in (
        (COMPOUND, compound_errors[None]),
        ("baseline_with_white_and_drift", baseline_samples),
        ("compound_with_white_and_drift", compound_samples),
    ):
        for target, name in enumerate(target_names):
            per_target_rows.append(
                {
                    "case": case,
                    "target": name,
                    "target_index": target + 1,
                    **base.summarize(values[:, target]),
                }
            )

    stochastic_row = next(
        row for row in summary_rows if row["case"] == "compound_with_white_and_drift"
    )
    stochastic_target_rows = [
        row for row in per_target_rows if row["case"] == "compound_with_white_and_drift"
    ]
    baseline_stochastic_target_rows = [
        row for row in per_target_rows if row["case"] == "baseline_with_white_and_drift"
    ]
    baseline_worst_target_rmse = max(
        float(row["rmse_2d_um"]) for row in baseline_stochastic_target_rows
    )
    worst_target_rmse = max(float(row["rmse_2d_um"]) for row in stochastic_target_rows)
    targets_above_required = [
        str(row["target"])
        for row in stochastic_target_rows
        if float(row["rmse_2d_um"]) >= args.required_rmse_um
    ]
    hard_gate = (
        float(stochastic_row["rmse_2d_um"]) < args.required_rmse_um
        and float(stochastic_row["p99_2d_um"]) < args.required_rmse_um
        and worst_target_rmse < args.required_rmse_um
    )
    preferred_gate = float(stochastic_row["rmse_2d_um"]) < args.preferred_rmse_um
    proxy_relative_rmse = (
        100.0 * float(stochastic_row["rmse_2d_um"]) / args.relative_error_scale_um
    )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "per_target_summary.csv", per_target_rows)
    np.save(output / "latent_bpm_gain_errors.npy", bpm_gain_errors)
    np.save(output / "baseline_deterministic_errors.npy", baseline_errors)
    np.save(output / "compound_deterministic_errors.npy", compound_errors)
    np.save(output / "linear_sum_errors.npy", linear_sum_errors)
    np.save(output / "nonlinear_interaction_errors.npy", interaction)
    np.save(output / "baseline_total_covariances.npy", baseline_total_cov)
    np.save(output / "baseline_white_covariances.npy", baseline_white_cov)
    np.save(output / "baseline_drift_covariances.npy", baseline_drift_cov)
    np.save(output / "compound_total_covariances.npy", compound_total_cov)
    np.save(output / "compound_white_covariances.npy", compound_white_cov)
    np.save(output / "compound_drift_covariances.npy", compound_drift_cov)
    np.savez_compressed(
        output / "center_error_samples.npz",
        baseline_with_white_and_drift=baseline_samples,
        compound_with_white_and_drift=compound_samples,
    )
    (output / "interaction_summary.json").write_text(
        json.dumps(interaction_summary, indent=2) + "\n", encoding="utf-8"
    )

    generation_seconds = sum(
        float(load_metadata(physical_root / case)["calculation_wall_seconds"])
        for case in required_directories
    )
    analysis_seconds = time.perf_counter() - started
    metadata_out = {
        "format": "cesr-eight-state-compound-nuisance-v1",
        "date": "2026-08-20",
        "lattice": metadata["lattice"],
        "engine": "latest-lattice exact SciBmad physical scans plus SciBmad/GTPSA inverse",
        "excluded_nuisance": "quadrupole_misalignment",
        "included_static_nuisances": list(COMPONENT_CASES),
        "target_count": target_count,
        "realizations_per_target": realization_count,
        "exact_states_per_physical_case": int(metadata["total_state_count"]),
        "summed_exact_generation_wall_seconds": generation_seconds,
        "analysis_wall_seconds": analysis_seconds,
        "bump_amplitude_m": bump_amplitude,
        "delta_k2_extrema_m3": [float(delta_k2.min()), float(delta_k2.max())],
        "repeats_per_signal_state": args.repeats,
        "reference_cycle_interval": args.reference_cycle_interval,
        "reference_calibration_reads_per_bump": args.reference_calibration_reads,
        "bpm_gain_rms": args.bpm_gain_rms,
        "bpm_noise_rms_per_read_m": args.bpm_noise_rms_m,
        "core_scan_drift_endpoint_rms_m": args.drift_endpoint_rms_m,
        "monte_carlo_center_draw_count": args.monte_carlo_seeds,
        "compound_combined_rmse_um": float(stochastic_row["rmse_2d_um"]),
        "compound_combined_p99_um": float(stochastic_row["p99_2d_um"]),
        "baseline_worst_target_rmse_um": baseline_worst_target_rmse,
        "compound_worst_target_rmse_um": worst_target_rmse,
        "targets_at_or_above_required_rmse_count": len(targets_above_required),
        "targets_at_or_above_required_rmse": targets_above_required,
        "proxy_relative_rmse_percent": proxy_relative_rmse,
        "baseline_stochastic_component_rmse_um": baseline_stochastic_rmse_um,
        "compound_white_component_rmse_um": compound_white_rmse_um,
        "compound_filtered_drift_component_rmse_um": compound_filtered_drift_rmse_um,
        "compound_total_stochastic_component_rmse_um": compound_stochastic_rmse_um,
        "hard_gate_passed": bool(hard_gate),
        "preferred_gate_passed": bool(preferred_gate),
        **interaction_summary,
    }
    (output / "result_metadata.json").write_text(
        json.dumps(metadata_out, indent=2) + "\n", encoding="utf-8"
    )

    table = "\n".join(
        f"| {row['label']} | {float(row['rmse_2d_um']):.3f} | "
        f"{float(row['increment_vs_clean_reference_rms_um']):.3f} | "
        f"{float(row['p90_2d_um']):.3f} | {float(row['p99_2d_um']):.3f} |"
        for row in summary_rows
    )
    report = f"""# Compound nuisance revalidation with the eight-state inverse

All deterministic states use the current latest-lattice SciBmad protocol:
signed +/-{bump_amplitude * 1e3:.3f} mm local bumps and delta-K2 extrema
{delta_k2.min():.3f}/{delta_k2.max():.3f} m^-3.  The paired benchmark covers
all {target_count} targets and {realization_count} hidden machines per target.
Quadrupole misalignment is deliberately excluded.  The fully combined case
simultaneously activates 1% BPM gain, 1% corrector gain, 1% K2 gain,
independent +/-1% quadrupole strength, 1 mrad RMS quadrupole roll, 5 um/read
BPM white noise, and the maintained 10 um endpoint random-walk drift.

| case | 2D RMSE [um] | increment vs clean reference RMS [um] | P90 [um] | P99 [um] |
|---|---:|---:|---:|---:|
{table}

The component-matched deterministic decomposition gives:

- quadrature of five component increments:
  {interaction_summary['quadrature_component_increment_rms_um']:.3f} um
- vector sum of paired component increments:
  {interaction_summary['linear_vector_sum_increment_rms_um']:.3f} um
- actual compound increment:
  {interaction_summary['actual_compound_increment_rms_um']:.3f} um
- nonlinear compound interaction RMS:
  {interaction_summary['nonlinear_interaction_rms_um']:.3f} um
- nonlinear interaction P99:
  {interaction_summary['nonlinear_interaction_p99_um']:.3f} um
- nonlinear interaction maximum:
  {interaction_summary['nonlinear_interaction_max_um']:.3f} um at
  {interaction_summary['nonlinear_interaction_max_target']}, realization
  {interaction_summary['nonlinear_interaction_max_realization_index']}
- nonlinear interaction / actual increment:
  {100 * interaction_summary['interaction_fraction_of_actual_increment']:.3f}%

For the complete stochastic case:

- aggregate 2D RMSE: {float(stochastic_row['rmse_2d_um']):.3f} um
- P99: {float(stochastic_row['p99_2d_um']):.3f} um
- worst target-level RMSE: {worst_target_rmse:.3f} um
- targets at or above {args.required_rmse_um:.0f} um RMSE:
  {len(targets_above_required)} ({', '.join(targets_above_required) or 'none'})
- compound white-noise component: {compound_white_rmse_um:.3f} um
- compound filtered-drift component: {compound_filtered_drift_rmse_um:.3f} um
- compound total stochastic component: {compound_stochastic_rmse_um:.3f} um
- clean-reference total stochastic component: {baseline_stochastic_rmse_um:.3f} um
- proxy relative RMSE on a {args.relative_error_scale_um:.0f} um scale:
  {proxy_relative_rmse:.3f}%
- hard aggregate/P99/all-target 50 um gate: {'PASS' if hard_gate else 'FAIL'}
- preferred aggregate 30 um gate: {'PASS' if preferred_gate else 'FAIL'}

BPM gain is applied to both signal and reference readbacks but is not supplied
to the nominal center template.  Corrector and K2 gains likewise alter the
physical SciBmad scan while the inverse uses commanded spans.  The compound
drift response is independently recovered from a paired exact SciBmad secant
with identical static nuisance realizations.  White noise and random-walk
histories are then propagated through the repeated acquisition sequence.

These remain synthetic sensitivities, not measured CESR error priors.  The
drift basis is one calibrated scalar mode, BPM white noise is temporally
independent, actuator hysteresis and polarity asymmetry are absent, and the
nuisance ensemble has only four hidden machines per target.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
