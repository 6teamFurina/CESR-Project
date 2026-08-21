#!/usr/bin/env python3
"""Independent checks for the eight-state time-series drift inverse."""

from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

import numpy as np

import analyze_stochastic_inverse as base
import analyze_time_series_inverse as series


HERE = Path(__file__).resolve().parent
PHYSICAL_ROOT = HERE / "results" / "exact_k5_b3"
MODEL_DIR = HERE.parent / "finite_bpm_inversion" / "results" / "local_orbit_model"
OUTPUT_DIR = HERE / "results" / "time_series_analysis"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def dense_functional_covariance(
    core_vectors: np.ndarray,
    schedule: list[tuple[int, int]],
    repeats: int,
    step_variance: float,
    reference_variances: np.ndarray | None = None,
    calibration_reads: int = 0,
    signal_only_schedule: list[tuple[int, int]] | None = None,
    reference_cycle_interval: int = 1,
) -> np.ndarray:
    """Conventional dense Kalman covariance for one target/realization."""
    reference_count = 0 if reference_variances is None else len(reference_variances)
    dimension = 3 + reference_count
    covariance = np.zeros((dimension, dimension))
    if reference_count:
        diagonal = np.arange(reference_count) + 3
        covariance[diagonal, diagonal] = reference_variances / calibration_reads
    time_index = 0
    for cycle in range(repeats):
        if signal_only_schedule is not None and reference_cycle_interval > 1:
            is_reference_cycle = (
                cycle % reference_cycle_interval == 0 or cycle == repeats - 1
            )
            cycle_schedule = schedule if is_reference_cycle else signal_only_schedule
        else:
            cycle_schedule = schedule
        for core, reference in cycle_schedule:
            if time_index > 0:
                covariance[0, 0] += step_variance
            if core >= 0:
                transition = np.eye(dimension)
                transition[1:3, 0] = core_vectors[core]
                covariance = transition @ covariance @ transition.T
            if reference >= 0 and reference_count:
                observation = np.zeros(dimension)
                observation[0] = 1.0
                observation[3 + reference] = 1.0
                innovation = float(
                    observation @ covariance @ observation
                    + reference_variances[reference]
                )
                gain_numerator = covariance @ observation
                covariance -= np.outer(gain_numerator, gain_numerator) / innovation
            time_index += 1
    return covariance[1:3, 1:3]


def main() -> int:
    baseline_dir = PHYSICAL_ROOT / "baseline"
    drift_dir = PHYSICAL_ROOT / "time_drift"
    with (baseline_dir / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    with (drift_dir / "scan_metadata.toml").open("rb") as stream:
        drift_metadata = tomllib.load(stream)
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta_k2 = levels * float(metadata["k2_step_m3"])
    bump_rows = base.read_rows(baseline_dir / "bump_points.csv")
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    bpm = np.asarray(np.load(baseline_dir / "bpm_orbits.npy", mmap_mode="r"))
    drift_scan = np.asarray(np.load(drift_dir / "bpm_orbits.npy", mmap_mode="r"))
    response = base.recover_drift_response(
        bpm, drift_scan, float(drift_metadata["drift_halfwidth_m"])
    )
    templates = base.source_templates(MODEL_DIR, 0.272)
    design = base.center_design(templates)
    channel_count = response.shape[-1]
    noise_rms = 5.0e-6
    channel_variance = np.full(channel_count, noise_rms**2)
    left = series.covariance_matched_left_inverses(design, channel_variance)
    np.testing.assert_allclose(
        left,
        np.asarray([np.linalg.pinv(matrix, rcond=1.0e-12) for matrix in design]),
        rtol=2.0e-13,
        atol=1.0e-15,
    )

    states = base.signed_state_indices(bump_commands, delta_k2)
    amplitude = float(np.max(np.abs(bump_commands)))
    normalization = float(np.ptp(delta_k2)) * 2.0 * amplitude
    short_repeats = 3
    white = series.white_center_covariances(
        left, states, channel_variance, short_repeats, normalization
    )
    old_white = base.white_center_covariances(
        design,
        noise_rms,
        short_repeats,
        float(np.ptp(delta_k2)),
        2.0 * amplitude,
    )
    np.testing.assert_allclose(white, old_white, rtol=2.0e-13, atol=1.0e-30)

    vectors = series.core_drift_vectors(
        left, response, states, short_repeats, normalization
    )
    step_variance = (1.0e-5) ** 2 / (short_repeats * len(states) - 1)
    compact_balanced = series.functional_drift_covariance(
        vectors,
        series.core_schedule(len(states)),
        short_repeats,
        step_variance,
    ).reshape(76, 4, 2, 2)
    old_balanced = base.random_walk_center_covariances(
        design,
        response,
        states,
        short_repeats,
        float(np.ptp(delta_k2)),
        2.0 * amplitude,
        1.0e-5,
    )
    np.testing.assert_allclose(
        compact_balanced, old_balanced, rtol=3.0e-12, atol=2.0e-26
    )

    reference_schedule, reference_bumps = series.interleaved_reference_schedule(states)
    zero_k2 = int(np.flatnonzero(levels == 0.0)[0])
    reference_noise = series.reference_variances(
        response, reference_bumps, zero_k2, channel_variance
    )
    compact_filtered = series.functional_drift_covariance(
        vectors,
        reference_schedule,
        short_repeats,
        step_variance,
        reference_noise,
        32,
        series.core_schedule(len(states)),
        2,
    )
    for target, realization in ((0, 0), (13, 2), (43, 1), (75, 3)):
        case = target * 4 + realization
        dense = dense_functional_covariance(
            vectors[:, case],
            reference_schedule,
            short_repeats,
            step_variance,
            reference_noise[case],
            32,
            series.core_schedule(len(states)),
            2,
        )
        np.testing.assert_allclose(
            compact_filtered[case], dense, rtol=3.0e-12, atol=2.0e-26
        )

    with (OUTPUT_DIR / "result_metadata.json").open(encoding="utf-8") as stream:
        result = json.load(stream)
    require(result["target_count"] == 76, "Time-series result does not cover 76 targets")
    require(result["signal_state_count"] == 8, "Signal protocol is not eight-state")
    require(result["reference_cycle_interval"] == 256, "Reference cadence changed")
    require(
        result["reference_calibration_reads_per_bump"] == 32,
        "Reference calibration count changed",
    )
    require(
        result["lattice"] == metadata["lattice"],
        "Result lattice does not match the exact SciBmad scan",
    )
    require(result["hard_gate_passed"] is True, "Saved hard 50 um gate failed")
    require(result["preferred_gate_passed"] is True, "Saved preferred 30 um gate failed")
    require(float(result["filtered_combined_rmse_um"]) < 30.0, "Combined RMSE >= 30 um")
    require(
        float(result["proxy_relative_rmse_percent"]) < 10.0,
        "Combined RMSE is not below 10% of the requested 300 um scale",
    )
    require(float(result["filtered_combined_p99_um"]) < 50.0, "Combined P99 >= 50 um")
    require(
        float(result["filtered_worst_target_rmse_um"]) < 50.0,
        "A target-level combined RMSE is >= 50 um",
    )

    summaries = {row["case"]: row for row in rows(OUTPUT_DIR / "summary.csv")}
    require(
        float(summaries["reference_filtered_combined"]["rmse_2d_um"])
        < float(summaries["balanced_8state_combined"]["rmse_2d_um"]),
        "Drift filtering did not reduce the combined RMSE",
    )
    schedule_rows = rows(OUTPUT_DIR / "protocol_schedule.csv")
    require(len(schedule_rows) == 20, "Expected a 20-acquisition interleaved cycle")
    require(
        sum(row["kind"] == "signal" for row in schedule_rows) == 8,
        "Interleaved schedule does not retain eight core states",
    )
    require(
        sum(row["kind"] == "K2_zero_reference" for row in schedule_rows) == 12,
        "Interleaved schedule does not contain twelve K2=0 references",
    )

    saved_white = np.load(OUTPUT_DIR / "white_center_covariances.npy")
    saved_balanced = np.load(OUTPUT_DIR / "balanced_drift_covariances.npy")
    saved_filtered = np.load(OUTPUT_DIR / "filtered_drift_covariances.npy")
    require(saved_white.shape == (76, 4, 2, 2), "Saved white covariance shape mismatch")
    require(saved_balanced.shape == saved_white.shape, "Balanced covariance shape mismatch")
    require(saved_filtered.shape == saved_white.shape, "Filtered covariance shape mismatch")
    for name, covariance in (
        ("white", saved_white),
        ("balanced", saved_balanced),
        ("filtered", saved_filtered),
    ):
        eigenvalues = np.linalg.eigvalsh(covariance)
        require(np.min(eigenvalues) >= -1.0e-24, f"{name} covariance is not PSD")

    print("PASS: iid covariance-matched GLS reduces exactly to the maintained OLS fit")
    print("PASS: white covariance matches the maintained analytic result")
    print("PASS: compact core-only recursion matches explicit read-by-read drift")
    print("PASS: partitioned reference filter matches an independent dense Kalman recursion")
    print("PASS: all 76 targets satisfy the saved 50 um gate and aggregate RMSE is below 30 um")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
