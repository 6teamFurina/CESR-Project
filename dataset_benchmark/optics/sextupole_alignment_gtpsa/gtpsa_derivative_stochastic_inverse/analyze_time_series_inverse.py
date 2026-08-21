#!/usr/bin/env python3
"""Eight-state time-series inverse with BPM whitening and drift filtering.

The exact SciBmad scan tensor supplies the deterministic orbit at the eight
signed signal states.  This analysis promotes the measurement model from
eight state averages to a repeated, time-ordered acquisition stream.  The
signal states remain the same K2-odd/bump-odd eight-state protocol.  Optional
same-bump K2=0 references are interleaved as ``0,+,0,-,0`` and a linear
Gaussian state-space filter estimates the random-walk drift contribution to
the final two-parameter center estimator.

The full BPM stream is not materialized.  Since the estimator, white noise,
and random walk are linear Gaussian, the code propagates the exact sufficient
statistics: the covariance-matched full-BPM center weights and a three-state
``[drift, accumulated_x_error, accumulated_y_error]`` functional filter.  A
finite calibration uncertainty for the four reference baselines is retained.
The companion validator compares this reduced recursion with a dense Kalman
covariance calculation and the existing explicit read-by-read random walk.
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


HERE = Path(__file__).resolve().parent
STUDY_ROOT = HERE.parent
DEFAULT_SCAN = HERE / "results" / "exact_k5_b3"
DEFAULT_MODEL = STUDY_ROOT / "finite_bpm_inversion" / "results" / "local_orbit_model"


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def covariance_matched_left_inverses(
    design: np.ndarray,
    channel_covariance: np.ndarray,
) -> np.ndarray:
    """Return GLS/matched-filter left inverses for every target.

    ``design`` stacks the x-bump and y-bump gradient blocks.  The same
    full-BPM covariance is used in both blocks.  A one-dimensional covariance
    argument denotes channel variances and avoids forming a dense inverse.
    """
    target_count, stacked_channels, parameter_count = design.shape
    if parameter_count != 2 or stacked_channels % 2:
        raise ValueError(f"Unexpected center design shape: {design.shape}")
    channel_count = stacked_channels // 2
    covariance = np.asarray(channel_covariance, dtype=float)
    if covariance.shape == (channel_count,):
        if np.any(covariance <= 0.0):
            raise ValueError("Every BPM channel variance must be positive")

        def precision_apply(matrix: np.ndarray) -> np.ndarray:
            return matrix / covariance[:, None]

    elif covariance.shape == (channel_count, channel_count):
        np.linalg.cholesky(covariance)

        def precision_apply(matrix: np.ndarray) -> np.ndarray:
            return np.linalg.solve(covariance, matrix)

    else:
        raise ValueError(
            f"Expected {channel_count} variances or a square covariance, got {covariance.shape}"
        )

    result = np.empty((target_count, 2, stacked_channels))
    for target in range(target_count):
        precision_design = np.concatenate(
            (
                precision_apply(design[target, :channel_count]),
                precision_apply(design[target, channel_count:]),
            ),
            axis=0,
        )
        information = design[target].T @ precision_design
        result[target] = np.linalg.solve(information, precision_design.T)
    return result


def apply_channel_covariance(
    weights: np.ndarray,
    channel_covariance: np.ndarray,
) -> np.ndarray:
    """Return W C W' for target-indexed 2 x channel weights."""
    covariance = np.asarray(channel_covariance, dtype=float)
    if covariance.ndim == 1:
        return np.einsum("tim,m,tjm->tij", weights, covariance, weights)
    return np.einsum("tim,mn,tjn->tij", weights, covariance, weights)


def white_center_covariances(
    left_inverses: np.ndarray,
    states: list[tuple[int, int, int, int]],
    channel_covariance: np.ndarray,
    repeats: int,
    normalization: float,
) -> np.ndarray:
    """Propagate per-read BPM covariance through the eight-state matched filter."""
    if repeats <= 0 or normalization <= 0.0:
        raise ValueError("Repeats and contrast normalization must be positive")
    channel_count = left_inverses.shape[-1] // 2
    covariance = np.zeros((left_inverses.shape[0], 2, 2))
    for block, sign, _bump, _k2 in states:
        block_weights = (
            sign
            * left_inverses[
                :, :, block * channel_count : (block + 1) * channel_count
            ]
            / normalization
        )
        covariance += apply_channel_covariance(block_weights, channel_covariance)
    return covariance / repeats


def core_drift_vectors(
    left_inverses: np.ndarray,
    drift_response: np.ndarray,
    states: list[tuple[int, int, int, int]],
    repeats: int,
    normalization: float,
) -> np.ndarray:
    """Map scalar drift at each signed read to final x/y center error."""
    target_count, realization_count, _nb, _nk, channel_count = drift_response.shape
    if left_inverses.shape != (target_count, 2, 2 * channel_count):
        raise ValueError("Left inverse and drift response shapes are inconsistent")
    vectors = np.empty((len(states), target_count * realization_count, 2))
    for state_index, (block, sign, bump, k2) in enumerate(states):
        weights = left_inverses[
            :, :, block * channel_count : (block + 1) * channel_count
        ]
        response = drift_response[:, :, bump, k2]
        projected = np.einsum("tim,trm->tri", weights, response)
        vectors[state_index] = (
            sign * projected.reshape(target_count * realization_count, 2)
            / (normalization * repeats)
        )
    return vectors


def reference_variances(
    drift_response: np.ndarray,
    reference_bumps: list[int],
    zero_k2: int,
    channel_covariance: np.ndarray,
) -> np.ndarray:
    """Matched-filter variance of one scalar drift observation per reference."""
    covariance = np.asarray(channel_covariance, dtype=float)
    target_count, realization_count = drift_response.shape[:2]
    result = np.empty((target_count * realization_count, len(reference_bumps)))
    for reference, bump in enumerate(reference_bumps):
        response = drift_response[:, :, bump, zero_k2]
        flat = response.reshape(target_count * realization_count, response.shape[-1])
        if covariance.ndim == 1:
            information = np.sum(flat * flat / covariance[None], axis=-1)
        else:
            solved = np.linalg.solve(covariance, flat.T).T
            information = np.sum(flat * solved, axis=-1)
        if np.any(information <= 0.0):
            raise ValueError("A reference drift template has zero information")
        result[:, reference] = 1.0 / information
    return result


def core_schedule(state_count: int) -> list[tuple[int, int]]:
    """Return ``(core state, reference type)`` entries for the eight signals."""
    return [(state, -1) for state in range(state_count)]


def interleaved_reference_schedule(
    states: list[tuple[int, int, int, int]],
) -> tuple[list[tuple[int, int]], list[int]]:
    """Build four same-bump ``0,+,0,-,0`` blocks around the eight core states."""
    if len(states) != 8:
        raise ValueError("The interleaved protocol requires exactly eight signal states")
    schedule: list[tuple[int, int]] = []
    reference_bumps: list[int] = []
    for pair_start in range(0, len(states), 2):
        first = states[pair_start]
        second = states[pair_start + 1]
        if first[2] != second[2]:
            raise ValueError("Each K+/K- signal pair must keep the same bump")
        reference = len(reference_bumps)
        reference_bumps.append(first[2])
        schedule.extend(
            (
                (-1, reference),
                (pair_start, -1),
                (-1, reference),
                (pair_start + 1, -1),
                (-1, reference),
            )
        )
    return schedule, reference_bumps


def functional_drift_covariance(
    core_vectors: np.ndarray,
    schedule: list[tuple[int, int]],
    repeats: int,
    step_variance: float,
    reference_noise_variance: np.ndarray | None = None,
    reference_calibration_reads: int = 0,
    signal_only_schedule: list[tuple[int, int]] | None = None,
    reference_cycle_interval: int = 1,
) -> np.ndarray:
    """Exact covariance of the residual drift contribution after filtering.

    The latent drift is a scalar random walk.  The filter state contains its
    current value, the two accumulated center-error functionals, and (when
    references are active) one static calibration error per reference bump.
    Later reference observations update the accumulated functionals through
    their covariance with the current drift, which is the state-space
    equivalent of smoothing the complete time series for the final center.
    """
    if repeats <= 0 or step_variance < 0.0:
        raise ValueError("Invalid repeat count or random-walk step variance")
    if reference_cycle_interval <= 0:
        raise ValueError("Reference cycle interval must be positive")
    case_count = core_vectors.shape[1]
    if core_vectors.shape[2] != 2:
        raise ValueError("Core drift vectors must end in two center components")
    use_references = reference_noise_variance is not None
    if use_references:
        reference_noise_variance = np.asarray(reference_noise_variance, dtype=float)
        if reference_noise_variance.shape[0] != case_count:
            raise ValueError("Reference variance case count mismatch")
        reference_count = reference_noise_variance.shape[1]
        if reference_calibration_reads <= 0:
            raise ValueError("Reference calibration reads must be positive")
    else:
        reference_count = 0

    # Partitioned covariance for [q, z_x, z_y, a_1, ..., a_n].
    pqq = np.zeros(case_count)
    pqz = np.zeros((case_count, 2))
    pzz = np.zeros((case_count, 2, 2))
    pqa = np.zeros((case_count, reference_count))
    pza = np.zeros((case_count, 2, reference_count))
    paa = np.zeros((case_count, reference_count, reference_count))
    if use_references:
        diagonal = np.arange(reference_count)
        paa[:, diagonal, diagonal] = (
            reference_noise_variance / reference_calibration_reads
        )

    time_index = 0
    for cycle in range(repeats):
        use_reference_cycle = (
            use_references
            and signal_only_schedule is not None
            and reference_cycle_interval > 1
        )
        if use_reference_cycle:
            is_reference_cycle = (
                cycle % reference_cycle_interval == 0 or cycle == repeats - 1
            )
            cycle_schedule = schedule if is_reference_cycle else signal_only_schedule
        else:
            cycle_schedule = schedule
        for core, reference in cycle_schedule:
            if time_index > 0:
                pqq += step_variance

            if core >= 0:
                vector = core_vectors[core]
                old_pqz = pqz.copy()
                pzz += (
                    vector[:, :, None] * old_pqz[:, None, :]
                    + old_pqz[:, :, None] * vector[:, None, :]
                    + pqq[:, None, None]
                    * vector[:, :, None]
                    * vector[:, None, :]
                )
                if reference_count:
                    pza += vector[:, :, None] * pqa[:, None, :]
                pqz += pqq[:, None] * vector

            if reference >= 0 and use_references:
                if reference >= reference_count:
                    raise ValueError("Schedule references an unknown baseline state")
                innovation_variance = (
                    pqq
                    + 2.0 * pqa[:, reference]
                    + paa[:, reference, reference]
                    + reference_noise_variance[:, reference]
                )
                cov_q = pqq + pqa[:, reference]
                cov_z = pqz + pza[:, :, reference]
                cov_a = pqa + paa[:, :, reference]
                inverse_innovation = 1.0 / innovation_variance
                pqq -= cov_q * cov_q * inverse_innovation
                pqz -= cov_q[:, None] * cov_z * inverse_innovation[:, None]
                pqa -= cov_q[:, None] * cov_a * inverse_innovation[:, None]
                pzz -= (
                    cov_z[:, :, None]
                    * cov_z[:, None, :]
                    * inverse_innovation[:, None, None]
                )
                pza -= (
                    cov_z[:, :, None]
                    * cov_a[:, None, :]
                    * inverse_innovation[:, None, None]
                )
                paa -= (
                    cov_a[:, :, None]
                    * cov_a[:, None, :]
                    * inverse_innovation[:, None, None]
                )
            time_index += 1

    return 0.5 * (pzz + np.swapaxes(pzz, -1, -2))


def reference_cycle_count(repeats: int, interval: int) -> int:
    """Number of periodic reference cycles, including a final endpoint anchor."""
    if repeats <= 0 or interval <= 0:
        raise ValueError("Repeats and reference interval must be positive")
    count = (repeats - 1) // interval + 1
    if (repeats - 1) % interval:
        count += 1
    return count


def expected_metrics(
    clean_errors: np.ndarray,
    covariance: np.ndarray,
) -> dict[str, float]:
    """Return aggregate and target-level expected RMSE for fixed clean bias."""
    squared = np.sum(clean_errors * clean_errors, axis=-1) + np.trace(
        covariance, axis1=-2, axis2=-1
    )
    return {
        "rmse_2d_um": float(np.sqrt(np.mean(squared)) * 1e6),
        "worst_target_rmse_2d_um": float(
            np.max(np.sqrt(np.mean(squared, axis=1))) * 1e6
        ),
        "stochastic_component_rmse_2d_um": float(
            np.sqrt(np.mean(np.trace(covariance, axis1=-2, axis2=-1))) * 1e6
        ),
    }


def sample_errors(
    clean_errors: np.ndarray,
    covariance: np.ndarray,
    sample_count: int,
    seed: int,
) -> np.ndarray:
    return clean_errors[None] + base.covariance_samples(
        covariance, sample_count, np.random.default_rng(seed)
    )


def summarize_case(case: str, errors: np.ndarray) -> dict[str, object]:
    return {
        "case": case,
        "fit_count": int(np.prod(errors.shape[:-1])),
        **base.summarize(errors),
    }


def protocol_rows(
    schedule: list[tuple[int, int]],
    states: list[tuple[int, int, int, int]],
    reference_bumps: list[int],
    bump_commands: np.ndarray,
    delta_k2: np.ndarray,
    zero_k2: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for slot, (core, reference) in enumerate(schedule, start=1):
        if core >= 0:
            block, sign, bump, k2 = states[core]
            kind = "signal"
            reference_label = ""
            core_label: object = core + 1
        else:
            block, sign = -1, 0
            bump = reference_bumps[reference]
            k2 = zero_k2
            kind = "K2_zero_reference"
            reference_label = reference + 1
            core_label = ""
        rows.append(
            {
                "slot_in_cycle": slot,
                "kind": kind,
                "core_state": core_label,
                "reference_type": reference_label,
                "gradient_block": "x" if block == 0 else ("y" if block == 1 else ""),
                "contrast_sign": sign,
                "bump_x_command_m": bump_commands[bump, 0],
                "bump_y_command_m": bump_commands[bump, 1],
                "delta_k2_m3": delta_k2[k2],
            }
        )
    return rows


def main() -> int:
    started = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument("--physical-root", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--sextupole-length-m", type=float, default=0.272)
    parser.add_argument("--repeats", type=int, default=3072)
    parser.add_argument(
        "--candidate-repeats", default="256,512,1024,1280,2048,3072,4096"
    )
    parser.add_argument("--bpm-noise-rms-m", type=float, default=5.0e-6)
    parser.add_argument("--core-scan-drift-endpoint-rms-m", type=float, default=1.0e-5)
    parser.add_argument("--reference-calibration-reads", type=int, default=32)
    parser.add_argument("--reference-cycle-interval", type=int, default=256)
    parser.add_argument("--monte-carlo-seeds", type=int, default=512)
    parser.add_argument("--measurement-seed", type=int, default=20260911)
    parser.add_argument("--required-rmse-um", type=float, default=50.0)
    parser.add_argument("--preferred-rmse-um", type=float, default=30.0)
    parser.add_argument("--relative-error-scale-um", type=float, default=300.0)
    parser.add_argument(
        "--output-dir", type=Path, default=HERE / "results" / "time_series_analysis"
    )
    args = parser.parse_args()
    if (
        args.repeats <= 0
        or args.reference_calibration_reads <= 0
        or args.reference_cycle_interval <= 0
        or args.relative_error_scale_um <= 0.0
    ):
        raise ValueError("Repeat and reference calibration counts must be positive")

    physical_root = args.physical_root.resolve()
    baseline_dir = physical_root / "baseline"
    drift_dir = physical_root / "time_drift"
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
    gradients, bump_amplitude = base.parity_gradients(bpm, delta_k2, bump_commands)
    templates = base.source_templates(args.model_dir.resolve(), args.sextupole_length_m)
    design = base.center_design(templates)
    target_count, realization_count = bpm.shape[:2]
    channel_count = bpm.shape[-2] * bpm.shape[-1]
    if templates.shape[:2] != (target_count, channel_count):
        raise ValueError("SciBmad scan and GTPSA template channel counts differ")

    channel_variance = np.full(channel_count, args.bpm_noise_rms_m**2)
    left_inverses = covariance_matched_left_inverses(design, channel_variance)
    right = np.concatenate((gradients[:, :, 0], gradients[:, :, 1]), axis=-1)
    clean_estimates = np.einsum("tic,trc->tri", left_inverses, right)
    zero_bump = int(np.flatnonzero(np.all(bump_commands == 0.0, axis=1))[0])
    zero_k2 = int(np.flatnonzero(levels == 0.0)[0])
    target_truth = np.load(baseline_dir / "target_truth.npy")
    target_orbits = np.asarray(np.load(baseline_dir / "target_orbits.npy", mmap_mode="r"))
    relative_truth = target_truth - target_orbits[:, :, zero_bump, zero_k2]
    clean_errors = clean_estimates - relative_truth

    drift_scan = np.asarray(np.load(drift_dir / "bpm_orbits.npy", mmap_mode="r"))
    drift_response = base.recover_drift_response(
        bpm, drift_scan, float(drift_metadata["drift_halfwidth_m"])
    )
    states = base.signed_state_indices(bump_commands, delta_k2)
    signal_schedule = core_schedule(len(states))
    reference_schedule, reference_bumps = interleaved_reference_schedule(states)
    reference_noise = reference_variances(
        drift_response,
        reference_bumps,
        zero_k2,
        channel_variance,
    )
    normalization = float(np.ptp(delta_k2)) * 2.0 * bump_amplitude
    case_shape = (target_count, realization_count, 2, 2)

    candidate_repeats = sorted(
        {int(value) for value in args.candidate_repeats.split(",")} | {args.repeats}
    )
    if any(value <= 0 for value in candidate_repeats):
        raise ValueError("Every candidate repeat count must be positive")
    tradeoff_rows: list[dict[str, object]] = []
    selected: dict[str, np.ndarray] = {}
    for repeats in candidate_repeats:
        white_target = white_center_covariances(
            left_inverses, states, channel_variance, repeats, normalization
        )
        white = np.broadcast_to(
            white_target[:, None], case_shape
        ).copy()
        vectors = core_drift_vectors(
            left_inverses, drift_response, states, repeats, normalization
        )
        baseline_read_count = repeats * len(signal_schedule)
        step_variance = args.core_scan_drift_endpoint_rms_m**2 / max(
            baseline_read_count - 1, 1
        )
        balanced_drift = functional_drift_covariance(
            vectors,
            signal_schedule,
            repeats,
            step_variance,
        ).reshape(case_shape)
        filtered_drift = functional_drift_covariance(
            vectors,
            reference_schedule,
            repeats,
            step_variance,
            reference_noise,
            args.reference_calibration_reads,
            signal_schedule,
            args.reference_cycle_interval,
        ).reshape(case_shape)
        balanced_metrics = expected_metrics(clean_errors, white + balanced_drift)
        filtered_metrics = expected_metrics(clean_errors, white + filtered_drift)
        filtered_drift_rms = float(
            np.sqrt(np.mean(np.trace(filtered_drift, axis1=-2, axis2=-1))) * 1e6
        )
        balanced_drift_rms = float(
            np.sqrt(np.mean(np.trace(balanced_drift, axis1=-2, axis2=-1))) * 1e6
        )
        reference_cycles = reference_cycle_count(
            repeats, args.reference_cycle_interval
        )
        reference_read_count = (
            repeats * len(signal_schedule)
            + reference_cycles
            * (len(reference_schedule) - len(signal_schedule))
        )
        reference_endpoint = args.core_scan_drift_endpoint_rms_m * np.sqrt(
            max(reference_read_count - 1, 0) / max(baseline_read_count - 1, 1)
        )
        tradeoff_rows.append(
            {
                "repeats_per_signal_state": repeats,
                "balanced_signal_acquisitions": baseline_read_count,
                "reference_protocol_acquisitions": reference_read_count,
                "interleaved_reference_cycle_count": reference_cycles,
                "reference_cycle_interval": args.reference_cycle_interval,
                "reference_calibration_acquisitions": (
                    len(reference_bumps) * args.reference_calibration_reads
                ),
                "balanced_drift_component_rmse_um": balanced_drift_rms,
                "filtered_drift_component_rmse_um": filtered_drift_rms,
                "balanced_combined_rmse_um": balanced_metrics["rmse_2d_um"],
                "filtered_combined_rmse_um": filtered_metrics["rmse_2d_um"],
                "balanced_worst_target_rmse_um": balanced_metrics[
                    "worst_target_rmse_2d_um"
                ],
                "filtered_worst_target_rmse_um": filtered_metrics[
                    "worst_target_rmse_2d_um"
                ],
                "reference_scan_endpoint_drift_rms_um": reference_endpoint * 1e6,
            }
        )
        print(
            f"R={repeats:5d} balanced={balanced_metrics['rmse_2d_um']:8.3f} um "
            f"filtered={filtered_metrics['rmse_2d_um']:8.3f} um "
            f"drift {balanced_drift_rms:7.3f}->{filtered_drift_rms:7.3f} um"
        )
        if repeats == args.repeats:
            selected = {
                "white": white,
                "balanced_drift": balanced_drift,
                "filtered_drift": filtered_drift,
            }

    if not selected:
        raise AssertionError("Selected repeat count was not evaluated")
    balanced_covariance = selected["white"] + selected["balanced_drift"]
    filtered_covariance = selected["white"] + selected["filtered_drift"]
    white_errors = sample_errors(
        clean_errors,
        selected["white"],
        args.monte_carlo_seeds,
        args.measurement_seed,
    )
    balanced_errors = sample_errors(
        clean_errors,
        balanced_covariance,
        args.monte_carlo_seeds,
        args.measurement_seed + 1,
    )
    filtered_drift_errors = sample_errors(
        clean_errors,
        selected["filtered_drift"],
        args.monte_carlo_seeds,
        args.measurement_seed + 2,
    )
    filtered_errors = sample_errors(
        clean_errors,
        filtered_covariance,
        args.monte_carlo_seeds,
        args.measurement_seed + 3,
    )
    cases = {
        "clean": clean_errors,
        "bpm_white_noise_matched_filter": white_errors,
        "balanced_8state_combined": balanced_errors,
        "reference_filtered_drift": filtered_drift_errors,
        "reference_filtered_combined": filtered_errors,
    }
    summary_rows = [summarize_case(case, errors) for case, errors in cases.items()]
    target_names = (baseline_dir / "target_names.txt").read_text(
        encoding="utf-8"
    ).splitlines()
    per_target_rows: list[dict[str, object]] = []
    for case in ("balanced_8state_combined", "reference_filtered_combined"):
        for target, name in enumerate(target_names):
            per_target_rows.append(
                {
                    "case": case,
                    "target": name,
                    "target_index": target + 1,
                    **base.summarize(cases[case][:, target]),
                }
            )

    filtered_row = next(
        row for row in summary_rows if row["case"] == "reference_filtered_combined"
    )
    balanced_row = next(
        row for row in summary_rows if row["case"] == "balanced_8state_combined"
    )
    filtered_targets = [
        row for row in per_target_rows if row["case"] == "reference_filtered_combined"
    ]
    worst_target = max(float(row["rmse_2d_um"]) for row in filtered_targets)
    hard_gate = (
        float(filtered_row["rmse_2d_um"]) < args.required_rmse_um
        and float(filtered_row["p99_2d_um"]) < args.required_rmse_um
        and worst_target < args.required_rmse_um
    )
    preferred_gate = float(filtered_row["rmse_2d_um"]) < args.preferred_rmse_um
    proxy_relative_rmse_percent = (
        100.0 * float(filtered_row["rmse_2d_um"]) / args.relative_error_scale_um
    )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "protocol_tradeoff.csv", tradeoff_rows)
    write_rows(output / "per_target_summary.csv", per_target_rows)
    write_rows(
        output / "protocol_schedule.csv",
        protocol_rows(
            reference_schedule,
            states,
            reference_bumps,
            bump_commands,
            delta_k2,
            zero_k2,
        ),
    )
    np.save(output / "center_design.npy", design)
    np.save(output / "matched_filter_left_inverses.npy", left_inverses)
    np.save(output / "white_center_covariances.npy", selected["white"])
    np.save(output / "balanced_drift_covariances.npy", selected["balanced_drift"])
    np.save(output / "filtered_drift_covariances.npy", selected["filtered_drift"])
    np.savez_compressed(output / "center_error_samples.npz", **cases)

    analysis_seconds = time.perf_counter() - started
    selected_tradeoff = next(
        row
        for row in tradeoff_rows
        if int(row["repeats_per_signal_state"]) == args.repeats
    )
    metadata_out = {
        "format": "cesr-eight-state-time-series-inverse-v1",
        "date": "2026-08-20",
        "lattice": metadata["lattice"],
        "engine": "SciBmad exact RF-on scan plus SciBmad/GTPSA fixed source templates",
        "target_count": target_count,
        "realizations_per_target": realization_count,
        "signal_state_count": len(states),
        "interleaved_cycle_acquisition_count": len(reference_schedule),
        "repeats_per_signal_state": args.repeats,
        "bpm_noise_rms_per_read_m": args.bpm_noise_rms_m,
        "core_only_endpoint_drift_rms_m": args.core_scan_drift_endpoint_rms_m,
        "reference_scan_endpoint_drift_rms_m": float(
            selected_tradeoff["reference_scan_endpoint_drift_rms_um"]
        )
        * 1e-6,
        "reference_calibration_reads_per_bump": args.reference_calibration_reads,
        "reference_cycle_interval": args.reference_cycle_interval,
        "reference_baseline_count": len(reference_bumps),
        "reference_protocol": "same-bump delta-K2 sequence 0,+,0,-,0",
        "drift_inverse": (
            "random-walk state-space functional filter with four finite-precision "
            "reference-baseline nuisance states"
        ),
        "white_noise_inverse": "full-BPM covariance-matched fixed-template GLS",
        "analysis_wall_seconds": analysis_seconds,
        "monte_carlo_center_draw_count": args.monte_carlo_seeds,
        "required_rmse_um": args.required_rmse_um,
        "preferred_rmse_um": args.preferred_rmse_um,
        "relative_error_reference_scale_um": args.relative_error_scale_um,
        "proxy_relative_rmse_percent": proxy_relative_rmse_percent,
        "balanced_combined_rmse_um": float(balanced_row["rmse_2d_um"]),
        "filtered_combined_rmse_um": float(filtered_row["rmse_2d_um"]),
        "filtered_combined_p99_um": float(filtered_row["p99_2d_um"]),
        "filtered_worst_target_rmse_um": worst_target,
        "hard_gate_passed": bool(hard_gate),
        "preferred_gate_passed": bool(preferred_gate),
        "time_series_semantics": (
            "drift evolves at every acquisition and is never reset at state or cycle boundaries"
        ),
        "reference_baseline_assumption": (
            "each same-bump K2=0 reference mean is independently calibrated with the "
            "declared finite read count; calibration-time drift is not yet modeled"
        ),
    }
    (output / "result_metadata.json").write_text(
        json.dumps(metadata_out, indent=2) + "\n", encoding="utf-8"
    )

    table = "\n".join(
        f"| {row['case']} | {float(row['rmse_2d_um']):.3f} | "
        f"{float(row['median_2d_um']):.3f} | {float(row['p90_2d_um']):.3f} | "
        f"{float(row['p99_2d_um']):.3f} | {float(row['max_2d_um']):.3f} |"
        for row in summary_rows
    )
    report = f"""# Eight-state time-series inverse result

The deterministic signal comes from the latest-lattice SciBmad scan for all
{target_count} active normal sextupoles and {realization_count} hidden
all-sextupole-offset realizations per target.  Every signal-state read has
{args.bpm_noise_rms_m * 1e6:.1f} um RMS BPM white noise.  A scalar physical
orbit drift evolves continuously as a random walk; its endpoint RMS over the
core-only eight-state scan is {args.core_scan_drift_endpoint_rms_m * 1e6:.1f}
um.  Adding references lengthens the equal-cadence scan, so its endpoint RMS
is conservatively increased to
{float(selected_tradeoff['reference_scan_endpoint_drift_rms_um']):.3f} um.

The eight signal states are unchanged.  Every
{args.reference_cycle_interval} cycles (and at the final endpoint), each fixed
bump uses an interleaved `K2=0,+,0,-,0` block.  Other cycles contain only the
eight signal states.  The four K2=0 baseline means are calibrated with
{args.reference_calibration_reads} reads each and retained as
finite-uncertainty nuisance states rather than treated as exact.

| case | 2D RMSE [um] | median [um] | P90 [um] | P99 [um] | maximum [um] |
|---|---:|---:|---:|---:|---:|
{table}

- selected signal reads/state: {args.repeats}
- balanced eight-state acquisitions/target:
  {int(selected_tradeoff['balanced_signal_acquisitions'])}
- interleaved time-series acquisitions/target:
  {int(selected_tradeoff['reference_protocol_acquisitions'])}
- interleaved reference cycles/target:
  {int(selected_tradeoff['interleaved_reference_cycle_count'])}
- separate reference-calibration acquisitions/target:
  {int(selected_tradeoff['reference_calibration_acquisitions'])}
- drift stochastic component:
  {float(selected_tradeoff['balanced_drift_component_rmse_um']):.3f} ->
  {float(selected_tradeoff['filtered_drift_component_rmse_um']):.3f} um
- filtered worst target-level RMSE: {worst_target:.3f} um
- hard gate (aggregate RMSE, P99, and every target RMSE <
  {args.required_rmse_um:.1f} um): {'PASS' if hard_gate else 'FAIL'}
- preferred aggregate RMSE < {args.preferred_rmse_um:.1f} um:
  {'PASS' if preferred_gate else 'FAIL'}
- RMSE relative to the requested {args.relative_error_scale_um:.0f} um scale:
  {proxy_relative_rmse_percent:.3f}%
- stochastic analysis wall time: {analysis_seconds:.3f} s

White-noise suppression is the fixed-template, covariance-matched full-BPM
GLS estimator.  With the present equal independent BPM noise this is
numerically the same matched filter as OLS, but the implementation accepts a
nonuniform diagonal or full measured BPM covariance.  Drift suppression uses
the actual acquisition order and elapsed-step random walk.  Later references
update the already accumulated center-error functional, which is equivalent
to smoothing the time series for the final center without storing the full
BPM tensor.

These are synthetic SciBmad sensitivity results, not demonstrated CESR
precision.  The drift has one calibrated spatial mode per latent machine,
the acquisition cadence is uniform, and the BPM/drift magnitudes are assumed
rather than measured.  Machine deployment still requires measured BPM
covariance and cadence, calibrated drift modes, K2/corrector readbacks,
settling masks, missing/outlier BPM handling, and tests with multidirectional
drift.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    if not hard_gate or not preferred_gate:
        raise RuntimeError("The time-series inverse did not pass all requested gates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
