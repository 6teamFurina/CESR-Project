#!/usr/bin/env python3
"""Run the corrected 76-target inverse with BPM/GTPSA local orbit estimates.

The machine-facing phase consumes only BPM observable readbacks, known bump
and K2 commands, and the nominal latest-lattice SciBmad/GTPSA response and
transport cache.  It reconstructs the relative local orbit at each target,
reconstructs the absolute zero-bump reference orbit, and fits the sextupole
center.  Exact target-local orbit and latent sextupole offsets are loaded only
after every machine-facing estimate has been persisted and are evaluation-only.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
import tomllib
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares


HERE = Path(__file__).resolve().parent
STUDY_ROOT = HERE.parent
DEFAULT_SCAN_ROOT = (
    STUDY_ROOT / "sequential_joint_inverse" / "results" / "exact_joint_machines"
)
DEFAULT_CASE = "with_quadrupole_misalignment_gtpsa_noisy_corrected"
DEFAULT_MODEL = HERE / "results" / "local_orbit_model"
DEFAULT_KNOBS = (
    STUDY_ROOT
    / "quadrupole_affinity"
    / "exact_11_triplet_validation"
    / "results"
    / "bump_knobs"
    / "local_bump_knobs.csv"
)
DEFAULT_OUTPUT = HERE / "results" / "sequential_bpm_gtpsa_inverse"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fieldnames = list(rows[0])
    for row in rows[1:]:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize_vectors(errors_m: np.ndarray) -> dict[str, float]:
    errors = np.asarray(errors_m, dtype=float)
    radial_um = np.linalg.norm(errors, axis=-1) * 1.0e6
    return {
        "x_rmse_um": float(np.sqrt(np.mean(errors[..., 0] ** 2)) * 1.0e6),
        "y_rmse_um": float(np.sqrt(np.mean(errors[..., 1] ** 2)) * 1.0e6),
        "rmse_2d_um": float(np.sqrt(np.mean(radial_um**2))),
        "median_2d_um": float(np.median(radial_um)),
        "p90_2d_um": float(np.percentile(radial_um, 90)),
        "p95_2d_um": float(np.percentile(radial_um, 95)),
        "p99_2d_um": float(np.percentile(radial_um, 99)),
        "max_2d_um": float(np.max(radial_um)),
    }


def source_matrix(local_xy_m: np.ndarray, center_m: np.ndarray) -> np.ndarray:
    x = local_xy_m[:, 0] - center_m[0]
    y = local_xy_m[:, 1] - center_m[1]
    return np.column_stack((0.5 * (x * x - y * y), x * y))


def fit_profiled_center(
    slopes: np.ndarray,
    local_xy_m: np.ndarray,
    center_bound_m: float,
) -> tuple[np.ndarray, float, bool]:
    """Fit a two-plane center while profiling out two BPM response vectors."""
    channel_scale = np.sqrt(np.mean(slopes * slopes, axis=0))
    positive = channel_scale[np.isfinite(channel_scale) & (channel_scale > 0)]
    floor = np.median(positive) * 1.0e-8 if positive.size else 1.0
    normalized = slopes / np.maximum(channel_scale, floor)

    def residual(center: np.ndarray) -> np.ndarray:
        source = source_matrix(local_xy_m, center)
        propagation = np.linalg.lstsq(source, normalized, rcond=1.0e-12)[0]
        return (source @ propagation - normalized).ravel()

    raw_starts = (
        np.zeros(2),
        np.mean(local_xy_m, axis=0),
        np.array((np.min(local_xy_m[:, 0]), 0.0)),
        np.array((np.max(local_xy_m[:, 0]), 0.0)),
        np.array((0.0, np.min(local_xy_m[:, 1]))),
        np.array((0.0, np.max(local_xy_m[:, 1]))),
    )
    margin = max(1.0e-15, center_bound_m * 1.0e-12)
    starts = [
        np.clip(start, -center_bound_m + margin, center_bound_m - margin)
        for start in raw_starts
    ]
    solutions = [
        least_squares(
            residual,
            start,
            bounds=(-center_bound_m, center_bound_m),
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=500,
        )
        for start in starts
    ]
    selected = min(solutions, key=lambda result: float(np.dot(result.fun, result.fun)))
    relative_residual = float(
        np.linalg.norm(selected.fun) / max(np.linalg.norm(normalized), 1.0e-30)
    )
    at_bound = bool(
        np.any(np.abs(selected.x) >= center_bound_m * (1.0 - 1.0e-7))
    )
    return np.asarray(selected.x), relative_residual, at_bound


def fit_all_centers(
    slopes: np.ndarray,
    local_orbits_m: np.ndarray,
    center_bound_m: float,
    progress_label: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit arrays shaped machine x target x bump x BPM-channel."""
    machine_count, target_count = slopes.shape[:2]
    estimates = np.zeros((machine_count, target_count, 2))
    residuals = np.zeros((machine_count, target_count))
    bound_hits = np.zeros((machine_count, target_count), dtype=bool)
    for machine in range(machine_count):
        for target in range(target_count):
            estimate, residual, at_bound = fit_profiled_center(
                slopes[machine, target],
                local_orbits_m[machine, target],
                center_bound_m,
            )
            estimates[machine, target] = estimate
            residuals[machine, target] = residual
            bound_hits[machine, target] = at_bound
        print(
            f"{progress_label}: machine {machine + 1}/{machine_count}",
            flush=True,
        )
    return estimates, residuals, bound_hits


@dataclass(frozen=True)
class OrbitModel:
    target_names: list[str]
    bpm_names: list[str]
    control_rows: list[dict[str, str]]
    model_bpm_bumps: np.ndarray
    model_target_bumps: np.ndarray
    two_sided_maps: np.ndarray
    neighbor_rows: list[dict[str, object]]
    nominal_bpm_orbits: np.ndarray
    nominal_target_orbits: np.ndarray


def load_orbit_model(
    model_dir: Path,
    knobs_path: Path,
    target_names: list[str],
    bpm_names: list[str],
    bump_commands: np.ndarray,
) -> OrbitModel:
    # This import is deliberately local: the machine-facing module reuses the
    # maintained target-specific two-sided map construction without importing
    # any target-orbit truth.
    sys.path.insert(0, str(HERE))
    from analyze_local_orbit_predictors import build_two_sided_maps

    bpm_rows = read_rows(model_dir / "bpm_locations.csv")
    target_rows = read_rows(model_dir / "target_locations.csv")
    control_rows = read_rows(model_dir / "control_inventory.csv")
    model_bpm_names = [row["bpm"] for row in bpm_rows]
    model_target_names = [row["target"] for row in target_rows]
    if bpm_names != model_bpm_names:
        raise ValueError("Scan and GTPSA model BPM inventories differ")
    if target_names != model_target_names:
        raise ValueError("Scan and GTPSA model target inventories differ")

    bpm_response = np.load(model_dir / "bpm_control_response.npy")
    target_response = np.load(model_dir / "target_control_response.npy").reshape(
        len(target_names), 2, -1
    )
    bpm_maps = np.load(model_dir / "bpm_cumulative_maps.npy")
    target_maps = np.load(model_dir / "target_cumulative_maps.npy")
    one_turn_map = np.load(model_dir / "one_turn_map.npy")
    nominal_bpm = np.load(model_dir / "nominal_bpm_orbits.npy")
    nominal_target = np.load(model_dir / "nominal_target_orbits.npy")
    detector_count = len(bpm_names)
    control_count = len(control_rows)
    if bpm_response.shape != (2 * detector_count, control_count):
        raise ValueError(f"Unexpected BPM response shape: {bpm_response.shape}")
    if target_response.shape != (len(target_names), 2, control_count):
        raise ValueError(f"Unexpected target response shape: {target_response.shape}")
    if nominal_bpm.shape != (detector_count, 2):
        raise ValueError(f"Unexpected nominal BPM orbit shape: {nominal_bpm.shape}")
    if nominal_target.shape != (len(target_names), 2):
        raise ValueError(
            f"Unexpected nominal target orbit shape: {nominal_target.shape}"
        )

    control_lookup = {
        (row["corrector"], row["field"]): index
        for index, row in enumerate(control_rows)
    }
    target_lookup = {name: index for index, name in enumerate(target_names)}
    knob_x = np.zeros((len(target_names), control_count))
    knob_y = np.zeros_like(knob_x)
    for row in read_rows(knobs_path):
        target_name = row["target_sextupole"]
        if target_name not in target_lookup:
            continue
        key = (row["corrector"], row["field"])
        if key not in control_lookup:
            raise ValueError(f"Bump knob uses an unknown control: {key}")
        target = target_lookup[target_name]
        control = control_lookup[key]
        knob_x[target, control] = float(row["field_per_x_bump_m"])
        knob_y[target, control] = float(row["field_per_y_bump_m"])
    command_vectors = (
        knob_x[:, None, :] * bump_commands[None, :, 0, None]
        + knob_y[:, None, :] * bump_commands[None, :, 1, None]
    )
    model_bpm = np.einsum("oc,tbc->tbo", bpm_response, command_vectors)
    model_target = np.einsum("toc,tbc->tbo", target_response, command_vectors)

    two_sided_maps, neighbor_rows = build_two_sided_maps(
        bpm_maps,
        target_maps,
        one_turn_map,
        np.asarray([int(row["line_index"]) for row in bpm_rows]),
        np.asarray([int(row["line_index"]) for row in target_rows]),
    )
    for row, target_name in zip(neighbor_rows, target_names):
        row["target"] = target_name
        row["upstream_bpm"] = bpm_names[int(row["upstream_bpm_index"]) - 1]
        row["downstream_bpm"] = bpm_names[int(row["downstream_bpm_index"]) - 1]
    return OrbitModel(
        target_names=target_names,
        bpm_names=bpm_names,
        control_rows=control_rows,
        model_bpm_bumps=model_bpm,
        model_target_bumps=model_target,
        two_sided_maps=two_sided_maps,
        neighbor_rows=neighbor_rows,
        nominal_bpm_orbits=nominal_bpm,
        nominal_target_orbits=nominal_target,
    )


@dataclass(frozen=True)
class MachineFacingResult:
    predicted_relative_local_orbits: np.ndarray
    predicted_reference_absolute_orbits: np.ndarray
    bpm_k2_slopes: np.ndarray
    relative_center_estimates: np.ndarray
    absolute_offset_estimates: np.ndarray
    fit_relative_residuals: np.ndarray
    fit_bound_hits: np.ndarray


def machine_facing_inverse(
    observable_readbacks: np.ndarray,
    delta_k2: np.ndarray,
    zero_bump: int,
    nominal_k2: int,
    model: OrbitModel,
    center_bound_m: float,
    progress_label: str,
) -> MachineFacingResult:
    """Infer local orbit and center without access to any target-local truth."""
    readbacks = np.asarray(observable_readbacks, dtype=float)
    machine_count, target_count, bump_count, k2_count, detector_count, planes = (
        readbacks.shape
    )
    expected = (
        target_count == len(model.target_names)
        and bump_count == model.model_bpm_bumps.shape[1]
        and k2_count == len(delta_k2)
        and detector_count == len(model.bpm_names)
        and planes == 2
    )
    if not expected:
        raise ValueError(f"Unexpected observable-readback shape: {readbacks.shape}")

    reference_readback = readbacks[:, :, zero_bump, nominal_k2, :, :]
    nominal_k2_readback = readbacks[:, :, :, nominal_k2, :, :]
    observed_relative = nominal_k2_readback - reference_readback[:, :, None, :, :]
    observed_relative_flat = observed_relative.reshape(
        machine_count, target_count, bump_count, 2 * detector_count
    )
    residual = observed_relative_flat - model.model_bpm_bumps[None, :, :, :]
    predicted_relative = np.broadcast_to(
        model.model_target_bumps[None, :, :, :],
        (machine_count, target_count, bump_count, 2),
    ).copy()
    predicted_reference = np.broadcast_to(
        model.nominal_target_orbits[None, :, :],
        (machine_count, target_count, 2),
    ).copy()
    nominal_bpm_flat = model.nominal_bpm_orbits.reshape(-1)
    reference_residual = reference_readback.reshape(
        machine_count, target_count, 2 * detector_count
    ) - nominal_bpm_flat
    for target in range(target_count):
        row = model.neighbor_rows[target]
        upstream = int(row["upstream_bpm_index"]) - 1
        downstream = int(row["downstream_bpm_index"]) - 1
        channels = np.asarray(
            (2 * upstream, 2 * upstream + 1, 2 * downstream, 2 * downstream + 1)
        )
        transport = model.two_sided_maps[target]
        predicted_relative[:, target] += (
            np.take(residual[:, target], channels, axis=-1) @ transport.T
        )
        predicted_reference[:, target] += (
            np.take(reference_residual[:, target], channels, axis=-1) @ transport.T
        )

    centered_k2 = delta_k2 - np.mean(delta_k2)
    slopes = np.einsum(
        "mtbkdp,k->mtbdp", readbacks, centered_k2
    ) / float(np.dot(centered_k2, centered_k2))
    slopes = slopes.reshape(
        machine_count, target_count, bump_count, 2 * detector_count
    )
    estimates, fit_residuals, bound_hits = fit_all_centers(
        slopes,
        predicted_relative,
        center_bound_m,
        progress_label,
    )
    return MachineFacingResult(
        predicted_relative_local_orbits=predicted_relative,
        predicted_reference_absolute_orbits=predicted_reference,
        bpm_k2_slopes=slopes,
        relative_center_estimates=estimates,
        absolute_offset_estimates=estimates + predicted_reference,
        fit_relative_residuals=fit_residuals,
        fit_bound_hits=bound_hits,
    )


def recover_drift_response(
    baseline_readbacks: np.ndarray,
    drift_readbacks: np.ndarray,
    drift_halfwidth_m: float,
) -> np.ndarray:
    """Recover the state-specific BPM derivative for the saved drift secant."""
    if baseline_readbacks.shape != drift_readbacks.shape:
        raise ValueError("Baseline and drift readback tensors differ")
    bump_count, k2_count = baseline_readbacks.shape[2:4]
    fractions = np.linspace(-1.0, 1.0, bump_count * k2_count).reshape(
        bump_count, k2_count
    )
    response = np.zeros_like(baseline_readbacks)
    for bump in range(bump_count):
        for k2 in range(k2_count):
            fraction = float(fractions[bump, k2])
            if fraction != 0.0:
                response[:, :, bump, k2] = (
                    drift_readbacks[:, :, bump, k2]
                    - baseline_readbacks[:, :, bump, k2]
                ) / (drift_halfwidth_m * fraction)
    zero = np.argwhere(fractions == 0.0)
    if zero.shape != (1, 2):
        raise ValueError("Expected one zero-drift scan state")
    bump, k2 = map(int, zero[0])
    response[:, :, bump, k2] = 0.5 * (
        response[:, :, bump, k2 - 1] + response[:, :, bump, k2 + 1]
    )
    return response


def state_mean_random_walk_covariance(
    repeats: int,
    state_count: int,
    endpoint_rms_m: float,
) -> np.ndarray:
    """Covariance of interleaved per-state means for a scalar random walk."""
    read_count = repeats * state_count
    step_variance = endpoint_rms_m**2 / max(read_count - 1, 1)
    covariance = np.zeros((state_count, state_count))
    for p in range(state_count):
        for q in range(state_count):
            total = 0.0
            for repeat in range(repeats):
                n = state_count * repeat + p
                last_before = min(repeats - 1, math.floor((n - q) / state_count))
                count_before = max(last_before + 1, 0)
                arithmetic = (
                    state_count * last_before * (last_before + 1) / 2.0
                    + count_before * (q + 1)
                    if count_before
                    else 0.0
                )
                total += arithmetic + (repeats - count_before) * (n + 1)
            covariance[p, q] = step_variance * total / repeats**2
    covariance = 0.5 * (covariance + covariance.T)
    values, vectors = np.linalg.eigh(covariance)
    return (vectors * np.maximum(values, 0.0)) @ vectors.T


def covariance_square_root(covariance: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    return vectors @ np.diag(np.sqrt(np.maximum(values, 0.0)))


def save_machine_facing_result(output: Path, prefix: str, result: MachineFacingResult) -> None:
    np.save(
        output / f"{prefix}_predicted_relative_local_orbits.npy",
        result.predicted_relative_local_orbits,
    )
    np.save(
        output / f"{prefix}_predicted_reference_absolute_orbits.npy",
        result.predicted_reference_absolute_orbits,
    )
    np.save(output / f"{prefix}_bpm_k2_slopes.npy", result.bpm_k2_slopes)
    np.save(
        output / f"{prefix}_relative_center_estimates.npy",
        result.relative_center_estimates,
    )
    np.save(
        output / f"{prefix}_absolute_offset_estimates.npy",
        result.absolute_offset_estimates,
    )
    np.save(
        output / f"{prefix}_fit_relative_residuals.npy",
        result.fit_relative_residuals,
    )
    np.save(output / f"{prefix}_fit_bound_hits.npy", result.fit_bound_hits)


def per_target_rows(
    acquisition: str,
    method: str,
    target_names: list[str],
    relative_errors_m: np.ndarray,
    absolute_errors_m: np.ndarray,
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for target, name in enumerate(target_names):
        relative = summarize_vectors(relative_errors_m[..., target, :])
        absolute = summarize_vectors(absolute_errors_m[..., target, :])
        result.append(
            {
                "acquisition": acquisition,
                "method": method,
                "target": name,
                "target_index": target + 1,
                **{f"relative_{key}": value for key, value in relative.items()},
                **{f"absolute_{key}": value for key, value in absolute.items()},
            }
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-root", type=Path, default=DEFAULT_SCAN_ROOT)
    parser.add_argument("--case", default=DEFAULT_CASE)
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--knobs", type=Path, default=DEFAULT_KNOBS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--center-bound-m", type=float, default=2.5e-3)
    parser.add_argument("--stochastic-augmentations", type=int, default=32)
    parser.add_argument("--stochastic-machine-count", type=int, default=3)
    parser.add_argument("--measurement-repeats", type=int, default=3072)
    parser.add_argument("--bpm-noise-rms-m", type=float, default=5.0e-6)
    parser.add_argument("--drift-endpoint-rms-m", type=float, default=1.0e-5)
    parser.add_argument("--measurement-seed", type=int, default=20261229)
    args = parser.parse_args()
    started = time.time()
    scan_root = args.scan_root.resolve()
    source = scan_root / args.case
    latent_root = scan_root / "paired_latents"
    model_dir = args.model_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.center_bound_m <= 0:
        raise ValueError("--center-bound-m must be positive")
    if args.stochastic_augmentations < 0:
        raise ValueError("--stochastic-augmentations cannot be negative")
    if args.stochastic_machine_count < 0:
        raise ValueError("--stochastic-machine-count cannot be negative")
    if args.measurement_repeats <= 0:
        raise ValueError("--measurement-repeats must be positive")

    with (source / "scan_metadata.toml").open("rb") as stream:
        scan_metadata = tomllib.load(stream)
    with (model_dir / "model_metadata.toml").open("rb") as stream:
        model_metadata = tomllib.load(stream)
    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    bump_rows = read_rows(source / "bump_points.csv")
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    zero_candidates = np.flatnonzero(np.all(bump_commands == 0.0, axis=1))
    if zero_candidates.size != 1:
        raise ValueError("Expected exactly one zero-bump state")
    zero_bump = int(zero_candidates[0])
    delta_k2 = np.asarray(scan_metadata["k2_delta_m3"], dtype=float)
    nominal_candidates = np.flatnonzero(delta_k2 == 0.0)
    if nominal_candidates.size != 1:
        raise ValueError("Expected exactly one nominal-K2 state")
    nominal_k2 = int(nominal_candidates[0])

    model = load_orbit_model(
        model_dir,
        args.knobs.resolve(),
        target_names,
        bpm_names,
        bump_commands,
    )
    write_rows(output / "two_sided_neighbors.csv", model.neighbor_rows)

    # Machine-facing inputs.  The stored SciBmad orbit is converted to the
    # observable readback by the fixed latent BPM gains.  Target-local arrays
    # and latent sextupole offsets are intentionally not loaded in this phase.
    physical_bpm = np.asarray(np.load(source / "bpm_orbits.npy", mmap_mode="r"))
    drift_physical_bpm = np.asarray(
        np.load(source / "drift_bpm_orbits.npy", mmap_mode="r")
    )
    bpm_gain_errors = np.asarray(np.load(latent_root / "bpm_gain_errors.npy"))
    gain_factors = 1.0 + bpm_gain_errors[:, None, None, None, :, :]
    measured_bpm = physical_bpm * gain_factors
    measured_drift_bpm = drift_physical_bpm * gain_factors
    deterministic = machine_facing_inverse(
        measured_bpm,
        delta_k2,
        zero_bump,
        nominal_k2,
        model,
        args.center_bound_m,
        "deterministic BPM+GTPSA inverse",
    )
    save_machine_facing_result(output, "deterministic", deterministic)

    machine_count, target_count, bump_count, k2_count, detector_count, _ = (
        measured_bpm.shape
    )
    stochastic_count = min(args.stochastic_machine_count, machine_count)
    stochastic_indices = np.arange(machine_count - stochastic_count, machine_count)
    stochastic_payload: dict[str, np.ndarray] = {}
    if args.stochastic_augmentations > 0 and stochastic_count > 0:
        drift_response = recover_drift_response(
            measured_bpm,
            measured_drift_bpm,
            float(scan_metadata["drift_halfwidth_m"]),
        )[stochastic_indices]
        state_count = bump_count * k2_count
        walk_covariance = state_mean_random_walk_covariance(
            args.measurement_repeats,
            state_count,
            args.drift_endpoint_rms_m,
        )
        walk_sqrt = covariance_square_root(walk_covariance)
        mean_noise_std = args.bpm_noise_rms_m / math.sqrt(args.measurement_repeats)
        rng = np.random.default_rng(args.measurement_seed)
        stochastic_relative = np.zeros(
            (args.stochastic_augmentations, stochastic_count, target_count, 2)
        )
        stochastic_absolute = np.zeros_like(stochastic_relative)
        stochastic_local = np.zeros(
            (
                args.stochastic_augmentations,
                stochastic_count,
                target_count,
                bump_count,
                2,
            )
        )
        stochastic_reference = np.zeros_like(stochastic_relative)
        stochastic_fit_residuals = np.zeros(
            (args.stochastic_augmentations, stochastic_count, target_count)
        )
        stochastic_bound_hits = np.zeros_like(
            stochastic_fit_residuals, dtype=bool
        )
        selected_baseline = measured_bpm[stochastic_indices]
        for augmentation in range(args.stochastic_augmentations):
            white = mean_noise_std * rng.standard_normal(selected_baseline.shape)
            state_means = rng.standard_normal(
                (stochastic_count, target_count, state_count)
            ) @ walk_sqrt.T
            state_means = state_means.reshape(
                stochastic_count, target_count, bump_count, k2_count, 1, 1
            )
            readbacks = selected_baseline + white + drift_response * state_means
            result = machine_facing_inverse(
                readbacks,
                delta_k2,
                zero_bump,
                nominal_k2,
                model,
                args.center_bound_m,
                f"stochastic augmentation {augmentation + 1}/{args.stochastic_augmentations}",
            )
            stochastic_relative[augmentation] = result.relative_center_estimates
            stochastic_absolute[augmentation] = result.absolute_offset_estimates
            stochastic_local[augmentation] = result.predicted_relative_local_orbits
            stochastic_reference[augmentation] = (
                result.predicted_reference_absolute_orbits
            )
            stochastic_fit_residuals[augmentation] = result.fit_relative_residuals
            stochastic_bound_hits[augmentation] = result.fit_bound_hits
        stochastic_payload = {
            "relative_center_estimates": stochastic_relative,
            "absolute_offset_estimates": stochastic_absolute,
            "predicted_relative_local_orbits": stochastic_local,
            "predicted_reference_absolute_orbits": stochastic_reference,
            "fit_relative_residuals": stochastic_fit_residuals,
            "fit_bound_hits": stochastic_bound_hits,
        }
        for name, array in stochastic_payload.items():
            np.save(output / f"stochastic_{name}.npy", array)
        np.save(output / "stochastic_machine_indices_zero_based.npy", stochastic_indices)

    # Evaluation-only boundary.  Nothing below this point contributes to the
    # persisted BPM+GTPSA predictions or their fit/model selection.
    exact_target_orbits = np.asarray(
        np.load(source / "target_orbits.npy", mmap_mode="r")
    )
    exact_reference_target = np.asarray(
        np.load(source / "reference_target_orbits.npy", mmap_mode="r")
    )
    latent_offsets = np.asarray(np.load(latent_root / "sextupole_offsets.npy"))
    exact_relative_local = (
        exact_target_orbits[:, :, :, nominal_k2, :]
        - exact_reference_target[:, :, None, :]
    )
    relative_truth = latent_offsets - exact_reference_target
    absolute_truth = latent_offsets
    nonzero_bumps = np.arange(bump_count) != zero_bump

    oracle_relative, oracle_fit_residuals, oracle_bound_hits = fit_all_centers(
        deterministic.bpm_k2_slopes,
        exact_relative_local,
        args.center_bound_m,
        "evaluation-only exact-local oracle",
    )
    oracle_absolute = oracle_relative + exact_reference_target
    np.save(output / "oracle_relative_center_estimates.npy", oracle_relative)
    np.save(output / "oracle_absolute_offset_estimates.npy", oracle_absolute)
    np.save(output / "oracle_fit_relative_residuals.npy", oracle_fit_residuals)
    np.save(output / "oracle_fit_bound_hits.npy", oracle_bound_hits)

    local_rows: list[dict[str, object]] = []
    local_rows.append(
        {
            "acquisition": "deterministic_static_readback",
            "quantity": "relative_local_orbit_nonzero_bumps",
            **summarize_vectors(
                deterministic.predicted_relative_local_orbits[:, :, nonzero_bumps, :]
                - exact_relative_local[:, :, nonzero_bumps, :]
            ),
        }
    )
    local_rows.append(
        {
            "acquisition": "deterministic_static_readback",
            "quantity": "absolute_reference_orbit",
            **summarize_vectors(
                deterministic.predicted_reference_absolute_orbits
                - exact_reference_target
            ),
        }
    )

    center_rows: list[dict[str, object]] = []
    deterministic_relative_error = (
        deterministic.relative_center_estimates - relative_truth
    )
    deterministic_absolute_error = (
        deterministic.absolute_offset_estimates - absolute_truth
    )
    oracle_relative_error = oracle_relative - relative_truth
    oracle_absolute_error = oracle_absolute - absolute_truth
    for method, relative_error, absolute_error, fit_residuals, bound_hits in (
        (
            "bpm_gtpsa_two_sided",
            deterministic_relative_error,
            deterministic_absolute_error,
            deterministic.fit_relative_residuals,
            deterministic.fit_bound_hits,
        ),
        (
            "exact_local_orbit_oracle",
            oracle_relative_error,
            oracle_absolute_error,
            oracle_fit_residuals,
            oracle_bound_hits,
        ),
    ):
        center_rows.append(
            {
                "acquisition": "deterministic_static_readback",
                "method": method,
                "fit_count": int(np.prod(relative_error.shape[:-1])),
                **{
                    f"relative_{key}": value
                    for key, value in summarize_vectors(relative_error).items()
                },
                **{
                    f"absolute_{key}": value
                    for key, value in summarize_vectors(absolute_error).items()
                },
                "fit_residual_median": float(np.median(fit_residuals)),
                "fit_bound_hit_count": int(np.count_nonzero(bound_hits)),
                "worst_target_relative_rmse_2d_um": float(
                    np.max(
                        np.sqrt(
                            np.mean(np.sum(relative_error**2, axis=-1), axis=0)
                        )
                        * 1.0e6
                    )
                ),
                "worst_target_absolute_rmse_2d_um": float(
                    np.max(
                        np.sqrt(
                            np.mean(np.sum(absolute_error**2, axis=-1), axis=0)
                        )
                        * 1.0e6
                    )
                ),
            }
        )

    target_summary_rows = per_target_rows(
        "deterministic_static_readback",
        "bpm_gtpsa_two_sided",
        target_names,
        deterministic_relative_error,
        deterministic_absolute_error,
    ) + per_target_rows(
        "deterministic_static_readback",
        "exact_local_orbit_oracle",
        target_names,
        oracle_relative_error,
        oracle_absolute_error,
    )

    if stochastic_payload:
        selected_exact_relative = exact_relative_local[stochastic_indices]
        selected_reference = exact_reference_target[stochastic_indices]
        selected_relative_truth = relative_truth[stochastic_indices]
        selected_absolute_truth = absolute_truth[stochastic_indices]
        stochastic_local_error = (
            stochastic_payload["predicted_relative_local_orbits"]
            - selected_exact_relative[None]
        )
        stochastic_reference_error = (
            stochastic_payload["predicted_reference_absolute_orbits"]
            - selected_reference[None]
        )
        stochastic_relative_error = (
            stochastic_payload["relative_center_estimates"]
            - selected_relative_truth[None]
        )
        stochastic_absolute_error = (
            stochastic_payload["absolute_offset_estimates"]
            - selected_absolute_truth[None]
        )
        local_rows.extend(
            (
                {
                    "acquisition": "stochastic_15_state_means",
                    "quantity": "relative_local_orbit_nonzero_bumps",
                    **summarize_vectors(
                        stochastic_local_error[:, :, :, nonzero_bumps, :]
                    ),
                },
                {
                    "acquisition": "stochastic_15_state_means",
                    "quantity": "absolute_reference_orbit",
                    **summarize_vectors(stochastic_reference_error),
                },
            )
        )
        center_rows.append(
            {
                "acquisition": "stochastic_15_state_means",
                "method": "bpm_gtpsa_two_sided",
                "fit_count": int(np.prod(stochastic_relative_error.shape[:-1])),
                **{
                    f"relative_{key}": value
                    for key, value in summarize_vectors(
                        stochastic_relative_error
                    ).items()
                },
                **{
                    f"absolute_{key}": value
                    for key, value in summarize_vectors(
                        stochastic_absolute_error
                    ).items()
                },
                "fit_residual_median": float(
                    np.median(stochastic_payload["fit_relative_residuals"])
                ),
                "fit_bound_hit_count": int(
                    np.count_nonzero(stochastic_payload["fit_bound_hits"])
                ),
                "worst_target_relative_rmse_2d_um": float(
                    np.max(
                        np.sqrt(
                            np.mean(
                                np.sum(stochastic_relative_error**2, axis=-1),
                                axis=(0, 1),
                            )
                        )
                        * 1.0e6
                    )
                ),
                "worst_target_absolute_rmse_2d_um": float(
                    np.max(
                        np.sqrt(
                            np.mean(
                                np.sum(stochastic_absolute_error**2, axis=-1),
                                axis=(0, 1),
                            )
                        )
                        * 1.0e6
                    )
                ),
            }
        )
        target_summary_rows.extend(
            per_target_rows(
                "stochastic_15_state_means",
                "bpm_gtpsa_two_sided",
                target_names,
                stochastic_relative_error,
                stochastic_absolute_error,
            )
        )

    write_rows(output / "local_orbit_summary.csv", local_rows)
    write_rows(output / "center_summary.csv", center_rows)
    write_rows(output / "per_target_summary.csv", target_summary_rows)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))
    deterministic_local_error = (
        deterministic.predicted_relative_local_orbits[:, :, nonzero_bumps, :]
        - exact_relative_local[:, :, nonzero_bumps, :]
    )
    per_target_local = (
        np.sqrt(np.mean(np.sum(deterministic_local_error**2, axis=-1), axis=(0, 2)))
        * 1.0e6
    )
    per_target_center = (
        np.sqrt(np.mean(np.sum(deterministic_relative_error**2, axis=-1), axis=0))
        * 1.0e6
    )
    x = np.arange(target_count) + 1
    axes[0].plot(x, per_target_local, marker="o", ms=2.5, lw=1.0)
    axes[0].set_title("BPM+GTPSA local-orbit reconstruction")
    axes[0].set_ylabel("Per-target 2D RMSE [um]")
    axes[0].set_xlabel("Sextupole inventory index")
    axes[0].grid(alpha=0.25)
    axes[1].plot(x, per_target_center, marker="o", ms=2.5, lw=1.0, color="#e15759")
    axes[1].set_title("Beam-relative center after local-orbit reconstruction")
    axes[1].set_xlabel("Sextupole inventory index")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "per_target_bpm_gtpsa_reconstruction.png", dpi=180)
    plt.close(fig)

    deterministic_row = next(
        row
        for row in center_rows
        if row["acquisition"] == "deterministic_static_readback"
        and row["method"] == "bpm_gtpsa_two_sided"
    )
    oracle_row = next(
        row
        for row in center_rows
        if row["acquisition"] == "deterministic_static_readback"
        and row["method"] == "exact_local_orbit_oracle"
    )
    stochastic_row = next(
        (
            row
            for row in center_rows
            if row["acquisition"] == "stochastic_15_state_means"
        ),
        None,
    )
    local_table = "\n".join(
        f"| {row['acquisition']} | {row['quantity']} | "
        f"{float(row['rmse_2d_um']):.3f} | {float(row['p90_2d_um']):.3f} | "
        f"{float(row['p99_2d_um']):.3f} | {float(row['max_2d_um']):.3f} |"
        for row in local_rows
    )
    center_table = "\n".join(
        f"| {row['acquisition']} | {row['method']} | "
        f"{float(row['relative_rmse_2d_um']):.3f} | "
        f"{float(row['relative_p99_2d_um']):.3f} | "
        f"{float(row['absolute_rmse_2d_um']):.3f} | "
        f"{float(row['absolute_p99_2d_um']):.3f} | "
        f"{int(row['fit_bound_hit_count'])} |"
        for row in center_rows
    )
    if stochastic_row is not None:
        matched_static_relative = summarize_vectors(
            deterministic_relative_error[stochastic_indices]
        )
        matched_static_absolute = summarize_vectors(
            deterministic_absolute_error[stochastic_indices]
        )
        stochastic_description = (
            f"The stochastic result uses the last {stochastic_count} machines, "
            f"{args.stochastic_augmentations} independent measurement realizations, "
            f"{args.bpm_noise_rms_m*1e6:.1f} micrometers RMS white noise per BPM "
            f"plane/read averaged over {args.measurement_repeats:,} reads, and a "
            f"{args.drift_endpoint_rms_m*1e6:.1f}-micrometer endpoint-RMS scalar "
            f"random walk over a repeated 15-state acquisition.  On the same "
            f"three-machine subset, the deterministic beam-relative and absolute "
            f"center RMSE values are "
            f"{matched_static_relative['rmse_2d_um']:.3f} and "
            f"{matched_static_absolute['rmse_2d_um']:.3f} micrometers; the full "
            f"16-machine deterministic rows are therefore not the paired "
            f"stochastic baseline."
        )
    else:
        stochastic_description = "The stochastic extension was disabled for this run."
    report = f"""# Sequential BPM/GTPSA local-orbit and sextupole-center inverse

## Result

This run replaces exact target-local orbit coordinates with a machine-facing
two-sided reconstruction.  For every one of {target_count} sextupoles in each
of {machine_count} paired latest-lattice machines, the estimator uses only the
nearest upstream/downstream BPM observable readbacks, known local-bump and K2
commands, and nominal order-one SciBmad/GTPSA cumulative and one-turn maps.
All target-local SciBmad orbit arrays and latent sextupole offsets are loaded
only after the BPM/GTPSA estimates have been saved, and are evaluation-only.

| acquisition | local-orbit quantity | 2D RMSE [um] | P90 [um] | P99 [um] | maximum [um] |
|---|---|---:|---:|---:|---:|
{local_table}

| acquisition | center method | beam-relative RMSE [um] | relative P99 [um] | absolute-offset RMSE [um] | absolute P99 [um] | bound hits |
|---|---|---:|---:|---:|---:|---:|
{center_table}

On deterministic static readbacks, replacing the exact local orbit changes the
beam-relative center RMSE from
`{float(oracle_row['relative_rmse_2d_um']):.3f} um` for the evaluation-only
oracle to `{float(deterministic_row['relative_rmse_2d_um']):.3f} um` for the
BPM/GTPSA reconstruction.  Adding the separately reconstructed absolute
reference orbit gives an absolute sextupole-offset RMSE of
`{float(deterministic_row['absolute_rmse_2d_um']):.3f} um`.

{stochastic_description}

## Method boundary

The local predictor first subtracts the zero-bump, nominal-K2 BPM readback.
Nominal control-to-BPM and control-to-target responses supply the commanded
bump prediction.  The measured-minus-model x/y residual at the nearest BPM on
each side supplies four position coordinates; the GTPSA transport maps infer
the missing upstream transverse momenta and transport the residual to the
target.  The same two-sided operator acts on the absolute BPM-minus-nominal
reference orbit to recover the target reference orbit.  The K2 slopes of all
111 BPMs are then fit against these reconstructed local coordinates while two
unknown propagation vectors are profiled out.

The exact nonlinear RF-on closed orbits remain the SciBmad forward simulation
used to synthesize BPM readbacks.  They are not inputs to the machine-facing
inverse.  The direct high-order K2/offset GTPSA-map limitation is also not
invoked: this route uses the stable order-one transport maps.

## Scope

The static physical ensemble includes all-sextupole offsets, 1% RMS local
corrector gains, 1% RMS K2 gains, 1% RMS BPM gains, independent quadrupole
strength errors within +/-1%, 1 mrad RMS quadrupole rolls, and 50-micrometer
per-plane quadrupole offsets after the recorded noisy-BPM GTPSA-ORM baseline
correction.  BPM offsets, rolls, missing/outlier channels, actuator hysteresis,
unknown model-to-machine calibration mismatch, and sim-to-real validation are
not included.  The 15-state random-walk schedule is a transparent sensitivity
model, not an optimized CESR acquisition order.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")

    metadata = {
        "format": "cesr-sequential-bpm-gtpsa-local-orbit-inverse-v1",
        "date": "2026-08-30",
        "scan_root": str(scan_root),
        "case": args.case,
        "model_dir": str(model_dir),
        "model_engine": model_metadata.get("engine", ""),
        "lattice": scan_metadata.get("lattice", ""),
        "machine_count": machine_count,
        "target_count": target_count,
        "bpm_count": detector_count,
        "bump_count": bump_count,
        "k2_count": k2_count,
        "center_bound_m": args.center_bound_m,
        "stochastic_machine_indices_zero_based": stochastic_indices.tolist(),
        "stochastic_augmentations": args.stochastic_augmentations,
        "measurement_repeats": args.measurement_repeats,
        "bpm_noise_rms_m_per_read": args.bpm_noise_rms_m,
        "bpm_mean_noise_std_m": args.bpm_noise_rms_m
        / math.sqrt(args.measurement_repeats),
        "drift_endpoint_rms_m": args.drift_endpoint_rms_m,
        "measurement_seed": args.measurement_seed,
        "truth_boundary": (
            "target_orbits, reference_target_orbits, and sextupole_offsets are "
            "loaded only after all BPM/GTPSA predictions are persisted"
        ),
        "machine_facing_inputs": (
            "BPM observable readbacks, known bump/K2 commands, nominal "
            "SciBmad/GTPSA control responses and order-one transport"
        ),
        "generated_seconds": time.time() - started,
    }
    (output / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
