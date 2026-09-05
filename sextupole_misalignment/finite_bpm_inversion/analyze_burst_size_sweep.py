#!/usr/bin/env python3
"""Paired burst-size sweep for the maintained full-error BPM/GTPSA inverse.

Only acquisition order changes. Every signed state retains 3,072 BPM turns,
the same 13 periodic-reference cycles, and the same finite calibration. Fixed
latent machines and standardized stochastic draws are shared across bursts.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_state_space_bpm_gtpsa_inverse as maintained


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "results" / "burst_size_sweep"
DEFAULT_PRODUCTION = HERE / "results" / "state_space_sequential_bpm_gtpsa_inverse"


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_bursts(text: str) -> list[int]:
    values = sorted({int(value.strip()) for value in text.split(",")})
    if not values or values[0] != 1 or any(value <= 0 for value in values):
        raise ValueError("Burst sizes must be positive and include the B=1 baseline")
    return values


def build_burst_protocol(
    bump_commands: np.ndarray,
    delta_k2: np.ndarray,
    repeats: int,
    reference_cycle_interval: int,
    burst_size: int,
) -> tuple[maintained.AcquisitionProtocol, dict[str, int]]:
    """Group consecutive turns without changing total signal/reference reads."""
    if repeats % burst_size or reference_cycle_interval % burst_size:
        raise ValueError(
            "Burst size must divide both repeats and the reference interval"
        )
    states = maintained.signed_core_states(bump_commands, delta_k2)
    reference_bumps: list[int] = []
    for pair_start in range(0, len(states), 2):
        first, second = states[pair_start : pair_start + 2]
        if first[2] != second[2]:
            raise ValueError("Each signed K2 pair must retain one bump")
        reference_bumps.append(first[2])

    signal_times: list[list[int]] = [[] for _ in states]
    reference_times: list[int] = []
    reference_types: list[int] = []
    physical_settings: list[tuple[int, int]] = []
    acquisition = 0

    def append(core: int = -1, reference: int = -1) -> None:
        nonlocal acquisition
        acquisition += 1
        if core >= 0:
            signal_times[core].append(acquisition)
            physical_settings.append((states[core][2], states[core][3]))
        else:
            reference_times.append(acquisition)
            reference_types.append(reference)
            zero_k2 = int(np.flatnonzero(delta_k2 == 0.0)[0])
            physical_settings.append((reference_bumps[reference], zero_k2))

    visit_count = repeats // burst_size
    reference_visit_interval = reference_cycle_interval // burst_size
    reference_cycles = 0
    for visit in range(visit_count):
        is_reference = (
            visit % reference_visit_interval == 0 or visit == visit_count - 1
        )
        reference_cycles += int(is_reference)
        if is_reference:
            for pair_start, reference in zip(
                range(0, len(states), 2), range(len(reference_bumps))
            ):
                append(reference=reference)
                for _ in range(burst_size):
                    append(core=pair_start)
                append(reference=reference)
                for _ in range(burst_size):
                    append(core=pair_start + 1)
                append(reference=reference)
        else:
            for core in range(len(states)):
                for _ in range(burst_size):
                    append(core=core)

    if any(len(times) != repeats for times in signal_times):
        raise AssertionError("A burst protocol changed the turns per signal state")
    protocol = maintained.AcquisitionProtocol(
        core_states=states,
        reference_bumps=tuple(reference_bumps),
        reference_times=np.asarray(reference_times, dtype=float),
        reference_types=np.asarray(reference_types, dtype=int),
        signal_times=tuple(np.asarray(times, dtype=float) for times in signal_times),
        total_acquisitions=acquisition,
        reference_cycle_count=reference_cycles,
    )
    state_visits = 1 + sum(
        left != right
        for left, right in zip(physical_settings[:-1], physical_settings[1:])
    )
    bump_changes = sum(
        left[0] != right[0]
        for left, right in zip(physical_settings[:-1], physical_settings[1:])
    )
    k2_changes = sum(
        left[1] != right[1]
        for left, right in zip(physical_settings[:-1], physical_settings[1:])
    )
    return protocol, {
        "burst_size": burst_size,
        "visits_per_signal_state": visit_count,
        "idealized_signal_state_visits": len(states) * visit_count,
        "reference_cycle_count": reference_cycles,
        "reference_event_count": len(reference_times),
        "total_signal_turns": repeats * len(states),
        "total_protocol_acquisitions": acquisition,
        "physical_state_visits_excluding_calibration": state_visits,
        "physical_state_switches_excluding_calibration": state_visits - 1,
        "bump_changes_excluding_calibration": bump_changes,
        "k2_changes_excluding_calibration": k2_changes,
    }


def error_metrics(errors_m: np.ndarray) -> dict[str, float]:
    summary = maintained.base.summarize_vectors(errors_m)
    radial_um = np.linalg.norm(errors_m, axis=-1) * 1.0e6
    per_target = np.sqrt(np.mean(radial_um**2, axis=(0, 1)))
    return {
        "rmse_2d_um": float(summary["rmse_2d_um"]),
        "median_2d_um": float(summary["median_2d_um"]),
        "p90_2d_um": float(summary["p90_2d_um"]),
        "p95_2d_um": float(summary["p95_2d_um"]),
        "p99_2d_um": float(summary["p99_2d_um"]),
        "maximum_2d_um": float(summary["max_2d_um"]),
        "worst_target_rmse_2d_um": float(np.max(per_target)),
    }


def prefixed(prefix: str, values: dict[str, float]) -> dict[str, float]:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def plot_summary(rows: list[dict[str, object]], output: Path) -> None:
    burst = np.asarray([int(row["burst_size"]) for row in rows])
    visits = np.asarray(
        [int(row["physical_state_visits_excluding_calibration"]) for row in rows]
    )
    figure, axes = plt.subplots(2, 1, figsize=(8.0, 7.2), sharex=True)
    for prefix, label, style in (
        ("filtered_absolute", "Filtered absolute", "o-"),
        ("unfiltered_absolute", "Unfiltered absolute", "s--"),
    ):
        axes[0].plot(
            burst,
            [float(row[f"{prefix}_rmse_2d_um"]) for row in rows],
            style,
            label=f"{label} RMSE",
        )
        axes[0].plot(
            burst,
            [float(row[f"{prefix}_p99_2d_um"]) for row in rows],
            style,
            alpha=0.55,
            label=f"{label} P99",
        )
    axes[0].axhline(30.0, color="0.35", linestyle=":", label="30 um aggregate gate")
    axes[0].axhline(50.0, color="0.15", linestyle="-.", label="50 um tail gate")
    axes[0].set_ylabel("Absolute center error [um]")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(ncol=2, fontsize=8)
    axes[1].plot(burst, visits, "o-", color="#7b3294")
    axes[1].set_xscale("log", base=2)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Consecutive BPM turns per state visit (burst size)")
    axes[1].set_ylabel("State visits / target")
    axes[1].grid(True, which="both", alpha=0.25)
    figure.suptitle("Full-error paired burst-size sweep")
    figure.tight_layout()
    figure.savefig(output / "burst_size_tradeoff.png", dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-root", type=Path, default=maintained.SCAN_ROOT)
    parser.add_argument("--case", default=maintained.DEFAULT_CASE)
    parser.add_argument("--model-dir", type=Path, default=maintained.DEFAULT_MODEL)
    parser.add_argument("--knobs", type=Path, default=maintained.DEFAULT_KNOBS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--production-results", type=Path, default=DEFAULT_PRODUCTION)
    parser.add_argument("--burst-sizes", default="1,2,4,8,16,32,64,128")
    parser.add_argument("--stochastic-augmentations", type=int, default=32)
    parser.add_argument("--measurement-repeats", type=int, default=3072)
    parser.add_argument("--reference-cycle-interval", type=int, default=256)
    parser.add_argument("--reference-calibration-reads", type=int, default=32)
    parser.add_argument("--bpm-noise-rms-m", type=float, default=5.0e-6)
    parser.add_argument("--drift-endpoint-rms-m", type=float, default=1.0e-5)
    parser.add_argument("--measurement-seed", type=int, default=20261230)
    parser.add_argument("--machine-limit", type=int, default=0)
    parser.add_argument("--target-limit", type=int, default=0)
    args = parser.parse_args()
    started = time.time()
    bursts = parse_bursts(args.burst_sizes)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = args.scan_root.resolve() / args.case
    with (source / "scan_metadata.toml").open("rb") as stream:
        scan_metadata = tomllib.load(stream)
    if (
        scan_metadata.get("baseline_orbit_correction_applied") is not True
        or scan_metadata.get("baseline_response_method") != "reference_gtpsa_orm"
        or scan_metadata.get("baseline_gtpsa_response_model") != "nominal"
        or scan_metadata.get("baseline_gtpsa_validation_enabled") is not False
    ):
        raise ValueError("Source is not the maintained nominal-GTPSA corrected case")

    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in maintained.read_rows(source / "bump_points.csv")
        ]
    )
    delta_k2 = np.asarray(scan_metadata["k2_delta_m3"], dtype=float)
    zero_bump = int(np.flatnonzero(np.all(bump_commands == 0.0, axis=1))[0])
    zero_k2 = int(np.flatnonzero(delta_k2 == 0.0)[0])
    machine_count = int(scan_metadata["machine_count"])
    if args.machine_limit > 0:
        machine_count = min(machine_count, args.machine_limit)
    if args.target_limit > 0:
        target_names = target_names[: args.target_limit]
    target_count = len(target_names)
    augmentations = args.stochastic_augmentations

    model = maintained.subset_orbit_model(
        args.model_dir.resolve(), args.knobs.resolve(), target_names, bpm_names,
        bump_commands,
    )
    all_templates = maintained.derivative.source_templates(
        args.model_dir.resolve(), 0.272
    )
    all_model_names = [
        row["target"]
        for row in maintained.read_rows(args.model_dir.resolve() / "target_locations.csv")
    ]
    template_lookup = {name: index for index, name in enumerate(all_model_names)}
    templates = all_templates[
        np.asarray([template_lookup[name] for name in target_names], dtype=int)
    ]
    observable = np.asarray(
        np.load(source / "observable_bpm_readbacks.npy", mmap_mode="r")
    )[:machine_count, :target_count]
    drift_observable = np.asarray(
        np.load(source / "observable_drift_bpm_readbacks.npy", mmap_mode="r")
    )[:machine_count, :target_count]
    forward_drift = maintained.recover_forward_drift_response(
        observable, drift_observable, float(scan_metadata["drift_halfwidth_m"])
    )
    observable = observable.reshape(
        machine_count, target_count, len(bump_commands), len(delta_k2), -1
    )
    scalar_step_variance = args.drift_endpoint_rms_m**2 / max(
        args.measurement_repeats * 8 - 1, 1
    )

    shape = (len(bursts), augmentations, machine_count, target_count, 2)
    unfiltered_relative = np.zeros(shape)
    unfiltered_absolute = np.zeros(shape)
    filtered_relative = np.zeros(shape)
    filtered_absolute = np.zeros(shape)
    bpm_shape = shape[:-1]
    unfiltered_bpm_error = np.zeros(bpm_shape)
    filtered_bpm_error = np.zeros(bpm_shape)
    protocol_rows: list[dict[str, object]] = []

    baseline_protocol = maintained.build_protocol(
        bump_commands,
        delta_k2,
        args.measurement_repeats,
        args.reference_cycle_interval,
    )
    for burst_index, burst in enumerate(bursts):
        protocol, protocol_stats = build_burst_protocol(
            bump_commands,
            delta_k2,
            args.measurement_repeats,
            args.reference_cycle_interval,
            burst,
        )
        if burst == 1:
            if (
                protocol.core_states != baseline_protocol.core_states
                or protocol.reference_bumps != baseline_protocol.reference_bumps
                or not np.array_equal(
                    protocol.reference_times, baseline_protocol.reference_times
                )
                or not np.array_equal(
                    protocol.reference_types, baseline_protocol.reference_types
                )
                or any(
                    not np.array_equal(left, right)
                    for left, right in zip(
                        protocol.signal_times, baseline_protocol.signal_times
                    )
                )
            ):
                raise AssertionError("B=1 does not exactly reproduce production order")
        protocol_rows.append(protocol_stats)
        brownian_sqrt = maintained.covariance_sqrt(
            maintained.brownian_functional_covariance(
                protocol, scalar_step_variance
            )
        )
        nonzero_bumps = np.asarray(protocol.reference_bumps, dtype=int)
        for target in range(target_count):
            operator = maintained.build_state_space_operator(
                maintained.nominal_drift_matrix(model, target, bump_commands),
                protocol,
                args.bpm_noise_rms_m,
                args.reference_calibration_reads,
                scalar_step_variance,
            )
            simulated = maintained.simulate_forward_target_observables(
                observable[:, target],
                forward_drift[:, target],
                operator,
                protocol,
                zero_k2,
                brownian_sqrt,
                args.bpm_noise_rms_m,
                args.reference_calibration_reads,
                augmentations,
                args.measurement_seed,
                target,
            )
            filtered_drift = maintained.hidden_state_filtered_averages(
                simulated.projected_reference_observations, operator
            )
            unfiltered_core = simulated.unfiltered_core_readbacks
            filtered_core = (
                unfiltered_core
                - filtered_drift @ operator.nominal_bpm_drift.T
            )
            local_all, reference = maintained.reconstruct_target_local_orbits(
                simulated.calibration_readbacks,
                target,
                zero_bump,
                model,
            )
            local = local_all[:, nonzero_bumps]
            unfiltered_slopes, _ = maintained.k2_slopes_from_core(
                unfiltered_core, protocol, delta_k2
            )
            filtered_slopes, _ = maintained.k2_slopes_from_core(
                filtered_core, protocol, delta_k2
            )
            unfiltered_centers = maintained.fixed_template_centers(
                unfiltered_slopes, local, templates[target]
            ).reshape(augmentations, machine_count, 2)
            filtered_centers = maintained.fixed_template_centers(
                filtered_slopes, local, templates[target]
            ).reshape(augmentations, machine_count, 2)
            reference = reference.reshape(augmentations, machine_count, 2)
            unfiltered_relative[burst_index, :, :, target] = unfiltered_centers
            filtered_relative[burst_index, :, :, target] = filtered_centers
            unfiltered_absolute[burst_index, :, :, target] = (
                unfiltered_centers + reference
            )
            filtered_absolute[burst_index, :, :, target] = filtered_centers + reference
            unfiltered_bpm_error[burst_index, :, :, target] = np.sqrt(
                np.mean(
                    (unfiltered_core - simulated.static_core_readbacks) ** 2,
                    axis=(1, 2),
                )
            ).reshape(augmentations, machine_count)
            filtered_bpm_error[burst_index, :, :, target] = np.sqrt(
                np.mean(
                    (filtered_core - simulated.static_core_readbacks) ** 2,
                    axis=(1, 2),
                )
            ).reshape(augmentations, machine_count)
        print(
            f"burst {burst:3d}: {target_count} targets, {machine_count} machines complete",
            flush=True,
        )

    # Persist every machine-facing estimate before opening evaluation truth.
    products = {
        "burst_sizes": np.asarray(bursts, dtype=int),
        "unfiltered_relative_center_estimates": unfiltered_relative,
        "unfiltered_absolute_offset_estimates": unfiltered_absolute,
        "filtered_relative_center_estimates": filtered_relative,
        "filtered_absolute_offset_estimates": filtered_absolute,
        "unfiltered_bpm_state_error_rms_m": unfiltered_bpm_error,
        "filtered_bpm_state_error_rms_m": filtered_bpm_error,
    }
    for name, values in products.items():
        np.save(output / f"{name}.npy", values)
    write_rows(output / "protocol_summary.csv", protocol_rows)

    # Evaluation-only boundary.
    latent_root = args.scan_root.resolve() / "paired_latents"
    exact_reference = np.asarray(
        np.load(source / "reference_target_orbits.npy", mmap_mode="r")
    )[:machine_count, :target_count]
    latent_offsets = np.asarray(
        np.load(latent_root / "sextupole_offsets.npy", mmap_mode="r")
    )[:machine_count, :target_count]
    relative_truth = latent_offsets - exact_reference
    summary_rows: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    for burst_index, burst in enumerate(bursts):
        errors = {
            "unfiltered_relative": (
                unfiltered_relative[burst_index] - relative_truth[None]
            ),
            "unfiltered_absolute": (
                unfiltered_absolute[burst_index] - latent_offsets[None]
            ),
            "filtered_relative": filtered_relative[burst_index] - relative_truth[None],
            "filtered_absolute": filtered_absolute[burst_index] - latent_offsets[None],
        }
        row: dict[str, object] = {**protocol_rows[burst_index]}
        for name, values in errors.items():
            row.update(prefixed(name, error_metrics(values)))
        row["unfiltered_bpm_state_rmse_um"] = float(
            np.sqrt(np.mean(unfiltered_bpm_error[burst_index] ** 2)) * 1.0e6
        )
        row["filtered_bpm_state_rmse_um"] = float(
            np.sqrt(np.mean(filtered_bpm_error[burst_index] ** 2)) * 1.0e6
        )
        summary_rows.append(row)
        for target, name in enumerate(target_names):
            target_row: dict[str, object] = {
                "burst_size": burst,
                "target_index": target + 1,
                "target": name,
            }
            for metric_name, values in errors.items():
                radial_um = np.linalg.norm(values[:, :, target], axis=-1) * 1.0e6
                target_row[f"{metric_name}_rmse_2d_um"] = float(
                    np.sqrt(np.mean(radial_um**2))
                )
                target_row[f"{metric_name}_p99_2d_um"] = float(
                    np.percentile(radial_um, 99)
                )
            per_target_rows.append(target_row)

    baseline = summary_rows[0]
    for row in summary_rows:
        for key in (
            "filtered_relative_rmse_2d_um",
            "filtered_relative_p99_2d_um",
            "filtered_absolute_rmse_2d_um",
            "filtered_absolute_p99_2d_um",
        ):
            row[f"delta_{key}_vs_b1_um"] = float(row[key]) - float(baseline[key])
    write_rows(output / "burst_summary.csv", summary_rows)
    write_rows(output / "per_target_summary.csv", per_target_rows)

    production_differences: dict[str, float] = {}
    production = args.production_results.resolve()
    if machine_count == 16 and target_count == 76 and augmentations == 32:
        comparisons = {
            "unfiltered_fixed_template_relative_center_estimates": (
                unfiltered_relative[0]
            ),
            "unfiltered_fixed_template_absolute_offset_estimates": (
                unfiltered_absolute[0]
            ),
            "filtered_fixed_template_relative_center_estimates": filtered_relative[0],
            "filtered_fixed_template_absolute_offset_estimates": filtered_absolute[0],
            "unfiltered_bpm_state_error_rms_m": unfiltered_bpm_error[0],
            "filtered_bpm_state_error_rms_m": filtered_bpm_error[0],
        }
        for name, current in comparisons.items():
            saved = np.load(production / f"{name}.npy")
            difference = float(np.max(np.abs(current - saved)))
            production_differences[name] = difference
            if difference > 2.0e-15:
                raise AssertionError(
                    f"B=1 changed production {name}: maximum difference {difference}"
                )
    plot_summary(summary_rows, output)

    metadata = {
        "format": "cesr-full-error-burst-size-sweep-v1",
        "date": "2026-09-05",
        "lattice": scan_metadata.get("lattice", ""),
        "source_case": args.case,
        "machine_count": machine_count,
        "target_count": target_count,
        "stochastic_augmentations": augmentations,
        "burst_sizes": bursts,
        "signal_reads_per_state": args.measurement_repeats,
        "reference_cycle_interval_signal_reads": args.reference_cycle_interval,
        "reference_calibration_reads_per_bump": args.reference_calibration_reads,
        "bpm_noise_rms_m_per_read": args.bpm_noise_rms_m,
        "core_only_drift_endpoint_rms_m": args.drift_endpoint_rms_m,
        "measurement_seed": args.measurement_seed,
        "paired_random_semantics": "identical latent machines and seeds; identical standardized Gaussian draws because every burst retains the same stochastic-array shapes",
        "changed_variable": "consecutive signal turns per fixed physical state visit",
        "fixed_quantities": "all maintained latent errors, exact SciBmad observations, nominal GTPSA models, total signal reads, reference count, calibration, noise priors, filter, and fixed-template inverse",
        "drift_time_semantics": "one random-walk step per BPM acquisition; total acquisitions are equal across bursts; hardware settling time is not yet included",
        "truth_boundary": "all burst-indexed machine-facing estimates persisted before exact reference or sextupole offsets were opened",
        "production_b1_max_abs_differences_m": production_differences,
        "analysis_seconds": time.time() - started,
    }
    (output / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    baseline_visits = float(
        summary_rows[0]["physical_state_visits_excluding_calibration"]
    )
    rows_by_burst = {int(row["burst_size"]): row for row in summary_rows}
    b16 = rows_by_burst[16]
    b32 = rows_by_burst[32]
    b64 = rows_by_burst[64]
    lines = [
        "# Full-error burst-size sweep",
        "",
        "This paired latest-lattice SciBmad study changes only the number of",
        "consecutive BPM turns acquired during one fixed magnet-state visit.",
        f"All {machine_count} latent machines, {target_count} sextupoles,",
        f"{augmentations} stochastic realizations, 3,072 turns per signed state,",
        "13 reference cycles, finite calibration, nominal GTPSA transport,",
        "state-space filtering, and the fixed-template inverse are retained.",
        "",
        "| burst | state visits / target | visit reduction | absolute RMSE [um] | RMSE delta vs B=1 [um] | absolute P99 [um] | worst-target RMSE [um] | BPM-state RMSE [um] |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        visit_reduction = baseline_visits / float(
            row["physical_state_visits_excluding_calibration"]
        )
        lines.append(
            "| {burst_size} | {physical_state_visits_excluding_calibration} | "
            f"{visit_reduction:.2f}x | "
            "{filtered_absolute_rmse_2d_um:.3f} | "
            "{delta_filtered_absolute_rmse_2d_um_vs_b1_um:+.3f} | "
            "{filtered_absolute_p99_2d_um:.3f} | "
            "{filtered_absolute_worst_target_rmse_2d_um:.3f} | "
            "{filtered_bpm_state_rmse_um:.3f} |".format(**row)
        )
    lines.extend(
        (
            "",
            "All burst rows contain the same 24,576 signal turns and 156 periodic",
            "reference observations per target. Any wall-time reduction therefore",
            "comes from fewer physical state visits. The current random-walk model",
            "advances per BPM acquisition and contains no corrector or sextupole",
            "settling interval; a facility-time conclusion requires measured settling",
            "and drift spectra.",
            "",
            f"`B=16` reduces visits by {baseline_visits / float(b16['physical_state_visits_excluding_calibration']):.2f}x while changing filtered absolute RMSE by only {float(b16['delta_filtered_absolute_rmse_2d_um_vs_b1_um']):+.3f} micrometers and P99 by {float(b16['delta_filtered_absolute_p99_2d_um_vs_b1_um']):+.3f} micrometers. It is the conservative modeled operating candidate. `B=32` reduces visits by {baseline_visits / float(b32['physical_state_visits_excluding_calibration']):.2f}x, with RMSE and P99 penalties of {float(b32['delta_filtered_absolute_rmse_2d_um_vs_b1_um']):+.3f} and {float(b32['delta_filtered_absolute_p99_2d_um_vs_b1_um']):+.3f} micrometers; it is an aggressive candidate requiring target-wise and machine-data checks. `B=64` reaches {float(b64['filtered_absolute_rmse_2d_um']):.3f} micrometers and crosses the 30-micrometer aggregate gate, so it is not the current default.",
            "",
            "The scalar filtered BPM-state RMSE decreases slightly with burst size,",
            "but final center error rises from `B=32` onward. Consecutive-state",
            "clustering changes how residual drift projects onto the signed parity",
            "contrasts used by the inverse, so BPM-state RMS alone is not a valid",
            "selection metric. The baseline P99 and worst-target RMSE already exceed",
            "the stricter 50-micrometer tail gate; burst acquisition does not resolve",
            "that pre-existing model/template limitation.",
            "",
            "The latest lattice emits the documented straight-multipole-in-curved-",
            "reference warning. This study does not vary girder pitch.",
            "",
        )
    )
    (output / "SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
