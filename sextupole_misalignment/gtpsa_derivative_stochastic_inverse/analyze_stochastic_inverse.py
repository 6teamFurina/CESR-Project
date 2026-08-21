#!/usr/bin/env python3
"""Fast GTPSA-transport dO/dK2 inverse under BPM noise and orbit drift.

The forward model combines the exact local normal-sextupole dO/dK2 source
polynomial with the latest-lattice order-1 SciBmad/GTPSA cumulative maps.  The
inverse uses fixed source templates, parity contrasts, and closed-form noise
propagation rather than fitting a separate scale to every noisy BPM channel.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
STUDY_ROOT = HERE.parent
DEFAULT_SCAN = HERE / "results" / "exact_k5_b3"
DEFAULT_MODEL = STUDY_ROOT / "finite_bpm_inversion" / "results" / "local_orbit_model"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(errors_m: np.ndarray) -> dict[str, float]:
    radial_um = np.linalg.norm(errors_m, axis=-1) * 1e6
    return {
        "x_rmse_um": float(np.sqrt(np.mean(errors_m[..., 0] ** 2)) * 1e6),
        "y_rmse_um": float(np.sqrt(np.mean(errors_m[..., 1] ** 2)) * 1e6),
        "rmse_2d_um": float(np.sqrt(np.mean(radial_um**2))),
        "median_2d_um": float(np.median(radial_um)),
        "p90_2d_um": float(np.percentile(radial_um, 90)),
        "p95_2d_um": float(np.percentile(radial_um, 95)),
        "p99_2d_um": float(np.percentile(radial_um, 99)),
        "max_2d_um": float(np.max(radial_um)),
    }


def source_templates(model_dir: Path, sextupole_length_m: float) -> np.ndarray:
    """Return BPM dO/dK2 source templates with shape target x channel x source.

    Source 0 multiplies 0.5*(x^2-y^2); source 1 multiplies x*y.  A thin normal
    sextupole gives integrated kicks -L*source0 in px and +L*source1 in py.
    The periodic response is reconstructed from the saved SciBmad/GTPSA
    cumulative maps without any new closed-orbit tracking.
    """
    bpm_maps = np.load(model_dir / "bpm_cumulative_maps.npy")
    target_maps = np.load(model_dir / "target_cumulative_maps.npy")
    one_turn = np.load(model_dir / "one_turn_map.npy")
    bpm_rows = read_rows(model_dir / "bpm_locations.csv")
    target_rows = read_rows(model_dir / "target_locations.csv")
    bpm_lines = np.asarray([int(row["line_index"]) for row in bpm_rows])
    target_lines = np.asarray([int(row["line_index"]) for row in target_rows])
    nt, nd = len(target_rows), len(bpm_rows)
    templates = np.zeros((nt, 2 * nd, 2))
    closure = np.linalg.inv(np.eye(6) - one_turn)

    kick_px = np.zeros(6)
    kick_px[1] = 1.0
    kick_py = np.zeros(6)
    kick_py[3] = 1.0
    for target in range(nt):
        inverse_target = np.linalg.inv(target_maps[target])
        for source, kick in enumerate((kick_px, kick_py)):
            one_turn_source = one_turn @ inverse_target @ kick
            start_closed = closure @ one_turn_source
            for bpm in range(nd):
                response = bpm_maps[bpm] @ start_closed
                if bpm_lines[bpm] > target_lines[target]:
                    response = response + bpm_maps[bpm] @ inverse_target @ kick
                scale = -sextupole_length_m if source == 0 else sextupole_length_m
                templates[target, 2 * bpm, source] = scale * response[0]
                templates[target, 2 * bpm + 1, source] = scale * response[2]
    return templates


def parity_gradients(
    bpm_orbits: np.ndarray,
    delta_k2: np.ndarray,
    bump_commands: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Extract K2-odd and then bump-odd x/y gradients from raw BPM states."""
    negative_k = int(np.argmin(delta_k2))
    positive_k = int(np.argmax(delta_k2))
    k_span = float(delta_k2[positive_k] - delta_k2[negative_k])
    slopes = (
        bpm_orbits[..., positive_k, :, :] - bpm_orbits[..., negative_k, :, :]
    ) / k_span
    slopes = slopes.reshape(*slopes.shape[:3], -1)

    def unique_state(x: float, y: float) -> int:
        selected = np.flatnonzero(
            np.isclose(bump_commands[:, 0], x, atol=1e-15)
            & np.isclose(bump_commands[:, 1], y, atol=1e-15)
        )
        if selected.size != 1:
            raise ValueError(f"Expected one bump state at {(x, y)}, found {selected}")
        return int(selected[0])

    amplitude = float(np.max(np.abs(bump_commands)))
    minus_x, plus_x = unique_state(-amplitude, 0.0), unique_state(amplitude, 0.0)
    minus_y, plus_y = unique_state(0.0, -amplitude), unique_state(0.0, amplitude)
    gradient_x = (slopes[:, :, plus_x] - slopes[:, :, minus_x]) / (2.0 * amplitude)
    gradient_y = (slopes[:, :, plus_y] - slopes[:, :, minus_y]) / (2.0 * amplitude)
    return np.stack((gradient_x, gradient_y), axis=2), amplitude


def center_design(templates: np.ndarray) -> np.ndarray:
    """Map [cx,cy] to stacked [dS/dbx,dS/dby] channel gradients."""
    normal = templates[:, :, 0]
    skew = templates[:, :, 1]
    design_x = np.stack((-normal, -skew), axis=-1)
    design_y = np.stack((-skew, normal), axis=-1)
    return np.concatenate((design_x, design_y), axis=1)


def solve_centers(gradients: np.ndarray, design: np.ndarray) -> np.ndarray:
    nt, nr, _, channels = gradients.shape
    right = np.concatenate((gradients[:, :, 0], gradients[:, :, 1]), axis=-1)
    estimates = np.zeros((nt, nr, 2))
    for target in range(nt):
        estimates[target] = np.linalg.lstsq(
            design[target], right[target].T, rcond=1e-12
        )[0].T
    return estimates


def white_noise_crlb_um(
    design: np.ndarray,
    bpm_noise_rms_m: float,
    repeats_per_signed_state: int,
    k2_span_m3: float,
    bump_span_m: float,
) -> float:
    gradient_variance = (
        2.0 * bpm_noise_rms_m
        / (np.sqrt(repeats_per_signed_state) * k2_span_m3 * bump_span_m)
    ) ** 2
    traces = np.asarray(
        [gradient_variance * np.trace(np.linalg.inv(matrix.T @ matrix)) for matrix in design]
    )
    return float(np.sqrt(np.mean(traces)) * 1e6)


def recover_drift_response(
    baseline: np.ndarray,
    drift_scan: np.ndarray,
    drift_halfwidth_m: float,
) -> np.ndarray:
    """Recover state-specific BPM orbit derivative with respect to scalar drift."""
    if baseline.shape != drift_scan.shape:
        raise ValueError("Baseline and drift scan tensors differ")
    nb, nk = baseline.shape[2:4]
    fractions = np.linspace(-1.0, 1.0, nb * nk).reshape(nb, nk)
    response = np.zeros_like(baseline)
    for bump in range(nb):
        for k2 in range(nk):
            fraction = fractions[bump, k2]
            if fraction != 0.0:
                response[:, :, bump, k2] = (
                    drift_scan[:, :, bump, k2] - baseline[:, :, bump, k2]
                ) / (drift_halfwidth_m * fraction)
    zero = np.argwhere(fractions == 0.0)
    if zero.shape != (1, 2):
        raise ValueError("Expected one zero-drift state")
    bump, k2 = map(int, zero[0])
    response[:, :, bump, k2] = 0.5 * (
        response[:, :, bump, k2 - 1] + response[:, :, bump, k2 + 1]
    )
    return response.reshape(*response.shape[:4], -1)


def signed_state_indices(
    bump_commands: np.ndarray,
    delta_k2: np.ndarray,
) -> list[tuple[int, int, int, int]]:
    """Return (gradient block, sign, bump index, K2 index) in a balanced order."""
    amplitude = float(np.max(np.abs(bump_commands)))

    def state(x: float, y: float) -> int:
        selected = np.flatnonzero(
            np.isclose(bump_commands[:, 0], x, atol=1e-15)
            & np.isclose(bump_commands[:, 1], y, atol=1e-15)
        )
        if selected.size != 1:
            raise ValueError(f"Missing signed bump state {(x, y)}")
        return int(selected[0])

    k_minus, k_plus = int(np.argmin(delta_k2)), int(np.argmax(delta_k2))
    x_minus, x_plus = state(-amplitude, 0.0), state(amplitude, 0.0)
    y_minus, y_plus = state(0.0, -amplitude), state(0.0, amplitude)
    # Within each four-state block the signs are +,-,-,+.  Consequently both
    # sum(w) and sum(t*w) vanish for a state-independent additive drift.
    return [
        (0, +1, x_plus, k_plus),
        (0, -1, x_plus, k_minus),
        (0, -1, x_minus, k_plus),
        (0, +1, x_minus, k_minus),
        (1, +1, y_plus, k_plus),
        (1, -1, y_plus, k_minus),
        (1, -1, y_minus, k_plus),
        (1, +1, y_minus, k_minus),
    ]


def random_walk_center_covariances(
    design: np.ndarray,
    drift_response: np.ndarray,
    states: list[tuple[int, int, int, int]],
    repeats: int,
    k2_span_m3: float,
    bump_span_m: float,
    endpoint_change_rms_m: float,
) -> np.ndarray:
    """Exact center covariance for repeated balanced contrasts and a random walk.

    The endpoint RMS is fixed as the repeat count changes.  The closed form
    sums the reverse cumulative contrast weights over repeated eight-state
    cycles, avoiding construction of an O(repeats) acquisition tensor.
    """
    nt, nr, _, _, channels = drift_response.shape
    covariances = np.zeros((nt, nr, 2, 2))
    read_count = repeats * len(states)
    step_variance = endpoint_change_rms_m**2 / max(read_count - 1, 1)
    q_sum = repeats * (repeats - 1) / 2.0
    q2_sum = repeats * (repeats - 1) * (2 * repeats - 1) / 6.0
    normalization = k2_span_m3 * bump_span_m
    for target in range(nt):
        inverse = np.linalg.inv(design[target].T @ design[target]) @ design[target].T
        blocks = (inverse[:, :channels], inverse[:, channels:])
        for realization in range(nr):
            base_vectors = []
            for block, sign, bump, k2 in states:
                response = drift_response[target, realization, bump, k2]
                base_vectors.append(sign * (blocks[block] @ response) / normalization)
            base_vectors = np.asarray(base_vectors)
            total = np.sum(base_vectors, axis=0)
            suffix = np.cumsum(base_vectors[::-1], axis=0)[::-1]
            tail_outer_sum = np.zeros((2, 2))
            for partial in suffix:
                tail_outer_sum += (
                    q2_sum * np.outer(total, total)
                    + q_sum * (np.outer(total, partial) + np.outer(partial, total))
                    + repeats * np.outer(partial, partial)
                ) / repeats**2
            # d[0] is fixed to zero; remove the nonexistent increment before
            # the first acquisition from the reverse-tail sum.
            first_tail = ((repeats - 1) * total + suffix[0]) / repeats
            tail_outer_sum -= np.outer(first_tail, first_tail)
            covariances[target, realization] = step_variance * tail_outer_sum
    return covariances


def white_center_covariances(
    design: np.ndarray,
    bpm_noise_rms_m: float,
    repeats: int,
    k2_span_m3: float,
    bump_span_m: float,
) -> np.ndarray:
    gradient_variance = (
        2.0 * bpm_noise_rms_m / (np.sqrt(repeats) * k2_span_m3 * bump_span_m)
    ) ** 2
    return np.asarray(
        [gradient_variance * np.linalg.inv(matrix.T @ matrix) for matrix in design]
    )


def expected_rmse_um(clean_errors: np.ndarray, covariance: np.ndarray) -> float:
    covariance = np.broadcast_to(covariance, clean_errors.shape[:2] + (2, 2))
    squared = np.sum(clean_errors * clean_errors, axis=-1) + np.trace(
        covariance, axis1=-2, axis2=-1
    )
    return float(np.sqrt(np.mean(squared)) * 1e6)


def covariance_samples(
    covariance: np.ndarray,
    sample_count: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw zero-mean Gaussian vectors for every target/realization covariance."""
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    if np.min(eigenvalues) < -1e-12 * max(float(np.max(eigenvalues)), 1.0):
        raise ValueError("Center covariance is not positive semidefinite")
    square_root = eigenvectors * np.sqrt(np.maximum(eigenvalues, 0.0))[..., None, :]
    standard = rng.standard_normal((sample_count,) + covariance.shape[:-2] + (2,))
    return np.einsum("...ij,s...j->s...i", square_root, standard)


def case_summary(case: str, errors_m: np.ndarray) -> dict[str, object]:
    return {"case": case, "fit_count": int(np.prod(errors_m.shape[:-1])), **summarize(errors_m)}


def main() -> int:
    analysis_start = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument("--physical-root", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--sextupole-length-m", type=float, default=0.272)
    parser.add_argument("--clean-only", action="store_true")
    parser.add_argument("--target-limit", type=int, default=0)
    parser.add_argument("--repeats-per-signed-state", type=int, default=4096)
    parser.add_argument(
        "--candidate-repeats", default="64,256,1024,1280,2048,4096,16384"
    )
    parser.add_argument("--bpm-noise-rms-m", type=float, default=5.0e-6)
    parser.add_argument("--random-walk-endpoint-rms-m", type=float, default=1.0e-5)
    parser.add_argument("--monte-carlo-seeds", type=int, default=512)
    parser.add_argument("--measurement-seed", type=int, default=20260901)
    parser.add_argument("--required-rmse-um", type=float, default=50.0)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results" / "analysis")
    args = parser.parse_args()
    baseline_dir = args.physical_root.resolve() / "baseline"
    model_dir = args.model_dir.resolve()
    with (baseline_dir / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta_k2 = levels * float(metadata["k2_step_m3"])
    bump_rows = read_rows(baseline_dir / "bump_points.csv")
    bump_commands = np.asarray(
        [(float(row["bump_x_command_m"]), float(row["bump_y_command_m"])) for row in bump_rows]
    )
    bpm = np.asarray(np.load(baseline_dir / "bpm_orbits.npy", mmap_mode="r"))
    if args.target_limit > 0:
        bpm = bpm[: args.target_limit]
    gradients, amplitude = parity_gradients(bpm, delta_k2, bump_commands)
    templates = source_templates(model_dir, args.sextupole_length_m)[: bpm.shape[0]]
    design = center_design(templates)
    estimates = solve_centers(gradients, design)

    zero_bump = int(np.flatnonzero(np.all(bump_commands == 0.0, axis=1))[0])
    zero_k2 = int(np.flatnonzero(levels == 0.0)[0])
    target_truth = np.load(baseline_dir / "target_truth.npy")[: bpm.shape[0]]
    target_orbits = np.load(baseline_dir / "target_orbits.npy", mmap_mode="r")
    relative_truth = target_truth - np.asarray(
        target_orbits[: bpm.shape[0], :, zero_bump, zero_k2, :]
    )
    clean_errors = estimates - relative_truth
    stats = summarize(clean_errors)
    residual = np.concatenate((gradients[:, :, 0], gradients[:, :, 1]), axis=-1)
    fitted = np.einsum("tci,tri->trc", design, estimates)
    relative_residual = float(np.linalg.norm(residual - fitted) / np.linalg.norm(residual))
    print(f"targets={len(templates)} realizations={bpm.shape[1]} bump_amplitude={amplitude:.6g} m")
    print(f"clean 2D RMSE={stats['rmse_2d_um']:.6f} um median={stats['median_2d_um']:.6f} um P90={stats['p90_2d_um']:.6f} um")
    print(f"clean maximum={stats['max_2d_um']:.6f} um parity-gradient relative residual={relative_residual:.6e}")
    if args.clean_only:
        return 0
    k2_span = float(delta_k2.max() - delta_k2.min())
    bump_span = 2.0 * amplitude
    drift_dir = args.physical_root.resolve() / "time_drift"
    with (drift_dir / "scan_metadata.toml").open("rb") as stream:
        drift_metadata = tomllib.load(stream)
    drift_scan = np.asarray(np.load(drift_dir / "bpm_orbits.npy", mmap_mode="r"))[
        : bpm.shape[0]
    ]
    drift_response = recover_drift_response(
        bpm,
        drift_scan,
        float(drift_metadata["drift_halfwidth_m"]),
    )
    states = signed_state_indices(bump_commands, delta_k2)
    candidate_repeats = sorted(
        {int(value) for value in args.candidate_repeats.split(",")}
        | {args.repeats_per_signed_state}
    )
    if any(value <= 0 for value in candidate_repeats):
        raise ValueError("Repeat counts must be positive")
    tradeoff_rows: list[dict[str, object]] = []
    selected_white_covariance = None
    selected_drift_covariance = None
    for repeats in candidate_repeats:
        white_covariance = white_center_covariances(
            design, args.bpm_noise_rms_m, repeats, k2_span, bump_span
        )
        drift_covariance = random_walk_center_covariances(
            design,
            drift_response,
            states,
            repeats,
            k2_span,
            bump_span,
            args.random_walk_endpoint_rms_m,
        )
        white_rmse = expected_rmse_um(clean_errors, white_covariance[:, None])
        drift_rmse = expected_rmse_um(clean_errors, drift_covariance)
        combined_rmse = expected_rmse_um(
            clean_errors,
            white_covariance[:, None] + drift_covariance,
        )
        combined_covariance = white_covariance[:, None] + drift_covariance
        expected_squared = np.sum(clean_errors**2, axis=-1) + np.trace(
            combined_covariance, axis1=-2, axis2=-1
        )
        worst_target_rmse = float(
            np.max(np.sqrt(np.mean(expected_squared, axis=1))) * 1e6
        )
        tradeoff_rows.append(
            {
                "repeats_per_signed_state": repeats,
                "acquisitions_per_target_scan": repeats * len(states),
                "analytic_white_noise_rmse_2d_um": white_rmse,
                "analytic_random_walk_drift_rmse_2d_um": drift_rmse,
                "analytic_combined_rmse_2d_um": combined_rmse,
                "analytic_combined_worst_target_rmse_2d_um": worst_target_rmse,
            }
        )
        print(
            f"R={repeats:6d}: white={white_rmse:9.3f} um "
            f"drift={drift_rmse:9.3f} um combined={combined_rmse:9.3f} um "
            f"worst-target={worst_target_rmse:9.3f} um"
        )
        if repeats == args.repeats_per_signed_state:
            selected_white_covariance = white_covariance
            selected_drift_covariance = drift_covariance

    if selected_white_covariance is None or selected_drift_covariance is None:
        raise AssertionError("Selected covariance was not calculated")
    nt, nr = clean_errors.shape[:2]
    white_covariance_full = np.broadcast_to(
        selected_white_covariance[:, None], (nt, nr, 2, 2)
    ).copy()
    white_draws = covariance_samples(
        white_covariance_full,
        args.monte_carlo_seeds,
        np.random.default_rng(args.measurement_seed),
    )
    drift_draws = covariance_samples(
        selected_drift_covariance,
        args.monte_carlo_seeds,
        np.random.default_rng(args.measurement_seed + 1),
    )
    case_errors = {
        "clean": clean_errors,
        "bpm_white_noise": clean_errors[None] + white_draws,
        "random_walk_drift": clean_errors[None] + drift_draws,
        "combined": clean_errors[None] + white_draws + drift_draws,
    }
    expected_by_case = {
        "clean": stats["rmse_2d_um"],
        "bpm_white_noise": expected_rmse_um(clean_errors, white_covariance_full),
        "random_walk_drift": expected_rmse_um(clean_errors, selected_drift_covariance),
        "combined": expected_rmse_um(
            clean_errors, white_covariance_full + selected_drift_covariance
        ),
    }
    summary_rows = []
    for case, errors in case_errors.items():
        row = case_summary(case, errors)
        row["analytic_expected_rmse_2d_um"] = expected_by_case[case]
        summary_rows.append(row)

    target_names = (baseline_dir / "target_names.txt").read_text(
        encoding="utf-8"
    ).splitlines()[:nt]
    per_target_rows: list[dict[str, object]] = []
    for case, errors in case_errors.items():
        sampled = errors if errors.ndim == 4 else errors[None]
        for target, name in enumerate(target_names):
            per_target_rows.append(
                {
                    "case": case,
                    "target": name,
                    "target_index": target + 1,
                    **summarize(sampled[:, target]),
                }
            )
    per_seed_rows: list[dict[str, object]] = []
    for seed_index in range(args.monte_carlo_seeds):
        for case in ("bpm_white_noise", "random_walk_drift", "combined"):
            per_seed_rows.append(
                {
                    "measurement_seed_index": seed_index + 1,
                    **case_summary(case, case_errors[case][seed_index]),
                }
            )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "protocol_tradeoff.csv", tradeoff_rows)
    write_rows(output / "per_target_summary.csv", per_target_rows)
    write_rows(output / "per_seed_summary.csv", per_seed_rows)
    np.save(output / "source_templates.npy", templates)
    np.save(output / "center_design.npy", design)
    np.save(output / "clean_relative_center_estimates.npy", estimates)
    np.save(output / "white_center_covariances.npy", white_covariance_full)
    np.save(output / "drift_center_covariances.npy", selected_drift_covariance)
    np.savez_compressed(output / "center_error_samples.npz", **case_errors)

    combined_row = next(row for row in summary_rows if row["case"] == "combined")
    stochastic_rows = [row for row in summary_rows if row["case"] != "clean"]
    combined_target_rows = [
        row for row in per_target_rows if row["case"] == "combined"
    ]
    worst_target_rmse = max(float(row["rmse_2d_um"]) for row in combined_target_rows)
    below_threshold_fraction = float(
        np.mean(
            np.linalg.norm(case_errors["combined"], axis=-1) * 1e6
            < args.required_rmse_um
        )
    )
    gate_passed = (
        all(float(row["rmse_2d_um"]) < args.required_rmse_um for row in stochastic_rows)
        and float(combined_row["p99_2d_um"]) < args.required_rmse_um
        and worst_target_rmse < args.required_rmse_um
    )
    analysis_seconds = time.perf_counter() - analysis_start
    generation_seconds = float(metadata["calculation_wall_seconds"]) + float(
        drift_metadata["calculation_wall_seconds"]
    )
    result_metadata = {
        "format": "cesr-gtpsa-derivative-stochastic-inverse-v1",
        "date": "2026-08-19",
        "lattice": metadata["lattice"],
        "forward_model": (
            "exact local normal-sextupole dO/dK2 polynomial propagated by "
            "latest-lattice SciBmad/GTPSA cumulative and one-turn maps"
        ),
        "only_machine_error": "fixed x/y misalignment of all 76 active normal sextupoles",
        "target_count": nt,
        "realizations_per_target": nr,
        "exact_forward_state_count": int(metadata["total_state_count"]) + int(
            drift_metadata["total_state_count"]
        ),
        "exact_generation_wall_seconds": generation_seconds,
        "analysis_wall_seconds": analysis_seconds,
        "bump_amplitude_m": amplitude,
        "delta_k2_extrema_m3": [float(delta_k2.min()), float(delta_k2.max())],
        "repeats_per_signed_state": args.repeats_per_signed_state,
        "signed_state_count": len(states),
        "acquisitions_per_target_scan": args.repeats_per_signed_state * len(states),
        "bpm_noise_rms_per_read_m": args.bpm_noise_rms_m,
        "random_walk_endpoint_change_rms_m": args.random_walk_endpoint_rms_m,
        "monte_carlo_measurement_seed_count": args.monte_carlo_seeds,
        "measurement_seed_base": args.measurement_seed,
        "required_rmse_2d_um": args.required_rmse_um,
        "additional_gate": "combined P99 and every target-level combined RMSE must also be below threshold",
        "combined_monte_carlo_rmse_2d_um": float(combined_row["rmse_2d_um"]),
        "combined_monte_carlo_p90_2d_um": float(combined_row["p90_2d_um"]),
        "combined_monte_carlo_p99_2d_um": float(combined_row["p99_2d_um"]),
        "combined_draw_fraction_below_threshold": below_threshold_fraction,
        "combined_worst_target_rmse_2d_um": worst_target_rmse,
        "acceptance_gate_passed": bool(gate_passed),
        "truth_semantics": "beam-relative target magnetic center; truth is loaded only after clean estimates are formed",
        "stochastic_simulation": (
            "independent Gaussian BPM noise and Gaussian random-walk increments are "
            "propagated exactly through the fixed linear parity/matched-filter estimator"
        ),
    }
    (output / "result_metadata.json").write_text(
        json.dumps(result_metadata, indent=2) + "\n", encoding="utf-8"
    )
    table = "\n".join(
        f"| {row['case']} | {float(row['rmse_2d_um']):.3f} | "
        f"{float(row['median_2d_um']):.3f} | {float(row['p90_2d_um']):.3f} | "
        f"{float(row['p99_2d_um']):.3f} | "
        f"{float(row['max_2d_um']):.3f} |"
        for row in summary_rows
    )
    report = f"""# GTPSA-derivative stochastic inverse result

The full latest-lattice benchmark uses all {nt} active normal sextupoles and
{nr} hidden all-sextupole-misalignment realizations per target.  The protocol
uses delta-K2 extrema {delta_k2.min():.3f}/{delta_k2.max():.3f} m^-3, signed
local bumps +/-{amplitude * 1e3:.3f} mm, and
{args.repeats_per_signed_state} reads per signed state.  Every read has
{args.bpm_noise_rms_m * 1e6:.1f} um independent BPM white noise; the physical
orbit drift is a random walk with {args.random_walk_endpoint_rms_m * 1e6:.1f}
um endpoint-change RMS.

| case | 2D RMSE [um] | median [um] | P90 [um] | P99 [um] | maximum [um] |
|---|---:|---:|---:|---:|---:|
{table}

- combined worst target-level RMSE: {worst_target_rmse:.3f} um
- combined draws below {args.required_rmse_um:.1f} um: {100 * below_threshold_fraction:.3f}%
- required threshold: {args.required_rmse_um:.1f} um
- acceptance gate: {'PASS' if gate_passed else 'FAIL'}
- exact SciBmad generation: {generation_seconds:.1f} s for
  {result_metadata['exact_forward_state_count']} paired states
- stochastic inverse and {args.monte_carlo_seeds} measurement seeds:
  {analysis_seconds:.3f} s

The estimator does not fit a separate propagation vector for every noisy BPM
channel.  It fixes the two local sextupole source templates using the validated
SciBmad/GTPSA transport maps, takes K2-odd/bump-odd contrasts in a time-balanced
`+,-,-,+` order, and solves only for x/y center with the known white-noise
covariance.  Even and bump-independent terms cancel from this contrast.  The
random-walk covariance is propagated with an exact reverse-cumulative closed
form, so runtime is independent of the repeat count.

This is a synthetic sensitivity benchmark.  The 5 um BPM noise and 10 um drift
span are assumed settings rather than measured CESR priors.  Gaussian white
noise is unbounded, so the threshold is imposed on aggregate RMSE, P99, and
every target-level RMSE rather than on the single largest Monte Carlo draw.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    if not gate_passed:
        raise RuntimeError("The stochastic inverse did not pass the required error gate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
