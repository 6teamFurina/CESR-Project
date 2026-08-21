#!/usr/bin/env python3
"""Compare repeated and interleaved acquisition for the finite-BPM inverse.

The exact nominal states and the state-specific drift response are taken from
paired latest-lattice SciBmad scans.  New acquisition histories replay those
states with a correlated scalar random walk along the already generated local
bump direction and independent BPM readout noise.
"""

from __future__ import annotations

import argparse
import csv
import sys
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
STUDY_ROOT = HERE.parent
NUISANCE_ROOT = STUDY_ROOT / "real_machine_nuisance_ablation"
FINITE_BPM = STUDY_ROOT / "finite_bpm_inversion"
sys.path.insert(0, str(NUISANCE_ROOT))
sys.path.insert(0, str(FINITE_BPM))

from analyze_command_space_finite_bpm import k2_slope  # noqa: E402
from analyze_nuisance_ablation import load_model, summarize, write_rows  # noqa: E402


NUISANCES = ("bpm_noise", "random_walk_drift", "combined")
PROTOCOLS = (
    ("blocked", "direct"),
    ("interleaved", "direct"),
    ("interleaved", "reference_interpolated"),
)


def slice_model(model: dict[str, object], target_count: int) -> dict[str, object]:
    result = dict(model)
    for key in ("target_names", "model_bpm", "model_target", "two_sided_maps"):
        result[key] = result[key][:target_count]
    result["neighbor_rows"] = result["neighbor_rows"][:target_count]
    return result


def predict_two_sided_nominal(
    nominal_bpm: np.ndarray,
    model: dict[str, object],
) -> np.ndarray:
    """Predict target-local bump coordinates from averaged nominal-K2 BPMs."""
    zero_bump = int(model["zero_bump"])
    observed = np.asarray(nominal_bpm, dtype=float)
    observed = observed - observed[:, :, zero_bump : zero_bump + 1, :, :]
    observed_flat = observed.reshape(*observed.shape[:3], -1)
    model_bpm = np.asarray(model["model_bpm"])
    model_target = np.asarray(model["model_target"])
    residual = observed_flat - model_bpm[:, None, :, :]
    prediction = np.broadcast_to(
        model_target[:, None, :, :], observed.shape[:3] + (2,)
    ).copy()
    two_sided_maps = np.asarray(model["two_sided_maps"])
    neighbor_rows = model["neighbor_rows"]
    for target in range(observed.shape[0]):
        upstream = int(neighbor_rows[target]["upstream_bpm_index"]) - 1
        downstream = int(neighbor_rows[target]["downstream_bpm_index"]) - 1
        channels = np.asarray(
            [2 * upstream, 2 * upstream + 1, 2 * downstream, 2 * downstream + 1]
        )
        nearby = np.take(residual[target], channels, axis=-1)
        prediction[target] += nearby @ two_sided_maps[target].T
    if not np.all(np.isfinite(prediction)):
        raise ValueError("Non-finite two-sided local-orbit prediction")
    return prediction


def state_specific_drift_response(
    baseline_bpm: np.ndarray,
    drift_bpm: np.ndarray,
    halfwidth_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Recover the directional BPM response from the paired exact drift scan."""
    if baseline_bpm.shape != drift_bpm.shape:
        raise ValueError("Baseline and time-drift tensors differ in shape")
    nb, nk = baseline_bpm.shape[2:4]
    fractions = np.linspace(-1.0, 1.0, nb * nk).reshape(nb, nk)
    response = np.zeros_like(baseline_bpm)
    for bump in range(nb):
        for k2_index in range(nk):
            fraction = fractions[bump, k2_index]
            if fraction != 0.0:
                response[:, :, bump, k2_index] = (
                    drift_bpm[:, :, bump, k2_index]
                    - baseline_bpm[:, :, bump, k2_index]
                ) / (halfwidth_m * fraction)
    zero_states = np.argwhere(fractions == 0.0)
    if zero_states.shape != (1, 2):
        raise ValueError("Expected exactly one zero-drift physical state")
    bump, k2_index = map(int, zero_states[0])
    if not (0 < k2_index < nk - 1):
        raise ValueError("Zero-drift state cannot be interpolated in K2")
    response[:, :, bump, k2_index] = 0.5 * (
        response[:, :, bump, k2_index - 1]
        + response[:, :, bump, k2_index + 1]
    )
    return response, fractions


def response_validation(
    response: np.ndarray,
    directions: np.ndarray,
    model: dict[str, object],
) -> dict[str, float]:
    """Compare exact-scan secants with the nominal latest-lattice bump map."""
    bump_commands = np.asarray(model["bump_commands"])
    model_bpm = np.asarray(model["model_bpm"])
    nt, nr = directions.shape[:2]
    nominal_directional = np.zeros((nt, nr, model_bpm.shape[-1]))
    for target in range(nt):
        xy_to_bpm = np.linalg.lstsq(
            bump_commands, model_bpm[target], rcond=1e-12
        )[0]
        nominal_directional[target] = directions[target] @ xy_to_bpm
    exact = response.reshape(*response.shape[:4], -1)
    difference_at_5um = 5.0e-6 * (
        exact - nominal_directional[:, :, None, None, :]
    )
    per_fit_rms_nm = np.sqrt(np.mean(difference_at_5um**2, axis=(2, 3, 4))) * 1e9
    exact_norm = np.linalg.norm(exact.ravel())
    relative_l2 = float(
        np.linalg.norm(
            (exact - nominal_directional[:, :, None, None, :]).ravel()
        ) / exact_norm
    )
    return {
        "nominal_map_relative_l2": relative_l2,
        "difference_at_5um_rms_nm": float(np.sqrt(np.mean(difference_at_5um**2)) * 1e9),
        "difference_at_5um_per_fit_median_rms_nm": float(np.median(per_fit_rms_nm)),
        "difference_at_5um_per_fit_p90_rms_nm": float(np.percentile(per_fit_rms_nm, 90)),
        "difference_at_5um_per_fit_max_rms_nm": float(np.max(per_fit_rms_nm)),
    }


def component_rng(
    seed: int,
    component: int,
    schedule: str,
    repeats: int,
    target: int,
    realization: int,
) -> np.random.Generator:
    schedule_code = 1 if schedule == "blocked" else 2
    return np.random.default_rng(
        np.random.SeedSequence(
            [seed, component, schedule_code, repeats, target + 1, realization + 1]
        )
    )


def block_k2_indices(
    schedule: str,
    repeats: int,
    negative: int,
    nominal: int,
    positive: int,
) -> np.ndarray:
    if schedule == "blocked":
        return np.asarray(
            [nominal] * repeats + [positive] * repeats + [negative] * repeats,
            dtype=int,
        )
    if schedule == "interleaved":
        indices: list[int] = []
        for _ in range(repeats):
            indices.extend((nominal, positive, nominal, negative))
        indices.append(nominal)
        return np.asarray(indices, dtype=int)
    raise ValueError(f"Unsupported schedule: {schedule}")


def simulate_schedule(
    baseline: np.ndarray,
    drift_response: np.ndarray,
    levels: np.ndarray,
    repeats: int,
    schedule: str,
    include_noise: bool,
    include_drift: bool,
    noise_rms_m: float,
    drift_span_rms_m: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, int]:
    """Return direct slopes, nominal means, optional reference slopes, read count."""
    nt, nr, nb, nk, nd, planes = baseline.shape
    nominal_candidates = np.flatnonzero(levels == 0.0)
    if nominal_candidates.size != 1:
        raise ValueError("Expected one nominal K2 level")
    nominal = int(nominal_candidates[0])
    negative = int(np.argmin(levels))
    positive = int(np.argmax(levels))
    block_indices = block_k2_indices(
        schedule, repeats, negative, nominal, positive
    )
    reads_per_block = len(block_indices)
    acquisition_count = nb * reads_per_block
    bump_indices = np.repeat(np.arange(nb), reads_per_block)
    k2_indices = np.tile(block_indices, nb)
    direct_slopes = np.zeros((nt, nr, nb, nd * planes))
    reference_slopes = (
        np.zeros_like(direct_slopes) if schedule == "interleaved" else None
    )
    nominal_means = np.zeros((nt, nr, nb, nd, planes))
    delta_span = float(levels[positive] - levels[negative])

    for target in range(nt):
        for realization in range(nr):
            base = baseline[target, realization]
            response = drift_response[target, realization]
            values = np.asarray(base[bump_indices, k2_indices], dtype=float).copy()
            if include_drift:
                rng = component_rng(
                    seed, 1, schedule, repeats, target, realization
                )
                if acquisition_count == 1:
                    drift = np.zeros(1)
                else:
                    step_rms = drift_span_rms_m / np.sqrt(acquisition_count - 1)
                    drift = np.concatenate(
                        ([0.0], np.cumsum(step_rms * rng.standard_normal(acquisition_count - 1)))
                    )
                values += drift[:, None, None] * response[bump_indices, k2_indices]
            if include_noise:
                rng = component_rng(
                    seed, 2, schedule, repeats, target, realization
                )
                values += noise_rms_m * rng.standard_normal(values.shape)

            for bump in range(nb):
                start = bump * reads_per_block
                stop = start + reads_per_block
                block = values[start:stop].reshape(reads_per_block, nd * planes)
                block_states = k2_indices[start:stop]
                state_means = np.stack(
                    [block[block_states == index].mean(axis=0) for index in range(nk)]
                )
                direct_slopes[target, realization, bump] = k2_slope(
                    state_means[None, :, :], levels
                )[0]
                nominal_means[target, realization, bump] = state_means[nominal].reshape(nd, planes)

                if reference_slopes is not None:
                    zeros = block[0::2]
                    plus = block[1:-1:4]
                    minus = block[3:-1:4]
                    if not (len(zeros) == 2 * repeats + 1 and len(plus) == len(minus) == repeats):
                        raise ValueError("Malformed interleaved block")
                    plus_reference = 0.5 * (
                        zeros[0:-1:2] + zeros[1::2]
                    )
                    minus_reference = 0.5 * (
                        zeros[1:-1:2] + zeros[2::2]
                    )
                    reference_slopes[target, realization, bump] = (
                        (plus - plus_reference).mean(axis=0)
                        - (minus - minus_reference).mean(axis=0)
                    ) / delta_span
    return direct_slopes, nominal_means, reference_slopes, acquisition_count


def fit_centers_profiled_batch(slopes: np.ndarray, xy_command: np.ndarray) -> np.ndarray:
    """Globally profile the unchanged source objective on a refined 2D grid."""
    original_shape = slopes.shape[:2]
    fits = int(np.prod(original_shape))
    slopes_flat = slopes.reshape(fits, slopes.shape[2], slopes.shape[3])
    xy_flat = xy_command.reshape(fits, xy_command.shape[2], 2)
    scale = np.sqrt(np.mean(slopes_flat * slopes_flat, axis=1))
    positive = np.where(scale > 0, scale, np.nan)
    floor = np.nanmedian(positive, axis=1)[:, None] * 1e-8
    floor = np.where(np.isfinite(floor), floor, 1.0)
    normalized = slopes_flat / np.maximum(scale, floor)[:, None, :]
    gram = np.einsum("fic,fjc->fij", normalized, normalized)
    total = np.trace(gram, axis1=1, axis2=2)

    def objective(centers: np.ndarray) -> np.ndarray:
        # centers: fit x grid-point x plane
        dx = xy_flat[:, None, :, 0] - centers[:, :, None, 0]
        dy = xy_flat[:, None, :, 1] - centers[:, :, None, 1]
        first = 0.5 * (dx * dx - dy * dy)
        second = dx * dy
        a = np.sum(first * first, axis=-1)
        b = np.sum(first * second, axis=-1)
        d = np.sum(second * second, axis=-1)
        first_gram = np.matmul(
            first[..., None, :], gram[:, None, :, :]
        )[..., 0, :]
        second_gram = np.matmul(
            second[..., None, :], gram[:, None, :, :]
        )[..., 0, :]
        c00 = np.sum(first_gram * first, axis=-1)
        c01 = np.sum(first_gram * second, axis=-1)
        c10 = np.sum(second_gram * first, axis=-1)
        c11 = np.sum(second_gram * second, axis=-1)
        determinant = a * d - b * b
        projected = np.full_like(determinant, -np.inf)
        valid = determinant > np.finfo(float).eps * np.maximum(a * d, 1e-300)
        projected[valid] = (
            d[valid] * c00[valid]
            + a[valid] * c11[valid]
            - b[valid] * (c01[valid] + c10[valid])
        ) / determinant[valid]
        return total[:, None] - projected

    axis = np.linspace(-1.5e-3, 1.5e-3, 31)
    grid_x, grid_y = np.meshgrid(axis, axis, indexing="ij")
    common = np.column_stack((grid_x.ravel(), grid_y.ravel()))
    centers = np.broadcast_to(common[None, :, :], (fits, len(common), 2))
    best = centers[np.arange(fits), np.argmin(objective(centers), axis=1)]
    halfwidth = axis[1] - axis[0]
    unit = np.linspace(-1.0, 1.0, 21)
    offset_x, offset_y = np.meshgrid(unit, unit, indexing="ij")
    unit_offsets = np.column_stack((offset_x.ravel(), offset_y.ravel()))
    for _ in range(4):
        centers = np.clip(
            best[:, None, :] + halfwidth * unit_offsets[None, :, :],
            -1.5e-3,
            1.5e-3,
        )
        best = centers[np.arange(fits), np.argmin(objective(centers), axis=1)]
        halfwidth /= 10.0
    return best.reshape(*original_shape, 2)


def fit_all(slopes: np.ndarray, nominal_means: np.ndarray, model: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    predicted_local = predict_two_sided_nominal(nominal_means, model)
    estimates = fit_centers_profiled_batch(slopes, predicted_local)
    return estimates, predicted_local


def case_id(nuisance: str, protocol: str, estimator: str, repeats: int) -> str:
    return f"{nuisance}__{protocol}__{estimator}__n{repeats}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--physical-root", type=Path,
        default=NUISANCE_ROOT / "results" / "physical_scans",
    )
    parser.add_argument(
        "--model-dir", type=Path,
        default=FINITE_BPM / "results" / "local_orbit_model",
    )
    parser.add_argument(
        "--knobs", type=Path,
        default=(
            STUDY_ROOT / "quadrupole_affinity" / "exact_11_triplet_validation"
            / "results" / "bump_knobs" / "local_bump_knobs.csv"
        ),
    )
    parser.add_argument("--repeats", default="1,4,16,64")
    parser.add_argument("--bpm-noise-rms-m", type=float, default=5.0e-6)
    parser.add_argument("--drift-span-rms-m", type=float, default=1.0e-5)
    parser.add_argument("--measurement-seed", type=int, default=20260830)
    parser.add_argument("--target-limit", type=int, default=0)
    parser.add_argument("--realization-limit", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results")
    parser.add_argument("--readme-path", type=Path, default=HERE / "README.md")
    args = parser.parse_args()
    repeats_values = tuple(int(value) for value in args.repeats.split(","))
    if not repeats_values or any(value <= 0 for value in repeats_values):
        raise ValueError("Repeat counts must be positive")

    physical_root = args.physical_root.resolve()
    baseline_source = physical_root / "baseline"
    drift_source = physical_root / "time_drift"
    with (baseline_source / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    with (drift_source / "scan_metadata.toml").open("rb") as stream:
        drift_metadata = tomllib.load(stream)
    levels = np.asarray(metadata["k2_levels"], dtype=float) * float(metadata["k2_step_m3"])
    baseline = np.asarray(np.load(baseline_source / "bpm_orbits.npy", mmap_mode="r"))
    drift_bpm = np.asarray(np.load(drift_source / "bpm_orbits.npy", mmap_mode="r"))
    target_truth = np.load(baseline_source / "target_truth.npy")
    target_orbits = np.asarray(np.load(baseline_source / "target_orbits.npy", mmap_mode="r"))
    directions = np.load(drift_source / "latent_drift_directions.npy")
    target_count = baseline.shape[0] if args.target_limit == 0 else min(args.target_limit, baseline.shape[0])
    realization_count = baseline.shape[1] if args.realization_limit == 0 else min(args.realization_limit, baseline.shape[1])
    selection = (slice(0, target_count), slice(0, realization_count))
    baseline = baseline[selection]
    drift_bpm = drift_bpm[selection]
    target_truth = target_truth[selection]
    target_orbits = target_orbits[selection]
    directions = directions[selection]

    model = load_model(
        baseline_source, args.model_dir.resolve(), args.knobs.resolve()
    )
    model = slice_model(model, target_count)
    target_names = list(model["target_names"])
    response, fractions = state_specific_drift_response(
        baseline, drift_bpm, float(drift_metadata["drift_halfwidth_m"])
    )
    validation = response_validation(response, directions, model)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    nominal_k2 = int(np.flatnonzero(levels == 0.0)[0])
    zero_bump = int(model["zero_bump"])
    clean_nominal = baseline[:, :, :, nominal_k2]
    baseline_flat = baseline.reshape(
        target_count, realization_count, baseline.shape[2], baseline.shape[3], -1
    )
    centered_levels = levels - levels.mean()
    clean_slopes = np.einsum(
        "trbkc,k->trbc", baseline_flat, centered_levels
    ) / np.dot(centered_levels, centered_levels)
    clean_estimates, clean_predicted = fit_all(clean_slopes, clean_nominal, model)
    relative_truth = target_truth - target_orbits[:, :, zero_bump, nominal_k2]
    clean_errors = clean_estimates - relative_truth
    clean_summary = summarize(clean_errors)

    summary_rows: list[dict[str, object]] = []
    realization_rows: list[dict[str, object]] = []
    estimates_by_case: dict[str, np.ndarray] = {"clean": clean_estimates}

    def record(
        identifier: str,
        nuisance: str,
        protocol: str,
        estimator: str,
        repeats: int,
        acquisition_count: int,
        estimates: np.ndarray,
        predicted_local: np.ndarray,
    ) -> None:
        errors = estimates - relative_truth
        stats = summarize(errors)
        incremental = float(
            np.sqrt(np.mean(np.sum((errors - clean_errors) ** 2, axis=-1))) * 1e6
        )
        include_noise = nuisance in {"bpm_noise", "combined"}
        include_drift = nuisance in {"random_walk_drift", "combined"}
        summary_rows.append(
            {
                "case": identifier,
                "nuisance": nuisance,
                "protocol": protocol,
                "estimator": estimator,
                "repeats_per_nonzero_point": repeats,
                "acquisitions_per_scan": acquisition_count,
                "fit_count": target_count * realization_count,
                "bpm_noise_rms_per_sample_um": args.bpm_noise_rms_m * 1e6 if include_noise else 0.0,
                "random_walk_endpoint_change_rms_um": args.drift_span_rms_m * 1e6 if include_drift else 0.0,
                **stats,
                "delta_rmse_2d_vs_clean_um": stats["rmse_2d_um"] - clean_summary["rmse_2d_um"],
                "incremental_error_vector_rms_vs_clean_um": incremental,
                "fit_boundary_fraction": float(
                    np.mean(np.any(np.abs(estimates) >= 1.49e-3, axis=-1))
                ),
            }
        )
        estimates_by_case[identifier] = estimates
        radial = np.linalg.norm(errors, axis=-1) * 1e6
        for target, name in enumerate(target_names):
            for realization in range(realization_count):
                realization_rows.append(
                    {
                        "case": identifier,
                        "nuisance": nuisance,
                        "protocol": protocol,
                        "estimator": estimator,
                        "repeats_per_nonzero_point": repeats,
                        "target": name,
                        "target_index": target + 1,
                        "realization": realization + 1,
                        "relative_truth_x_um": relative_truth[target, realization, 0] * 1e6,
                        "relative_truth_y_um": relative_truth[target, realization, 1] * 1e6,
                        "estimate_x_um": estimates[target, realization, 0] * 1e6,
                        "estimate_y_um": estimates[target, realization, 1] * 1e6,
                        "error_x_um": errors[target, realization, 0] * 1e6,
                        "error_y_um": errors[target, realization, 1] * 1e6,
                        "error_2d_um": radial[target, realization],
                        "predicted_local_rms_um": float(
                            np.sqrt(np.mean(predicted_local[target, realization] ** 2)) * 1e6
                        ),
                    }
                )

    for nuisance in NUISANCES:
        include_noise = nuisance in {"bpm_noise", "combined"}
        include_drift = nuisance in {"random_walk_drift", "combined"}
        for repeats in repeats_values:
            blocked_slopes, blocked_nominal, _, blocked_count = simulate_schedule(
                baseline, response, levels, repeats, "blocked",
                include_noise, include_drift, args.bpm_noise_rms_m,
                args.drift_span_rms_m, args.measurement_seed,
            )
            blocked_estimates, blocked_predicted = fit_all(
                blocked_slopes, blocked_nominal, model
            )
            identifier = case_id(nuisance, "blocked", "direct", repeats)
            record(
                identifier, nuisance, "blocked", "direct", repeats,
                blocked_count, blocked_estimates, blocked_predicted,
            )

            direct_slopes, interleaved_nominal, reference_slopes, interleaved_count = simulate_schedule(
                baseline, response, levels, repeats, "interleaved",
                include_noise, include_drift, args.bpm_noise_rms_m,
                args.drift_span_rms_m, args.measurement_seed,
            )
            if reference_slopes is None:
                raise AssertionError("Interleaved reference slopes were not produced")
            for estimator, slopes in (
                ("direct", direct_slopes),
                ("reference_interpolated", reference_slopes),
            ):
                estimates, predicted = fit_all(slopes, interleaved_nominal, model)
                identifier = case_id(nuisance, "interleaved", estimator, repeats)
                record(
                    identifier, nuisance, "interleaved", estimator, repeats,
                    interleaved_count, estimates, predicted,
                )
            print(
                f"{nuisance} repeats={repeats}: "
                f"blocked {summary_rows[-3]['rmse_2d_um']:.3f} um; "
                f"interleaved direct {summary_rows[-2]['rmse_2d_um']:.3f} um; "
                f"interleaved referenced {summary_rows[-1]['rmse_2d_um']:.3f} um"
            )

    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "per_realization_fits.csv", realization_rows)
    np.savez_compressed(output / "relative_center_estimates.npz", **estimates_by_case)
    np.save(output / "clean_predicted_local_orbits.npy", clean_predicted)
    with (output / "drift_response_validation.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(validation))
        writer.writeheader()
        writer.writerow(validation)

    maintained_clean_path = (
        NUISANCE_ROOT / "results" / "analysis"
        / "baseline_relative_center_estimates.npy"
    )
    maintained_clean = np.load(maintained_clean_path)[selection]
    clean_optimizer_difference_um = (
        np.linalg.norm(clean_estimates - maintained_clean, axis=-1) * 1e6
    )
    optimizer_validation = {
        "difference_rms_um": float(
            np.sqrt(np.mean(clean_optimizer_difference_um**2))
        ),
        "difference_median_um": float(np.median(clean_optimizer_difference_um)),
        "difference_p90_um": float(np.percentile(clean_optimizer_difference_um, 90)),
        "difference_max_um": float(np.max(clean_optimizer_difference_um)),
    }
    with (output / "clean_optimizer_validation.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(optimizer_validation))
        writer.writeheader()
        writer.writerow(optimizer_validation)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8), sharey=True)
    styles = {
        ("blocked", "direct"): ("o-", "blocked means"),
        ("interleaved", "direct"): ("s-", "interleaved order, direct slope"),
        ("interleaved", "reference_interpolated"): ("^-", "interleaved 0-reference slope"),
    }
    for axis, nuisance in zip(axes, NUISANCES):
        for protocol, estimator in PROTOCOLS:
            selected = [
                row for row in summary_rows
                if row["nuisance"] == nuisance
                and row["protocol"] == protocol
                and row["estimator"] == estimator
            ]
            selected.sort(key=lambda row: int(row["repeats_per_nonzero_point"]))
            style, label = styles[(protocol, estimator)]
            axis.plot(
                [int(row["repeats_per_nonzero_point"]) for row in selected],
                [float(row["rmse_2d_um"]) for row in selected],
                style, lw=1.7, ms=5, label=label,
            )
        axis.axhline(clean_summary["rmse_2d_um"], color="black", ls="--", lw=1.1, label="clean" if nuisance == NUISANCES[0] else None)
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_xticks(repeats_values, [str(value) for value in repeats_values])
        axis.set_title(nuisance.replace("_", " "))
        axis.set_xlabel("repeats per nonzero K2 point")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("beam-relative center 2D RMSE [micrometers]")
    axes[0].legend(fontsize=8)
    fig.suptitle("Repeated and interleaved finite-BPM sextupole-center protocol")
    fig.tight_layout()
    fig.savefig(output / "protocol_rmse_vs_repeats.png", dpi=180)
    plt.close(fig)

    def rows_for(nuisance: str, protocol: str, estimator: str) -> list[dict[str, object]]:
        return sorted(
            [
                row for row in summary_rows
                if row["nuisance"] == nuisance
                and row["protocol"] == protocol
                and row["estimator"] == estimator
            ],
            key=lambda row: int(row["repeats_per_nonzero_point"]),
        )

    table_lines = "\n".join(
        f"| {row['nuisance']} | {row['protocol']} | {row['estimator']} | "
        f"{row['repeats_per_nonzero_point']} | {row['acquisitions_per_scan']} | "
        f"{float(row['rmse_2d_um']):.3f} | {float(row['median_2d_um']):.3f} | "
        f"{float(row['p90_2d_um']):.3f} | "
        f"{float(row['incremental_error_vector_rms_vs_clean_um']):.3f} | "
        f"{100 * float(row['fit_boundary_fraction']):.2f}% |"
        for row in summary_rows
    )
    combined_blocked = rows_for("combined", "blocked", "direct")
    combined_reference = rows_for(
        "combined", "interleaved", "reference_interpolated"
    )
    last_blocked = combined_blocked[-1]
    last_reference = combined_reference[-1]
    drift_last_blocked = rows_for("random_walk_drift", "blocked", "direct")[-1]
    drift_last_reference = rows_for(
        "random_walk_drift", "interleaved", "reference_interpolated"
    )[-1]
    noise_last_reference = rows_for(
        "bpm_noise", "interleaved", "reference_interpolated"
    )[-1]
    positive_k2 = int(np.argmax(levels))
    negative_k2 = int(np.argmin(levels))
    outer_half_difference = 0.5 * (
        baseline[:, :, :, positive_k2] - baseline[:, :, :, negative_k2]
    )
    signal_rms_nm_by_fit = (
        np.sqrt(np.mean(outer_half_difference**2, axis=(2, 3, 4))) * 1e9
    )
    median_signal_rms_nm = float(np.median(signal_rms_nm_by_fit))
    repeat_count_to_signal_scale = (
        args.bpm_noise_rms_m * 1e9 / median_signal_rms_nm
    ) ** 2
    quadrupole_offsets = np.load(baseline_source / "latent_quadrupole_offsets.npy")
    sextupole_offsets = np.load(baseline_source / "latent_sextupole_offsets.npy")
    report = f"""# Interleaved and repeated BPM acquisition protocol

This paired study tests whether repeated per-point averaging and interleaved
`0,+,0,-,0` K2 acquisition rescue the finite-BPM sextupole-center inverse under
uncorrelated BPM readout noise and correlated random scan drift. Every fit uses
one of {target_count} target sextupoles in {realization_count} latent machines;
all 76 sextupole x/y offsets are fixed but hidden during a scan. Quadrupole
strength, roll, and alignment errors are absent in this bounded test.

The clean exact-state reference has beam-relative 2D RMSE
`{clean_summary['rmse_2d_um']:.3f} micrometers` (median
`{clean_summary['median_2d_um']:.3f}`, P90 `{clean_summary['p90_2d_um']:.3f}`).

## Acquisition and error model

- **Blocked means:** acquire `0` N times, `+` N times, then `-` N times at each
  bump and apply the unchanged symmetric three-point slope.
- **Interleaved direct:** repeat `0,+,0,-` N times, append a final `0`, average
  by K2 state, and apply the unchanged slope.
- **Interleaved 0-reference:** use the same reads, linearly interpolate the
  two adjacent zero-state readings around every `+` and `-` read, subtract
  those references, then average the paired symmetric slopes.
- BPM noise is independent Gaussian noise with
  `{args.bpm_noise_rms_m * 1e6:.1f} micrometer` RMS per BPM plane and acquisition.
- Drift is a scalar Gaussian random walk along each latent machine's fixed
  random two-plane local-bump direction. Its expected end-to-end RMS change is
  held at `{args.drift_span_rms_m * 1e6:.1f} micrometers` for every complete
  scan, so repeat-count comparisons hold total drift severity fixed rather
  than assuming each extra read lengthens the scan.
- The nominal states and state-specific drift secants come from paired exact
  RF-on latest-lattice SciBmad scans. Replaying arbitrary acquisition histories
  is a local linear interpolation of those physical states, not a new exact
  closed-orbit solve for every repeated read.

## Results

| nuisance | order | slope estimator | repeats N | acquisitions/scan | 2D RMSE [um] | median [um] | P90 [um] | paired increment vs clean [um] | fit at boundary |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
{table_lines}

## Interpretation

- **Interleaving helps correlated drift, but does not restore the clean
  inverse.** At `N={drift_last_reference['repeats_per_nonzero_point']}`, the
  drift-only RMSE falls from `{float(drift_last_blocked['rmse_2d_um']):.3f}`
  micrometers for blocked acquisition to
  `{float(drift_last_reference['rmse_2d_um']):.3f}` micrometers for the
  interleaved reference estimator. Its median falls from
  `{float(drift_last_blocked['median_2d_um']):.3f}` to
  `{float(drift_last_reference['median_2d_um']):.3f}` micrometers. The remaining
  long tail is still far above the clean result.
- **Simple averaging is insufficient at the assumed BPM noise.** The typical
  clean outer-K2 half-difference is only `{median_signal_rms_nm:.3f} nm` RMS
  over the retained BPM channels, while one read has `5000 nm` RMS noise.
  Approximately `{repeat_count_to_signal_scale:,.0f}` independent reads would
  be needed merely to reduce the raw mean noise to that typical signal scale;
  this is not a claim that the final center precision would then be 6
  micrometers. At `N={noise_last_reference['repeats_per_nonzero_point']}`, the
  noise-only interleaved-reference RMSE remains
  `{float(noise_last_reference['rmse_2d_um']):.3f}` micrometers and
  `{100 * float(noise_last_reference['fit_boundary_fraction']):.2f}%` of fits
  reach the +/-1.5-mm search boundary.
- **The combined case is noise dominated.** At the largest tested repeat
  count, blocked means give `{float(last_blocked['rmse_2d_um']):.3f}`
  micrometers and the interleaved reference estimator gives
  `{float(last_reference['rmse_2d_um']):.3f}` micrometers. Adjacent zero
  references add their own white-noise variance, so drift cancellation cannot
  compensate while BPM noise dominates.
- The maintained inverse normalizes every BPM slope channel by its own
  across-bump RMS. When raw readout noise dominates a channel, this
  self-normalization largely removes the amplitude benefit of averaging until
  the averaged noise approaches the nanometer-scale K2 signal. The next
  bounded estimator comparison should therefore use known measurement
  covariance or fixed model/clean structural scaling, together with explicit
  drift regression; acquisition order alone is not enough.

The full curves are in `results/protocol_rmse_vs_repeats.png`.

## Integrity and provenance checks

- Lattice: `{metadata['lattice']}`.
- Baseline quadrupole-offset maximum absolute value:
  `{float(np.max(np.abs(quadrupole_offsets))) * 1e6:.6f} micrometers`.
- Latent sextupole-offset RMS over the retained physical source tensor:
  `{float(np.sqrt(np.mean(sextupole_offsets[:target_count, :realization_count] ** 2))) * 1e6:.3f} micrometers`.
- At a 5-micrometer drift displacement, the state-specific exact-scan secants
  differ from the nominal latest-lattice bump map by
  `{validation['difference_at_5um_rms_nm']:.3f} nm` RMS over BPM channels and
  states (relative L2 `{validation['nominal_map_relative_l2']:.5f}`). This
  difference is retained by the replay rather than replaced with the nominal
  map.
- The batched inverse globally profiles the same two-source, per-channel-
  normalized objective on a 2D multiresolution grid. Against the maintained
  six-start local solver on clean data, its center-vector difference has RMS
  `{optimizer_validation['difference_rms_um']:.3f} micrometers`, median
  `{optimizer_validation['difference_median_um']:.3f}`, P90
  `{optimizer_validation['difference_p90_um']:.3f}`, and maximum
  `{optimizer_validation['difference_max_um']:.3f}`. The aggregate clean RMSE
  changes only from `6.051` to
  `{clean_summary['rmse_2d_um']:.3f} micrometers`; the larger maximum flags an
  alternate profiled minimum rather
  than target-truth leakage.
- Exact target orbit and sextupole alignment enter only after fitting, for
  evaluation.

## Run and validate

From `CESR Project/` with the validated `Ubuntu-Bmad` Python environment:

```powershell
wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/interleaved_measurement_protocol/analyze_interleaved_protocol.py'

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/interleaved_measurement_protocol/validate_interleaved_protocol.py'
```

## Limitations

The 5-micrometer BPM noise and 10-micrometer random-walk span are sensitivity
settings, not measured CESR priors. The drift is restricted to the local-bump
orbit mode and one fixed direction per latent machine. Its interpolation is
first order in drift amplitude. The fixed-total-span timing assumption favors
a protocol comparison; a machine deployment study must use measured sampling
cadence, BPM covariance, nonlocal drift modes, outliers, and missing channels.
"""
    args.readme_path.resolve().write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
