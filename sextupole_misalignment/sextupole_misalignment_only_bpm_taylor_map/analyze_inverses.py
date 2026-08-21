#!/usr/bin/env python3
"""Compare BPM-local-orbit, derivative, and high-order Taylor-map inverses."""

from __future__ import annotations

import argparse
import csv
import sys
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares


HERE = Path(__file__).resolve().parent
STUDY_ROOT = HERE.parent
FINITE_BPM = STUDY_ROOT / "finite_bpm_inversion"
sys.path.insert(0, str(FINITE_BPM))

from analyze_command_space_finite_bpm import fit_center  # noqa: E402
from analyze_local_orbit_predictors import build_two_sided_maps  # noqa: E402


DEFAULT_MODEL = FINITE_BPM / "results" / "local_orbit_model"
DEFAULT_KNOBS = (
    STUDY_ROOT
    / "quadrupole_affinity"
    / "exact_11_triplet_validation"
    / "results"
    / "bump_knobs"
    / "local_bump_knobs.csv"
)


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


def summarize_vectors(errors_m: np.ndarray) -> dict[str, float]:
    radial_um = np.linalg.norm(errors_m, axis=-1) * 1e6
    return {
        "x_rmse_um": float(np.sqrt(np.mean(errors_m[..., 0] ** 2)) * 1e6),
        "y_rmse_um": float(np.sqrt(np.mean(errors_m[..., 1] ** 2)) * 1e6),
        "rmse_2d_um": float(np.sqrt(np.mean(radial_um**2))),
        "median_2d_um": float(np.median(radial_um)),
        "p90_2d_um": float(np.percentile(radial_um, 90)),
        "p99_2d_um": float(np.percentile(radial_um, 99)),
        "max_2d_um": float(np.max(radial_um)),
    }


def polynomial_k2_slope(values: np.ndarray, delta_k2: np.ndarray, degree: int) -> np.ndarray:
    """Derivative at zero from a scaled polynomial fit along the K2 axis."""
    scale = float(np.max(np.abs(delta_k2)))
    if not scale > 0:
        raise ValueError("K2 scan has zero span")
    q = delta_k2 / scale
    design = np.column_stack([q**power for power in range(degree + 1)])
    if np.linalg.matrix_rank(design) != degree + 1:
        raise ValueError("K2 polynomial design is rank deficient")
    moved = np.moveaxis(values, -2, 0)
    coefficients = np.linalg.lstsq(design, moved.reshape(len(q), -1), rcond=1e-12)[0]
    derivative = coefficients[1].reshape(moved.shape[1:]) / scale
    return derivative


def linear_k2_slope(values: np.ndarray, delta_k2: np.ndarray) -> np.ndarray:
    centered = delta_k2 - np.mean(delta_k2)
    return np.einsum("...kc,k->...c", values, centered) / np.dot(centered, centered)


def quadratic_design(xy_scaled: np.ndarray) -> np.ndarray:
    x = xy_scaled[:, 0]
    y = xy_scaled[:, 1]
    return np.column_stack(
        (np.ones_like(x), x, y, 0.5 * x * x, x * y, 0.5 * y * y)
    )


def channel_scales(slopes: np.ndarray) -> np.ndarray:
    scale = np.sqrt(np.mean(slopes * slopes, axis=0))
    positive = scale[np.isfinite(scale) & (scale > 0)]
    floor = np.median(positive) * 1e-8 if positive.size else 1.0
    return np.maximum(scale, floor)


def solve_stacked_center(
    gradients: np.ndarray,
    hessians: np.ndarray,
    scales: np.ndarray,
    coordinate_scale_m: float,
) -> tuple[np.ndarray, float, int]:
    rows = (hessians / scales[:, None, None]).reshape(-1, 2)
    right = (-gradients / scales[:, None]).reshape(-1)
    sensitivity = np.linalg.norm(rows, axis=1)
    finite = np.isfinite(right) & np.all(np.isfinite(rows), axis=1)
    if np.any(finite):
        threshold = np.max(sensitivity[finite]) * 1e-9
        finite &= sensitivity > threshold
    rows = rows[finite]
    right = right[finite]
    if rows.shape[0] < 2:
        return np.full(2, np.nan), np.inf, rows.shape[0]
    solution = np.linalg.lstsq(rows, right, rcond=1e-12)[0]
    singular = np.linalg.svd(rows, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else np.inf
    return solution * coordinate_scale_m, condition, rows.shape[0]


def direct_quadratic_derivative_center(
    slopes: np.ndarray,
    local_xy_m: np.ndarray,
    coordinate_scale_m: float,
) -> tuple[np.ndarray, float, int, float]:
    design = quadratic_design(local_xy_m / coordinate_scale_m)
    coefficients = np.linalg.lstsq(design, slopes, rcond=1e-12)[0]
    fitted = design @ coefficients
    gradients = coefficients[1:3].T
    hessians = np.empty((slopes.shape[1], 2, 2))
    hessians[:, 0, 0] = coefficients[3]
    hessians[:, 0, 1] = coefficients[4]
    hessians[:, 1, 0] = coefficients[4]
    hessians[:, 1, 1] = coefficients[5]
    center, condition, rows = solve_stacked_center(
        gradients,
        hessians,
        channel_scales(slopes),
        coordinate_scale_m,
    )
    residual = float(np.linalg.norm(slopes - fitted) / max(np.linalg.norm(slopes), 1e-30))
    return center, condition, rows, residual


def chain_rule_derivative_center(
    slopes: np.ndarray,
    bump_commands_m: np.ndarray,
    predicted_local_xy_m: np.ndarray,
    coordinate_scale_m: float,
) -> tuple[np.ndarray, float, int, float, float, float]:
    bump_scaled = bump_commands_m / coordinate_scale_m
    design = quadratic_design(bump_scaled)
    local_coefficients = np.linalg.lstsq(
        design,
        predicted_local_xy_m / coordinate_scale_m,
        rcond=1e-12,
    )[0]
    local_fitted_m = design @ local_coefficients * coordinate_scale_m
    local_residual_um = float(
        np.sqrt(np.mean(np.sum((local_fitted_m - predicted_local_xy_m) ** 2, axis=1)))
        * 1e6
    )
    jacobian = local_coefficients[1:3].T
    jacobian_condition = float(np.linalg.cond(jacobian))
    jacobian_inverse = np.linalg.pinv(jacobian, rcond=1e-12)

    coefficients = np.linalg.lstsq(design, slopes, rcond=1e-12)[0]
    fitted = design @ coefficients
    gradient_b = coefficients[1:3].T
    hessian_b = np.empty((slopes.shape[1], 2, 2))
    hessian_b[:, 0, 0] = coefficients[3]
    hessian_b[:, 0, 1] = coefficients[4]
    hessian_b[:, 1, 0] = coefficients[4]
    hessian_b[:, 1, 1] = coefficients[5]
    gradient_z = np.einsum("ab,mb->ma", jacobian_inverse.T, gradient_b)
    hessian_z = np.einsum(
        "ab,mbc,cd->mad",
        jacobian_inverse.T,
        hessian_b,
        jacobian_inverse,
    )
    center, condition, rows = solve_stacked_center(
        gradient_z,
        hessian_z,
        channel_scales(slopes),
        coordinate_scale_m,
    )
    residual = float(np.linalg.norm(slopes - fitted) / max(np.linalg.norm(slopes), 1e-30))
    return center, condition, rows, residual, jacobian_condition, local_residual_um


def total_degree_exponents(
    order: int, *, require_k2_dependence: bool = False
) -> list[tuple[int, int, int]]:
    return [
        (x_power, y_power, k_power)
        for degree in range(order + 1)
        for x_power in range(degree + 1)
        for y_power in range(degree - x_power + 1)
        for k_power in [degree - x_power - y_power]
        if not require_k2_dependence or k_power >= 1
    ]


def polynomial_design(
    x_scaled: np.ndarray,
    y_scaled: np.ndarray,
    k_scaled: np.ndarray,
    exponents: list[tuple[int, int, int]],
) -> np.ndarray:
    return np.column_stack(
        [
            x_scaled**x_power * y_scaled**y_power * k_scaled**k_power
            for x_power, y_power, k_power in exponents
        ]
    )


def raw_taylor_map_center(
    bpm_values: np.ndarray,
    local_xy_by_state_m: np.ndarray,
    delta_k2: np.ndarray,
    nominal_k2: int,
    coordinate_scale_m: float,
    order: int,
    initial_center_m: np.ndarray,
) -> tuple[np.ndarray, float, int, float, float]:
    """Fit O(x,y,K2), then find the common root of dO/dK2 at K2=0."""
    nb, nk, channel_count = bpm_values.shape
    if local_xy_by_state_m.shape != (nb, nk, 2):
        raise ValueError("Taylor-map local-coordinate tensor has the wrong shape")
    k_scale = float(np.max(np.abs(delta_k2)))
    x_scaled = local_xy_by_state_m[..., 0].reshape(-1) / coordinate_scale_m
    y_scaled = local_xy_by_state_m[..., 1].reshape(-1) / coordinate_scale_m
    k_scaled = np.broadcast_to(delta_k2[None, :], (nb, nk)).reshape(-1) / k_scale
    # The response is explicitly referenced to the same bump's K2=0 state,
    # so all K2-independent monomials are exactly absent.  Omitting them avoids
    # asking a finite grid to identify irrelevant pure-x/y powers and makes the
    # retained map a K2-dependent Taylor map by construction.
    exponents = total_degree_exponents(order, require_k2_dependence=True)
    # With five K2 levels (including zero), only powers K2^1 through K2^4 are
    # independently identifiable after removal of the K2-independent block.
    exponents = [item for item in exponents if item[2] <= nk - 1]
    design = polynomial_design(x_scaled, y_scaled, k_scaled, exponents)
    rank = int(np.linalg.matrix_rank(design, tol=np.max(design.shape) * np.finfo(float).eps))
    if rank != len(exponents):
        return np.full(2, np.nan), np.inf, rank, np.inf, np.inf

    # Remove each bump's nominal-K2 orbit.  This preserves every nonzero-K2
    # raw state while avoiding a large K2-independent bump orbit in the fit.
    centered = bpm_values - bpm_values[:, nominal_k2 : nominal_k2 + 1, :]
    coefficients = np.linalg.lstsq(design, centered.reshape(-1, channel_count), rcond=1e-12)[0]
    fitted = design @ coefficients
    map_residual = float(
        np.linalg.norm(fitted - centered.reshape(-1, channel_count))
        / max(np.linalg.norm(centered), 1e-30)
    )
    derivative_terms = [
        (index, x_power, y_power)
        for index, (x_power, y_power, k_power) in enumerate(exponents)
        if k_power == 1
    ]
    if not derivative_terms:
        raise ValueError("Taylor map contains no first K2 derivative terms")
    slopes = polynomial_k2_slope(bpm_values, delta_k2, min(4, nk - 1))
    scale = channel_scales(slopes) * k_scale

    def derivative(center_scaled: np.ndarray) -> np.ndarray:
        x, y = center_scaled
        result = np.zeros(channel_count)
        for index, x_power, y_power in derivative_terms:
            result += coefficients[index] * x**x_power * y**y_power
        return result

    def residual(center_scaled: np.ndarray) -> np.ndarray:
        return derivative(center_scaled) / scale

    bound = 1.5e-3 / coordinate_scale_m
    initial_scaled = np.clip(initial_center_m / coordinate_scale_m, -bound, bound)
    starts = [
        initial_scaled,
        np.zeros(2),
        np.mean(local_xy_by_state_m[:, nominal_k2, :], axis=0) / coordinate_scale_m,
        np.array([0.7, 0.0]),
        np.array([-0.7, 0.0]),
        np.array([0.0, 0.7]),
        np.array([0.0, -0.7]),
    ]
    solutions = [
        least_squares(
            residual,
            np.clip(start, -bound, bound),
            bounds=(-bound, bound),
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
            max_nfev=500,
        )
        for start in starts
    ]
    best = min(solutions, key=lambda result: float(np.dot(result.fun, result.fun)))
    singular = np.linalg.svd(design, compute_uv=False)
    condition = float(singular[0] / singular[-1])
    root_rms = float(np.linalg.norm(best.fun) / np.sqrt(channel_count))
    return best.x * coordinate_scale_m, condition, rank, map_residual, root_rms


def load_bump_commands(path: Path) -> np.ndarray:
    rows = read_rows(path)
    rows.sort(key=lambda row: int(row["bump_index"]))
    return np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in rows
        ],
        dtype=float,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=HERE / "results" / "exact_scans")
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--knobs", type=Path, default=DEFAULT_KNOBS)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results" / "analysis")
    args = parser.parse_args()
    source = args.input_dir.resolve()
    model_dir = args.model_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with (source / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    if metadata.get("only_machine_error") != "fixed x/y misalignment of all 76 active normal sextupoles":
        raise ValueError("Input tensor is not the declared sextupole-misalignment-only benchmark")
    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    bump_commands = load_bump_commands(source / "bump_points.csv")
    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta_k2 = levels * float(metadata["k2_step_m3"])
    nominal_k2_candidates = np.flatnonzero(levels == 0.0)
    zero_bump_candidates = np.flatnonzero(np.all(bump_commands == 0.0, axis=1))
    if nominal_k2_candidates.size != 1 or zero_bump_candidates.size != 1:
        raise ValueError("Expected unique zero K2 and zero bump states")
    nominal_k2 = int(nominal_k2_candidates[0])
    zero_bump = int(zero_bump_candidates[0])
    coordinate_scale = float(metadata["bump_amplitude_m"])

    bpm_rows = read_rows(model_dir / "bpm_locations.csv")
    model_target_rows = read_rows(model_dir / "target_locations.csv")
    control_rows = read_rows(model_dir / "control_inventory.csv")
    if bpm_names != [row["bpm"] for row in bpm_rows]:
        raise ValueError("Scan and local-orbit model BPM inventories differ")
    model_target_names = [row["target"] for row in model_target_rows]
    target_model_indices = np.asarray([model_target_names.index(name) for name in target_names])

    bpm_response = np.load(model_dir / "bpm_control_response.npy")
    all_target_response = np.load(model_dir / "target_control_response.npy").reshape(
        len(model_target_names), 2, -1
    )
    target_response = all_target_response[target_model_indices]
    bpm_maps = np.load(model_dir / "bpm_cumulative_maps.npy")
    all_target_maps = np.load(model_dir / "target_cumulative_maps.npy")
    target_maps = all_target_maps[target_model_indices]
    one_turn = np.load(model_dir / "one_turn_map.npy")
    nt = len(target_names)
    nd = len(bpm_names)
    nc = len(control_rows)

    control_lookup = {
        (row["corrector"], row["field"]): index
        for index, row in enumerate(control_rows)
    }
    target_lookup = {name: index for index, name in enumerate(target_names)}
    knob_x = np.zeros((nt, nc))
    knob_y = np.zeros((nt, nc))
    for row in read_rows(args.knobs.resolve()):
        name = row["target_sextupole"]
        if name not in target_lookup:
            continue
        target = target_lookup[name]
        control = control_lookup[(row["corrector"], row["field"])]
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
        one_turn,
        np.asarray([int(row["line_index"]) for row in bpm_rows]),
        np.asarray([int(model_target_rows[index]["line_index"]) for index in target_model_indices]),
    )
    for row, target_name in zip(neighbor_rows, target_names):
        row["target"] = target_name
        row["upstream_bpm"] = bpm_names[int(row["upstream_bpm_index"]) - 1]
        row["downstream_bpm"] = bpm_names[int(row["downstream_bpm_index"]) - 1]
    write_rows(output / "two_sided_neighbors.csv", neighbor_rows)

    # Machine-facing local-orbit prediction.  No target-local truth or target
    # alignment is loaded before these arrays are constructed and persisted.
    bpm_orbits = np.load(source / "bpm_orbits.npy", mmap_mode="r")
    nt_shape, nr, nb, nk, nd_shape, planes = bpm_orbits.shape
    if (nt_shape, nb, nk, nd_shape, planes) != (nt, len(bump_commands), len(levels), nd, 2):
        raise ValueError(f"Unexpected BPM tensor shape: {bpm_orbits.shape}")
    reference_bpm = np.asarray(
        bpm_orbits[:, :, zero_bump, nominal_k2, :, :], dtype=float
    )
    observed_relative = np.asarray(bpm_orbits, dtype=float) - reference_bpm[:, :, None, None, :, :]
    observed_relative_flat = observed_relative.reshape(nt, nr, nb, nk, 2 * nd)
    model_bpm_states = model_bpm[:, None, :, None, :]
    residual = observed_relative_flat - model_bpm_states
    predicted_relative = np.broadcast_to(
        model_target[:, None, :, None, :],
        (nt, nr, nb, nk, 2),
    ).copy()
    for target in range(nt):
        upstream = int(neighbor_rows[target]["upstream_bpm_index"]) - 1
        downstream = int(neighbor_rows[target]["downstream_bpm_index"]) - 1
        channels = np.array(
            [2 * upstream, 2 * upstream + 1, 2 * downstream, 2 * downstream + 1]
        )
        predicted_relative[target] += (
            np.take(residual[target], channels, axis=-1) @ two_sided_maps[target].T
        )
    np.save(output / "predicted_relative_local_orbits.npy", predicted_relative)

    nominal_bpm = np.load(source / "nominal_bpm_orbits.npy").reshape(2 * nd)
    nominal_target = np.load(source / "nominal_target_orbits.npy")
    reference_residual = reference_bpm.reshape(nt, nr, 2 * nd) - nominal_bpm
    predicted_reference_absolute = np.broadcast_to(
        nominal_target[:, None, :], (nt, nr, 2)
    ).copy()
    for target in range(nt):
        upstream = int(neighbor_rows[target]["upstream_bpm_index"]) - 1
        downstream = int(neighbor_rows[target]["downstream_bpm_index"]) - 1
        channels = np.array(
            [2 * upstream, 2 * upstream + 1, 2 * downstream, 2 * downstream + 1]
        )
        predicted_reference_absolute[target] += (
            np.take(reference_residual[target], channels, axis=-1)
            @ two_sided_maps[target].T
        )
    np.save(output / "predicted_reference_absolute_orbits.npy", predicted_reference_absolute)

    estimates: dict[str, np.ndarray] = {
        "fd_linear_source_predicted": np.zeros((nt, nr, 2)),
        "fd_quartic_source_predicted": np.zeros((nt, nr, 2)),
        "quadratic_o_derivative_predicted": np.zeros((nt, nr, 2)),
        "chain_rule_o_derivative_predicted": np.zeros((nt, nr, 2)),
        "o_taylor_order3_nominal_local": np.zeros((nt, nr, 2)),
        "o_taylor_order4_nominal_local": np.zeros((nt, nr, 2)),
        "o_taylor_order5_nominal_local": np.zeros((nt, nr, 2)),
        "o_taylor_order4_all_state_local": np.zeros((nt, nr, 2)),
    }
    diagnostics: list[dict[str, object]] = []
    bpm_flat = np.asarray(bpm_orbits, dtype=float).reshape(nt, nr, nb, nk, 2 * nd)
    for target in range(nt):
        for realization in range(nr):
            values = bpm_flat[target, realization]
            slopes_linear = linear_k2_slope(values, delta_k2)
            slopes_quartic = polynomial_k2_slope(values, delta_k2, 4)
            local_nominal = predicted_relative[target, realization, :, nominal_k2, :]
            estimates["fd_linear_source_predicted"][target, realization] = fit_center(
                slopes_linear, local_nominal
            )
            estimates["fd_quartic_source_predicted"][target, realization] = fit_center(
                slopes_quartic, local_nominal
            )

            direct = direct_quadratic_derivative_center(
                slopes_quartic, local_nominal, coordinate_scale
            )
            estimates["quadratic_o_derivative_predicted"][target, realization] = direct[0]
            diagnostics.append(
                {
                    "target": target_names[target],
                    "target_index": target + 1,
                    "realization": realization + 1,
                    "method": "quadratic_o_derivative_predicted",
                    "inverse_condition": direct[1],
                    "retained_rows": direct[2],
                    "surface_relative_residual": direct[3],
                    "local_map_jacobian_condition": "",
                    "local_map_fit_rmse_um": "",
                    "taylor_design_rank": "",
                    "taylor_design_condition": "",
                    "taylor_map_relative_residual": "",
                    "taylor_root_rms_scaled": "",
                }
            )
            chain = chain_rule_derivative_center(
                slopes_quartic, bump_commands, local_nominal, coordinate_scale
            )
            estimates["chain_rule_o_derivative_predicted"][target, realization] = chain[0]
            diagnostics.append(
                {
                    "target": target_names[target],
                    "target_index": target + 1,
                    "realization": realization + 1,
                    "method": "chain_rule_o_derivative_predicted",
                    "inverse_condition": chain[1],
                    "retained_rows": chain[2],
                    "surface_relative_residual": chain[3],
                    "local_map_jacobian_condition": chain[4],
                    "local_map_fit_rmse_um": chain[5],
                    "taylor_design_rank": "",
                    "taylor_design_condition": "",
                    "taylor_map_relative_residual": "",
                    "taylor_root_rms_scaled": "",
                }
            )

            nominal_local_states = np.broadcast_to(
                local_nominal[:, None, :], (nb, nk, 2)
            )
            initial = estimates["fd_quartic_source_predicted"][target, realization]
            for order in (3, 4, 5):
                method = f"o_taylor_order{order}_nominal_local"
                result = raw_taylor_map_center(
                    values,
                    nominal_local_states,
                    delta_k2,
                    nominal_k2,
                    coordinate_scale,
                    order,
                    initial,
                )
                estimates[method][target, realization] = result[0]
                diagnostics.append(
                    {
                        "target": target_names[target],
                        "target_index": target + 1,
                        "realization": realization + 1,
                        "method": method,
                        "inverse_condition": "",
                        "retained_rows": "",
                        "surface_relative_residual": "",
                        "local_map_jacobian_condition": "",
                        "local_map_fit_rmse_um": "",
                        "taylor_design_rank": result[2],
                        "taylor_design_condition": result[1],
                        "taylor_map_relative_residual": result[3],
                        "taylor_root_rms_scaled": result[4],
                    }
                )
            all_state = raw_taylor_map_center(
                values,
                predicted_relative[target, realization],
                delta_k2,
                nominal_k2,
                coordinate_scale,
                4,
                initial,
            )
            estimates["o_taylor_order4_all_state_local"][target, realization] = all_state[0]
            diagnostics.append(
                {
                    "target": target_names[target],
                    "target_index": target + 1,
                    "realization": realization + 1,
                    "method": "o_taylor_order4_all_state_local",
                    "inverse_condition": "",
                    "retained_rows": "",
                    "surface_relative_residual": "",
                    "local_map_jacobian_condition": "",
                    "local_map_fit_rmse_um": "",
                    "taylor_design_rank": all_state[2],
                    "taylor_design_condition": all_state[1],
                    "taylor_map_relative_residual": all_state[3],
                    "taylor_root_rms_scaled": all_state[4],
                }
            )
        print(f"inverse {target + 1}/{nt}: {target_names[target]}")

    for method, values in estimates.items():
        np.save(output / f"{method}_relative_center_estimates.npy", values)
    write_rows(output / "fit_diagnostics.csv", diagnostics)

    # Evaluation-only truth begins here.
    target_orbits = np.load(source / "target_orbits.npy", mmap_mode="r")
    target_truth_increment = np.load(source / "target_truth.npy")
    nominal_target_centers = np.load(source / "nominal_target_centers.npy")
    exact_reference = np.asarray(
        target_orbits[:, :, zero_bump, nominal_k2, :], dtype=float
    )
    true_total_center = nominal_target_centers[:, None, :] + target_truth_increment
    relative_truth = true_total_center - exact_reference
    exact_relative = np.asarray(target_orbits, dtype=float) - exact_reference[:, :, None, None, :]

    # Add two evaluation-only oracle variants after every machine-facing fit.
    oracle_source = np.zeros((nt, nr, 2))
    oracle_taylor = np.zeros((nt, nr, 2))
    for target in range(nt):
        for realization in range(nr):
            values = bpm_flat[target, realization]
            slopes = polynomial_k2_slope(values, delta_k2, 4)
            exact_nominal = exact_relative[target, realization, :, nominal_k2, :]
            oracle_source[target, realization] = fit_center(slopes, exact_nominal)
            oracle_taylor[target, realization] = raw_taylor_map_center(
                values,
                exact_relative[target, realization],
                delta_k2,
                nominal_k2,
                coordinate_scale,
                4,
                oracle_source[target, realization],
            )[0]
    estimates["fd_quartic_source_oracle_local"] = oracle_source
    estimates["o_taylor_order4_oracle_all_state_local"] = oracle_taylor
    np.save(output / "fd_quartic_source_oracle_local_relative_center_estimates.npy", oracle_source)
    np.save(
        output / "o_taylor_order4_oracle_all_state_local_relative_center_estimates.npy",
        oracle_taylor,
    )

    local_nominal_errors = (
        predicted_relative[:, :, :, nominal_k2, :]
        - exact_relative[:, :, :, nominal_k2, :]
    )
    nonzero_bumps = np.arange(nb) != zero_bump
    local_all_errors = predicted_relative - exact_relative
    reference_errors = predicted_reference_absolute - exact_reference
    local_rows = [
        {
            "quantity": "relative_local_orbit_nominal_k2_nonzero_bumps",
            **summarize_vectors(local_nominal_errors[:, :, nonzero_bumps, :]),
        },
        {
            "quantity": "relative_local_orbit_all_states",
            **summarize_vectors(local_all_errors),
        },
        {
            "quantity": "absolute_zero_bump_reference_orbit",
            **summarize_vectors(reference_errors),
        },
    ]
    write_rows(output / "local_orbit_summary.csv", local_rows)

    summary_rows: list[dict[str, object]] = []
    per_case_rows: list[dict[str, object]] = []
    for method, relative_estimate in estimates.items():
        relative_error = relative_estimate - relative_truth
        absolute_increment_estimate = (
            relative_estimate
            + predicted_reference_absolute
            - nominal_target_centers[:, None, :]
        )
        absolute_increment_error = absolute_increment_estimate - target_truth_increment
        row = {
            "method": method,
            **{f"relative_{key}": value for key, value in summarize_vectors(relative_error).items()},
            **{
                f"absolute_increment_{key}": value
                for key, value in summarize_vectors(absolute_increment_error).items()
            },
        }
        summary_rows.append(row)
        for target, name in enumerate(target_names):
            for realization in range(nr):
                relative_vector_um = relative_error[target, realization] * 1e6
                absolute_vector_um = absolute_increment_error[target, realization] * 1e6
                per_case_rows.append(
                    {
                        "method": method,
                        "target": name,
                        "target_index": target + 1,
                        "realization": realization + 1,
                        "relative_truth_x_um": relative_truth[target, realization, 0] * 1e6,
                        "relative_truth_y_um": relative_truth[target, realization, 1] * 1e6,
                        "relative_estimate_x_um": relative_estimate[target, realization, 0] * 1e6,
                        "relative_estimate_y_um": relative_estimate[target, realization, 1] * 1e6,
                        "relative_error_x_um": relative_vector_um[0],
                        "relative_error_y_um": relative_vector_um[1],
                        "relative_error_2d_um": np.linalg.norm(relative_vector_um),
                        "absolute_increment_truth_x_um": target_truth_increment[target, realization, 0] * 1e6,
                        "absolute_increment_truth_y_um": target_truth_increment[target, realization, 1] * 1e6,
                        "absolute_increment_estimate_x_um": absolute_increment_estimate[target, realization, 0] * 1e6,
                        "absolute_increment_estimate_y_um": absolute_increment_estimate[target, realization, 1] * 1e6,
                        "absolute_increment_error_2d_um": np.linalg.norm(absolute_vector_um),
                    }
                )
    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "per_case_estimates.csv", per_case_rows)

    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    methods = [row["method"] for row in summary_rows]
    values = [float(row["relative_rmse_2d_um"]) for row in summary_rows]
    ax.barh(np.arange(len(methods)), values, color="#4472c4")
    ax.set_yticks(np.arange(len(methods)), [name.replace("_", " ") for name in methods])
    ax.invert_yaxis()
    ax.set_xlabel("Beam-relative center 2D RMSE [micrometers]")
    ax.set_title("Only-sextupole-misalignment inverse benchmark")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "relative_center_rmse_by_method.png", dpi=180)
    plt.close(fig)

    summary_lines = "\n".join(
        f"| {row['method']} | {row['relative_rmse_2d_um']:.6f} | "
        f"{row['relative_median_2d_um']:.6f} | {row['relative_p90_2d_um']:.6f} | "
        f"{row['absolute_increment_rmse_2d_um']:.6f} |"
        for row in summary_rows
    )
    local_lines = "\n".join(
        f"| {row['quantity']} | {row['rmse_2d_um']:.6f} | {row['median_2d_um']:.6f} | "
        f"{row['p90_2d_um']:.6f} | {row['max_2d_um']:.6f} |"
        for row in local_rows
    )
    report = f"""# Only-sextupole-misalignment BPM/Taylor-map results

All forward states use the validated latest repaired SciBmad lattice.  The only
machine error is fixed x/y misalignment on all 76 active normal sextupoles.
BPM errors, time drift, corrector/K2 calibration errors, and quadrupole errors
are absent.

- targets / latent realizations: {nt} / {nr} per target
- exact states: {nt * nr * nb * nk}
- bump grid: {nb} points, amplitude +/-{coordinate_scale * 1e3:.3f} mm per plane
- K2 grid: {nk} points, delta K2 range {delta_k2.min():.6g} to {delta_k2.max():.6g} m^-3
- BPM count: {nd}

## BPM-predicted local orbit

| quantity | 2D RMSE [um] | median [um] | P90 [um] | maximum [um] |
|---|---:|---:|---:|---:|
{local_lines}

## Center inverse

| method | beam-relative 2D RMSE [um] | median [um] | P90 [um] | absolute-increment 2D RMSE [um] |
|---|---:|---:|---:|---:|
{summary_lines}

`fd_*_source_predicted` is the maintained physical two-source inverse after a
linear or quartic-in-K2 derivative extraction.  `quadratic_o_derivative` fits
the observable K2 derivative directly as a quadratic function of the
BPM-predicted target orbit.  `chain_rule_o_derivative` instead obtains the
local orbit Jacobian with respect to the commanded bump and transforms the
observable derivatives by the chain rule.  `o_taylor_orderN` fits all retained
raw nonzero-K2 states to a total-order-N polynomial in local x, y, and K2 and
finds the common zero of its analytic K2 derivative.

Oracle-local rows use exact target coordinates only after all machine-facing
fits and are evaluation diagnostics, not deployable methods.
"""
    (output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
