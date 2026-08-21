#!/usr/bin/env python3
"""Compare local-orbit Taylor maps with exact SciBmad amplitude rays."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DEFAULT_EXACT = HERE / "results" / "exact_validation"
DEFAULT_MAPS = HERE / "results" / "gtpsa_maps"
DEFAULT_OUTPUT = HERE / "results" / "analysis"
COORDINATES = ("x", "px", "y", "py", "z", "pz")
LOCATIONS = ("entry", "exit")
POSITION_COLUMNS = ("entry_x", "entry_y", "exit_x", "exit_y")
ANGLE_COLUMNS = ("entry_px", "entry_py", "exit_px", "exit_py")


def powers_through(order: int) -> list[tuple[int, int, int]]:
    return [
        (px, py, total - px - py)
        for total in range(order + 1)
        for px in range(total + 1)
        for py in range(total - px + 1)
    ]


def design_matrix(q: np.ndarray, powers: list[tuple[int, int, int]]) -> np.ndarray:
    return np.column_stack(
        [q[:, 0] ** px * q[:, 1] ** py * q[:, 2] ** pk for px, py, pk in powers]
    )


def load_direct_coefficients(path: Path, powers: list[tuple[int, int, int]]) -> np.ndarray:
    table = pd.read_csv(path)
    coefficients = np.full((len(LOCATIONS) * len(COORDINATES), len(powers)), np.nan)
    power_index = {power: index for index, power in enumerate(powers)}
    output_index = {
        (location, coordinate): i
        for i, (location, coordinate) in enumerate(
            (pair for location in LOCATIONS for pair in ((location, c) for c in COORDINATES))
        )
    }
    for row in table.itertuples(index=False):
        power = (int(row.qx_power), int(row.qy_power), int(row.qk_power))
        if power in power_index:
            coefficients[output_index[(row.location, row.coordinate)], power_index[power]] = (
                float(row.coefficient)
            )
    if not np.all(np.isfinite(coefficients)):
        raise ValueError(f"Incomplete direct coefficient table: {path}")
    return coefficients


def fit_empirical_coefficients(
    target: pd.DataFrame,
    powers: list[tuple[int, int, int]],
    training_radius: float,
) -> tuple[np.ndarray, dict[str, float]]:
    train = target[target["converged"] & (target["radius_scale"] <= training_radius)].copy()
    q = train[["qx", "qy", "qk"]].to_numpy(float)
    design = design_matrix(q, powers)
    rank = int(np.linalg.matrix_rank(design))
    if rank != len(powers):
        raise ValueError(
            f"Empirical Taylor design rank {rank}/{len(powers)} for {target.iloc[0]['target']}"
        )
    outputs = train[[f"{location}_{coordinate}" for location in LOCATIONS for coordinate in COORDINATES]].to_numpy(float)
    coefficients, *_ = np.linalg.lstsq(design, outputs, rcond=None)
    prediction = design @ coefficients
    residual = prediction - outputs
    return coefficients.T, {
        "training_state_count": len(train),
        "training_design_rank": rank,
        "training_condition_number": float(np.linalg.cond(design)),
        "training_position_rmse_um": float(
            np.sqrt(np.mean(residual[:, [0, 2, 6, 8]] ** 2)) * 1e6
        ),
        "training_angle_rmse_urad": float(
            np.sqrt(np.mean(residual[:, [1, 3, 7, 9]] ** 2)) * 1e6
        ),
    }


def prediction(coefficients: np.ndarray, q: np.ndarray, powers, order: int) -> np.ndarray:
    selected = [index for index, power in enumerate(powers) if sum(power) <= order]
    return design_matrix(q, [powers[index] for index in selected]) @ coefficients[:, selected].T


def finite_max(values: pd.Series, scale: float = 1.0) -> float:
    """Return a scaled finite maximum without warning on an all-NaN ray."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(scale * np.max(finite)) if finite.size else math.nan


def prefix_limit(group: pd.DataFrame, pass_column: str) -> dict[str, float | bool]:
    by_radius = (
        group.groupby("radius_scale", sort=True)
        .agg(
            passed=(pass_column, "all"),
            all_converged=("converged", "all"),
            worst_position_error_um=("position_error_um", "max"),
            worst_angle_error_urad=("angle_error_urad", "max"),
            maximum_ring_abs_x_mm=("maximum_ring_abs_x_m", lambda x: finite_max(x, 1e3)),
            maximum_ring_abs_y_mm=("maximum_ring_abs_y_m", lambda x: finite_max(x, 1e3)),
            maximum_abs_corrector_delta=("maximum_abs_corrector_delta", "max"),
        )
        .reset_index()
    )
    passing_radius = 0.0
    first_failure = math.nan
    for row in by_radius.itertuples(index=False):
        if bool(row.passed):
            passing_radius = float(row.radius_scale)
        else:
            first_failure = float(row.radius_scale)
            break
    censored = math.isnan(first_failure)
    return {
        "last_passing_radius": passing_radius,
        "first_failing_radius": first_failure,
        "scan_cap_lower_bound": censored,
        "maximum_tested_radius": float(by_radius["radius_scale"].max()),
    }


def physical_limit(family: str, radius: float, bump_scale: float, k2_scale: float) -> float:
    if family in {"x_axis", "y_axis", "x_bump_fixed_k", "y_bump_fixed_k"}:
        return radius * bump_scale * 1e3
    if family in {"k_axis", "k_fixed_x_bump", "k_fixed_y_bump"}:
        return radius * k2_scale
    return radius


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-dir", type=Path, default=DEFAULT_EXACT)
    parser.add_argument("--maps-root", type=Path, default=DEFAULT_MAPS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--order", type=int, default=4)
    parser.add_argument("--training-radius", type=float, default=1.5)
    parser.add_argument("--position-tolerance-um", type=float, default=1.0)
    parser.add_argument("--angle-tolerance-urad", type=float, default=1.0)
    args = parser.parse_args()
    if args.order < 2:
        raise ValueError("Taylor order must be at least two")

    exact = pd.read_csv(args.exact_dir / "exact_local_orbit_states.csv")
    exact["converged"] = exact["converged"].astype(str).str.lower().eq("true")
    metadata_text = (args.exact_dir / "scan_metadata.toml").read_text(encoding="utf-8")
    # These scales are intentionally parsed without adding a TOML dependency.
    bump_scale = float(
        next(line.split("=", 1)[1] for line in metadata_text.splitlines() if line.startswith("bump_scale_m"))
    )
    k2_scale = float(
        next(line.split("=", 1)[1] for line in metadata_text.splitlines() if line.startswith("k2_scale_m3"))
    )
    powers = powers_through(args.order)
    state_outputs: list[pd.DataFrame] = []
    family_outputs: list[dict] = []
    target_outputs: list[dict] = []
    map_diagnostics: list[dict] = []

    for target_name, target in exact.groupby("target", sort=False):
        target = target.copy().reset_index(drop=True)
        direct_path = (
            args.maps_root / target_name.lower() / "local_orbit_taylor_coefficients.csv"
        )
        fit_coefficients, fit_diagnostics = fit_empirical_coefficients(
            target, powers, args.training_radius
        )
        if direct_path.exists():
            coefficients = load_direct_coefficients(direct_path, powers)
            source = "direct_gtpsa"
        else:
            coefficients = fit_coefficients
            source = "exact_scan_fitted_fallback"

        q = target[["qx", "qy", "qk"]].to_numpy(float)
        exact_outputs = target[
            [f"{location}_{coordinate}" for location in LOCATIONS for coordinate in COORDINATES]
        ].to_numpy(float)
        diagnostic = {"target": target_name, "map_source": source, **fit_diagnostics}
        if direct_path.exists():
            coefficient_difference = coefficients - fit_coefficients
            diagnostic["direct_vs_fit_coefficient_relative_l2"] = float(
                np.linalg.norm(coefficient_difference) / max(np.linalg.norm(coefficients), 1e-30)
            )
        else:
            diagnostic["direct_vs_fit_coefficient_relative_l2"] = math.nan
        map_diagnostics.append(diagnostic)

        order_limits: dict[int, dict[str, dict]] = {}
        for order in range(1, args.order + 1):
            predicted = prediction(coefficients, q, powers, order)
            residual = predicted - exact_outputs
            target_order = target.copy()
            target_order["map_source"] = source
            target_order["taylor_order"] = order
            converged = target_order["converged"].to_numpy(bool)
            position_error = np.full(len(target_order), np.inf)
            angle_error = np.full(len(target_order), np.inf)
            position_error[converged] = (
                np.max(np.abs(residual[converged][:, [0, 2, 6, 8]]), axis=1) * 1e6
            )
            angle_error[converged] = (
                np.max(np.abs(residual[converged][:, [1, 3, 7, 9]]), axis=1) * 1e6
            )
            target_order["position_error_um"] = position_error
            target_order["angle_error_urad"] = angle_error
            target_order["passes_taylor_gate"] = (
                target_order["converged"]
                & (target_order["position_error_um"] <= args.position_tolerance_um)
                & (target_order["angle_error_urad"] <= args.angle_tolerance_urad)
            )
            state_outputs.append(target_order)
            order_limits[order] = {}
            for family, family_table in target_order.groupby("direction_family", sort=False):
                if family == "origin":
                    continue
                limit = prefix_limit(family_table, "passes_taylor_gate")
                limit["target"] = target_name
                limit["map_source"] = source
                limit["taylor_order"] = order
                limit["direction_family"] = family
                limit["last_passing_physical_limit"] = physical_limit(
                    family, float(limit["last_passing_radius"]), bump_scale, k2_scale
                )
                family_outputs.append(limit)
                order_limits[order][family] = limit

        def limit(order: int, family: str) -> float:
            return float(order_limits[order][family]["last_passing_physical_limit"])

        def radius(order: int, family: str) -> float:
            return float(order_limits[order][family]["last_passing_radius"])

        def first_failing_physical(order: int, family: str) -> float:
            failing_radius = float(order_limits[order][family]["first_failing_radius"])
            if math.isnan(failing_radius):
                return math.nan
            return physical_limit(family, failing_radius, bump_scale, k2_scale)

        def maximum_corrector_demand(families: tuple[str, ...], at_radius: float) -> float:
            if at_radius <= 0:
                return 0.0
            selected = target[
                target["direction_family"].isin(families)
                & np.isclose(target["radius_scale"], at_radius)
            ]["maximum_abs_corrector_delta"].to_numpy(float)
            finite = selected[np.isfinite(selected)]
            return float(np.max(finite)) if finite.size else math.nan

        summary = {
            "target": target_name,
            "target_s_m": float(target.iloc[0]["target_s_m"]),
            "nominal_k2_m3": float(target.iloc[0]["target_nominal_kn2_m3"]),
            "map_source": source,
            "all_exact_states_converged": bool(target["converged"].all()),
            "failed_exact_state_count": int((~target["converged"]).sum()),
            "order1_x_bump_zero_k2_limit_mm": limit(1, "x_axis"),
            "order1_y_bump_zero_k2_limit_mm": limit(1, "y_axis"),
            "order1_x_bump_at_abs_k2_0p1_limit_mm": limit(1, "x_bump_fixed_k"),
            "order1_y_bump_at_abs_k2_0p1_limit_mm": limit(1, "y_bump_fixed_k"),
            "order1_k2_at_abs_x_bump_1p5mm_limit_m3": limit(1, "k_fixed_x_bump"),
            "order1_k2_at_abs_y_bump_1p5mm_limit_m3": limit(1, "k_fixed_y_bump"),
            "order2_x_bump_zero_k2_limit_mm": limit(2, "x_axis"),
            "order2_y_bump_zero_k2_limit_mm": limit(2, "y_axis"),
            "order2_x_bump_at_abs_k2_0p1_limit_mm": limit(2, "x_bump_fixed_k"),
            "order2_y_bump_at_abs_k2_0p1_limit_mm": limit(2, "y_bump_fixed_k"),
            "order2_k2_at_abs_x_bump_1p5mm_limit_m3": limit(2, "k_fixed_x_bump"),
            "order2_k2_at_abs_y_bump_1p5mm_limit_m3": limit(2, "k_fixed_y_bump"),
            "order2_xk_joint_scale_limit": radius(2, "xk_protocol"),
            "order2_yk_joint_scale_limit": radius(2, "yk_protocol"),
            "order2_common_joint_scale_limit": min(
                radius(2, "xk_protocol"), radius(2, "yk_protocol")
            ),
            "order4_x_bump_zero_k2_limit_mm": limit(args.order, "x_axis"),
            "order4_y_bump_zero_k2_limit_mm": limit(args.order, "y_axis"),
            "order4_x_bump_at_abs_k2_0p1_limit_mm": limit(args.order, "x_bump_fixed_k"),
            "order4_y_bump_at_abs_k2_0p1_limit_mm": limit(args.order, "y_bump_fixed_k"),
            "order4_k2_at_abs_x_bump_1p5mm_limit_m3": limit(args.order, "k_fixed_x_bump"),
            "order4_k2_at_abs_y_bump_1p5mm_limit_m3": limit(args.order, "k_fixed_y_bump"),
            "order4_xk_joint_scale_limit": radius(args.order, "xk_protocol"),
            "order4_yk_joint_scale_limit": radius(args.order, "yk_protocol"),
            "order4_common_joint_scale_limit": min(
                radius(args.order, "xk_protocol"), radius(args.order, "yk_protocol")
            ),
            "order4_x_bump_at_abs_k2_0p1_first_fail_mm": first_failing_physical(
                args.order, "x_bump_fixed_k"
            ),
            "order4_y_bump_at_abs_k2_0p1_first_fail_mm": first_failing_physical(
                args.order, "y_bump_fixed_k"
            ),
            "order4_k2_at_abs_x_bump_1p5mm_first_fail_m3": first_failing_physical(
                args.order, "k_fixed_x_bump"
            ),
            "order4_k2_at_abs_y_bump_1p5mm_first_fail_m3": first_failing_physical(
                args.order, "k_fixed_y_bump"
            ),
            "order4_x_bump_limit_is_scan_cap_lower_bound": bool(
                order_limits[args.order]["x_bump_fixed_k"]["scan_cap_lower_bound"]
            ),
            "order4_y_bump_limit_is_scan_cap_lower_bound": bool(
                order_limits[args.order]["y_bump_fixed_k"]["scan_cap_lower_bound"]
            ),
            "order4_k2_at_x_bump_limit_is_scan_cap_lower_bound": bool(
                order_limits[args.order]["k_fixed_x_bump"]["scan_cap_lower_bound"]
            ),
            "order4_k2_at_y_bump_limit_is_scan_cap_lower_bound": bool(
                order_limits[args.order]["k_fixed_y_bump"]["scan_cap_lower_bound"]
            ),
        }
        joint_failures = [
            float(order_limits[args.order][family]["first_failing_radius"])
            for family in ("xk_protocol", "yk_protocol")
        ]
        finite_joint_failures = [value for value in joint_failures if math.isfinite(value)]
        summary["order4_common_joint_first_failing_scale"] = (
            min(finite_joint_failures) if finite_joint_failures else math.nan
        )
        summary["order4_common_joint_bump_limit_mm"] = (
            1.5 * summary["order4_common_joint_scale_limit"]
        )
        summary["order4_common_joint_k2_limit_m3"] = (
            k2_scale * summary["order4_common_joint_scale_limit"]
        )
        summary["order4_common_joint_first_failing_bump_mm"] = (
            1.5 * summary["order4_common_joint_first_failing_scale"]
        )
        summary["order4_common_joint_first_failing_k2_m3"] = (
            k2_scale * summary["order4_common_joint_first_failing_scale"]
        )
        summary["order4_x_bump_limit_max_abs_corrector_delta"] = maximum_corrector_demand(
            ("x_bump_fixed_k",), radius(args.order, "x_bump_fixed_k")
        )
        summary["order4_y_bump_limit_max_abs_corrector_delta"] = maximum_corrector_demand(
            ("y_bump_fixed_k",), radius(args.order, "y_bump_fixed_k")
        )
        summary["order4_joint_limit_max_abs_corrector_delta"] = maximum_corrector_demand(
            ("xk_protocol", "yk_protocol"), summary["order4_common_joint_scale_limit"]
        )
        summary["order2_common_joint_bump_limit_mm"] = (
            1.5 * summary["order2_common_joint_scale_limit"]
        )
        summary["order2_common_joint_k2_limit_m3"] = (
            k2_scale * summary["order2_common_joint_scale_limit"]
        )
        target_outputs.append(summary)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    states = pd.concat(state_outputs, ignore_index=True)
    families = pd.DataFrame(family_outputs)
    targets = pd.DataFrame(target_outputs)
    diagnostics = pd.DataFrame(map_diagnostics)
    states.to_csv(args.output_dir / "state_taylor_errors.csv", index=False)
    families.to_csv(args.output_dir / "family_validity_limits.csv", index=False)
    targets.to_csv(args.output_dir / "per_sextupole_limits.csv", index=False)
    compact_order4_columns = {
        "target": "target",
        "target_s_m": "target_s_m",
        "nominal_k2_m3": "nominal_k2_m3",
        "map_source": "map_source",
        "failed_exact_state_count": "failed_outer_scan_state_count",
        "order4_x_bump_at_abs_k2_0p1_limit_mm": "x_bump_last_pass_mm_at_abs_delta_k2_0p1",
        "order4_x_bump_at_abs_k2_0p1_first_fail_mm": "x_bump_first_fail_mm_at_abs_delta_k2_0p1",
        "order4_x_bump_limit_is_scan_cap_lower_bound": "x_bump_last_pass_is_scan_cap_lower_bound",
        "order4_y_bump_at_abs_k2_0p1_limit_mm": "y_bump_last_pass_mm_at_abs_delta_k2_0p1",
        "order4_y_bump_at_abs_k2_0p1_first_fail_mm": "y_bump_first_fail_mm_at_abs_delta_k2_0p1",
        "order4_y_bump_limit_is_scan_cap_lower_bound": "y_bump_last_pass_is_scan_cap_lower_bound",
        "order4_k2_at_abs_x_bump_1p5mm_limit_m3": "abs_delta_k2_last_pass_m3_at_abs_x_bump_1p5mm",
        "order4_k2_at_abs_x_bump_1p5mm_first_fail_m3": "abs_delta_k2_first_fail_m3_at_abs_x_bump_1p5mm",
        "order4_k2_at_x_bump_limit_is_scan_cap_lower_bound": "k2_at_x_bump_last_pass_is_scan_cap_lower_bound",
        "order4_k2_at_abs_y_bump_1p5mm_limit_m3": "abs_delta_k2_last_pass_m3_at_abs_y_bump_1p5mm",
        "order4_k2_at_abs_y_bump_1p5mm_first_fail_m3": "abs_delta_k2_first_fail_m3_at_abs_y_bump_1p5mm",
        "order4_k2_at_y_bump_limit_is_scan_cap_lower_bound": "k2_at_y_bump_last_pass_is_scan_cap_lower_bound",
        "order4_common_joint_scale_limit": "joint_last_pass_scale_of_1p5mm_0p1m3",
        "order4_common_joint_first_failing_scale": "joint_first_fail_scale_of_1p5mm_0p1m3",
        "order4_common_joint_bump_limit_mm": "joint_last_pass_bump_mm",
        "order4_common_joint_k2_limit_m3": "joint_last_pass_abs_delta_k2_m3",
        "order4_common_joint_first_failing_bump_mm": "joint_first_fail_bump_mm",
        "order4_common_joint_first_failing_k2_m3": "joint_first_fail_abs_delta_k2_m3",
        "order4_x_bump_limit_max_abs_corrector_delta": "x_limit_max_abs_corrector_field_delta",
        "order4_y_bump_limit_max_abs_corrector_delta": "y_limit_max_abs_corrector_field_delta",
        "order4_joint_limit_max_abs_corrector_delta": "joint_limit_max_abs_corrector_field_delta",
    }
    targets[list(compact_order4_columns)].rename(columns=compact_order4_columns).to_csv(
        args.output_dir / "per_sextupole_order4_envelope.csv", index=False
    )
    diagnostics.to_csv(args.output_dir / "map_diagnostics.csv", index=False)

    numeric_columns = [
        "order2_x_bump_at_abs_k2_0p1_limit_mm",
        "order2_y_bump_at_abs_k2_0p1_limit_mm",
        "order2_k2_at_abs_x_bump_1p5mm_limit_m3",
        "order2_k2_at_abs_y_bump_1p5mm_limit_m3",
        "order2_common_joint_scale_limit",
        "order4_x_bump_at_abs_k2_0p1_limit_mm",
        "order4_y_bump_at_abs_k2_0p1_limit_mm",
        "order4_k2_at_abs_x_bump_1p5mm_limit_m3",
        "order4_k2_at_abs_y_bump_1p5mm_limit_m3",
        "order4_common_joint_scale_limit",
    ]
    maintained_protocol = states[
        states["direction_family"].isin(["xk_protocol", "yk_protocol"])
        & np.isclose(states["radius_scale"], 1.0)
    ]
    order2_protocol = maintained_protocol[maintained_protocol["taylor_order"] == 2]
    order4_protocol = maintained_protocol[maintained_protocol["taylor_order"] == 4]
    direct_order4_protocol = order4_protocol[order4_protocol["map_source"] == "direct_gtpsa"]
    direct_coefficient_differences = diagnostics.loc[
        diagnostics["map_source"] == "direct_gtpsa",
        "direct_vs_fit_coefficient_relative_l2",
    ].to_numpy(float)
    aggregate = {
        "format": "cesr-sextupole-local-orbit-validity-envelope-analysis-v1",
        "engine": "direct SciBmad/GTPSA maps with exact SciBmad validation; exact-scan Taylor fallback when direct GTPSA is unavailable",
        "target_count": len(targets),
        "exact_state_count": int(len(exact)),
        "converged_exact_state_count": int(exact["converged"].sum()),
        "failed_exact_state_count": int((~exact["converged"]).sum()),
        "targets_with_any_exact_failure": int(
            exact.groupby("target")["converged"].all().eq(False).sum()
        ),
        "direct_gtpsa_target_count": int((targets["map_source"] == "direct_gtpsa").sum()),
        "fallback_target_count": int((targets["map_source"] != "direct_gtpsa").sum()),
        "direct_gtpsa_target_inventory": targets.loc[
            targets["map_source"] == "direct_gtpsa", "target"
        ].tolist(),
        "fallback_target_inventory": targets.loc[
            targets["map_source"] != "direct_gtpsa", "target"
        ].tolist(),
        "position_tolerance_um": args.position_tolerance_um,
        "angle_tolerance_urad": args.angle_tolerance_urad,
        "maximum_taylor_order": args.order,
        "bump_scale_m": bump_scale,
        "k2_scale_m3": k2_scale,
        "fallback_training_radius_scale": args.training_radius,
        "gate_semantics": (
            "every signed state at a radius must converge and have maximum target-entry/exit "
            "transverse position error <= tolerance and slope error <= tolerance"
        ),
        "maintained_protocol": {
            "signed_state_count": int(len(order4_protocol)),
            "order2_passing_state_count": int(order2_protocol["passes_taylor_gate"].sum()),
            "order4_passing_state_count": int(order4_protocol["passes_taylor_gate"].sum()),
            "direct_order4_signed_state_count": int(len(direct_order4_protocol)),
            "direct_order4_maximum_position_error_um": float(
                direct_order4_protocol["position_error_um"].max()
            ),
            "direct_order4_maximum_angle_error_urad": float(
                direct_order4_protocol["angle_error_urad"].max()
            ),
        },
        "map_diagnostics": {
            "direct_vs_fit_coefficient_relative_l2_minimum": float(
                np.min(direct_coefficient_differences)
            ),
            "direct_vs_fit_coefficient_relative_l2_median": float(
                np.median(direct_coefficient_differences)
            ),
            "direct_vs_fit_coefficient_relative_l2_maximum": float(
                np.max(direct_coefficient_differences)
            ),
            "fit_training_position_rmse_um_maximum": float(
                diagnostics["training_position_rmse_um"].max()
            ),
            "fit_training_angle_rmse_urad_maximum": float(
                diagnostics["training_angle_rmse_urad"].max()
            ),
        },
        "aggregate_limits": {
            column: {
                "minimum": float(targets[column].min()),
                "median": float(targets[column].median()),
                "p90": float(targets[column].quantile(0.9)),
                "maximum": float(targets[column].max()),
            }
            for column in numeric_columns
        },
        "interpretation_boundary": (
            "Taylor/model-validity limits only. No CESR corrector, sextupole-power-supply, "
            "aperture, lifetime, interlock, or operator limit is asserted. A value at the "
            "largest scanned amplitude is a lower bound, not a located boundary."
        ),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(aggregate, indent=2), encoding="utf-8"
    )
    print(json.dumps(aggregate, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
