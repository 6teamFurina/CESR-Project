#!/usr/bin/env python3
"""Independent inventory and consistency checks for the validity envelope."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
LOCAL_OUTPUT_COLUMNS = [
    f"{location}_{coordinate}"
    for location in ("entry", "exit")
    for coordinate in ("x", "px", "y", "py", "z", "pz")
]


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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-dir", type=Path, default=HERE / "results" / "exact_validation")
    parser.add_argument("--analysis-dir", type=Path, default=HERE / "results" / "analysis")
    parser.add_argument("--maps-root", type=Path, default=HERE / "results" / "gtpsa_maps")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    exact = pd.read_csv(args.exact_dir / "exact_local_orbit_states.csv")
    exact["converged"] = exact["converged"].astype(str).str.lower().eq("true")
    states = pd.read_csv(args.analysis_dir / "state_taylor_errors.csv")
    states["converged"] = states["converged"].astype(str).str.lower().eq("true")
    states["passes_taylor_gate"] = (
        states["passes_taylor_gate"].astype(str).str.lower().eq("true")
    )
    families = pd.read_csv(args.analysis_dir / "family_validity_limits.csv")
    limits = pd.read_csv(args.analysis_dir / "per_sextupole_limits.csv")
    compact = pd.read_csv(args.analysis_dir / "per_sextupole_order4_envelope.csv")
    diagnostics = pd.read_csv(args.analysis_dir / "map_diagnostics.csv")
    summary = json.loads((args.analysis_dir / "summary.json").read_text(encoding="utf-8"))
    map_source_by_target = dict(zip(limits["target"], limits["map_source"]))

    targets = list(dict.fromkeys(exact["target"]))
    require(len(targets) == 76, f"Expected 76 targets, found {len(targets)}")
    counts = exact.groupby("target").size()
    require((counts == 651).all(), f"Unexpected exact state counts: {counts.value_counts().to_dict()}")
    require(len(exact) == 76 * 651, "Exact row inventory mismatch")
    require(np.isfinite(exact.loc[exact["converged"], [
        "entry_x", "entry_px", "entry_y", "entry_py",
        "exit_x", "exit_px", "exit_y", "exit_py",
    ]].to_numpy(float)).all(), "Non-finite converged local orbit")

    require(len(limits) == 76, "Per-target limit inventory mismatch")
    require(set(limits["target"]) == set(targets), "Limit targets differ from exact targets")
    require(len(compact) == 76, "Compact order-four inventory mismatch")
    require(set(compact["target"]) == set(targets), "Compact targets differ from exact targets")
    require(len(diagnostics) == 76, "Map diagnostic inventory mismatch")
    require(summary["target_count"] == 76, "Summary target count mismatch")
    require(
        summary["direct_gtpsa_target_count"] + summary["fallback_target_count"] == 76,
        "Direct/fallback map counts do not sum to 76",
    )
    direct_directories = {
        path.name.upper()
        for path in args.maps_root.iterdir()
        if path.is_dir() and (path / "map_metadata.toml").exists()
    }
    direct_targets = set(limits.loc[limits["map_source"] == "direct_gtpsa", "target"])
    require(
        direct_targets == direct_directories,
        "Direct-map directory/source mismatch: "
        f"analysis-only={sorted(direct_targets - direct_directories)}, "
        f"map-only={sorted(direct_directories - direct_targets)}",
    )

    require(set(states["taylor_order"]) == {1, 2, 3, 4}, "Taylor-order inventory mismatch")
    require(len(states) == 4 * len(exact), "State-error inventory mismatch")
    require((states["position_error_um"] >= 0).all(), "Negative position error")
    require((states["angle_error_urad"] >= 0).all(), "Negative angle error")

    compact_brackets = [
        (
            "x_bump_last_pass_mm_at_abs_delta_k2_0p1",
            "x_bump_first_fail_mm_at_abs_delta_k2_0p1",
        ),
        (
            "y_bump_last_pass_mm_at_abs_delta_k2_0p1",
            "y_bump_first_fail_mm_at_abs_delta_k2_0p1",
        ),
        (
            "abs_delta_k2_last_pass_m3_at_abs_x_bump_1p5mm",
            "abs_delta_k2_first_fail_m3_at_abs_x_bump_1p5mm",
        ),
        (
            "abs_delta_k2_last_pass_m3_at_abs_y_bump_1p5mm",
            "abs_delta_k2_first_fail_m3_at_abs_y_bump_1p5mm",
        ),
        (
            "joint_last_pass_scale_of_1p5mm_0p1m3",
            "joint_first_fail_scale_of_1p5mm_0p1m3",
        ),
    ]
    for passing_column, failing_column in compact_brackets:
        located = compact[failing_column].notna()
        require(
            (
                compact.loc[located, failing_column]
                > compact.loc[located, passing_column]
            ).all(),
            f"Non-increasing compact bracket: {passing_column}/{failing_column}",
        )

    # Recheck the maintained eight signed states at radius scale one.
    protocol = states[
        states["direction_family"].isin(["xk_protocol", "yk_protocol"])
        & np.isclose(states["radius_scale"], 1.0)
        & (states["taylor_order"] == 4)
    ]
    require(len(protocol) == 76 * 8, "Maintained protocol state inventory mismatch")
    require(protocol["converged"].all(), "At least one maintained protocol state did not converge")
    require(protocol["passes_taylor_gate"].all(), "Maintained protocol exceeds order-four gate")

    # Independent fallback check: fit only inside a small normalized cube and
    # predict the maintained protocol outside that cube.  None of the tested
    # q_bump=+/-1.5 maintained-protocol points enters this fit.
    powers = powers_through(4)
    cross_validation_position_errors = []
    cross_validation_angle_errors = []
    cross_validation_training_counts = []
    cross_validation_rows = []
    for target_name, target in exact.groupby("target", sort=False):
        q = target[["qx", "qy", "qk"]].to_numpy(float)
        training_mask = target["converged"].to_numpy(bool) & (
            np.max(np.abs(q), axis=1) <= 1.0 + 1e-12
        )
        training_design = design_matrix(q[training_mask], powers)
        require(
            np.linalg.matrix_rank(training_design) == len(powers),
            f"Small-cube cross-validation design is rank deficient for {target_name}",
        )
        coefficients, *_ = np.linalg.lstsq(
            training_design,
            target.loc[training_mask, LOCAL_OUTPUT_COLUMNS].to_numpy(float),
            rcond=None,
        )
        test_mask = (
            target["direction_family"].isin(["xk_protocol", "yk_protocol"]).to_numpy(bool)
            & np.isclose(target["radius_scale"].to_numpy(float), 1.0)
        )
        require(target.loc[test_mask, "converged"].all(), f"Cross-validation test failed for {target_name}")
        residual = (
            design_matrix(q[test_mask], powers) @ coefficients
            - target.loc[test_mask, LOCAL_OUTPUT_COLUMNS].to_numpy(float)
        )
        target_position_errors = (
            np.max(np.abs(residual[:, [0, 2, 6, 8]]), axis=1) * 1e6
        )
        target_angle_errors = (
            np.max(np.abs(residual[:, [1, 3, 7, 9]]), axis=1) * 1e6
        )
        cross_validation_position_errors.extend(target_position_errors)
        cross_validation_angle_errors.extend(target_angle_errors)
        cross_validation_training_counts.append(int(training_mask.sum()))
        cross_validation_rows.append(
            {
                "target": target_name,
                "map_source": map_source_by_target[target_name],
                "normalized_training_cube_bound": 1.0,
                "training_state_count": int(training_mask.sum()),
                "training_design_rank": int(np.linalg.matrix_rank(training_design)),
                "held_out_protocol_state_count": int(test_mask.sum()),
                "maximum_position_error_um": float(np.max(target_position_errors)),
                "maximum_angle_error_urad": float(np.max(target_angle_errors)),
                "passes_1um_1urad_gate": bool(
                    np.max(target_position_errors) <= 1.0
                    and np.max(target_angle_errors) <= 1.0
                ),
            }
        )

    cv_position_max = float(np.max(cross_validation_position_errors))
    cv_angle_max = float(np.max(cross_validation_angle_errors))
    require(cv_position_max <= 1.0, "Small-cube cross-validation exceeds 1 um")
    require(cv_angle_max <= 1.0, "Small-cube cross-validation exceeds 1 urad")

    # The stored family limit must equal the prefix of radii passing the gate.
    order4 = states[states["taylor_order"] == 4]
    for row in families[families["taylor_order"] == 4].itertuples(index=False):
        group = order4[
            (order4["target"] == row.target)
            & (order4["direction_family"] == row.direction_family)
        ]
        by_radius = group.groupby("radius_scale")["passes_taylor_gate"].all().sort_index()
        expected = 0.0
        for radius, passed in by_radius.items():
            if not passed:
                break
            expected = float(radius)
        require(
            np.isclose(expected, float(row.last_passing_radius)),
            f"Family prefix limit mismatch for {row.target} {row.direction_family}",
        )

    validation_summary = {
        "validated_targets": 76,
        "exact_states": len(exact),
        "converged_exact_states": int(exact["converged"].sum()),
        "direct_gtpsa_targets": summary["direct_gtpsa_target_count"],
        "fallback_targets": summary["fallback_target_count"],
        "maintained_order4_protocol_states": len(protocol),
        "small_cube_cv_bound": 1.0,
        "small_cube_cv_training_states_per_target_min": min(
            cross_validation_training_counts
        ),
        "small_cube_cv_training_states_per_target_max": max(
            cross_validation_training_counts
        ),
        "small_cube_cv_protocol_position_max_um": cv_position_max,
        "small_cube_cv_protocol_angle_max_urad": cv_angle_max,
    }
    pd.DataFrame(cross_validation_rows).to_csv(
        args.analysis_dir / "small_cube_cross_validation.csv", index=False
    )
    output_path = args.output or (args.analysis_dir / "validation_summary.json")
    output_path.write_text(json.dumps(validation_summary, indent=2), encoding="utf-8")
    print(json.dumps(validation_summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
