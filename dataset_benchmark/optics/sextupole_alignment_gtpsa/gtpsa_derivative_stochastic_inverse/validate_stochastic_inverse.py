#!/usr/bin/env python3
"""Independent consistency checks for the stochastic dO/dK2 inverse.

The validation deliberately reconstructs the explicit read-by-read random-walk
covariance for a short acquisition sequence.  This checks the closed form used
by the production analysis without allocating the 32,768-read default scan.
"""

from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

import numpy as np

import analyze_stochastic_inverse as inverse


HERE = Path(__file__).resolve().parent
PHYSICAL_ROOT = HERE / "results" / "exact_k5_b3"
MODEL_DIR = HERE.parent / "finite_bpm_inversion" / "results" / "local_orbit_model"
OUTPUT_DIR = HERE / "results" / "analysis"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def explicit_random_walk_covariance(
    design: np.ndarray,
    response: np.ndarray,
    states: list[tuple[int, int, int, int]],
    repeats: int,
    normalization: float,
    endpoint_rms_m: float,
) -> np.ndarray:
    """Build one center covariance by explicitly enumerating acquisition reads."""
    channels = response.shape[-1]
    left_inverse = np.linalg.inv(design.T @ design) @ design.T
    blocks = (left_inverse[:, :channels], left_inverse[:, channels:])
    one_cycle = np.asarray(
        [
            sign * (blocks[block] @ response[bump, k2]) / normalization
            for block, sign, bump, k2 in states
        ]
    )
    read_weights = np.tile(one_cycle, (repeats, 1)) / repeats
    reverse_tails = np.cumsum(read_weights[::-1], axis=0)[::-1]
    step_variance = endpoint_rms_m**2 / (len(read_weights) - 1)
    # The initial drift value is fixed to zero, so there is no increment before
    # the first read.  Increment j contributes the tail beginning at read j.
    return step_variance * reverse_tails[1:].T @ reverse_tails[1:]


def main() -> int:
    baseline_dir = PHYSICAL_ROOT / "baseline"
    drift_dir = PHYSICAL_ROOT / "time_drift"
    with (baseline_dir / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    with (drift_dir / "scan_metadata.toml").open("rb") as stream:
        drift_metadata = tomllib.load(stream)

    levels = np.asarray(metadata["k2_levels"], dtype=float)
    delta_k2 = levels * float(metadata["k2_step_m3"])
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in rows(baseline_dir / "bump_points.csv")
        ]
    )
    bump_amplitude = float(np.max(np.abs(bump_commands)))
    normalization = float(np.ptp(delta_k2)) * 2.0 * bump_amplitude

    baseline = np.asarray(np.load(baseline_dir / "bpm_orbits.npy", mmap_mode="r"))
    drift_scan = np.asarray(np.load(drift_dir / "bpm_orbits.npy", mmap_mode="r"))
    require(baseline.shape[:2] == (76, 4), "Expected the full 76 x 4 exact benchmark")
    require(drift_scan.shape == baseline.shape, "Paired drift tensor shape mismatch")
    require(np.all(np.isfinite(baseline)), "Non-finite baseline orbit")
    require(np.all(np.isfinite(drift_scan)), "Non-finite drift orbit")

    saved_templates = np.load(OUTPUT_DIR / "source_templates.npy")
    regenerated_templates = inverse.source_templates(MODEL_DIR, 0.272)
    np.testing.assert_allclose(saved_templates, regenerated_templates, rtol=0.0, atol=0.0)
    saved_design = np.load(OUTPUT_DIR / "center_design.npy")
    regenerated_design = inverse.center_design(regenerated_templates)
    np.testing.assert_allclose(saved_design, regenerated_design, rtol=0.0, atol=0.0)
    require(saved_templates.shape[0] == 76, "Source templates do not cover all targets")
    require(saved_design.shape == (76, 2 * saved_templates.shape[1], 2), "Bad design shape")

    response = inverse.recover_drift_response(
        baseline,
        drift_scan,
        float(drift_metadata["drift_halfwidth_m"]),
    )
    states = inverse.signed_state_indices(bump_commands, delta_k2)
    short_repeats = 3
    endpoint_rms = 1.0e-5
    closed_form = inverse.random_walk_center_covariances(
        saved_design,
        response,
        states,
        short_repeats,
        float(np.ptp(delta_k2)),
        2.0 * bump_amplitude,
        endpoint_rms,
    )
    for target, realization in ((0, 0), (13, 2), (43, 1), (75, 3)):
        explicit = explicit_random_walk_covariance(
            saved_design[target],
            response[target, realization],
            states,
            short_repeats,
            normalization,
            endpoint_rms,
        )
        np.testing.assert_allclose(
            closed_form[target, realization], explicit, rtol=2.0e-13, atol=1.0e-30
        )

    selected_repeats = 4096
    bpm_noise = 5.0e-6
    saved_white = np.load(OUTPUT_DIR / "white_center_covariances.npy")
    analytic_white = inverse.white_center_covariances(
        saved_design,
        bpm_noise,
        selected_repeats,
        float(np.ptp(delta_k2)),
        2.0 * bump_amplitude,
    )
    np.testing.assert_allclose(
        saved_white,
        np.broadcast_to(analytic_white[:, None], saved_white.shape),
        rtol=2.0e-14,
        atol=1.0e-30,
    )
    saved_drift = np.load(OUTPUT_DIR / "drift_center_covariances.npy")
    analytic_drift = inverse.random_walk_center_covariances(
        saved_design,
        response,
        states,
        selected_repeats,
        float(np.ptp(delta_k2)),
        2.0 * bump_amplitude,
        endpoint_rms,
    )
    np.testing.assert_allclose(saved_drift, analytic_drift, rtol=2.0e-14, atol=1.0e-30)

    samples = np.load(OUTPUT_DIR / "center_error_samples.npz")
    summary = {row["case"]: row for row in rows(OUTPUT_DIR / "summary.csv")}
    expected_shapes = {
        "clean": (76, 4, 2),
        "bpm_white_noise": (512, 76, 4, 2),
        "random_walk_drift": (512, 76, 4, 2),
        "combined": (512, 76, 4, 2),
    }
    for case, shape in expected_shapes.items():
        errors = samples[case]
        require(errors.shape == shape, f"Unexpected {case} sample shape: {errors.shape}")
        require(np.all(np.isfinite(errors)), f"Non-finite {case} errors")
        recalculated = inverse.summarize(errors)
        for key, value in recalculated.items():
            np.testing.assert_allclose(float(summary[case][key]), value, rtol=2.0e-14)

    with (OUTPUT_DIR / "result_metadata.json").open(encoding="utf-8") as stream:
        result_metadata = json.load(stream)
    require(result_metadata["target_count"] == 76, "Metadata target count mismatch")
    require(result_metadata["exact_forward_state_count"] == 9120, "Exact state count mismatch")
    require(result_metadata["repeats_per_signed_state"] == selected_repeats, "Repeat mismatch")
    require(result_metadata["acquisitions_per_target_scan"] == 32768, "Read count mismatch")
    require(result_metadata["acceptance_gate_passed"] is True, "Saved acceptance gate failed")
    require(float(summary["combined"]["rmse_2d_um"]) < 50.0, "Combined RMSE gate failed")
    require(float(summary["combined"]["p99_2d_um"]) < 50.0, "Combined P99 gate failed")
    combined_targets = [
        row for row in rows(OUTPUT_DIR / "per_target_summary.csv") if row["case"] == "combined"
    ]
    require(len(combined_targets) == 76, "Missing combined target rows")
    require(
        max(float(row["rmse_2d_um"]) for row in combined_targets) < 50.0,
        "Worst-target combined RMSE gate failed",
    )

    # A Gaussian draw check is statistical, so use a deliberately loose bound;
    # deterministic formula/file checks above carry the strict validation.
    for case in ("bpm_white_noise", "random_walk_drift", "combined"):
        measured = float(summary[case]["rmse_2d_um"])
        expected = float(summary[case]["analytic_expected_rmse_2d_um"])
        require(abs(measured / expected - 1.0) < 0.03, f"{case} Monte Carlo mismatch")

    print("PASS: saved summaries and metadata reproduce exactly")
    print("PASS: white-noise covariance matches the analytic least-squares covariance")
    print("PASS: random-walk closed form matches explicit read-by-read propagation")
    print("PASS: combined P99 and all 76 target-level RMSEs satisfy the 50 um gate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
