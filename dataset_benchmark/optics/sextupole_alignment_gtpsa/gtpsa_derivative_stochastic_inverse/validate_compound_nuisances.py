#!/usr/bin/env python3
"""Validate the current-amplitude compound nuisance reanalysis."""

from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

import numpy as np

import analyze_compound_nuisances as compound
import analyze_stochastic_inverse as base


HERE = Path(__file__).resolve().parent
PHYSICAL_ROOT = HERE / "results" / "exact_k5_b3"
OUTPUT = HERE / "results" / "compound_nuisance_analysis"
TIME_SERIES_OUTPUT = HERE / "results" / "time_series_analysis"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def metadata(path: Path) -> dict[str, object]:
    with (path / "scan_metadata.toml").open("rb") as stream:
        return tomllib.load(stream)


def main() -> int:
    physical_cases = (
        compound.BASELINE,
        "time_drift",
        *compound.PHYSICAL_COMPONENTS,
        compound.COMPOUND,
        compound.COMPOUND_DRIFT,
    )
    baseline_dir = PHYSICAL_ROOT / compound.BASELINE
    baseline_meta = metadata(baseline_dir)
    baseline_truth = np.load(baseline_dir / "target_truth.npy")
    baseline_sextupoles = np.load(baseline_dir / "latent_sextupole_offsets.npy")
    baseline_targets = (baseline_dir / "target_names.txt").read_text(
        encoding="utf-8"
    )
    baseline_bpms = (baseline_dir / "bpm_names.txt").read_text(encoding="utf-8")
    for case in physical_cases:
        source = PHYSICAL_ROOT / case
        case_meta = metadata(source)
        require(case_meta["lattice"] == baseline_meta["lattice"], f"{case} lattice differs")
        require(case_meta["target_count"] == 76, f"{case} does not cover 76 targets")
        require(
            case_meta["realization_count_per_target"] == 4,
            f"{case} does not contain four realizations",
        )
        require(
            case_meta["bump_amplitude_m"] == 0.0015,
            f"{case} is not the current 1.5 mm protocol",
        )
        require(case_meta["k2_step_m3"] == 0.05, f"{case} has the old K2 span")
        require(
            case_meta["k2_levels"] == [-2.0, 0.0, 2.0],
            f"{case} does not reach the current +/-0.10 m^-3 extrema",
        )
        require(
            (source / "target_names.txt").read_text(encoding="utf-8") == baseline_targets,
            f"{case} target inventory differs",
        )
        require(
            (source / "bpm_names.txt").read_text(encoding="utf-8") == baseline_bpms,
            f"{case} BPM inventory differs",
        )
        np.testing.assert_array_equal(np.load(source / "target_truth.npy"), baseline_truth)
        np.testing.assert_array_equal(
            np.load(source / "latent_sextupole_offsets.npy"), baseline_sextupoles
        )

    component_arrays = {
        "corrector_gain": "latent_corrector_gain_errors.npy",
        "k2_calibration": "latent_k2_gain_errors.npy",
        "quadrupole_strength": "latent_quadrupole_relative_errors.npy",
        "quadrupole_roll": "latent_quadrupole_rolls.npy",
    }
    combined_dir = PHYSICAL_ROOT / compound.COMPOUND
    combined_drift_dir = PHYSICAL_ROOT / compound.COMPOUND_DRIFT
    for case, filename in component_arrays.items():
        individual = np.load(PHYSICAL_ROOT / case / filename)
        combined = np.load(combined_dir / filename)
        combined_drift = np.load(combined_drift_dir / filename)
        require(np.any(individual != 0.0), f"{case} nuisance is unexpectedly zero")
        np.testing.assert_array_equal(combined, individual)
        np.testing.assert_array_equal(combined_drift, individual)
    for filename in (
        "target_truth.npy",
        "latent_sextupole_offsets.npy",
        "latent_corrector_gain_errors.npy",
        "latent_k2_gain_errors.npy",
        "latent_quadrupole_relative_errors.npy",
        "latent_quadrupole_rolls.npy",
        "latent_quadrupole_offsets.npy",
    ):
        np.testing.assert_array_equal(
            np.load(combined_dir / filename), np.load(combined_drift_dir / filename)
        )
    require(
        not np.any(np.load(combined_dir / "latent_quadrupole_offsets.npy")),
        "Combined case accidentally contains quadrupole misalignment",
    )
    np.testing.assert_array_equal(
        np.load(combined_drift_dir / "latent_drift_directions.npy"),
        np.load(PHYSICAL_ROOT / "time_drift" / "latent_drift_directions.npy"),
    )

    baseline_errors = np.load(OUTPUT / "baseline_deterministic_errors.npy")
    compound_errors = np.load(OUTPUT / "compound_deterministic_errors.npy")
    linear_errors = np.load(OUTPUT / "linear_sum_errors.npy")
    interaction = np.load(OUTPUT / "nonlinear_interaction_errors.npy")
    require(baseline_errors.shape == (76, 4, 2), "Baseline error shape mismatch")
    require(compound_errors.shape == baseline_errors.shape, "Compound error shape mismatch")
    np.testing.assert_allclose(
        compound_errors - linear_errors, interaction, rtol=0.0, atol=2.0e-22
    )

    baseline_covariance = np.load(OUTPUT / "baseline_total_covariances.npy")
    maintained_covariance = (
        np.load(TIME_SERIES_OUTPUT / "white_center_covariances.npy")
        + np.load(TIME_SERIES_OUTPUT / "filtered_drift_covariances.npy")
    )
    np.testing.assert_allclose(
        baseline_covariance, maintained_covariance, rtol=3.0e-13, atol=2.0e-28
    )
    for filename in (
        "baseline_total_covariances.npy",
        "baseline_white_covariances.npy",
        "baseline_drift_covariances.npy",
        "compound_total_covariances.npy",
        "compound_white_covariances.npy",
        "compound_drift_covariances.npy",
    ):
        covariance = np.load(OUTPUT / filename)
        require(covariance.shape == (76, 4, 2, 2), f"{filename} shape mismatch")
        require(np.all(np.isfinite(covariance)), f"{filename} contains non-finite values")
        require(
            np.min(np.linalg.eigvalsh(covariance)) >= -1.0e-24,
            f"{filename} is not positive semidefinite",
        )

    samples = np.load(OUTPUT / "center_error_samples.npz")
    compound_samples = samples["compound_with_white_and_drift"]
    require(compound_samples.shape == (512, 76, 4, 2), "Compound sample shape mismatch")
    recalculated = base.summarize(compound_samples)
    summary = {row["case"]: row for row in rows(OUTPUT / "summary.csv")}
    saved = summary["compound_with_white_and_drift"]
    for key, value in recalculated.items():
        np.testing.assert_allclose(float(saved[key]), value, rtol=2.0e-14)

    with (OUTPUT / "result_metadata.json").open(encoding="utf-8") as stream:
        result = json.load(stream)
    require(result["excluded_nuisance"] == "quadrupole_misalignment", "Wrong exclusion")
    require(result["target_count"] == 76, "Saved result target count mismatch")
    require(result["realizations_per_target"] == 4, "Saved realization count mismatch")
    np.testing.assert_allclose(
        result["compound_combined_rmse_um"], recalculated["rmse_2d_um"], rtol=2.0e-14
    )
    np.testing.assert_allclose(
        result["compound_combined_p99_um"], recalculated["p99_2d_um"], rtol=2.0e-14
    )
    target_rows = [
        row
        for row in rows(OUTPUT / "per_target_summary.csv")
        if row["case"] == "compound_with_white_and_drift"
    ]
    baseline_target_rows = [
        row
        for row in rows(OUTPUT / "per_target_summary.csv")
        if row["case"] == "baseline_with_white_and_drift"
    ]
    require(len(target_rows) == 76, "Compound per-target summary is incomplete")
    require(len(baseline_target_rows) == 76, "Baseline per-target summary is incomplete")
    worst_target = max(float(row["rmse_2d_um"]) for row in target_rows)
    baseline_worst_target = max(
        float(row["rmse_2d_um"]) for row in baseline_target_rows
    )
    np.testing.assert_allclose(
        result["baseline_worst_target_rmse_um"],
        baseline_worst_target,
        rtol=2.0e-14,
    )
    targets_above_required = [
        row["target"] for row in target_rows if float(row["rmse_2d_um"]) >= 50.0
    ]
    np.testing.assert_allclose(
        result["compound_worst_target_rmse_um"], worst_target, rtol=2.0e-14
    )
    require(
        result["targets_at_or_above_required_rmse"] == targets_above_required,
        "Saved over-threshold target list mismatch",
    )
    expected_hard_gate = (
        recalculated["rmse_2d_um"] < 50.0
        and recalculated["p99_2d_um"] < 50.0
        and worst_target < 50.0
    )
    require(result["hard_gate_passed"] is expected_hard_gate, "Hard gate flag mismatch")
    require(
        result["preferred_gate_passed"] is (recalculated["rmse_2d_um"] < 30.0),
        "Preferred gate flag mismatch",
    )

    with (OUTPUT / "interaction_summary.json").open(encoding="utf-8") as stream:
        interaction_saved = json.load(stream)
    interaction_rms = float(
        np.sqrt(np.mean(np.sum(interaction * interaction, axis=-1))) * 1e6
    )
    np.testing.assert_allclose(
        interaction_saved["nonlinear_interaction_rms_um"], interaction_rms, rtol=2.0e-14
    )
    interaction_norm = np.linalg.norm(interaction, axis=-1) * 1e6
    interaction_peak = np.unravel_index(
        int(np.argmax(interaction_norm)), interaction_norm.shape
    )
    np.testing.assert_allclose(
        interaction_saved["nonlinear_interaction_max_um"],
        interaction_norm[interaction_peak],
        rtol=2.0e-14,
    )
    require(
        interaction_saved["nonlinear_interaction_max_target"]
        == baseline_targets.splitlines()[interaction_peak[0]],
        "Saved peak-interaction target mismatch",
    )

    covariance_components = {
        "baseline_stochastic_component_rmse_um": "baseline_total_covariances.npy",
        "compound_white_component_rmse_um": "compound_white_covariances.npy",
        "compound_filtered_drift_component_rmse_um": "compound_drift_covariances.npy",
        "compound_total_stochastic_component_rmse_um": "compound_total_covariances.npy",
    }
    for key, filename in covariance_components.items():
        covariance = np.load(OUTPUT / filename)
        value = float(
            np.sqrt(np.mean(np.trace(covariance, axis1=-2, axis2=-1))) * 1e6
        )
        np.testing.assert_allclose(result[key], value, rtol=2.0e-14)

    print("PASS: every physical scan uses the current 1.5 mm / 0.10 m^-3 protocol")
    print("PASS: individual and compound nuisance tensors use exactly paired component seeds")
    print("PASS: quadrupole misalignment is zero in the combined case")
    print("PASS: compound = paired linear sum + independently saved interaction residual")
    print("PASS: baseline stochastic covariance reproduces the maintained time-series result")
    print(
        "PASS: compound statistics and gate flags reproduce; "
        f"RMSE={recalculated['rmse_2d_um']:.3f} um, P99={recalculated['p99_2d_um']:.3f} um"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
