#!/usr/bin/env python3
"""Run the paired finite-BPM center inverse for one nuisance at a time."""

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
FINITE_BPM = STUDY_ROOT / "finite_bpm_inversion"
sys.path.insert(0, str(FINITE_BPM))

from analyze_command_space_finite_bpm import fit_center, k2_slope  # noqa: E402
from analyze_local_orbit_predictors import build_two_sided_maps, read_rows  # noqa: E402


PHYSICAL_CASES = (
    "baseline",
    "corrector_gain",
    "k2_calibration",
    "quadrupole_strength",
    "quadrupole_roll",
    "quadrupole_misalignment",
    "time_drift",
)
CASE_ORDER = (
    "baseline",
    "bpm_gain",
    "corrector_gain",
    "k2_calibration",
    "quadrupole_strength",
    "quadrupole_roll",
    "quadrupole_misalignment",
    "time_drift",
    "bpm_noise",
)
CASE_LABELS = {
    "baseline": "Reference: no added nuisance",
    "bpm_gain": "BPM gain",
    "corrector_gain": "Corrector gain",
    "k2_calibration": "K2 calibration gain",
    "quadrupole_strength": "Quadrupole strength",
    "quadrupole_roll": "Quadrupole roll",
    "quadrupole_misalignment": "Quadrupole misalignment",
    "time_drift": "Time drift",
    "bpm_noise": "BPM noise",
}
CASE_MAGNITUDES = {
    "baseline": "none",
    "bpm_gain": "1% RMS per BPM/plane, fixed in scan",
    "corrector_gain": "1% RMS per corrector, fixed in scan",
    "k2_calibration": "1% RMS intervention gain per target scan",
    "quadrupole_strength": "independent uniform +/-1%",
    "quadrupole_roll": "1 mrad RMS",
    "quadrupole_misalignment": "100 um RMS in x and y",
    "time_drift": "+/-5 um target-local linear drift across scan",
    "bpm_noise": "5 um RMS per BPM/plane/state",
}


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
        "p99_2d_um": float(np.percentile(radial_um, 99)),
        "max_2d_um": float(np.max(radial_um)),
    }


def load_model(source: Path, model_dir: Path, knobs_path: Path) -> dict[str, object]:
    bpm_names = (source / "bpm_names.txt").read_text(encoding="utf-8").splitlines()
    target_names = (source / "target_names.txt").read_text(encoding="utf-8").splitlines()
    bpm_rows = read_rows(model_dir / "bpm_locations.csv")
    all_target_rows = read_rows(model_dir / "target_locations.csv")
    control_rows = read_rows(model_dir / "control_inventory.csv")
    if bpm_names != [row["bpm"] for row in bpm_rows]:
        raise ValueError("Scan and model BPM inventories differ")
    model_target_lookup = {
        row["target"]: index for index, row in enumerate(all_target_rows)
    }
    if any(name not in model_target_lookup for name in target_names):
        raise ValueError("At least one scan target is absent from the model inventory")
    target_model_indices = np.asarray(
        [model_target_lookup[name] for name in target_names], dtype=int
    )
    target_rows = [all_target_rows[index] for index in target_model_indices]

    bpm_response = np.load(model_dir / "bpm_control_response.npy")
    target_response_flat = np.load(model_dir / "target_control_response.npy")
    bpm_maps = np.load(model_dir / "bpm_cumulative_maps.npy")
    target_maps = np.load(model_dir / "target_cumulative_maps.npy")[target_model_indices]
    one_turn = np.load(model_dir / "one_turn_map.npy")
    nt, nd, nc = len(target_names), len(bpm_names), len(control_rows)
    target_response = target_response_flat.reshape(len(all_target_rows), 2, nc)[
        target_model_indices
    ]

    bump_rows = read_rows(source / "bump_points.csv")
    bump_commands = np.asarray(
        [
            (float(row["bump_x_command_m"]), float(row["bump_y_command_m"]))
            for row in bump_rows
        ]
    )
    zero_candidates = np.flatnonzero(np.all(bump_commands == 0.0, axis=1))
    if zero_candidates.size != 1:
        raise ValueError("Expected exactly one zero bump")
    zero_bump = int(zero_candidates[0])

    control_lookup = {
        (row["corrector"], row["field"]): index
        for index, row in enumerate(control_rows)
    }
    target_lookup = {name: index for index, name in enumerate(target_names)}
    knob_x = np.zeros((nt, nc))
    knob_y = np.zeros((nt, nc))
    for row in read_rows(knobs_path):
        if row["target_sextupole"] not in target_lookup:
            continue
        target = target_lookup[row["target_sextupole"]]
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
        np.asarray([int(row["line_index"]) for row in target_rows]),
    )
    return {
        "bpm_names": bpm_names,
        "target_names": target_names,
        "bump_commands": bump_commands,
        "zero_bump": zero_bump,
        "model_bpm": model_bpm,
        "model_target": model_target,
        "two_sided_maps": two_sided_maps,
        "neighbor_rows": neighbor_rows,
    }


def predict_two_sided(
    measured_bpm: np.ndarray,
    nominal_k2: int,
    model: dict[str, object],
) -> np.ndarray:
    """Leakage-safe local-orbit prediction using BPM data and nominal model only."""
    zero_bump = int(model["zero_bump"])
    model_bpm = np.asarray(model["model_bpm"])
    model_target = np.asarray(model["model_target"])
    two_sided_maps = np.asarray(model["two_sided_maps"])
    neighbor_rows = model["neighbor_rows"]
    observed = np.asarray(measured_bpm[:, :, :, nominal_k2, :, :], dtype=float)
    observed = observed - observed[:, :, zero_bump : zero_bump + 1, :, :]
    observed_flat = observed.reshape(*observed.shape[:3], -1)
    residual = observed_flat - model_bpm[:, None, :, :]
    nt, nr, nb = observed.shape[:3]
    prediction = np.broadcast_to(model_target[:, None, :, :], (nt, nr, nb, 2)).copy()
    for target in range(nt):
        upstream = int(neighbor_rows[target]["upstream_bpm_index"]) - 1
        downstream = int(neighbor_rows[target]["downstream_bpm_index"]) - 1
        channels = np.asarray(
            [2 * upstream, 2 * upstream + 1, 2 * downstream, 2 * downstream + 1]
        )
        nearby_residual = np.take(residual[target], channels, axis=-1)
        prediction[target] += nearby_residual @ two_sided_maps[target].T
    if not np.all(np.isfinite(prediction)):
        raise ValueError("Non-finite two-sided local-orbit prediction")
    return prediction


def measured_case(
    case_name: str,
    physical_root: Path,
    gain_rms: float,
    noise_rms_m: float,
    measurement_seed: int,
    output: Path,
) -> tuple[Path, np.ndarray]:
    source = physical_root / ("baseline" if case_name in {"bpm_gain", "bpm_noise"} else case_name)
    bpm = np.array(np.load(source / "bpm_orbits.npy", mmap_mode="r"), dtype=float, copy=True)
    if case_name == "bpm_gain":
        rng = np.random.default_rng(measurement_seed)
        gains = gain_rms * rng.standard_normal((bpm.shape[0], bpm.shape[1], bpm.shape[4], 2))
        bpm *= 1.0 + gains[:, :, None, None, :, :]
        np.save(output / "latent_bpm_gain_errors.npy", gains)
    elif case_name == "bpm_noise":
        rng = np.random.default_rng(measurement_seed + 1)
        bpm += noise_rms_m * rng.standard_normal(bpm.shape)
    return source, bpm


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--physical-root", type=Path,
        default=HERE / "results" / "physical_scans",
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
    parser.add_argument("--bpm-gain-rms", type=float, default=0.01)
    parser.add_argument("--bpm-noise-rms-m", type=float, default=5.0e-6)
    parser.add_argument("--measurement-seed", type=int, default=20260829)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results" / "analysis")
    parser.add_argument("--readme-path", type=Path, default=HERE / "README.md")
    args = parser.parse_args()
    physical_root = args.physical_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    baseline_source = physical_root / "baseline"
    with (baseline_source / "scan_metadata.toml").open("rb") as stream:
        baseline_metadata = tomllib.load(stream)
    model = load_model(baseline_source, args.model_dir.resolve(), args.knobs.resolve())
    target_names = list(model["target_names"])
    zero_bump = int(model["zero_bump"])
    summary_rows: list[dict[str, object]] = []
    realization_rows: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    errors_by_case: dict[str, np.ndarray] = {}
    estimates_by_case: dict[str, np.ndarray] = {}

    for case_index, case_name in enumerate(CASE_ORDER):
        source, measured_bpm = measured_case(
            case_name,
            physical_root,
            args.bpm_gain_rms,
            args.bpm_noise_rms_m,
            args.measurement_seed,
            output,
        )
        with (source / "scan_metadata.toml").open("rb") as stream:
            metadata = tomllib.load(stream)
        levels = np.asarray(metadata["k2_levels"], dtype=float)
        nominal_candidates = np.flatnonzero(levels == 0.0)
        if nominal_candidates.size != 1:
            raise ValueError("Expected exactly one nominal K2 level")
        nominal_k2 = int(nominal_candidates[0])
        commanded_delta_k2 = levels * float(metadata["k2_step_m3"])

        predicted_local = predict_two_sided(measured_bpm, nominal_k2, model)
        nt, nr, nb, nk, nd, planes = measured_bpm.shape
        if (nt, nb, nk, planes) != (len(target_names), 5, 3, 2):
            raise ValueError(f"Unexpected measured tensor shape: {measured_bpm.shape}")
        estimates = np.zeros((nt, nr, 2), dtype=float)
        for target in range(nt):
            for realization in range(nr):
                slopes = k2_slope(
                    measured_bpm[target, realization].reshape(nb, nk, 2 * nd),
                    commanded_delta_k2,
                )
                estimates[target, realization] = fit_center(
                    slopes, predicted_local[target, realization]
                )
        np.save(output / f"{case_name}_relative_center_estimates.npy", estimates)
        np.save(output / f"{case_name}_predicted_local_orbits.npy", predicted_local)

        # Evaluation-only truth begins here.
        target_truth = np.load(source / "target_truth.npy")
        target_orbits = np.load(source / "target_orbits.npy", mmap_mode="r")
        zero_orbit = np.asarray(target_orbits[:, :, zero_bump, nominal_k2, :], dtype=float)
        relative_truth = target_truth - zero_orbit
        exact_relative_local = (
            np.asarray(target_orbits[:, :, :, nominal_k2, :], dtype=float)
            - zero_orbit[:, :, None, :]
        )
        center_errors = estimates - relative_truth
        errors_by_case[case_name] = center_errors
        estimates_by_case[case_name] = estimates
        nonzero_bumps = np.arange(nb) != zero_bump
        local_errors = predicted_local - exact_relative_local
        local_rmse_um = float(
            np.sqrt(
                np.mean(
                    np.sum(local_errors[:, :, nonzero_bumps, :] ** 2, axis=-1)
                )
            )
            * 1e6
        )
        overall = summarize(center_errors)
        baseline_error = errors_by_case.get("baseline")
        incremental_rms_um = 0.0 if case_name == "baseline" else float(
            np.sqrt(np.mean(np.sum((center_errors - baseline_error) ** 2, axis=-1)))
            * 1e6
        )
        reference_rmse_um = (
            overall["rmse_2d_um"]
            if case_name == "baseline"
            else float(summary_rows[0]["rmse_2d_um"])
        )
        summary_rows.append(
            {
                "case": case_name,
                "label": CASE_LABELS[case_name],
                "nuisance_magnitude": CASE_MAGNITUDES[case_name],
                "fit_count": nt * nr,
                "local_orbit_rmse_2d_um": local_rmse_um,
                **overall,
                "delta_rmse_2d_vs_baseline_um": overall["rmse_2d_um"]
                - reference_rmse_um,
                "incremental_error_vector_rms_um": incremental_rms_um,
            }
        )
        radial_um = np.linalg.norm(center_errors, axis=-1) * 1e6
        for target, name in enumerate(target_names):
            target_summary = summarize(center_errors[target])
            per_target_rows.append(
                {
                    "case": case_name,
                    "target": name,
                    "target_index": target + 1,
                    **target_summary,
                }
            )
            for realization in range(nr):
                realization_rows.append(
                    {
                        "case": case_name,
                        "target": name,
                        "target_index": target + 1,
                        "realization": realization + 1,
                        "relative_truth_x_um": relative_truth[target, realization, 0] * 1e6,
                        "relative_truth_y_um": relative_truth[target, realization, 1] * 1e6,
                        "estimate_x_um": estimates[target, realization, 0] * 1e6,
                        "estimate_y_um": estimates[target, realization, 1] * 1e6,
                        "error_x_um": center_errors[target, realization, 0] * 1e6,
                        "error_y_um": center_errors[target, realization, 1] * 1e6,
                        "error_2d_um": radial_um[target, realization],
                    }
                )
        print(
            f"{case_index + 1}/{len(CASE_ORDER)} {case_name}: "
            f"center RMSE {overall['rmse_2d_um']:.3f} um; "
            f"local RMSE {local_rmse_um:.3f} um"
        )

    write_rows(output / "summary.csv", summary_rows)
    write_rows(output / "per_realization_fits.csv", realization_rows)
    write_rows(output / "per_target_summary.csv", per_target_rows)

    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    x = np.arange(len(summary_rows))
    values = np.asarray([float(row["rmse_2d_um"]) for row in summary_rows])
    ax.bar(x, values, color=["#6c757d"] + ["#4472c4"] * (len(x) - 1))
    ax.axhline(values[0], color="#c00000", ls="--", lw=1.3, label="reference")
    ax.set_xticks(x, [str(row["label"]).replace("Reference: ", "") for row in summary_rows], rotation=35, ha="right")
    ax.set_ylabel("Beam-relative center 2D RMSE [micrometers]")
    ax.set_title("One-at-a-time real-machine nuisance ablation")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "center_rmse_by_nuisance.png", dpi=180)
    plt.close(fig)

    table_lines = "\n".join(
        f"| {row['label']} | {row['nuisance_magnitude']} | "
        f"{float(row['local_orbit_rmse_2d_um']):.3f} | "
        f"{float(row['rmse_2d_um']):.3f} | "
        f"{float(row['delta_rmse_2d_vs_baseline_um']):+.3f} | "
        f"{float(row['incremental_error_vector_rms_um']):.3f} | "
        f"{float(row['median_2d_um']):.3f} | {float(row['p90_2d_um']):.3f} | "
        f"{float(row['max_2d_um']):.3f} |"
        for row in summary_rows
    )
    quadrupole_misalignment_rows = [
        row for row in realization_rows if row["case"] == "quadrupole_misalignment"
    ]
    quadrupole_misalignment_truth = np.asarray(
        [
            (float(row["relative_truth_x_um"]), float(row["relative_truth_y_um"]))
            for row in quadrupole_misalignment_rows
        ]
    )
    outside_fit_box = int(
        np.sum(np.any(np.abs(quadrupole_misalignment_truth) > 1500.0, axis=1))
    )
    outside_bump_radius = int(
        np.sum(np.linalg.norm(quadrupole_misalignment_truth, axis=1) > 500.0)
    )
    physical_seconds = sum(
        float(tomllib.load((physical_root / case_name / "scan_metadata.toml").open("rb"))["calculation_wall_seconds"])
        for case_name in PHYSICAL_CASES
    )
    report = f"""# Real-machine nuisance ablation for finite-BPM sextupole alignment

This paired study starts from a clean reference that retains the central hard
condition: the target alignment is unknown and the other 75 sextupoles carry
independent `300 micrometer` RMS x/y offsets. Each non-reference row adds
exactly one nuisance. Physical lattice/actuator cases were regenerated with
the validated latest repaired SciBmad lattice; BPM gain and BPM noise were
applied only to simulated readbacks. The inverse always receives nominal bump
and K2 commands plus BPM readings, never latent nuisance values or exact
target-local orbit.

The magnitudes below are representative sensitivity-test settings, not measured
CESR calibration distributions. There are {baseline_metadata['target_count']}
targets and {baseline_metadata['realization_count_per_target']} paired latent
machines per target, or {baseline_metadata['target_count'] * baseline_metadata['realization_count_per_target']}
fits per row.

## Result table

| added machine error | test magnitude | local-orbit 2D RMSE [um] | center 2D RMSE [um] | aggregate change [um] | paired error-vector increment RMS [um] | median [um] | P90 [um] | max [um] |
|---|---|---:|---:|---:|---:|---:|---:|---:|
{table_lines}

The paired error-vector increment compares each nuisance realization with its
matched reference and is the most direct one-at-a-time impact measure. It must
not be confused with the difference between aggregate RMSE values: a nuisance
can partially cancel the reference fit error and lower aggregate RMSE while
still changing individual estimates.

## Interpretation checks

- A fixed multiplicative K2 calibration gain has negligible effect here because
  the maintained fit normalizes every BPM K2-slope channel and fits a free
  propagation matrix. This does not cover K2 hysteresis, polarity asymmetry, or
  point-to-point calibration drift.
- Quadrupole strength and roll slightly lower aggregate RMSE in this finite
  paired sample, but their nonzero paired increment quantifies the actual
  estimate change; they are not beneficial corrections.
- The quadrupole-misalignment row is an **uncorrected-orbit stress test**, not a
  transfer-matrix-only result. Its 100-micrometer RMS offsets push
  {outside_bump_radius}/{len(quadrupole_misalignment_rows)} beam-relative truths
  outside the 0.5-mm bump radius and {outside_fit_box}/{len(quadrupole_misalignment_rows)}
  outside the current +/-1.5-mm-per-plane fit box. A separate orbit-corrected
  misalignment study is needed to isolate residual matrix mismatch.
- Time drift and BPM noise are intentionally passed to the unchanged three-point
  slope estimator without drift regression, repeated-read averaging, or
  covariance weighting. Their large errors diagnose protocol sensitivity, not
  the best achievable calibrated-machine performance.

## Nuisance definitions

- **BPM gain:** independent multiplicative x/y calibration error, fixed for a
  BPM throughout one scan.
- **Corrector gain:** independent multiplicative error on each corrector's bump
  increment, fixed throughout one scan; the base corrector setting is unchanged.
- **K2 calibration:** one multiplicative error on the target's commanded K2
  intervention, shared by every bump and K2 point in the scan.
- **Quadrupole strength:** independent physical Kn1 errors uniformly bounded by
  `+/-1%`.
- **Quadrupole roll:** independent roll added coherently to every tracking slice
  belonging to the same physical quadrupole.
- **Quadrupole misalignment:** independent x/y displacement added coherently to
  every tracking slice belonging to the same physical quadrupole.
- **Time drift:** a random transverse direction with a target-local command that
  changes linearly from `-5` to `+5 micrometers` in acquisition order. It is zero
  at the zero-bump, nominal-K2 reference state and is propagated physically by
  the same local-bump correctors.
- **BPM noise:** independent Gaussian readout noise for every BPM plane and scan
  state.

The reference intentionally omits quadrupole strength error, unlike the earlier
maintained all-76 result. This makes every row a one-at-a-time ablation; its
absolute reference RMSE therefore need not equal the earlier `5.864 micrometer`
mixed-nuisance result.

## Method and provenance

- lattice: `{baseline_metadata['lattice']}`
- exact physical SciBmad states per physical case: `{baseline_metadata['total_state_count']}`
- summed physical-generation wall time: `{physical_seconds:.1f} s`
- local orbit: nominal-model command prediction corrected by the nearest
  upstream/downstream BPM pair
- center inverse: all-111-BPM symmetric three-point K2 slope and the maintained
  shared thin-sextupole source fit
- exact target orbit and target alignment: evaluation only

Run from `CESR Project/`:

```powershell
julia --project=. sextupole_misalignment/real_machine_nuisance_ablation/generate_physical_nuisance_scans.jl

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/real_machine_nuisance_ablation/analyze_nuisance_ablation.py'

wsl.exe -d Ubuntu-Bmad -- /home/joeyfurina/miniforge3/envs/bmad/bin/python `
  '/mnt/d/Ring_Design_Development/CESR Project/sextupole_misalignment/real_machine_nuisance_ablation/validate_nuisance_ablation.py'
```

## Limitations

This is still synthetic and noise magnitudes are assumed. Each error is tested
alone, so the table is a sensitivity decomposition rather than a prediction of
the fully combined real-machine error. A later combined study should use
measured calibration priors, correlated girder/family errors, interleaved scan
timing, and covariance-aware or joint-nuisance inference.
"""
    readme_path = args.readme_path.resolve()
    readme_path.parent.mkdir(parents=True, exist_ok=True)
    readme_path.write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
