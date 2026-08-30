#!/usr/bin/env python3
"""Compare uncorrected, finite-difference, and GTPSA/noisy correction runs."""

from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SEQUENTIAL = HERE.parent / "sequential_joint_inverse"
RESULTS = SEQUENTIAL / "results"
SCANS = RESULTS / "exact_joint_machines"

WITHOUT = "without_quadrupole_misalignment"
UNCORRECTED = "with_quadrupole_misalignment"
FD_CORRECTED = "with_quadrupole_misalignment_corrected"
GTPSA_NOISY_CORRECTED = (
    "with_quadrupole_misalignment_gtpsa_noisy_corrected"
)

MODELS = (
    "physics_gls",
    "shared_target_local_ridge",
    "shared_joint_ridge",
    "shared_joint_random_feature",
)
LEARNED_MODELS = MODELS[1:]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def same_distribution(
    rows: list[dict[str, str]], case: str, model: str
) -> dict[str, str]:
    candidates = [
        row
        for row in rows
        if row["evaluation_case"] == case
        and row["model"] == model
        and (
            row["training_case"] == "fixed_nominal_physics"
            or row["training_case"] == case
        )
    ]
    if len(candidates) != 1:
        raise ValueError(f"Expected one matched row for {case}/{model}")
    return candidates[0]


def metric(row: dict[str, str], name: str) -> float:
    value = float(row[name])
    if not np.isfinite(value):
        raise ValueError(f"Non-finite {name}")
    return value


def best_learned(
    rows: list[dict[str, str]], case: str
) -> dict[str, str]:
    return min(
        (same_distribution(rows, case, model) for model in LEARNED_MODELS),
        key=lambda row: metric(row, "rmse_2d_um"),
    )


def ood_joint(rows: list[dict[str, str]], case: str) -> dict[str, str]:
    candidates = [
        row
        for row in rows
        if row["training_case"] == WITHOUT
        and row["evaluation_case"] == case
        and row["model"] == "shared_joint_ridge"
    ]
    if len(candidates) != 1:
        raise ValueError(f"Expected one zero-trained joint row for {case}")
    return candidates[0]


def rms(array: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(array, dtype=np.float64) ** 2)))


def main() -> int:
    uncorrected_dir = RESULTS / "joint_inverse_analysis"
    fd_dir = RESULTS / "joint_inverse_analysis_corrected"
    noisy_dir = RESULTS / "joint_inverse_analysis_gtpsa_noisy_corrected"
    noisy_dir.mkdir(parents=True, exist_ok=True)

    uncorrected_rows = read_rows(uncorrected_dir / "summary.csv")
    fd_rows = read_rows(fd_dir / "summary.csv")
    noisy_rows = read_rows(noisy_dir / "summary.csv")

    comparison_rows: list[dict[str, object]] = []
    for model in MODELS:
        zero = same_distribution(noisy_rows, WITHOUT, model)
        uncorrected = same_distribution(uncorrected_rows, UNCORRECTED, model)
        fd = same_distribution(fd_rows, FD_CORRECTED, model)
        noisy = same_distribution(noisy_rows, GTPSA_NOISY_CORRECTED, model)
        zero_rmse = metric(zero, "rmse_2d_um")
        uncorrected_rmse = metric(uncorrected, "rmse_2d_um")
        fd_rmse = metric(fd, "rmse_2d_um")
        noisy_rmse = metric(noisy, "rmse_2d_um")
        comparison_rows.append(
            {
                "model": model,
                "zero_offset_rmse_2d_um": zero_rmse,
                "uncorrected_offset_rmse_2d_um": uncorrected_rmse,
                "fd_noiseless_corrected_rmse_2d_um": fd_rmse,
                "gtpsa_noisy_corrected_rmse_2d_um": noisy_rmse,
                "gtpsa_noisy_reduction_percent_vs_uncorrected": (
                    100.0 * (1.0 - noisy_rmse / uncorrected_rmse)
                ),
                "gtpsa_noisy_excess_percent_vs_zero": (
                    100.0 * (noisy_rmse / zero_rmse - 1.0)
                ),
                "gtpsa_noisy_change_percent_vs_fd_noiseless": (
                    100.0 * (noisy_rmse / fd_rmse - 1.0)
                ),
                "fd_noiseless_p99_2d_um": metric(fd, "p99_2d_um"),
                "gtpsa_noisy_p99_2d_um": metric(noisy, "p99_2d_um"),
                "fd_noiseless_worst_target_rmse_2d_um": metric(
                    fd, "worst_target_rmse_2d_um"
                ),
                "gtpsa_noisy_worst_target_rmse_2d_um": metric(
                    noisy, "worst_target_rmse_2d_um"
                ),
            }
        )
    write_rows(noisy_dir / "gtpsa_noisy_protocol_comparison.csv", comparison_rows)

    best_zero = best_learned(noisy_rows, WITHOUT)
    best_uncorrected = best_learned(uncorrected_rows, UNCORRECTED)
    best_fd = best_learned(fd_rows, FD_CORRECTED)
    best_noisy = best_learned(noisy_rows, GTPSA_NOISY_CORRECTED)
    zero_rmse = metric(best_zero, "rmse_2d_um")
    uncorrected_rmse = metric(best_uncorrected, "rmse_2d_um")
    fd_rmse = metric(best_fd, "rmse_2d_um")
    noisy_rmse = metric(best_noisy, "rmse_2d_um")

    fd_scan = SCANS / FD_CORRECTED
    noisy_scan = SCANS / GTPSA_NOISY_CORRECTED
    with (noisy_scan / "scan_metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    fd_commands = np.load(fd_scan / "baseline_corrector_commands.npy")
    noisy_commands = np.load(noisy_scan / "baseline_corrector_commands.npy")
    fd_bpm = np.load(fd_scan / "reference_bpm_orbits.npy")
    noisy_bpm = np.load(noisy_scan / "reference_bpm_orbits.npy")
    fd_target = np.load(fd_scan / "reference_target_orbits.npy")
    noisy_target = np.load(noisy_scan / "reference_target_orbits.npy")
    reference_noise = np.load(noisy_scan / "baseline_reference_bpm_noise_m.npy")

    command_delta = noisy_commands - fd_commands
    bpm_delta = noisy_bpm - fd_bpm
    target_delta = noisy_target - fd_target
    target_delta_2d_rms_um = float(
        np.sqrt(np.mean(np.sum(target_delta**2, axis=-1))) * 1.0e6
    )
    removed_excess = 100.0 * (
        1.0 - (noisy_rmse - zero_rmse) / (uncorrected_rmse - zero_rmse)
    )
    summary: dict[str, object] = {
        "format": "cesr-gtpsa-noisy-correction-comparison-v1",
        "zero_offset_best_model": best_zero["model"],
        "zero_offset_best_rmse_2d_um": zero_rmse,
        "uncorrected_best_model": best_uncorrected["model"],
        "uncorrected_best_rmse_2d_um": uncorrected_rmse,
        "finite_difference_noiseless_best_model": best_fd["model"],
        "finite_difference_noiseless_best_rmse_2d_um": fd_rmse,
        "gtpsa_noisy_best_model": best_noisy["model"],
        "gtpsa_noisy_best_rmse_2d_um": noisy_rmse,
        "gtpsa_noisy_reduction_percent_vs_uncorrected": (
            100.0 * (1.0 - noisy_rmse / uncorrected_rmse)
        ),
        "gtpsa_noisy_excess_percent_vs_zero": (
            100.0 * (noisy_rmse / zero_rmse - 1.0)
        ),
        "gtpsa_noisy_change_percent_vs_fd_noiseless": (
            100.0 * (noisy_rmse / fd_rmse - 1.0)
        ),
        "gtpsa_noisy_removed_excess_rmse_percent": removed_excess,
        "uncorrected_zero_trained_joint_rmse_2d_um": metric(
            ood_joint(uncorrected_rows, UNCORRECTED), "rmse_2d_um"
        ),
        "fd_noiseless_zero_trained_joint_rmse_2d_um": metric(
            ood_joint(fd_rows, FD_CORRECTED), "rmse_2d_um"
        ),
        "gtpsa_noisy_zero_trained_joint_rmse_2d_um": metric(
            ood_joint(noisy_rows, GTPSA_NOISY_CORRECTED), "rmse_2d_um"
        ),
        "bpm_noise_rms_um_per_read": (
            1.0e6 * float(metadata["baseline_bpm_noise_rms_m_per_read"])
        ),
        "correction_measurement_repeats": int(
            metadata["baseline_measurement_repeats"]
        ),
        "expected_bpm_mean_noise_std_um": (
            1.0e6 * float(metadata["baseline_bpm_mean_noise_std_m"])
        ),
        "realized_reference_bpm_mean_noise_rms_um": rms(reference_noise) * 1.0e6,
        "gtpsa_vs_finite_difference_relative_l2_max": float(
            metadata["baseline_gtpsa_vs_finite_difference_relative_l2_max"]
        ),
        "gtpsa_vs_finite_difference_max_abs": float(
            metadata["baseline_gtpsa_vs_finite_difference_max_abs"]
        ),
        "gtpsa_response_closure_norm_max": float(
            metadata["baseline_response_closure_norm_max"]
        ),
        "baseline_command_delta_rms_urad_vs_fd_noiseless": (
            rms(command_delta) * 1.0e6
        ),
        "baseline_command_delta_max_abs_urad_vs_fd_noiseless": (
            float(np.max(np.abs(command_delta))) * 1.0e6
        ),
        "corrected_bpm_delta_rms_um_vs_fd_noiseless": rms(bpm_delta) * 1.0e6,
        "corrected_bpm_delta_max_abs_um_vs_fd_noiseless": (
            float(np.max(np.abs(bpm_delta))) * 1.0e6
        ),
        "corrected_target_delta_2d_rms_um_vs_fd_noiseless": (
            target_delta_2d_rms_um
        ),
        "corrected_target_delta_max_abs_um_vs_fd_noiseless": (
            float(np.max(np.abs(target_delta))) * 1.0e6
        ),
    }
    (noisy_dir / "gtpsa_noisy_protocol_comparison.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    table_rows = []
    for row in comparison_rows:
        table_rows.append(
            f"| {row['model']} | {row['zero_offset_rmse_2d_um']:.3f} | "
            f"{row['uncorrected_offset_rmse_2d_um']:.3f} | "
            f"{row['fd_noiseless_corrected_rmse_2d_um']:.3f} | "
            f"{row['gtpsa_noisy_corrected_rmse_2d_um']:.3f} | "
            f"{row['gtpsa_noisy_change_percent_vs_fd_noiseless']:+.3f}% |"
        )

    report = f"""# GTPSA-ORM plus noisy-reference correction comparison

All four columns use the same deterministic 16-machine latest-lattice SciBmad
ensemble, 10/3/3 machine split, static nuisance draws, sextupole scans, and
inverse definitions. The new protocol uses the zero-quadrupole-offset
SciBmad/GTPSA closed-orbit Jacobian, with fixed BPM and corrector gains, as its
103-control ORM. The stored reference and each current-orbit correction input
receive independent Gaussian BPM-noise means. The solved baseline command is
then frozen during all 76 sextupole scans.

| inverse | zero offset [um] | uncorrected offset [um] | FD ORM, noiseless correction [um] | GTPSA ORM, noisy correction [um] | GTPSA/noisy change vs FD/noiseless |
|---|---:|---:|---:|---:|---:|
{chr(10).join(table_rows)}

The correction inputs use `{summary['bpm_noise_rms_um_per_read']:.1f} um` RMS
noise per BPM plane and read, averaged over
`{summary['correction_measurement_repeats']:,}` reads. The expected noise of
each mean is only `{summary['expected_bpm_mean_noise_std_um']:.3f} um`; the
realized reference-noise RMS is
`{summary['realized_reference_bpm_mean_noise_rms_um']:.3f} um`.

The GTPSA ORM agrees with the central finite-difference check to maximum
relative L2 difference
`{summary['gtpsa_vs_finite_difference_relative_l2_max']:.3e}` across the 16
machines; its maximum periodic-response closure norm is
`{summary['gtpsa_response_closure_norm_max']:.3e}`. Relative to the finite-
difference/noiseless baseline, the new commands differ by only
`{summary['baseline_command_delta_rms_urad_vs_fd_noiseless']:.4f} urad` RMS,
the corrected BPM orbit by
`{summary['corrected_bpm_delta_rms_um_vs_fd_noiseless']:.4f} um` RMS, and the
corrected target orbit by
`{summary['corrected_target_delta_2d_rms_um_vs_fd_noiseless']:.4f} um` 2D RMS.

The best learned GTPSA/noisy result is
`{best_noisy['model']}` at `{noisy_rmse:.3f} um` 2D RMSE. It is
`{summary['gtpsa_noisy_reduction_percent_vs_uncorrected']:.3f}%` below the best
uncorrected result, `{summary['gtpsa_noisy_excess_percent_vs_zero']:.3f}%`
above the best zero-offset result, and
`{summary['gtpsa_noisy_change_percent_vs_fd_noiseless']:+.3f}%` relative to
the best finite-difference/noiseless corrected result. The last difference is
numerically negligible and is not evidence that adding noise improves the
inverse.

The strict tail gate still fails: the best GTPSA/noisy P99 is
`{metric(best_noisy, 'p99_2d_um'):.3f} um` and its worst-target RMSE is
`{metric(best_noisy, 'worst_target_rmse_2d_um'):.3f} um`. This experiment
therefore validates compatibility of the requested workflow at the tested
3,072-read averaging level. It does not establish robustness to single-shot
or low-repeat BPM noise, correlated BPM noise, missing channels, response
measurement error, or hardware effects. A repeat-count/covariance sweep is
required before selecting a CESR correction acquisition protocol.

The ORM is scaled with the realized simulated BPM and baseline-corrector gains,
so this is an exact-calibration/model-conditioned response test rather than an
unknown gain-mismatch test. The 103 baseline controls and 62 local-bump
controls also retain separate deterministic gain registries. That convention
is shared by the two corrected protocols and preserves their paired numerical
comparison, but it must be unified at the physical-device level before a
facility-facing calibration conclusion.
"""
    (noisy_dir / "GTPSA_NOISY_COMPARISON.md").write_text(
        report, encoding="utf-8"
    )
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
