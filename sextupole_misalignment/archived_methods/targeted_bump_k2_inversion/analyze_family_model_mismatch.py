#!/usr/bin/env python3
"""Measure family-wise mismatch of the nominal local inverse at known truth."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

import invert_scan as inv


HERE = Path(__file__).resolve().parent


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze_case(case: str) -> list[dict[str, object]]:
    root = HERE / "results" / case
    observations = inv.read_rows(root / "scan_observations.csv")
    states = inv.read_rows(root / "scan_states.csv")
    target = states[0]["target_sextupole"]
    truth = np.asarray(
        [float(states[0]["true_x_offset_m"]), float(states[0]["true_y_offset_m"])],
        dtype=float,
    )
    slopes, keys = inv.fit_slopes(observations)
    bumps = inv.achieved_bumps(states)
    response = inv.load_response_map(
        inv.RESPONSE_MAP / "alignment_coefficients.csv", target
    )
    scales = inv.load_scales(inv.RESPONSE_MAP / "local_response_svd_scales.csv")
    rows = []
    for family in ("orbit", "phase", "beta", "alpha", "coupling", "tune"):
        selected = [
            key for key in keys
            if key in response and inv.observable_family(key[2]) == family
        ]
        a = np.asarray(
            [[response[key][1], response[key][2]] for key in selected], dtype=float
        )
        baseline = np.asarray([response[key][0] for key in selected], dtype=float)
        scale = np.asarray([scales[response[key][3]] for key in selected], dtype=float)
        design = []
        measured = []
        for bump_index in sorted(slopes):
            g = np.asarray([slopes[bump_index][key] for key in selected], dtype=float)
            design.append(a / scale[:, None])
            measured.append((g - baseline + a @ bumps[bump_index]) / scale)
        design_matrix = np.vstack(design)
        measured_vector = np.concatenate(measured)
        truth_prediction = design_matrix @ truth
        residual = measured_vector - truth_prediction
        estimate = np.linalg.lstsq(design_matrix, measured_vector, rcond=None)[0]
        signal_rms = float(np.linalg.norm(truth_prediction) / np.sqrt(len(residual)))
        mismatch_rms = float(np.linalg.norm(residual) / np.sqrt(len(residual)))
        rows.append(
            {
                "case": case,
                "observable_family": family,
                "observation_count_per_bump": len(selected),
                "truth_signal_rms_structural_units": signal_rms,
                "truth_mismatch_rms_structural_units": mismatch_rms,
                "mismatch_over_truth_signal": mismatch_rms / signal_rms,
                "family_estimated_x_offset_um": 1e6 * float(estimate[0]),
                "family_estimated_y_offset_um": 1e6 * float(estimate[1]),
                "family_error_2d_um": 1e6 * float(np.linalg.norm(estimate - truth)),
            }
        )
    return rows


def main() -> None:
    rows = analyze_case("smoke_exact") + analyze_case("smoke_background")
    path = HERE / "results" / "family_model_mismatch.csv"
    write_csv(path, rows)
    for row in rows:
        print(
            f"{row['case']:16s} {row['observable_family']:8s} "
            f"mismatch/signal={row['mismatch_over_truth_signal']:.3f} "
            f"family error={row['family_error_2d_um']:.1f} um"
        )
    print(f"Output: {path}")


if __name__ == "__main__":
    main()

