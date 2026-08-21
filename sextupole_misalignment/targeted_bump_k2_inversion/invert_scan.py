#!/usr/bin/env python3
"""Invert one exact bump-by-K2 scan and compare the estimate with truth.

The baseline estimator uses the completed per-target GTPSA response map. For
each bump it first fits the full-ring K2 slope, reconstructs the two effective
local alignment/feeddown coordinates, adds the measured beam bump at the
target, and combines the bump estimates. Full-slope and two-mode views are
reported for several observable ablations.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESPONSE_MAP = HERE.parent / "response_map" / "results" / "full"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(source))


def observable_family(observable: str) -> str:
    if observable.startswith("orbit_"):
        return "orbit"
    if observable.startswith("phi_"):
        return "phase"
    if observable.startswith("beta_"):
        return "beta"
    if observable.startswith("alpha_"):
        return "alpha"
    if observable.startswith("c"):
        return "coupling"
    if observable.startswith("tune_"):
        return "tune"
    raise ValueError(observable)


ABLATIONS = {
    "orbit_only": {"orbit"},
    "orbit_phase_tune": {"orbit", "phase", "tune"},
    "orbit_phase_coupling_tune": {"orbit", "phase", "coupling", "tune"},
    "all_saved": {"orbit", "phase", "beta", "alpha", "coupling", "tune"},
}


def fit_slopes(
    observation_rows: list[dict[str, str]],
) -> tuple[dict[int, dict[tuple[str, str, str], float]], list[tuple[str, str, str]]]:
    grouped: dict[tuple[int, tuple[str, str, str]], list[tuple[float, float]]] = defaultdict(list)
    for row in observation_rows:
        bump = int(row["bump_index"])
        key = (row["observation_scope"], row["observation_name"], row["observable"])
        grouped[(bump, key)].append(
            (float(row["delta_k2_m3"]), float(row["observable_readback"]))
        )
    keys = sorted({key for _, key in grouped})
    slopes: dict[int, dict[tuple[str, str, str], float]] = defaultdict(dict)
    for (bump, key), points in grouped.items():
        points.sort()
        x = np.asarray([point[0] for point in points], dtype=float)
        y = np.asarray([point[1] for point in points], dtype=float)
        degree = min(2, len(points) - 1)
        design_columns = [np.ones_like(x), x]
        if degree >= 2:
            design_columns.append(0.5 * x * x)
        coefficients = np.linalg.lstsq(np.column_stack(design_columns), y, rcond=None)[0]
        slopes[bump][key] = float(coefficients[1])
    return dict(slopes), keys


def load_response_map(
    path: Path, target: str
) -> dict[tuple[str, str, str], tuple[float, float, float, str]]:
    response = {}
    for row in read_rows(path):
        if row["sextupole"].upper() != target.upper():
            continue
        key = (row["observation_scope"], row["observation_name"], row["observable"])
        response[key] = (
            float(row["d_k2"]),
            float(row["d2_k2_x"]),
            float(row["d2_k2_y"]),
            row["observable"],
        )
    if len(response) != 1191:
        raise ValueError(f"Expected 1191 response rows for {target}, found {len(response)}")
    return response


def load_scales(path: Path) -> dict[str, float]:
    return {
        row["observable"]: float(row["global_mixed_response_rms"])
        for row in read_rows(path)
    }


def achieved_bumps(state_rows: list[dict[str, str]]) -> dict[int, np.ndarray]:
    zero_k2 = [row for row in state_rows if abs(float(row["delta_k2_m3"])) < 1e-18]
    if not zero_k2:
        raise ValueError("Scan has no zero-K2 state")
    center = next(
        row for row in zero_k2
        if abs(float(row["bump_x_command_m"])) < 1e-18
        and abs(float(row["bump_y_command_m"])) < 1e-18
    )
    reference = np.asarray(
        [float(center["target_orbit_x_m"]), float(center["target_orbit_y_m"])], dtype=float
    )
    return {
        int(row["bump_index"]): np.asarray(
            [float(row["target_orbit_x_m"]), float(row["target_orbit_y_m"])], dtype=float
        )
        - reference
        for row in zero_k2
    }


def estimate_case(
    keys: list[tuple[str, str, str]],
    slopes: dict[int, dict[tuple[str, str, str], float]],
    bumps: dict[int, np.ndarray],
    response: dict[tuple[str, str, str], tuple[float, float, float, str]],
    scales: dict[str, float],
    families: set[str],
    representation: str,
) -> tuple[np.ndarray, float, int, list[dict[str, float | int]]]:
    selected = [
        key for key in keys
        if key in response and observable_family(key[2]) in families
    ]
    if not selected:
        raise ValueError("No observables selected")
    a = np.asarray([[response[key][1], response[key][2]] for key in selected], dtype=float)
    baseline = np.asarray([response[key][0] for key in selected], dtype=float)
    scale = np.asarray([scales[response[key][3]] for key in selected], dtype=float)
    a_white = a / scale[:, None]
    if representation == "modes2":
        u, _, _ = np.linalg.svd(a_white, full_matrices=False)
        projector = u[:, :2].T
    elif representation == "full_slopes":
        projector = np.eye(len(selected))
    else:
        raise ValueError(representation)

    design_blocks = []
    right_blocks = []
    per_bump = []
    for bump_index in sorted(slopes):
        g = np.asarray([slopes[bump_index][key] for key in selected], dtype=float)
        right = (g - baseline + a @ bumps[bump_index]) / scale
        design_blocks.append(projector @ a_white)
        right_blocks.append(projector @ right)

        local_effective = np.linalg.lstsq(a_white, (g - baseline) / scale, rcond=None)[0]
        local_center = local_effective + bumps[bump_index]
        per_bump.append(
            {
                "bump_index": bump_index,
                "effective_x_offset_m": float(local_effective[0]),
                "effective_y_offset_m": float(local_effective[1]),
                "center_x_estimate_m": float(local_center[0]),
                "center_y_estimate_m": float(local_center[1]),
            }
        )
    design = np.vstack(design_blocks)
    right = np.concatenate(right_blocks)
    estimate = np.linalg.lstsq(design, right, rcond=None)[0]
    residual = float(np.linalg.norm(design @ estimate - right) / np.sqrt(len(right)))
    return estimate, residual, len(selected), per_bump


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-dir", type=Path, default=HERE / "results" / "smoke_exact")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    scan_dir = args.scan_dir.resolve()
    output_dir = (args.output_dir or scan_dir / "inversion").resolve()
    observation_rows = read_rows(scan_dir / "scan_observations.csv")
    state_rows = read_rows(scan_dir / "scan_states.csv")
    target = state_rows[0]["target_sextupole"]
    truth = np.asarray(
        [float(state_rows[0]["true_x_offset_m"]), float(state_rows[0]["true_y_offset_m"])],
        dtype=float,
    )
    slopes, keys = fit_slopes(observation_rows)
    bumps = achieved_bumps(state_rows)
    response = load_response_map(RESPONSE_MAP / "alignment_coefficients.csv", target)
    scales = load_scales(RESPONSE_MAP / "local_response_svd_scales.csv")

    estimates = []
    per_bump_rows = []
    for ablation, families in ABLATIONS.items():
        for representation in ("full_slopes", "modes2"):
            estimate, residual, count, per_bump = estimate_case(
                keys, slopes, bumps, response, scales, families, representation
            )
            estimates.append(
                {
                    "target_sextupole": target,
                    "observable_ablation": ablation,
                    "input_representation": representation,
                    "observable_count": count,
                    "true_x_offset_m": truth[0],
                    "true_y_offset_m": truth[1],
                    "estimated_x_offset_m": estimate[0],
                    "estimated_y_offset_m": estimate[1],
                    "error_x_m": estimate[0] - truth[0],
                    "error_y_m": estimate[1] - truth[1],
                    "absolute_error_2d_m": float(np.linalg.norm(estimate - truth)),
                    "whitened_residual_rms": residual,
                }
            )
            for row in per_bump:
                per_bump_rows.append(
                    {
                        "target_sextupole": target,
                        "observable_ablation": ablation,
                        "input_representation": representation,
                        **row,
                    }
                )
    write_csv(output_dir / "offset_estimates.csv", estimates)
    write_csv(output_dir / "per_bump_feeddown_estimates.csv", per_bump_rows)
    summary = {
        "format": "cesr-targeted-sextupole-inversion-v1",
        "target_sextupole": target,
        "truth_m": {"x": float(truth[0]), "y": float(truth[1])},
        "best_case": min(estimates, key=lambda row: float(row["absolute_error_2d_m"])),
        "interpretation_boundary": (
            "observable_rms scaling is structural; experimental precision requires measured covariance"
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "inversion_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Estimates: {output_dir / 'offset_estimates.csv'}")
    print(f"Summary: {output_dir / 'inversion_summary.json'}")


if __name__ == "__main__":
    main()
