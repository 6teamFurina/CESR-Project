#!/usr/bin/env python3
"""Compare paired P0, P1, P2a, and P2b inverses on one saved exact scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import invert_scan as base


HERE = Path(__file__).resolve().parent
DEFAULT_SCAN = HERE / "results" / "smoke_background"
DEFAULT_MODELS = HERE / "results" / "paired_benchmark" / "conditioned_models"
DEFAULT_OUTPUT = HERE / "results" / "paired_benchmark"


def row_key(row: dict[str, str]) -> tuple[str, str, str]:
    return row["observation_scope"], row["observation_name"], row["observable"]


def conditioned_by_bump(rows: list[dict[str, str]]) -> dict[int, dict[tuple[str, str, str], dict[str, str]]]:
    result: dict[int, dict[tuple[str, str, str], dict[str, str]]] = {}
    for row in rows:
        result.setdefault(int(row["bump_index"]), {})[row_key(row)] = row
    return result


def estimate_p1(
    slopes: dict[int, dict[tuple[str, str, str], float]],
    conditioned: dict[int, dict[tuple[str, str, str], dict[str, str]]],
    scales: dict[str, float],
    families: set[str],
) -> tuple[np.ndarray, float, int]:
    design_blocks, right_blocks = [], []
    count = 0
    for bump in sorted(slopes):
        rows = conditioned[bump]
        keys = [
            key for key in sorted(rows)
            if base.observable_family(key[2]) in families
        ]
        scale = np.asarray([scales[key[2]] for key in keys])
        design_blocks.append(np.asarray([
            [float(rows[key]["d2_k2_x"]), float(rows[key]["d2_k2_y"])]
            for key in keys
        ]) / scale[:, None])
        right_blocks.append(np.asarray([
            slopes[bump][key] - float(rows[key]["reference_k2_slope"])
            for key in keys
        ]) / scale)
        count = len(keys)
    design = np.vstack(design_blocks)
    right = np.concatenate(right_blocks)
    estimate = np.linalg.lstsq(design, right, rcond=None)[0]
    residual = float(np.linalg.norm(design @ estimate - right) / np.sqrt(len(right)))
    return estimate, residual, count


def estimate_p2(
    slopes: dict[int, dict[tuple[str, str, str], float]],
    conditioned: dict[int, dict[tuple[str, str, str], dict[str, str]]],
    sources: dict[int, dict[tuple[str, str, str], dict[str, str]]],
    scales: dict[str, float],
    families: set[str],
    source_columns: tuple[str, ...],
) -> tuple[np.ndarray, float, int, list[dict[str, object]]]:
    design_blocks, right_blocks = [], []
    kick_rows: list[dict[str, object]] = []
    count = 0
    for bump in sorted(slopes):
        response_rows = conditioned[bump]
        source_rows = sources[bump]
        keys = [
            key for key in sorted(response_rows)
            if base.observable_family(key[2]) in families
        ]
        scale = np.asarray([scales[key[2]] for key in keys])
        a = np.asarray([
            [float(response_rows[key]["d2_k2_x"]), float(response_rows[key]["d2_k2_y"])]
            for key in keys
        ]) / scale[:, None]
        right = np.asarray([
            slopes[bump][key] - float(response_rows[key]["reference_k2_slope"])
            for key in keys
        ]) / scale
        source = np.asarray([
            [float(source_rows[key][column]) for column in source_columns]
            for key in keys
        ]) / scale[:, None]
        u, singular, vt = np.linalg.svd(source, full_matrices=False)
        threshold = np.finfo(float).eps * max(source.shape) * singular[0]
        rank = int(np.sum(singular > threshold))
        u = u[:, :rank]
        singular = singular[:rank]
        vt = vt[:rank]
        q = vt.T @ ((u.T @ right) / singular)
        b = vt.T @ ((u.T @ a) / singular[:, None])
        # Whiten the reconstructed local-source covariance before the offset fit.
        whitener = singular[:, None] * vt
        design_blocks.append(whitener @ b)
        right_blocks.append(whitener @ q)
        kick_rows.append({
            "bump_index": bump,
            "source_rank": rank,
            **{f"estimated_{name}_per_k2": float(value) for name, value in zip(source_columns, q)},
        })
        count = len(keys)
    design = np.vstack(design_blocks)
    right = np.concatenate(right_blocks)
    estimate = np.linalg.lstsq(design, right, rcond=None)[0]
    residual = float(np.linalg.norm(design @ estimate - right) / np.sqrt(len(right)))
    return estimate, residual, count, kick_rows


def result_row(
    method: str,
    view: str,
    estimate: np.ndarray,
    truth: np.ndarray,
    residual: float,
    count: int,
    conditioning: str,
) -> dict[str, object]:
    error = estimate - truth
    return {
        "method": method,
        "observable_view": view,
        "conditioning": conditioning,
        "observable_count_per_bump": count,
        "true_x_offset_m": float(truth[0]),
        "true_y_offset_m": float(truth[1]),
        "estimated_x_offset_m": float(estimate[0]),
        "estimated_y_offset_m": float(estimate[1]),
        "error_x_m": float(error[0]),
        "error_y_m": float(error[1]),
        "absolute_error_2d_m": float(np.linalg.norm(error)),
        "fit_residual_rms": residual,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-dir", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--models-dir", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--conditioning-label", default="other-offset oracle")
    parser.add_argument("--p1-method", default="P1_background_conditioned_mixed")
    args = parser.parse_args()
    scan_dir, models_dir, output_dir = map(Path.resolve, (args.scan_dir, args.models_dir, args.output_dir))
    observations = base.read_rows(scan_dir / "scan_observations.csv")
    states = base.read_rows(scan_dir / "scan_states.csv")
    target = states[0]["target_sextupole"]
    truth = np.asarray([float(states[0]["true_x_offset_m"]), float(states[0]["true_y_offset_m"])])
    slopes, keys = base.fit_slopes(observations)
    bumps = base.achieved_bumps(states)
    nominal = base.load_response_map(base.RESPONSE_MAP / "alignment_coefficients.csv", target)
    scales = base.load_scales(base.RESPONSE_MAP / "local_response_svd_scales.csv")
    conditioned = conditioned_by_bump(base.read_rows(models_dir / "conditioned_mixed_response.csv"))
    sources = conditioned_by_bump(base.read_rows(models_dir / "local_source_response.csv"))

    results: list[dict[str, object]] = []
    for view, families in {
        "orbit_only": {"orbit"},
        "orbit_phase_coupling_tune": {"orbit", "phase", "coupling", "tune"},
    }.items():
        estimate, residual, count, _ = base.estimate_case(
            keys, slopes, bumps, nominal, scales, families, "full_slopes"
        )
        results.append(result_row("P0_nominal_mixed_GTPSA", view, estimate, truth, residual, count, "nominal"))
        estimate, residual, count = estimate_p1(slopes, conditioned, scales, families)
        results.append(result_row(args.p1_method, view, estimate, truth, residual, count, args.conditioning_label))

    p2a, residual, count, kicks2 = estimate_p2(
        slopes, conditioned, sources, scales, {"orbit"}, ("d_kn0l", "d_ks0l")
    )
    results.append(result_row("P2a_two_local_dipole_kicks", "orbit_only", p2a, truth, residual, count, args.conditioning_label))
    p2b, residual, count, kicks4 = estimate_p2(
        slopes, conditioned, sources, scales,
        {"orbit", "phase", "coupling", "tune"},
        ("d_kn0l", "d_ks0l", "d_kn1l", "d_ks1l"),
    )
    results.append(result_row("P2b_four_local_kicks", "orbit_phase_coupling_tune", p2b, truth, residual, count, args.conditioning_label))

    output_dir.mkdir(parents=True, exist_ok=True)
    base.write_csv(output_dir / "p0_p2_offset_estimates.csv", results)
    base.write_csv(output_dir / "p2a_local_kicks.csv", kicks2)
    base.write_csv(output_dir / "p2b_local_kicks.csv", kicks4)
    summary = {
        "format": "cesr-paired-p0-p2-benchmark-v1",
        "target": target,
        "scan_dir": str(scan_dir),
        "paired_input_semantics": "all methods consume the same exact scan and saved offset realization",
        "results": results,
    }
    (output_dir / "p0_p2_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    for row in results:
        print(f"{row['method']:38s} {row['observable_view']:28s} error={1e6*float(row['absolute_error_2d_m']):9.3f} um")


if __name__ == "__main__":
    main()
