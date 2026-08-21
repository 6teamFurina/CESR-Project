#!/usr/bin/env python3
"""Compute the thin SVD of each sextupole's 2-by-n mixed-response matrix.

For one sextupole, the two rows are d2O/(dKn2 dx_offset) and
d2O/(dKn2 dy_offset).  The columns are the 1,191 saved detector/ring
observables.  This is an observable-compression analysis; it is not a fit to
position samples.

Two scalings are retained:

* ``raw``: the saved coefficients without combining-unit normalization.  This
  is an algebraic reference only.
* ``observable_rms``: each named observable (for example orbit_x or beta_1)
  is divided by one global RMS response scale computed across all 76 magnets,
  locations, and both offset rows.  This gives a useful structural comparison,
  but is not a substitute for whitening with measured covariance.
"""

from __future__ import annotations

import csv
import json
import math
from collections import OrderedDict, defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results" / "full"
COEFFICIENTS = RESULTS / "alignment_coefficients.csv"
SUMMARY_CSV = RESULTS / "local_response_svd_summary.csv"
MODES_CSV = RESULTS / "local_response_svd_modes.csv"
SCALES_CSV = RESULTS / "local_response_svd_scales.csv"
SUMMARY_JSON = RESULTS / "local_response_svd_summary.json"

EXPECTED_SEXTUPOLES = 76
EXPECTED_OBSERVATIONS = 1191
SCHEMES = ("raw", "observable_rms")


def percentile(values: list[float], probability: float) -> float:
    return float(np.quantile(np.asarray(values, dtype=float), probability))


def canonicalize_svd(
    parameter_rotation: np.ndarray, mode_weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Fix arbitrary SVD signs using the largest mode-weight magnitude."""
    parameter_rotation = parameter_rotation.copy()
    mode_weights = mode_weights.copy()
    for mode in range(mode_weights.shape[0]):
        pivot = int(np.argmax(np.abs(mode_weights[mode])))
        if mode_weights[mode, pivot] < 0:
            mode_weights[mode] *= -1
            parameter_rotation[:, mode] *= -1
    return parameter_rotation, mode_weights


def main() -> None:
    if not COEFFICIENTS.is_file():
        raise FileNotFoundError(COEFFICIENTS)

    rows_by_sextupole: OrderedDict[str, list[dict[str, str]]] = OrderedDict()
    observable_squares: dict[str, float] = defaultdict(float)
    observable_value_counts: dict[str, int] = defaultdict(int)

    with COEFFICIENTS.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        required = {
            "sextupole",
            "sextupole_index",
            "sextupole_s_m",
            "observation_scope",
            "observation_name",
            "observable",
            "d2_k2_x",
            "d2_k2_y",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}")
        for row in reader:
            name = row["sextupole"]
            rows_by_sextupole.setdefault(name, []).append(row)
            observable = row["observable"]
            x_value = float(row["d2_k2_x"])
            y_value = float(row["d2_k2_y"])
            if not math.isfinite(x_value) or not math.isfinite(y_value):
                raise ValueError(f"Non-finite mixed response in {name}: {row}")
            observable_squares[observable] += x_value * x_value + y_value * y_value
            observable_value_counts[observable] += 2

    if len(rows_by_sextupole) != EXPECTED_SEXTUPOLES:
        raise ValueError(
            f"Expected {EXPECTED_SEXTUPOLES} sextupoles, found {len(rows_by_sextupole)}"
        )
    counts = {len(rows) for rows in rows_by_sextupole.values()}
    if counts != {EXPECTED_OBSERVATIONS}:
        raise ValueError(
            f"Expected {EXPECTED_OBSERVATIONS} observations per sextupole, found {sorted(counts)}"
        )

    reference_rows = next(iter(rows_by_sextupole.values()))
    observation_keys = [
        (row["observation_scope"], row["observation_name"], row["observable"])
        for row in reference_rows
    ]
    if len(set(observation_keys)) != EXPECTED_OBSERVATIONS:
        raise ValueError("Observation keys are not unique within the reference sextupole")
    for name, rows in rows_by_sextupole.items():
        keys = [
            (row["observation_scope"], row["observation_name"], row["observable"])
            for row in rows
        ]
        if keys != observation_keys:
            raise ValueError(f"Observation order mismatch for {name}")

    observable_scales = {
        observable: math.sqrt(observable_squares[observable] / observable_value_counts[observable])
        for observable in observable_squares
    }
    zero_scales = [name for name, value in observable_scales.items() if not value > 0]
    if zero_scales:
        raise ValueError(f"Zero global response scale for observables: {zero_scales}")

    with SCALES_CSV.open("w", newline="", encoding="utf-8") as output:
        fieldnames = ["observable", "global_mixed_response_rms"]
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        for observable in sorted(observable_scales):
            writer.writerow(
                {
                    "observable": observable,
                    "global_mixed_response_rms": f"{observable_scales[observable]:.17g}",
                }
            )

    summary_fieldnames = [
        "sextupole",
        "sextupole_index",
        "sextupole_s_m",
        "scaling",
        "observation_count",
        "sigma_1",
        "sigma_2",
        "condition_number",
        "effective_rank",
        "response_x_norm",
        "response_y_norm",
        "response_x_y_cosine",
        "parameter_rotation_x_mode1",
        "parameter_rotation_y_mode1",
        "parameter_rotation_x_mode2",
        "parameter_rotation_y_mode2",
        "relative_reconstruction_error",
        "parameter_rotation_orthogonality_error",
        "mode_weight_orthogonality_error",
    ]
    mode_fieldnames = [
        "sextupole",
        "scaling",
        "mode",
        "observation_index",
        "observation_scope",
        "observation_name",
        "observable",
        "mode_weight",
    ]

    summary_rows: list[dict[str, object]] = []
    with MODES_CSV.open("w", newline="", encoding="utf-8") as modes_output:
        modes_writer = csv.DictWriter(modes_output, fieldnames=mode_fieldnames)
        modes_writer.writeheader()
        for sextupole, rows in rows_by_sextupole.items():
            base_response = np.asarray(
                [
                    [float(row["d2_k2_x"]) for row in rows],
                    [float(row["d2_k2_y"]) for row in rows],
                ],
                dtype=float,
            )
            for scheme in SCHEMES:
                if scheme == "raw":
                    response = base_response
                else:
                    scales = np.asarray(
                        [observable_scales[row["observable"]] for row in rows], dtype=float
                    )
                    response = base_response / scales[np.newaxis, :]

                parameter_rotation, singular_values, mode_weights = np.linalg.svd(
                    response, full_matrices=False
                )
                parameter_rotation, mode_weights = canonicalize_svd(
                    parameter_rotation, mode_weights
                )
                reconstruction = (
                    parameter_rotation @ np.diag(singular_values) @ mode_weights
                )
                response_norm = float(np.linalg.norm(response))
                reconstruction_error = float(
                    np.linalg.norm(reconstruction - response) / response_norm
                )
                parameter_orthogonality = float(
                    np.linalg.norm(parameter_rotation.T @ parameter_rotation - np.eye(2))
                )
                mode_orthogonality = float(
                    np.linalg.norm(mode_weights @ mode_weights.T - np.eye(2))
                )
                tolerance = max(response.shape) * np.finfo(float).eps * singular_values[0]
                effective_rank = int(np.sum(singular_values > tolerance))
                x_norm = float(np.linalg.norm(response[0]))
                y_norm = float(np.linalg.norm(response[1]))
                cosine = float(np.dot(response[0], response[1]) / (x_norm * y_norm))
                condition = float(singular_values[0] / singular_values[1])

                if reconstruction_error > 1e-12:
                    raise ValueError(
                        f"SVD reconstruction failed for {sextupole}/{scheme}: "
                        f"{reconstruction_error:.3e}"
                    )
                if parameter_orthogonality > 1e-12 or mode_orthogonality > 1e-12:
                    raise ValueError(
                        f"SVD orthogonality failed for {sextupole}/{scheme}: "
                        f"{parameter_orthogonality:.3e}, {mode_orthogonality:.3e}"
                    )

                first = rows[0]
                summary_rows.append(
                    {
                        "sextupole": sextupole,
                        "sextupole_index": first["sextupole_index"],
                        "sextupole_s_m": first["sextupole_s_m"],
                        "scaling": scheme,
                        "observation_count": len(rows),
                        "sigma_1": f"{singular_values[0]:.17g}",
                        "sigma_2": f"{singular_values[1]:.17g}",
                        "condition_number": f"{condition:.17g}",
                        "effective_rank": effective_rank,
                        "response_x_norm": f"{x_norm:.17g}",
                        "response_y_norm": f"{y_norm:.17g}",
                        "response_x_y_cosine": f"{cosine:.17g}",
                        "parameter_rotation_x_mode1": f"{parameter_rotation[0, 0]:.17g}",
                        "parameter_rotation_y_mode1": f"{parameter_rotation[1, 0]:.17g}",
                        "parameter_rotation_x_mode2": f"{parameter_rotation[0, 1]:.17g}",
                        "parameter_rotation_y_mode2": f"{parameter_rotation[1, 1]:.17g}",
                        "relative_reconstruction_error": f"{reconstruction_error:.17g}",
                        "parameter_rotation_orthogonality_error": f"{parameter_orthogonality:.17g}",
                        "mode_weight_orthogonality_error": f"{mode_orthogonality:.17g}",
                    }
                )

                for mode in range(2):
                    for observation_index, (row, weight) in enumerate(
                        zip(rows, mode_weights[mode], strict=True), start=1
                    ):
                        modes_writer.writerow(
                            {
                                "sextupole": sextupole,
                                "scaling": scheme,
                                "mode": mode + 1,
                                "observation_index": observation_index,
                                "observation_scope": row["observation_scope"],
                                "observation_name": row["observation_name"],
                                "observable": row["observable"],
                                "mode_weight": f"{weight:.17g}",
                            }
                        )

    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=summary_fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    scheme_summaries: dict[str, dict[str, object]] = {}
    for scheme in SCHEMES:
        selected = [row for row in summary_rows if row["scaling"] == scheme]
        conditions = [float(row["condition_number"]) for row in selected]
        sigma_ratios = [float(row["sigma_2"]) / float(row["sigma_1"]) for row in selected]
        worst = sorted(
            (
                {
                    "sextupole": str(row["sextupole"]),
                    "condition_number": float(row["condition_number"]),
                    "sigma_2_over_sigma_1": float(row["sigma_2"]) / float(row["sigma_1"]),
                }
                for row in selected
            ),
            key=lambda item: item["condition_number"],
            reverse=True,
        )[:10]
        scheme_summaries[scheme] = {
            "condition_number": {
                "minimum": min(conditions),
                "p10": percentile(conditions, 0.10),
                "median": float(np.median(conditions)),
                "p90": percentile(conditions, 0.90),
                "maximum": max(conditions),
            },
            "sigma_2_over_sigma_1": {
                "minimum": min(sigma_ratios),
                "p10": percentile(sigma_ratios, 0.10),
                "median": float(np.median(sigma_ratios)),
                "p90": percentile(sigma_ratios, 0.90),
                "maximum": max(sigma_ratios),
            },
            "effective_rank_counts": {
                str(rank): sum(int(row["effective_rank"]) == rank for row in selected)
                for rank in sorted({int(row["effective_rank"]) for row in selected})
            },
            "largest_condition_numbers": worst,
        }

    result = {
        "format": "cesr-sextupole-local-response-svd-v1",
        "matrix_convention": (
            "R_j has shape 2 x 1191; rows are d2O/(dKn2 dx_offset) and "
            "d2O/(dKn2 dy_offset); columns are saved observables"
        ),
        "svd_convention": "R_j = P_j diag(sigma_1,sigma_2) Q_j^T",
        "sextupoles": len(rows_by_sextupole),
        "observations_per_sextupole": EXPECTED_OBSERVATIONS,
        "scalings": {
            "raw": "No scaling; mixed physical units make cross-observable interpretation invalid.",
            "observable_rms": (
                "Each named observable is divided by its global RMS mixed response across all "
                "magnets, locations, and x/y rows. Structural only; not measurement-noise whitening."
            ),
        },
        "projection": (
            "For a consistently scaled n-vector y_j, compressed amplitudes are "
            "z_j = Q_j^T y_j and the local model is z_j = diag(sigma) P_j^T [dx_j,dy_j]."
        ),
        "scheme_summaries": scheme_summaries,
        "artifacts": {
            "summary_csv": str(SUMMARY_CSV),
            "modes_csv": str(MODES_CSV),
            "scales_csv": str(SCALES_CSV),
        },
        "interpretation_boundary": (
            "Experimental alignment precision requires replacing observable_rms scaling with "
            "measured covariance whitening and adding nuisance-response columns."
        ),
    }
    SUMMARY_JSON.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(f"Computed {len(summary_rows)} local SVD summaries")
    print(f"Summary: {SUMMARY_CSV}")
    print(f"Modes: {MODES_CSV}")
    print(f"Scales: {SCALES_CSV}")
    print(f"Machine-readable summary: {SUMMARY_JSON}")


if __name__ == "__main__":
    main()
