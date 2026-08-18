#!/usr/bin/env python3
"""Relate existing per-sextupole detector-vector contributions to beta and phase."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import tomllib
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
ERROR_ANALYSIS = HERE.parent
DEFAULT_RESULTS = HERE / "results"
DEFAULT_ORBITS = (
    ERROR_ANALYSIS
    / "quadratic_x_attribution"
    / "element_results"
    / "element_exposure_directions.csv"
)
DEFAULT_THICK = ERROR_ANALYSIS / "thick_element_sextupole_sourcing"
PREDICTORS = (
    ("source_only", "predictor_source_only_rad"),
    ("nominal_beta_envelope", "predictor_nominal_beta_envelope_m"),
    ("direction_beta_envelope", "predictor_direction_beta_envelope_m"),
    ("nominal_beta_phase", "predictor_nominal_beta_phase_m"),
    ("direction_beta_phase", "predictor_direction_beta_phase_m"),
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--optics", type=Path, default=DEFAULT_RESULTS / "nominal_optics_points.csv")
    parser.add_argument(
        "--optics-metadata",
        type=Path,
        default=DEFAULT_RESULTS / "nominal_optics_metadata.toml",
    )
    parser.add_argument(
        "--direction-optics",
        type=Path,
        default=DEFAULT_RESULTS / "direction_optics_points.csv",
    )
    parser.add_argument(
        "--direction-tunes",
        type=Path,
        default=DEFAULT_RESULTS / "direction_optics_tunes.csv",
    )
    parser.add_argument("--orbits", type=Path, default=DEFAULT_ORBITS)
    parser.add_argument(
        "--horizontal-contributions",
        type=Path,
        default=DEFAULT_THICK / "horizontal_results" / "thick_sextupole_direction_contributions.csv",
    )
    parser.add_argument(
        "--horizontal-directions",
        type=Path,
        default=DEFAULT_THICK / "horizontal_results" / "direction_closure.csv",
    )
    parser.add_argument(
        "--vertical-contributions",
        type=Path,
        default=DEFAULT_THICK / "vertical_results" / "thick_sextupole_direction_contributions.csv",
    )
    parser.add_argument(
        "--vertical-directions",
        type=Path,
        default=DEFAULT_THICK / "vertical_results" / "direction_closure.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        return math.nan
    mean_left = statistics.fmean(left)
    mean_right = statistics.fmean(right)
    dl = [value - mean_left for value in left]
    dr = [value - mean_right for value in right]
    denominator = math.sqrt(sum(value * value for value in dl) * sum(value * value for value in dr))
    return sum(a * b for a, b in zip(dl, dr)) / denominator if denominator else math.nan


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        average_rank = 0.5 * (cursor + end - 1) + 1.0
        for position in range(cursor, end):
            result[order[position]] = average_rank
        cursor = end
    return result


def spearman(left: list[float], right: list[float]) -> float:
    return pearson(ranks(left), ranks(right))


def log_pearson(left: list[float], right: list[float]) -> float:
    pairs = [(a, b) for a, b in zip(left, right) if a > 0.0 and b > 0.0]
    if len(pairs) < 2:
        return math.nan
    return pearson(
        [math.log10(a) for a, _ in pairs],
        [math.log10(b) for _, b in pairs],
    )


def quantile(values: list[float], probability: float) -> float:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return math.nan
    position = probability * (len(clean) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return clean[lower]
    fraction = position - lower
    return clean[lower] * (1.0 - fraction) + clean[upper] * fraction


def optics_products(
    optics_path: Path,
    metadata_path: Path,
) -> tuple[list[dict[str, object]], dict[str, dict[str, float]]]:
    points = read_rows(optics_path)
    with metadata_path.open("rb") as stream:
        metadata = tomllib.load(stream)
    tunes = {
        "x": float(metadata["full_tune_1_turn"]),
        "y": float(metadata["full_tune_2_turn"]),
    }
    sextupoles = [row for row in points if row["point_type"] == "sextupole_exit"]
    detectors = [row for row in points if row["point_type"] == "detector"]
    if len(sextupoles) != 76 or len(detectors) != 99:
        raise ValueError(f"Expected 76 sextupoles and 99 detectors, got {len(sextupoles)} and {len(detectors)}")

    pair_rows: list[dict[str, object]] = []
    sums: dict[str, dict[str, float]] = {
        row["element_name"].lower(): {
            "envelope_x_sq": 0.0,
            "phase_x_sq": 0.0,
            "envelope_y_sq": 0.0,
            "phase_y_sq": 0.0,
        }
        for row in sextupoles
    }
    for sextupole in sextupoles:
        sext_name = sextupole["element_name"].lower()
        for detector in detectors:
            output: dict[str, object] = {
                "sextupole_name": sextupole["element_name"],
                "sextupole_index": int(sextupole["element_index"]),
                "sextupole_exit_s_m": float(sextupole["s_m"]),
                "k2l_m2": float(sextupole["k2l_m2"]),
                "detector_name": detector["element_name"],
                "detector_index": int(detector["element_index"]),
                "detector_s_m": float(detector["s_m"]),
            }
            for plane, mode in (("x", "1"), ("y", "2")):
                beta_s = float(sextupole[f"beta_{mode}_m"])
                beta_d = float(detector[f"beta_{mode}_m"])
                phi_s = float(sextupole[f"phi_{mode}_turn"])
                phi_d = float(detector[f"phi_{mode}_turn"])
                tune = tunes[plane]
                denominator_signed = 2.0 * math.sin(math.pi * tune)
                if denominator_signed == 0.0:
                    raise ValueError(f"Mode-{mode} tune lies on an integer: {tune}")
                beta_product_sqrt = math.sqrt(beta_s * beta_d)
                phase_advance = abs(phi_d - phi_s)
                envelope = beta_product_sqrt / abs(denominator_signed)
                response = beta_product_sqrt / denominator_signed * math.cos(
                    2.0 * math.pi * phase_advance - math.pi * tune
                )
                output.update(
                    {
                        f"beta_{plane}_sext_m": beta_s,
                        f"beta_{plane}_detector_m": beta_d,
                        f"phi_{plane}_sext_turn": phi_s,
                        f"phi_{plane}_detector_turn": phi_d,
                        f"abs_phase_advance_{plane}_turn": phase_advance,
                        f"closed_orbit_envelope_{plane}_m": envelope,
                        f"closed_orbit_response_{plane}_m": response,
                    }
                )
                sums[sext_name][f"envelope_{plane}_sq"] += envelope * envelope
                sums[sext_name][f"phase_{plane}_sq"] += response * response
            pair_rows.append(output)
    factors: dict[str, dict[str, float]] = {}
    for sextupole in sextupoles:
        name = sextupole["element_name"].lower()
        factors[name] = {
            "k2l_m2": float(sextupole["k2l_m2"]),
            "beta_x_sext_m": float(sextupole["beta_1_m"]),
            "beta_y_sext_m": float(sextupole["beta_2_m"]),
            "envelope_x_l2_m": math.sqrt(sums[name]["envelope_x_sq"]),
            "phase_x_l2_m": math.sqrt(sums[name]["phase_x_sq"]),
            "envelope_y_l2_m": math.sqrt(sums[name]["envelope_y_sq"]),
            "phase_y_l2_m": math.sqrt(sums[name]["phase_y_sq"]),
        }
    return pair_rows, factors


def direction_optics_products(
    optics_path: Path,
    tunes_path: Path,
) -> tuple[list[dict[str, object]], dict[tuple[int, str], dict[str, float]]]:
    points_by_trial: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in read_rows(optics_path):
        points_by_trial[int(row["trial"])].append(row)
    tunes = {
        int(row["trial"]): {
            "x": float(row["full_tune_1_turn"]),
            "y": float(row["full_tune_2_turn"]),
        }
        for row in read_rows(tunes_path)
    }
    if len(points_by_trial) != 100 or set(points_by_trial) != set(tunes):
        raise ValueError("Direction optics and tune tables must contain the same 100 trials")

    rows: list[dict[str, object]] = []
    factors: dict[tuple[int, str], dict[str, float]] = {}
    for trial in sorted(points_by_trial):
        points = points_by_trial[trial]
        sextupoles = [row for row in points if row["point_type"] == "sextupole_exit"]
        detectors = [row for row in points if row["point_type"] == "detector"]
        if len(sextupoles) != 76 or len(detectors) != 99:
            raise ValueError(
                f"Trial {trial}: expected 76 sextupoles and 99 detectors, got "
                f"{len(sextupoles)} and {len(detectors)}"
            )
        for sextupole in sextupoles:
            name = sextupole["element_name"].lower()
            result: dict[str, object] = {
                "trial": trial,
                "element_index": int(sextupole["element_index"]),
                "element_name": sextupole["element_name"],
                "s_m": float(sextupole["s_m"]),
                "k2l_m2": float(sextupole["k2l_m2"]),
                "beta_x_sext_m": float(sextupole["beta_1_m"]),
                "beta_y_sext_m": float(sextupole["beta_2_m"]),
                "full_tune_x_turn": tunes[trial]["x"],
                "full_tune_y_turn": tunes[trial]["y"],
            }
            factor: dict[str, float] = {
                "k2l_m2": float(sextupole["k2l_m2"]),
                "beta_x_sext_m": float(sextupole["beta_1_m"]),
                "beta_y_sext_m": float(sextupole["beta_2_m"]),
            }
            for plane, mode in (("x", "1"), ("y", "2")):
                tune = tunes[trial][plane]
                denominator_signed = 2.0 * math.sin(math.pi * tune)
                beta_s = float(sextupole[f"beta_{mode}_m"])
                phi_s = float(sextupole[f"phi_{mode}_turn"])
                envelope_sq = 0.0
                phase_sq = 0.0
                for detector in detectors:
                    beta_d = float(detector[f"beta_{mode}_m"])
                    phi_d = float(detector[f"phi_{mode}_turn"])
                    beta_product_sqrt = math.sqrt(beta_s * beta_d)
                    envelope = beta_product_sqrt / abs(denominator_signed)
                    response = beta_product_sqrt / denominator_signed * math.cos(
                        2.0 * math.pi * abs(phi_d - phi_s) - math.pi * tune
                    )
                    envelope_sq += envelope * envelope
                    phase_sq += response * response
                factor[f"envelope_{plane}_l2_m"] = math.sqrt(envelope_sq)
                factor[f"phase_{plane}_l2_m"] = math.sqrt(phase_sq)
                result[f"transport_envelope_{plane}_l2_m"] = factor[f"envelope_{plane}_l2_m"]
                result[f"transport_phase_{plane}_l2_m"] = factor[f"phase_{plane}_l2_m"]
            rows.append(result)
            factors[(trial, name)] = factor
    return rows, factors


def optics_variation_summary(
    nominal_points_path: Path,
    direction_points_path: Path,
    nominal_metadata_path: Path,
    direction_tunes_path: Path,
    nominal_factors: dict[str, dict[str, float]],
    direction_factors: dict[tuple[int, str], dict[str, float]],
) -> list[dict[str, object]]:
    nominal_points = {
        (row["point_type"], int(row["element_index"])): row
        for row in read_rows(nominal_points_path)
    }
    direction_points = read_rows(direction_points_path)
    changes: dict[str, list[float]] = defaultdict(list)
    signed_changes: dict[str, list[float]] = defaultdict(list)
    for row in direction_points:
        nominal = nominal_points[(row["point_type"], int(row["element_index"]))]
        for mode in ("1", "2"):
            baseline = float(nominal[f"beta_{mode}_m"])
            change = float(row[f"beta_{mode}_m"]) / baseline - 1.0
            scope = "sextupole" if row["point_type"] == "sextupole_exit" else "detector"
            quantity = f"{scope}_beta_{mode}"
            changes[quantity].append(abs(change))
            signed_changes[quantity].append(change)
    for (trial, name), factor in direction_factors.items():
        nominal = nominal_factors[name]
        for plane in ("x", "y"):
            for kind in ("envelope", "phase"):
                key = f"{kind}_{plane}_l2_m"
                change = factor[key] / nominal[key] - 1.0
                quantity = f"transport_{kind}_{plane}_l2"
                changes[quantity].append(abs(change))
                signed_changes[quantity].append(change)
    with nominal_metadata_path.open("rb") as stream:
        nominal_metadata = tomllib.load(stream)
    tune_rows = read_rows(direction_tunes_path)
    for plane, mode in (("x", "1"), ("y", "2")):
        baseline = float(nominal_metadata[f"full_tune_{mode}_turn"])
        quantity = f"full_tune_{plane}"
        for row in tune_rows:
            change = float(row[f"full_tune_{mode}_turn"]) / baseline - 1.0
            changes[quantity].append(abs(change))
            signed_changes[quantity].append(change)
    output: list[dict[str, object]] = []
    for quantity in changes:
        output.append(
            {
                "quantity": quantity,
                "samples": len(changes[quantity]),
                "mean_signed_relative_change": statistics.fmean(signed_changes[quantity]),
                "median_abs_relative_change": quantile(changes[quantity], 0.50),
                "p90_abs_relative_change": quantile(changes[quantity], 0.90),
                "maximum_abs_relative_change": max(changes[quantity]),
            }
        )
    return output


def contribution_lookup(path: Path) -> dict[tuple[int, str], dict[str, float]]:
    result: dict[tuple[int, str], dict[str, float]] = {}
    for row in read_rows(path):
        key = (int(row["trial"]), row["element_name"].lower())
        result[key] = {
            "contribution_norm_m": float(row["contribution_norm_m"]),
            "projection_numerator": float(row["projection_numerator"]),
        }
    return result


def direction_norms(path: Path) -> dict[int, float]:
    return {int(row["trial"]): float(row["q_total_norm_m"]) for row in read_rows(path)}


def build_analysis_rows(
    orbit_path: Path,
    nominal_factors: dict[str, dict[str, float]],
    direction_factors: dict[tuple[int, str], dict[str, float]],
    contribution_paths: dict[str, Path],
    direction_paths: dict[str, Path],
) -> list[dict[str, object]]:
    orbits = read_rows(orbit_path)
    if len(orbits) != 7600:
        raise ValueError(f"Expected 7600 direction-element orbit rows, found {len(orbits)}")
    contributions = {plane: contribution_lookup(path) for plane, path in contribution_paths.items()}
    norms = {plane: direction_norms(path) for plane, path in direction_paths.items()}
    rows: list[dict[str, object]] = []
    for orbit in orbits:
        trial = int(orbit["trial"])
        name = orbit["element_name"].lower()
        if name not in nominal_factors:
            raise ValueError(f"Orbit element missing from optics factors: {name}")
        nominal_factor = nominal_factors[name]
        direction_key = (trial, name)
        if direction_key not in direction_factors:
            raise ValueError(f"Direction optics missing trial={trial}, element={name}")
        direction_factor = direction_factors[direction_key]
        k2l = float(orbit["k2l_m2"])
        if not math.isclose(k2l, nominal_factor["k2l_m2"], rel_tol=1e-11, abs_tol=1e-14):
            raise ValueError(f"K2L mismatch for {name}: {k2l} versus {nominal_factor['k2l_m2']}")
        x_h = float(orbit["x_h_m"])
        y_v = float(orbit["y_v_m"])
        source_kicks = {
            "x": abs(-0.5 * k2l * (x_h * x_h - y_v * y_v)),
            "y": abs(k2l * x_h * y_v),
        }
        for plane in ("x", "y"):
            key = (trial, name)
            if key not in contributions[plane]:
                raise ValueError(f"Missing {plane} contribution for trial={trial}, element={name}")
            actual_norm = contributions[plane][key]["contribution_norm_m"]
            q_norm = norms[plane][trial]
            source = source_kicks[plane]
            rows.append(
                {
                    "plane": plane,
                    "trial": trial,
                    "element_order": int(orbit["element_order"]),
                    "element_index": int(orbit["element_index"]),
                    "element_name": orbit["element_name"],
                    "s_m": float(orbit["s_m"]),
                    "k2l_m2": k2l,
                    "x_h_m": x_h,
                    "y_v_m": y_v,
                    "nominal_beta_sext_m": nominal_factor[f"beta_{plane}_sext_m"],
                    "direction_beta_sext_m": direction_factor[f"beta_{plane}_sext_m"],
                    "local_source_kick_rad": source,
                    "nominal_transport_envelope_l2_m": nominal_factor[f"envelope_{plane}_l2_m"],
                    "nominal_transport_phase_l2_m": nominal_factor[f"phase_{plane}_l2_m"],
                    "direction_transport_envelope_l2_m": direction_factor[f"envelope_{plane}_l2_m"],
                    "direction_transport_phase_l2_m": direction_factor[f"phase_{plane}_l2_m"],
                    "predictor_source_only_rad": source,
                    "predictor_nominal_beta_envelope_m": source * nominal_factor[f"envelope_{plane}_l2_m"],
                    "predictor_direction_beta_envelope_m": source * direction_factor[f"envelope_{plane}_l2_m"],
                    "predictor_nominal_beta_phase_m": source * nominal_factor[f"phase_{plane}_l2_m"],
                    "predictor_direction_beta_phase_m": source * direction_factor[f"phase_{plane}_l2_m"],
                    "actual_contribution_norm_m": actual_norm,
                    "q_total_norm_m": q_norm,
                    "actual_relative_magnitude": actual_norm / q_norm,
                    "projection_numerator": contributions[plane][key]["projection_numerator"],
                }
            )

    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["plane"]), int(row["trial"]))].append(row)
    predictor_columns = tuple(column for _, column in PREDICTORS)
    for group in grouped.values():
        actual_sum = sum(float(row["actual_contribution_norm_m"]) for row in group)
        predictor_sums = {
            column: sum(float(row[column]) for row in group) for column in predictor_columns
        }
        for row in group:
            row["actual_magnitude_share"] = float(row["actual_contribution_norm_m"]) / actual_sum
            for column in predictor_columns:
                row[f"{column}_share"] = float(row[column]) / predictor_sums[column]
    return rows


def aggregate_elements(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["plane"]), str(row["element_name"]).lower())].append(row)
    output: list[dict[str, object]] = []
    predictors = tuple(column for _, column in PREDICTORS)
    for (plane, _), group in grouped.items():
        denominator = sum(float(row["q_total_norm_m"]) ** 2 for row in group)
        result: dict[str, object] = {
            "plane": plane,
            "element_order": group[0]["element_order"],
            "element_index": group[0]["element_index"],
            "element_name": group[0]["element_name"],
            "s_m": group[0]["s_m"],
            "k2l_m2": group[0]["k2l_m2"],
            "nominal_beta_sext_m": group[0]["nominal_beta_sext_m"],
            "actual_rms_relative_magnitude": math.sqrt(
                sum(float(row["actual_contribution_norm_m"]) ** 2 for row in group) / denominator
            ),
        }
        for predictor in predictors:
            result[f"{predictor}_rms"] = math.sqrt(
                statistics.fmean(float(row[predictor]) ** 2 for row in group)
            )
        output.append(result)
    for plane in ("x", "y"):
        plane_rows = [row for row in output if row["plane"] == plane]
        actual_sum = sum(float(row["actual_rms_relative_magnitude"]) for row in plane_rows)
        predictor_sums = {
            predictor: sum(float(row[f"{predictor}_rms"]) for row in plane_rows)
            for predictor in predictors
        }
        for row in plane_rows:
            row["actual_element_magnitude_share"] = float(row["actual_rms_relative_magnitude"]) / actual_sum
            for predictor in predictors:
                row[f"{predictor}_element_share"] = (
                    float(row[f"{predictor}_rms"]) / predictor_sums[predictor]
                )
    return sorted(output, key=lambda row: (str(row["plane"]), int(row["element_order"])))


def correlation_summary(
    rows: list[dict[str, object]],
    elements: list[dict[str, object]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    predictors = PREDICTORS
    for plane in ("x", "y"):
        plane_rows = [row for row in rows if row["plane"] == plane]
        plane_elements = [row for row in elements if row["plane"] == plane]
        trials = sorted({int(row["trial"]) for row in plane_rows})
        for label, predictor in predictors:
            actual = [float(row["actual_magnitude_share"]) for row in plane_rows]
            predicted = [float(row[f"{predictor}_share"]) for row in plane_rows]
            direction_spearman: list[float] = []
            direction_log_pearson: list[float] = []
            for trial in trials:
                trial_rows = [row for row in plane_rows if int(row["trial"]) == trial]
                trial_actual = [float(row["actual_magnitude_share"]) for row in trial_rows]
                trial_predicted = [float(row[f"{predictor}_share"]) for row in trial_rows]
                direction_spearman.append(spearman(trial_actual, trial_predicted))
                direction_log_pearson.append(log_pearson(trial_actual, trial_predicted))
            element_actual = [float(row["actual_element_magnitude_share"]) for row in plane_elements]
            element_predicted = [float(row[f"{predictor}_element_share"]) for row in plane_elements]
            output.append(
                {
                    "plane": plane,
                    "predictor": label,
                    "direction_element_rows": len(plane_rows),
                    "directions": len(trials),
                    "elements": len(plane_elements),
                    "pooled_share_spearman": spearman(actual, predicted),
                    "pooled_share_log_pearson": log_pearson(actual, predicted),
                    "direction_spearman_p10": quantile(direction_spearman, 0.10),
                    "direction_spearman_median": quantile(direction_spearman, 0.50),
                    "direction_spearman_p90": quantile(direction_spearman, 0.90),
                    "direction_log_pearson_p10": quantile(direction_log_pearson, 0.10),
                    "direction_log_pearson_median": quantile(direction_log_pearson, 0.50),
                    "direction_log_pearson_p90": quantile(direction_log_pearson, 0.90),
                    "element_share_spearman": spearman(element_actual, element_predicted),
                    "element_share_log_pearson": log_pearson(element_actual, element_predicted),
                }
            )
    return output


def svg_escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def render_svg(path: Path, plane: str, elements: list[dict[str, object]], summary: list[dict[str, object]]) -> None:
    plane_rows = [row for row in elements if row["plane"] == plane]
    predictors = (
        ("source_only", "Source only", "predictor_source_only_rad_element_share"),
        ("direction_beta_envelope", "Source × direction β", "predictor_direction_beta_envelope_m_element_share"),
        ("direction_beta_phase", "Source × direction β/phase", "predictor_direction_beta_phase_m_element_share"),
    )
    width, height = 1120, 410
    panel_width, panel_height = 300, 280
    top, gap = 72, 55
    left0 = 75
    colors = {"x": "#2563eb", "y": "#d946ef"}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.axis{stroke:#334155;stroke-width:1}.grid{stroke:#dbe3ec;stroke-width:1}.point{fill-opacity:.68}</style>',
        f'<text x="{width/2}" y="31" text-anchor="middle" font-size="21" font-weight="600">Detector-{plane}: sextupole contribution versus optics predictor</text>',
    ]
    for panel, (label, title, column) in enumerate(predictors):
        left = left0 + panel * (panel_width + gap)
        actual = [float(row["actual_element_magnitude_share"]) for row in plane_rows]
        predicted = [float(row[column]) for row in plane_rows]
        logs = [math.log10(value) for value in actual + predicted if value > 0]
        low = math.floor(min(logs))
        high = math.ceil(max(logs))
        if high <= low:
            high = low + 1

        def coordinate(value: float, horizontal: bool) -> float:
            fraction = (math.log10(value) - low) / (high - low)
            return left + fraction * panel_width if horizontal else top + (1.0 - fraction) * panel_height

        for tick in range(low, high + 1):
            x = coordinate(10.0**tick, True)
            y = coordinate(10.0**tick, False)
            lines.append(f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+panel_height}"/>')
            lines.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left+panel_width}" y2="{y:.2f}"/>')
            lines.append(f'<text x="{x:.2f}" y="{top+panel_height+19}" text-anchor="middle" font-size="11">10^{tick}</text>')
            lines.append(f'<text x="{left-8}" y="{y+4:.2f}" text-anchor="end" font-size="11">10^{tick}</text>')
        lines.extend(
            [
                f'<line class="axis" x1="{left}" y1="{top+panel_height}" x2="{left+panel_width}" y2="{top+panel_height}"/>',
                f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+panel_height}"/>',
                f'<line x1="{left}" y1="{top+panel_height}" x2="{left+panel_width}" y2="{top}" stroke="#94a3b8" stroke-dasharray="4 4"/>',
                f'<text x="{left+panel_width/2}" y="{top-23}" text-anchor="middle" font-size="15" font-weight="600">{svg_escape(title)}</text>',
                f'<text x="{left+panel_width/2}" y="{height-12}" text-anchor="middle" font-size="12">Predicted element magnitude share</text>',
            ]
        )
        stat = next(row for row in summary if row["plane"] == plane and row["predictor"] == label)
        lines.append(
            f'<text x="{left+8}" y="{top+18}" font-size="12">Spearman ρ = {float(stat["element_share_spearman"]):.3f}</text>'
        )
        for row in plane_rows:
            value_x = float(row[column])
            value_y = float(row["actual_element_magnitude_share"])
            lines.append(
                f'<circle class="point" cx="{coordinate(value_x, True):.2f}" cy="{coordinate(value_y, False):.2f}" r="3.2" fill="{colors[plane]}"><title>{svg_escape(str(row["element_name"]))}</title></circle>'
            )
    lines.append(
        f'<text transform="translate(18 {top+panel_height/2}) rotate(-90)" text-anchor="middle" font-size="12">Actual element magnitude share</text>'
    )
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_report(
    path: Path,
    summary: list[dict[str, object]],
    variation: list[dict[str, object]],
    inputs: argparse.Namespace,
) -> None:
    lookup = {(str(row["plane"]), str(row["predictor"])): row for row in summary}
    lines = [
        "# Sextupole contribution correlation with beta and phase",
        "",
        "## Scope",
        "",
        "This experiment reuses the 100 fixed corrector-direction pairs, the saved first-order sextupole orbit exposure, and the maintained thick-element Hessian contribution norms. RF-on Twiss functions and sextupole-to-detector phase advances are calculated both at the nominal state and separately at every simultaneous h+v direction state.",
        "",
        "The saved thick-element contribution is a 99-detector vector norm for each sextupole and direction. Therefore the beta/phase predictors are also assembled as 99-detector vectors and reduced with the same Euclidean norm. No existing contribution is recomputed.",
        "",
        "## Predictors",
        "",
        "- `source_only`: the dominant-plane thin-sextupole local source, `|K2L (x_h^2-y_v^2)/2|` for detector-x and `|K2L x_h y_v|` for detector-y.",
        "- `nominal_beta_envelope` and `direction_beta_envelope`: the source multiplied by the L2 norm over detectors of `sqrt(beta_i beta_j)/(2 |sin(pi Q)|)`, using nominal or direction-matched optics.",
        "- `nominal_beta_phase` and `direction_beta_phase`: the corresponding envelope including `cos(2 pi |phi_j-phi_i| - pi Q)` before the detector-vector norm.",
        "",
        "The mode-1 Twiss functions are used as the x-like predictor and mode-2 as the y-like predictor. These are uncoupled-style proxies; the exact thick contribution retains the full coupled six-dimensional transport.",
        "",
        "## Correlations",
        "",
        "| plane | predictor | pooled Spearman | pooled log Pearson | direction Spearman median [P10, P90] | element Spearman | element log Pearson |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for plane in ("x", "y"):
        for predictor, _ in PREDICTORS:
            row = lookup[(plane, predictor)]
            lines.append(
                f"| {plane} | `{predictor}` | {float(row['pooled_share_spearman']):.4f} | "
                f"{float(row['pooled_share_log_pearson']):.4f} | "
                f"{float(row['direction_spearman_median']):.4f} "
                f"[{float(row['direction_spearman_p10']):.4f}, {float(row['direction_spearman_p90']):.4f}] | "
                f"{float(row['element_share_spearman']):.4f} | "
                f"{float(row['element_share_log_pearson']):.4f} |"
            )
    variation_lookup = {str(row["quantity"]): row for row in variation}
    lines.extend(
        [
            "",
            "## Direction-optics variation from nominal",
            "",
            "| quantity | median absolute relative change | P90 | maximum |",
            "|---|---:|---:|---:|",
        ]
    )
    for quantity in (
        "sextupole_beta_1",
        "detector_beta_1",
        "sextupole_beta_2",
        "detector_beta_2",
        "transport_envelope_x_l2",
        "transport_phase_x_l2",
        "transport_envelope_y_l2",
        "transport_phase_y_l2",
    ):
        row = variation_lookup[quantity]
        lines.append(
            f"| `{quantity}` | {float(row['median_abs_relative_change']):.3e} | "
            f"{float(row['p90_abs_relative_change']):.3e} | "
            f"{float(row['maximum_abs_relative_change']):.3e} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "A correlation increase from `source_only` to a beta-envelope predictor measures the additional ranking information supplied by the beta functions. A further increase for the phase-aware predictor measures the value of phase advance for predicting the detector-vector magnitude. Comparing nominal and direction-matched rows tests whether orbit-dependent feed-down optics materially changes that relationship. This is still an association study, not a controlled beta-beating scan.",
            "",
            "The orbit file retains the dominant `x_h` and `y_v` local responses but not the smaller cross-plane `x_v` and `y_h` responses. The proxy therefore does not attempt to reproduce the exact coupled local Hessian source. Residual disagreement may come from those cross-plane terms, solenoidal coupling, finite-length sourcing, and non-sextupole complete-element sources.",
            "",
            "![Detector-x correlation](sextupole_beta_phase_correlation_x.svg)",
            "",
            "![Detector-y correlation](sextupole_beta_phase_correlation_y.svg)",
            "",
            "## Reused inputs",
            "",
            f"- Orbit exposure: `{inputs.orbits.resolve()}`",
            f"- Horizontal contribution: `{inputs.horizontal_contributions.resolve()}`",
            f"- Vertical contribution: `{inputs.vertical_contributions.resolve()}`",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = arguments()
    pair_rows, nominal_factors = optics_products(args.optics, args.optics_metadata)
    direction_factor_rows, direction_factors = direction_optics_products(
        args.direction_optics, args.direction_tunes
    )
    variation = optics_variation_summary(
        args.optics,
        args.direction_optics,
        args.optics_metadata,
        args.direction_tunes,
        nominal_factors,
        direction_factors,
    )
    rows = build_analysis_rows(
        args.orbits,
        nominal_factors,
        direction_factors,
        {"x": args.horizontal_contributions, "y": args.vertical_contributions},
        {"x": args.horizontal_directions, "y": args.vertical_directions},
    )
    elements = aggregate_elements(rows)
    summary = correlation_summary(rows, elements)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.output_dir / "sextupole_detector_optics.csv", pair_rows)
    write_rows(args.output_dir / "direction_sextupole_transport_factors.csv", direction_factor_rows)
    write_rows(args.output_dir / "direction_element_correlation_data.csv", rows)
    write_rows(args.output_dir / "element_correlation_data.csv", elements)
    write_rows(args.output_dir / "correlation_summary.csv", summary)
    write_rows(args.output_dir / "direction_optics_variation_summary.csv", variation)
    render_svg(args.output_dir / "sextupole_beta_phase_correlation_x.svg", "x", elements, summary)
    render_svg(args.output_dir / "sextupole_beta_phase_correlation_y.svg", "y", elements, summary)
    write_report(args.output_dir / "RESULTS.md", summary, variation, args)
    print(f"Results: {args.output_dir / 'RESULTS.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
