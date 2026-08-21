#!/usr/bin/env python3
"""Predict signed normal-sextupole detector contributions from source and optics.

The exact contribution file is produced by
``sextupole_detector_contributions/run_sextupole_detector_contributions.jl``.
This analysis deliberately uses the total detector vector only.  It compares
the per-element contribution magnitude/ranking with three physical
predictors: the local normal-sextupole source alone, source times a beta
Green-function envelope, and source times the phase-aware Green function.
"""

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
DEFAULT_RESULTS = HERE / "results" / "latest_cesr"
DEFAULT_CONTRIBUTIONS = (
    ERROR_ANALYSIS
    / "sextupole_detector_contributions"
    / "results"
    / "latest_cesr"
    / "sextupole_direction_contributions.csv"
)
DEFAULT_CLOSURE = (
    ERROR_ANALYSIS
    / "sextupole_detector_contributions"
    / "results"
    / "latest_cesr"
    / "direction_closure.csv"
)
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
    parser.add_argument("--optics-metadata", type=Path, default=DEFAULT_RESULTS / "nominal_optics_metadata.toml")
    parser.add_argument("--direction-optics", type=Path, default=DEFAULT_RESULTS / "direction_optics_points.csv")
    parser.add_argument("--direction-tunes", type=Path, default=DEFAULT_RESULTS / "direction_optics_tunes.csv")
    parser.add_argument("--contributions", type=Path, default=DEFAULT_CONTRIBUTIONS)
    parser.add_argument("--closure", type=Path, default=DEFAULT_CLOSURE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS)
    # Retain the old spellings as compatibility aliases.  New latest runs use
    # the one combined total-vector contribution file above.
    parser.add_argument("--horizontal-contributions", type=Path, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--vertical-contributions", type=Path, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--horizontal-directions", type=Path, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--vertical-directions", type=Path, default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.contributions == DEFAULT_CONTRIBUTIONS and args.horizontal_contributions is not None:
        args.contributions = args.horizontal_contributions
    if args.closure == DEFAULT_CLOSURE and args.horizontal_directions is not None:
        args.closure = args.horizontal_directions
    return args


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
    return pearson([math.log10(a) for a, _ in pairs], [math.log10(b) for _, b in pairs])


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
    tunes = {"x": float(metadata["full_tune_1_turn"]), "y": float(metadata["full_tune_2_turn"])}
    sextupoles = [row for row in points if row["point_type"] == "sextupole_exit"]
    detectors = [row for row in points if row["point_type"] == "detector"]
    if not sextupoles or not detectors:
        raise ValueError(f"Optics table must contain active sextupoles and detectors; got {len(sextupoles)} and {len(detectors)}")

    pair_rows: list[dict[str, object]] = []
    sums: dict[str, dict[str, float]] = {
        row["element_name"].lower(): {
            "envelope_x_sq": 0.0, "phase_x_sq": 0.0,
            "envelope_y_sq": 0.0, "phase_y_sq": 0.0,
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
                output.update({
                    f"beta_{plane}_sext_m": beta_s,
                    f"beta_{plane}_detector_m": beta_d,
                    f"phi_{plane}_sext_turn": phi_s,
                    f"phi_{plane}_detector_turn": phi_d,
                    f"abs_phase_advance_{plane}_turn": phase_advance,
                    f"closed_orbit_envelope_{plane}_m": envelope,
                    f"closed_orbit_response_{plane}_m": response,
                })
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
    if not points_by_trial or set(points_by_trial) != set(tunes):
        raise ValueError("Direction optics and tune tables must contain the same trials")

    rows: list[dict[str, object]] = []
    factors: dict[tuple[int, str], dict[str, float]] = {}
    for trial in sorted(points_by_trial):
        points = points_by_trial[trial]
        sextupoles = [row for row in points if row["point_type"] == "sextupole_exit"]
        detectors = [row for row in points if row["point_type"] == "detector"]
        if not sextupoles or not detectors:
            raise ValueError(f"Trial {trial}: expected active sextupoles and detectors")
        for sextupole in sextupoles:
            name = sextupole["element_name"].lower()
            result: dict[str, object] = {
                "trial": trial,
                "element_index": int(sextupole["element_index"]),
                "element_name": sextupole["element_name"],
                "s_m": float(sextupole["s_m"]),
                "k2l_m2": float(sextupole["k2l_m2"]),
            }
            sums = {"envelope_x_sq": 0.0, "phase_x_sq": 0.0, "envelope_y_sq": 0.0, "phase_y_sq": 0.0}
            for detector in detectors:
                for plane, mode in (("x", "1"), ("y", "2")):
                    beta_s = float(sextupole[f"beta_{mode}_m"])
                    beta_d = float(detector[f"beta_{mode}_m"])
                    phi_s = float(sextupole[f"phi_{mode}_turn"])
                    phi_d = float(detector[f"phi_{mode}_turn"])
                    tune = tunes[trial][plane]
                    denominator_signed = 2.0 * math.sin(math.pi * tune)
                    if denominator_signed == 0.0:
                        raise ValueError(f"Trial {trial}: mode-{mode} tune lies on an integer")
                    root_beta = math.sqrt(beta_s * beta_d)
                    envelope = root_beta / abs(denominator_signed)
                    response = root_beta / denominator_signed * math.cos(
                        2.0 * math.pi * abs(phi_d - phi_s) - math.pi * tune
                    )
                    sums[f"envelope_{plane}_sq"] += envelope * envelope
                    sums[f"phase_{plane}_sq"] += response * response
            result.update({
                "beta_x_sext_m": float(sextupole["beta_1_m"]),
                "beta_y_sext_m": float(sextupole["beta_2_m"]),
                "envelope_x_l2_m": math.sqrt(sums["envelope_x_sq"]),
                "phase_x_l2_m": math.sqrt(sums["phase_x_sq"]),
                "envelope_y_l2_m": math.sqrt(sums["envelope_y_sq"]),
                "phase_y_l2_m": math.sqrt(sums["phase_y_sq"]),
            })
            rows.append(result)
            factors[(trial, name)] = result
    return rows, factors


def optics_variation_summary(
    nominal: dict[str, dict[str, float]],
    direction: dict[tuple[int, str], dict[str, float]],
) -> list[dict[str, object]]:
    quantities = (
        "beta_x_sext_m", "beta_y_sext_m", "envelope_x_l2_m", "phase_x_l2_m",
        "envelope_y_l2_m", "phase_y_l2_m",
    )
    output: list[dict[str, object]] = []
    for quantity in quantities:
        changes: list[float] = []
        for (trial, name), row in direction.items():
            if name not in nominal:
                raise ValueError(f"Direction optics element {name} is absent from nominal optics")
            baseline = nominal[name][quantity]
            value = float(row[quantity])
            changes.append(abs(value - baseline) / max(abs(baseline), 1.0e-30))
        output.append({
            "quantity": quantity,
            "direction_rows": len(changes),
            "median_abs_relative_change": quantile(changes, 0.50),
            "p90_abs_relative_change": quantile(changes, 0.90),
            "maximum_abs_relative_change": max(changes),
        })
    return output


def contribution_rows(
    path: Path,
    closure_path: Path,
) -> tuple[list[dict[str, str]], dict[int, dict[str, float]]]:
    rows = read_rows(path)
    if not rows:
        raise ValueError(f"Contribution table is empty: {path}")
    required = {
        "trial", "element_name", "k2l_m2", "source_kick_x_rad", "source_kick_y_rad",
        "x_contribution_norm_m", "y_contribution_norm_m", "total_contribution_norm_m",
    }
    missing = sorted(required.difference(rows[0]))
    if missing:
        raise ValueError(f"Latest contribution table is missing columns: {', '.join(missing)}")
    closure_rows = read_rows(closure_path)
    closure: dict[int, dict[str, float]] = {}
    for row in closure_rows:
        trial = int(row["trial"])
        closure[trial] = {
            "x": float(row["q_x_norm_m"]),
            "y": float(row["q_y_norm_m"]),
            "total": float(row["q_total_norm_m"]),
        }
    if not closure:
        raise ValueError(f"Closure table is empty: {closure_path}")
    for row in rows:
        trial = int(row["trial"])
        if trial not in closure:
            raise ValueError(f"Contribution row trial={trial} has no closure row")
        for key in required - {"trial", "element_name"}:
            if not math.isfinite(float(row[key])):
                raise ValueError(f"Non-finite contribution value in {key}, trial={trial}")
    return rows, closure


def build_analysis_rows(
    contributions: list[dict[str, str]],
    closure: dict[int, dict[str, float]],
    nominal_factors: dict[str, dict[str, float]],
    direction_factors: dict[tuple[int, str], dict[str, float]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for contribution in contributions:
        trial = int(contribution["trial"])
        name = contribution["element_name"].lower()
        if name not in nominal_factors:
            raise ValueError(f"Contribution element {name} is absent from nominal optics")
        direction_key = (trial, name)
        if direction_key not in direction_factors:
            raise ValueError(f"Direction optics is missing trial={trial}, element={name}")
        nominal = nominal_factors[name]
        direction = direction_factors[direction_key]
        k2l = float(contribution["k2l_m2"])
        if not math.isclose(k2l, nominal["k2l_m2"], rel_tol=1e-11, abs_tol=1e-14):
            raise ValueError(f"K2L mismatch for {name}: {k2l} versus {nominal['k2l_m2']}")
        for plane in ("x", "y"):
            source_signed = float(contribution[f"source_kick_{plane}_rad"])
            source = abs(source_signed)
            actual = float(contribution[f"{plane}_contribution_norm_m"])
            q_norm = closure[trial][plane]
            if q_norm <= 0.0:
                raise ValueError(f"Non-positive Q norm for trial={trial}, plane={plane}")
            rows = {
                "plane": plane,
                "trial": trial,
                "element_order": int(contribution["element_order"]),
                "element_index": int(contribution["element_index"]),
                "element_name": contribution["element_name"],
                "s_m": float(contribution["s_m"]),
                "k2l_m2": k2l,
                "source_kick_signed_rad": source_signed,
                "local_source_kick_rad": source,
                "nominal_beta_sext_m": nominal[f"beta_{plane}_sext_m"],
                "direction_beta_sext_m": direction[f"beta_{plane}_sext_m"],
                "nominal_transport_envelope_l2_m": nominal[f"envelope_{plane}_l2_m"],
                "nominal_transport_phase_l2_m": nominal[f"phase_{plane}_l2_m"],
                "direction_transport_envelope_l2_m": direction[f"envelope_{plane}_l2_m"],
                "direction_transport_phase_l2_m": direction[f"phase_{plane}_l2_m"],
                "predictor_source_only_rad": source,
                "predictor_nominal_beta_envelope_m": source * nominal[f"envelope_{plane}_l2_m"],
                "predictor_direction_beta_envelope_m": source * direction[f"envelope_{plane}_l2_m"],
                "predictor_nominal_beta_phase_m": source * nominal[f"phase_{plane}_l2_m"],
                "predictor_direction_beta_phase_m": source * direction[f"phase_{plane}_l2_m"],
                "actual_contribution_norm_m": actual,
                "q_total_norm_m": q_norm,
                "actual_relative_magnitude": actual / q_norm,
            }
            output.append(rows)

    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for row in output:
        grouped[(str(row["plane"]), int(row["trial"]))].append(row)
    predictor_columns = tuple(column for _, column in PREDICTORS)
    for group in grouped.values():
        actual_sum = sum(float(row["actual_contribution_norm_m"]) for row in group)
        predictor_sums = {column: sum(float(row[column]) for row in group) for column in predictor_columns}
        if actual_sum <= 0.0 or any(value <= 0.0 for value in predictor_sums.values()):
            raise ValueError("Contribution/predictor sum is not positive")
        for row in group:
            row["actual_magnitude_share"] = float(row["actual_contribution_norm_m"]) / actual_sum
            for column in predictor_columns:
                row[f"{column}_share"] = float(row[column]) / predictor_sums[column]
    return output


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
                row[f"{predictor}_element_share"] = float(row[f"{predictor}_rms"]) / predictor_sums[predictor]
    return sorted(output, key=lambda row: (str(row["plane"]), int(row["element_order"])))


def correlation_summary(rows: list[dict[str, object]], elements: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for plane in ("x", "y"):
        plane_rows = [row for row in rows if row["plane"] == plane]
        plane_elements = [row for row in elements if row["plane"] == plane]
        trials = sorted({int(row["trial"]) for row in plane_rows})
        for label, predictor in PREDICTORS:
            actual = [float(row["actual_magnitude_share"]) for row in plane_rows]
            predicted = [float(row[f"{predictor}_share"]) for row in plane_rows]
            direction_spearman: list[float] = []
            direction_log_pearson: list[float] = []
            for trial in trials:
                trial_rows = [row for row in plane_rows if int(row["trial"]) == trial]
                direction_spearman.append(spearman(
                    [float(row["actual_magnitude_share"]) for row in trial_rows],
                    [float(row[f"{predictor}_share"]) for row in trial_rows],
                ))
                direction_log_pearson.append(log_pearson(
                    [float(row["actual_magnitude_share"]) for row in trial_rows],
                    [float(row[f"{predictor}_share"]) for row in trial_rows],
                ))
            element_actual = [float(row["actual_element_magnitude_share"]) for row in plane_elements]
            element_predicted = [float(row[f"{predictor}_element_share"]) for row in plane_elements]
            output.append({
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
            })
    return output


def svg_escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def render_svg(path: Path, plane: str, elements: list[dict[str, object]], summary: list[dict[str, object]]) -> None:
    plane_rows = [row for row in elements if row["plane"] == plane]
    predictors = (
        ("source_only", "Source only", "predictor_source_only_rad_element_share"),
        ("direction_beta_envelope", "Source × direction β Green envelope", "predictor_direction_beta_envelope_m_element_share"),
        ("direction_beta_phase", "Source × direction β/phase Green", "predictor_direction_beta_phase_m_element_share"),
    )
    width, height = 1120, 410
    panel_width, panel_height = 300, 280
    top, gap, left0 = 72, 55, 75
    colors = {"x": "#2563eb", "y": "#d946ef"}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#172033}.axis{stroke:#334155;stroke-width:1}.grid{stroke:#dbe3ec;stroke-width:1}.point{fill-opacity:.68}</style>',
        f'<text x="{width/2}" y="31" text-anchor="middle" font-size="21" font-weight="600">Total detector-{plane} contribution versus source/Green-function predictor</text>',
    ]
    for panel, (_, title, column) in enumerate(predictors):
        left = left0 + panel * (panel_width + gap)
        actual = [max(float(row["actual_element_magnitude_share"]), 1.0e-300) for row in plane_rows]
        predicted = [max(float(row[column]), 1.0e-300) for row in plane_rows]
        logs = [math.log10(value) for value in actual + predicted]
        low = math.floor(min(logs))
        high = max(low + 1, math.ceil(max(logs)))

        def coordinate(value: float, horizontal: bool) -> float:
            fraction = (math.log10(max(value, 1.0e-300)) - low) / (high - low)
            return left + fraction * panel_width if horizontal else top + (1.0 - fraction) * panel_height

        for tick in range(low, high + 1):
            x = coordinate(10.0**tick, True)
            y = coordinate(10.0**tick, False)
            lines.append(f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+panel_height}"/>')
            lines.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left+panel_width}" y2="{y:.2f}"/>')
            lines.append(f'<text x="{x:.2f}" y="{top+panel_height+19}" text-anchor="middle" font-size="11">10^{tick}</text>')
            lines.append(f'<text x="{left-8}" y="{y+4:.2f}" text-anchor="end" font-size="11">10^{tick}</text>')
        lines.extend([
            f'<line class="axis" x1="{left}" y1="{top+panel_height}" x2="{left+panel_width}" y2="{top+panel_height}"/>',
            f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+panel_height}"/>',
            f'<line x1="{left}" y1="{top+panel_height}" x2="{left+panel_width}" y2="{top}" stroke="#94a3b8" stroke-dasharray="4 4"/>',
            f'<text x="{left+panel_width/2}" y="{top-23}" text-anchor="middle" font-size="15" font-weight="600">{svg_escape(title)}</text>',
            f'<text x="{left+panel_width/2}" y="{height-12}" text-anchor="middle" font-size="12">Predicted element magnitude share</text>',
        ])
        stat = next(row for row in summary if row["plane"] == plane and row["predictor"] == _)
        lines.append(f'<text x="{left+8}" y="{top+18}" font-size="12">Spearman ρ = {float(stat["element_share_spearman"]):.3f}</text>')
        for row in plane_rows:
            value_x = float(row[column])
            value_y = float(row["actual_element_magnitude_share"])
            lines.append(
                f'<circle class="point" cx="{coordinate(value_x, True):.2f}" cy="{coordinate(value_y, False):.2f}" r="3.2" fill="{colors[plane]}"><title>{svg_escape(str(row["element_name"]))}</title></circle>'
            )
    lines.append(f'<text transform="translate(18 {top+panel_height/2}) rotate(-90)" text-anchor="middle" font-size="12">Actual element magnitude share</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_report(path: Path, summary: list[dict[str, object]], variation: list[dict[str, object]], args: argparse.Namespace) -> None:
    lookup = {(str(row["plane"]), str(row["predictor"])): row for row in summary}
    lines = [
        "# Latest-CESR normal-sextupole source--beta--phase predictor",
        "",
        "## Scope",
        "",
        "The target is the total second-order nonlinear detector vector `(Qx, Qy)`. Exact complete-element SciBmad Hessian sources are used for the signed normal-sextupole contribution; no hh/hv/vv block shares or third-order terms are reported.",
        "",
        "The predictors use the same sextupole source kick and detector set as the contribution table. They are physical ranking predictors, not replacements for the exact coupled six-dimensional source transport.",
        "",
        "## Green-function predictors",
        "",
        "For each transverse plane, the uncoupled reference Green function is `G_ij = sqrt(beta_i beta_j)/(2 sin(pi Q)) cos(2 pi |phi_i-phi_j| - pi Q)`. The envelope predictor uses `|source| ||G_env||_2`; the phase predictor uses `|source| ||G_ij||_2` over the configured detector registry. The source-only predictor is the absolute signed local normal-sextupole source kick reconstructed from the direction-matched first-order orbit; that kick already includes `K2L`.",
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
    lines.extend(["", "## Direction-optics variation from nominal", "", "| quantity | median absolute relative change | P90 | maximum |", "|---|---:|---:|---:|"])
    for quantity in variation_lookup:
        row = variation_lookup[quantity]
        lines.append(
            f"| `{quantity}` | {float(row['median_abs_relative_change']):.3e} | "
            f"{float(row['p90_abs_relative_change']):.3e} | {float(row['maximum_abs_relative_change']):.3e} |"
        )
    lines.extend([
        "",
        "## Interpretation boundary",
        "",
        "A correlation increase from source-only to the beta envelope measures the ranking information supplied by the beta response amplitude. A further increase for the phase-aware Green function measures the value of phase advance. Direction-matched optics retain the orbit-dependent operating point, while nominal optics are a fixed reference. This is an association/predictor study, not a controlled beta-beating scan.",
        "",
        "The exact contribution retains coupled six-dimensional transport and finite-element source terms. The predictor uses same-plane uncoupled Twiss quantities, so residual disagreement is expected from cross-plane coupling, finite-length sourcing, and non-sextupole sources.",
        "",
        f"- Contributions: `{args.contributions.resolve()}`",
        f"- Closure: `{args.closure.resolve()}`",
        f"- Nominal optics: `{args.optics.resolve()}`",
        f"- Direction optics: `{args.direction_optics.resolve()}`",
        "",
        "![Detector-x correlation](sextupole_beta_phase_correlation_x.svg)",
        "",
        "![Detector-y correlation](sextupole_beta_phase_correlation_y.svg)",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_metadata(path: Path, label: str) -> dict[str, object]:
    if not path.is_file():
        raise ValueError(f"Missing {label} metadata: {path}")
    with path.open("rb") as stream:
        return tomllib.load(stream)


def _normalized_metadata_path(value: object) -> str:
    return str(value).replace("\\", "/").rstrip("/").casefold()


def _consistent_metadata_field(
    sources: list[tuple[str, dict[str, object]]],
    field: str,
    *,
    normalizer=lambda value: value,
) -> object:
    values: list[tuple[str, object]] = []
    for label, metadata in sources:
        if field not in metadata:
            raise ValueError(f"{label} metadata is missing required field {field}")
        values.append((label, metadata[field]))
    reference = normalizer(values[0][1])
    mismatches = [
        f"{label}={value!r}"
        for label, value in values[1:]
        if normalizer(value) != reference
    ]
    if mismatches:
        details = ", ".join([f"{values[0][0]}={values[0][1]!r}"] + mismatches)
        raise ValueError(f"Metadata mismatch for {field}: {details}")
    return values[0][1]


def _optional_consistent_metadata_field(
    sources: list[tuple[str, dict[str, object]]],
    field: str,
    *,
    required_labels: tuple[str, ...],
    normalizer=lambda value: value,
) -> object:
    values = [(label, metadata[field]) for label, metadata in sources if field in metadata]
    present_labels = {label for label, _ in values}
    missing = [label for label in required_labels if label not in present_labels]
    if missing:
        raise ValueError(f"Metadata is missing required {field} in: {', '.join(missing)}")
    reference = normalizer(values[0][1])
    mismatches = [
        f"{label}={value!r}"
        for label, value in values[1:]
        if normalizer(value) != reference
    ]
    if mismatches:
        details = ", ".join([f"{values[0][0]}={values[0][1]!r}"] + mismatches)
        raise ValueError(f"Metadata mismatch for {field}: {details}")
    return values[0][1]


def _metadata_list(
    label: str,
    metadata: dict[str, object],
    *fields: str,
) -> list[str]:
    present = [metadata[field] for field in fields if field in metadata]
    if not present:
        raise ValueError(f"{label} metadata is missing one of: {', '.join(fields)}")
    values = [[str(item) for item in value] for value in present]
    if any(value != values[0] for value in values[1:]):
        raise ValueError(f"Metadata aliases disagree in {label}: {', '.join(fields)}")
    return values[0]


def _equal_metadata_lists(
    named_lists: list[tuple[str, list[str]]],
    field: str,
) -> list[str]:
    reference = named_lists[0][1]
    mismatches = [label for label, values in named_lists[1:] if values != reference]
    if mismatches:
        raise ValueError(f"Metadata mismatch for {field}: {', '.join(mismatches)}")
    return reference


def _check_metadata_count(
    sources: list[tuple[str, dict[str, object]]],
    fields: tuple[str, ...],
    expected: int,
    label: str,
    *,
    required: bool = True,
) -> None:
    for source_label, metadata in sources:
        values = [int(metadata[field]) for field in fields if field in metadata]
        if not values:
            if required:
                raise ValueError(f"{source_label} metadata is missing count for {label}")
            continue
        if any(value != values[0] for value in values[1:]) or values[0] != expected:
            raise ValueError(
                f"Metadata mismatch for {label}: {source_label} reports {values}, expected {expected}"
            )


def _toml_string(value: object) -> str:
    escaped = str(value).replace("\\", "/").replace('"', '\\"')
    return f'"{escaped}"'


def _toml_array(values: list[object]) -> str:
    return "[" + ", ".join(_toml_string(value) for value in values) + "]"


def write_metadata(path: Path, args: argparse.Namespace, contribution_rows: list[dict[str, str]], closure: dict[int, dict[str, float]], nominal: dict[str, dict[str, float]], direction: dict[tuple[int, str], dict[str, float]]) -> None:
    contribution_metadata_path = args.contributions.with_name("metadata.toml").resolve()
    nominal_metadata_path = args.optics_metadata.resolve()
    direction_metadata_path = args.direction_optics.with_name("direction_optics_metadata.toml").resolve()
    sources = [
        ("contribution", _read_metadata(contribution_metadata_path, "contribution")),
        ("nominal optics", _read_metadata(nominal_metadata_path, "nominal optics")),
        ("direction optics", _read_metadata(direction_metadata_path, "direction optics")),
    ]
    contribution_metadata = sources[0][1]
    nominal_metadata = sources[1][1]
    direction_metadata = sources[2][1]

    ring_id = str(_consistent_metadata_field(sources, "ring_id", normalizer=lambda value: str(value).casefold()))
    lattice_path_value = _consistent_metadata_field(sources, "lattice_path", normalizer=_normalized_metadata_path)
    lattice_path = str(lattice_path_value).replace("\\", "/")
    scibmad_version = str(_consistent_metadata_field(sources, "scibmad_version", normalizer=str))
    branch = int(_consistent_metadata_field(sources, "branch", normalizer=int))
    rf_on = bool(_consistent_metadata_field(sources, "rf_on", normalizer=bool))
    rf_voltage = float(
        _optional_consistent_metadata_field(
            sources,
            "rf_voltage",
            required_labels=("contribution", "nominal optics", "direction optics"),
            normalizer=float,
        )
    )
    state_dimension = int(_consistent_metadata_field(sources, "state_dimension", normalizer=int))

    control_lists = _equal_metadata_lists(
        [
            ("contribution", _metadata_list("contribution", contribution_metadata, "control_names", "control_names_from_config")),
            ("nominal optics", _metadata_list("nominal optics", nominal_metadata, "control_names", "control_names_from_config")),
            ("direction optics", _metadata_list("direction optics", direction_metadata, "control_names", "control_names_from_config")),
        ],
        "control_names",
    )
    detector_lists = _equal_metadata_lists(
        [(label, _metadata_list(label, metadata, "detector_names")) for label, metadata in sources],
        "detector_names",
    )
    sextupole_lists = _equal_metadata_lists(
        [
            ("contribution", _metadata_list("contribution", contribution_metadata, "active_normal_sextupole_names", "normal_sextupole_names")),
            ("nominal optics", _metadata_list("nominal optics", nominal_metadata, "active_normal_sextupole_names", "normal_sextupole_names")),
            ("direction optics", _metadata_list("direction optics", direction_metadata, "active_normal_sextupole_names", "normal_sextupole_names")),
        ],
        "normal_sextupole_names",
    )
    observable_lists = _equal_metadata_lists(
        [
            ("contribution", _metadata_list("contribution", contribution_metadata, "observable_labels", "observable_labels_from_config")),
            ("nominal optics", _metadata_list("nominal optics", nominal_metadata, "observable_labels", "observable_labels_from_config")),
            ("direction optics", _metadata_list("direction optics", direction_metadata, "observable_labels", "observable_labels_from_config")),
        ],
        "observable_labels",
    )

    control_count = len(control_lists)
    detector_count = len(detector_lists)
    active_normal_sextupoles = len(sextupole_lists)
    observable_count = len(observable_lists)
    _check_metadata_count(sources, ("control_count",), control_count, "control_count", required=False)
    _check_metadata_count(sources, ("detector_count", "detectors"), detector_count, "detector_count")
    _check_metadata_count(sources, ("active_normal_sextupoles", "sextupole_count"), active_normal_sextupoles, "active_normal_sextupoles")
    if observable_count != 2 * detector_count:
        raise ValueError(
            f"Metadata mismatch for observable_count: {observable_count} labels for {detector_count} detectors"
        )

    direction_count = len(closure)
    contribution_trials = {int(row["trial"]) for row in contribution_rows}
    direction_trials = {trial for trial, _ in direction}
    if contribution_trials != set(closure) or direction_trials != set(closure):
        raise ValueError("Contribution, closure, and direction-optics trial labels disagree")
    trials = int(_consistent_metadata_field([sources[0], sources[2]], "trials", normalizer=int))
    if trials != direction_count:
        raise ValueError(f"Metadata mismatch for trials: {trials} versus {direction_count} closure rows")
    seed = int(
        _optional_consistent_metadata_field(
            sources,
            "seed",
            required_labels=("contribution", "direction optics"),
            normalizer=int,
        )
    )
    base_kick_rad = float(
        _optional_consistent_metadata_field(
            sources,
            "base_kick_rad",
            required_labels=("contribution", "direction optics"),
            normalizer=float,
        )
    )

    contribution_path = str(args.contributions.resolve()).replace("\\", "/")
    closure_path = str(args.closure.resolve()).replace("\\", "/")
    metadata_paths = [
        str(contribution_metadata_path).replace("\\", "/"),
        str(nominal_metadata_path).replace("\\", "/"),
        str(direction_metadata_path).replace("\\", "/"),
    ]
    lines = [
        "format = \"cesr-latest-sextupole-source-beta-phase-predictor-v1\"",
        f"ring_id = {_toml_string(ring_id)}",
        f"lattice_path = {_toml_string(lattice_path)}",
        "engine = \"SciBmad\"",
        f"scibmad_version = {_toml_string(scibmad_version)}",
        f"branch = {branch}",
        f"rf_on = {'true' if rf_on else 'false'}",
        f"rf_voltage = {rf_voltage:.17g}",
        f"state_dimension = {state_dimension}",
        f"control_count = {control_count}",
        f"detector_count = {detector_count}",
        f"observable_count = {observable_count}",
        f"active_normal_sextupoles = {active_normal_sextupoles}",
        f"direction_count = {direction_count}",
        f"trials = {trials}",
        f"seed = {seed}",
        f"base_kick_rad = {base_kick_rad:.17g}",
        f"control_names = {_toml_array(control_lists)}",
        f"detector_names = {_toml_array(detector_lists)}",
        f"normal_sextupole_names = {_toml_array(sextupole_lists)}",
        f"observable_labels = {_toml_array(observable_lists)}",
        f"input_metadata_paths = {_toml_array(metadata_paths)}",
        f"contribution_metadata_toml = {_toml_string(metadata_paths[0])}",
        f"nominal_optics_metadata_toml = {_toml_string(metadata_paths[1])}",
        f"direction_optics_metadata_toml = {_toml_string(metadata_paths[2])}",
        f"contribution_rows = {len(contribution_rows)}",
        f"direction_optics_rows = {len(direction)}",
        "target = \"total nonlinear detector vector (x and y)\"",
        "source_boundary = \"exact complete-element exit boundary; predictor source is same-plane thin kick proxy\"",
        "phase_units = \"turn\"",
        "green_function = \"sqrt(beta_s beta_d)/(2 sin(pi Q)) times cosine phase factor\"",
        "predictors = [\"source_only\", \"nominal_beta_envelope\", \"direction_beta_envelope\", \"nominal_beta_phase\", \"direction_beta_phase\"]",
        f"contributions_csv = {_toml_string(contribution_path)}",
        f"closure_csv = {_toml_string(closure_path)}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = arguments()
    pair_rows, nominal = optics_products(args.optics, args.optics_metadata)
    direction_rows, direction = direction_optics_products(args.direction_optics, args.direction_tunes)
    variation = optics_variation_summary(nominal, direction)
    contributions, closure = contribution_rows(args.contributions, args.closure)
    rows = build_analysis_rows(contributions, closure, nominal, direction)
    elements = aggregate_elements(rows)
    summary = correlation_summary(rows, elements)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.output_dir / "sextupole_detector_optics.csv", pair_rows)
    write_rows(args.output_dir / "direction_sextupole_transport_factors.csv", direction_rows)
    write_rows(args.output_dir / "direction_element_correlation_data.csv", rows)
    write_rows(args.output_dir / "element_correlation_data.csv", elements)
    write_rows(args.output_dir / "correlation_summary.csv", summary)
    write_rows(args.output_dir / "direction_optics_variation_summary.csv", variation)
    render_svg(args.output_dir / "sextupole_beta_phase_correlation_x.svg", "x", elements, summary)
    render_svg(args.output_dir / "sextupole_beta_phase_correlation_y.svg", "y", elements, summary)
    write_report(args.output_dir / "RESULTS.md", summary, variation, args)
    write_metadata(args.output_dir / "metadata.toml", args, contributions, closure, nominal, direction)
    print(f"Results: {args.output_dir / 'RESULTS.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
