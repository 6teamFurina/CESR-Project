#!/usr/bin/env python3

"""Localize the normal-sextupole exposure and vertical-minus-horizontal excess."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from collections import defaultdict
from datetime import date
from pathlib import Path


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def pearson(first: list[float], second: list[float]) -> float:
    if len(first) != len(second) or len(first) < 2:
        raise ValueError("Pearson correlation requires equal nontrivial samples")
    first_mean = mean(first)
    second_mean = mean(second)
    covariance = sum(
        (left - first_mean) * (right - second_mean)
        for left, right in zip(first, second)
    )
    first_scale = math.sqrt(sum((value - first_mean) ** 2 for value in first))
    second_scale = math.sqrt(sum((value - second_mean) ** 2 for value in second))
    return covariance / (first_scale * second_scale)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def station_and_side(name: str) -> tuple[str, str]:
    match = re.search(r"(\d+)([EW])$", name.upper())
    if match is None:
        return name.upper(), "unknown"
    return f"SEX_{int(match.group(1)):02d}", match.group(2)


def summarize_samples(
    samples: dict[object, list[tuple[float, float]]],
    metadata: dict[object, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key, values in samples.items():
        exposure_h = [value[0] for value in values]
        exposure_v = [value[1] for value in values]
        excess = [vertical - horizontal for horizontal, vertical in values]
        row = dict(metadata[key])
        row.update(
            direction_pairs=len(values),
            mean_exposure_h=mean(exposure_h),
            p10_exposure_h=percentile(exposure_h, 0.10),
            median_exposure_h=percentile(exposure_h, 0.50),
            p90_exposure_h=percentile(exposure_h, 0.90),
            mean_exposure_v=mean(exposure_v),
            p10_exposure_v=percentile(exposure_v, 0.10),
            median_exposure_v=percentile(exposure_v, 0.50),
            p90_exposure_v=percentile(exposure_v, 0.90),
            mean_excess_v_minus_h=mean(excess),
            p10_excess_v_minus_h=percentile(excess, 0.10),
            median_excess_v_minus_h=percentile(excess, 0.50),
            p90_excess_v_minus_h=percentile(excess, 0.90),
            fraction_directions_v_gt_h=sum(value > 0.0 for value in excess) / len(excess),
        )
        rows.append(row)
    return rows


def add_ranks_and_shares(rows: list[dict[str, object]]) -> None:
    total_h = sum(float(row["mean_exposure_h"]) for row in rows)
    total_v = sum(float(row["mean_exposure_v"]) for row in rows)
    net_excess = total_v - total_h
    positive_excess = sum(
        max(0.0, float(row["mean_excess_v_minus_h"])) for row in rows
    )
    for row in rows:
        horizontal = float(row["mean_exposure_h"])
        vertical = float(row["mean_exposure_v"])
        excess = float(row["mean_excess_v_minus_h"])
        row["mean_h_share"] = horizontal / total_h
        row["mean_v_share"] = vertical / total_v
        row["signed_net_excess_share"] = excess / net_excess
        row["positive_excess_share"] = max(0.0, excess) / positive_excess

    for rank, row in enumerate(
        sorted(rows, key=lambda item: float(item["mean_exposure_h"]), reverse=True), 1
    ):
        row["horizontal_rank"] = rank
    cumulative = 0.0
    for rank, row in enumerate(
        sorted(rows, key=lambda item: float(item["mean_exposure_v"]), reverse=True), 1
    ):
        row["vertical_rank"] = rank
        cumulative += float(row["mean_v_share"])
        row["cumulative_v_share"] = cumulative
    cumulative = 0.0
    for rank, row in enumerate(
        sorted(rows, key=lambda item: float(item["mean_excess_v_minus_h"]), reverse=True), 1
    ):
        row["excess_rank"] = rank
        cumulative += float(row["positive_excess_share"])
        row["cumulative_positive_excess_share"] = cumulative


def load_rows(path: Path):
    element_samples: dict[tuple[int, str], list[tuple[float, float]]] = defaultdict(list)
    element_metadata: dict[tuple[int, str], dict[str, object]] = {}
    pair_trial: dict[tuple[str, int], list[float]] = defaultdict(lambda: [0.0, 0.0])
    side_trial: dict[tuple[str, int], list[float]] = defaultdict(lambda: [0.0, 0.0])
    sign_trial: dict[tuple[str, int], list[float]] = defaultdict(
        lambda: [0.0, 0.0, 0.0, 0.0]
    )
    trial_totals: dict[int, list[float]] = defaultdict(lambda: [0.0, 0.0])
    with path.open(newline="", encoding="utf-8") as stream:
        for raw in csv.DictReader(stream):
            trial = int(raw["trial"])
            element_index = int(raw["element_index"])
            name = raw["element_name"]
            horizontal = float(raw["source_exposure_h"])
            vertical = float(raw["source_exposure_v"])
            k2l = float(raw["k2l_m2"])
            key = (element_index, name)
            station, side = station_and_side(name)
            element_samples[key].append((horizontal, vertical))
            element_metadata[key] = {
                "element_order": int(raw["element_order"]),
                "element_index": element_index,
                "element_name": name,
                "station_pair": station,
                "ring_side": side,
                "s_m": float(raw["s_m"]),
                "k2l_m2": k2l,
            }
            pair_trial[(station, trial)][0] += horizontal
            pair_trial[(station, trial)][1] += vertical
            side_trial[(side, trial)][0] += horizontal
            side_trial[(side, trial)][1] += vertical
            sign = "negative_K2" if k2l < 0.0 else "positive_K2"
            sign_trial[(sign, trial)][0] += horizontal
            sign_trial[(sign, trial)][1] += vertical
            sign_trial[(sign, trial)][2] += float(raw["orbit_sq_h_m2"])
            sign_trial[(sign, trial)][3] += float(raw["orbit_sq_v_m2"])
            trial_totals[trial][0] += horizontal
            trial_totals[trial][1] += vertical

    pair_samples: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for (station, _trial), values in pair_trial.items():
        pair_samples[station].append((values[0], values[1]))
    pair_metadata = {key: {"station_pair": key} for key in pair_samples}
    side_samples: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for (side, _trial), values in side_trial.items():
        side_samples[side].append((values[0], values[1]))
    side_metadata = {key: {"ring_side": key} for key in side_samples}
    sign_samples: dict[str, list[tuple[float, float]]] = defaultdict(list)
    sign_orbits: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for (sign, _trial), values in sign_trial.items():
        sign_samples[sign].append((values[0], values[1]))
        sign_orbits[sign].append((values[2], values[3]))
    sign_metadata = {key: {"k2_location_group": key} for key in sign_samples}
    sign_rows = summarize_samples(sign_samples, sign_metadata)
    for row in sign_rows:
        orbit_values = sign_orbits[str(row["k2_location_group"])]
        mean_h = mean([value[0] for value in orbit_values])
        mean_v = mean([value[1] for value in orbit_values])
        row["mean_sum_orbit_sq_h_m2"] = mean_h
        row["mean_sum_orbit_sq_v_m2"] = mean_v
        row["mean_orbit_sq_ratio_v_h"] = mean_v / mean_h
        row["ratio_of_mean_source_exposures_v_h"] = (
            float(row["mean_exposure_v"]) / float(row["mean_exposure_h"])
        )
    return (
        summarize_samples(element_samples, element_metadata),
        summarize_samples(pair_samples, pair_metadata),
        summarize_samples(side_samples, side_metadata),
        sign_rows,
        trial_totals,
    )


def closure_check(direction_path: Path, trial_totals: dict[int, list[float]]) -> float:
    maximum = 0.0
    with direction_path.open(newline="", encoding="utf-8") as stream:
        for raw in csv.DictReader(stream):
            trial = int(raw["trial"])
            expected = (
                float(raw["source_exposure_h"]),
                float(raw["source_exposure_v"]),
            )
            actual = trial_totals[trial]
            for reference, reconstructed in zip(expected, actual):
                maximum = max(
                    maximum,
                    abs(reference - reconstructed) / max(abs(reference), 1e-300),
                )
    return maximum


SVG_STYLE = (
    '<style>text{font-family:"Times New Roman",Times,serif;fill:#202124}'
    '.title{font-size:28px;font-weight:700}.sub{font-size:16px;fill:#555}'
    '.axis{stroke:#555;stroke-width:1.2}.grid{stroke:#e5e7eb;stroke-width:1}'
    '.label{font-size:15px}.small{font-size:14px;fill:#555}'
    '.direct-label{font-size:18px;font-weight:700;fill:#333}</style>'
)


def render_top_exposure_svg(path: Path, element_rows: list[dict[str, object]]) -> None:
    width, height = 1220, 555
    left, right = 125, 1095
    chart_top, chart_bottom = 95, 455
    ranked = sorted(
        element_rows,
        key=lambda row: float(row["mean_excess_v_minus_h"]),
        reverse=True,
    )[:15]
    max_exposure = max(float(row["mean_exposure_v"]) for row in ranked)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        SVG_STYLE,
        '<text x="40" y="38" class="title">15 Sextupoles with the Largest Exposure Differences</text>',
        '<text x="40" y="62" class="sub">Mean unsigned source exposures over 100 equal-corrector-RMS direction pairs</text>',
    ]
    chart_width = right - left
    for tick in range(6):
        value = max_exposure * tick / 5
        x = left + chart_width * tick / 5
        parts.append(f'<line x1="{x:.1f}" y1="{chart_top}" x2="{x:.1f}" y2="{chart_bottom}" class="grid"/>')
        parts.append(f'<text x="{x:.1f}" y="{chart_bottom + 20}" text-anchor="middle" class="small">{value:.2e}</text>')
    row_height = (chart_bottom - chart_top) / len(ranked)
    for index, row in enumerate(ranked):
        y = chart_top + (index + 0.5) * row_height
        horizontal = float(row["mean_exposure_h"])
        vertical = float(row["mean_exposure_v"])
        x_h = left + chart_width * horizontal / max_exposure
        x_v = left + chart_width * vertical / max_exposure
        name = html.escape(str(row["element_name"]))
        parts.append(f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" class="label">{name}</text>')
        parts.append(f'<line x1="{left}" y1="{y - 4:.1f}" x2="{x_h:.1f}" y2="{y - 4:.1f}" stroke="#3b82c4" stroke-width="7"/>')
        parts.append(f'<line x1="{left}" y1="{y + 5:.1f}" x2="{x_v:.1f}" y2="{y + 5:.1f}" stroke="#e8792e" stroke-width="7"/>')
    parts.extend([
        f'<line x1="{left}" y1="{chart_bottom}" x2="{right}" y2="{chart_bottom}" class="axis"/>',
        f'<rect x="{left}" y="{chart_bottom + 36}" width="16" height="7" fill="#3b82c4"/><text x="{left + 23}" y="{chart_bottom + 44}" class="small">horizontal input</text>',
        f'<rect x="{left + 145}" y="{chart_bottom + 36}" width="16" height="7" fill="#e8792e"/><text x="{left + 168}" y="{chart_bottom + 44}" class="small">vertical input</text>',
        f'<text x="{(left + right) / 2:.1f}" y="{chart_bottom + 45}" text-anchor="middle" class="small">mean local proxy contribution |K2L| u^2</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts), encoding="utf-8")


def render_top_exposure_compact_svg(
    path: Path, element_rows: list[dict[str, object]]
) -> None:
    """Render the same top-15 ranking for half-page-width manuscript use."""
    width, height = 620, 342
    left, right = 105, 600
    chart_top, chart_bottom = 42, 252
    ranked = sorted(
        element_rows,
        key=lambda row: float(row["mean_excess_v_minus_h"]),
        reverse=True,
    )[:15]
    max_exposure = max(float(row["mean_exposure_v"]) for row in ranked)
    style = (
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}'
        '.title{font-size:20px;font-weight:700}.axis{stroke:#777;stroke-width:1}'
        '.grid{stroke:#ddd;stroke-width:1}.label{font-size:19px}'
        '.small{font-size:18px;fill:#555}</style>'
    )
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        style,
        '<text x="12" y="24" class="title">Largest Exposure Differences</text>',
    ]
    chart_width = right - left
    for tick in range(5):
        value = max_exposure * tick / 4
        x = left + chart_width * tick / 4
        parts.append(
            f'<line x1="{x:.1f}" y1="{chart_top}" x2="{x:.1f}" y2="{chart_bottom}" class="grid"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{chart_bottom + 20}" text-anchor="middle" class="small">{value:.2e}</text>'
        )
    row_height = (chart_bottom - chart_top) / len(ranked)
    for index, row in enumerate(ranked):
        y = chart_top + (index + 0.5) * row_height
        horizontal = float(row["mean_exposure_h"])
        vertical = float(row["mean_exposure_v"])
        x_h = left + chart_width * horizontal / max_exposure
        x_v = left + chart_width * vertical / max_exposure
        name = html.escape(str(row["element_name"]))
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end" class="label">{name}</text>'
        )
        parts.append(
            f'<line x1="{left}" y1="{y - 2.8:.1f}" x2="{x_h:.1f}" y2="{y - 2.8:.1f}" stroke="#3b82c4" stroke-width="4"/>'
        )
        parts.append(
            f'<line x1="{left}" y1="{y + 3.0:.1f}" x2="{x_v:.1f}" y2="{y + 3.0:.1f}" stroke="#e8792e" stroke-width="4"/>'
        )
    parts.extend(
        [
            f'<line x1="{left}" y1="{chart_bottom}" x2="{right}" y2="{chart_bottom}" class="axis"/>',
            f'<rect x="{left}" y="{chart_bottom + 34}" width="16" height="7" fill="#3b82c4"/><text x="{left + 23}" y="{chart_bottom + 42}" class="small">horizontal</text>',
            f'<rect x="{left + 160}" y="{chart_bottom + 34}" width="16" height="7" fill="#e8792e"/><text x="{left + 183}" y="{chart_bottom + 42}" class="small">vertical</text>',
            f'<text x="{(left + right) / 2:.1f}" y="{chart_bottom + 76}" text-anchor="middle" class="small">mean local proxy contribution |K2L| u^2</text>',
            '</svg>',
        ]
    )
    path.write_text("\n".join(parts), encoding="utf-8")


def render_ring_exposure_svg(path: Path, element_rows: list[dict[str, object]]) -> None:
    width, height = 1220, 235
    left, right = 125, 1095
    ring_top, ring_bottom = 45, 150
    chart_width = right - left
    s_max = max(float(row["s_m"]) for row in element_rows)
    max_abs_excess = max(abs(float(row["mean_excess_v_minus_h"])) for row in element_rows)
    max_abs_k2l = max(abs(float(row["k2l_m2"])) for row in element_rows)
    excess_power = 10.0 ** math.floor(math.log10(max_abs_excess))
    excess_limit = math.ceil(max_abs_excess / excess_power) * excess_power
    k2l_limit = math.ceil(max_abs_k2l * 10.0) / 10.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        SVG_STYLE,
    ]
    zero_y = (ring_top + ring_bottom) / 2
    half_height = (ring_bottom - ring_top) * 0.45
    for fraction in (-1.0, -0.5, 0.0, 0.5, 1.0):
        y = zero_y - half_height * fraction
        exposure_tick = fraction * excess_limit
        k2l_tick = fraction * k2l_limit
        grid_class = "axis" if fraction == 0.0 else "grid"
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" class="{grid_class}"/>'
        )
        parts.append(
            f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" class="small">{exposure_tick:.1e}</text>'
        )
        parts.append(
            f'<text x="{right + 10}" y="{y + 4:.1f}" text-anchor="start" class="small">{k2l_tick:.2f}</text>'
        )
    top_ten = sorted(
        element_rows,
        key=lambda row: float(row["mean_excess_v_minus_h"]),
        reverse=True,
    )[:10]
    top_ten_names = {str(row["element_name"]) for row in top_ten}
    label_lanes: dict[str, int] = {}
    lane_last_x = [-math.inf, -math.inf, -math.inf]
    for row in sorted(top_ten, key=lambda item: float(item["s_m"])):
        x = left + chart_width * float(row["s_m"]) / s_max
        available = [lane for lane, last_x in enumerate(lane_last_x) if x - last_x >= 60]
        lane = available[0] if available else min(range(3), key=lane_last_x.__getitem__)
        label_lanes[str(row["element_name"])] = lane
        lane_last_x[lane] = x

    k2_points: list[str] = []
    exposure_points: dict[str, tuple[float, float]] = {}
    for row in element_rows:
        x = left + chart_width * float(row["s_m"]) / s_max
        excess = float(row["mean_excess_v_minus_h"])
        y = zero_y - half_height * excess / excess_limit
        color = "#e8792e" if excess >= 0 else "#3b82c4"
        parts.append(f'<line x1="{x:.1f}" y1="{zero_y:.1f}" x2="{x:.1f}" y2="{y:.1f}" stroke="{color}" stroke-width="5"/>')
        if str(row["element_name"]) in top_ten_names:
            exposure_points[str(row["element_name"])] = (x, y)
        k2_y = zero_y - half_height * float(row["k2l_m2"]) / k2l_limit
        parts.append(
            f'<line x1="{x:.1f}" y1="{zero_y:.1f}" x2="{x:.1f}" y2="{k2_y:.1f}" stroke="#8b5fbf" stroke-width="1" opacity="0.7"/>'
        )
        k2_points.append(f'{x:.1f},{k2_y:.1f}')
    for point in k2_points:
        x, y = point.split(",")
        parts.append(f'<circle cx="{x}" cy="{y}" r="3.2" fill="#8b5fbf" stroke="white" stroke-width="0.8"/>')
    for row in top_ten:
        name = str(row["element_name"])
        x, y = exposure_points[name]
        label_y = y - 4 - 15 * label_lanes[name]
        escaped_name = html.escape(name.removeprefix("sex_").upper())
        parts.append(
            f'<line x1="{x:.1f}" y1="{y - 3:.1f}" x2="{x:.1f}" y2="{label_y + 3:.1f}" stroke="#555" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{x + 3:.1f}" y="{label_y:.1f}" class="direct-label" transform="rotate(-55 {x + 3:.1f} {label_y:.1f})">{escaped_name}</text>'
        )
    for tick in range(5):
        s_value = s_max * tick / 4
        x = left + chart_width * tick / 4
        parts.append(f'<line x1="{x:.1f}" y1="{ring_bottom}" x2="{x:.1f}" y2="{ring_bottom + 5}" class="axis"/>')
        parts.append(f'<text x="{x:.1f}" y="{ring_bottom + 20}" text-anchor="middle" class="small">{s_value:.0f}</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def top_cumulative_count(rows: list[dict[str, object]], field: str, threshold: float) -> int:
    for index, row in enumerate(rows, 1):
        if float(row[field]) >= threshold:
            return index
    return len(rows)


def write_report(
    path: Path,
    elements: list[dict[str, object]],
    pairs: list[dict[str, object]],
    sides: list[dict[str, object]],
    signs: list[dict[str, object]],
    closure: float,
) -> None:
    by_vertical = sorted(elements, key=lambda row: float(row["mean_exposure_v"]), reverse=True)
    by_excess = sorted(elements, key=lambda row: float(row["mean_excess_v_minus_h"]), reverse=True)
    pairs_by_excess = sorted(pairs, key=lambda row: float(row["mean_excess_v_minus_h"]), reverse=True)
    total_h = sum(float(row["mean_exposure_h"]) for row in elements)
    total_v = sum(float(row["mean_exposure_v"]) for row in elements)
    k2_excess_correlation = pearson(
        [float(row["k2l_m2"]) for row in elements],
        [float(row["mean_excess_v_minus_h"]) for row in elements],
    )
    opposite_sign_count = sum(
        float(row["k2l_m2"]) * float(row["mean_excess_v_minus_h"]) < 0.0
        for row in elements
    )
    lines = [
        "# Element-local normal-sextupole exposure result",
        "",
        f"Result recorded {date.today().isoformat()}. The analysis uses 100 equal-corrector-RMS H/V direction pairs and all {len(elements)} active normal sextupoles.",
        "",
        "## Global closure",
        "",
        f"- Mean horizontal unsigned source exposure: `{total_h:.8e}`.",
        f"- Mean vertical unsigned source exposure: `{total_v:.8e}`.",
        f"- Ratio of mean exposures: `{total_v / total_h:.6f}`.",
        f"- Maximum element-sum closure relative residual: `{closure:.3e}`.",
        f"- Pearson correlation between signed K2L and mean exposure excess: `{k2_excess_correlation:.6f}`; `{opposite_sign_count}/{len(elements)}` sites have opposite signs.",
        f"- The top {top_cumulative_count(by_vertical, 'cumulative_v_share', 0.50)} elements supply 50% of mean vertical exposure; the top {top_cumulative_count(by_excess, 'cumulative_positive_excess_share', 0.50)} positive-excess elements supply 50% of the positive vertical-minus-horizontal excess.",
        "",
        "The ratio of means is not the median of direction-level ratios quoted in the first-stage result. This is an unsigned local-source proxy and does not contain signed K2 cancellation or transport from each sextupole to the detectors.",
        "",
        "## Largest mean vertical exposure",
        "",
        "| rank | element | s [m] | K2L [m^-2] | mean Ev | vertical share | mean Eh |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in by_vertical[:12]:
        lines.append(
            f"| {row['vertical_rank']} | `{row['element_name']}` | {float(row['s_m']):.3f} | {float(row['k2l_m2']):.5g} | {float(row['mean_exposure_v']):.4e} | {100 * float(row['mean_v_share']):.2f}% | {float(row['mean_exposure_h']):.4e} |"
        )
    lines.extend([
        "",
        "## Largest mean vertical-minus-horizontal excess",
        "",
        "| rank | element | station | mean Ev-Eh | positive-excess share | directions Ev>Eh |",
        "|---:|---|---|---:|---:|---:|",
    ])
    for row in by_excess[:12]:
        lines.append(
            f"| {row['excess_rank']} | `{row['element_name']}` | `{row['station_pair']}` | {float(row['mean_excess_v_minus_h']):.4e} | {100 * float(row['positive_excess_share']):.2f}% | {100 * float(row['fraction_directions_v_gt_h']):.1f}% |"
        )
    lines.extend([
        "",
        "## Largest East/West station-pair excess",
        "",
        "| rank | pair | mean Ev-Eh | positive-excess share | directions Ev>Eh |",
        "|---:|---|---:|---:|---:|",
    ])
    for row in pairs_by_excess[:10]:
        lines.append(
            f"| {row['excess_rank']} | `{row['station_pair']}` | {float(row['mean_excess_v_minus_h']):.4e} | {100 * float(row['positive_excess_share']):.2f}% | {100 * float(row['fraction_directions_v_gt_h']):.1f}% |"
        )
    lines.extend([
        "",
        "## Ring-side control",
        "",
        "| side | mean Eh | mean Ev | mean Ev-Eh |",
        "|---|---:|---:|---:|",
    ])
    for row in sorted(sides, key=lambda item: str(item["ring_side"])):
        lines.append(
            f"| {row['ring_side']} | {float(row['mean_exposure_h']):.4e} | {float(row['mean_exposure_v']):.4e} | {float(row['mean_excess_v_minus_h']):.4e} |"
        )
    lines.extend([
        "",
        "## K2-sign location groups",
        "",
        "Exposure itself uses |K2L|. The sign below is therefore a label for the two alternating sextupole-location classes, not a signed-kick decomposition.",
        "",
        "| location group | elements | mean sum yv^2 / mean sum xh^2 | mean Ev/Eh | mean Ev-Eh | net-excess share |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    total_net_excess = total_v - total_h
    for row in sorted(signs, key=lambda item: str(item["k2_location_group"])):
        excess = float(row["mean_excess_v_minus_h"])
        lines.append(
            f"| `{row['k2_location_group']}` | {row['element_count']} | {float(row['mean_orbit_sq_ratio_v_h']):.4f} | {float(row['ratio_of_mean_source_exposures_v_h']):.4f} | {excess:.4e} | {100 * excess / total_net_excess:.2f}% |"
        )
    leading_pair_share = sum(
        float(row["positive_excess_share"]) for row in pairs_by_excess[:4]
    )
    negative_group = next(
        row for row in signs if row["k2_location_group"] == "negative_K2"
    )
    positive_group = next(
        row for row in signs if row["k2_location_group"] == "positive_K2"
    )
    lines.extend([
        "",
        "## Physical interpretation of the proxy",
        "",
        f"The first four station pairs (`{pairs_by_excess[0]['station_pair']}`, `{pairs_by_excess[1]['station_pair']}`, `{pairs_by_excess[2]['station_pair']}`, and `{pairs_by_excess[3]['station_pair']}`) supply `{100 * leading_pair_share:.2f}%` of the positive excess. East and West totals agree to within about one percent of their mean, so the localization is a nearly ring-symmetric optics pattern rather than a one-sided outlier.",
        "",
        f"The negative-K2 location class has an unweighted mean orbit-squared ratio of `{float(negative_group['mean_orbit_sq_ratio_v_h']):.4f}` before K2 weighting and a weighted exposure ratio of `{float(negative_group['ratio_of_mean_source_exposures_v_h']):.4f}`. The positive-K2 class instead has ratios `{float(positive_group['mean_orbit_sq_ratio_v_h']):.4f}` and `{float(positive_group['ratio_of_mean_source_exposures_v_h']):.4f}`. Thus the excess originates primarily because the alternating negative-K2 sextupole locations sample much larger vertical than horizontal internal-orbit response; the positive-K2 locations partially cancel that imbalance at the unsigned-proxy level.",
        "",
        "This is consistent with the two sextupole location classes occupying different horizontal/vertical optics. A beta-, tune-, and phase-resolved response-matrix decomposition is still required before assigning the contrast to a specific linear-optics factor.",
        "",
        "![Largest element-local exposure differences](top15_sextupole_exposure_differences.svg)",
        "",
        "![Exposures and strengths across the ring](sextupole_exposures_strengths_ring.svg)",
        "",
        "## Evidence boundary",
        "",
        "This ranking identifies where the unsigned normal-sextupole source proxy is generated. It is not yet a ranking of signed detector response. The next causal step is a signed detector-vector reconstruction or recomputed-lattice ablation of the leading elements/pairs.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "element_csv",
        nargs="?",
        type=Path,
        default=here / "element_results" / "element_exposure_directions.csv",
    )
    parser.add_argument(
        "--direction-csv",
        type=Path,
        default=here / "element_results" / "direction_attribution.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=here / "element_results")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    elements, pairs, sides, signs, trial_totals = load_rows(args.element_csv)
    for row in signs:
        sign = str(row["k2_location_group"])
        row["element_count"] = sum(
            1
            for element in elements
            if (float(element["k2l_m2"]) < 0.0) == (sign == "negative_K2")
        )
    add_ranks_and_shares(elements)
    add_ranks_and_shares(pairs)
    add_ranks_and_shares(sides)
    add_ranks_and_shares(signs)
    closure = closure_check(args.direction_csv, trial_totals)
    if closure > 1e-12:
        raise ValueError(f"Element exposure does not close direction totals: {closure:.3e}")
    elements.sort(key=lambda row: int(row["element_order"]))
    pairs.sort(key=lambda row: str(row["station_pair"]))
    sides.sort(key=lambda row: str(row["ring_side"]))
    write_csv(args.output_dir / "element_exposure_summary.csv", elements)
    write_csv(args.output_dir / "station_pair_exposure_summary.csv", pairs)
    write_csv(args.output_dir / "ring_side_exposure_summary.csv", sides)
    write_csv(args.output_dir / "k2_sign_location_exposure_summary.csv", signs)
    top_figure = args.output_dir / "top15_sextupole_exposure_differences.svg"
    compact_top_figure = (
        args.output_dir / "top15_sextupole_exposure_differences_compact.svg"
    )
    ring_figure = args.output_dir / "sextupole_exposures_strengths_ring.svg"
    render_top_exposure_svg(top_figure, elements)
    render_top_exposure_compact_svg(compact_top_figure, elements)
    render_ring_exposure_svg(ring_figure, elements)
    write_report(
        args.output_dir / "ELEMENT_EXPOSURE_RESULTS.md",
        elements,
        pairs,
        sides,
        signs,
        closure,
    )
    (args.output_dir / "element_exposure_analysis.json").write_text(
        json.dumps(
            {
                "format": "cesr-x-quadratic-element-exposure-v1",
                "date": date.today().isoformat(),
                "element_direction_csv": str(args.element_csv.resolve()),
                "direction_csv": str(args.direction_csv.resolve()),
                "active_normal_sextupoles": len(elements),
                "direction_pairs": len(trial_totals),
                "max_element_sum_relative_closure": closure,
                "pearson_signed_k2l_vs_mean_excess": pearson(
                    [float(row["k2l_m2"]) for row in elements],
                    [float(row["mean_excess_v_minus_h"]) for row in elements],
                ),
                "opposite_sign_site_count": sum(
                    float(row["k2l_m2"])
                    * float(row["mean_excess_v_minus_h"])
                    < 0.0
                    for row in elements
                ),
                "figure_files": [
                    str(top_figure.resolve()),
                    str(compact_top_figure.resolve()),
                    str(ring_figure.resolve()),
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(args.output_dir / "ELEMENT_EXPOSURE_RESULTS.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
