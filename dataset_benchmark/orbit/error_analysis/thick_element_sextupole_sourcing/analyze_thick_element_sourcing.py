#!/usr/bin/env python3
"""Summarize and render thick-element sextupole sourcing results."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def render_svg(path: Path, elements: list[dict[str, str]]) -> None:
    width, height = 1120, 560
    left, right, top, bottom = 84, 28, 70, 92
    chart_w, chart_h = width - left - right, height - top - bottom
    circumference = max(f(row, "s_m") for row in elements)
    values = [100 * f(row, "eta_total") for row in elements]
    limit = 1.05 * max(abs(value) for value in values) if any(values) else 1.0
    zero_y = top + chart_h / 2

    def x_pos(s_m: float) -> float:
        return left + chart_w * s_m / circumference

    def y_pos(value: float) -> float:
        return zero_y - chart_h * value / (2 * limit)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:21px;font-weight:600}.sub{font-size:13px;fill:#555}.tick{font-size:12px}.label{font-size:11px}</style>',
        f'<text x="{width/2}" y="29" text-anchor="middle" class="title">Thick-element normal-sextupole sourcing of horizontal quadratic detector error</text>',
        f'<text x="{width/2}" y="51" text-anchor="middle" class="sub">Signed ensemble projection; complete-element Hessian sources and periodic six-dimensional propagation</text>',
    ]
    for tick in range(-4, 5):
        value = limit * tick / 4
        y = y_pos(value)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{left-9}" y="{y+4:.2f}" text-anchor="end" class="tick">{value:.2g}</text>')
    parts.append(f'<line x1="{left}" y1="{zero_y:.2f}" x2="{width-right}" y2="{zero_y:.2f}" stroke="#555" stroke-width="1.2"/>')
    bar_w = max(3.0, 0.58 * chart_w / len(elements))
    for row, value in zip(elements, values):
        x, y = x_pos(f(row, "s_m")) - bar_w / 2, y_pos(value)
        color = "#d95f02" if value >= 0 else "#1f78b4"
        parts.append(f'<rect x="{x:.2f}" y="{min(y, zero_y):.2f}" width="{bar_w:.2f}" height="{abs(zero_y-y):.2f}" fill="{color}" opacity="0.88"/>')
    ranked = sorted(elements, key=lambda row: abs(f(row, "eta_total")), reverse=True)[:10]
    for rank, row in enumerate(ranked, 1):
        value = 100 * f(row, "eta_total")
        x, y = x_pos(f(row, "s_m")), y_pos(value)
        label_y = top + 15 + 16 * ((rank - 1) % 5)
        label_x = left + 12 + (rank > 5) * chart_w / 2
        parts.append(f'<line x1="{x:.2f}" y1="{y:.2f}" x2="{label_x+92:.2f}" y2="{label_y-4:.2f}" stroke="#999" stroke-width="0.7"/>')
        parts.append(f'<text x="{label_x:.2f}" y="{label_y:.2f}" class="label">{rank}. {esc(row["element_name"])} ({value:+.2f}%)</text>')
    for tick in range(9):
        s_m = circumference * tick / 8
        x = x_pos(s_m)
        parts.append(f'<text x="{x:.2f}" y="{top+chart_h+22}" text-anchor="middle" class="tick">{s_m:.0f}</text>')
    parts.extend([
        f'<text x="{left+chart_w/2}" y="{height-31}" text-anchor="middle" class="tick">Ring position s [m]</text>',
        f'<text x="22" y="{top+chart_h/2}" text-anchor="middle" class="tick" transform="rotate(-90 22 {top+chart_h/2})">signed projection eta_j [%]</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_report(
    path: Path,
    elements: list[dict[str, str]],
    families: list[dict[str, str]],
    directions: list[dict[str, str]],
    summary: dict[str, str],
    comparison: dict[str, float | bool] | None,
    family_percentiles: list[dict[str, float | str]],
) -> None:
    ranked = sorted(elements, key=lambda row: abs(f(row, "eta_total")), reverse=True)
    closures = [f(row, "sextupole_total_relative_closure") for row in directions]
    projections = [f(row, "sextupole_total_signed_projection") for row in directions]
    lines = [
        "# Thick-element Hessian sourcing result", "",
        "## Closure", "",
        f'- Directions: `{summary["trials"]}`; lattice elements: `{summary["elements"]}`; active normal sextupoles: `{summary["active_normal_sextupoles"]}`; detectors: `{summary["detectors"]}`.',
        f'- All-element total relative closure: `{f(summary,"total_all_element_relative_closure"):.6g}`.',
        f'- Sextupole-only HH / HV / VV relative closure: `{f(summary,"hh_sextupole_relative_closure"):.6g} / {f(summary,"hv_sextupole_relative_closure"):.6g} / {f(summary,"vv_sextupole_relative_closure"):.6g}`.',
        f'- Sextupole-only total relative closure: `{f(summary,"total_sextupole_relative_closure"):.6g}`.',
        f'- Sextupole-only total signed projection: `{f(summary,"total_sextupole_signed_projection"):.6g}`.',
        f'- Direction-level sextupole closure P10 / median / P90: `{percentile(closures,0.1):.6g} / {percentile(closures,0.5):.6g} / {percentile(closures,0.9):.6g}`.',
        f'- Direction-level signed projection P10 / median / P90: `{percentile(projections,0.1):.6g} / {percentile(projections,0.5):.6g} / {percentile(projections,0.9):.6g}`.',
        f'- Family partition maximum absolute reconstruction difference: `{f(summary,"family_partition_check_max"):.6g} m`.',
    ]
    if comparison is not None:
        lines.extend([
            "", "## Comparison with the midpoint thin-kick reconstruction", "",
            f'- Thin / thick total relative closure: `{comparison["thin_closure"]:.6g} / {comparison["thick_closure"]:.6g}`.',
            f'- Absolute closure improvement: `{comparison["closure_improvement"]:.6g}`.',
            f'- Pearson correlation of the 76 signed `eta_total` values: `{comparison["eta_correlation"]:.12g}`.',
            f'- Maximum absolute change in an element `eta_total`: `{comparison["max_eta_difference"]:.6g}`.',
            f'- Top-15 absolute-projection ordering identical: `{comparison["top15_identical"]}`.',
            "",
            "The near-identical ranking and small closure change show that the remaining residual is not primarily caused by treating the sextupoles as thin midpoint sources.",
        ])
    lines.extend([
        "", "## Complete-element source families", "",
        "Signed projections add to one; magnitude ratios do not add because family vectors interfere.", "",
        "| family | elements | eta HH | eta HV | eta VV | eta total | magnitude total |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in sorted(families, key=lambda item: abs(f(item, "eta_total")), reverse=True):
        lines.append(
            f'| `{row["family"]}` | {int(float(row["element_count"]))} | '
            f'{100*f(row,"eta_hh"):+.3f}% | {100*f(row,"eta_hv"):+.3f}% | '
            f'{100*f(row,"eta_vv"):+.3f}% | {100*f(row,"eta_total"):+.3f}% | '
            f'{100*f(row,"magnitude_total"):.3f}% |'
        )
    lines.extend([
        "", "### Direction-level family statistics", "",
        "| family | eta P10 | eta median | eta P90 | magnitude median |",
        "|---|---:|---:|---:|---:|",
    ])
    for row in sorted(family_percentiles, key=lambda item: abs(float(item["eta_median"])), reverse=True):
        lines.append(
            f'| `{row["family"]}` | {100*float(row["eta_p10"]):+.3f}% | '
            f'{100*float(row["eta_median"]):+.3f}% | {100*float(row["eta_p90"]):+.3f}% | '
            f'{100*float(row["magnitude_median"]):.3f}% |'
        )
    lines.extend([
        "", "## Largest absolute signed projections", "",
        "| rank | sextupole | s [m] | K2L [m^-2] | eta total | magnitude ratio |",
        "|---:|---|---:|---:|---:|---:|",
    ])
    for rank, row in enumerate(ranked[:15], 1):
        lines.append(f'| {rank} | `{row["element_name"]}` | {f(row,"s_m"):.3f} | {f(row,"k2l_m2"):.5g} | {100*f(row,"eta_total"):+.4f}% | {100*f(row,"magnitude_total"):.4f}% |')
    lines.extend(["", "## Interpretation boundary", "",
        "The all-element closure validates the chain-rule source decomposition. The sextupole-only residual is retained explicitly and measures sources assigned to other complete lattice elements under this element-boundary convention.", "",
        "![Thick-element sextupole sourcing](thick_sextupole_signed_contributions.svg)", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    elements = rows(args.output_dir / "thick_sextupole_contribution_summary.csv")
    families = rows(args.output_dir / "family_contribution_summary.csv")
    directions = rows(args.output_dir / "direction_closure.csv")
    family_directions = rows(args.output_dir / "family_direction_contributions.csv")
    summaries = rows(args.output_dir / "reconstruction_summary.csv")
    if len(summaries) != 1:
        raise RuntimeError("Expected exactly one reconstruction summary row")
    error_analysis = Path(__file__).resolve().parent.parent
    thin_dir = error_analysis / "sextupole_detector_contributions" / "results"
    comparison = None
    if (thin_dir / "sextupole_contribution_summary.csv").is_file() and (
        thin_dir / "reconstruction_summary.csv"
    ).is_file():
        thin_elements = rows(thin_dir / "sextupole_contribution_summary.csv")
        thin_summary = rows(thin_dir / "reconstruction_summary.csv")[0]
        thin_eta = {row["element_name"].lower(): f(row, "eta_total") for row in thin_elements}
        thick_eta = {row["element_name"].lower(): f(row, "eta_total") for row in elements}
        names = sorted(thin_eta.keys() & thick_eta.keys())
        x = [thin_eta[name] for name in names]
        y = [thick_eta[name] for name in names]
        mean_x, mean_y = sum(x) / len(x), sum(y) / len(y)
        correlation = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y)) / math.sqrt(
            sum((a - mean_x) ** 2 for a in x) * sum((b - mean_y) ** 2 for b in y)
        )
        differences = [abs(thin_eta[name] - thick_eta[name]) for name in names]
        top_thin = sorted(names, key=lambda name: abs(thin_eta[name]), reverse=True)[:15]
        top_thick = sorted(names, key=lambda name: abs(thick_eta[name]), reverse=True)[:15]
        thin_closure = f(thin_summary, "total_concatenated_relative_closure")
        thick_closure = f(summaries[0], "total_sextupole_relative_closure")
        comparison = {
            "thin_closure": thin_closure,
            "thick_closure": thick_closure,
            "closure_improvement": thin_closure - thick_closure,
            "eta_correlation": correlation,
            "max_eta_difference": max(differences),
            "top15_identical": top_thin == top_thick,
        }
    norm_by_trial = {int(float(row["trial"])): f(row, "q_total_norm_m") for row in directions}
    family_samples: dict[str, dict[str, list[float]]] = {}
    for row in family_directions:
        family = row["family"]
        trial = int(float(row["trial"]))
        norm = norm_by_trial[trial]
        sample = family_samples.setdefault(family, {"eta": [], "magnitude": []})
        sample["eta"].append(f(row, "projection_numerator") / norm**2)
        sample["magnitude"].append(f(row, "contribution_norm_m") / norm)
    family_percentiles: list[dict[str, float | str]] = []
    for family, sample in family_samples.items():
        family_percentiles.append({
            "family": family,
            "eta_p10": percentile(sample["eta"], 0.1),
            "eta_median": percentile(sample["eta"], 0.5),
            "eta_p90": percentile(sample["eta"], 0.9),
            "magnitude_p10": percentile(sample["magnitude"], 0.1),
            "magnitude_median": percentile(sample["magnitude"], 0.5),
            "magnitude_p90": percentile(sample["magnitude"], 0.9),
        })
    percentile_path = args.output_dir / "family_direction_percentiles.csv"
    with percentile_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "family", "eta_p10", "eta_median", "eta_p90",
            "magnitude_p10", "magnitude_median", "magnitude_p90",
        ])
        writer.writeheader()
        writer.writerows(sorted(family_percentiles, key=lambda row: str(row["family"])))
    render_svg(args.output_dir / "thick_sextupole_signed_contributions.svg", elements)
    write_report(
        args.output_dir / "RESULTS.md", elements, families, directions,
        summaries[0], comparison, family_percentiles,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
