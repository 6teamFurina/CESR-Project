#!/usr/bin/env python3
"""Render the orbit-response rho sweep as a dependency-free SVG chart."""

from __future__ import annotations

import argparse
import csv
import math
from html import escape
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_FIGURES = HERE / "error_analysis" / "response_rho_sweep_600" / "figures"
RESULT_CANDIDATES = (
    HERE / "error_analysis" / "response_rho_sweep_600" / "combined",
    HERE / "results" / "response_rho_sweep_600" / "combined",
    HERE / "results" / "response_rho_sweep",
)


def default_results() -> Path:
    return next(
        (path for path in RESULT_CANDIDATES if (path / "rho_sweep_summary.csv").is_file()),
        RESULT_CANDIDATES[0],
    )


def arguments() -> argparse.Namespace:
    results = default_results()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        type=Path,
        default=results / "rho_sweep_summary.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--normalize-rho-squared",
        action="store_true",
        help="Plot RMSE divided by rho squared to expose departure from quadratic scaling.",
    )
    parser.add_argument(
        "--base-kick-urad",
        type=float,
        default=5.0,
        help="Active-corrector RMS kick represented by rho=1 (default: 5 microrad).",
    )
    args = parser.parse_args()
    if args.output is None:
        filename = (
            "scibmad_orbit_response_error_rho2_normalized.svg"
            if args.normalize_rho_squared
            else "scibmad_orbit_response_error.svg"
        )
        args.output = DEFAULT_FIGURES / filename
    return args


def read_summary(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def text(x: float, y: float, value: str, **attributes: object) -> str:
    properties = " ".join(
        f'{name.replace("_", "-")}="{item}"' for name, item in attributes.items()
    )
    return f'<text x="{x:.2f}" y="{y:.2f}" {properties}>{escape(value)}</text>'


def power_label(x: float, y: float, power: int) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="end" font-size="12" fill="#4A4A4A">'
        f'10<tspan baseline-shift="super" font-size="9">{power}</tspan></text>'
    )


def marker(shape: str, x: float, y: float, color: str) -> str:
    common = f'fill="white" stroke="{color}" stroke-width="1.7"'
    if shape == "circle":
        return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.5" {common}/>'
    if shape == "square":
        return f'<rect x="{x - 3.5:.2f}" y="{y - 3.5:.2f}" width="7" height="7" {common}/>'
    if shape == "triangle":
        points = f"{x:.2f},{y - 4.2:.2f} {x - 4:.2f},{y + 3.5:.2f} {x + 4:.2f},{y + 3.5:.2f}"
        return f'<polygon points="{points}" {common}/>'
    raise ValueError(f"Unknown marker: {shape}")


def cross_marker(
    x: float,
    y: float,
    color: str,
    *,
    size: float = 5.0,
    stroke_width: float = 2.0,
) -> str:
    return "\n".join(
        (
            f'<line x1="{x - size:.2f}" y1="{y - size:.2f}" x2="{x + size:.2f}" y2="{y + size:.2f}" stroke="{color}" stroke-width="{stroke_width:g}"/>',
            f'<line x1="{x - size:.2f}" y1="{y + size:.2f}" x2="{x + size:.2f}" y2="{y - size:.2f}" stroke="{color}" stroke-width="{stroke_width:g}"/>',
        )
    )


def log_ticks(lower: float, upper: float) -> tuple[list[float], list[float]]:
    major = [
        10.0**power
        for power in range(math.floor(math.log10(lower)), math.ceil(math.log10(upper)) + 1)
        if lower <= 10.0**power <= upper
    ]
    for endpoint in (lower, upper):
        if not any(math.isclose(endpoint, tick, rel_tol=1.0e-10) for tick in major):
            major.append(endpoint)
    minor = [
        multiplier * 10.0**power
        for power in range(math.floor(math.log10(lower)), math.ceil(math.log10(upper)) + 1)
        for multiplier in (2.0, 5.0)
        if lower < multiplier * 10.0**power < upper
    ]
    return sorted(major), minor


def format_log_tick(value: float) -> str:
    if value >= 10.0:
        return f"{value:g}"
    if value >= 1.0:
        return f"{value:.0f}"
    return f"{value:.2g}"


def main() -> int:
    args = arguments()
    if not math.isfinite(args.base_kick_urad) or args.base_kick_urad <= 0:
        raise ValueError("--base-kick-urad must be finite and positive")
    summary = args.summary.expanduser().resolve()
    output = args.output.expanduser().resolve()
    rows = read_summary(summary)
    trial_counts = sorted({int(row["trials"]) for row in rows})
    trial_note = (
        f"{trial_counts[0]} random directions per input scale"
        if len(trial_counts) == 1
        else "trial count varies by input scale"
    )
    chart_title = (
        "Quadratic-normalized orbit-response residual"
        if args.normalize_rho_squared
        else "First-order orbit-response error versus corrector input scale"
    )
    y_label = (
        "Detector-orbit RMSE / rho squared [micrometre]"
        if args.normalize_rho_squared
        else "Detector-orbit response residual, RMSE [micrometre]"
    )

    def plotted_value(row: dict[str, str], field: str) -> float:
        value = 1.0e6 * float(row[field])
        if args.normalize_rho_squared:
            rho = float(row["rho"])
            return value / (rho * rho)
        return value
    scenarios = (
        ("all", "All correctors (119)", "#0072B2", "", "circle"),
        ("horizontal", "Horizontal only (58)", "#D55E00", "9 4", "square"),
        ("vertical", "Vertical only (61)", "#009E73", "2 3", "triangle"),
    )
    width, height = 1260, 620
    top, bottom = 128, 132
    panel_width, panel_gap = 530, 92
    lefts = (94, 94 + panel_width + panel_gap)
    plot_width = panel_width - 34
    plot_height = height - top - bottom
    rho_values = sorted({float(row["rho"]) for row in rows})
    positive_rhos = [rho for rho in rho_values if rho > 0]
    if not positive_rhos:
        raise RuntimeError("At least one positive rho is required for the plot")
    rho_min, rho_max = min(positive_rhos), max(positive_rhos)
    major_x_ticks, minor_x_ticks = log_ticks(rho_min, rho_max)

    positive_errors_um: list[float] = []
    for row in rows:
        if float(row["rho"]) <= 0:
            continue
        for plane in ("x", "y"):
            for prefix in ("mean", "max_trial"):
                value = plotted_value(row, f"{prefix}_{plane}_rmse_m")
                if math.isfinite(value) and value > 0:
                    positive_errors_um.append(value)
    if not positive_errors_um:
        raise RuntimeError("No positive finite orbit errors are available")
    lower_power = math.floor(math.log10(min(positive_errors_um)))
    upper_power = math.ceil(math.log10(max(positive_errors_um)))
    if lower_power == upper_power:
        upper_power += 1

    def x_position(rho: float, left: float) -> float:
        if rho_min == rho_max:
            return left + plot_width / 2
        fraction = (math.log10(rho) - math.log10(rho_min)) / (
            math.log10(rho_max) - math.log10(rho_min)
        )
        return left + plot_width * fraction

    def y_position(value_um: float) -> float:
        fraction = (math.log10(value_um) - lower_power) / (upper_power - lower_power)
        return top + plot_height * (1.0 - fraction)

    subtitle = (
        f"CESR RF-on closed orbit; {trial_note}; rho = 1 corresponds to "
        f"{args.base_kick_urad:g} microrad active-corrector RMS"
    )
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title description">',
        f'<title id="title">{escape(chart_title)}</title>',
        '<desc id="description">Horizontal and vertical detector-orbit response residuals versus normalized active-corrector RMS input radius. Curves show mean trial RMSE, upper whiskers show maximum trial RMSE, and crosses mark incomplete exact-reference samples.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#202020" stroke-linecap="round" stroke-linejoin="round">',
        text(width / 2, 30, chart_title, text_anchor="middle", font_size="20", font_weight="600"),
        text(width / 2, 54, subtitle, text_anchor="middle", font_size="12.5", fill="#555555"),
    ]

    legend_y = 88
    for start, (_, label, color, dash, shape) in zip((70, 300, 535), scenarios):
        dash_attribute = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(f'<line x1="{start}" y1="{legend_y}" x2="{start + 34}" y2="{legend_y}" stroke="{color}" stroke-width="2.2"{dash_attribute}/>')
        parts.append(marker(shape, start + 17, legend_y, color))
        parts.append(text(start + 43, legend_y + 4, label, font_size="12.5"))
    guide_x = 770
    parts.append(f'<line x1="{guide_x}" y1="{legend_y}" x2="{guide_x + 34}" y2="{legend_y}" stroke="#6B6B6B" stroke-width="1.5" stroke-dasharray="7 5" opacity="0.8"/>')
    parts.append(text(guide_x + 43, legend_y + 4, "Pure quadratic-error reference", font_size="12.5"))
    incomplete_x = 1030
    parts.append(cross_marker(incomplete_x + 5, legend_y, "#555555"))
    parts.append(text(incomplete_x + 18, legend_y + 4, "Incomplete exact reference", font_size="12.5"))

    for panel, (plane, title) in enumerate((("x", "Horizontal detector orbit"), ("y", "Vertical detector orbit"))):
        left = lefts[panel]
        right = left + plot_width
        bottom_y = top + plot_height
        for power in range(lower_power, upper_power + 1):
            y = y_position(10.0**power)
            parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#D4D7DA" stroke-width="1"/>')
            if panel == 0:
                parts.append(power_label(left - 12, y + 4, power))
        for power in range(lower_power, upper_power):
            for multiplier in (2.0, 5.0):
                y = y_position(multiplier * 10.0**power)
                parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#ECEEEF" stroke-width="0.8"/>')
        for rho in minor_x_ticks:
            x = x_position(rho, left)
            parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom_y}" stroke="#F0F1F2" stroke-width="0.8"/>')
        for rho in major_x_ticks:
            x = x_position(rho, left)
            parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom_y}" stroke="#E0E2E4" stroke-width="1"/>')
            parts.append(text(x, bottom_y + 23, format_log_tick(rho), text_anchor="middle", font_size="12", fill="#4A4A4A"))
        parts.extend((
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#333333" stroke-width="1.15"/>',
            text(left, top - 17, f"({chr(97 + panel)})  {title}", font_size="15", font_weight="600"),
            text((left + right) / 2, bottom_y + 52, "Normalized active-corrector input radius, rho", text_anchor="middle", font_size="13.5"),
        ))

        reference_rows = sorted(
            (row for row in rows if row["scenario"] == "all" and float(row["rho"]) > 0 and int(row["converged_trials"]) == int(row["trials"])),
            key=lambda row: float(row["rho"]),
        )
        if reference_rows:
            anchor = reference_rows[0]
            anchor_rho = float(anchor["rho"])
            anchor_value = plotted_value(anchor, f"mean_{plane}_rmse_m")
            guide = []
            for rho in positive_rhos:
                value = anchor_value if args.normalize_rho_squared else anchor_value * (rho / anchor_rho) ** 2
                if 10.0**lower_power <= value <= 10.0**upper_power:
                    guide.append((rho, value))
            if len(guide) >= 2:
                guide_path = " ".join(("M" if index == 0 else "L") + f" {x_position(rho, left):.2f} {y_position(value):.2f}" for index, (rho, value) in enumerate(guide))
                parts.append(f'<path d="{guide_path}" fill="none" stroke="#6B6B6B" stroke-width="1.5" stroke-dasharray="7 5" opacity="0.8"/>')

        series_data: list[
            tuple[str, str, str, str, list[tuple[float, float, float, bool]]]
        ] = []
        for scenario, _, color, dash, shape in scenarios:
            selected = sorted(
                (
                    row
                    for row in rows
                    if row["scenario"] == scenario and float(row["rho"]) > 0
                ),
                key=lambda row: float(row["rho"]),
            )
            points: list[tuple[float, float, float, bool]] = []
            for row in selected:
                mean_um = plotted_value(row, f"mean_{plane}_rmse_m")
                maximum_um = plotted_value(row, f"max_trial_{plane}_rmse_m")
                if mean_um > 0 and math.isfinite(mean_um) and math.isfinite(maximum_um):
                    complete = int(row["converged_trials"]) == int(row["trials"])
                    points.append((float(row["rho"]), mean_um, maximum_um, complete))
            series_data.append((scenario, color, dash, shape, points))

        # Draw all mean curves first so every whisker remains visible above them.
        for _, color, dash, _, points in series_data:
            if len(points) >= 2:
                path = " ".join(
                    ("M" if index == 0 else "L")
                    + f" {x_position(rho, left):.2f} {y_position(mean_um):.2f}"
                    for index, (rho, mean_um, _, _) in enumerate(points)
                )
                dash_attribute = f' stroke-dasharray="{dash}"' if dash else ""
                parts.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2.15"{dash_attribute}/>')

        # The horizontal-only series is rendered last and slightly wider because
        # its shorter orange whiskers otherwise disappear beneath the other curves.
        whisker_order = {"all": 0, "vertical": 1, "horizontal": 2}
        for scenario, color, _, _, points in sorted(
            series_data, key=lambda item: whisker_order[item[0]]
        ):
            whisker_width = 2.15 if scenario == "horizontal" else 1.35
            whisker_opacity = 0.88 if scenario == "horizontal" else 0.62
            cap_half_width = 5.5 if scenario == "horizontal" else 4.0
            for rho, mean_um, maximum_um, complete in points:
                x = x_position(rho, left)
                y_mean = y_position(mean_um)
                y_maximum = y_position(maximum_um)
                parts.append(f'<line x1="{x:.2f}" y1="{y_mean:.2f}" x2="{x:.2f}" y2="{y_maximum:.2f}" stroke="{color}" stroke-width="{whisker_width:g}" opacity="{whisker_opacity:g}"/>')
                parts.append(f'<line x1="{x - cap_half_width:.2f}" y1="{y_maximum:.2f}" x2="{x + cap_half_width:.2f}" y2="{y_maximum:.2f}" stroke="{color}" stroke-width="{whisker_width:g}" opacity="{whisker_opacity:g}"/>')

        # Complete samples use the scenario marker. Incomplete samples use only
        # a same-color cross, avoiding a marker/cross overlay at the same point.
        for _, color, _, shape, points in series_data:
            for rho, mean_um, _, complete in points:
                x = x_position(rho, left)
                y_mean = y_position(mean_um)
                if not complete:
                    parts.append(cross_marker(x, y_mean, color))
                else:
                    parts.append(marker(shape, x, y_mean, color))

    center_y = top + plot_height / 2
    parts.append(f'<text x="24" y="{center_y:.2f}" transform="rotate(-90 24 {center_y:.2f})" text-anchor="middle" font-size="13.5">{escape(y_label)}</text>')
    parts.append(text(width / 2, height - 29, "Markers and lines: mean over trials; one-sided whiskers: maximum over trials.", text_anchor="middle", font_size="11.5", fill="#555555"))
    parts.extend(("</g>", "</svg>"))

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
