#!/usr/bin/env python3
"""Analyze and plot this directory's paired vertical-corrector experiment."""

from __future__ import annotations

import argparse
import csv
import math
from html import escape
from pathlib import Path


ORIGINAL_RHOS = (
    0.05, 0.075, 0.1, 0.14, 0.2, 0.28, 0.4, 0.57, 0.8, 1.13,
    1.6, 2.26, 3.2, 4.53, 6.4,
)
PERCENTILE_RHOS = ORIGINAL_RHOS + (7.5, 8.8, 10.05)

def local_slopes(rhos: list[float], values: list[float]) -> list[float]:
    slopes = [math.nan]
    for index in range(1, len(rhos)):
        slopes.append(
            math.log(values[index] / values[index - 1])
            / math.log(rhos[index] / rhos[index - 1])
        )
    return slopes


def load_summary(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    result: list[dict[str, float]] = []
    for row in rows:
        parsed = {key: float(value) for key, value in row.items()}
        result.append(parsed)
    return result


def load_pairs(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return [
        {
            "rho": float(row["rho"]),
            "trial": float(row["trial"]),
            "y_even_rmse_m": float(row["y_even_rmse_m"]),
            "y_odd_nl_rmse_m": float(row["y_odd_nl_rmse_m"]),
        }
        for row in rows
        if row["plus_converged"].lower() == "true"
        and row["minus_converged"].lower() == "true"
    ]


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def direction_percentiles(pairs: list[dict[str, float]]) -> list[dict[str, float]]:
    result = []
    for rho in PERCENTILE_RHOS:
        selected = [row for row in pairs if math.isclose(row["rho"], rho)]
        if not selected:
            raise ValueError(f"No pair-level rows found for original rho={rho:g}")
        ratios = [row["y_odd_nl_rmse_m"] / row["y_even_rmse_m"] for row in selected]
        fractions = [ratio**2 / (1.0 + ratio**2) for ratio in ratios]
        row = {"rho": rho, "directions": float(len(selected))}
        for label, probability in (("p10", 0.1), ("median", 0.5), ("p90", 0.9)):
            row[f"ratio_{label}"] = percentile(ratios, probability)
            row[f"fraction_{label}"] = percentile(fractions, probability)
        result.append(row)
    return result


def write_direction_percentiles(path: Path, rows: list[dict[str, float]]) -> None:
    fieldnames = [
        "rho", "directions", "ratio_p10", "ratio_median", "ratio_p90",
        "fraction_p10", "fraction_median", "fraction_p90",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_analysis(path: Path, rows: list[dict[str, float]]) -> None:
    rhos = [row["rho"] for row in rows]
    x_even = [row["mean_x_even_rmse_m"] for row in rows]
    x_odd = [row["mean_x_odd_nl_rmse_m"] for row in rows]
    y_even = [row["mean_y_even_rmse_m"] for row in rows]
    y_odd = [row["mean_y_odd_nl_rmse_m"] for row in rows]
    x_even_slopes = local_slopes(rhos, x_even)
    x_odd_slopes = local_slopes(rhos, x_odd)
    y_even_slopes = local_slopes(rhos, y_even)
    y_odd_slopes = local_slopes(rhos, y_odd)

    lines = [
        "# Vertical-corrector signed-parity experiment",
        "",
        "The same random vertical-corrector direction was evaluated at both signs.",
        "`even` contains even Taylor orders; `odd_nl` is the odd part after subtracting",
        "the nominal first-order detector response.",
        "",
        "| rho | kick RMS (urad) | X even (um) | X odd-nl (um) | X p_even | X p_odd | Y even (um) | Y odd-nl (um) | Y p_even | Y p_odd | Y odd/even |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for index, rho in enumerate(rhos):
        def slope_text(value: float) -> str:
            return "n/a" if math.isnan(value) else f"{value:.3f}"

        lines.append(
            f"| {rho:g} | {5 * rho:.3g} | {1e6*x_even[index]:.6g} | "
            f"{1e6*x_odd[index]:.6g} | {slope_text(x_even_slopes[index])} | "
            f"{slope_text(x_odd_slopes[index])} | {1e6*y_even[index]:.6g} | "
            f"{1e6*y_odd[index]:.6g} | {slope_text(y_even_slopes[index])} | "
            f"{slope_text(y_odd_slopes[index])} | {y_odd[index]/y_even[index]:.4g} |"
        )

    lines.extend(
        [
            "",
            "## Compact checks",
            "",
            f"- Smallest-rho X even/rho^2: `{1e6*x_even[0]/rhos[0]**2:.8g} um`.",
            f"- Smallest-rho Y even/rho^2: `{1e6*y_even[0]/rhos[0]**2:.8g} um`.",
            f"- Largest-rho X odd/even ratio: `{x_odd[-1]/x_even[-1]:.8g}`.",
            f"- Largest-rho Y odd/even ratio: `{y_odd[-1]/y_even[-1]:.8g}`.",
            f"- Final Y even local slope: `{y_even_slopes[-1]:.8g}`.",
            f"- Final Y odd-nonlinear local slope: `{y_odd_slopes[-1]:.8g}`.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def render(path: Path, rows: list[dict[str, float]]) -> None:
    rhos = [row["rho"] for row in rows]
    plane_values = {}
    all_values: list[float] = []
    for plane in ("x", "y"):
        even = [1e6 * row[f"mean_{plane}_even_rmse_m"] for row in rows]
        odd = [1e6 * row[f"mean_{plane}_odd_nl_rmse_m"] for row in rows]
        plane_values[plane] = (even, odd)
        all_values.extend(even)
        all_values.extend(odd)

    width, height = 1160, 570
    top, plot_height, plot_width, gap = 125, 350, 465, 90
    lefts = (85, 85 + plot_width + gap)
    rho_lo, rho_hi = min(rhos), max(rhos)
    y_lo = 10.0 ** math.floor(math.log10(min(all_values)))
    y_hi = 10.0 ** math.ceil(math.log10(max(all_values)))

    def xp(rho: float, left: float) -> float:
        return left + plot_width * (
            (math.log10(rho) - math.log10(rho_lo))
            / (math.log10(rho_hi) - math.log10(rho_lo))
        )

    def yp(value: float) -> float:
        return top + plot_height * (
            1.0
            - (math.log10(value) - math.log10(y_lo))
            / (math.log10(y_hi) - math.log10(y_lo))
        )

    def polyline(values: list[float], left: float, color: str, dash: str = "") -> str:
        points = " ".join(f"{xp(rho,left):.2f},{yp(value):.2f}" for rho, value in zip(rhos, values))
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        return f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.2"{dash_attr}/>'

    def marker(shape: str, x: float, y: float, color: str) -> str:
        common = f'fill="white" stroke="{color}" stroke-width="1.5"'
        if shape == "square":
            return f'<rect x="{x-3:.2f}" y="{y-3:.2f}" width="6" height="6" {common}/>'
        points = f"{x:.2f},{y-3.8:.2f} {x-3.6:.2f},{y+3:.2f} {x+3.6:.2f},{y+3:.2f}"
        return f'<polygon points="{points}" {common}/>'

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#202020">',
        '<text x="580" y="32" text-anchor="middle" font-size="21" font-weight="600">Vertical-corrector signed-parity decomposition</text>',
        '<text x="580" y="56" text-anchor="middle" font-size="13" fill="#555">100 paired random directions; rho = 1 is 5 microrad active-corrector RMS</text>',
        '<line x1="95" y1="88" x2="129" y2="88" stroke="#D55E00" stroke-width="2.4"/>',
        '<rect x="109" y="85" width="6" height="6" fill="white" stroke="#D55E00" stroke-width="1.5"/>',
        '<text x="138" y="93" font-size="13">Even component</text>',
        '<line x1="315" y1="88" x2="349" y2="88" stroke="#009E73" stroke-width="2.4"/>',
        '<polygon points="332,84.2 328.4,91 335.6,91" fill="white" stroke="#009E73" stroke-width="1.5"/>',
        '<text x="358" y="93" font-size="13">Odd nonlinear component</text>',
        '<line x1="610" y1="88" x2="644" y2="88" stroke="#D55E00" stroke-width="1.8" stroke-dasharray="7 5"/>',
        '<text x="653" y="93" font-size="13">Quadratic guide (rho^2)</text>',
        '<line x1="865" y1="88" x2="899" y2="88" stroke="#009E73" stroke-width="1.8" stroke-dasharray="7 5"/>',
        '<text x="908" y="93" font-size="13">Cubic guide (rho^3)</text>',
    ]

    for panel, plane in enumerate(("x", "y")):
        left = lefts[panel]
        right = left + plot_width
        bottom = top + plot_height
        for power in range(int(math.log10(y_lo)), int(math.log10(y_hi)) + 1):
            value = 10.0**power
            y = yp(value)
            parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#E0E3E5"/>')
            if panel == 0:
                parts.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" font-size="11">10^{power}</text>')
        for tick in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0):
            if rho_lo <= tick <= rho_hi:
                x = xp(tick, left)
                parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#ECEEEF"/>')
                parts.append(f'<text x="{x:.2f}" y="{bottom+20}" text-anchor="middle" font-size="11">{tick:g}</text>')
        parts.append(f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#555"/>')
        parts.append(f'<text x="{(left+right)/2:.2f}" y="{top-12}" text-anchor="middle" font-size="15" font-weight="600">{plane.upper()} detector orbit</text>')
        even, odd = plane_values[plane]
        even_guide = [even[0] * (rho / rhos[0]) ** 2 for rho in rhos]
        odd_guide = [odd[0] * (rho / rhos[0]) ** 3 for rho in rhos]
        parts.append(polyline(even_guide, left, "#D55E00", "7 5"))
        parts.append(polyline(odd_guide, left, "#009E73", "7 5"))
        parts.append(polyline(even, left, "#D55E00"))
        parts.append(polyline(odd, left, "#009E73"))
        for values, color, shape in (
            (even, "#D55E00", "square"),
            (odd, "#009E73", "triangle"),
        ):
            for rho, value in zip(rhos, values):
                parts.append(marker(shape, xp(rho, left), yp(value), color))
        parts.append(f'<text x="{(left+right)/2:.2f}" y="{bottom+48}" text-anchor="middle" font-size="13">Normalized vertical-corrector radius, rho</text>')

    parts.extend(
        [
            '<text x="20" y="300" text-anchor="middle" font-size="13" transform="rotate(-90 20 300)">Mean detector RMSE [micrometre]</text>',
            '</g>',
            '</svg>',
        ]
    )
    path.with_suffix(".svg").write_text("\n".join(parts), encoding="utf-8")
    render_linear_growth(path.with_name("vertical_parity_growth_linear"), rows)


def render_linear_growth(path: Path, rows: list[dict[str, float]]) -> None:
    rhos = [row["rho"] for row in rows]
    even = [1e6 * row["mean_y_even_rmse_m"] for row in rows]
    odd = [1e6 * row["mean_y_odd_nl_rmse_m"] for row in rows]
    rho_hi = max(rhos)
    y_hi = math.ceil(4 * max(even + odd) * 1.06) / 4
    even_coefficients = sorted(value / rho**2 for rho, value in zip(rhos, even))
    odd_coefficients = sorted(value / rho**3 for rho, value in zip(rhos, odd))
    middle = len(rhos) // 2
    crossover = even_coefficients[middle] / odd_coefficients[middle]

    width, height = 1040, 580
    left, right, top, bottom = 92, 952, 108, 480
    plot_width, plot_height = right - left, bottom - top

    def xp(rho: float) -> float:
        return left + plot_width * rho / rho_hi

    def yp(value: float) -> float:
        return bottom - plot_height * value / y_hi

    def polyline(values: list[float], color: str) -> str:
        points = " ".join(f"{xp(rho):.2f},{yp(value):.2f}" for rho, value in zip(rhos, values))
        return f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.6"/>'

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="linear-title linear-desc">',
        '<title id="linear-title">Vertical detector even and odd nonlinear growth on a linear scale</title>',
        '<desc id="linear-desc">The odd cubic component overtakes the even quadratic component near rho 1.31 and grows to 4.91 times the even component at rho 6.4.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#202020">',
        '<text x="520" y="32" text-anchor="middle" font-size="21" font-weight="600">Vertical detector: cubic growth versus quadratic growth</text>',
        '<text x="520" y="56" text-anchor="middle" font-size="13" fill="#555">Linear y-axis; 100 paired vertical-corrector directions</text>',
        '<line x1="322" y1="82" x2="358" y2="82" stroke="#D55E00" stroke-width="2.6"/>',
        '<rect x="337" y="79" width="6" height="6" fill="white" stroke="#D55E00" stroke-width="1.5"/>',
        '<text x="368" y="87" font-size="13">Even quadratic component</text>',
        '<line x1="590" y1="82" x2="626" y2="82" stroke="#009E73" stroke-width="2.6"/>',
        '<polygon points="608,78.2 604.4,85 611.6,85" fill="white" stroke="#009E73" stroke-width="1.5"/>',
        '<text x="636" y="87" font-size="13">Odd cubic component</text>',
    ]

    y_step = 0.25
    for index in range(int(round(y_hi / y_step)) + 1):
        value = index * y_step
        y = yp(value)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#E1E4E6"/>')
        parts.append(f'<text x="{left-12}" y="{y+4:.2f}" text-anchor="end" font-size="11">{value:g}</text>')
    for value in (0, 1, 2, 3, 4, 5, 6, rho_hi):
        x = xp(value)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#ECEEEF"/>')
        parts.append(f'<text x="{x:.2f}" y="{bottom+20}" text-anchor="middle" font-size="11">{value:g}</text>')

    cross_x = xp(crossover)
    parts.append(f'<line x1="{cross_x:.2f}" y1="{top}" x2="{cross_x:.2f}" y2="{bottom}" stroke="#666" stroke-width="1.4" stroke-dasharray="5 5"/>')
    parts.append(f'<text x="{cross_x+8:.2f}" y="{top+18}" font-size="12" fill="#555">Mean crossover: rho = {crossover:.2f}</text>')
    parts.append(polyline(even, "#D55E00"))
    parts.append(polyline(odd, "#009E73"))

    for rho, value in zip(rhos, even):
        x, y = xp(rho), yp(value)
        parts.append(f'<rect x="{x-3:.2f}" y="{y-3:.2f}" width="6" height="6" fill="white" stroke="#D55E00" stroke-width="1.5"/>')
    for rho, value in zip(rhos, odd):
        x, y = xp(rho), yp(value)
        points = f"{x:.2f},{y-3.8:.2f} {x-3.6:.2f},{y+3:.2f} {x+3.6:.2f},{y+3:.2f}"
        parts.append(f'<polygon points="{points}" fill="white" stroke="#009E73" stroke-width="1.5"/>')

    ratio = odd[-1] / even[-1]
    parts.append(f'<text x="{right-10}" y="{yp(odd[-1])-12:.2f}" text-anchor="end" font-size="12">Odd: {odd[-1]:.3f} um ({ratio:.2f}x even)</text>')
    parts.append(f'<text x="{right-10}" y="{yp(even[-1])-12:.2f}" text-anchor="end" font-size="12">Even: {even[-1]:.3f} um</text>')
    parts.extend(
        [
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#555"/>',
            f'<text x="{(left+right)/2:.2f}" y="{bottom+50}" text-anchor="middle" font-size="13">Normalized vertical-corrector radius, rho</text>',
            '<text x="22" y="294" text-anchor="middle" font-size="13" transform="rotate(-90 22 294)">Mean detector RMSE [micrometre]</text>',
            '</g>',
            '</svg>',
        ]
    )
    path.with_suffix(".svg").write_text("\n".join(parts), encoding="utf-8")


def render_crossover_loglog(path: Path, rows: list[dict[str, float]]) -> None:
    focused = [row for row in rows if 0.4 <= row["rho"] <= 3.6]
    rhos = [row["rho"] for row in focused]
    even = [1e6 * row["mean_y_even_rmse_m"] for row in focused]
    odd = [1e6 * row["mean_y_odd_nl_rmse_m"] for row in focused]
    ratios = [odd_value / even_value for even_value, odd_value in zip(even, odd)]
    crossover = math.nan
    for index in range(1, len(rhos)):
        if ratios[index - 1] <= 1.0 <= ratios[index]:
            fraction = -math.log(ratios[index - 1]) / math.log(
                ratios[index] / ratios[index - 1]
            )
            crossover = math.exp(
                math.log(rhos[index - 1])
                + fraction * math.log(rhos[index] / rhos[index - 1])
            )
            break
    if math.isnan(crossover):
        raise ValueError("Focused rho range does not bracket the even/odd crossover")

    width, height = 960, 600
    left, right, top, bottom = 105, 875, 105, 500
    plot_width, plot_height = right - left, bottom - top
    rho_lo, rho_hi = min(rhos), max(rhos)
    all_values = even + odd
    y_lo = 10.0 ** math.floor(math.log10(min(all_values)))
    y_hi = 10.0 ** math.ceil(math.log10(max(all_values)))

    def xp(rho: float) -> float:
        return left + plot_width * (
            (math.log10(rho) - math.log10(rho_lo))
            / (math.log10(rho_hi) - math.log10(rho_lo))
        )

    def yp(value: float) -> float:
        return bottom - plot_height * (
            (math.log10(value) - math.log10(y_lo))
            / (math.log10(y_hi) - math.log10(y_lo))
        )

    def points(values: list[float]) -> str:
        return " ".join(
            f"{xp(rho):.2f},{yp(value):.2f}"
            for rho, value in zip(rhos, values)
        )

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="cross-title cross-desc">',
        '<title id="cross-title">Vertical even and odd response near their crossover</title>',
        '<desc id="cross-desc">A focused log-log plot showing the quadratic even and cubic odd nonlinear mean responses crossing near rho 1.31.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#202020">',
        '<text x="480" y="32" text-anchor="middle" font-size="21" font-weight="600">Vertical response near the even–odd crossover</text>',
        '<text x="480" y="56" text-anchor="middle" font-size="13" fill="#555">Focused log–log view; mean across the same 100 corrector directions</text>',
        '<line x1="246" y1="82" x2="282" y2="82" stroke="#D55E00" stroke-width="2.6"/>',
        '<rect x="261" y="79" width="6" height="6" fill="white" stroke="#D55E00" stroke-width="1.5"/>',
        '<text x="292" y="87" font-size="13">Even component (rho²)</text>',
        '<line x1="510" y1="82" x2="546" y2="82" stroke="#009E73" stroke-width="2.6"/>',
        '<polygon points="528,78.2 524.4,85 531.6,85" fill="white" stroke="#009E73" stroke-width="1.5"/>',
        '<text x="556" y="87" font-size="13">Odd nonlinear component (rho³)</text>',
    ]
    for power in range(int(math.log10(y_lo)), int(math.log10(y_hi)) + 1):
        value = 10.0**power
        y = yp(value)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#E1E4E6"/>')
        parts.append(f'<text x="{left-12}" y="{y+4:.2f}" text-anchor="end" font-size="11">10^{power}</text>')
    for tick in (0.4, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 3.6):
        x = xp(tick)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#ECEEEF"/>')
        parts.append(f'<text x="{x:.2f}" y="{bottom+20}" text-anchor="middle" font-size="11">{tick:g}</text>')

    cross_x = xp(crossover)
    parts.extend(
        [
            f'<line x1="{cross_x:.2f}" y1="{top}" x2="{cross_x:.2f}" y2="{bottom}" stroke="#666" stroke-width="1.5" stroke-dasharray="6 5"/>',
            f'<text x="{cross_x+9:.2f}" y="{top+20}" font-size="12" fill="#555">Crossover: rho = {crossover:.2f}</text>',
            f'<polyline points="{points(even)}" fill="none" stroke="#D55E00" stroke-width="2.6"/>',
            f'<polyline points="{points(odd)}" fill="none" stroke="#009E73" stroke-width="2.6"/>',
        ]
    )
    for rho, value in zip(rhos, even):
        x, y = xp(rho), yp(value)
        parts.append(f'<rect x="{x-3:.2f}" y="{y-3:.2f}" width="6" height="6" fill="white" stroke="#D55E00" stroke-width="1.5"/>')
    for rho, value in zip(rhos, odd):
        x, y = xp(rho), yp(value)
        triangle = f"{x:.2f},{y-3.8:.2f} {x-3.6:.2f},{y+3:.2f} {x+3.6:.2f},{y+3:.2f}"
        parts.append(f'<polygon points="{triangle}" fill="white" stroke="#009E73" stroke-width="1.5"/>')
    parts.extend(
        [
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#555"/>',
            f'<text x="{(left+right)/2:.2f}" y="{bottom+52}" text-anchor="middle" font-size="13">Normalized vertical-corrector radius, rho (log scale)</text>',
            '<text x="25" y="302" text-anchor="middle" font-size="13" transform="rotate(-90 25 302)">Mean detector RMSE [micrometre] (log scale)</text>',
            '</g>',
            '</svg>',
        ]
    )
    path.with_suffix(".svg").write_text("\n".join(parts), encoding="utf-8")


def render_direction_percentiles(
    path: Path,
    rows: list[dict[str, float]],
    metric: str,
) -> None:
    is_fraction = metric == "fraction"
    rhos = [row["rho"] for row in rows]
    p10 = [row[f"{metric}_p10"] for row in rows]
    median = [row[f"{metric}_median"] for row in rows]
    p90 = [row[f"{metric}_p90"] for row in rows]

    width, height = 1080, 620
    left, right, top, bottom = 105, 965, 120, 520
    plot_width, plot_height = right - left, bottom - top
    rho_lo, rho_hi = min(rhos), max(rhos)

    def xp(rho: float) -> float:
        return left + plot_width * (
            (math.log10(rho) - math.log10(rho_lo))
            / (math.log10(rho_hi) - math.log10(rho_lo))
        )

    if is_fraction:
        y_lo, y_hi = 0.0, 1.0

        def yp(value: float) -> float:
            return bottom - plot_height * (value - y_lo) / (y_hi - y_lo)

        y_ticks = [(index / 10, f"{10 * index}%") for index in range(11)]
        reference, reference_label = 0.5, "Equal odd/even squared-RMSE share (50%)"
        title = "Direction dependence of odd squared-error share"
        subtitle = "P10–P90 across the same 100 directions; original scan plus 3 high-rho points"
        y_label = "Odd share  f = E_odd^2 / (E_odd^2 + E_even^2)"
        desc = "The median and tenth-to-ninetieth percentile band of the bounded odd squared-error share."
    else:
        y_lo = 10.0 ** math.floor(math.log10(min(p10)))
        y_hi = 10.0 ** math.ceil(math.log10(max(p90)))

        def yp(value: float) -> float:
            return bottom - plot_height * (
                (math.log10(value) - math.log10(y_lo))
                / (math.log10(y_hi) - math.log10(y_lo))
            )

        y_ticks = [
            (10.0**power, f"10^{power}")
            for power in range(
                int(math.log10(y_lo)), int(math.log10(y_hi)) + 1
            )
        ]
        reference, reference_label = 1.0, "Equal odd/even RMSE (r = 1)"
        title = "Direction dependence of odd/even RMSE ratio"
        subtitle = "P10–P90 across the same 100 directions; original scan plus 3 high-rho points"
        y_label = "Odd/even ratio  r = E_odd / E_even"
        desc = "The median and tenth-to-ninetieth percentile band of the unbounded odd-to-even RMSE ratio."

    def points(values: list[float]) -> str:
        return " ".join(
            f"{xp(rho):.2f},{yp(value):.2f}"
            for rho, value in zip(rhos, values)
        )

    upper_points = [(xp(rho), yp(value)) for rho, value in zip(rhos, p90)]
    lower_points = [(xp(rho), yp(value)) for rho, value in zip(rhos, p10)]
    band = " ".join(
        f"{x:.2f},{y:.2f}" for x, y in upper_points + list(reversed(lower_points))
    )

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="percentile-title percentile-desc">',
        f'<title id="percentile-title">{escape(title)}</title>',
        f'<desc id="percentile-desc">{escape(desc)}</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#202020">',
        f'<text x="540" y="34" text-anchor="middle" font-size="21" font-weight="600">{escape(title)}</text>',
        f'<text x="540" y="58" text-anchor="middle" font-size="13" fill="#555">{escape(subtitle)}</text>',
        '<rect x="284" y="78" width="36" height="12" fill="#009E73" fill-opacity="0.18"/>',
        '<line x1="284" y1="84" x2="320" y2="84" stroke="#007A59" stroke-width="1.8" stroke-dasharray="6 4"/>',
        '<text x="330" y="89" font-size="13">P10–P90 band</text>',
        '<line x1="488" y1="84" x2="524" y2="84" stroke="#007A59" stroke-width="2.8"/>',
        '<circle cx="506" cy="84" r="3.5" fill="white" stroke="#007A59" stroke-width="1.5"/>',
        '<text x="534" y="89" font-size="13">Median</text>',
        '<line x1="646" y1="84" x2="682" y2="84" stroke="#666" stroke-width="1.5" stroke-dasharray="7 5"/>',
        f'<text x="692" y="89" font-size="13">{escape(reference_label)}</text>',
    ]

    for value, label in y_ticks:
        y = yp(value)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#E1E4E6"/>')
        parts.append(f'<text x="{left-12}" y="{y+4:.2f}" text-anchor="end" font-size="11">{label}</text>')
    for tick in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0):
        if rho_lo <= tick <= rho_hi:
            x = xp(tick)
            parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#ECEEEF"/>')
            parts.append(f'<text x="{x:.2f}" y="{bottom+20}" text-anchor="middle" font-size="11">{tick:g}</text>')

    reference_y = yp(reference)
    parts.extend(
        [
            f'<polygon points="{band}" fill="#009E73" fill-opacity="0.18" stroke="none"/>',
            f'<line x1="{left}" y1="{reference_y:.2f}" x2="{right}" y2="{reference_y:.2f}" stroke="#666" stroke-width="1.5" stroke-dasharray="7 5"/>',
            f'<polyline points="{points(p10)}" fill="none" stroke="#007A59" stroke-width="1.7" stroke-dasharray="6 4"/>',
            f'<polyline points="{points(p90)}" fill="none" stroke="#007A59" stroke-width="1.7" stroke-dasharray="6 4"/>',
            f'<polyline points="{points(median)}" fill="none" stroke="#007A59" stroke-width="2.8"/>',
        ]
    )
    for rho, value in zip(rhos, median):
        parts.append(f'<circle cx="{xp(rho):.2f}" cy="{yp(value):.2f}" r="3.4" fill="white" stroke="#007A59" stroke-width="1.5"/>')

    parts.extend(
        [
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#555"/>',
            f'<text x="{(left+right)/2:.2f}" y="{bottom+52}" text-anchor="middle" font-size="13">Normalized vertical-corrector radius, rho (log scale)</text>',
            f'<text x="24" y="320" text-anchor="middle" font-size="13" transform="rotate(-90 24 320)">{escape(y_label)}</text>',
            '</g>',
            '</svg>',
        ]
    )
    path.with_suffix(".svg").write_text("\n".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=Path)
    parser.add_argument("--pairs", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir or args.summary.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = load_summary(args.summary)
    pairs = load_pairs(args.pairs or args.summary.parent / "vertical_parity_pairs.csv")
    percentiles = direction_percentiles(pairs)
    write_analysis(output_dir / "VERTICAL_PARITY_RESULTS.md", rows)
    render_crossover_loglog(output_dir / "vertical_parity_crossover_loglog", rows)
    render_linear_growth(output_dir / "vertical_parity_growth_linear", rows)
    write_direction_percentiles(
        output_dir / "vertical_parity_direction_percentiles.csv", percentiles
    )
    render_direction_percentiles(
        output_dir / "vertical_parity_odd_fraction_percentiles", percentiles, "fraction"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
