#!/usr/bin/env python3
"""Analyze the all-corrector H/V mixed-term finite-difference experiment."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_csv(path: Path) -> list[dict[str, float | bool]]:
    with path.open(newline="", encoding="utf-8") as stream:
        source = list(csv.DictReader(stream))
    rows: list[dict[str, float | bool]] = []
    for source_row in source:
        row: dict[str, float | bool] = {}
        for key, value in source_row.items():
            if value.lower() in ("true", "false"):
                row[key] = value.lower() == "true"
            else:
                row[key] = float(value)
        rows.append(row)
    return rows


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def local_slopes(rhos: list[float], values: list[float]) -> list[float]:
    result = [math.nan]
    for index in range(1, len(rhos)):
        result.append(
            math.log(values[index] / values[index - 1])
            / math.log(rhos[index] / rhos[index - 1])
        )
    return result


def direction_percentiles(rows: list[dict[str, float | bool]]) -> list[dict[str, float]]:
    complete = [row for row in rows if row["all_states_converged"]]
    rhos = sorted({float(row["rho"]) for row in complete})
    result = []
    for rho in rhos:
        selected = [row for row in complete if math.isclose(float(row["rho"]), rho)]
        output: dict[str, float] = {"rho": rho, "directions": float(len(selected))}
        for plane in ("x", "y"):
            for metric in (
                "mixed_to_pure_ratio", "mixed_energy_share",
                "reconstruction_improvement", "relative_reconstruction_remainder",
            ):
                values = [float(row[f"{plane}_{metric}"]) for row in selected]
                for label, probability in (("p10", 0.1), ("median", 0.5), ("p90", 0.9)):
                    output[f"{plane}_{metric}_{label}"] = percentile(values, probability)

        result.append(output)
    return result


def component_share_percentiles(rows: list[dict[str, float | bool]]) -> list[dict[str, float]]:
    """Percentiles of ||Q_block||^2 / sum_block ||Q_block||^2 per orbit plane."""
    complete = [row for row in rows if row["all_states_converged"]]
    rhos = sorted({float(row["rho"]) for row in complete})
    result = []
    for rho in rhos:
        selected = [row for row in complete if math.isclose(float(row["rho"]), rho)]
        output: dict[str, float] = {"rho": rho, "directions": float(len(selected))}
        for plane in ("x", "y"):
            shares = {term: [] for term in ("hh", "hv", "vv")}
            for row in selected:
                squared = {
                    term: float(row[f"{plane}_q{term}_rmse_m"]) ** 2
                    for term in ("hh", "hv", "vv")
                }
                denominator = sum(squared.values())
                for term, value in squared.items():
                    shares[term].append(value / denominator)
            for term, values in shares.items():
                for label, probability in (("p10", 0.1), ("median", 0.5), ("p90", 0.9)):
                    output[f"{plane}_{term}_{label}"] = percentile(values, probability)
        result.append(output)
    return result


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


SHARE_DIGITS = {
    ("x", "hh"): 2,
    ("x", "hv"): 4,
    ("x", "vv"): 2,
    ("y", "hh"): 4,
    ("y", "hv"): 3,
    ("y", "vv"): 4,
}


def share_cell(row: dict[str, float], plane: str, term: str) -> str:
    digits = SHARE_DIGITS[(plane, term)]
    values = [100 * row[f"{plane}_{term}_{label}"] for label in ("median", "p10", "p90")]
    return f"{values[0]:.{digits}f} [{values[1]:.{digits}f}, {values[2]:.{digits}f}]"


def write_component_share_tex(path: Path, rows: list[dict[str, float]]) -> None:
    lines = [
        r"\begin{table*}[!t]",
        r"  \centering",
        r"  \caption{Four-sign finite-difference validation of the GTPSA direction-resolved quadratic-block shares. Each entry is median [P10, P90] in percent across the same 100 fixed H/V direction pairs.}",
        r"  \resizebox{\textwidth}{!}{%",
        r"  \begin{tabular}{@{}crrrrrr@{}}",
        r"    \toprule",
        r"    & \multicolumn{3}{c}{\textbf{X orbit component}} & \multicolumn{3}{c}{\textbf{Y orbit component}} \\",
        r"    $\rho$ & $hh$ & $hv$ & $vv$ & $hh$ & $hv$ & $vv$ \\",
        r"    \midrule",
    ]
    for row in rows:
        cells = [share_cell(row, plane, term) for plane in ("x", "y") for term in ("hh", "hv", "vv")]
        lines.append(f"    {row['rho']:g} & " + " & ".join(cells) + r" \\")
    lines.extend([
        r"    \bottomrule",
        r"  \end{tabular}%",
        r"  }",
        r"  \label{tab:mixedshares}",
        r"\end{table*}",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_report(path: Path, summary, percentiles, component_shares) -> None:
    rhos = [float(row["rho"]) for row in summary]
    slopes = {
        plane: local_slopes(rhos, [float(row[f"mean_{plane}_qhv_rmse_m"]) for row in summary])
        for plane in ("x", "y")
    }
    lines = [
        "# Four-sign validation of the GTPSA mixed-term response", "",
        "The GTPSA direction-contracted values in",
        "`../gtpsa_results/GTPSA_RESULTS.md` are the adopted final results.",
        "This report preserves the independent finite-difference validation and",
        "its exact nonlinear reconstruction checks.", "",
        "The four-sign finite difference directly separates the pure horizontal",
        "`Q_hh`, pure vertical `Q_vv`, and mixed `Q_hv` response vectors. All",
        "reconstructions and remainders are computed as vectors before detector RMS.", "",
        "## Mean decomposition", "",
        "| ρ | H/V kick RMS (µrad) | plane | Q_hh (µm) | Q_vv (µm) | Q_hv (µm) | Q_hv slope | mixed energy share | full remainder / exact |",
        "|---:|---:|:---:|---:|---:|---:|---:|---:|---:|",
    ]
    for index, (rho, row) in enumerate(zip(rhos, summary)):
        for plane in ("x", "y"):
            slope = slopes[plane][index]
            slope_text = "n/a" if math.isnan(slope) else f"{slope:.4f}"
            lines.append(
                f"| {rho:g} | {5*rho:.3g} | {plane.upper()} | "
                f"{1e6*float(row[f'mean_{plane}_qhh_rmse_m']):.6g} | "
                f"{1e6*float(row[f'mean_{plane}_qvv_rmse_m']):.6g} | "
                f"{1e6*float(row[f'mean_{plane}_qhv_rmse_m']):.6g} | {slope_text} | "
                f"{100*float(row[f'mean_{plane}_mixed_energy_share']):.4f}% | "
                f"{100*float(row[f'mean_{plane}_relative_reconstruction_remainder']):.4f}% |"
            )

    last, last_pct, first = summary[-1], percentiles[-1], summary[0]
    lines.extend([
        "", "## Direction-resolved result at the largest fitted radius", "",
        f"At `ρ = {rhos[-1]:g}` ({5*rhos[-1]:.3g} µrad RMS in each family):", "",
        "| plane | mixed/pure P10 | median | P90 | mixed-energy share P10 | median | P90 | reconstruction improvement P10 | median | P90 |",
        "|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for plane in ("x", "y"):
        label = plane.upper()
        lines.append(
            f"| {label} | {float(last_pct[f'{plane}_mixed_to_pure_ratio_p10']):.4g} | "
            f"{float(last_pct[f'{plane}_mixed_to_pure_ratio_median']):.4g} | "
            f"{float(last_pct[f'{plane}_mixed_to_pure_ratio_p90']):.4g} | "
            f"{100*float(last_pct[f'{plane}_mixed_energy_share_p10']):.3f}% | "
            f"{100*float(last_pct[f'{plane}_mixed_energy_share_median']):.3f}% | "
            f"{100*float(last_pct[f'{plane}_mixed_energy_share_p90']):.3f}% | "
            f"{100*float(last_pct[f'{plane}_reconstruction_improvement_p10']):.3f}% | "
            f"{100*float(last_pct[f'{plane}_reconstruction_improvement_median']):.3f}% | "
            f"{100*float(last_pct[f'{plane}_reconstruction_improvement_p90']):.3f}% |"
        )

    lines.extend([
        "", "## Direction-resolved block squared-norm shares", "",
        "For each orbit plane, the three shares use the common denominator",
        "`||Q_hh||² + ||Q_hv||² + ||Q_vv||²` and therefore sum to one for",
        "every direction before percentiles are taken. Entries are",
        "`median [P10, P90]` in percent.", "",
        "| ρ | X hh | X hv | X vv | Y hh | Y hv | Y vv |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in component_shares:
        cells = [share_cell(row, plane, term) for plane in ("x", "y") for term in ("hh", "hv", "vv")]
        lines.append(f"| {row['rho']:g} | " + " | ".join(cells) + " |")

    y_coefficients = [
        1e6 * float(row["mean_y_qhv_rmse_m"]) / float(row["rho"]) ** 2
        for row in summary
    ]
    lines.extend([
        "", "## Compact checks", "",
        f"- Mean Y `Q_hv/ρ²` changes by `{100*(max(y_coefficients)/min(y_coefficients)-1):.5g}%` over the fitted interval.",
        f"- At the smallest radius, mean Y mixed energy share is `{100*float(first['mean_y_mixed_energy_share']):.6g}%`.",
        f"- At the largest radius, mean Y mixed energy share is `{100*float(last['mean_y_mixed_energy_share']):.6g}%`.",
        f"- At the largest radius, the full signed reconstruction leaves `{100*float(last['mean_y_relative_reconstruction_remainder']):.6g}%` of the exact Y residual RMS.",
        f"- At the largest radius, adding `Q_hv` reduces the mean squared Y reconstruction error by `{100*float(last['mean_y_reconstruction_improvement']):.6g}%` relative to the pure-only reconstruction.",
        "", "The cross-term hypothesis is supported only if `Q_hv` remains quadratic,",
        "dominates the Y direction distribution, and the signed reconstruction closes",
        "the four exact joint responses with a much smaller higher-order remainder.", "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def render_component_share_svg(path: Path, rows) -> None:
    rhos = [float(row["rho"]) for row in rows]
    width, height = 1160, 560
    top, bottom, panel_width, gap = 125, 445, 470, 95
    lefts = (85, 85 + panel_width + gap)
    colors = {"hh": "#0072B2", "hv": "#009E73", "vv": "#D55E00"}
    titles = {"x": "X orbit component: vv-dominated", "y": "Y orbit component: hv-dominated"}
    xlo, xhi = math.log10(min(rhos)), math.log10(max(rhos))

    def xp(value: float, left: float) -> float:
        return left + 4 + (panel_width - 8) * (math.log10(value) - xlo) / (xhi - xlo)

    def yp(value: float) -> float:
        return bottom - 4 - (bottom - top - 8) * value

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="share-title share-desc">',
        '<title id="share-title">Quadratic orbit-response block squared-norm shares</title>',
        '<desc id="share-desc">Two side-by-side panels show median and P10 to P90 bands across 100 directions. The X orbit component is dominated by vv and the Y orbit component by hv.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="580" y="32" text-anchor="middle" font-size="21" font-weight="600">Quadratic block squared-norm shares</text>',
        '<text x="580" y="57" text-anchor="middle" font-size="13" fill="#555">Median and P10–P90 across 100 fixed H/V direction pairs</text>',
    ]
    for index, term in enumerate(("hh", "hv", "vv")):
        x = 415 + 115 * index
        parts.append(f'<line x1="{x}" y1="87" x2="{x+28}" y2="87" stroke="{colors[term]}" stroke-width="2.7"/>')
        parts.append(f'<text x="{x+35}" y="91" font-size="13">{term}</text>')

    tick_rhos = (0.1, 0.2, 0.4, 0.8, 1.13)
    for plane, left in zip(("x", "y"), lefts):
        right = left + panel_width
        parts.append(f'<text x="{(left+right)/2:.2f}" y="112" text-anchor="middle" font-size="15" font-weight="600">{titles[plane]}</text>')
        for percentage in range(0, 101, 20):
            y = yp(percentage / 100)
            parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#dddddd" stroke-width="1"/>')
            parts.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" font-size="11">{percentage}</text>')
        for rho in tick_rhos:
            x = xp(rho, left)
            parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#eeeeee" stroke-width="1"/>')
            parts.append(f'<text x="{x:.2f}" y="{bottom+21}" text-anchor="middle" font-size="11">{rho:g}</text>')

        for term in ("hh", "hv", "vv"):
            p10 = [float(row[f"{plane}_{term}_p10"]) for row in rows]
            median = [float(row[f"{plane}_{term}_median"]) for row in rows]
            p90 = [float(row[f"{plane}_{term}_p90"]) for row in rows]
            upper = [(xp(rho, left), yp(value)) for rho, value in zip(rhos, p90)]
            lower = [(xp(rho, left), yp(value)) for rho, value in zip(rhos, p10)]
            band = " ".join(f"{x:.2f},{y:.2f}" for x, y in upper + list(reversed(lower)))
            line = " ".join(f"{xp(rho,left):.2f},{yp(value):.2f}" for rho, value in zip(rhos, median))
            parts.append(f'<polygon points="{band}" fill="{colors[term]}" fill-opacity="0.14" stroke="none"/>')
            parts.append(f'<polyline points="{line}" fill="none" stroke="{colors[term]}" stroke-width="2.7"/>')
            for rho, value in zip(rhos, median):
                parts.append(f'<circle cx="{xp(rho,left):.2f}" cy="{yp(value):.2f}" r="3.0" fill="white" stroke="{colors[term]}" stroke-width="1.5"/>')
        parts.append(f'<rect x="{left}" y="{top}" width="{panel_width}" height="{bottom-top}" fill="none" stroke="#555" stroke-width="1"/>')

    parts.extend([
        f'<text x="580" y="{bottom+61}" text-anchor="middle" font-size="13.5">Normalized active-corrector input radius, ρ</text>',
        f'<text x="24" y="{(top+bottom)/2}" text-anchor="middle" font-size="13.5" transform="rotate(-90 24 {(top+bottom)/2})">Squared-norm share (%)</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).parent / "results")
    args = parser.parse_args()
    summary = load_csv(args.results_dir / "mixed_term_summary.csv")
    directions = load_csv(args.results_dir / "mixed_term_directions.csv")
    percentiles = direction_percentiles(directions)
    component_shares = component_share_percentiles(directions)
    write_csv(args.results_dir / "mixed_term_direction_percentiles.csv", percentiles)
    write_csv(args.results_dir / "mixed_term_component_share_percentiles.csv", component_shares)
    write_component_share_tex(args.results_dir / "mixed_term_component_share_table.tex", component_shares)
    write_report(args.results_dir / "MIXED_TERM_RESULTS.md", summary, percentiles, component_shares)


if __name__ == "__main__":
    main()
