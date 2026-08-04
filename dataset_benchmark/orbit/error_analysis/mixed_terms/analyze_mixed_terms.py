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

        combined_values = {
            "mixed_to_pure_ratio": [],
            "mixed_energy_share": [],
            "reconstruction_improvement": [],
            "relative_reconstruction_remainder": [],
        }
        pooled_mixed_squared = 0.0
        pooled_pure_squared = 0.0
        for row in selected:
            mixed_squared = sum(float(row[f"{plane}_qhv_rmse_m"]) ** 2 for plane in ("x", "y"))
            pure_squared = sum(float(row[f"{plane}_pure_rmse_m"]) ** 2 for plane in ("x", "y"))
            exact_squared = sum(float(row[f"{plane}_exact_residual_rmse_m"]) ** 2 for plane in ("x", "y"))
            pure_remainder_squared = sum(
                float(row[f"{plane}_pure_only_remainder_rmse_m"]) ** 2
                for plane in ("x", "y")
            )
            full_remainder_squared = sum(
                float(row[f"{plane}_full_reconstruction_remainder_rmse_m"]) ** 2
                for plane in ("x", "y")
            )
            combined_values["mixed_to_pure_ratio"].append(math.sqrt(mixed_squared / pure_squared))
            combined_values["mixed_energy_share"].append(mixed_squared / (mixed_squared + pure_squared))
            combined_values["reconstruction_improvement"].append(
                1.0 - full_remainder_squared / pure_remainder_squared
            )
            combined_values["relative_reconstruction_remainder"].append(
                math.sqrt(full_remainder_squared / exact_squared)
            )
            pooled_mixed_squared += mixed_squared
            pooled_pure_squared += pure_squared
        for metric, values in combined_values.items():
            for label, probability in (("p10", 0.1), ("median", 0.5), ("p90", 0.9)):
                output[f"combined_{metric}_{label}"] = percentile(values, probability)
        output["combined_pooled_mixed_energy_share"] = (
            pooled_mixed_squared / (pooled_mixed_squared + pooled_pure_squared)
        )
        result.append(output)
    return result


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, summary, percentiles) -> None:
    rhos = [float(row["rho"]) for row in summary]
    slopes = {
        plane: local_slopes(rhos, [float(row[f"mean_{plane}_qhv_rmse_m"]) for row in summary])
        for plane in ("x", "y")
    }
    lines = [
        "# All-corrector horizontal--vertical mixed-term experiment", "",
        "The four-sign finite difference directly separates the pure horizontal",
        "`Q_hh`, pure vertical `Q_vv`, and mixed `Q_hv` response vectors. All",
        "reconstructions and remainders are computed as vectors before detector RMS.", "",
        "## Mean decomposition", "",
        "| rho | H/V kick RMS (urad) | plane | Q_hh (um) | Q_vv (um) | Q_hv (um) | Q_hv slope | mixed energy share | full remainder / exact |",
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
        f"At `rho = {rhos[-1]:g}` ({5*rhos[-1]:.3g} urad RMS in each family):", "",
        "| plane | mixed/pure P10 | median | P90 | mixed-energy share P10 | median | P90 | reconstruction improvement P10 | median | P90 |",
        "|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for plane in ("x", "y", "combined"):
        label = "X+Y" if plane == "combined" else plane.upper()
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

    y_coefficients = [
        1e6 * float(row["mean_y_qhv_rmse_m"]) / float(row["rho"]) ** 2
        for row in summary
    ]
    lines.extend([
        "", "## Compact checks", "",
        f"- Mean Y `Q_hv/rho^2` changes by `{100*(max(y_coefficients)/min(y_coefficients)-1):.5g}%` over the fitted interval.",
        f"- At the smallest radius, mean Y mixed energy share is `{100*float(first['mean_y_mixed_energy_share']):.6g}%`.",
        f"- At the largest radius, mean Y mixed energy share is `{100*float(last['mean_y_mixed_energy_share']):.6g}%`.",
        f"- At the largest radius, the full signed reconstruction leaves `{100*float(last['mean_y_relative_reconstruction_remainder']):.6g}%` of the exact Y residual RMS.",
        f"- At the largest radius, adding `Q_hv` reduces the mean squared Y reconstruction error by `{100*float(last['mean_y_reconstruction_improvement']):.6g}%` relative to the pure-only reconstruction.",
        f"- For the combined 198-dimensional X+Y vector at the largest radius, the pooled mixed squared-norm share is `{100*float(last_pct['combined_pooled_mixed_energy_share']):.6g}%`; its direction-resolved P10/median/P90 shares are `{100*float(last_pct['combined_mixed_energy_share_p10']):.4g}% / {100*float(last_pct['combined_mixed_energy_share_median']):.4g}% / {100*float(last_pct['combined_mixed_energy_share_p90']):.4g}%`.",
        "", "The cross-term hypothesis is supported only if `Q_hv` remains quadratic,",
        "dominates the Y direction distribution, and the signed reconstruction closes",
        "the four exact joint responses with a much smaller higher-order remainder.", "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def render_svg(path: Path, summary) -> None:
    rhos = [float(row["rho"]) for row in summary]
    width, height = 1040, 500
    left, right, top, bottom = 90, 970, 85, 410
    colors = {"qhh": "#2369A7", "qvv": "#D97904", "qhv": "#8B3FA6"}
    values = [1e6*float(row[f"mean_y_{term}_rmse_m"]) for row in summary for term in colors]
    xlo, xhi = math.log10(min(rhos)), math.log10(max(rhos))
    ylo, yhi = math.floor(math.log10(min(values))), math.ceil(math.log10(max(values)))

    def xp(value: float) -> float:
        return left + (right-left)*(math.log10(value)-xlo)/(xhi-xlo)

    def yp(value: float) -> float:
        return bottom - (bottom-top)*(math.log10(value)-ylo)/(yhi-ylo)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="520" y="32" text-anchor="middle" font-size="21" font-weight="600">All-corrector Y quadratic block decomposition</text>',
        '<text x="520" y="55" text-anchor="middle" font-size="13" fill="#555">100 fixed H/V direction pairs; vector RMS at 99 detectors</text>',
    ]
    for exponent in range(ylo, yhi + 1):
        y = yp(10.0**exponent)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#ddd"/>')
        parts.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end" font-size="12">10^{exponent}</text>')
    for rho in rhos:
        x = xp(rho)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#eee"/>')
        parts.append(f'<text x="{x:.2f}" y="{bottom+22}" text-anchor="middle" font-size="11">{rho:g}</text>')
    for term, color in colors.items():
        term_values = [1e6*float(row[f"mean_y_{term}_rmse_m"]) for row in summary]
        points = " ".join(f"{xp(rho):.2f},{yp(value):.2f}" for rho, value in zip(rhos, term_values))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.3"/>')
        for rho, value in zip(rhos, term_values):
            parts.append(f'<circle cx="{xp(rho):.2f}" cy="{yp(value):.2f}" r="3.2" fill="white" stroke="{color}" stroke-width="1.6"/>')
    for index, (term, label) in enumerate((("qhh", "Q_hh"), ("qvv", "Q_vv"), ("qhv", "Q_hv"))):
        x = 350 + 130*index
        parts.append(f'<line x1="{x}" y1="75" x2="{x+25}" y2="75" stroke="{colors[term]}" stroke-width="2.3"/>')
        parts.append(f'<text x="{x+31}" y="79" font-size="12">{label}</text>')
    parts.extend([
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#333"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#333"/>',
        f'<text x="{(left+right)/2}" y="{bottom+55}" text-anchor="middle" font-size="13">Normalized radius rho (log scale)</text>',
        f'<text x="24" y="{(top+bottom)/2}" text-anchor="middle" font-size="13" transform="rotate(-90 24 {(top+bottom)/2})">Mean detector RMSE (um, log scale)</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).parent / "results")
    args = parser.parse_args()
    summary = load_csv(args.results_dir / "mixed_term_summary.csv")
    directions = load_csv(args.results_dir / "mixed_term_directions.csv")
    percentiles = direction_percentiles(directions)
    write_csv(args.results_dir / "mixed_term_direction_percentiles.csv", percentiles)
    write_report(args.results_dir / "MIXED_TERM_RESULTS.md", summary, percentiles)
    render_svg(args.results_dir / "mixed_term_decomposition.svg", summary)


if __name__ == "__main__":
    main()
