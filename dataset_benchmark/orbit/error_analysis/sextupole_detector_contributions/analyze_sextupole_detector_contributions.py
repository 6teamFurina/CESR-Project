#!/usr/bin/env python3
"""Summarize and render the signed sextupole detector contributions."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def value(row: dict[str, str], key: str) -> float:
    return float(row[key])


def escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def render_svg(path: Path, rows: list[dict[str, str]]) -> None:
    width, height = 1120, 560
    left, right, top, bottom = 84, 28, 70, 92
    chart_w, chart_h = width - left - right, height - top - bottom
    circumference = max(value(row, "s_m") for row in rows)
    eta = [100.0 * value(row, "eta_total") for row in rows]
    limit = 1.05 * max(abs(item) for item in eta) if any(eta) else 1.0
    zero_y = top + chart_h / 2

    def x_pos(s_m: float) -> float:
        return left + chart_w * s_m / circumference

    def y_pos(item: float) -> float:
        return zero_y - (chart_h / 2) * item / limit

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:21px;font-weight:600}.sub{font-size:13px;fill:#555}.tick{font-size:12px}.label{font-size:11px}</style>',
        f'<text x="{width/2}" y="29" text-anchor="middle" class="title">Signed normal-sextupole reconstruction of horizontal quadratic detector error</text>',
        f'<text x="{width/2}" y="51" text-anchor="middle" class="sub">Ensemble projection onto Qx = Qhh,x + Qhv,x + Qvv,x; positive reinforces, negative cancels</text>',
    ]
    for tick in range(-4, 5):
        item = limit * tick / 4
        y = y_pos(item)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{left-9}" y="{y+4:.2f}" text-anchor="end" class="tick">{item:.2g}</text>')
    parts.append(f'<line x1="{left}" y1="{zero_y:.2f}" x2="{width-right}" y2="{zero_y:.2f}" stroke="#555" stroke-width="1.2"/>')
    bar_w = max(3.0, 0.58 * chart_w / len(rows))
    for row, item in zip(rows, eta):
        x, y = x_pos(value(row, "s_m")) - bar_w / 2, y_pos(item)
        color = "#d95f02" if item >= 0 else "#1f78b4"
        parts.append(f'<rect x="{x:.2f}" y="{min(y, zero_y):.2f}" width="{bar_w:.2f}" height="{abs(zero_y-y):.2f}" fill="{color}" opacity="0.88"/>')
    for tick in range(9):
        s_m = circumference * tick / 8
        x = x_pos(s_m)
        parts.append(f'<line x1="{x:.2f}" y1="{top+chart_h}" x2="{x:.2f}" y2="{top+chart_h+5}" stroke="#444"/>')
        parts.append(f'<text x="{x:.2f}" y="{top+chart_h+22}" text-anchor="middle" class="tick">{s_m:.0f}</text>')
    ranked = sorted(rows, key=lambda row: abs(value(row, "eta_total")), reverse=True)[:10]
    for rank, row in enumerate(ranked, 1):
        item = 100 * value(row, "eta_total")
        x, y = x_pos(value(row, "s_m")), y_pos(item)
        label_y = top + 15 + 16 * ((rank - 1) % 5)
        label_x = left + 12 + (rank > 5) * (chart_w / 2)
        parts.append(f'<line x1="{x:.2f}" y1="{y:.2f}" x2="{label_x+92:.2f}" y2="{label_y-4:.2f}" stroke="#999" stroke-width="0.7"/>')
        parts.append(f'<text x="{label_x:.2f}" y="{label_y:.2f}" class="label">{rank}. {escape(row["element_name"])} ({item:+.2f}%)</text>')
    parts.extend([
        f'<text x="{left+chart_w/2}" y="{height-31}" text-anchor="middle" class="tick">Ring position s [m]</text>',
        f'<text x="22" y="{top+chart_h/2}" text-anchor="middle" class="tick" transform="rotate(-90 22 {top+chart_h/2})">signed projection eta_j [%]</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_report(
    path: Path,
    rows: list[dict[str, str]],
    direction_rows: list[dict[str, str]],
    summary: dict[str, str],
) -> None:
    ranked = sorted(rows, key=lambda row: abs(value(row, "eta_total")), reverse=True)
    closures = [value(row, "total_relative_closure") for row in direction_rows]
    projections = [value(row, "total_signed_projection") for row in direction_rows]
    total_closure = float(summary["total_concatenated_relative_closure"])
    lines = [
        "# Signed normal-sextupole detector-contribution result", "",
        "Positive signed projection reinforces the target and negative projection cancels it.", "",
        "## Reconstruction closure", "",
        f'- Directions: `{summary["trials"]}`.',
        f'- Active normal sextupoles: `{summary["active_normal_sextupoles"]}`.',
        f'- HH concatenated relative closure: `{float(summary["hh_concatenated_relative_closure"]):.6g}`.',
        f'- HV concatenated relative closure: `{float(summary["hv_concatenated_relative_closure"]):.6g}`.',
        f'- VV concatenated relative closure: `{float(summary["vv_concatenated_relative_closure"]):.6g}`.',
        f'- Total concatenated relative closure: `{float(summary["total_concatenated_relative_closure"]):.6g}`.',
        f'- Total reconstruction signed projection: `{float(summary["total_reconstruction_signed_projection"]):.6g}`.',
        f'- Direction-level total closure P10 / median / P90: `{percentile(closures, 0.1):.6g} / {percentile(closures, 0.5):.6g} / {percentile(closures, 0.9):.6g}`.',
        f'- Direction-level signed projection P10 / median / P90: `{percentile(projections, 0.1):.6g} / {percentile(projections, 0.5):.6g} / {percentile(projections, 0.9):.6g}`.',
        "",
        "## Interpretation", "",
        (
            f"The total concatenated residual is {100 * total_closure:.2f}%, so this "
            "calculation does not close tightly enough to be presented as a complete "
            "per-sextupole attribution of the final detector error. The ranking below "
            "is a signed projection of the leading thin normal-sextupole reconstruction; "
            "the residual must remain explicit."
        ),
        "", "## Largest absolute signed projections", "",
        "| rank | sextupole | s [m] | K2L [m^-2] | eta total | magnitude ratio |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for rank, row in enumerate(ranked[:15], 1):
        lines.append(f'| {rank} | `{row["element_name"]}` | {value(row,"s_m"):.3f} | {value(row,"k2l_m2"):.5g} | {100*value(row,"eta_total"):+.4f}% | {100*value(row,"magnitude_total"):.4f}% |')
    lines.extend(["", "## Evidence boundary", "",
        "The signed projection is additive only to the extent that the propagated local-source vectors reconstruct the GTPSA target. The magnitude ratio is not additive. If closure is not small, retain this as a leading thin-kick normal-sextupole reconstruction and keep the residual explicit.", "",
        "![Signed sextupole contributions](sextupole_signed_detector_contributions.svg)", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", type=Path, default=Path(__file__).resolve().parent / "results" / "latest_cesr")
    args = parser.parse_args()
    rows = read_rows(args.output_dir / "sextupole_contribution_summary.csv")
    direction_rows = read_rows(args.output_dir / "direction_closure.csv")
    summaries = read_rows(args.output_dir / "reconstruction_summary.csv")
    if len(summaries) != 1:
        raise RuntimeError("Expected exactly one reconstruction summary row")
    render_svg(args.output_dir / "sextupole_signed_detector_contributions.svg", rows)
    write_report(args.output_dir / "RESULTS.md", rows, direction_rows, summaries[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
