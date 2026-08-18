#!/usr/bin/env python3
"""Render sparse sextupole--quadrupole affinity heatmaps from scored pairs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
LATEST_RESULTS = HERE / "results" / "scibmad_latest"


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--scores",
        type=Path,
        default=LATEST_RESULTS / "affinity" / "quadrupole_affinity_scores.csv",
    )
    result.add_argument(
        "--screen",
        type=Path,
        default=LATEST_RESULTS / "responses" / "quadrupole_optics_screen.csv",
    )
    result.add_argument("--output-dir", type=Path, default=LATEST_RESULTS / "affinity")
    return result


def quadrupole_key(name: str) -> tuple[int, str]:
    match = re.fullmatch(r"Q(\d+)([EW])", name.upper())
    if match is None:
        return (10_000, name)
    return (int(match.group(1)), match.group(2))


def finite_float(value: str, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {label}: {value}")
    return result


def render(
    rows: list[dict[str, str]],
    sextupoles: list[str],
    quadrupoles: list[str],
    metric: str,
    title: str,
    colorbar_label: str,
    stem: str,
    output_dir: Path,
) -> dict[str, Any]:
    s_index = {name: index for index, name in enumerate(sextupoles)}
    q_index = {name: index for index, name in enumerate(quadrupoles)}
    values = np.full((len(quadrupoles), len(sextupoles)), np.nan)
    for row in rows:
        values[q_index[row["quadrupole"]], s_index[row["sextupole"]]] = finite_float(
            row[metric], metric
        )
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError(f"No finite values for {metric}")

    cmap = plt.colormaps["viridis"].copy()
    cmap.set_bad("#f3f4f6")
    figure, axis = plt.subplots(figsize=(18.0, 22.0), constrained_layout=True)
    image = axis.imshow(values, aspect="auto", interpolation="nearest", cmap=cmap)
    axis.set_title(title, fontsize=17, pad=14)
    axis.set_xlabel(
        "Target sextupole (ring s order; blank cells were not evaluated)", fontsize=13
    )
    axis.set_ylabel("Quadrupole (ring s order)", fontsize=13)
    axis.set_xticks(np.arange(len(sextupoles)))
    axis.set_xticklabels(sextupoles, rotation=90, fontsize=7.0)
    axis.set_yticks(np.arange(len(quadrupoles)))
    axis.set_yticklabels(quadrupoles, fontsize=6.5)
    axis.tick_params(length=0)

    for column_index in range(values.shape[1]):
        if np.any(np.isfinite(values[:, column_index])):
            row_index = int(np.nanargmax(values[:, column_index]))
            axis.add_patch(
                Rectangle(
                    (column_index - 0.5, row_index - 0.5),
                    1.0,
                    1.0,
                    fill=False,
                    edgecolor="#111827",
                    linewidth=0.75,
                )
            )

    colorbar = figure.colorbar(image, ax=axis, shrink=0.82, pad=0.012)
    colorbar.set_label(colorbar_label, fontsize=11)
    output_dir.mkdir(parents=True, exist_ok=True)
    svg = output_dir / f"{stem}.svg"
    png = output_dir / f"{stem}.png"
    figure.savefig(svg, metadata={"Date": None})
    figure.savefig(png, dpi=180, metadata={"Software": "matplotlib"})
    plt.close(figure)
    return {
        "metric": metric,
        "minimum": float(np.min(finite)),
        "median": float(np.median(finite)),
        "maximum": float(np.max(finite)),
        "svg": svg.name,
        "png": png.name,
    }


def main() -> int:
    args = parser().parse_args()
    scores = args.scores.expanduser().resolve()
    screen = args.screen.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    with scores.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"No rows in {scores}")

    sextupole_s = {
        row["sextupole"]: finite_float(row["sextupole_s_m"], "sextupole_s_m") for row in rows
    }
    sextupoles = sorted(sextupole_s, key=lambda name: sextupole_s[name])
    with screen.open(encoding="utf-8", newline="") as stream:
        screen_rows = list(csv.DictReader(stream))
    retained_quadrupoles = {row["quadrupole"] for row in rows}
    if screen_rows and all(row.get("quadrupole_s_m", "").strip() for row in screen_rows):
        quadrupole_s = {
            row["quadrupole"]: finite_float(row["quadrupole_s_m"], "quadrupole_s_m")
            for row in screen_rows
            if row["quadrupole"] in retained_quadrupoles
        }
        quadrupoles = sorted(quadrupole_s, key=lambda name: quadrupole_s[name])
    else:
        quadrupoles = sorted(retained_quadrupoles, key=quadrupole_key)
    summaries = [
        render(
            rows,
            sextupoles,
            quadrupoles,
            "information_gain_logdet",
            "Nuisance-marginalized information gain",
            "log det information gain",
            "information_gain_heatmap",
            output_dir,
        ),
        render(
            rows,
            sextupoles,
            quadrupoles,
            "precision_improvement_worst_axis",
            "Nuisance-marginalized center precision improvement",
            "worst-axis precision ratio",
            "precision_improvement_heatmap",
            output_dir,
        ),
    ]
    (output_dir / "heatmap_summary.json").write_text(
        json.dumps(
            {
                "retained_pairs": len(rows),
                "sextupoles": len(sextupoles),
                "independent_quadrupoles_screened": len(
                    {row["quadrupole"] for row in screen_rows}
                ),
                "quadrupoles_shown": len(quadrupoles),
                "quadrupoles_with_at_least_one_retained_pair": len(
                    {row["quadrupole"] for row in rows}
                ),
                "blank_cells_are_not_evaluated": True,
                "plots": summaries,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    for summary in summaries:
        print(
            f"{summary['metric']}: {summary['minimum']:.6g} .. "
            f"{summary['maximum']:.6g}; {summary['svg']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
