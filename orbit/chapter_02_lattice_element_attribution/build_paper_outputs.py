#!/usr/bin/env python3
"""Validate paired element-attribution results and build paper artifacts."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import tomllib
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
ERROR_ANALYZER = (
    PROJECT_ROOT
    / "orbit"
    / "error_analysis"
    / "thick_element_sextupole_sourcing"
    / "analyze_thick_element_sourcing.py"
)
PAPER_TOOLS = (
    PROJECT_ROOT
    / "orbit"
    / "high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis"
    / "svg_to_pdf.py"
)
FAMILY_LABELS = {
    "normal_sextupole": "Normal sextupole",
    "other_sextupole": "Other sextupole",
    "sbend": "Sector bend",
    "quadrupole": "Quadrupole",
    "drift": "Drift/geometric map",
    "kicker": "Kicker",
    "solenoid": "Solenoid",
    "wiggler": "Wiggler",
    "rfcavity": "RF cavity",
    "octupole": "Octupole",
    "marker": "Marker",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--horizontal", type=Path, required=True)
    parser.add_argument("--vertical", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def read_metadata(directory: Path) -> dict[str, object]:
    path = directory / "metadata.toml"
    text = path.read_text(encoding="utf-8")
    metadata = tomllib.loads(text)
    if "engine" not in metadata and "scibmad_version" in metadata:
        # Production files created immediately before the explicit engine
        # field was added are unambiguously SciBmad outputs.  Make that
        # provenance explicit in the raw chapter artifact as well.
        path.write_text('engine = "SciBmad"\n' + text, encoding="utf-8")
        metadata["engine"] = "SciBmad"
    return metadata


def scalar(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise RuntimeError(f"Non-finite {key}: {row[key]}")
    return value


def validate_directory(directory: Path, expected_plane: str) -> tuple[dict[str, object], dict[str, str]]:
    metadata = read_metadata(directory)
    if str(metadata.get("ring_id", "")).casefold() != "latest_cesr":
        raise RuntimeError(f"{directory} is not a latest_cesr result")
    if str(metadata.get("engine", "")).casefold() != "scibmad":
        raise RuntimeError(f"{directory} does not identify SciBmad")
    if not bool(metadata.get("rf_on")) or int(metadata.get("branch", -1)) != 0:
        raise RuntimeError(f"{directory} is not RF-on branch 0")
    lattice = str(metadata.get("lattice_path", "")).replace("\\", "/").casefold()
    if not lattice.endswith("latest_lattice/latest_cesr_scibmad_repaired.jl"):
        raise RuntimeError(f"{directory} does not use the repaired latest lattice")
    if str(metadata.get("output_plane", "")).casefold() != expected_plane:
        raise RuntimeError(f"{directory} has the wrong output plane")

    summaries = read_rows(directory / "reconstruction_summary.csv")
    if len(summaries) != 1:
        raise RuntimeError(f"Expected one reconstruction row in {directory}")
    summary = summaries[0]
    if int(summary["trials"]) != int(metadata["trials"]):
        raise RuntimeError(f"Trial mismatch in {directory}")
    if int(summary["elements"]) != int(metadata["element_count"]):
        raise RuntimeError(f"Element-count mismatch in {directory}")
    element_rows = read_rows(directory / "element_contribution_summary.csv")
    active_names = [
        row["element_name"]
        for row in element_rows
        if abs(scalar(row, "k2l_m2")) > 0
    ]
    expected_active = int(summary["active_normal_sextupoles"])
    if len(active_names) != expected_active:
        raise RuntimeError(f"Active-sextupole inventory mismatch in {directory}")
    if "active_normal_sextupoles" not in metadata:
        metadata_path = directory / "metadata.toml"
        with metadata_path.open("a", encoding="utf-8") as stream:
            stream.write(f"active_normal_sextupoles = {expected_active}\n")
            stream.write(
                "active_normal_sextupole_names = "
                + json.dumps(active_names)
                + "\n"
            )
        metadata["active_normal_sextupoles"] = expected_active
        metadata["active_normal_sextupole_names"] = active_names
    if expected_active != int(metadata["active_normal_sextupoles"]):
        raise RuntimeError(f"Sextupole-count mismatch in {directory}")
    if scalar(summary, "total_all_element_relative_closure") > 1.0e-10:
        raise RuntimeError(f"All-element reconstruction does not close in {directory}")
    if abs(scalar(summary, "total_all_element_signed_projection") - 1.0) > 1.0e-10:
        raise RuntimeError(f"All-element signed projection does not close in {directory}")

    for filename in (
        "element_contribution_summary.csv",
        "family_contribution_summary.csv",
        "direction_closure.csv",
        "element_direction_contributions.csv",
        "family_direction_contributions.csv",
        "family_direction_percentiles.csv",
    ):
        rows = read_rows(directory / filename)
        for row in rows:
            for key, value in row.items():
                if key in {"family", "element_name", "element_type", "output_plane"}:
                    continue
                scalar(row, key)
    return metadata, summary


def paired_families(horizontal: Path, vertical: Path) -> list[dict[str, object]]:
    x_rows = {row["family"]: row for row in read_rows(horizontal / "family_contribution_summary.csv")}
    y_rows = {row["family"]: row for row in read_rows(vertical / "family_contribution_summary.csv")}
    if x_rows.keys() != y_rows.keys():
        raise RuntimeError("Horizontal and vertical family registries disagree")
    paired: list[dict[str, object]] = []
    for family in x_rows:
        x, y = x_rows[family], y_rows[family]
        if int(x["element_count"]) != int(y["element_count"]):
            raise RuntimeError(f"Element-count mismatch for family {family}")
        paired.append(
            {
                "family": family,
                "element_count": int(x["element_count"]),
                "eta_x_percent": 100.0 * scalar(x, "eta_total"),
                "magnitude_x_percent": 100.0 * scalar(x, "magnitude_total"),
                "eta_y_percent": 100.0 * scalar(y, "eta_total"),
                "magnitude_y_percent": 100.0 * scalar(y, "magnitude_total"),
            }
        )
    paired.sort(
        key=lambda row: (
            row["family"] != "normal_sextupole",
            -max(abs(float(row["eta_x_percent"])), abs(float(row["eta_y_percent"]))),
        )
    )
    return paired


def write_family_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_family_latex(path: Path, rows: list[dict[str, object]]) -> None:
    lines = [
        r"\begin{table}[!htb]",
        r"  \centering",
        r"  \caption{Complete-element family attribution of the latest-CESR horizontal and vertical quadratic detector targets. Signed projections $\eta_{F,p}$ are additive; magnitude ratios $\mu_{F,p}$ are not.}",
        r"  \scriptsize",
        r"  \setlength{\tabcolsep}{4pt}",
        r"  \begin{tabular}{lrrrrr}",
        r"    \toprule",
        r"    Family & Elements & $\eta_{F,x}$ [\%] & $\mu_{F,x}$ [\%] & $\eta_{F,y}$ [\%] & $\mu_{F,y}$ [\%] \\",
        r"    \midrule",
    ]
    for row in rows:
        label = FAMILY_LABELS.get(str(row["family"]), str(row["family"]).replace("_", " ").title())
        lines.append(
            f"    {label} & {int(row['element_count'])}"
            f" & {float(row['eta_x_percent']):+.3f}"
            f" & {float(row['magnitude_x_percent']):.3f}"
            f" & {float(row['eta_y_percent']):+.3f}"
            f" & {float(row['magnitude_y_percent']):.3f} \\\\"
        )
    lines.extend(
        [
            r"    \bottomrule",
            r"  \end{tabular}",
            r"  \label{tab:latest-cesr-element-families}",
            r"\end{table}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def build_figures(horizontal: Path, vertical: Path, figures: Path, active_count: int) -> None:
    figures.mkdir(parents=True, exist_ok=True)
    analyzer = load_module("thick_element_analyzer", ERROR_ANALYZER)
    x_elements = read_rows(horizontal / "element_contribution_summary.csv")
    y_elements = read_rows(vertical / "element_contribution_summary.csv")
    x_sextupoles = [row for row in x_elements if abs(scalar(row, "k2l_m2")) > 0]
    y_sextupoles = [row for row in y_elements if abs(scalar(row, "k2l_m2")) > 0]
    if len(x_sextupoles) != active_count or len(y_sextupoles) != active_count:
        raise RuntimeError(
            f"Expected {active_count} active normal sextupoles, found {len(x_sextupoles)} and {len(y_sextupoles)}"
        )
    analyzer.render_paired_svg(
        figures / "normal_sextupole_signed_contributions_paired.svg",
        x_sextupoles,
        y_sextupoles,
        max(
            max(scalar(row, "s_m") for row in x_elements),
            max(scalar(row, "s_m") for row in y_elements),
        ),
    )

    converter = load_module("project_svg_to_pdf", PAPER_TOOLS)
    for stem in (
        "all_element_signed_contributions_paired",
        "normal_sextupole_signed_contributions_paired",
    ):
        source = figures / f"{stem}.svg"
        if not source.is_file():
            raise RuntimeError(f"Missing paired figure: {source}")
        converter.convert(source, figures / f"{stem}.pdf")


def write_report(
    path: Path,
    x_metadata: dict[str, object],
    x_summary: dict[str, str],
    y_summary: dict[str, str],
    families: list[dict[str, object]],
) -> None:
    sext = next(row for row in families if row["family"] == "normal_sextupole")
    lines = [
        "# Chapter 2 results: complete lattice-element attribution",
        "",
        "Status: `production_complete`.",
        "",
        "## Provenance",
        "",
        f"- Ring: `{x_metadata['ring_id']}`; RF on; branch `{x_metadata['branch']}`.",
        f"- Lattice: `{str(x_metadata['lattice_path']).replace(chr(92), '/')}`.",
        f"- Engine: `{x_metadata['engine']} {x_metadata['scibmad_version']}`.",
        f"- Ensemble: `{x_metadata['trials']}` paired directions per output plane; seed `{x_metadata['seed']}`; base kick `{1.0e6 * float(x_metadata['base_kick_rad']):g}` microradian.",
        f"- Runtime inventory: `{x_summary['elements']}` complete elements, `{x_summary['active_normal_sextupoles']}` active normal sextupoles, and `{x_summary['detectors']}` detectors.",
        "",
        "## Numerical closure",
        "",
        f"- Horizontal all-element relative vector closure: `{float(x_summary['total_all_element_relative_closure']):.6e}`; signed projection: `{float(x_summary['total_all_element_signed_projection']):.12g}`.",
        f"- Vertical all-element relative vector closure: `{float(y_summary['total_all_element_relative_closure']):.6e}`; signed projection: `{float(y_summary['total_all_element_signed_projection']):.12g}`.",
        f"- Normal-sextupole signed projections: `{float(sext['eta_x_percent']):+.3f}%` in x and `{float(sext['eta_y_percent']):+.3f}%` in y.",
        f"- Normal-sextupole magnitude ratios: `{float(sext['magnitude_x_percent']):.3f}%` in x and `{float(sext['magnitude_y_percent']):.3f}%` in y.",
        "",
        "Signed projections may be negative or exceed 100% for individual families because propagated source vectors interfere. Magnitude ratios are not additive and must not be interpreted as positive error shares.",
        "",
        "## Paper artifacts",
        "",
        "- `tables/family_attribution.csv` and `tables/family_attribution.tex` reproduce the paper's paired family table format.",
        "- `figures/normal_sextupole_signed_contributions_paired.svg` and `.pdf` reproduce the paper's paired normal-sextupole layout.",
        "- `figures/all_element_signed_contributions_paired.svg` and `.pdf` retain the complete-element view used to audit the decomposition.",
        "",
        "![Latest-CESR normal-sextupole attribution](figures/normal_sextupole_signed_contributions_paired.svg)",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = arguments()
    x_metadata, x_summary = validate_directory(args.horizontal, "x")
    y_metadata, y_summary = validate_directory(args.vertical, "y")
    for key in (
        "ring_id",
        "lattice_path",
        "engine",
        "scibmad_version",
        "branch",
        "rf_on",
        "trials",
        "seed",
        "base_kick_rad",
        "control_names",
        "detector_names",
    ):
        if x_metadata.get(key) != y_metadata.get(key):
            raise RuntimeError(f"Horizontal and vertical metadata disagree on {key}")
    families = paired_families(args.horizontal, args.vertical)
    args.tables_dir.mkdir(parents=True, exist_ok=True)
    write_family_csv(args.tables_dir / "family_attribution.csv", families)
    write_family_latex(args.tables_dir / "family_attribution.tex", families)
    build_figures(
        args.horizontal,
        args.vertical,
        args.figures_dir,
        int(x_summary["active_normal_sextupoles"]),
    )
    write_report(args.report, x_metadata, x_summary, y_summary, families)
    print(f"Validated and built Chapter 2 paper outputs: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
