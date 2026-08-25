#!/usr/bin/env python3
"""Validate the production rho sweep and build paper-format artifacts."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
PAPER_TOOLS = (
    PROJECT_ROOT
    / "orbit"
    / "high_throughput_nonlinear_closed_orbit_calculation_and_response_error_analysis"
    / "svg_to_pdf.py"
)
PAPER_RADII = (1.13, 3.2, 4.53, 6.4, 9.05)
SCENARIOS = ("all", "horizontal", "vertical")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def finite(value: object, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise RuntimeError(f"Non-finite {label}: {value}")
    return number


def validate(rows: list[dict[str, str]], metadata: dict[str, object]) -> None:
    if str(metadata.get("ring_id", "")).casefold() != "latest_cesr":
        raise RuntimeError("Production metadata is not for latest_cesr")
    if str(metadata.get("engine", "")).casefold() != "scibmad":
        raise RuntimeError("Production metadata does not identify SciBmad")
    if not bool(metadata.get("rf_on")) or int(metadata.get("branch", -1)) != 0:
        raise RuntimeError("Production result is not RF-on branch 0")
    lattice = str(metadata.get("lattice_path", "")).replace("\\", "/").casefold()
    if not lattice.endswith("latest_lattice/latest_cesr_scibmad_repaired.jl"):
        raise RuntimeError("Production metadata does not use the repaired latest lattice")
    expected_trials = int(metadata["trials_per_positive_rho_scenario"])
    if expected_trials < 2:
        raise RuntimeError("Production direction count is not valid")
    for row in rows:
        scenario = row["scenario"]
        if scenario not in SCENARIOS:
            raise RuntimeError(f"Unexpected scenario: {scenario}")
        rho = finite(row["rho"], "rho")
        trials = int(row["trials"])
        converged = int(row["converged_trials"])
        if rho > 0 and trials != expected_trials:
            raise RuntimeError(
                f"{scenario} rho={rho:g} has {trials} trials, expected {expected_trials}"
            )
        if converged > trials:
            raise RuntimeError("Converged count exceeds trial count")
        for key, value in row.items():
            if key in {"scenario", "source_chunk"}:
                continue
            finite(value, f"{scenario}/{rho:g}/{key}")


def selected_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    for scenario in SCENARIOS:
        for target in PAPER_RADII:
            matches = [
                row
                for row in rows
                if row["scenario"] == scenario
                and math.isclose(float(row["rho"]), target, rel_tol=0.0, abs_tol=1.0e-12)
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    f"Expected one {scenario} row at rho={target:g}, found {len(matches)}"
                )
            selected.append(matches[0])
    return selected


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = (
        "scenario",
        "rho",
        "trials",
        "converged_trials",
        "mean_x_rmse_um",
        "max_x_rmse_um",
        "mean_y_rmse_um",
        "max_y_rmse_um",
        "max_closure_norm",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "scenario": row["scenario"],
                    "rho": f"{float(row['rho']):g}",
                    "trials": row["trials"],
                    "converged_trials": row["converged_trials"],
                    "mean_x_rmse_um": f"{1.0e6 * float(row['mean_x_rmse_m']):.9g}",
                    "max_x_rmse_um": f"{1.0e6 * float(row['max_trial_x_rmse_m']):.9g}",
                    "mean_y_rmse_um": f"{1.0e6 * float(row['mean_y_rmse_m']):.9g}",
                    "max_y_rmse_um": f"{1.0e6 * float(row['max_trial_y_rmse_m']):.9g}",
                    "max_closure_norm": f"{float(row['max_closure_norm']):.9g}",
                }
            )


def write_latex(path: Path, rows: list[dict[str, str]]) -> None:
    labels = {"all": "All", "horizontal": "Horizontal", "vertical": "Vertical"}
    lines = [
        r"\begin{table*}[!t]",
        r"  \centering",
        r"  \caption{Latest-CESR nonlinear orbit-response residuals at the five radii used by the matched paper benchmark. Mean and maximum values are over the fixed direction ensemble.}",
        r"  \scriptsize",
        r"  \setlength{\tabcolsep}{3.2pt}",
        r"  \begin{tabular}{lrrrrrr}",
        r"    \toprule",
        r"    Input & $\rho$ & Conv. & $\bar e_x$ [$\mu$m] & $e_{x,\max}$ [$\mu$m] & $\bar e_y$ [$\mu$m] & $e_{y,\max}$ [$\mu$m] \\",
        r"    \midrule",
    ]
    previous = None
    for row in rows:
        scenario = row["scenario"]
        label = labels[scenario] if scenario != previous else ""
        previous = scenario
        lines.append(
            "    "
            + f"{label} & {float(row['rho']):g} & {int(row['converged_trials'])}/{int(row['trials'])}"
            + f" & {1.0e6 * float(row['mean_x_rmse_m']):.4g}"
            + f" & {1.0e6 * float(row['max_trial_x_rmse_m']):.4g}"
            + f" & {1.0e6 * float(row['mean_y_rmse_m']):.4g}"
            + f" & {1.0e6 * float(row['max_trial_y_rmse_m']):.4g} \\\\"
        )
    lines.extend(
        [
            r"    \bottomrule",
            r"  \end{tabular}",
            r"  \label{tab:latest-cesr-rho-response}",
            r"\end{table*}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def load_converter():
    spec = importlib.util.spec_from_file_location("project_svg_to_pdf", PAPER_TOOLS)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load SVG converter: {PAPER_TOOLS}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def convert_figures(figures: Path) -> None:
    converter = load_converter()
    for stem in (
        "scibmad_orbit_response_error_stacked",
        "scibmad_orbit_response_error_rho2_normalized_stacked",
    ):
        source = figures / f"{stem}.svg"
        destination = figures / f"{stem}.pdf"
        if not source.is_file():
            raise RuntimeError(f"Missing figure: {source}")
        converter.convert(source, destination)


def first_crossing(rows: list[dict[str, str]], key: str, threshold_um: float = 1.0) -> float | None:
    candidates = sorted(
        (
            (float(row["rho"]), 1.0e6 * float(row[key]))
            for row in rows
            if row["scenario"] == "all" and float(row["rho"]) > 0
        ),
        key=lambda item: item[0],
    )
    return next((rho for rho, value in candidates if value >= threshold_um), None)


def write_report(
    path: Path,
    rows: list[dict[str, str]],
    metadata: dict[str, object],
    trial_rows: list[dict[str, str]],
) -> None:
    failed = int(metadata.get("failed_count", 0))
    total = int(metadata["total_unique_samples"])
    converged = int(metadata["converged_count"])
    fallback = int(metadata.get("fallback_count", 0))
    max_closure = max(float(row["max_closure_norm"]) for row in rows)
    x_cross = first_crossing(rows, "mean_x_rmse_m")
    y_cross = first_crossing(rows, "mean_y_rmse_m")
    incomplete = [
        row for row in trial_rows if str(row.get("converged", "")).casefold() != "true"
    ]
    lines = [
        "# Chapter 1 results: nonlinear closed-orbit response-radius sweep",
        "",
        "Status: `production_complete`.",
        "",
        "## Provenance",
        "",
        f"- Ring: `{metadata['ring_id']}`; RF on; branch `{metadata['branch']}`.",
        f"- Lattice: `{str(metadata['lattice_path']).replace(chr(92), '/')}`.",
        f"- Engine: `{metadata['engine']} {metadata['scibmad_version']}`.",
        f"- Response: `{metadata['response_method']}`; `{metadata['control_count']}` controls and `{metadata['observable_count']}` ordered detector observables.",
        f"- Directions: `{metadata['trials_per_positive_rho_scenario']}` per positive radius and input scenario; seed `{metadata['seed']}`; `rho=1` is `{1.0e6 * float(metadata['base_kick_rad']):g}` microradian active-control RMS.",
        "",
        "## Numerical completion",
        "",
        f"- Unique exact states: `{total}`; converged: `{converged}`; failed: `{failed}`.",
        f"- Recorded fallback count across chunks: `{fallback}`.",
        f"- Maximum reported closed-orbit closure norm: `{max_closure:.6e}`.",
        f"- First sampled all-control radius with mean horizontal residual at least 1 micrometre: `{x_cross if x_cross is not None else 'not reached'}`.",
        f"- First sampled all-control radius with mean vertical residual at least 1 micrometre: `{y_cross if y_cross is not None else 'not reached'}`.",
        "",
        "The crossing statements refer to the discrete sampled radii, not interpolated operating limits. Corrector power-supply, aperture, lifetime, and machine-protection constraints are not represented.",
        "",
        "## Paper artifacts",
        "",
        "- `tables/rho_response_paper_radii.csv` and `tables/rho_response_paper_radii.tex` contain the five matched-benchmark radii.",
        "- `figures/scibmad_orbit_response_error_stacked.svg` and `.pdf` reproduce the paper's stacked response-error layout.",
        "- The rho-squared-normalized companion is a model-validity diagnostic and is not required in the main manuscript.",
        "",
        "![Latest-CESR nonlinear orbit-response error](figures/scibmad_orbit_response_error_stacked.svg)",
        "",
    ]
    if incomplete:
        numerical_end = lines.index(
            "The crossing statements refer to the discrete sampled radii, not interpolated operating limits. Corrector power-supply, aperture, lifetime, and machine-protection constraints are not represented."
        )
        details = "; ".join(
            f"sample {row.get('sample_id')} ({row.get('scenario')}, rho={float(row.get('rho', 0)):g}, trial={row.get('trial_id')}, closure={float(row.get('closure_norm', 0)):.6g})"
            for row in incomplete
        )
        lines.insert(
            numerical_end,
            f"- Incomplete exact reference retained in the survival boundary: {details}.",
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = arguments()
    rows = read_rows(args.summary)
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    chunks = metadata.get("chunk_metadata")
    if isinstance(chunks, list) and chunks and isinstance(chunks[0], dict):
        for key, value in chunks[0].items():
            metadata.setdefault(key, value)
    validate(rows, metadata)
    selected = selected_rows(rows)
    trial_rows = read_rows(args.summary.with_name("rho_sweep_trial_errors.csv"))
    args.tables_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.tables_dir / "rho_response_paper_radii.csv", selected)
    write_latex(args.tables_dir / "rho_response_paper_radii.tex", selected)
    convert_figures(args.figures_dir)
    write_report(args.report, rows, metadata, trial_rows)
    print(f"Validated and built Chapter 1 paper outputs: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
