#!/usr/bin/env python3
"""Compare checkpointed SciBmad and Bmad nonlinear-rho optics results."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

PHASE_COLUMNS = {
    "phi_1",
    "dphi_1_ddelta",
    "phi_2",
    "dphi_2_ddelta",
    "phi_3",
    "dphi_3_ddelta",
}
DETECTOR_KEYS = ["sample_id", "name"]
DETECTOR_NONVALUES = {"sample_id", "s", "beamline_index", "name"}
RING_COMMON = {
    "Qx_fractional",
    "Qy_fractional",
    "slip_tps_constant",
    "xi_1",
    "xi_2",
    "slip_factor",
}


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scibmad-dir", type=Path, required=True)
    parser.add_argument("--bmad-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def truth(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def finite_max(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    return float(np.max(values)) if len(values) else math.nan


def convergence_summary(scibmad_dir: Path, bmad_dir: Path) -> pd.DataFrame:
    sci = pd.read_csv(scibmad_dir / "scibmad_sample_status.csv")
    bmad = pd.read_csv(bmad_dir / "bmad_sample_status.csv")
    for column in ("rf_on_converged", "coasting_converged", "twiss_converged"):
        sci[column] = truth(sci[column])
    bmad["optics_converged"] = truth(bmad["optics_converged"])
    merged = sci.merge(
        bmad[["sample_id", "optics_converged", "closure_norm"]],
        on="sample_id",
        how="outer",
        suffixes=("_scibmad", "_bmad"),
        validate="one_to_one",
    )
    merged["paired_converged"] = (
        merged["twiss_converged"].fillna(False)
        & merged["optics_converged"].fillna(False)
    )
    rows = []
    for (scenario, rho), group in merged.groupby(["scenario", "rho"], sort=False):
        rows.append(
            {
                "scenario": scenario,
                "rho": rho,
                "samples": len(group),
                "scibmad_rf_on_converged": int(group["rf_on_converged"].sum()),
                "scibmad_coasting_converged": int(group["coasting_converged"].sum()),
                "scibmad_twiss_converged": int(group["twiss_converged"].sum()),
                "bmad_optics_converged": int(group["optics_converged"].sum()),
                "paired_converged": int(group["paired_converged"].sum()),
                "scibmad_maximum_coasting_closure_norm": finite_max(
                    group["coasting_closure_norm"]
                ),
                "bmad_maximum_transverse_closure_norm": finite_max(
                    group["closure_norm"]
                ),
            }
        )
    return pd.DataFrame(rows)


def category(scope: str, column: str) -> str:
    if scope == "ring":
        return "ring_tune" if column.startswith("Q") else "ring_chromatic"
    if column.startswith("orbit_") or column.startswith("dorbit_"):
        return "orbit_derived"
    if column.startswith("d") and column.endswith("_ddelta"):
        return "chromatic_derivative"
    if column in {"gamma_c", "c11", "c12", "c21", "c22"}:
        return "coupled_optics"
    return "ordinary_twiss"


def remove_phase_origins(frame: pd.DataFrame) -> None:
    for column in PHASE_COLUMNS & set(frame.columns):
        frame[column] -= frame.groupby("sample_id")[column].transform("first")


def align_bmad_conventions(detector: pd.DataFrame, ring: pd.DataFrame) -> None:
    if "dphi_3_ddelta" in detector:
        detector["dphi_3_ddelta"] *= -1.0
    if "slip_factor" in ring:
        ring["slip_factor"] *= -1.0


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or np.std(a) == 0.0 or np.std(b) == 0.0:
        return math.nan
    return float(np.corrcoef(a, b)[0, 1])


def metric_row(
    scenario: str,
    rho: float,
    scope: str,
    column: str,
    reference: np.ndarray,
    candidate: np.ndarray,
    reference_response: np.ndarray | None,
    candidate_response: np.ndarray | None,
) -> dict[str, object]:
    error = candidate - reference
    direct_rmse = float(np.sqrt(np.mean(error * error)))
    direct_scale = float(np.sqrt(np.mean(reference * reference)))
    result: dict[str, object] = {
        "scenario": scenario,
        "rho": rho,
        "scope": scope,
        "category": category(scope, column),
        "column": column,
        "value_count": len(reference),
        "direct_rmse": direct_rmse,
        "direct_reference_rms": direct_scale,
        "direct_nrmse_percent": (
            100.0 * direct_rmse / direct_scale if direct_scale > 1e-14 else math.nan
        ),
        "direct_maximum_absolute_error": float(np.max(np.abs(error))),
        "direct_correlation": correlation(reference, candidate),
        "response_rmse": math.nan,
        "response_reference_rms": math.nan,
        "response_nrmse_percent": math.nan,
        "response_maximum_absolute_error": math.nan,
        "response_correlation": math.nan,
    }
    if reference_response is not None and candidate_response is not None:
        response_error = candidate_response - reference_response
        response_rmse = float(np.sqrt(np.mean(response_error * response_error)))
        response_scale = float(np.sqrt(np.mean(reference_response * reference_response)))
        result.update(
            {
                "response_rmse": response_rmse,
                "response_reference_rms": response_scale,
                "response_nrmse_percent": (
                    100.0 * response_rmse / response_scale
                    if response_scale > 1e-14
                    else math.nan
                ),
                "response_maximum_absolute_error": float(np.max(np.abs(response_error))),
                "response_correlation": correlation(reference_response, candidate_response),
            }
        )
    return result


def read_cell(directory: Path, prefix: str, kind: str) -> pd.DataFrame:
    path = directory / f"{prefix}_{kind}.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def validate_cell_tables(detector: pd.DataFrame, ring: pd.DataFrame, label: str) -> None:
    if detector[DETECTOR_KEYS].duplicated().any():
        raise ValueError(f"duplicate detector keys in {label}")
    counts = detector.groupby("sample_id").size()
    if counts.empty or not counts.eq(99).all():
        raise ValueError(f"detector rows per sample are not exactly 99 in {label}")
    if ring["sample_id"].duplicated().any():
        raise ValueError(f"duplicate ring sample IDs in {label}")
    if set(counts.index) != set(ring["sample_id"]):
        raise ValueError(f"detector and ring sample IDs differ in {label}")
    for table_name, frame in (("detector", detector), ("ring", ring)):
        numeric = frame.select_dtypes(include=[np.number]).to_numpy(dtype=float)
        if not np.isfinite(numeric).all():
            raise ValueError(f"non-finite {table_name} values in {label}")


def cell_metrics(
    sci_dir: Path,
    bmad_dir: Path,
    scenario: str,
    rho: float,
    sci_detector_zero: pd.DataFrame,
    bmad_detector_zero: pd.DataFrame,
    sci_ring_zero: pd.DataFrame,
    bmad_ring_zero: pd.DataFrame,
) -> list[dict[str, object]]:
    sci_detector = read_cell(sci_dir, "scibmad", "detector_chromatic_twiss")
    bmad_detector = read_cell(bmad_dir, "bmad", "detector_chromatic_twiss")
    sci_ring = read_cell(sci_dir, "scibmad", "ring_chromatic_twiss")
    bmad_ring = read_cell(bmad_dir, "bmad", "ring_chromatic_twiss")
    if any(frame.empty for frame in (sci_detector, bmad_detector, sci_ring, bmad_ring)):
        return []
    validate_cell_tables(sci_detector, sci_ring, f"SciBmad {scenario} rho={rho}")
    validate_cell_tables(bmad_detector, bmad_ring, f"Bmad {scenario} rho={rho}")
    sci_detector["name"] = sci_detector["name"].str.lower()
    bmad_detector["name"] = bmad_detector["name"].str.lower()
    remove_phase_origins(sci_detector)
    remove_phase_origins(bmad_detector)
    align_bmad_conventions(bmad_detector, bmad_ring)

    detector = sci_detector.merge(
        bmad_detector,
        on=DETECTOR_KEYS,
        suffixes=("_sci", "_bmad"),
        validate="one_to_one",
    )
    ring = sci_ring.merge(
        bmad_ring,
        on="sample_id",
        suffixes=("_sci", "_bmad"),
        validate="one_to_one",
    )
    detector_columns = sorted(
        (set(sci_detector.columns) & set(bmad_detector.columns)) - DETECTOR_NONVALUES
    )
    ring_columns = sorted(RING_COMMON & set(sci_ring.columns) & set(bmad_ring.columns))
    rows: list[dict[str, object]] = []
    baseline = scenario == "baseline"

    sci_zero_by_name = sci_detector_zero.set_index("name")
    bmad_zero_by_name = bmad_detector_zero.set_index("name")
    for column in detector_columns:
        reference = detector[f"{column}_sci"].to_numpy(dtype=float)
        candidate = detector[f"{column}_bmad"].to_numpy(dtype=float)
        if baseline:
            reference_response = candidate_response = None
        else:
            reference_response = reference - detector["name"].map(
                sci_zero_by_name[column]
            ).to_numpy(dtype=float)
            candidate_response = candidate - detector["name"].map(
                bmad_zero_by_name[column]
            ).to_numpy(dtype=float)
        rows.append(
            metric_row(
                scenario,
                rho,
                "detector",
                column,
                reference,
                candidate,
                reference_response,
                candidate_response,
            )
        )

    for column in ring_columns:
        reference = ring[f"{column}_sci"].to_numpy(dtype=float)
        candidate = ring[f"{column}_bmad"].to_numpy(dtype=float)
        if baseline:
            reference_response = candidate_response = None
        else:
            reference_response = reference - float(sci_ring_zero[column].iloc[0])
            candidate_response = candidate - float(bmad_ring_zero[column].iloc[0])
        rows.append(
            metric_row(
                scenario,
                rho,
                "ring",
                column,
                reference,
                candidate,
                reference_response,
                candidate_response,
            )
        )
    return rows


def markdown(frame: pd.DataFrame) -> str:
    def render(value: object) -> str:
        if pd.isna(value):
            return "--"
        if isinstance(value, (float, np.floating)):
            return f"{value:.6g}"
        return str(value)

    header = "| " + " | ".join(frame.columns) + " |"
    rule = "|" + "|".join("---" for _ in frame.columns) + "|"
    rows = [
        "| " + " | ".join(render(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, rule, *rows])


def main() -> int:
    args = make_parser().parse_args()
    sci_root = args.scibmad_dir.resolve()
    bmad_root = args.bmad_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    convergence = convergence_summary(sci_root, bmad_root)
    convergence.to_csv(output_dir / "convergence_summary.csv", index=False)

    sci_chunks = {path.name: path for path in (sci_root / "chunks").iterdir() if path.is_dir()}
    bmad_chunks = {path.name: path for path in (bmad_root / "chunks").iterdir() if path.is_dir()}
    common_chunks = sorted(set(sci_chunks) & set(bmad_chunks))
    baseline_name = next(name for name in common_chunks if name.endswith("baseline"))
    sci_detector_zero = read_cell(
        sci_chunks[baseline_name], "scibmad", "detector_chromatic_twiss"
    )
    bmad_detector_zero = read_cell(
        bmad_chunks[baseline_name], "bmad", "detector_chromatic_twiss"
    )
    sci_ring_zero = read_cell(sci_chunks[baseline_name], "scibmad", "ring_chromatic_twiss")
    bmad_ring_zero = read_cell(bmad_chunks[baseline_name], "bmad", "ring_chromatic_twiss")
    sci_detector_zero["name"] = sci_detector_zero["name"].str.lower()
    bmad_detector_zero["name"] = bmad_detector_zero["name"].str.lower()
    remove_phase_origins(sci_detector_zero)
    remove_phase_origins(bmad_detector_zero)
    align_bmad_conventions(bmad_detector_zero, bmad_ring_zero)

    metrics: list[dict[str, object]] = []
    for chunk_name in common_chunks:
        row = convergence.iloc[int(chunk_name.split("_", 1)[0])]
        metrics.extend(
            cell_metrics(
                sci_chunks[chunk_name],
                bmad_chunks[chunk_name],
                str(row["scenario"]),
                float(row["rho"]),
                sci_detector_zero,
                bmad_detector_zero,
                sci_ring_zero,
                bmad_ring_zero,
            )
        )
    details = pd.DataFrame(metrics)
    details.to_csv(output_dir / "per_cell_quantity_metrics.csv", index=False)

    nonbaseline = details[details["scenario"] != "baseline"].copy()
    category_rows = []
    for keys, group in nonbaseline.groupby(["scenario", "rho", "category"], sort=False):
        finite = group[np.isfinite(group["response_nrmse_percent"])]
        category_rows.append(
            {
                "scenario": keys[0],
                "rho": keys[1],
                "category": keys[2],
                "quantity_count": len(group),
                "finite_response_scale_count": len(finite),
                "median_response_nrmse_percent": finite["response_nrmse_percent"].median(),
                "p90_response_nrmse_percent": finite["response_nrmse_percent"].quantile(0.9),
                "maximum_response_nrmse_percent": finite["response_nrmse_percent"].max(),
                "minimum_response_correlation": finite["response_correlation"].min(),
            }
        )
    category_summary = pd.DataFrame(category_rows)
    category_summary.to_csv(output_dir / "per_cell_category_summary.csv", index=False)

    overall_rows = []
    for category_name, group in nonbaseline.groupby("category", sort=False):
        finite = group[np.isfinite(group["response_nrmse_percent"])]
        worst = finite.sort_values("response_nrmse_percent", ascending=False).iloc[0]
        overall_rows.append(
            {
                "category": category_name,
                "quantity_cells": len(finite),
                "median_response_nrmse_percent": finite["response_nrmse_percent"].median(),
                "p90_response_nrmse_percent": finite["response_nrmse_percent"].quantile(0.9),
                "maximum_response_nrmse_percent": worst["response_nrmse_percent"],
                "worst_case": f"{worst['scenario']}, rho={worst['rho']:.2f}, {worst['column']}",
                "minimum_response_correlation": finite["response_correlation"].min(),
            }
        )
    overall = pd.DataFrame(overall_rows)
    overall.to_csv(output_dir / "overall_category_summary.csv", index=False)

    report = (
        "# Nonlinear-rho RF-off optics comparison\n\n"
        "All comparisons use paired corrector inputs. Detector phase origins are removed "
        "per sample. Bmad `dphi_3/ddelta` and ring `slip_factor` are sign-aligned to "
        "the SciBmad convention. Response errors subtract each engine's zero-input "
        "baseline before comparison.\n\n"
        "## Convergence\n\n"
        + markdown(convergence)
        + "\n\n## Baseline-subtracted response agreement by category\n\n"
        + markdown(overall)
        + "\n"
    )
    (output_dir / "RESULTS.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
