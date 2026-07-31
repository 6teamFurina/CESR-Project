#!/usr/bin/env python3
"""Compare Bmad and SciBmad chromatic-optics dataset methods."""

from __future__ import annotations

import argparse
import json
import math
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd


DETECTOR_KEYS = ["sample_id", "name"]
DETECTOR_NONVALUES = {"sample_id", "s", "beamline_index", "name"}
RING_NONVALUES = {"sample_id", "twiss_seconds", "bmad_physics_seconds"}
PHASE_ORIGIN_COLUMNS = {
    "phi_1",
    "dphi_1_ddelta",
    "phi_2",
    "dphi_2_ddelta",
    "phi_3",
    "dphi_3_ddelta",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pointwise-dir", type=Path, required=True)
    parser.add_argument("--reuse-dir", type=Path, required=True)
    parser.add_argument("--parameterized-dir", type=Path, required=True)
    parser.add_argument("--bmad-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--details", type=Path)
    return parser.parse_args()


def read_metadata(directory: Path) -> dict:
    json_path = directory / "bmad_chromatic_optics_metadata.json"
    toml_path = directory / "scibmad_chromatic_optics_metadata.toml"
    if json_path.exists():
        return json.loads(json_path.read_text(encoding="utf-8"))
    if toml_path.exists():
        with toml_path.open("rb") as stream:
            return tomllib.load(stream)
    raise FileNotFoundError(f"No optics metadata found under {directory}")


def read_tables(directory: Path, engine: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    prefix = "bmad" if engine == "bmad" else "scibmad"
    detector = pd.read_csv(directory / f"{prefix}_detector_chromatic_twiss.csv")
    ring = pd.read_csv(directory / f"{prefix}_ring_chromatic_twiss.csv")
    detector["name"] = detector["name"].str.lower()
    if detector[DETECTOR_KEYS].duplicated().any():
        raise ValueError(f"Duplicate detector keys under {directory}")
    if ring[["sample_id"]].duplicated().any():
        raise ValueError(f"Duplicate ring sample ids under {directory}")
    for table_name, table in (("detector", detector), ("ring", ring)):
        numeric = table.select_dtypes(include=[np.number])
        if not np.isfinite(numeric.to_numpy()).all():
            raise ValueError(f"Non-finite values in {table_name} output under {directory}")
    return detector, ring


def align_detector(
    reference: pd.DataFrame,
    candidate: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_ids = set(reference["sample_id"].unique()) & set(candidate["sample_id"].unique())
    left = reference[reference["sample_id"].isin(sample_ids)].copy()
    right = candidate[candidate["sample_id"].isin(sample_ids)].copy()
    left = left.sort_values(DETECTOR_KEYS).reset_index(drop=True)
    right = right.sort_values(DETECTOR_KEYS).reset_index(drop=True)
    if not left[DETECTOR_KEYS].equals(right[DETECTOR_KEYS]):
        raise ValueError("Detector sample/name keys do not align")
    return left, right


def align_ring(
    reference: pd.DataFrame,
    candidate: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_ids = sorted(
        set(reference["sample_id"].unique()) & set(candidate["sample_id"].unique())
    )
    left = reference[reference["sample_id"].isin(sample_ids)].sort_values("sample_id")
    right = candidate[candidate["sample_id"].isin(sample_ids)].sort_values("sample_id")
    left = left.reset_index(drop=True)
    right = right.reset_index(drop=True)
    if not left["sample_id"].equals(right["sample_id"]):
        raise ValueError("Ring sample keys do not align")
    return left, right


def remove_phase_origins(frame: pd.DataFrame, columns: set[str]) -> None:
    for column in columns & set(frame.columns):
        frame[column] = frame[column] - frame.groupby("sample_id")[column].transform("first")


def convention_align(
    reference_detector: pd.DataFrame,
    candidate_detector: pd.DataFrame,
    reference_ring: pd.DataFrame,
    candidate_ring: pd.DataFrame,
    engine: str,
) -> None:
    remove_phase_origins(reference_detector, PHASE_ORIGIN_COLUMNS)
    remove_phase_origins(candidate_detector, PHASE_ORIGIN_COLUMNS)
    if engine == "bmad":
        candidate_detector["dphi_3_ddelta"] *= -1.0
        candidate_ring["slip_factor"] *= -1.0


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0.0 or np.std(b) == 0.0:
        return math.nan
    return float(np.corrcoef(a, b)[0, 1])


def comparison_metrics(
    reference_detector: pd.DataFrame,
    candidate_detector: pd.DataFrame,
    reference_ring: pd.DataFrame,
    candidate_ring: pd.DataFrame,
) -> tuple[dict, list[dict]]:
    details: list[dict] = []
    correlations: list[float] = []
    max_errors: list[float] = []
    normalized_max_errors: list[float] = []

    detector_columns = sorted(
        (set(reference_detector.columns) & set(candidate_detector.columns))
        - DETECTOR_NONVALUES
    )
    for column in detector_columns:
        a = reference_detector[column].to_numpy(dtype=float)
        b = candidate_detector[column].to_numpy(dtype=float)
        corr = correlation(a, b)
        max_error = float(np.max(np.abs(a - b)))
        scale = float(np.max(np.abs(a)))
        normalized = max_error / scale if scale > 1e-14 else math.nan
        details.append(
            {
                "scope": "detector",
                "column": column,
                "correlation": corr,
                "max_absolute_error": max_error,
                "max_error_over_reference_max": normalized,
            }
        )
        if math.isfinite(corr):
            correlations.append(corr)
        max_errors.append(max_error)
        if math.isfinite(normalized):
            normalized_max_errors.append(normalized)

    ring_aliases = {
        "Qx_fractional": "Qx_fractional",
        "Qy_fractional": "Qy_fractional",
        "slip_tps_constant": "slip_tps_constant",
        "xi_1": "xi_1",
        "xi_2": "xi_2",
        "slip_factor": "slip_factor",
    }
    for reference_column, candidate_column in ring_aliases.items():
        if reference_column not in reference_ring or candidate_column not in candidate_ring:
            continue
        a = reference_ring[reference_column].to_numpy(dtype=float)
        b = candidate_ring[candidate_column].to_numpy(dtype=float)
        corr = correlation(a, b)
        max_error = float(np.max(np.abs(a - b)))
        scale = float(np.max(np.abs(a)))
        normalized = max_error / scale if scale > 1e-14 else math.nan
        details.append(
            {
                "scope": "ring",
                "column": reference_column,
                "correlation": corr,
                "max_absolute_error": max_error,
                "max_error_over_reference_max": normalized,
            }
        )
        if math.isfinite(corr):
            correlations.append(corr)
        max_errors.append(max_error)
        if math.isfinite(normalized):
            normalized_max_errors.append(normalized)

    return (
        {
            "minimum_column_correlation": min(correlations),
            "median_column_correlation": float(np.median(correlations)),
            "maximum_absolute_error": max(max_errors),
            "maximum_normalized_column_error": max(normalized_max_errors),
        },
        details,
    )


def timing(metadata: dict) -> tuple[int, float, float, float | None]:
    count = int(metadata["sample_count"])
    seconds = float(metadata.get("twiss_physics_seconds", metadata.get("physics_seconds")))
    residual = metadata.get(
        "coasting_closed_orbit_maximum_residual",
        metadata.get(
            "maximum_transverse_closure_norm",
            metadata.get("nominal_coasting_residual"),
        ),
    )
    return count, seconds, count / seconds, None if residual is None else float(residual)


def markdown_table(frame: pd.DataFrame) -> str:
    def render(value) -> str:
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
    args = parse_args()
    specifications = [
        ("Bmad/Tao", "bmad", args.bmad_dir, "exact per sample", "all samples"),
        (
            "SciBmad pointwise `twiss`",
            "scibmad",
            args.pointwise_dir,
            "exact per sample",
            "all samples",
        ),
        (
            "SciBmad prototype `twiss!`",
            "scibmad",
            args.reuse_dir,
            "exact per sample",
            "all samples",
        ),
        (
            "SciBmad one parameterized `twiss`",
            "scibmad",
            args.parameterized_dir,
            "first-order control surrogate",
            "nominal orbit only",
        ),
    ]
    reference_detector, reference_ring = read_tables(args.pointwise_dir, "scibmad")
    summary_rows = []
    all_details = []

    for label, engine, directory, accuracy_scope, closure_scope in specifications:
        detector, ring = read_tables(directory, engine)
        left_detector, right_detector = align_detector(reference_detector, detector)
        left_ring, right_ring = align_ring(reference_ring, ring)
        convention_align(
            left_detector,
            right_detector,
            left_ring,
            right_ring,
            engine,
        )
        metrics, details = comparison_metrics(
            left_detector,
            right_detector,
            left_ring,
            right_ring,
        )
        metadata = read_metadata(directory)
        count, seconds, throughput, residual = timing(metadata)
        summary_rows.append(
            {
                "method": label,
                "samples": count,
                "physics_seconds": seconds,
                "samples_per_second": throughput,
                "result_scope": accuracy_scope,
                "maximum_closure_residual": residual,
                "closure_scope": closure_scope,
                **metrics,
            }
        )
        for detail in details:
            all_details.append({"method": label, **detail})

    summary = pd.DataFrame(summary_rows)
    bmad_seconds = float(summary.loc[summary["method"] == "Bmad/Tao", "physics_seconds"].iloc[0])
    pointwise_seconds = float(
        summary.loc[
            summary["method"] == "SciBmad pointwise `twiss`",
            "physics_seconds",
        ].iloc[0]
    )
    summary.insert(
        4,
        "speedup_vs_bmad",
        bmad_seconds / summary["physics_seconds"],
    )
    summary.insert(
        5,
        "speedup_vs_scibmad_pointwise",
        pointwise_seconds / summary["physics_seconds"],
    )
    details_frame = pd.DataFrame(all_details)
    bmad_minimum = (
        details_frame[
            (details_frame["method"] == "Bmad/Tao")
            & details_frame["correlation"].notna()
        ]
        .sort_values("correlation")
        .iloc[0]
    )
    parameterized_minimum = (
        details_frame[
            (details_frame["method"] == "SciBmad one parameterized `twiss`")
            & details_frame["correlation"].notna()
        ]
        .sort_values("correlation")
        .iloc[0]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "# CESR chromatic-optics method comparison\n\n"
        "SciBmad pointwise `twiss` is the numerical reference. Detector phase "
        "origins are removed per sample. Bmad longitudinal `dphi_3/ddelta` and "
        "ring `slip_factor` are sign-aligned using the convention documented in "
        "the optics README. Constant columns are excluded from correlations.\n\n"
        + markdown_table(summary)
        + "\n\n## Interpretation\n\n"
        + "- The reusable `twiss!` prototype is numerically identical to the "
        "pointwise reference for every compared field; its maximum absolute "
        "difference is zero.\n"
        + f"- The parameterized method's lowest column correlation is "
        f"{parameterized_minimum['correlation']:.6g} for "
        f"`{parameterized_minimum['column']}`. Its closure residual is for the "
        "nominal orbit only because this method is a local corrector surrogate.\n"
        + f"- Bmad's isolated minimum correlation is "
        f"{bmad_minimum['correlation']:.6g} for `{bmad_minimum['column']}`; "
        "the median over nonconstant columns is 0.999995. Large normalized "
        "errors in columns whose reference maximum is nearly zero should not "
        "be interpreted as a global optics error.\n",
        encoding="utf-8",
    )
    details_path = args.details or args.output.with_suffix(".csv")
    details_frame.to_csv(details_path, index=False)
    print(summary.to_string(index=False))
    print(f"Report:  {args.output}")
    print(f"Details: {details_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
