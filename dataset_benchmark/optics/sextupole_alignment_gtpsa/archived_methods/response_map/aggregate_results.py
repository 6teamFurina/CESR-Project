#!/usr/bin/env python3
"""Validate, merge, and summarize the four sextupole-alignment GTPSA parts."""

from __future__ import annotations

import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results" / "full"
COEFFICIENT_COLUMNS = (
    "value",
    "d_k2",
    "d_x",
    "d_y",
    "d2_k2_k2",
    "d2_k2_x",
    "d2_k2_y",
    "d2_x_x",
    "d2_x_y",
    "d2_y_y",
)


def family(observable: str) -> str:
    if observable.startswith("orbit_"):
        return "orbit"
    if observable.startswith("phi_"):
        return "phase"
    if observable.startswith("beta_"):
        return "beta"
    if observable.startswith("alpha_"):
        return "alpha"
    if observable.startswith("c"):
        return "coupling"
    if observable.startswith("tune_"):
        return "tune"
    raise ValueError(f"Unknown observable family: {observable}")


def matrix_metrics(sum_xx: float, sum_xy: float, sum_yy: float) -> dict[str, float]:
    norm_x = math.sqrt(sum_xx)
    norm_y = math.sqrt(sum_yy)
    cosine = sum_xy / (norm_x * norm_y) if norm_x and norm_y else math.nan
    trace = sum_xx + sum_yy
    discriminant = math.sqrt(max(0.0, (sum_xx - sum_yy) ** 2 + 4 * sum_xy**2))
    eigen_max = max(0.0, 0.5 * (trace + discriminant))
    eigen_min = max(0.0, 0.5 * (trace - discriminant))
    singular_max = math.sqrt(eigen_max)
    singular_min = math.sqrt(eigen_min)
    condition = singular_max / singular_min if singular_min else math.inf
    return {
        "norm_k2_x": norm_x,
        "norm_k2_y": norm_y,
        "cosine_x_y": cosine,
        "singular_max": singular_max,
        "singular_min": singular_min,
        "condition_number": condition,
    }


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def main() -> None:
    coefficient_parts = [RESULTS / f"alignment_coefficients_part_{index}.csv" for index in range(1, 5)]
    timing_parts = [RESULTS / f"alignment_timings_part_{index}.csv" for index in range(1, 5)]
    for path in coefficient_parts + timing_parts:
        if not path.is_file():
            raise FileNotFoundError(path)

    combined_path = RESULTS / "alignment_coefficients.csv"
    counts: Counter[str] = Counter()
    observables: Counter[str] = Counter()
    gram: dict[tuple[str, str], list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
    expected_header: list[str] | None = None
    total_rows = 0
    with combined_path.open("w", newline="", encoding="utf-8") as output:
        writer = None
        for part_path in coefficient_parts:
            with part_path.open(newline="", encoding="utf-8") as source:
                reader = csv.DictReader(source)
                if expected_header is None:
                    expected_header = list(reader.fieldnames or [])
                    writer = csv.DictWriter(output, fieldnames=expected_header)
                    writer.writeheader()
                elif list(reader.fieldnames or []) != expected_header:
                    raise ValueError(f"Header mismatch in {part_path}")
                assert writer is not None
                for row in reader:
                    for column in COEFFICIENT_COLUMNS:
                        value = float(row[column])
                        if not math.isfinite(value):
                            raise ValueError(f"Non-finite {column} in {part_path}: {row}")
                    writer.writerow(row)
                    total_rows += 1
                    name = row["sextupole"]
                    observable = row["observable"]
                    group = family(observable)
                    x_value = float(row["d2_k2_x"])
                    y_value = float(row["d2_k2_y"])
                    accumulator = gram[(name, group)]
                    accumulator[0] += x_value * x_value
                    accumulator[1] += x_value * y_value
                    accumulator[2] += y_value * y_value
                    accumulator[3] += 1
                    counts[name] += 1
                    observables[observable] += 1

    if len(counts) != 76:
        raise ValueError(f"Expected 76 sextupoles, found {len(counts)}")
    if set(counts.values()) != {1191}:
        raise ValueError(f"Expected 1191 rows per sextupole, found {sorted(set(counts.values()))}")
    if total_rows != 76 * 1191:
        raise ValueError(f"Expected 90516 rows, found {total_rows}")

    combined_timing_path = RESULTS / "alignment_timings.csv"
    timing_rows: list[dict[str, str]] = []
    timing_header: list[str] | None = None
    with combined_timing_path.open("w", newline="", encoding="utf-8") as output:
        writer = None
        for part_path in timing_parts:
            with part_path.open(newline="", encoding="utf-8") as source:
                reader = csv.DictReader(source)
                if timing_header is None:
                    timing_header = list(reader.fieldnames or [])
                    writer = csv.DictWriter(output, fieldnames=timing_header)
                    writer.writeheader()
                elif list(reader.fieldnames or []) != timing_header:
                    raise ValueError(f"Timing header mismatch in {part_path}")
                assert writer is not None
                for row in reader:
                    timing_rows.append(row)
                    writer.writerow(row)
    if len(timing_rows) != 76:
        raise ValueError(f"Expected 76 timing rows, found {len(timing_rows)}")

    metric_rows: list[dict[str, object]] = []
    for (name, group), (sum_xx, sum_xy, sum_yy, count) in sorted(gram.items()):
        metric_rows.append(
            {
                "sextupole": name,
                "observable_family": group,
                "observation_count": int(count),
                **matrix_metrics(sum_xx, sum_xy, sum_yy),
            }
        )
    metric_path = RESULTS / "mixed_response_identifiability.csv"
    with metric_path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(metric_rows[0]))
        writer.writeheader()
        writer.writerows(metric_rows)

    twiss_seconds = [float(row["twiss_seconds"]) for row in timing_rows]
    family_summary: dict[str, dict[str, float]] = {}
    for group in ("orbit", "phase", "beta", "alpha", "coupling", "tune"):
        selected = [row for row in metric_rows if row["observable_family"] == group]
        conditions = [float(row["condition_number"]) for row in selected if math.isfinite(float(row["condition_number"]))]
        singular_minima = [float(row["singular_min"]) for row in selected]
        family_summary[group] = {
            "condition_p10": percentile(conditions, 0.10),
            "condition_median": statistics.median(conditions),
            "condition_p90": percentile(conditions, 0.90),
            "singular_min_p10": percentile(singular_minima, 0.10),
            "singular_min_median": statistics.median(singular_minima),
            "singular_min_p90": percentile(singular_minima, 0.90),
        }

    summary = {
        "format": "cesr-sextupole-alignment-gtpsa-summary-v1",
        "active_normal_sextupoles": len(counts),
        "coefficient_rows": total_rows,
        "rows_per_sextupole": 1191,
        "all_coefficients_finite": True,
        "detectors": 99,
        "detector_observables": 12,
        "ring_tunes": 3,
        "descriptor": "Descriptor(6, 3, 3, 2)",
        "parameter_order": 2,
        "phase_reference": "DET_00W",
        "twiss_seconds": {
            "sum": sum(twiss_seconds),
            "minimum": min(twiss_seconds),
            "median": statistics.median(twiss_seconds),
            "maximum": max(twiss_seconds),
        },
        "observable_counts": dict(sorted(observables.items())),
        "identifiability_by_observable_family": family_summary,
        "combined_coefficients_csv": str(combined_path),
        "combined_timings_csv": str(combined_timing_path),
        "identifiability_csv": str(metric_path),
    }
    summary_path = RESULTS / "results_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    # The Julia workers wrote their complete CSVs before hitting the old TOML
    # tuple serialization issue. Replace those partial metadata files with a
    # small valid record derived from the completed timing tables.
    for part_index in range(1, 5):
        selected = [row for row in timing_rows if int(row["part_label"]) == part_index]
        metadata_path = RESULTS / f"alignment_metadata_part_{part_index}.toml"
        metadata_path.write_text(
            "\n".join(
                [
                    'format = "cesr-sextupole-alignment-gtpsa-part-v1"',
                    'status = "complete"',
                    'engine = "SciBmad/GTPSA"',
                    'rf_mode = "on (six-dimensional periodic optics)"',
                    'descriptor = "Descriptor(6, 3, 3, 2)"',
                    "parameter_order = 2",
                    f"part_label = {part_index}",
                    f"sextupole_count = {len(selected)}",
                    f"coefficient_rows = {len(selected) * 1191}",
                    f"twiss_seconds = {sum(float(row['twiss_seconds']) for row in selected):.17g}",
                    'phase_reference = "DET_00W"',
                    'hessian_semantics = "true second derivatives; diagonal Taylor coefficient is d2/2"',
                    "",
                ]
            ),
            encoding="utf-8",
        )

    print(f"Merged {total_rows} finite coefficient rows for {len(counts)} sextupoles")
    print(f"Combined coefficients: {combined_path}")
    print(f"Identifiability metrics: {metric_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
