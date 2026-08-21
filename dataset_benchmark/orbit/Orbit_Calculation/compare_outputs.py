#!/usr/bin/env python3
"""Compare matched Bmad and SciBmad CESR dataset benchmark outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import tomllib
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read_output(
    path: Path,
) -> tuple[list[str], list[int], list[bool], list[list[float]]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    if len(rows) < 2 or rows[0][:2] != ["sample_id", "converged"]:
        raise RuntimeError(f"Invalid benchmark output: {path}")
    labels = rows[0][2:]
    sample_ids = [int(row[0]) for row in rows[1:]]
    converged = [row[1].lower() == "true" for row in rows[1:]]
    values = [[float(value) for value in row[2:]] for row in rows[1:]]
    return labels, sample_ids, converged, values


def vector_norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def correlation(left: list[float], right: list[float]) -> float:
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    dl = [value - left_mean for value in left]
    dr = [value - right_mean for value in right]
    denominator = vector_norm(dl) * vector_norm(dr)
    return sum(a * b for a, b in zip(dl, dr)) / denominator


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ring",
        choices=("latest", "legacy"),
        default="latest",
        help="Compare latest ring-scoped outputs, or explicitly compare archived legacy outputs",
    )
    parser.add_argument(
        "--bmad",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--scibmad",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--bmad-metadata",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--scibmad-metadata",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
    )
    return parser


def main() -> int:
    args = make_parser().parse_args()
    latest = args.ring == "latest"
    default_dir = HERE / "results" / ("latest_cesr" if latest else "formal_1000")
    defaults = {
        "bmad": (
            HERE / "results" / "latest_cesr" / "bmad_reference" / "bmad_rf_on_samples.csv"
            if latest
            else default_dir / "bmad" / "bmad_rf_on_samples.csv"
        ),
        "scibmad": (
            default_dir
            / "formal_1000" if latest else default_dir
        )
        / "scibmad_response_initial_frozen_fallback_bmad_tolerance"
        / "scibmad_rf_on_samples.csv",
        "bmad_metadata": (
            HERE / "results" / "latest_cesr" / "bmad_reference" / "bmad_rf_on_metadata.json"
            if latest
            else default_dir / "bmad" / "bmad_rf_on_metadata.json"
        ),
        "scibmad_metadata": (
            default_dir
            / "formal_1000" if latest else default_dir
        )
        / "scibmad_response_initial_frozen_fallback_bmad_tolerance"
        / "scibmad_rf_on_metadata.toml",
        "report": HERE / "results" / ("latest_cesr" if latest else "formal_1000") / "bmad_scibmad_cross_machine_comparison.md",
    }
    bmad_path = (args.bmad or defaults["bmad"]).resolve()
    scibmad_path = (args.scibmad or defaults["scibmad"]).resolve()
    bmad_metadata_path = (args.bmad_metadata or defaults["bmad_metadata"]).resolve()
    scibmad_metadata_path = (args.scibmad_metadata or defaults["scibmad_metadata"]).resolve()
    report_path = (args.report or defaults["report"]).resolve()
    labels_b, ids_b, good_b, values_b = read_output(bmad_path)
    labels_s, ids_s, good_s, values_s = read_output(scibmad_path)
    if labels_b != labels_s:
        raise RuntimeError("Observable labels differ")
    if ids_b != ids_s:
        raise RuntimeError("Sample IDs differ")

    usable = [
        index
        for index in range(len(ids_b))
        if good_b[index] and good_s[index]
    ]
    if not usable:
        raise RuntimeError("No jointly converged samples")

    per_sample_max: list[float] = []
    per_sample_relative: list[float] = []
    flat_b: list[float] = []
    flat_s: list[float] = []
    for index in usable:
        reference = values_b[index]
        candidate = values_s[index]
        difference = [sci - bmad for bmad, sci in zip(reference, candidate)]
        per_sample_max.append(max(abs(value) for value in difference))
        denominator = vector_norm(reference)
        per_sample_relative.append(vector_norm(difference) / denominator)
        flat_b.extend(reference)
        flat_s.extend(candidate)

    difference_flat = [sci - bmad for bmad, sci in zip(flat_b, flat_s)]
    bmad_metadata = json.loads(bmad_metadata_path.read_text(encoding="utf-8"))
    with scibmad_metadata_path.open("rb") as stream:
        scibmad_metadata = tomllib.load(stream)
    bmad_rate = float(bmad_metadata["samples_per_second"])
    scibmad_rate = float(scibmad_metadata["samples_per_second"])
    bmad_setup = float(bmad_metadata["initialization_seconds"])
    scibmad_setup = float(scibmad_metadata["model_setup_seconds"])
    bmad_warmup = float(bmad_metadata["warmup_seconds"])
    scibmad_warmup = float(scibmad_metadata["warmup_seconds"])
    bmad_physics = float(bmad_metadata["physics_seconds"])
    scibmad_physics = float(scibmad_metadata["physics_seconds"])
    bmad_write = float(bmad_metadata["write_seconds"])
    scibmad_write = float(scibmad_metadata["write_seconds"])
    bmad_total = bmad_setup + bmad_warmup + bmad_physics + bmad_write
    scibmad_total = scibmad_setup + scibmad_warmup + scibmad_physics + scibmad_write

    report = f"""# CESR matched-dataset Bmad-SciBmad benchmark

## Reproducibility

- Matched samples: `{len(ids_b)}`
- Jointly converged samples: `{len(usable)}`
- Controls per sample: `{bmad_metadata["control_count"]}`
- Observables per sample: `{len(labels_b)}`
- Bmad timed region: `{bmad_metadata["timed_region"]}`
- SciBmad timed region: `{scibmad_metadata["timed_region"]}`

## Numerical agreement

| Metric | Value |
|---|---:|
| Global RMSE (m) | `{math.sqrt(statistics.fmean(value * value for value in difference_flat)):.9e}` |
| Global maximum absolute difference (m) | `{max(abs(value) for value in difference_flat):.9e}` |
| Global correlation | `{correlation(flat_b, flat_s):.12f}` |
| Median per-sample relative 2-norm difference | `{statistics.median(per_sample_relative):.9e}` |
| Maximum per-sample relative 2-norm difference | `{max(per_sample_relative):.9e}` |

## Timing

| Engine | Init/setup (s) | Warmup (s) | Physics (s) | Write (s) | Recorded total (s) | Physics samples/s |
|---|---:|---:|---:|---:|---:|---:|
| Bmad/Tao/PyTao | `{bmad_setup:.6f}` | `{bmad_warmup:.6f}` | `{bmad_physics:.6f}` | `{bmad_write:.6f}` | `{bmad_total:.6f}` | `{bmad_rate:.6f}` |
| SciBmad | `{scibmad_setup:.6f}` | `{scibmad_warmup:.6f}` | `{scibmad_physics:.6f}` | `{scibmad_write:.6f}` | `{scibmad_total:.6f}` | `{scibmad_rate:.6f}` |

Steady-state SciBmad/Bmad physics-throughput ratio: `{scibmad_rate / bmad_rate:.6f}`.
The recorded one-shot totals include initialization or model setup, warmup or
compilation, the timed physics region, and file writing.

The throughput ratio is meaningful only when both runs use the same hardware,
CPU thread count, convergence tolerances, input file, and output schema. The
present Bmad and SciBmad results were measured on different machines and must
not be presented as a controlled same-hardware speedup.
"""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
    print(report)
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
