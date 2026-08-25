#!/usr/bin/env python3
"""Compare two SciBmad nonlinear-rho result directories sample by sample."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_rows(path: Path) -> tuple[list[str], dict[int, list[str]]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        return header, {int(row[0]): row for row in reader}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    reference_header, reference = read_rows(
        args.reference_dir / "scibmad_samples.csv"
    )
    candidate_header, candidate = read_rows(
        args.candidate_dir / "scibmad_samples.csv"
    )
    if reference_header != candidate_header:
        raise RuntimeError("SciBmad observable columns differ")
    if set(reference) != set(candidate):
        raise RuntimeError("SciBmad sample IDs differ")

    maximum = 0.0
    sum_square_difference = 0.0
    sum_square_reference = 0.0
    convergence_mismatches = 0
    value_count = 0
    for sample_id in sorted(reference):
        reference_row = reference[sample_id]
        candidate_row = candidate[sample_id]
        convergence_mismatches += reference_row[1].lower() != candidate_row[1].lower()
        for reference_text, candidate_text in zip(reference_row[2:], candidate_row[2:]):
            reference_value = float(reference_text)
            candidate_value = float(candidate_text)
            difference = candidate_value - reference_value
            maximum = max(maximum, abs(difference))
            sum_square_difference += difference * difference
            sum_square_reference += reference_value * reference_value
            value_count += 1

    _, reference_diagnostics = read_rows(
        args.reference_dir / "scibmad_sample_diagnostics.csv"
    )
    _, candidate_diagnostics = read_rows(
        args.candidate_dir / "scibmad_sample_diagnostics.csv"
    )
    if set(reference_diagnostics) != set(candidate_diagnostics):
        raise RuntimeError("SciBmad diagnostic sample IDs differ")
    iteration_mismatches = 0
    maximum_closure_difference = 0.0
    for sample_id in sorted(reference_diagnostics):
        reference_row = reference_diagnostics[sample_id]
        candidate_row = candidate_diagnostics[sample_id]
        iteration_mismatches += reference_row[5] != candidate_row[5]
        maximum_closure_difference = max(
            maximum_closure_difference,
            abs(float(candidate_row[6]) - float(reference_row[6])),
        )

    result = {
        "sample_count": len(reference),
        "observable_count_per_sample": len(reference_header) - 2,
        "value_count": value_count,
        "convergence_mismatches": convergence_mismatches,
        "iteration_mismatches": iteration_mismatches,
        "maximum_observable_difference_m": maximum,
        "observable_difference_relative_l2": (
            math.sqrt(sum_square_difference / sum_square_reference)
            if sum_square_reference
            else math.nan
        ),
        "maximum_closure_difference": maximum_closure_difference,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
