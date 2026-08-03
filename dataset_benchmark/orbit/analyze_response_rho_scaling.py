#!/usr/bin/env python3
"""Calculate local log slopes and rho-squared-normalized orbit errors."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = arguments()
    with args.summary.expanduser().resolve().open(encoding="utf-8", newline="") as stream:
        source = list(csv.DictReader(stream))
    output_rows: list[dict[str, object]] = []
    for scenario in ("all", "horizontal", "vertical"):
        rows = sorted(
            (row for row in source if row["scenario"] == scenario and float(row["rho"]) > 0),
            key=lambda row: float(row["rho"]),
        )
        for plane in ("x", "y"):
            previous_rho: float | None = None
            previous_error: float | None = None
            for row in rows:
                rho = float(row["rho"])
                error = float(row[f"mean_{plane}_rmse_m"])
                slope = math.nan
                if previous_rho is not None and previous_error is not None:
                    slope = math.log(error / previous_error) / math.log(rho / previous_rho)
                output_rows.append(
                    {
                        "scenario": scenario,
                        "plane": plane,
                        "rho": f"{rho:.17g}",
                        "mean_rmse_m": f"{error:.17g}",
                        "mean_rmse_over_rho_squared_m": f"{error / (rho * rho):.17g}",
                        "local_log_slope": f"{slope:.17g}",
                        "trials": row["trials"],
                        "converged_trials": row["converged_trials"],
                        "complete_reference": str(
                            int(row["trials"]) == int(row["converged_trials"])
                        ).lower(),
                    }
                )
                previous_rho = rho
                previous_error = error
    fields = list(output_rows[0])
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
