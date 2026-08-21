#!/usr/bin/env python3
"""Compare paired SciBmad and Bmad orbit outputs by scenario and rho."""

from __future__ import annotations

import csv
import argparse
import math
import tomllib
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULT_DIR = HERE / "results"


def ring_paths(ring: str) -> dict[str, Path]:
    artifact = "latest_cesr" if ring == "latest" else "legacy"
    bmad_dir = "bmad_reference" if ring == "latest" else "bmad"
    bmad_name = "bmad_rf_on_samples.csv" if ring == "latest" else "bmad_samples.csv"
    return {
        "scibmad": RESULT_DIR / artifact / "scibmad" / "scibmad_samples.csv",
        "bmad": RESULT_DIR / artifact / bmad_dir / bmad_name,
        "manifest": HERE / "shared_input" / artifact / "sample_manifest.csv",
        "scibmad_timing": RESULT_DIR / artifact / "scibmad" / "scibmad_group_timings.csv",
        "bmad_timing": RESULT_DIR / artifact / bmad_dir / "bmad_group_timings.csv",
        "scibmad_metadata": RESULT_DIR / artifact / "scibmad" / "scibmad_metadata.toml",
        "comparison": RESULT_DIR / artifact / "comparison",
    }


@dataclass
class Accumulator:
    count: int = 0
    sum_square: float = 0.0
    maximum: float = 0.0
    sum_a: float = 0.0
    sum_b: float = 0.0
    sum_aa: float = 0.0
    sum_bb: float = 0.0
    sum_ab: float = 0.0

    def add(self, a: float, b: float) -> None:
        difference = a - b
        self.count += 1
        self.sum_square += difference * difference
        self.maximum = max(self.maximum, abs(difference))
        self.sum_a += a
        self.sum_b += b
        self.sum_aa += a * a
        self.sum_bb += b * b
        self.sum_ab += a * b

    @property
    def rmse(self) -> float:
        return math.sqrt(self.sum_square / self.count) if self.count else math.nan

    @property
    def correlation(self) -> float:
        if not self.count:
            return math.nan
        covariance = self.sum_ab - self.sum_a * self.sum_b / self.count
        variance_a = self.sum_aa - self.sum_a * self.sum_a / self.count
        variance_b = self.sum_bb - self.sum_b * self.sum_b / self.count
        denominator = math.sqrt(max(0.0, variance_a) * max(0.0, variance_b))
        return covariance / denominator if denominator else math.nan

    @property
    def reference_rms(self) -> float:
        return math.sqrt(self.sum_bb / self.count) if self.count else math.nan

    @property
    def relative_rmse(self) -> float:
        return self.rmse / self.reference_rms if self.reference_rms else math.nan


def read_keyed(path: Path) -> tuple[list[str], dict[int, list[str]]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        return header, {int(row[0]): row for row in reader}


def read_manifest(path: Path) -> tuple[list[int], dict[int, tuple[str, float]]]:
    order: list[int] = []
    result: dict[int, tuple[str, float]] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            sample_id = int(row["sample_id"])
            order.append(sample_id)
            result[sample_id] = (row["scenario"], float(row["rho"]))
    return order, result


def read_group_timings(path: Path) -> dict[tuple[str, float], dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return {
            (row["scenario"], float(row["rho"])): row
            for row in csv.DictReader(stream)
        }


def format_number(value: float) -> str:
    return "nan" if math.isnan(value) else f"{value:.8g}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ring",
        choices=("latest", "legacy"),
        default="latest",
        help="Compare ring-scoped latest results, or explicitly compare archived legacy results",
    )
    args = parser.parse_args()
    paths = ring_paths(args.ring)
    sample_order, manifest = read_manifest(paths["manifest"])
    sci_header, sci_rows = read_keyed(paths["scibmad"])
    bmad_header, bmad_rows = read_keyed(paths["bmad"])
    if sci_header != bmad_header:
        raise RuntimeError("SciBmad and Bmad output columns differ")
    if set(sci_rows) != set(sample_order) or set(bmad_rows) != set(sample_order):
        raise RuntimeError("One or both engines did not emit every shared sample")

    x_count = sum(label.endswith(":x") for label in sci_header[2:])
    if x_count * 2 != len(sci_header) - 2:
        raise RuntimeError("Expected matching x/y observable blocks")

    baseline_ids = [sample_id for sample_id in sample_order if manifest[sample_id][0] == "baseline"]
    if len(baseline_ids) != 1:
        raise RuntimeError("Expected exactly one shared baseline")
    baseline_id = baseline_ids[0]
    if sci_rows[baseline_id][1].lower() != "true" or bmad_rows[baseline_id][1].lower() != "true":
        raise RuntimeError("The baseline must converge in both engines")
    sci_baseline = [float(value) for value in sci_rows[baseline_id][2:]]
    bmad_baseline = [float(value) for value in bmad_rows[baseline_id][2:]]
    baseline_x = Accumulator()
    baseline_y = Accumulator()
    for index, (a, b) in enumerate(zip(sci_baseline, bmad_baseline)):
        (baseline_x if index < x_count else baseline_y).add(a, b)

    groups: dict[tuple[str, float], dict[str, object]] = {}
    for sample_id in sample_order:
        scenario, rho = manifest[sample_id]
        if scenario == "baseline":
            continue
        key = (scenario, rho)
        group = groups.setdefault(
            key,
            {
                "samples": 0,
                "sci_converged": 0,
                "bmad_converged": 0,
                "paired": 0,
                "x": Accumulator(),
                "y": Accumulator(),
                "delta_x": Accumulator(),
                "delta_y": Accumulator(),
            },
        )
        group["samples"] = int(group["samples"]) + 1
        sci = sci_rows[sample_id]
        bmad = bmad_rows[sample_id]
        sci_good = sci[1].lower() == "true"
        bmad_good = bmad[1].lower() == "true"
        group["sci_converged"] = int(group["sci_converged"]) + int(sci_good)
        group["bmad_converged"] = int(group["bmad_converged"]) + int(bmad_good)
        if not (sci_good and bmad_good):
            continue
        group["paired"] = int(group["paired"]) + 1
        sci_values = [float(value) for value in sci[2:]]
        bmad_values = [float(value) for value in bmad[2:]]
        x_accumulator = group["x"]
        y_accumulator = group["y"]
        delta_x_accumulator = group["delta_x"]
        delta_y_accumulator = group["delta_y"]
        assert isinstance(x_accumulator, Accumulator)
        assert isinstance(y_accumulator, Accumulator)
        assert isinstance(delta_x_accumulator, Accumulator)
        assert isinstance(delta_y_accumulator, Accumulator)
        for index, (a, b) in enumerate(zip(sci_values[:x_count], bmad_values[:x_count])):
            x_accumulator.add(a, b)
            delta_x_accumulator.add(
                a - sci_baseline[index], b - bmad_baseline[index]
            )
        for local_index, (a, b) in enumerate(
            zip(sci_values[x_count:], bmad_values[x_count:])
        ):
            y_accumulator.add(a, b)
            index = x_count + local_index
            delta_y_accumulator.add(
                a - sci_baseline[index], b - bmad_baseline[index]
            )

    sci_timing = read_group_timings(paths["scibmad_timing"])
    bmad_timing = read_group_timings(paths["bmad_timing"])
    with paths["scibmad_metadata"].open("rb") as stream:
        sci_metadata = tomllib.load(stream)
    output_rows: list[dict[str, str | int]] = []
    for key, group in groups.items():
        x_result = group["x"]
        y_result = group["y"]
        delta_x_result = group["delta_x"]
        delta_y_result = group["delta_y"]
        assert isinstance(x_result, Accumulator)
        assert isinstance(y_result, Accumulator)
        assert isinstance(delta_x_result, Accumulator)
        assert isinstance(delta_y_result, Accumulator)
        sci_seconds = float(sci_timing[key]["physics_seconds"])
        sci_total_seconds = sci_seconds + float(sci_timing[key]["model_setup_seconds"])
        bmad_seconds = float(bmad_timing[key]["physics_seconds"])
        samples = int(group["samples"])
        output_rows.append(
            {
                "scenario": key[0],
                "rho": format_number(key[1]),
                "samples": samples,
                "scibmad_converged": int(group["sci_converged"]),
                "bmad_converged": int(group["bmad_converged"]),
                "paired_converged": int(group["paired"]),
                "x_rmse_m": format_number(x_result.rmse),
                "x_max_abs_m": format_number(x_result.maximum),
                "x_correlation": format_number(x_result.correlation),
                "y_rmse_m": format_number(y_result.rmse),
                "y_max_abs_m": format_number(y_result.maximum),
                "y_correlation": format_number(y_result.correlation),
                "baseline_subtracted_x_rmse_m": format_number(delta_x_result.rmse),
                "baseline_subtracted_x_max_abs_m": format_number(delta_x_result.maximum),
                "baseline_subtracted_x_correlation": format_number(delta_x_result.correlation),
                "baseline_subtracted_bmad_x_rms_m": format_number(delta_x_result.reference_rms),
                "baseline_subtracted_x_relative_rmse": format_number(delta_x_result.relative_rmse),
                "baseline_subtracted_y_rmse_m": format_number(delta_y_result.rmse),
                "baseline_subtracted_y_max_abs_m": format_number(delta_y_result.maximum),
                "baseline_subtracted_y_correlation": format_number(delta_y_result.correlation),
                "baseline_subtracted_bmad_y_rms_m": format_number(delta_y_result.reference_rms),
                "baseline_subtracted_y_relative_rmse": format_number(delta_y_result.relative_rmse),
                "scibmad_physics_seconds": format_number(sci_seconds),
                "scibmad_setup_plus_physics_seconds": format_number(sci_total_seconds),
                "bmad_physics_seconds": format_number(bmad_seconds),
                "scibmad_physics_samples_per_second": format_number(samples / sci_seconds),
                "scibmad_end_to_end_samples_per_second": format_number(samples / sci_total_seconds),
                "bmad_samples_per_second": format_number(samples / bmad_seconds),
                "speedup_physics_only": format_number(bmad_seconds / sci_seconds),
                "speedup_including_scibmad_setup": format_number(bmad_seconds / sci_total_seconds),
            }
        )

    comparison_dir = paths["comparison"]
    comparison_dir.mkdir(parents=True, exist_ok=True)
    summary_path = comparison_dir / "comparison_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    total_samples = sum(int(row["samples"]) for row in output_rows)
    sci_seconds = sum(float(row["scibmad_physics_seconds"]) for row in output_rows)
    sci_total_seconds = sum(
        float(row["scibmad_setup_plus_physics_seconds"]) for row in output_rows
    )
    bmad_seconds = sum(float(row["bmad_physics_seconds"]) for row in output_rows)
    sci_converged = sum(int(row["scibmad_converged"]) for row in output_rows)
    bmad_converged = sum(int(row["bmad_converged"]) for row in output_rows)
    paired = sum(int(row["paired_converged"]) for row in output_rows)
    sci_initial_guess_seconds = float(sci_metadata["shared_initial_guess_setup_seconds"])
    sci_all_runtime_seconds = sci_total_seconds + sci_initial_guess_seconds

    report_path = comparison_dir / "RESULTS.md"
    lines = [
        "# Nonlinear-rho orbit benchmark results",
        "",
        f"The comparison contains {total_samples:,} nonzero shared inputs "
        f"(up to {max(int(row['samples']) for row in output_rows)} for any scenario/rho cell), plus one shared baseline.",
        "",
        "| Metric | SciBmad | Bmad/Tao |",
        "|---|---:|---:|",
        f"| Converged nonzero inputs | {sci_converged}/{total_samples} | {bmad_converged}/{total_samples} |",
        f"| Paired converged inputs | {paired}/{total_samples} | {paired}/{total_samples} |",
        f"| Summed physics time | {sci_seconds:.3f} s | {bmad_seconds:.3f} s |",
        f"| Throughput | {total_samples / sci_seconds:.3f} samples/s | {total_samples / bmad_seconds:.3f} samples/s |",
        f"| SciBmad speedup, physics only | {bmad_seconds / sci_seconds:.3f}x | 1x |",
        f"| SciBmad setup + physics time | {sci_total_seconds:.3f} s | n/a |",
        f"| SciBmad speedup including per-group model setup | {bmad_seconds / sci_total_seconds:.3f}x | 1x |",
        f"| SciBmad initial-guess + model setup + physics time | {sci_all_runtime_seconds:.3f} s | n/a |",
        f"| SciBmad speedup including all runtime setup | {bmad_seconds / sci_all_runtime_seconds:.3f}x | 1x |",
        "",
        "## Zero-input baseline agreement",
        "",
        "| Plane | RMSE [m] | Maximum absolute difference [m] | Correlation |",
        "|---|---:|---:|---:|",
        f"| x | {format_number(baseline_x.rmse)} | {format_number(baseline_x.maximum)} | {format_number(baseline_x.correlation)} |",
        f"| y | {format_number(baseline_y.rmse)} | {format_number(baseline_y.maximum)} | {format_number(baseline_y.correlation)} |",
        "",
        "## Per-cell results",
        "",
        "The response RMSE compares each engine's orbit after subtracting that engine's own zero-input baseline.",
        "",
        "| Scenario | rho | Sci conv. | Bmad conv. | response x RMSE [m] | x relative RMSE | response y RMSE [m] | y relative RMSE | Sci/Bmad speedup (physics) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in output_rows:
        lines.append(
            f"| {row['scenario']} | {row['rho']} | {row['scibmad_converged']}/{row['samples']} | "
            f"{row['bmad_converged']}/{row['samples']} | {row['baseline_subtracted_x_rmse_m']} | "
            f"{100 * float(row['baseline_subtracted_x_relative_rmse']):.4f}% | "
            f"{row['baseline_subtracted_y_rmse_m']} | "
            f"{100 * float(row['baseline_subtracted_y_relative_rmse']):.4f}% | "
            f"{row['speedup_physics_only']}x |"
        )
    lines += [
        "",
        "## Timing interpretation",
        "",
        "SciBmad was run in one Julia process, with the configured simultaneous TPSA lanes per cell. "
        "Its physics-only number includes frozen-Jacobian iterations, explicit closure checks, "
        "tracking, and any full-AD fallback. The end-to-end variant additionally includes "
        "construction of each batch model, but excludes compilation warmup, CSV I/O, and "
        "the one-time first-order initial-guess preparation reported in metadata.",
        "",
        "Bmad was run sequentially in one persistent Tao/PyTao process in Ubuntu-Bmad. Its "
        "timed region includes corrector updates, one Tao model recalculation, and observable "
        "reads. Its convergence flag is based on Tao good_model flags and finite data; this "
        "path does not expose the explicit one-turn closure norm used by SciBmad.",
        "",
        "The two engines ran on the same physical machine but in different host runtimes "
        "(SciBmad on Windows Julia and Bmad inside WSL Ubuntu), so timing is an application-level "
        "comparison rather than a microarchitectural kernel benchmark.",
    ]
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
