#!/usr/bin/env python3
"""Generate matched CESR closed-orbit samples with Bmad through PyTao."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent.parent
TEST_CODES = PROJECT_ROOT / "test_codes"
sys.path.insert(0, str(TEST_CODES))

import test_control_response_tao as response_tools  # noqa: E402


def read_samples(path: Path) -> tuple[list[int], list[str], list[list[float]]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    if len(rows) < 2 or rows[0][0] != "sample_id":
        raise RuntimeError(f"Invalid sample CSV: {path}")
    names = rows[0][1:]
    if len(names) != 119 or len(set(names)) != 119:
        raise RuntimeError(f"Expected 119 unique controls, found {len(names)}")

    sample_ids: list[int] = []
    values: list[list[float]] = []
    for row_index, row in enumerate(rows[1:], start=1):
        if len(row) != len(rows[0]):
            raise RuntimeError(f"Input row {row_index} has the wrong width")
        sample_ids.append(int(row[0]))
        values.append([float(value) for value in row[1:]])
    return sample_ids, names, values


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        type=Path,
        default=HERE / "inputs" / "cesr_corrector_samples_1000.csv",
    )
    parser.add_argument(
        "--lattice",
        type=Path,
        default=HERE / "reference" / "cesr_bmad_compatible.bmad",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "results" / "formal_1000" / "bmad" / "bmad_rf_on_samples.csv",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=HERE / "results" / "formal_1000" / "bmad" / "bmad_rf_on_metadata.json",
    )
    parser.add_argument("--mode", choices=("rf_on",), default="rf_on")
    parser.add_argument("--warmup-samples", type=int, default=2)
    return parser


def activate_benchmark_data(tao: Any) -> None:
    tao.cmds(
        [
            "veto var *",
            "veto dat *",
            "use var h_steer[*]",
            "use var v_steer[*]",
            "set dat orbit.x[*]|meas = 0",
            "set dat orbit.y[*]|meas = 0",
            "use dat orbit.x[*]",
            "use dat orbit.y[*]",
        ],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )


def enable_rf(tao: Any) -> None:
    tao.cmds(
        [
            f"set ele {name} is_on = T"
            for name in ("RF_W1", "RF_W2", "RF_E1", "RF_E2")
        ],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )


def variable_references(horizontal: list[str], vertical: list[str]) -> list[str]:
    references = [f"h_steer[{index}]" for index in range(len(horizontal))]
    references += [f"v_steer[{index}]" for index in range(len(vertical))]
    return references


def apply_sample(tao: Any, references: list[str], values: list[float]) -> None:
    commands = [
        f"set var {reference}|model = {value:.17g}"
        for reference, value in zip(references, values)
    ]
    tao.cmds(
        commands,
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )


def read_observables(tao: Any) -> tuple[list[float], bool]:
    rows = response_tools.active_data(tao, "x")
    rows += response_tools.active_data(tao, "y")
    values = [float(row["model_value"]) for row in rows]
    good = all(row.get("good_model", True) for row in rows)
    good = good and all(math.isfinite(value) for value in values)
    return values, good


def write_outputs(
    path: Path,
    sample_ids: list[int],
    labels: list[str],
    observables: list[list[float]],
    converged: list[bool],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["sample_id", "converged", *labels])
        for sample_id, good, values in zip(sample_ids, converged, observables):
            writer.writerow(
                [sample_id, str(good).lower(), *(f"{value:.17g}" for value in values)]
            )


def main() -> int:
    args = make_parser().parse_args()
    inputs = args.inputs.expanduser().resolve()
    lattice = args.lattice.expanduser().resolve()
    output = args.output.expanduser().resolve()
    metadata_path = args.metadata.expanduser().resolve()
    sample_ids, names, samples = read_samples(inputs)
    detectors, horizontal, vertical = response_tools.parse_cesr_layout(lattice)
    expected_names = horizontal + vertical
    if names != expected_names:
        raise RuntimeError("Input control order differs from the Bmad lattice order")
    references = variable_references(horizontal, vertical)

    init_path = output.parent / "tao_dataset_benchmark_rf_on.init"
    init_path.parent.mkdir(parents=True, exist_ok=True)
    init_path.write_text(
        response_tools.build_tao_init(
            lattice,
            detectors,
            horizontal,
            vertical,
            derivative_step=1.0e-6,
        ),
        encoding="utf-8",
    )

    try:
        from pytao import Tao
    except ImportError as exc:
        raise RuntimeError(
            "Run this script in the Linux Bmad environment containing PyTao"
        ) from exc

    initialize_start = time.perf_counter()
    tao = Tao(init_file=str(init_path), noplot=True)
    enable_rf(tao)
    activate_benchmark_data(tao)
    initialization_seconds = time.perf_counter() - initialize_start

    x_rows = response_tools.active_data(tao, "x")
    y_rows = response_tools.active_data(tao, "y")
    labels = [
        f"{row['ele_name'].upper()}:x" for row in x_rows
    ] + [
        f"{row['ele_name'].upper()}:y" for row in y_rows
    ]
    if len(labels) != 198:
        raise RuntimeError(f"Expected 198 observables, found {len(labels)}")

    warmup_count = min(args.warmup_samples, len(samples))
    if warmup_count < 1:
        raise ValueError("--warmup-samples must be positive")
    warmup_start = time.perf_counter()
    for sample in samples[1 : 1 + warmup_count]:
        apply_sample(tao, references, sample)
        read_observables(tao)
    apply_sample(tao, references, samples[0])
    read_observables(tao)
    warmup_seconds = time.perf_counter() - warmup_start

    observables: list[list[float]] = []
    converged: list[bool] = []
    update_seconds = 0.0
    read_seconds = 0.0
    physics_start = time.perf_counter()
    for sample_index, sample in enumerate(samples):
        update_start = time.perf_counter()
        apply_sample(tao, references, sample)
        update_seconds += time.perf_counter() - update_start

        read_start = time.perf_counter()
        values, good = read_observables(tao)
        read_seconds += time.perf_counter() - read_start
        observables.append(values)
        converged.append(good)
        if (sample_index + 1) % 100 == 0:
            print(f"Completed {sample_index + 1}/{len(samples)} samples")
    physics_seconds = time.perf_counter() - physics_start

    write_start = time.perf_counter()
    write_outputs(output, sample_ids, labels, observables, converged)
    write_seconds = time.perf_counter() - write_start

    raw_version = tao.cmd("show version", raises=False)
    metadata = {
        "format": "cesr-dataset-benchmark-v1",
        "engine": "Bmad/Tao/PyTao",
        "mode": args.mode,
        "input_csv": str(inputs),
        "lattice": str(lattice),
        "output_csv": str(output),
        "sample_count": len(samples),
        "control_count": len(names),
        "observable_count": len(labels),
        "detector_count": len(detectors),
        "converged_count": sum(converged),
        "failed_count": len(converged) - sum(converged),
        "initialization_seconds": initialization_seconds,
        "warmup_sample_count": warmup_count,
        "warmup_seconds": warmup_seconds,
        "physics_seconds": physics_seconds,
        "variable_update_and_lattice_seconds": update_seconds,
        "observable_read_seconds": read_seconds,
        "samples_per_second": len(samples) / physics_seconds,
        "write_seconds": write_seconds,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "tao_version_raw": raw_version,
        "execution_model": (
            "one persistent Tao instance; 119 variable commands batched with "
            "suppress_lattice_calc=True; one model recalculation per sample"
        ),
        "timed_region": "variable update + lattice recalculation + observable read",
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"Physics: {physics_seconds:.3f} s, "
        f"{len(samples) / physics_seconds:.3f} samples/s, "
        f"converged {sum(converged)}/{len(converged)}"
    )
    print(f"Output:   {output}")
    print(f"Metadata: {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
