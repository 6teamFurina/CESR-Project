#!/usr/bin/env python3
"""Finite-difference the RF-off CESR optics with respect to 119 correctors.

One persistent Tao instance is kept at the nominal machine.  Each corrector is
changed to +step and -step in turn, so the timed region contains 238 exact
Bmad periodic-optics recalculations and the corresponding detector queries.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any, Sequence

HERE = Path(__file__).resolve().parent
OPTICS_DIR = HERE.parent
PROJECT_DIR = OPTICS_DIR.parent
ORBIT_DIR = PROJECT_DIR / "orbit"
ORBIT_CALCULATION_DIR = ORBIT_DIR / "Orbit_Calculation"
sys.path.insert(0, str(OPTICS_DIR))
sys.path.insert(0, str(ORBIT_CALCULATION_DIR))

from benchmark_bmad import (  # noqa: E402
    activate_benchmark_data,
    apply_sample,
    read_samples,
    response_tools,
    variable_references,
)
from benchmark_bmad_optics import (  # noqa: E402
    ORBIT_COLUMNS,
    TWISS_COLUMNS,
    collect_state,
    detector_inventory,
)

ALL_COLUMNS = (*TWISS_COLUMNS, *ORBIT_COLUMNS)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        type=Path,
        default=ORBIT_CALCULATION_DIR / "inputs" / "cesr_corrector_samples_1000.csv",
    )
    parser.add_argument(
        "--lattice",
        type=Path,
        default=ORBIT_DIR / "reference" / "cesr_bmad_compatible.bmad",
    )
    parser.add_argument("--step", type=float, default=1.0e-6)
    parser.add_argument("--warmup", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results" / "bmad")
    return parser


def set_one_control(tao: Any, reference: str, value: float) -> None:
    tao.cmds(
        [f"set var {reference}|model = {value:.17g}"],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )


def central_difference(high: float, low: float, step: float, phase: bool) -> float:
    difference = high - low
    if phase:
        difference = (difference + 0.5) % 1.0 - 0.5
    return difference / (2.0 * step)


def write_matrix(
    path: Path,
    row_labels: Sequence[str],
    control_names: Sequence[str],
    values: Sequence[Sequence[float]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["observable", *control_names])
        for label, row in zip(row_labels, values):
            writer.writerow([label, *(f"{value:.17g}" for value in row)])


def detector_row_labels(inventory: dict[str, Any]) -> list[str]:
    return [
        f"{name.lower()}:{column}"
        for name in inventory["names"]
        for column in ALL_COLUMNS
    ]


def flatten_state(state: dict[str, Any]) -> list[float]:
    return [
        state["columns"][column][detector]
        for detector in range(99)
        for column in ALL_COLUMNS
    ]


def main() -> int:
    args = make_parser().parse_args()
    if not math.isfinite(args.step) or args.step <= 0.0:
        raise ValueError("--step must be a positive finite number")

    inputs = args.inputs.expanduser().resolve()
    lattice = args.lattice.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    _, control_names, _ = read_samples(inputs)
    expected_detectors, horizontal, vertical = response_tools.parse_cesr_layout(lattice)
    if control_names != horizontal + vertical:
        raise RuntimeError("Input control order differs from the Bmad lattice order")
    references = variable_references(horizontal, vertical)

    output_dir.mkdir(parents=True, exist_ok=True)
    init_path = output_dir / "tao_bmad_corrector_jacobian_rf_off.init"
    init_path.write_text(
        response_tools.build_tao_init(
            lattice,
            expected_detectors,
            horizontal,
            vertical,
            derivative_step=args.step,
        ),
        encoding="utf-8",
    )

    try:
        from pytao import Tao
    except ImportError as exc:
        raise RuntimeError("Run in the Linux Bmad environment containing PyTao") from exc

    initialization_start = time.perf_counter()
    tao = Tao(init_file=str(init_path), noplot=True)
    tao.cmds(
        ["set global rf_on = F", "set particle_start pz = 0"],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )
    activate_benchmark_data(tao)
    inventory = detector_inventory(tao, expected_detectors)
    apply_sample(tao, references, [0.0] * len(references))
    collect_state(tao, inventory, 0.0)
    initialization_seconds = time.perf_counter() - initialization_start

    warmup_seconds = 0.0
    if args.warmup:
        warmup_start = time.perf_counter()
        set_one_control(tao, references[0], args.step)
        collect_state(tao, inventory, 0.0)
        set_one_control(tao, references[0], -args.step)
        collect_state(tao, inventory, 0.0)
        set_one_control(tao, references[0], 0.0)
        collect_state(tao, inventory, 0.0)
        warmup_seconds = time.perf_counter() - warmup_start

    row_labels = detector_row_labels(inventory)
    detector_jacobian = [
        [0.0] * len(control_names) for _ in range(len(row_labels))
    ]
    closed_orbit_jacobian = [[0.0] * len(control_names) for _ in range(6)]
    ring_jacobian = [[0.0] * len(control_names) for _ in range(2)]
    update_seconds = 0.0
    recalc_seconds = 0.0
    query_seconds = 0.0
    seconds_per_control: list[float] = []

    physics_start = time.perf_counter()
    for control, reference in enumerate(references):
        control_start = time.perf_counter()
        update_start = time.perf_counter()
        set_one_control(tao, reference, args.step)
        update_seconds += time.perf_counter() - update_start
        plus, recalc, query = collect_state(tao, inventory, 0.0)
        recalc_seconds += recalc
        query_seconds += query

        update_start = time.perf_counter()
        set_one_control(tao, reference, -args.step)
        update_seconds += time.perf_counter() - update_start
        minus, recalc, query = collect_state(tao, inventory, 0.0)
        recalc_seconds += recalc
        query_seconds += query

        plus_flat = flatten_state(plus)
        minus_flat = flatten_state(minus)
        for row, (high, low) in enumerate(zip(plus_flat, minus_flat)):
            column = ALL_COLUMNS[row % len(ALL_COLUMNS)]
            detector_jacobian[row][control] = central_difference(
                high, low, args.step, column in ("phi_1", "phi_2")
            )
        for coordinate in range(6):
            closed_orbit_jacobian[coordinate][control] = central_difference(
                plus["start_orbit"][coordinate],
                minus["start_orbit"][coordinate],
                args.step,
                False,
            )
        ring_jacobian[0][control] = central_difference(
            plus["Q1"], minus["Q1"], args.step, True
        )
        ring_jacobian[1][control] = central_difference(
            plus["Q2"], minus["Q2"], args.step, True
        )

        update_start = time.perf_counter()
        set_one_control(tao, reference, 0.0)
        update_seconds += time.perf_counter() - update_start
        elapsed = time.perf_counter() - control_start
        seconds_per_control.append(elapsed)
        print(
            f"Bmad corrector Jacobian {control + 1}/{len(references)} "
            f"({control_names[control]}): {elapsed:.3f} s"
        )
    physics_seconds = time.perf_counter() - physics_start

    detector_path = output_dir / "bmad_detector_optics_jacobian.csv"
    orbit_path = output_dir / "bmad_closed_orbit_jacobian_6x119.csv"
    ring_path = output_dir / "bmad_ring_tune_jacobian.csv"
    metadata_path = output_dir / "bmad_corrector_jacobian_metadata.json"
    write_start = time.perf_counter()
    write_matrix(detector_path, row_labels, control_names, detector_jacobian)
    write_matrix(
        orbit_path,
        ["start:x", "start:px", "start:y", "start:py", "start:z", "start:pz"],
        control_names,
        closed_orbit_jacobian,
    )
    write_matrix(ring_path, ["ring:Q1", "ring:Q2"], control_names, ring_jacobian)
    write_seconds = time.perf_counter() - write_start

    metadata = {
        "format": "cesr-corrector-optics-jacobian-v1",
        "engine": "Bmad/Tao/PyTao",
        "method": "119 independent symmetric corrector finite differences",
        "rf_mode": "off (4D coasting)",
        "corrector_step_rad": args.step,
        "control_count": len(control_names),
        "detector_count": len(inventory["names"]),
        "detector_quantity_count": len(ALL_COLUMNS),
        "detector_jacobian_shape": [len(row_labels), len(control_names)],
        "closed_orbit_jacobian_shape": [6, len(control_names)],
        "exact_optics_solve_count": 2 * len(control_names),
        "initialization_seconds": initialization_seconds,
        "warmup_enabled": args.warmup,
        "warmup_seconds": warmup_seconds,
        "physics_seconds": physics_seconds,
        "seconds_per_control": seconds_per_control,
        "variable_update_seconds": update_seconds,
        "lattice_recalculation_seconds": recalc_seconds,
        "output_query_seconds": query_seconds,
        "write_seconds": write_seconds,
        "timed_region": "238 exact Bmad optics recalculations, detector/ring queries, and single-variable updates",
        "input_csv": str(inputs),
        "lattice": str(lattice),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "tao_version_raw": tao.cmd("show version", raises=False),
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(f"Bmad finite-difference Jacobian physics: {physics_seconds:.3f} s")
    print(f"Detector Jacobian: {detector_path}")
    print(f"Closed-orbit Jacobian: {orbit_path}")
    print(f"Metadata: {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
