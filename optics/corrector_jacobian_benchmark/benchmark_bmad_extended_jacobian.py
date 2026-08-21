#!/usr/bin/env python3
"""Benchmark Bmad optics Jacobians for correctors plus K1/K2 strengths."""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OPTICS_DIR = HERE.parent
PROJECT_DIR = OPTICS_DIR.parent
ORBIT_DIR = PROJECT_DIR / "orbit"
ORBIT_CALCULATION_DIR = ORBIT_DIR / "Orbit_Calculation"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(OPTICS_DIR))
sys.path.insert(0, str(ORBIT_CALCULATION_DIR))

from audit_bmad_extended_inventory import (  # noqa: E402
    ELEMENT_RE,
    dictionary_float,
)
from benchmark_bmad import (  # noqa: E402
    activate_benchmark_data,
    apply_sample,
    read_samples,
    response_tools,
    variable_references,
)
from benchmark_bmad_corrector_jacobian import (  # noqa: E402
    ALL_COLUMNS,
    central_difference,
    detector_row_labels,
    flatten_state,
    write_matrix,
)
from benchmark_bmad_optics import collect_state, detector_inventory  # noqa: E402


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--inputs",
        type=Path,
        default=ORBIT_CALCULATION_DIR / "inputs" / "cesr_corrector_samples_1000.csv",
    )
    result.add_argument(
        "--lattice",
        type=Path,
        default=ORBIT_DIR / "reference" / "cesr_bmad_compatible.bmad",
    )
    result.add_argument(
        "--case",
        choices=("correctors_quads", "correctors_quads_sextupoles", "both"),
        default="both",
    )
    result.add_argument("--corrector-step", type=float, default=1.0e-6)
    result.add_argument("--quadrupole-step", type=float, default=1.0e-6)
    result.add_argument("--sextupole-step", type=float, default=1.0e-4)
    result.add_argument("--warmup", action=argparse.BooleanOptionalAction, default=True)
    result.add_argument("--output-root", type=Path, default=HERE / "results" / "extended")
    return result


def strength_inventory(tao: Any, lattice: Path) -> list[dict[str, Any]]:
    text = response_tools.strip_bmad_comments(lattice.read_text(encoding="utf-8"))
    candidates = [
        (match.group(1).upper(), match.group(2).lower())
        for match in ELEMENT_RE.finditer(text)
    ]
    result: list[dict[str, Any]] = []
    for name, kind in candidates:
        attribute = "K1" if kind == "quadrupole" else "K2"
        strength = dictionary_float(tao.ele_gen_attribs(name), attribute)
        if strength == 0.0:
            continue
        result.append(
            {"name": name, "kind": kind, "attribute": attribute, "baseline": strength}
        )
    return sorted(result, key=lambda item: (item["kind"], item["name"]))


def parameter_inventory(
    control_names: list[str],
    references: list[str],
    magnets: list[dict[str, Any]],
    case: str,
    corrector_step: float,
    quadrupole_step: float,
    sextupole_step: float,
) -> list[dict[str, Any]]:
    result = [
        {
            "label": f"COR:{name}",
            "family": "corrector",
            "reference": reference,
            "baseline": 0.0,
            "step": corrector_step,
        }
        for name, reference in zip(control_names, references)
    ]
    for magnet in magnets:
        if magnet["kind"] == "sextupole" and case == "correctors_quads":
            continue
        result.append(
            {
                **magnet,
                "label": f"{magnet['attribute']}:{magnet['name']}",
                "family": magnet["kind"],
                "step": quadrupole_step if magnet["kind"] == "quadrupole" else sextupole_step,
            }
        )
    return result


def set_parameter(tao: Any, parameter: dict[str, Any], value: float) -> None:
    if parameter["family"] == "corrector":
        command = f"set var {parameter['reference']}|model = {value:.17g}"
    else:
        command = (
            f"set ele {parameter['name']} {parameter['attribute']} = {value:.17g}"
        )
    tao.cmds([command], suppress_lattice_calc=True, suppress_plotting=True)


def run_case(
    tao: Any,
    inventory: dict[str, Any],
    parameters: list[dict[str, Any]],
    output_dir: Path,
    initialization_seconds: float,
    run_warmup: bool,
    lattice: Path,
) -> dict[str, Any]:
    warmup_seconds = 0.0
    if run_warmup:
        warmup_start = time.perf_counter()
        first = parameters[0]
        set_parameter(tao, first, first["baseline"] + first["step"])
        collect_state(tao, inventory, 0.0)
        set_parameter(tao, first, first["baseline"] - first["step"])
        collect_state(tao, inventory, 0.0)
        set_parameter(tao, first, first["baseline"])
        collect_state(tao, inventory, 0.0)
        warmup_seconds = time.perf_counter() - warmup_start

    row_labels = detector_row_labels(inventory)
    parameter_count = len(parameters)
    detector_jacobian = [[0.0] * parameter_count for _ in row_labels]
    closed_orbit_jacobian = [[0.0] * parameter_count for _ in range(6)]
    ring_jacobian = [[0.0] * parameter_count for _ in range(2)]
    update_seconds = 0.0
    recalc_seconds = 0.0
    query_seconds = 0.0
    seconds_per_parameter: list[float] = []

    physics_start = time.perf_counter()
    for column, parameter in enumerate(parameters):
        parameter_start = time.perf_counter()
        update_start = time.perf_counter()
        set_parameter(tao, parameter, parameter["baseline"] + parameter["step"])
        update_seconds += time.perf_counter() - update_start
        plus, recalc, query = collect_state(tao, inventory, 0.0)
        recalc_seconds += recalc
        query_seconds += query

        update_start = time.perf_counter()
        set_parameter(tao, parameter, parameter["baseline"] - parameter["step"])
        update_seconds += time.perf_counter() - update_start
        minus, recalc, query = collect_state(tao, inventory, 0.0)
        recalc_seconds += recalc
        query_seconds += query

        plus_flat = flatten_state(plus)
        minus_flat = flatten_state(minus)
        for row, (high, low) in enumerate(zip(plus_flat, minus_flat)):
            quantity = ALL_COLUMNS[row % len(ALL_COLUMNS)]
            detector_jacobian[row][column] = central_difference(
                high,
                low,
                parameter["step"],
                quantity in ("phi_1", "phi_2"),
            )
        for coordinate in range(6):
            closed_orbit_jacobian[coordinate][column] = central_difference(
                plus["start_orbit"][coordinate],
                minus["start_orbit"][coordinate],
                parameter["step"],
                False,
            )
        ring_jacobian[0][column] = central_difference(
            plus["Q1"], minus["Q1"], parameter["step"], True
        )
        ring_jacobian[1][column] = central_difference(
            plus["Q2"], minus["Q2"], parameter["step"], True
        )

        update_start = time.perf_counter()
        set_parameter(tao, parameter, parameter["baseline"])
        update_seconds += time.perf_counter() - update_start
        elapsed = time.perf_counter() - parameter_start
        seconds_per_parameter.append(elapsed)
        if (column + 1) % 25 == 0 or column + 1 == parameter_count:
            print(
                f"Bmad extended Jacobian {column + 1}/{parameter_count} "
                f"({parameter['label']}): {elapsed:.3f} s"
            )
    physics_seconds = time.perf_counter() - physics_start

    labels = [parameter["label"] for parameter in parameters]
    output_dir.mkdir(parents=True, exist_ok=True)
    write_start = time.perf_counter()
    write_matrix(output_dir / "bmad_detector_optics_jacobian.csv", row_labels, labels, detector_jacobian)
    write_matrix(
        output_dir / "bmad_closed_orbit_jacobian.csv",
        ["start:x", "start:px", "start:y", "start:py", "start:z", "start:pz"],
        labels,
        closed_orbit_jacobian,
    )
    write_matrix(
        output_dir / "bmad_ring_tune_jacobian.csv",
        ["ring:Q1", "ring:Q2"],
        labels,
        ring_jacobian,
    )
    write_seconds = time.perf_counter() - write_start

    family_counts = {
        family: sum(parameter["family"] == family for parameter in parameters)
        for family in ("corrector", "quadrupole", "sextupole")
    }
    metadata = {
        "format": "cesr-extended-optics-jacobian-v1",
        "engine": "Bmad/Tao/PyTao",
        "method": "symmetric finite differences of exact scalar periodic optics",
        "rf_mode": "off (4D coasting)",
        "parameter_count": parameter_count,
        "family_counts": family_counts,
        "detector_count": 99,
        "detector_quantity_count": len(ALL_COLUMNS),
        "detector_jacobian_shape": [len(row_labels), parameter_count],
        "exact_optics_solve_count": 2 * parameter_count,
        "initialization_seconds": initialization_seconds,
        "warmup_seconds": warmup_seconds,
        "physics_seconds": physics_seconds,
        "seconds_per_parameter": seconds_per_parameter,
        "variable_update_seconds": update_seconds,
        "lattice_recalculation_seconds": recalc_seconds,
        "output_query_seconds": query_seconds,
        "write_seconds": write_seconds,
        "lattice": str(lattice),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "tao_version_raw": tao.cmd("show version", raises=False),
        "parameters": parameters,
    }
    metadata_path = output_dir / "bmad_extended_jacobian_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(f"Bmad P={parameter_count} physics: {physics_seconds:.3f} s")
    return metadata


def main() -> int:
    args = parser().parse_args()
    for value, name in (
        (args.corrector_step, "--corrector-step"),
        (args.quadrupole_step, "--quadrupole-step"),
        (args.sextupole_step, "--sextupole-step"),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite number")

    inputs = args.inputs.expanduser().resolve()
    lattice = args.lattice.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    _, control_names, _ = read_samples(inputs)
    detectors, horizontal, vertical = response_tools.parse_cesr_layout(lattice)
    if control_names != horizontal + vertical:
        raise RuntimeError("Input control order differs from the Bmad lattice order")
    references = variable_references(horizontal, vertical)

    output_root.mkdir(parents=True, exist_ok=True)
    init_path = output_root / "tao_bmad_extended_jacobian_rf_off.init"
    init_path.write_text(
        response_tools.build_tao_init(lattice, detectors, horizontal, vertical, args.corrector_step),
        encoding="utf-8",
    )
    from pytao import Tao

    initialization_start = time.perf_counter()
    tao = Tao(init_file=str(init_path), noplot=True)
    tao.cmds(
        ["set global rf_on = F", "set particle_start pz = 0"],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )
    activate_benchmark_data(tao)
    detector_data = detector_inventory(tao, detectors)
    apply_sample(tao, references, [0.0] * len(references))
    collect_state(tao, detector_data, 0.0)
    magnets = strength_inventory(tao, lattice)
    initialization_seconds = time.perf_counter() - initialization_start
    quad_count = sum(item["kind"] == "quadrupole" for item in magnets)
    sext_count = sum(item["kind"] == "sextupole" for item in magnets)
    if (quad_count, sext_count) != (106, 76):
        raise RuntimeError(
            f"Expected 106 quadrupoles and 76 sextupoles, got {quad_count} and {sext_count}"
        )

    cases = (
        ("correctors_quads", "correctors_quads_sextupoles")
        if args.case == "both"
        else (args.case,)
    )
    for case in cases:
        parameters = parameter_inventory(
            control_names,
            references,
            magnets,
            case,
            args.corrector_step,
            args.quadrupole_step,
            args.sextupole_step,
        )
        run_case(
            tao,
            detector_data,
            parameters,
            output_root / case / "bmad",
            initialization_seconds,
            args.warmup,
            lattice,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
