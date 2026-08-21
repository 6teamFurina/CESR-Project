#!/usr/bin/env python3
"""Run checkpointed Bmad/Tao RF-off optics on the nonlinear-rho inputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

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
    collect_chromatic_sample,
    detector_inventory,
    write_detector_output,
    write_ring_output,
    write_start_orbits,
)

DEFAULT_INPUTS = (
    ORBIT_CALCULATION_DIR
    / "nonlinear_rho_benchmark"
    / "shared_input"
    / "nonlinear_rho_correctors.csv"
)
DEFAULT_MANIFEST = (
    ORBIT_CALCULATION_DIR
    / "nonlinear_rho_benchmark"
    / "shared_input"
    / "sample_manifest.csv"
)
DEFAULT_LATTICE = ORBIT_DIR / "reference" / "cesr_bmad_compatible.bmad"


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--lattice", type=Path, default=DEFAULT_LATTICE)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results" / "bmad")
    parser.add_argument("--delta-step", type=float, default=1.0e-5)
    parser.add_argument("--max-groups", type=int, default=16)
    return parser


def read_manifest(path: Path) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        row["sample_id"] = int(row["sample_id"])
        row["rho"] = float(row["rho"])
        row["trial_id"] = int(row["trial_id"])
    return rows


def manifest_groups(rows: list[dict[str, object]]) -> list[tuple[str, float, list[int]]]:
    groups: list[tuple[str, float, list[int]]] = []
    keys: list[tuple[str, float]] = [("baseline", 0.0)]
    for row in rows:
        key = (str(row["scenario"]), float(row["rho"]))
        if key not in keys:
            keys.append(key)
    for scenario, rho in keys:
        indices = [
            index
            for index, row in enumerate(rows)
            if row["scenario"] == scenario and row["rho"] == rho
        ]
        groups.append((scenario, rho, indices))
    return groups


def group_label(index: int, scenario: str, rho: float) -> str:
    if scenario == "baseline":
        return f"{index:02d}_baseline"
    return f"{index:02d}_{scenario}_rho_{rho:.2f}".replace(".", "p")


def sanitize_error(exception: Exception) -> str:
    return str(exception).replace("\r", " ").replace("\n", " ")


def main() -> int:
    args = make_parser().parse_args()
    inputs = args.inputs.expanduser().resolve()
    manifest_path = args.manifest.expanduser().resolve()
    lattice = args.lattice.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not math.isfinite(args.delta_step) or args.delta_step <= 0.0:
        raise ValueError("--delta-step must be positive and finite")
    if not 1 <= args.max_groups <= 16:
        raise ValueError("--max-groups must be between 1 and 16")

    sample_ids, control_names, samples = read_samples(inputs)
    manifest = read_manifest(manifest_path)
    if sample_ids != [row["sample_id"] for row in manifest]:
        raise RuntimeError("Input and manifest sample IDs differ")
    groups = manifest_groups(manifest)[: args.max_groups]

    expected_detectors, horizontal, vertical = response_tools.parse_cesr_layout(lattice)
    if control_names != horizontal + vertical:
        raise RuntimeError("Input control order differs from the Bmad lattice order")
    references = variable_references(horizontal, vertical)
    output_dir.mkdir(parents=True, exist_ok=True)
    init_path = output_dir / "tao_bmad_chromatic_optics_rf_off.init"
    init_path.write_text(
        response_tools.build_tao_init(
            lattice,
            expected_detectors,
            horizontal,
            vertical,
            derivative_step=args.delta_step,
        ),
        encoding="utf-8",
    )

    from pytao import Tao

    initialization_start = time.perf_counter()
    tao = Tao(init_file=str(init_path), noplot=True)
    tao.cmds(
        [
            "set global rf_on = F",
            f"set global delta_e_chrom = {args.delta_step:.17g}",
            "set particle_start pz = 0",
        ],
        suppress_lattice_calc=True,
        suppress_plotting=True,
    )
    activate_benchmark_data(tao)
    inventory = detector_inventory(tao, expected_detectors)
    initialization_seconds = time.perf_counter() - initialization_start

    status_path = output_dir / "bmad_sample_status.csv"
    group_path = output_dir / "bmad_group_summary.csv"
    with status_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "sample_id",
                "scenario",
                "rho",
                "trial_id",
                "optics_converged",
                "closure_norm",
                "physics_seconds",
                "error_message",
            ),
        )
        writer.writeheader()
    with group_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "group_index",
                "scenario",
                "rho",
                "samples",
                "optics_converged",
                "physics_seconds",
                "maximum_closure_norm",
            ),
        )
        writer.writeheader()

    total_valid = 0
    total_physics_seconds = 0.0
    selected_fields = None
    derivative_sources = None
    native_derivative_query_errors = None

    for group_index, (scenario, rho, indices) in enumerate(groups):
        results: list[dict[str, object]] = []
        valid_ids: list[int] = []
        sample_seconds: list[float] = []
        status_rows: list[dict[str, object]] = []
        closure_norms: list[float] = []
        group_start = time.perf_counter()

        for index in indices:
            sample_start = time.perf_counter()
            result = None
            error_message = ""
            try:
                apply_sample(tao, references, samples[index])
                result, _, _ = collect_chromatic_sample(tao, inventory, args.delta_step)
                elapsed = time.perf_counter() - sample_start
                closure = float(result["transverse_closure_norm"])
                if not math.isfinite(closure):
                    raise RuntimeError("non-finite transverse closure norm")
                valid_ids.append(sample_ids[index])
                results.append(result)
                sample_seconds.append(elapsed)
                closure_norms.append(closure)
                selected_fields = result["selected_fields"]
                derivative_sources = result["derivative_sources"]
                native_derivative_query_errors = result["native_derivative_query_errors"]
                converged = True
            except Exception as exception:  # PyTao exception types vary by release.
                elapsed = time.perf_counter() - sample_start
                closure = math.inf
                converged = False
                error_message = sanitize_error(exception)
            status_rows.append(
                {
                    "sample_id": sample_ids[index],
                    "scenario": scenario,
                    "rho": f"{rho:.17g}",
                    "trial_id": manifest[index]["trial_id"],
                    "optics_converged": str(converged).lower(),
                    "closure_norm": f"{closure:.17g}",
                    "physics_seconds": f"{elapsed:.17g}",
                    "error_message": error_message,
                }
            )

        group_seconds = time.perf_counter() - group_start
        chunk_dir = output_dir / "chunks" / group_label(group_index, scenario, rho)
        chunk_dir.mkdir(parents=True, exist_ok=True)
        if results:
            write_detector_output(
                chunk_dir / "bmad_detector_chromatic_twiss.csv",
                valid_ids,
                inventory,
                results,
            )
            write_ring_output(
                chunk_dir / "bmad_ring_chromatic_twiss.csv",
                valid_ids,
                results,
                sample_seconds,
            )
            write_start_orbits(
                chunk_dir / "bmad_start_closed_orbits.csv",
                valid_ids,
                results,
            )
        with status_path.open("a", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=status_rows[0].keys())
            writer.writerows(status_rows)
        with group_path.open("a", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=(
                    "group_index",
                    "scenario",
                    "rho",
                    "samples",
                    "optics_converged",
                    "physics_seconds",
                    "maximum_closure_norm",
                ),
            )
            writer.writerow(
                {
                    "group_index": group_index,
                    "scenario": scenario,
                    "rho": f"{rho:.17g}",
                    "samples": len(indices),
                    "optics_converged": len(results),
                    "physics_seconds": f"{group_seconds:.17g}",
                    "maximum_closure_norm": (
                        f"{max(closure_norms):.17g}" if closure_norms else "inf"
                    ),
                }
            )
        total_valid += len(results)
        total_physics_seconds += group_seconds
        print(
            f"Bmad optics {scenario:10s} rho={rho:4.2f}: "
            f"{len(results)}/{len(indices)} in {group_seconds:.3f} s",
            flush=True,
        )

    sample_count = sum(len(indices) for _, _, indices in groups)
    metadata = {
        "format": "cesr-nonlinear-rho-optics-bmad-v1",
        "engine": "Bmad/Tao/PyTao",
        "input_csv": str(inputs),
        "manifest_csv": str(manifest_path),
        "lattice": str(lattice),
        "output_directory": str(output_dir),
        "group_count": len(groups),
        "sample_count": sample_count,
        "valid_count": total_valid,
        "failed_count": sample_count - total_valid,
        "initialization_seconds": initialization_seconds,
        "physics_seconds": total_physics_seconds,
        "delta_step": args.delta_step,
        "rf_mode": "off (coasting)",
        "lattice_recalculations_per_sample": 3,
        "selected_tao_fields": selected_fields,
        "local_derivative_sources": derivative_sources,
        "native_derivative_query_errors": native_derivative_query_errors,
        "tao_version_raw": tao.cmd("show version", raises=False),
        "status_csv": str(status_path),
        "group_summary_csv": str(group_path),
    }
    (output_dir / "bmad_metadata.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return 0 if total_valid == sample_count else 2


if __name__ == "__main__":
    raise SystemExit(main())
